#!/usr/bin/env python3
"""Insert or update the `intuition:` field in every catalog entry.

    uv run python tools/add_intuition.py <texts.json>
    uv run python tools/add_intuition.py <texts.json> --check

Reads a JSON mapping of chemistry id -> explanation and writes each one into
`catalog/chemistries/<id>.yaml` as a folded block scalar, immediately after the
`name:` line.

The insertion is textual on purpose. Round-tripping the hand-curated YAML files
through a parser would reformat them and drop their comments, so this edits only
the lines it owns and leaves everything else byte-for-byte. It is idempotent:
an existing `intuition:` block is replaced, not duplicated.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CATALOG = REPO / "catalog" / "chemistries"

NAME_RE = re.compile(r"^(\s*)name:\s", re.M)


def strip_existing(lines: list[str], indent: str) -> list[str]:
    """Drop any previous `intuition:` key and its indented continuation lines."""
    out, skipping = [], False
    for line in lines:
        if skipping:
            # Continuation lines are indented deeper than the key, or blank.
            if line.strip() == "" or line.startswith(indent + " "):
                continue
            skipping = False
        if re.match(rf"^{indent}intuition:", line):
            skipping = True
            continue
        out.append(line)
    return out


def render(text: str, indent: str) -> list[str]:
    """A folded block scalar, wrapped so the YAML stays readable.

    A folded scalar turns each line break into a space, so the wrapping must
    never split inside a word. textwrap breaks on hyphens by default, which
    silently turned "lambda-expressions" into "lambda- expressions"; long words
    must not be broken either, for the same reason. Overflowing the width is
    the lesser evil.
    """
    body = " ".join(text.split())
    wrapped = textwrap.wrap(body, width=74, break_on_hyphens=False,
                            break_long_words=False)
    return [f"{indent}intuition: >\n"] + [f"{indent}  {line}\n" for line in wrapped]


def apply(cid: str, text: str) -> str:
    path = CATALOG / f"{cid}.yaml"
    if not path.exists():
        return f"no such entry: {cid}"
    original = path.read_text()
    match = NAME_RE.search(original)
    if not match:
        return f"{cid}: no `name:` line to anchor on"
    indent = match.group(1)

    lines = strip_existing(original.splitlines(keepends=True), indent)
    # Find the anchor again in the stripped copy.
    for i, line in enumerate(lines):
        if re.match(rf"^{indent}name:\s", line):
            break
    else:
        return f"{cid}: anchor vanished after stripping"

    updated = lines[: i + 1] + render(text, indent) + lines[i + 1:]
    path.write_text("".join(updated))
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("texts", type=Path, help="JSON: {chemistry id: explanation}")
    ap.add_argument("--check", action="store_true",
                    help="report coverage without writing anything")
    args = ap.parse_args(argv)

    texts: dict[str, str] = json.loads(args.texts.read_text())

    from chemart.catalog import load
    ids = {c.id for c in load()}

    missing = sorted(ids - set(texts))
    unknown = sorted(set(texts) - ids)
    if unknown:
        print(f"error: {len(unknown)} unknown ids: {unknown[:5]}", file=sys.stderr)
        return 2

    if args.check:
        print(f"{len(texts)}/{len(ids)} entries have an explanation")
        if missing:
            print(f"missing ({len(missing)}): {', '.join(missing)}")
        return 1 if missing else 0

    problems = [p for cid, text in sorted(texts.items()) if (p := apply(cid, text))]
    for p in problems:
        print(f"error: {p}", file=sys.stderr)

    # Validate by reloading: the YAML must still parse and the text must match.
    from chemart import catalog
    catalog.load.cache_clear()
    entries = {c.id: c for c in catalog.load()}
    bad = [cid for cid, text in texts.items()
           if " ".join((entries[cid].intuition or "").split()) != " ".join(text.split())]
    if bad:
        print(f"error: {len(bad)} entries did not round-trip: {bad[:5]}", file=sys.stderr)
        return 2

    print(f"wrote intuition for {len(texts) - len(problems)} entries; "
          f"{len(ids) - len(texts)} still missing")
    return 2 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
