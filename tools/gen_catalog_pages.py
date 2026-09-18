#!/usr/bin/env python3
"""Generate one documentation page per chemistry, straight from the catalog.

    uv run python tools/gen_catalog_pages.py              # with live network summaries
    uv run python tools/gen_catalog_pages.py --fast       # catalog only, no generation

Writing these pages from the library itself means the site cannot drift from
the code: parameters, provenance and capabilities are whatever the catalog says
today, and the "what you get" block is a network actually generated at build
time.

Output: docs/catalog/index.md plus docs/catalog/<id>.md for all 98.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "docs" / "catalog"

import chemart                                    # noqa: E402
from chemart.catalog import load                   # noqa: E402

FIDELITY_NOTE = {
    "book": "implemented exactly as the book specifies",
    "book+decisions": "the book left gaps; each filled choice is listed below",
    "reconstructed": "built from the original papers listed below",
}


def esc(text) -> str:
    """Collapse whitespace and neutralise pipes so text survives a table cell."""
    return " ".join(str(text or "").split()).replace("|", "\\|")


def block(text) -> str:
    """Collapse whitespace but keep it prose, for paragraphs."""
    return " ".join(str(text or "").split())


def param_rows(c) -> list[str]:
    """One row per parameter.

    Bounds, choices and the published range live in the same cell as the
    meaning rather than in columns of their own: an enum with fourteen choices
    is 500+ characters, and as a separate column it forced the table wider than
    the page. Folding them in gives the long text room to wrap.
    """
    rows = []
    for p in c.params:
        # Some defaults are whole rule sets; the full value is in the YAML and
        # in describe_chemistry, so the table only needs to be recognisable.
        default = esc(p.default)
        if len(default) > 48:
            default = default[:47] + "…"

        notes = []
        if p.choices:
            # A long enum inlined in full is the single worst offender for cell
            # height (l-systems lists 14 named systems, ~550 characters). Show
            # enough to recognise the shape and point at the authoritative list.
            shown = [f"`{v}`" for v in p.choices[:8]]
            more = len(p.choices) - len(shown)
            notes.append("one of " + ", ".join(shown)
                         + (f" … and {more} more ({len(p.choices)} total)" if more else ""))
        if p.min is not None and p.max is not None:
            notes.append(f"`{p.min}` … `{p.max}`")
        elif p.min is not None:
            notes.append(f"≥ `{p.min}`")
        elif p.max is not None:
            notes.append(f"≤ `{p.max}`")
        if p.range:
            notes.append(f"*range:* {esc(p.range)}")

        what = esc(p.meaning) or "—"
        if notes:
            what += " <br>" + " · ".join(notes)
        rows.append(f"| `{p.name}` | `{p.type}` | `{default}` | {p.role} | {what} |")
    return rows


def page(c, net, error) -> str:
    L: list[str] = []
    L.append(f"# {c.name}")
    L.append("")
    L.append(f"`{c.id}`{' · *' + esc(c.origin) + '*' if c.origin else ''}")
    L.append("")
    if c.aliases:
        L.append("*Also known as:* " + ", ".join(f"*{esc(a)}*" for a in c.aliases))
        L.append("")

    # The explanation comes first, before any table: someone meeting this
    # chemistry needs to know what the idea is before the (S, R, A) detail
    # means anything.
    if c.intuition:
        L.append(block(c.intuition))
        L.append("")

    L.append("| | |")
    L.append("|---|---|")
    L.append(f"| **family** | {c.family} |")
    L.append(f"| **kind** | {c.kind} |")
    L.append(f"| **constructive** | {'yes — the species set grows at run time' if c.constructive else 'no — fixed species set'} |")
    L.append(f"| **fidelity** | `{c.fidelity}` — {FIDELITY_NOTE.get(c.fidelity, '')} |")
    L.append(f"| **book** | {esc(c.book)} |")
    if c.refs:
        L.append(f"| **refs** | {esc(', '.join(c.refs))} |")
    L.append(f"| **provides** | {', '.join(f'`{p}`' for p in c.provides) or '—'} |")
    L.append("")

    # --- the (S, R, A) triple
    L.append("## Molecules, reactions, reactor")
    L.append("")
    L.append(f"**S — molecules** ({c.S.get('definition', '?')}): {block(c.S.get('repr'))}")
    L.append("")
    arity = c.R.get("arity", "?")
    L.append(f"**R — reactions** ({c.R.get('definition', '?')}, arity {arity}): {block(c.R.get('scheme'))}")
    L.append("")
    reactor = c.A.get("reactor", "?")
    if isinstance(reactor, list):
        reactor = ", ".join(map(str, reactor))
    L.append(f"**A — reactor**: {reactor}")
    if c.A.get("dilution"):
        L.append(f" · *dilution:* {block(c.A.get('dilution'))}")
    L.append("")

    # --- what you actually get
    L.append("## What you get")
    L.append("")
    if error:
        L.append("!!! failure \"This page was built without generating a network\"")
        L.append(f"    `{error}`")
        L.append("")
    elif net is not None:
        L.append("```python")
        L.append(f'net = chemart.generate_network("{c.id}", seed=1)')
        L.append("```")
        L.append("")
        L.append("```")
        L.append(net.summary())
        L.append("```")
        L.append("")
        shown = net.reactions[:8]
        if shown:
            L.append("First reactions:")
            L.append("")
            L.append("```")
            for r in shown:
                L.append(r.to_text())
            if len(net.reactions) > len(shown):
                L.append(f"… and {len(net.reactions) - len(shown)} more")
            L.append("```")
            L.append("")
        else:
            L.append("This chemistry defines no transformational reactions; the law is in "
                     "`extras['interaction_law']`.")
            L.append("")

    # --- parameters
    L.append("## Parameters")
    L.append("")
    if c.params:
        L.append("| name | type | default | role | what it does |")
        L.append("|---|---|---|---|---|")
        L += param_rows(c)
    else:
        L.append("This chemistry takes no parameters.")
    L.append("")

    # --- phenomena
    if c.phenomena:
        L.append("## Published phenomena")
        L.append("")
        L.append("What the literature reports this model produces. Whether the "
                 "generator reproduces each one is recorded in the decisions below.")
        L.append("")
        for p in c.phenomena:
            L.append(f"- {block(p)}")
        L.append("")

    # --- provenance
    if c.sources:
        L.append("## Sources")
        L.append("")
        for s in c.sources:
            L.append(f"- {block(s)}")
        L.append("")

    if c.decisions:
        L.append("## Decisions")
        L.append("")
        L.append("Every gap, ambiguity or erratum in the sources, and how Chemart "
                 "resolved it. Read this before quoting a number from this entry.")
        L.append("")
        for d in c.decisions:
            L.append(f"- {block(d)}")
        L.append("")

    if c.notes:
        L.append("## Notes")
        L.append("")
        L.append(block(c.notes))
        L.append("")

    L.append("---")
    L.append("")
    L.append(f"*Specification: `catalog/chemistries/{c.id}.yaml` · "
             f"generator: `chemart/chemistries/{c.id.replace('-', '_')}.py` · "
             f"tests: `tests/chemistries/test_{c.id.replace('-', '_')}.py`*")
    L.append("")
    return "\n".join(L)


#: Capabilities worth putting in a table. `topology` and `stoichiometry` are on
#: all 98 entries and `initial-state` on 84, so listing them in every row costs
#: width and tells you nothing; the full set is on each chemistry's own page.
TABLE_CAPS = [
    ("rate-constants", "kinetics"),
    ("rate-law", "rate law"),
    ("energies", "energies"),
    ("thermodynamic-consistency", "thermo"),
    ("mass-conservation", "conservation"),
    ("flow", "flow"),
    ("space", "space"),
    ("compartments", "compartments"),
]


def short_origin(origin, limit: int = 26) -> str:
    """First attribution only, truncated — the full string is on the page."""
    text = block(origin)
    if not text:
        return "—"
    text = text.split(";")[0].strip()
    if len(text) > limit:
        text = text[: limit - 1].rstrip(" ,&") + "…"
    return esc(text)


def table_caps(c) -> str:
    have = set(c.provides)
    labels = [label for tag, label in TABLE_CAPS if tag in have]
    return ", ".join(labels) if labels else "—"


def index_page(entries, nets) -> str:
    L: list[str] = []
    L.append("# Catalog")
    L.append("")
    L.append(f"All **{len(entries)}** chemistries, one page each, generated from the "
             "catalog so these pages cannot drift from the library.")
    L.append("")
    L.append("Columns: **grows** is whether the species set is open and expands at "
             "run time; **S/R** is the species and reaction count at *default* "
             "parameters with `seed=1`, which for most chemistries scales up "
             "considerably; **beyond topology** lists only the capabilities that "
             "distinguish entries, since every one of them supplies topology and "
             "stoichiometry. **fidelity** says how close the implementation is to a "
             "published specification — see [Fidelity and trust](../trust.md).")
    L.append("")
    L.append("Each chemistry's own page has the full capability list, the complete "
             "attribution, its parameters and its provenance.")
    L.append("")

    by_family: dict[str, list] = {}
    for c in entries:
        by_family.setdefault(c.family, []).append(c)

    for family in sorted(by_family):
        group = sorted(by_family[family], key=lambda c: c.id)
        L.append(f"## {family}")
        L.append("")
        L.append("| chemistry | origin | grows | fidelity | S/R | beyond topology |")
        L.append("|---|---|:--:|---|--:|---|")
        for c in group:
            net = nets.get(c.id)
            size = f"{len(net.species)}/{len(net.reactions)}" if net else "—"
            L.append(
                f"| [{c.name}]({c.id}.md) | {short_origin(c.origin)} | "
                f"{'yes' if c.constructive else '·'} | {c.fidelity} | {size} | "
                f"{table_caps(c)} |"
            )
        L.append("")

    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fast", action="store_true",
                    help="skip generating networks (much faster, less useful pages)")
    args = ap.parse_args(argv)

    entries = sorted(load(), key=lambda c: c.id)
    OUT.mkdir(parents=True, exist_ok=True)
    for stale in OUT.glob("*.md"):
        stale.unlink()

    nets, errors = {}, {}
    if not args.fast:
        t0 = time.perf_counter()
        for i, c in enumerate(entries, 1):
            print(f"\r[{i:>3}/{len(entries)}] {c.id:<32}", end="", file=sys.stderr, flush=True)
            try:
                nets[c.id] = chemart.generate_network(c.id, seed=1)
            except Exception as err:                        # noqa: BLE001
                errors[c.id] = f"{type(err).__name__}: {err}"
        print(f"\rgenerated {len(nets)}/{len(entries)} networks in "
              f"{time.perf_counter() - t0:.0f}s{' ' * 20}", file=sys.stderr)

    for c in entries:
        (OUT / f"{c.id}.md").write_text(page(c, nets.get(c.id), errors.get(c.id)))
    (OUT / "index.md").write_text(index_page(entries, nets))

    print(f"wrote {len(entries) + 1} pages to {OUT.relative_to(REPO)}")
    if errors:
        print(f"{len(errors)} chemistries failed to generate:", file=sys.stderr)
        for cid, err in errors.items():
            print(f"  {cid}: {err}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
