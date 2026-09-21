#!/usr/bin/env python3
"""Generate one documentation page per chemistry, from the catalog and its explainers.

    uv run python tools/gen_catalog_pages.py              # with live network summaries
    uv run python tools/gen_catalog_pages.py --fast       # catalog only, no generation
    uv run python tools/gen_catalog_pages.py --only ccm   # rebuild one page, leave the rest

Every page has the same four parts, in reading order:

1. **Introduction** — what the chemistry is, who built it and why.
2. **How it works** — the mechanism, the formal (S, R, A) specification, how
   to use it in Chemart, its parameters, and the implementation decisions.
3. **Results** — what has been shown with it, and what Chemart reproduces.
4. **References** — the book section and the sources.

The prose comes from `catalog/explainers/<id>.md` (see the README there); the
rest is generated from `catalog/chemistries/<id>.yaml` and from a network
actually generated at build time, so parameters, provenance and capabilities
cannot drift from the library. An entry without an explainer falls back to its
YAML `intuition`, `notes` and `phenomena`.

Output: docs/catalog/index.md, which lists the chemistry catalog by family and
then the archive, plus docs/catalog/<id>.md for every entry, archived ones
included (with a banner saying so).
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "docs" / "catalog"
EXPLAINERS = REPO / "catalog" / "explainers"

import chemart                                    # noqa: E402
from chemart.catalog import ARCHIVES, active, load  # noqa: E402

FIDELITY_NOTE = {
    "book": "implemented exactly as the book specifies",
    "book+decisions": "the book left gaps; each filled choice is listed under *Implementation decisions*",
    "reconstructed": "built from the original papers listed under *References*",
}

NETWORK_NOTE = {
    "given": "given — the chemistry is a reaction network, instantiated from its parameters",
    "generated": "generated — the network is the output of the chemistry's algorithm",
}

#: A reader meets these reactor names in the formal specification; each gets a
#: plain-language gloss so the term is never left unexplained.
REACTOR_NOTE = {
    "well-stirred-multiset": "a well-stirred pot of discrete molecules that meet at random, with no notion of position",
    "ode": "deterministic rate equations: concentrations change continuously in time",
    "ssa": "stochastic simulation (Gillespie): individual reaction events happen at random times",
    "lattice-2d": "a two-dimensional grid on which only neighbours interact",
    "continuous-space": "particles that move in continuous space",
    "compartments": "nested compartments (membranes), each with its own contents",
    "graph-rewrite": "graph rewriting: molecules are graphs and reactions rewrite them",
    "maximally-parallel": "every rule that can fire does so at once, in synchronous steps",
    "sequential-vm": "a virtual machine that executes molecules as programs, one instruction at a time",
}

DEFINITION_NOTE = {
    ("S", "explicit"): "listed up front",
    ("S", "implicit"): "defined by a construction rule rather than listed",
    ("R", "explicit"): "listed reaction by reaction",
    ("R", "implicit"): "computed by an algorithm from the molecules that meet",
}

#: The explainer's sections, in order; True marks the required ones.
SECTIONS = {
    "Introduction": True,
    "How it works": True,
    "Using it": True,
    "Results": True,
    "Further reading": False,
}


def angles(text: str) -> str:
    """Escape angle brackets, so notation such as `q<row>=<column>` in the YAML is
    shown rather than swallowed by the browser as unknown HTML tags."""
    return text.replace("<", "&lt;").replace(">", "&gt;")


def esc(text, code: bool = False) -> str:
    """Collapse whitespace and neutralise pipes so text survives a table cell.

    With `code=True` the text goes inside backticks, where Markdown shows it
    verbatim, so angle brackets must not be escaped.
    """
    # `text or ""` would turn a default of 0, 0.0 or False into an empty cell.
    text = " ".join(("" if text is None else str(text)).split()).replace("|", "\\|")
    return text if code else angles(text)


def block(text) -> str:
    """Collapse whitespace but keep it prose, for paragraphs."""
    return angles(" ".join(("" if text is None else str(text)).split()))


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
        default = esc(p.default, code=True)
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



def load_explainer(cid: str, directory: Path = EXPLAINERS) -> dict[str, str] | None:
    """The hand-written prose for one chemistry, split into its `## ` sections.

    Returns None when the chemistry has no explainer yet. A malformed file is an
    error rather than a silently odd page: unknown or repeated sections, text
    before the first section, a `# ` title (the generator writes the title) or a
    missing required section all raise ValueError.
    """
    path = directory / f"{cid}.md"
    if not path.exists():
        return None
    where = path.relative_to(REPO) if path.is_relative_to(REPO) else path
    sections: dict[str, list[str]] = {}
    current = None
    fence = None
    for n, line in enumerate(path.read_text().split("\n"), 1):
        # Headings inside fenced code are code, not structure.
        m = re.match(r"^\s*(`{3,}|~{3,})", line)
        if m:
            if fence is None:
                fence = m.group(1)
            elif m.group(1).startswith(fence):
                fence = None
        if fence is None and not m:
            if line.startswith("# "):
                raise ValueError(f"{where}:{n}: no '# ' title; the generator writes it")
            h = re.match(r"^## (.+?)\s*$", line)
            if h:
                current = h.group(1)
                if current not in SECTIONS:
                    raise ValueError(f"{where}:{n}: unknown section {current!r}; "
                                     f"expected one of {list(SECTIONS)}")
                if current in sections:
                    raise ValueError(f"{where}:{n}: section {current!r} appears twice")
                sections[current] = []
                continue
        if current is None:
            if line.strip():
                raise ValueError(f"{where}:{n}: text before the first '## ' section")
            continue
        sections[current].append(line)

    out = {name: "\n".join(lines).strip() for name, lines in sections.items()}
    missing = [s for s, required in SECTIONS.items() if required and not out.get(s)]
    if missing:
        raise ValueError(f"{where}: missing section(s) {missing}")
    order = [s for s in SECTIONS if s in out]
    if list(out) != order:
        raise ValueError(f"{where}: sections must come in the order {order}")
    return out


def glance(c) -> list[str]:
    """The classification table that closes the introduction."""
    L = ["| at a glance | |", "|---|---|"]
    L.append(f"| **family** | {c.family} |")
    L.append(f"| **kind** | {c.kind} |")
    if c.network:
        L.append(f"| **network** | {NETWORK_NOTE.get(c.network, c.network)} |")
    L.append(f"| **constructive** | {'yes — the species set grows at run time' if c.constructive else 'no — fixed species set'} |")
    L.append(f"| **fidelity** | `{c.fidelity}` — {FIDELITY_NOTE.get(c.fidelity, '')} |")
    L.append(f"| **book** | {esc(c.book)} |")
    L.append(f"| **provides** | {', '.join(f'`{p}`' for p in c.provides) or '—'} |")
    return L + [""]


def paragraphs(text) -> list[str]:
    """YAML prose as Markdown paragraphs: whitespace collapsed, blank lines kept."""
    return [block(p) for p in re.split(r"\n\s*\n", str(text or "")) if p.strip()]


def specification(c) -> list[str]:
    """The (S, R, A) triple, with the notes that say what the one-liners leave out."""
    L = ["### Formal specification", ""]
    L.append("Banzhaf and Yamamoto describe every artificial chemistry by three things "
             "(book §2.3): the set of possible molecules **S**, the reaction rules **R** "
             "that transform them, and the reactor algorithm **A** that decides which "
             "reactions happen, and when. For this chemistry:")
    L.append("")

    def definition(part):
        d = getattr(c, part).get("definition", "?")
        note = DEFINITION_NOTE.get((part, d))
        return f"*{d}*" + (f" — {note}" if note else "")

    L.append(f"**S — molecules** ({definition('S')}). {block(c.S.get('repr'))}")
    L.append("")
    L += [p + "\n" for p in paragraphs(c.S.get("notes"))]

    arity = c.R.get("arity", "?")
    if isinstance(arity, list):
        arity = ", ".join(map(str, arity))
    L.append(f"**R — reactions** ({definition('R')}; molecules per reaction: {arity}). "
             f"{block(c.R.get('scheme'))}")
    L.append("")
    L += [p + "\n" for p in paragraphs(c.R.get("notes"))]

    reactors = c.A.get("reactor", "?")
    reactors = reactors if isinstance(reactors, list) else [reactors]
    glossed = [f"`{r}` ({REACTOR_NOTE[r]})" if r in REACTOR_NOTE else f"`{r}`" for r in reactors]
    L.append(f"**A — reactor**: {'; or '.join(glossed)}.")
    L.append("")
    if c.A.get("dilution"):
        L.append(f"*How the population is bounded:* {block(c.A.get('dilution'))}")
        L.append("")
    L += [p + "\n" for p in paragraphs(c.A.get("notes"))]
    return L


def default_network(c, net, error) -> list[str]:
    L = []
    if error:
        L.append("!!! failure \"This page was built without generating a network\"")
        L.append(f"    `{error}`")
        L.append("")
    elif net is not None:
        L.append("Every parameter at its default:")
        L.append("")
        L.append("```python")
        L.append("import chemart")
        L.append(f'net = chemart.generate_network("{c.id}", seed=1)')
        L.append("print(net.summary())")
        L.append("```")
        L.append("")
        L.append("```")
        L.append(net.summary())
        L.append("```")
        L.append("")
        shown = net.reactions[:8]
        if shown:
            L.append("Its first reactions (`net.reactions`):")
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
    return L


def parameters(c) -> list[str]:
    L = ["#### Parameters", ""]
    if not c.params:
        return L + ["This chemistry takes no parameters.", ""]
    L.append("Pass any of these as keyword arguments to `generate_network`. The *role* "
             "column says what a parameter controls: `structural` (which molecules and "
             "reactions exist), `kinetic` (rates), `thermodynamic` (energies, "
             "temperature), `population` (sizes, budgets, initial state), `spatial`, "
             "`stochastic` or `selection`. *range* gives the values used in the "
             "published work.")
    L.append("")
    L.append("| name | type | default | role | what it does |")
    L.append("|---|---|---|---|---|")
    L += param_rows(c)
    return L + [""]


def decisions(c) -> list[str]:
    if not c.decisions:
        return []
    L = ["### Implementation decisions", ""]
    L.append("The sources leave gaps, and sometimes contradict each other or the book. "
             "Each such case, and how Chemart resolved it, is listed here: read these "
             "before quoting a number from this page.")
    L.append("")
    n = len(c.decisions)
    L.append(f"??? note \"{n} decision{'s' if n != 1 else ''}\"")
    L.append("")
    for d in c.decisions:
        L.append(f"    - {block(d)}")
    return L + [""]


def link(ref: str) -> str:
    """A reference as a clickable link where it is one."""
    ref = ref.strip()
    if ref.startswith("doi:"):
        return f"[{ref}](https://doi.org/{ref[4:].strip()})"
    if ref.startswith(("http://", "https://")):
        return f"<{ref}>"
    return esc(ref)


def references(c, ex) -> list[str]:
    L = ["## References", ""]
    numbers = [r for r in c.refs if re.fullmatch(r"\[\d+\]", r.strip())]
    others = [r for r in c.refs if r not in numbers]
    book = (f"Banzhaf, W. & Yamamoto, L. (2015). *Artificial Chemistries*. MIT Press. "
            f"Section {esc(c.book)}")
    if numbers:
        book += f"; the book's bibliography entries {', '.join(numbers)}"
    L.append(book + ".")
    L.append("")
    if others:
        L += [f"- {link(r)}" for r in others] + [""]
    if c.sources:
        L.append("### Sources used for the implementation")
        L.append("")
        L += [f"- {block(s)}" for s in c.sources] + [""]
    if ex and ex.get("Further reading"):
        L.append("### Further reading")
        L.append("")
        L.append(ex["Further reading"])
        L.append("")
    return L


def test_files(cid: str) -> list[str]:
    """The chemistry's own test file, or else the shared ones that test it by id."""
    tests = REPO / "tests" / "chemistries"
    own = tests / f"test_{cid.replace('-', '_')}.py"
    if own.exists():
        return [str(own.relative_to(REPO))]
    return [str(p.relative_to(REPO)) for p in sorted(tests.glob("test_*.py"))
            if f'"{cid}"' in p.read_text()]


def page(c, net, error, ex=None) -> str:
    """One chemistry's page: Introduction > How it works > Results > References."""
    source = f"catalog/chemistries/{c.id}.yaml"
    if ex:
        source += f" and catalog/explainers/{c.id}.md"
    L: list[str] = [f"<!-- Generated by tools/gen_catalog_pages.py from {source}. "
                    "Edit those, not this file. -->", ""]
    L.append(f"# {c.name}")
    L.append("")
    L.append(f"`{c.id}`{' · *' + esc(c.origin) + '*' if c.origin else ''}")
    L.append("")
    if c.aliases:
        L.append("*Also known as:* " + ", ".join(f"*{esc(a)}*" for a in c.aliases))
        L.append("")

    if c.archived:
        L.append(f'!!! warning "Archived: {c.archived}"')
        L.append(f"    This entry is not part of the chemistry catalog "
                 f"({ARCHIVES.get(c.archived, 'archived')}). It keeps its specification, "
                 "generator and tests, and `generate_network` still runs it by id, but "
                 "listings and the LLM tools leave it out. See the "
                 "[catalog's archive](index.md#archive).")
        L.append("")

    # 1. What the chemistry is, before any table or formalism.
    L.append("## Introduction")
    L.append("")
    if ex:
        L.append(ex["Introduction"])
    else:
        L.append(block(c.intuition))
    L.append("")
    L += glance(c)

    # 2. The mechanism, then the precise model, then how to drive it.
    L.append("## How it works")
    L.append("")
    if ex:
        L.append(ex["How it works"])
        L.append("")
    elif c.notes:
        L += [p + "\n" for p in paragraphs(c.notes)]
    L += specification(c)
    L.append("### Using it in Chemart")
    L.append("")
    L += default_network(c, net, error)
    if ex:
        L.append(ex["Using it"])
        L.append("")
    L += parameters(c)
    L += decisions(c)

    # 3. What has been shown with it.
    L.append("## Results")
    L.append("")
    if ex:
        L.append(ex["Results"])
        L.append("")
    elif c.phenomena:
        L.append("What the literature reports this model produces:")
        L.append("")
        L += [f"- {block(p)}" for p in c.phenomena] + [""]

    # 4. Where it all comes from.
    L += references(c, ex)
    L.append("---")
    L.append("")
    L.append(f"*Specification: `catalog/chemistries/{c.id}.yaml` · "
             f"explainer: `catalog/explainers/{c.id}.md` · "
             f"generator: `chemart/chemistries/{c.id.replace('-', '_')}.py` · "
             f"tests: {', '.join(f'`{t}`' for t in test_files(c.id)) or 'none'}*")
    L.append("")
    return "\n".join(L)


#: Capabilities worth putting in a table. `topology` and `stoichiometry` are on
#: every entry and `initial-state` on most, so listing them in every row costs
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


def _index_rows(group, nets) -> list[str]:
    L = ["| chemistry | origin | network | grows | fidelity | S/R | beyond topology |",
         "|---|---|---|:--:|---|--:|---|"]
    for c in sorted(group, key=lambda c: c.id):
        net = nets.get(c.id)
        size = f"{len(net.species)}/{len(net.reactions)}" if net else "—"
        L.append(
            f"| [{c.name}]({c.id}.md) | {short_origin(c.origin)} | {c.network or '—'} | "
            f"{'yes' if c.constructive else '·'} | {c.fidelity} | {size} | "
            f"{table_caps(c)} |"
        )
    return L + [""]


def glance_counts(entries) -> list[str]:
    """What the catalog holds, counted here so no document states it by hand."""
    def having(*tags):
        return sum(1 for c in entries if set(tags) & set(c.provides))

    kinds = Counter(c.kind for c in entries)
    L = ["| | |", "|---|---|"]
    L.append(f"| network given (written down, built by a formula, or drawn at random) | {sum(1 for c in entries if c.network == 'given')} |")
    L.append(f"| network generated (the output of the chemistry's algorithm) | {sum(1 for c in entries if c.network == 'generated')} |")
    L.append(f"| constructive (open, growing species set) | {sum(1 for c in entries if c.constructive)} |")
    L.append(f"| carry their own rate constants or rate law | {having('rate-constants', 'rate-law')} |")
    L.append(f"| carry energetics or thermodynamic consistency | {having('energies', 'thermodynamic-consistency')} |")
    L.append(f"| declare a conservation law | {having('mass-conservation')} |")
    L.append(f"| define space | {having('space')} |")
    L.append(f"| define compartments | {having('compartments')} |")
    L.append("")
    L.append("By kind: " + ", ".join(f"{n} {k}" for k, n in kinds.most_common()) + ".")
    return L + [""]


def index_page(entries, nets) -> str:
    main = active(entries)
    archived = [c for c in entries if c.archived is not None]
    L: list[str] = []
    L.append("# Catalog")
    L.append("")
    L.append(f"**{len(main)}** chemistries, one page each, generated from the "
             "catalog so these pages cannot drift from the library."
             + (f" A further {len(archived)} archived entries are listed "
                "[at the end](#archive)." if archived else ""))
    L.append("")
    L += glance_counts(main)
    L.append("Columns: **network** is *given* when the chemistry is a reaction network "
             "Chemart instantiates from its parameters, and *generated* when the network "
             "is the output of the chemistry's algorithm (see `catalog/NETWORKS.md`); "
             "**grows** is whether the species set is open and expands at "
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
    for c in main:
        by_family.setdefault(c.family, []).append(c)
    for family in sorted(by_family):
        L.append(f"## {family}")
        L.append("")
        L += _index_rows(by_family[family], nets)

    if archived:
        L.append("## Archive")
        L.append("")
        L.append("These entries are not part of the chemistry catalog. `list_chemistries()` "
                 "and the LLM tools leave them out, but each keeps its specification, "
                 "generator, tests and page, and `generate_network(id)` still runs it.")
        L.append("")
        for group in sorted({c.archived for c in archived}):
            members = [c for c in archived if c.archived == group]
            L.append(f"### {group}")
            L.append("")
            L.append(f"{len(members)} entries: {ARCHIVES.get(group, '')}.")
            L.append("")
            L += _index_rows(members, nets)

    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fast", action="store_true",
                    help="skip generating networks (much faster, less useful pages)")
    ap.add_argument("--only", nargs="+", metavar="ID",
                    help="rebuild just these pages; the index and every other page are left alone")
    args = ap.parse_args(argv)

    entries = sorted(load(), key=lambda c: c.id)
    if args.only:
        known = {c.id for c in entries}
        unknown = [i for i in args.only if i not in known]
        if unknown:
            ap.error(f"unknown chemistry id(s): {', '.join(unknown)}")
        entries = [c for c in entries if c.id in args.only]

    # Read every explainer before touching the output, so a malformed one stops
    # the build instead of leaving half the pages regenerated.
    explainers = {c.id: load_explainer(c.id) for c in entries}

    OUT.mkdir(parents=True, exist_ok=True)
    if not args.only:
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
        (OUT / f"{c.id}.md").write_text(page(c, nets.get(c.id), errors.get(c.id), explainers[c.id]))
    if not args.only:
        (OUT / "index.md").write_text(index_page(entries, nets))

    written = len(entries) + (0 if args.only else 1)
    with_prose = sum(1 for e in explainers.values() if e)
    print(f"wrote {written} page(s) to {OUT.relative_to(REPO)} "
          f"({with_prose}/{len(entries)} with an explainer)")
    if errors:
        print(f"{len(errors)} chemistries failed to generate:", file=sys.stderr)
        for cid, err in errors.items():
            print(f"  {cid}: {err}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
