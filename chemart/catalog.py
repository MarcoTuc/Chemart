"""Loader, validator and index generator for the Chemart chemistry catalog.

The catalog is the ground truth for what Chemart can generate and the only
specification of each generator's parameters. It is kept as YAML (one file
per chemistry, `catalog/chemistries/<id>.yaml`, see catalog/SCHEMA.md) so
that it stays reviewable and diffable;
this module turns it into typed Python objects and enforces the invariants
that the rest of the library relies on.

An entry is on schema v2 exactly when its generator module
`chemart/chemistries/<id with underscores>.py` exists; v2 entries must have
JSON-typed parameters with defaults and provenance fields (see SCHEMA.md).

    python -m chemart.catalog validate
    python -m chemart.catalog index            # regenerates docs/CATALOG.md
    python -m chemart.catalog status           # implementation progress
    python -m chemart.catalog show matrix-chemistry
    python -m chemart.catalog query --provides rate-constants --ready yes
"""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import MISSING, dataclass, field, fields
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import yaml

ROOT = Path(__file__).resolve().parent.parent
#: A built wheel carries the catalog inside the package (hatch force-include);
#: a source checkout reads it from catalog/chemistries/.
_PACKAGED_CATALOG = Path(__file__).resolve().parent / "_catalog"
CATALOG_DIR = _PACKAGED_CATALOG if _PACKAGED_CATALOG.is_dir() else ROOT / "catalog" / "chemistries"
MODULES_DIR = Path(__file__).resolve().parent / "chemistries"

FAMILIES = {
    "core", "rewriting", "automata", "bio-inspired", "origin-of-life",
    "evolutionary-dynamics", "network", "spatial", "application",
    "systems-biology", "wet", "non-chemical",
}
KINDS = {"generator", "formalism", "framework", "analysis", "wet"}
ROLES = {
    "structural", "kinetic", "thermodynamic", "population", "spatial",
    "stochastic", "selection",
}
PROVIDES = {
    "topology", "stoichiometry", "catalysts", "rate-constants", "rate-law",
    "energies", "thermodynamic-consistency", "mass-conservation", "flow",
    "space", "compartments", "initial-state", "sequence-structure-function",
}
READINESS = {"yes", "partial", "no"}
FIDELITY = {"book", "book+decisions", "reconstructed"}
#: Hub entries may also be new chemistries that come from no publication.
HUB_FIDELITY = FIDELITY | {"original"}

#: Keyword arguments of generate_network; a parameter may not shadow them.
RESERVED_PARAMS = {"seed", "revision", "trust_remote_code"}

#: v2 parameter types: JSON values only.
PARAM_TYPES = {"int", "float", "bool", "str", "enum", "list", "dict"}
_JSON_TYPE = {
    "int": "integer", "float": "number", "bool": "boolean", "str": "string",
    "list": "array", "dict": "object",
}

#: Capability tiers, in increasing order of what a generator must supply.
TIERS = {
    "topology": {"topology"},
    "kinetics": {"rate-constants", "rate-law"},
    "thermodynamics": {"energies", "thermodynamic-consistency"},
}


@dataclass
class Param:
    name: str
    role: str
    type: str | None = None
    default: Any = None
    range: str | None = None
    meaning: str | None = None
    min: float | None = None
    max: float | None = None
    choices: list | None = None

    def coerce(self, value: Any) -> Any:
        """Return `value` checked (and int->float widened); raise ValueError otherwise."""
        t = self.type
        is_number = isinstance(value, (int, float)) and not isinstance(value, bool)
        if t == "int":
            ok = isinstance(value, int) and not isinstance(value, bool)
        elif t == "float":
            ok = is_number
            if ok:
                value = float(value)
        elif t == "bool":
            ok = isinstance(value, bool)
        elif t == "str":
            ok = isinstance(value, str)
        elif t == "enum":
            ok = value in (self.choices or [])
        elif t == "list":
            if isinstance(value, tuple):
                value = list(value)
            ok = isinstance(value, list)
        elif t == "dict":
            ok = isinstance(value, dict)
        else:
            raise ValueError(f"parameter {self.name!r} has non-JSON catalog type {t!r}")
        if not ok:
            raise ValueError(f"{self.name}={value!r} is invalid: expected {self.expected()}")
        if self.min is not None and value < self.min:
            raise ValueError(f"{self.name}={value!r} is invalid: must be >= {self.min}")
        if self.max is not None and value > self.max:
            raise ValueError(f"{self.name}={value!r} is invalid: must be <= {self.max}")
        return value

    def expected(self) -> str:
        if self.type == "enum":
            return f"one of {self.choices}"
        return {"int": "an integer", "float": "a number", "bool": "true or false",
                "str": "a string", "list": "a list", "dict": "an object"}.get(self.type, str(self.type))

    def json_schema(self) -> dict[str, Any]:
        schema: dict[str, Any] = {}
        if self.type == "enum":
            schema["enum"] = self.choices
        elif self.type in _JSON_TYPE:
            schema["type"] = _JSON_TYPE[self.type]
        if self.default is not None:
            schema["default"] = self.default
        if self.min is not None:
            schema["minimum"] = self.min
        if self.max is not None:
            schema["maximum"] = self.max
        description = " ".join(str(self.meaning or "").split())
        if self.range:
            description += f" (range: {self.range})"
        if self.type not in PARAM_TYPES:
            description += f" [catalog type {self.type!r}, not yet normalised]"
        schema["description"] = description.strip()
        return schema


@dataclass
class Chemistry:
    id: str
    name: str
    family: str
    kind: str
    constructive: bool
    S: dict
    R: dict
    A: dict
    #: Section(s) of Banzhaf & Yamamoto; required in the built-in catalog,
    #: optional for chemistries shared on the hub.
    book: str | None = None
    params: list[Param] = field(default_factory=list)
    provides: list[str] = field(default_factory=list)
    generator_ready: str | None = None
    fidelity: str | None = None
    sources: list[str] = field(default_factory=list)
    decisions: list[str] = field(default_factory=list)
    aliases: list[str] = field(default_factory=list)
    origin: str | None = None
    refs: list[str] = field(default_factory=list)
    reference_impl: str | None = None
    phenomena: list[str] = field(default_factory=list)
    #: Plain-language explanation of how the chemistry works and why, for
    #: readers meeting it for the first time. The (S, R, A) fields say what it
    #: is; this says what the idea is.
    intuition: str | None = None
    notes: str | None = None
    source_file: str = ""

    @property
    def module(self) -> str:
        return f"chemart.chemistries.{self.id.replace('-', '_')}"

    @property
    def implemented(self) -> bool:
        return (MODULES_DIR / f"{self.id.replace('-', '_')}.py").exists()

    @property
    def tiers(self) -> list[str]:
        """Which capability tiers this chemistry can populate."""
        have = set(self.provides)
        return [t for t, need in TIERS.items() if have & need]

    def params_by_role(self, role: str) -> list[Param]:
        return [p for p in self.params if p.role == role]


_PARAM_FIELDS = {f.name for f in fields(Param)}
_ENTRY_FIELDS = {f.name for f in fields(Chemistry)} - {"source_file"}


def _required(cls) -> list[str]:
    return [f.name for f in fields(cls) if f.default is MISSING and f.default_factory is MISSING]


def _check_keys(raw: Any, cls, allowed: set[str], what: str) -> None:
    if not isinstance(raw, dict):
        raise ValueError(f"{what} must be a mapping, got {type(raw).__name__}")
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"{what}: unknown field(s) {', '.join(map(str, unknown))}; "
                         f"valid fields: {', '.join(sorted(allowed))}")
    missing = [name for name in _required(cls) if name not in raw]
    if missing:
        raise ValueError(f"{what}: missing required field(s) {', '.join(missing)}")


def parse_entry(raw: Any, source_file: str = "") -> Chemistry:
    """Build one `Chemistry` from its YAML mapping; raise ValueError if malformed."""
    _check_keys(raw, Chemistry, _ENTRY_FIELDS, "catalog entry")
    raw = dict(raw)
    # YAML 1.1 turns bare yes/no into booleans; generator_ready is a
    # three-valued string, so put it back.
    ready = raw.get("generator_ready")
    if isinstance(ready, bool):
        raw["generator_ready"] = "yes" if ready else "no"
    params = raw.get("params") or []
    if not isinstance(params, list):
        raise ValueError(f"{raw['id']}: params must be a list")
    parsed = []
    for i, p in enumerate(params):
        _check_keys(p, Param, _PARAM_FIELDS, f"{raw['id']}: params[{i}]")
        parsed.append(Param(**p))
    raw["params"] = parsed
    raw["source_file"] = source_file
    return Chemistry(**raw)


@lru_cache(maxsize=None)
def _parse(catalog_dir: Path) -> tuple[Chemistry, ...]:
    entries: list[Chemistry] = []
    for path in sorted(catalog_dir.glob("*.yaml")):
        doc = yaml.safe_load(path.read_text())
        for raw in doc["chemistries"]:
            entries.append(parse_entry(raw, path.name))
    if not entries:
        raise RuntimeError(f"no catalog entries found in {catalog_dir}")
    return tuple(entries)


def load(catalog_dir: Path = CATALOG_DIR) -> list[Chemistry]:
    """Return the catalog entries, parsing the YAML at most once per directory.

    Callers each get their own list, so nothing is shared but the entries
    themselves, which no caller mutates. `Chemistry.implemented` stats the
    filesystem on every access, so a cached entry still reports fresh
    implementation status as generator modules appear.

    A process that writes a catalog file and re-reads it must call
    `load.cache_clear()` first; the CLI is a fresh process each time.
    """
    return list(_parse(catalog_dir))


load.cache_clear = _parse.cache_clear  # type: ignore[attr-defined]


def validate(entries: Iterable[Chemistry]) -> list[str]:
    """Return a list of problems; empty means the catalog is consistent."""
    problems: list[str] = []
    seen: dict[str, str] = {}

    for c in entries:
        where = f"{c.source_file}:{c.id}"
        if c.id in seen:
            problems.append(f"{where}: duplicate id, also in {seen[c.id]}")
        seen[c.id] = c.source_file
        problems += [f"{where}: {m}" for m in entry_problems(c)]

    for path in sorted(MODULES_DIR.glob("*.py")):
        if path.stem != "__init__" and path.stem.replace("_", "-") not in seen:
            problems.append(f"chemart/chemistries/{path.name}: no catalog entry with this id")
    return problems


def entry_problems(c: Chemistry, *, hub: bool = False) -> list[str]:
    """Problems with one entry on its own. `hub=True` applies the rules for
    chemistries shared on the hub: every entry there has a generator (so the
    v2 rules always apply), `book` is optional and fidelity may be `original`."""
    problems: list[str] = []
    if c.family not in FAMILIES:
        problems.append(f"unknown family {c.family!r}")
    if c.kind not in KINDS:
        problems.append(f"unknown kind {c.kind!r}")
    if not hub and not c.book:
        problems.append("missing book section")
    for p in c.provides:
        if p not in PROVIDES:
            problems.append(f"unknown capability {p!r}")
    for name, n in Counter(p.name for p in c.params).items():
        if n > 1:
            problems.append(f"duplicate param {name!r}")
    for p in c.params:
        if p.role not in ROLES:
            problems.append(f"param {p.name!r} has unknown role {p.role!r}")
        if p.name in RESERVED_PARAMS:
            problems.append(f"param {p.name!r} is reserved: it is an argument of generate_network")

    # Semantic invariants.
    if "stoichiometry" in c.provides and "topology" not in c.provides:
        problems.append("stoichiometry implies topology")
    if "thermodynamic-consistency" in c.provides and "rate-constants" not in c.provides:
        problems.append(
            "thermodynamic consistency constrains reverse rates, "
            "so rate-constants must also be provided"
        )
    for section in ("S", "R", "A"):
        if not getattr(c, section):
            problems.append(f"missing {section} section")

    if hub or c.implemented:
        problems += _v2_problems(c, hub=hub)
    else:
        if c.generator_ready not in READINESS:
            problems.append(f"generator_ready must be one of {sorted(READINESS)}")
        if c.fidelity is not None:
            problems.append(f"fidelity is set but {c.module} does not exist")
        if c.kind == "generator" and c.generator_ready == "yes" and not c.provides:
            problems.append("a ready generator must provide something")
    return problems


def _v2_problems(c: Chemistry, *, hub: bool = False) -> list[str]:
    """Schema v2 invariants, enforced for entries that have a generator module."""
    out: list[str] = []
    if c.generator_ready is not None:
        out.append("implemented entries replace generator_ready with fidelity")
    allowed = HUB_FIDELITY if hub else FIDELITY
    if c.fidelity not in allowed:
        out.append(f"fidelity must be one of {sorted(allowed)}")
    if c.fidelity == "reconstructed" and not c.sources:
        out.append("reconstructed entries must list their sources")
    if c.fidelity == "book+decisions" and not c.decisions:
        out.append("book+decisions entries must list their decisions")
    if not c.provides:
        out.append("implemented entries must declare provides")
    if not (c.intuition or "").strip():
        out.append("missing intuition: implemented entries need a plain-language "
                   "explanation of how the chemistry works and why")
    for p in c.params:
        if p.type not in PARAM_TYPES:
            hint = " (seed is an argument of generate_network, not a param)" if p.type == "seed" else ""
            out.append(f"param {p.name!r}: type {p.type!r} not in {sorted(PARAM_TYPES)}{hint}")
            continue
        if p.type == "enum" and not p.choices:
            out.append(f"param {p.name!r}: enum needs choices")
        if (p.min is not None or p.max is not None) and p.type not in ("int", "float"):
            out.append(f"param {p.name!r}: min/max only apply to int and float")
        if not p.meaning:
            out.append(f"param {p.name!r}: missing meaning")
        if p.default is None:
            out.append(f"param {p.name!r}: missing default")
        else:
            try:
                p.coerce(p.default)
            except ValueError as err:
                out.append(f"param {p.name!r}: default violates its own spec: {err}")
    return out


# --------------------------------------------------------------------------
# Index generation
# --------------------------------------------------------------------------

_TIER_MARK = {"topology": "T", "kinetics": "K", "thermodynamics": "H"}


def _status(c: Chemistry) -> str:
    return c.fidelity if c.implemented else f"todo ({c.generator_ready})"


def render_index(entries: list[Chemistry]) -> str:
    by_family: dict[str, list[Chemistry]] = {}
    for c in entries:
        by_family.setdefault(c.family, []).append(c)

    out: list[str] = []
    out.append("<!-- generated by `python -m chemart.catalog index` - do not edit -->")
    out.append("# Chemart catalog index\n")
    out.append(
        f"{len(entries)} chemistries collected from Banzhaf & Yamamoto, "
        "*Artificial Chemistries* (MIT Press, 2015).\n"
    )
    out.append(
        "Tier column: **T** topology, **K** kinetics (the chemistry itself "
        "prescribes rates or a rate law), **H** thermodynamics (energies "
        "and/or detailed balance). Status column: the fidelity of the "
        "implemented generator, or `todo (<book readiness>)`.\n"
    )

    done = [c for c in entries if c.implemented]
    fidelity = Counter(c.fidelity for c in done)
    kin = sum(1 for c in entries if "kinetics" in c.tiers)
    thermo = sum(1 for c in entries if "thermodynamics" in c.tiers)
    constructive = sum(1 for c in entries if c.constructive)
    out.append(
        f"- implemented: **{len(done)}** of {len(entries)}"
        + (" (" + ", ".join(f"{k}: {v}" for k, v in sorted(fidelity.items())) + ")" if done else "")
        + "\n"
        f"- constructive (open, growing species set): **{constructive}**\n"
        f"- carry their own kinetics: **{kin}**; carry energetics: **{thermo}**\n"
    )

    for family in sorted(by_family):
        out.append(f"\n## {family}\n")
        out.append("| id | name | kind | constructive | tier | status | book |")
        out.append("|---|---|---|---|---|---|---|")
        for c in sorted(by_family[family], key=lambda c: c.id):
            tier = "".join(_TIER_MARK[t] for t in ("topology", "kinetics", "thermodynamics") if t in c.tiers) or "-"
            out.append(
                f"| `{c.id}` | {c.name} | {c.kind} | "
                f"{'yes' if c.constructive else 'no'} | {tier} | "
                f"{_status(c)} | {(c.book or '-').split(';')[0]} |"
            )
    return "\n".join(out) + "\n"


def render_status(entries: list[Chemistry]) -> str:
    done = [c for c in entries if c.implemented]
    lines = [f"implemented {len(done)}/{len(entries)}"]
    for fid, n in sorted(Counter(c.fidelity for c in done).items()):
        lines.append(f"  {fid}: {n}")
    todo: dict[str, list[str]] = {}
    for c in entries:
        if not c.implemented:
            todo.setdefault(c.family, []).append(c.id)
    if todo:
        lines.append("remaining:")
        for family in sorted(todo):
            lines.append(f"  {family}: {', '.join(sorted(todo[family]))}")
    return "\n".join(lines)


def _show(c: Chemistry) -> str:
    lines = [f"# {c.name}  (`{c.id}`)", ""]
    if c.aliases:
        lines.append(f"aliases: {', '.join(c.aliases)}")
    if c.origin:
        lines.append(f"origin: {c.origin}")
    if c.book:
        lines.append(f"book: {c.book}")
    lines += [
        f"family/kind: {c.family} / {c.kind}"
        f"{'  (constructive)' if c.constructive else ''}",
        f"provides: {', '.join(c.provides) or '-'}",
        f"status: {_status(c)}",
        "",
        f"S  {c.S.get('definition', '?')}: {c.S.get('repr', '')}",
        f"R  {c.R.get('definition', '?')}, arity {c.R.get('arity', '?')}: "
        f"{c.R.get('scheme', '')}",
        f"A  reactor {c.A.get('reactor', '?')}; dilution: {c.A.get('dilution', '?')}",
        "",
        "parameters:",
    ]
    for role in sorted({p.role for p in c.params}):
        lines.append(f"  [{role}]")
        for p in c.params_by_role(role):
            default = "" if p.default is None else f" = {p.default}"
            bounds = []
            if p.choices:
                bounds.append(f"choices {p.choices}")
            if p.min is not None or p.max is not None:
                bounds.append(f"[{p.min}, {p.max}]")
            if p.range:
                bounds.append(p.range)
            extra = f"  ({'; '.join(bounds)})" if bounds else ""
            lines.append(f"    {p.name}: {p.type}{default}{extra}")
            if p.meaning:
                lines.append(f"        {p.meaning}")
    for title, items in (("sources", c.sources), ("decisions", c.decisions), ("phenomena", c.phenomena)):
        if items:
            lines += ["", f"{title}:"] + [f"  - {i}" for i in items]
    if c.notes:
        lines += ["", "notes:", "  " + c.notes.strip().replace("\n", "\n  ")]
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    cmd = argv[1] if len(argv) > 1 else "validate"
    entries = load()

    if cmd == "validate":
        problems = validate(entries)
        if "--only" in argv:
            wanted = argv[argv.index("--only") + 1]
            module = f"chemart/chemistries/{wanted.replace('-', '_')}.py"
            problems = [p for p in problems if p.startswith(f"{wanted}.yaml:{wanted}:") or p.startswith(module)]
        for p in problems:
            print(p, file=sys.stderr)
        print(f"{len(entries)} chemistries, {len(problems)} problems")
        return 1 if problems else 0

    if cmd == "index":
        target = ROOT / "docs" / "CATALOG.md"
        target.write_text(render_index(entries))
        print(f"wrote {target.relative_to(ROOT)}")
        return 0

    if cmd == "status":
        print(render_status(entries))
        return 0

    if cmd == "show":
        wanted = argv[2]
        for c in entries:
            if c.id == wanted:
                print(_show(c))
                return 0
        print(f"no such chemistry: {wanted}", file=sys.stderr)
        return 1

    if cmd == "query":
        args = argv[2:]
        provides = [args[i + 1] for i, a in enumerate(args) if a == "--provides"]
        ready = next((args[i + 1] for i, a in enumerate(args) if a == "--ready"), None)
        family = next((args[i + 1] for i, a in enumerate(args) if a == "--family"), None)
        for c in entries:
            if provides and not set(provides) <= set(c.provides):
                continue
            if ready and c.generator_ready != ready:
                continue
            if family and c.family != family:
                continue
            print(f"{c.id:34s} {c.name}")
        return 0

    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
