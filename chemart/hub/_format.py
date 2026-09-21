"""What a hub repo may contain, checked without running any of it.

The same function, `inspect`, runs on the uploader's machine before a push
and on the server when a commit arrives, so there is one set of rules. It
only parses: YAML with anchors and aliases refused, JSON with NaN refused,
Python through `ast.parse`. Uploaded code is never imported here.

A **generator repo** holds a chemistry:

    chemart.yaml    `hub:` block + `chemistries: [one catalog entry]`   (required)
    generator.py    defines generate(p, rng) -> Network                 (required)
    *.py            helpers, imported relatively (from .rules import x)
    preview.json    Network.to_dict() at default parameters, made by the uploader
    README.md, LICENSE

A **network repo** holds one reaction network:

    network.json    Network.to_dict()                                   (required)
    chemart.yaml    `hub:` block only                                   (optional)
    README.md, LICENSE, other .json/.txt data
"""

from __future__ import annotations

import ast
import json
import re
from dataclasses import dataclass, field
from typing import Any, Collection, Mapping

import yaml

from chemart import catalog
from chemart.hub import _config
from chemart.hub._ids import RepoId, check_path
from chemart.network import Network

FORMAT = 1
REPO_TYPES = ("generator", "network")
CHEMART_YAML = "chemart.yaml"

_HUB_KEYS = {
    "format", "repo_type", "title", "description", "license", "tags",
    "requires_chemart", "requires", "builtin",
}
_TAG_RE = re.compile(r"^[a-z0-9]+(?:[-.][a-z0-9]+)*$")
_MODULE_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*$")
_SPEC_RE = re.compile(r"^\s*(?:>=|<=|==|!=|>|<)\s*\d+(?:\.\d+)*\s*(?:,\s*(?:>=|<=|==|!=|>|<)\s*\d+(?:\.\d+)*\s*)*$")


class FormatError(ValueError):
    """A repo's files break the format; `problems` lists every reason."""

    def __init__(self, problems: list[str]):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))


@dataclass
class RepoCard:
    """Everything the hub needs to know about one snapshot of a repo."""

    repo_type: str
    hub: dict[str, Any]
    entry: catalog.Chemistry | None = None
    network: Network | None = None          # network.json, or a generator's preview.json
    readme: str | None = None
    imports: list[str] = field(default_factory=list)

    @property
    def has_code(self) -> bool:
        return self.repo_type == "generator" and not self.builtin

    @property
    def builtin(self) -> str | None:
        return self.hub.get("builtin")

    @property
    def title(self) -> str:
        if self.hub.get("title"):
            return self.hub["title"]
        if self.entry is not None:
            return self.entry.name
        return (self.network.chemistry if self.network else "") or ""

    @property
    def summary(self) -> str:
        if self.hub.get("description"):
            return " ".join(self.hub["description"].split())
        if self.entry is not None:
            text = " ".join((self.entry.intuition or self.entry.S.get("repr") or "").split())
            return text if len(text) <= 280 else text[:279] + "…"
        if self.network is not None:
            return self.network.summary().splitlines()[0]
        return ""

    @property
    def provides(self) -> list[str]:
        if self.entry is not None:
            return list(self.entry.provides)
        return self.network.provides if self.network is not None else []


# --------------------------------------------------------------------------
# Safe parsing
# --------------------------------------------------------------------------

def _text(data: bytes, what: str) -> str:
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        raise FormatError([f"{what} is not UTF-8 text"]) from None


def load_yaml(data: bytes, what: str = CHEMART_YAML) -> Any:
    """`yaml.safe_load`, refusing anchors and aliases (the "billion laughs" bomb)."""
    text = _text(data, what)
    try:
        for event in yaml.parse(text, Loader=yaml.SafeLoader):
            if isinstance(event, yaml.AliasEvent) or getattr(event, "anchor", None):
                raise FormatError([f"{what}: YAML anchors and aliases are not allowed"])
        return yaml.safe_load(text)
    except yaml.YAMLError as err:
        raise FormatError([f"{what} is not valid YAML: {err}"]) from None


def _no_constants(name: str):
    raise ValueError(f"{name} is not valid JSON")


def load_json(data: bytes, what: str) -> Any:
    """`json.loads`, refusing NaN and Infinity (they are not JSON)."""
    text = _text(data, what)
    try:
        return json.loads(text, parse_constant=_no_constants)
    except (ValueError, RecursionError) as err:
        raise FormatError([f"{what} is not valid JSON: {err}"]) from None


def load_network(data: bytes, what: str) -> Network:
    raw = load_json(data, what)
    try:
        return Network.from_dict(raw)
    except ValueError as err:
        raise FormatError([f"{what}: {err}"]) from None


# --------------------------------------------------------------------------
# chemart.yaml
# --------------------------------------------------------------------------

def parse_hub_block(raw: Any) -> tuple[dict[str, Any], list[str]]:
    """Validate the `hub:` block; return (block, problems)."""
    if raw is None:
        return {}, []
    if not isinstance(raw, dict):
        return {}, ["hub: must be a mapping"]
    problems = [f"hub: unknown field {k!r}" for k in sorted(set(raw) - _HUB_KEYS, key=str)]
    hub = {k: v for k, v in raw.items() if k in _HUB_KEYS}
    if hub.setdefault("format", FORMAT) != FORMAT:
        problems.append(f"hub.format must be {FORMAT}")
    if "repo_type" in hub and hub["repo_type"] not in REPO_TYPES:
        problems.append(f"hub.repo_type must be one of {list(REPO_TYPES)}")
    for key, limit in (("title", 120), ("description", 1000), ("license", 64)):
        value = hub.get(key)
        if value is not None and (not isinstance(value, str) or len(value) > limit):
            problems.append(f"hub.{key} must be a string of at most {limit} characters")
    tags = hub.get("tags", [])
    if not isinstance(tags, list) or len(tags) > 20 or not all(
        isinstance(t, str) and len(t) <= 40 and _TAG_RE.match(t) for t in tags
    ):
        problems.append("hub.tags must be a list of at most 20 lowercase slugs")
    requires = hub.get("requires", [])
    if not isinstance(requires, list) or not all(isinstance(m, str) and _MODULE_RE.match(m) for m in requires):
        problems.append("hub.requires must be a list of importable module names")
    spec = hub.get("requires_chemart")
    if spec is not None and (not isinstance(spec, str) or not _SPEC_RE.match(spec)):
        problems.append("hub.requires_chemart must look like '>=0.2' or '>=0.2, <1'")
    builtin = hub.get("builtin")
    if builtin is not None and not isinstance(builtin, str):
        problems.append("hub.builtin must be a catalog id")
    return hub, problems


def parse_chemart_yaml(data: bytes) -> tuple[dict[str, Any], catalog.Chemistry | None]:
    """Parse chemart.yaml into its hub block and (for generators) its entry."""
    doc = load_yaml(data)
    if not isinstance(doc, dict):
        raise FormatError([f"{CHEMART_YAML} must be a mapping with 'hub' and/or 'chemistries'"])
    problems = [f"{CHEMART_YAML}: unknown top-level key {k!r}"
                for k in sorted(set(doc) - {"hub", "chemistries"}, key=str)]
    hub, hub_problems = parse_hub_block(doc.get("hub"))
    problems += hub_problems
    entry = None
    if "chemistries" in doc:
        entries = doc["chemistries"]
        if not isinstance(entries, list) or len(entries) != 1:
            problems.append(f"{CHEMART_YAML}: 'chemistries' must be a list with exactly one entry")
        else:
            try:
                entry = catalog.parse_entry(entries[0], CHEMART_YAML)
            except ValueError as err:
                problems.append(f"{CHEMART_YAML}: {err}")
    if problems:
        raise FormatError(problems)
    return hub, entry


# --------------------------------------------------------------------------
# Python sources: parsed, never run
# --------------------------------------------------------------------------

def _python(data: bytes, path: str) -> tuple[ast.Module | None, list[str]]:
    try:
        return ast.parse(_text(data, path), filename=path), []
    except FormatError as err:
        return None, err.problems
    except (SyntaxError, ValueError, RecursionError, MemoryError) as err:
        return None, [f"{path} is not valid Python: {err}"]


def _imports(tree: ast.Module) -> set[str]:
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            out.add(node.module.split(".")[0])
    return out


def _defines_generate(tree: ast.Module) -> bool:
    return any(isinstance(n, ast.FunctionDef) and n.name == "generate" for n in tree.body)


# --------------------------------------------------------------------------
# The whole snapshot
# --------------------------------------------------------------------------

def _size_limit(path: str) -> int:
    if path == "preview.json":
        return _config.MAX_PREVIEW
    if path.endswith(".json"):
        return _config.MAX_NETWORK
    return _config.MAX_TEXT_FILE


def inspect(
    files: Mapping[str, bytes],
    repo: RepoId,
    repo_type: str | None = None,
    *,
    official_namespace: str = "chemart",
    builtin_ids: Collection[str] | None = None,
) -> RepoCard:
    """Check a complete snapshot; return its card or raise FormatError.

    `repo_type` is what the repo was created as (the server knows it); a
    `hub.repo_type` in chemart.yaml must agree. `builtin_ids`, when given,
    is the set of catalog ids a builtin-backed repo may point at.
    """
    problems: list[str] = []
    if not files:
        raise FormatError(["the snapshot is empty"])
    if len(files) > _config.MAX_FILES:
        problems.append(f"too many files: {len(files)} > {_config.MAX_FILES}")
    total = sum(len(b) for b in files.values())
    if total > _config.MAX_REPO:
        problems.append(f"repo too large: {total} bytes > {_config.MAX_REPO}")
    for path, data in files.items():
        try:
            check_path(path)
        except ValueError as err:
            problems.append(str(err))
            continue
        if len(data) > _size_limit(path):
            problems.append(f"{path} is too large: {len(data)} bytes > {_size_limit(path)}")
    if problems:
        raise FormatError(problems)

    hub: dict[str, Any] = {}
    entry = None
    if CHEMART_YAML in files:
        hub, entry = parse_chemart_yaml(files[CHEMART_YAML])
    declared = hub.get("repo_type")
    if repo_type is not None and declared is not None and declared != repo_type:
        raise FormatError([f"chemart.yaml says repo_type {declared!r} but the repo is a {repo_type} repo"])
    repo_type = repo_type or declared
    if repo_type not in REPO_TYPES:
        raise FormatError([f"cannot tell the repo type: set hub.repo_type to one of {list(REPO_TYPES)}"])
    hub["repo_type"] = repo_type

    card = RepoCard(repo_type=repo_type, hub=hub, entry=entry)
    if "README.md" in files:
        try:
            card.readme = _text(files["README.md"], "README.md")
        except FormatError as err:
            problems += err.problems

    sources = {p: d for p, d in files.items() if p.endswith(".py")}
    trees: dict[str, ast.Module] = {}
    for path, data in sources.items():
        tree, errs = _python(data, path)
        problems += errs
        if tree is not None:
            trees[path] = tree
    card.imports = sorted(set().union(*(_imports(t) for t in trees.values())) - {"chemart"}) if trees else []

    if repo_type == "generator":
        problems += _generator_problems(files, repo, card, trees, official_namespace, builtin_ids)
    else:
        problems += _network_problems(files, card, sources)
    if problems:
        raise FormatError(problems)
    return card


def _generator_problems(files, repo, card, trees, official_namespace, builtin_ids) -> list[str]:
    from chemart.api import resolve_params

    problems: list[str] = []
    entry = card.entry
    if entry is None:
        return [f"a generator repo needs {CHEMART_YAML} with one entry under 'chemistries'"]
    problems += [f"{CHEMART_YAML}: {p}" for p in catalog.entry_problems(entry, hub=True)]
    if entry.id != repo.name:
        problems.append(f"the entry id {entry.id!r} must equal the repo name {repo.name!r}")

    builtin = card.builtin
    if builtin is not None:
        if repo.namespace != official_namespace:
            problems.append(f"hub.builtin is reserved for the official {official_namespace}/ repos")
        if builtin != entry.id:
            problems.append("hub.builtin must equal the entry id")
        if builtin_ids is not None and builtin not in builtin_ids:
            problems.append(f"hub.builtin {builtin!r} is not a catalog id of this Chemart")
        if trees:
            problems.append("a builtin-backed repo carries no code: its generator ships with Chemart")
    elif "generator.py" not in files:
        problems.append("a generator repo needs generator.py defining generate(p, rng)")
    elif "generator.py" in trees and not _defines_generate(trees["generator.py"]):
        problems.append("generator.py must define a top-level function generate(p, rng)")

    if "preview.json" not in files:
        problems.append("a generator repo needs preview.json (made by `chemart push`)")
    else:
        try:
            card.network = load_network(files["preview.json"], "preview.json")
        except FormatError as err:
            problems += err.problems
    if card.network is not None and not problems:
        unclaimed = sorted(set(card.network.provides) - set(entry.provides))
        if unclaimed:
            problems.append(f"preview.json contains {unclaimed}, which the entry does not claim in provides")
        try:
            defaults = resolve_params(entry, {})
        except ValueError as err:
            problems.append(f"{CHEMART_YAML}: {err}")
        else:
            if card.network.params != defaults:
                problems.append("preview.json was not made with the default parameters")
    return problems


def _network_problems(files, card, sources) -> list[str]:
    problems: list[str] = []
    if card.entry is not None:
        problems.append(f"a network repo's {CHEMART_YAML} holds a 'hub' block only, no 'chemistries'")
    if card.builtin is not None:
        problems.append("hub.builtin only applies to generator repos")
    if sources:
        problems.append("a network repo carries no code; share code as a generator repo")
    if "network.json" not in files:
        problems.append("a network repo needs network.json (Network.to_dict())")
    else:
        try:
            card.network = load_network(files["network.json"], "network.json")
        except FormatError as err:
            problems += err.problems
    return problems


# --------------------------------------------------------------------------
# Version requirements, without a dependency on `packaging`
# --------------------------------------------------------------------------

def _version_tuple(text: str) -> tuple[int, ...]:
    parts = []
    for piece in text.split("."):
        digits = re.match(r"\d+", piece)
        if not digits:
            break
        parts.append(int(digits.group()))
    return tuple(parts)


def version_satisfies(version: str, spec: str) -> bool:
    """`version_satisfies("0.2.1", ">=0.2, <1")` -> True."""
    have = _version_tuple(version)
    for clause in spec.split(","):
        clause = clause.strip()
        op = re.match(r"(>=|<=|==|!=|>|<)", clause).group()
        want = _version_tuple(clause[len(op):].strip())
        width = max(len(have), len(want))
        a, b = have + (0,) * (width - len(have)), want + (0,) * (width - len(want))
        ok = {">=": a >= b, "<=": a <= b, "==": a == b, "!=": a != b, ">": a > b, "<": a < b}[op]
        if not ok:
            return False
    return True
