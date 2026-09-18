"""Chemart's public interface: three functions, also exposed as LLM tools.

    list_chemistries()                      every catalogued chemistry
    describe_chemistry(id)                  metadata + JSON Schema of its parameters
    generate_network(id, seed=None, **p)    run its generator -> Network

The catalog YAML is the only parameter specification. A chemistry module
`chemart.chemistries.<id with underscores>` defines ``generate(p, rng)``,
receives the validated parameters as attributes of `p`, and returns a
`Network`; this module fills in the provenance.
"""

from __future__ import annotations

import copy
import difflib
from functools import lru_cache
from importlib import import_module
from types import SimpleNamespace
from typing import Any

import numpy as np

from chemart import catalog
from chemart.network import Network


@lru_cache(maxsize=1)
def _entries() -> dict[str, catalog.Chemistry]:
    return {c.id: c for c in catalog.load()}


def _entry(chemistry_id: str) -> catalog.Chemistry:
    entries = _entries()
    if chemistry_id in entries:
        return entries[chemistry_id]
    close = difflib.get_close_matches(str(chemistry_id), entries, n=3)
    hint = f" Did you mean: {', '.join(close)}?" if close else " Call list_chemistries() for all ids."
    raise ValueError(f"unknown chemistry {chemistry_id!r}.{hint}")


def _oneline(value: Any) -> str:
    return " ".join(str(value or "").split())


def _summary(c: catalog.Chemistry, limit: int = 240) -> str:
    text = f"{_oneline(c.S.get('repr'))}. Reaction: {_oneline(c.R.get('scheme'))}"
    return text if len(text) <= limit else text[: limit - 1] + "…"


# ----------------------------------------------------------------------------
def list_chemistries() -> list[dict[str, Any]]:
    """Every chemistry in the catalog: id, name, one-line summary, status."""
    return [
        {
            "id": c.id,
            "name": c.name,
            "summary": _summary(c),
            "implemented": c.implemented,
            "fidelity": c.fidelity,
        }
        for c in _entries().values()
    ]


def params_schema(c: catalog.Chemistry) -> dict[str, Any]:
    return {
        "type": "object",
        "properties": {p.name: p.json_schema() for p in c.params},
        "additionalProperties": False,
    }


def describe_chemistry(chemistry: str) -> dict[str, Any]:
    """Everything known about one chemistry, including its parameter schema."""
    c = _entry(chemistry)
    return {
        "id": c.id,
        "name": c.name,
        "intuition": _oneline(c.intuition) or None,
        "aliases": c.aliases,
        "origin": c.origin,
        "family": c.family,
        "kind": c.kind,
        "constructive": c.constructive,
        "implemented": c.implemented,
        "fidelity": c.fidelity,
        "book": c.book,
        "refs": c.refs,
        "sources": c.sources,
        "decisions": c.decisions,
        "molecules": c.S,
        "reactions": c.R,
        "reactor": c.A,
        "provides": c.provides,
        "phenomena": c.phenomena,
        "notes": _oneline(c.notes) or None,
        "params": params_schema(c),
    }


def resolve_params(c: catalog.Chemistry, given: dict[str, Any]) -> dict[str, Any]:
    """Validate `given` against the catalog spec and fill in defaults."""
    known = {p.name: p for p in c.params}
    for name in given:
        if name not in known:
            close = difflib.get_close_matches(name, known, n=1)
            hint = f" Did you mean {close[0]!r}?" if close else ""
            raise ValueError(
                f"{c.id}: unknown parameter {name!r}.{hint} "
                f"Valid parameters: {', '.join(known) or '(none)'}."
            )
    values: dict[str, Any] = {}
    for name, p in known.items():
        raw = given[name] if name in given else copy.deepcopy(p.default)
        try:
            values[name] = p.coerce(raw)
        except ValueError as err:
            raise ValueError(f"{c.id}: {err}. {name}: {p.meaning}") from None
    return values


def generate_network(chemistry: str, seed: int | None = None, **params: Any) -> Network:
    """Generate a reaction network. Omitted parameters take their defaults."""
    c = _entry(chemistry)
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ValueError(f"seed must be an integer or None, got {seed!r}")
    if not c.implemented:
        raise NotImplementedError(f"{c.id!r} is catalogued but has no generator yet.")
    values = resolve_params(c, params)
    rng = np.random.default_rng(seed)
    net = import_module(c.module).generate(SimpleNamespace(**copy.deepcopy(values)), rng)
    if not isinstance(net, Network):
        raise TypeError(f"{c.module}.generate returned {type(net).__name__}, expected Network")
    net.chemistry, net.params, net.seed = c.id, values, seed
    return net


# ----------------------------------------------------------------------------
def tool_definitions() -> list[dict[str, Any]]:
    """The three functions as JSON-Schema tool specs (name, description, input_schema)."""
    implemented = [c.id for c in _entries().values() if c.implemented]
    chemistry_arg: dict[str, Any] = {"type": "string", "description": "chemistry id from list_chemistries"}
    if implemented:
        chemistry_arg = {**chemistry_arg, "enum": implemented}
    return [
        {
            "name": "list_chemistries",
            "description": "List every artificial chemistry in the Chemart catalog with its id, "
                           "name, a one-line summary and whether a generator is implemented.",
            "input_schema": {"type": "object", "properties": {}},
        },
        {
            "name": "describe_chemistry",
            "description": "Describe one artificial chemistry: molecules, reaction scheme, book and "
                           "paper sources, and the JSON Schema of its parameters.",
            "input_schema": {
                "type": "object",
                "properties": {"chemistry": {"type": "string", "description": "chemistry id"}},
                "required": ["chemistry"],
            },
        },
        {
            "name": "generate_network",
            "description": "Generate the chemical reaction network of one artificial chemistry. "
                           "Returns species, reactions (reactant and product stoichiometry, optional "
                           "rate law, optional firing count), status and the parameters used. "
                           "All parameters are optional; see describe_chemistry.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "chemistry": chemistry_arg,
                    "seed": {"type": "integer", "description": "random seed for reproducibility"},
                    "params": {"type": "object", "description": "parameter values by name"},
                },
                "required": ["chemistry"],
            },
        },
    ]


def call_tool(name: str, arguments: dict[str, Any] | None = None) -> Any:
    """Execute a tool call from `tool_definitions` and return JSON-ready data."""
    arguments = arguments or {}
    if name == "list_chemistries":
        return list_chemistries()
    if name == "describe_chemistry":
        return describe_chemistry(arguments["chemistry"])
    if name == "generate_network":
        net = generate_network(
            arguments["chemistry"], arguments.get("seed"), **arguments.get("params", {})
        )
        return net.to_dict()
    raise ValueError(f"unknown tool {name!r}")
