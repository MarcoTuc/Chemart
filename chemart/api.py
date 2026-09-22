"""Chemart's public interface, also exposed as LLM tools.

    list_chemistries()                      every catalogued chemistry
    describe_chemistry(id)                  metadata + JSON Schema of its parameters
    generate_network(id, seed=None, **p)    one network -> Network
    evolve(id, seed=None, **p)              a run of its process -> Trajectory

The catalog YAML is the only parameter specification. A chemistry module
`chemart.chemistries.<id with underscores>` has one or two faces:
``generate(p, rng)`` returns a `Network`, and ``evolve(p, rng)`` is a
generator that yields `Frame`s of a process and returns the observed
`Network` at the end. Both receive the validated parameters as attributes of
`p`; this module fills in the provenance. A parameter used by one face only
says so with `face` in the catalog.

An id of the form ``namespace/name`` (optionally ``@revision``) names a
chemistry shared on the Chemart Hub instead; see `chemart.hub`. Code from
the hub runs only when the caller passes ``trust_remote_code=True``.
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
from chemart.trajectory import Frame, Trajectory


@lru_cache(maxsize=1)
def _entries() -> dict[str, catalog.Chemistry]:
    return {c.id: c for c in catalog.load()}


def is_hub_id(chemistry_id: str) -> bool:
    """Built-in ids are bare slugs; hub ids are ``namespace/name[@revision]``."""
    return isinstance(chemistry_id, str) and ("/" in chemistry_id or "@" in chemistry_id)


def _entry(chemistry_id: str, revision: str | None = None) -> catalog.Chemistry:
    if is_hub_id(chemistry_id):
        from chemart import hub

        return hub.resolve_chemistry(chemistry_id, revision)
    if revision is not None:
        raise ValueError(
            f"revision={revision!r} only applies to hub ids (namespace/name); "
            f"{chemistry_id!r} is a built-in chemistry"
        )
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
def list_chemistries(include_archived: bool = False) -> list[dict[str, Any]]:
    """Every chemistry in the catalog: id, name, one-line summary, status.

    Archived entries (see `catalog.ARCHIVES`) are left out unless
    `include_archived` is true; they stay available by id either way.
    """
    rows = []
    for c in _entries().values():
        if c.archived is not None and not include_archived:
            continue
        row = {
            "id": c.id,
            "name": c.name,
            "summary": _summary(c),
            "implemented": c.implemented,
            "fidelity": c.fidelity,
            "type": c.type,
        }
        if include_archived:
            row["archived"] = c.archived
        rows.append(row)
    return rows


def params_schema(c: catalog.Chemistry, face: str | None = None) -> dict[str, Any]:
    """The JSON Schema of a chemistry's parameters, or of one face's."""
    return {
        "type": "object",
        "properties": {p.name: p.json_schema() for p in c.params if face is None or p.face in (None, face)},
        "additionalProperties": False,
    }


def describe_chemistry(chemistry: str, revision: str | None = None) -> dict[str, Any]:
    """Everything known about one chemistry, including its parameter schema.

    For a hub id this reads the repo's metadata only; no code is downloaded
    or run, so it needs no ``trust_remote_code``.
    """
    c = _entry(chemistry, revision)
    info = {
        "id": c.id,
        "name": c.name,
        "intuition": _oneline(c.intuition) or None,
        "aliases": c.aliases,
        "origin": c.origin,
        "family": c.family,
        "archived": c.archived,
        "kind": c.kind,
        "type": c.type,
        "faces": faces(c),
        "clock": c.clock,
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
        "params": params_schema(c, _face_of_generate(c)) if c.implemented else params_schema(c),
    }
    if "evolve" in info["faces"]:
        info["evolve_params"] = params_schema(c, "evolve")
    hub_info = getattr(c, "hub_info", None)
    if hub_info is not None:
        info["hub"] = hub_info()
    return info


#: Parameters that chose between a closure and a soup before the two became faces.
_REMOVED_MODES = ("method", "mode")


def faces(c: catalog.Chemistry) -> list[str]:
    """How the chemistry can be run: "generate" and/or "evolve", from its module."""
    if getattr(c, "load_generator", None) is not None:
        return ["generate"]                  # hub chemistries have a generate face only
    if not c.implemented:
        return []
    module = import_module(c.module)
    return [f for f in ("generate", "evolve") if callable(getattr(module, f, None))]


def _face_of_generate(c: catalog.Chemistry) -> str:
    """generate_network runs the generate face, or the evolve face of a gas without one."""
    return "evolve" if faces(c) == ["evolve"] else "generate"


def resolve_params(c: catalog.Chemistry, given: dict[str, Any], face: str | None = None) -> dict[str, Any]:
    """Validate `given` against the catalog spec and fill in defaults.

    With `face`, only the parameters of that face (and the shared ones) are
    accepted and returned.
    """
    known = {p.name: p for p in c.params if face is None or p.face in (None, face)}
    for name in given:
        if name not in known:
            other = next((p for p in c.params if p.name == name), None)
            if other is not None:
                call = "chemart.evolve" if other.face == "evolve" else "chemart.generate_network"
                raise ValueError(f"{c.id}: parameter {name!r} belongs to the {other.face} face; "
                                 f"pass it to {call}")
            if name in _REMOVED_MODES and "evolve" in faces(c):
                raise ValueError(f"{c.id}: parameter {name!r} is gone: generate_network returns the "
                                 "network (the closure) and chemart.evolve runs the process (the soup)")
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


def generate_network(
    chemistry: str,
    seed: int | None = None,
    *,
    revision: str | None = None,
    trust_remote_code: bool = False,
    **params: Any,
) -> Network:
    """Generate a reaction network. Omitted parameters take their defaults.

    `chemistry` is a catalog id (``"brusselator"``) or a hub id
    (``"alice/my-chem"``, ``"alice/my-chem@3f2a9c1"``). A hub chemistry that
    ships its own code runs only with ``trust_remote_code=True``; pin
    `revision` to the commit you reviewed.
    """
    c = _entry(chemistry, revision)
    if not c.implemented:
        raise NotImplementedError(f"{c.id!r} is catalogued but has no generator yet.")
    return run_generator(c, generator_for(c, trust_remote_code), seed, params)


def generator_for(c: catalog.Chemistry, trust_remote_code: bool = False):
    """The ``generate(p, rng)`` function of an entry, built-in or from the hub.

    For a gas with no generate face it runs ``evolve`` to the end and returns
    the observed network.
    """
    load = getattr(c, "load_generator", None)
    if load is not None:
        return load(trust_remote_code)
    module = import_module(c.module)
    if callable(getattr(module, "generate", None)):
        return module.generate

    def generate(p, rng):
        frames = module.evolve(p, rng)
        while True:
            try:
                next(frames)
            except StopIteration as stop:
                return stop.value

    return generate


def _check_seed(seed) -> None:
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ValueError(f"seed must be an integer or None, got {seed!r}")


def run_generator(c: catalog.Chemistry, generate, seed: int | None, params: dict[str, Any]) -> Network:
    """Validate `params` against `c`, call `generate(p, rng)`, record provenance."""
    _check_seed(seed)
    values = resolve_params(c, params, _face_of_generate(c))
    rng = np.random.default_rng(seed)
    net = generate(SimpleNamespace(**copy.deepcopy(values)), rng)
    if not isinstance(net, Network):
        raise TypeError(f"{c.id}: generate returned {type(net).__name__}, expected Network")
    net.chemistry, net.params, net.seed = getattr(c, "ref", c.id), values, seed
    return net


def evolver_for(c: catalog.Chemistry):
    """The ``evolve(p, rng)`` function of a built-in chemistry that has one."""
    if getattr(c, "load_generator", None) is not None:
        raise ValueError(f"{c.id}: chemistries from the hub can only be generated, not evolved")
    if not c.implemented:
        raise NotImplementedError(f"{c.id!r} is catalogued but has no generator yet.")
    evolve_ = getattr(import_module(c.module), "evolve", None)
    if not callable(evolve_):
        raise ValueError(f"{c.id} is a {c.type} with no process to evolve; "
                         "use generate_network (and chemart.simulate for its dynamics)")
    return evolve_


def run_evolver(c: catalog.Chemistry, evolve_, seed: int | None, params: dict[str, Any], every: int = 1):
    """Validate `params`, run `evolve(p, rng)` and yield its frames, `every`
    frames merged into one; return the observed network with its provenance."""
    _check_seed(seed)
    if isinstance(every, bool) or not isinstance(every, int) or every < 1:
        raise ValueError(f"every must be a positive integer, got {every!r}")
    values = resolve_params(c, params, "evolve")
    rng = np.random.default_rng(seed)
    frames = evolve_(SimpleNamespace(**copy.deepcopy(values)), rng)
    pending: list[Frame] = []
    first = True
    while True:
        try:
            frame = next(frames)
        except StopIteration as stop:
            net = stop.value
            break
        if first:
            first = False
            yield frame
            continue
        pending.append(frame)
        if len(pending) == every:
            yield _merge(pending)
            pending = []
    if pending:
        yield _merge(pending)
    if not isinstance(net, Network):
        raise TypeError(f"{c.id}: evolve returned {type(net).__name__}, expected Network")
    net.chemistry, net.params, net.seed = getattr(c, "ref", c.id), values, seed
    return net


def _merge(frames: list[Frame]) -> Frame:
    """Consecutive frames as one: the reactions add up, the last state stands."""
    if len(frames) == 1:
        return frames[0]
    from collections import Counter

    fired: dict[tuple, list] = {}
    for f in frames:
        for lhs, rhs, n in f.fired:
            key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
            entry = fired.get(key)
            if entry is None:
                fired[key] = [list(lhs), list(rhs), n]
            else:
                entry[2] += n
    last = frames[-1]
    return Frame(t=last.t, state=last.state, fired=list(fired.values()), observables=last.observables)


def evolve_frames(chemistry: str, seed: int | None = None, *, every: int = 1, **params: Any):
    """Run a chemistry's process, yielding its frames as they come.

    A generator: iterate it for live frames; its return value (``yield from``
    or ``StopIteration.value``) is the observed network.
    """
    c = _entry(chemistry)
    return run_evolver(c, evolver_for(c), seed, params, every)


def evolve(chemistry: str, seed: int | None = None, *, every: int = 1, **params: Any) -> Trajectory:
    """Run a chemistry's process to the end and return the whole `Trajectory`.

    Every chemistry with an evolve face can be evolved: the Turing gases, and
    the chemistries that define a reactor for their network (such as a
    lattice). Omitted parameters take their defaults; `every` keeps one frame
    in `every`, adding up the reactions fired in between.
    """
    c = _entry(chemistry)
    run = run_evolver(c, evolver_for(c), seed, params, every)
    frames: list[Frame] = []
    while True:
        try:
            frames.append(next(run))
        except StopIteration as stop:
            net = stop.value
            break
    return Trajectory(network=net, frames=frames, method="evolve", clock=c.clock or "steps",
                      settings={"seed": seed, "params": net.params, "every": every})


# ----------------------------------------------------------------------------
def tool_definitions() -> list[dict[str, Any]]:
    """The three functions as JSON-Schema tool specs (name, description, input_schema)."""
    implemented = [c.id for c in _entries().values() if c.implemented and c.archived is None]
    chemistry_arg: dict[str, Any] = {"type": "string", "description": "chemistry id from list_chemistries"}
    if implemented:
        chemistry_arg = {**chemistry_arg, "enum": implemented}
    return [
        {
            "name": "list_chemistries",
            "description": "List every artificial chemistry in the Chemart catalog with its id, "
                           "name, a one-line summary and whether a generator is implemented.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "include_archived": {
                        "type": "boolean",
                        "description": "also list archived entries, which are not part of the "
                                       "chemistry catalog (default false)",
                    },
                },
            },
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
        return list_chemistries(bool(arguments.get("include_archived", False)))
    if name == "describe_chemistry":
        return describe_chemistry(arguments["chemistry"])
    if name == "generate_network":
        # A tool call never runs code from the hub: a model that was talked
        # into it must not be able to switch the trust gate on via params.
        params = dict(arguments.get("params") or {})
        reserved = sorted(set(params) & catalog.RESERVED_PARAMS)
        if reserved:
            raise ValueError(f"{', '.join(reserved)} cannot be passed as chemistry parameters")
        net = generate_network(arguments["chemistry"], arguments.get("seed"), **params)
        return net.to_dict()
    raise ValueError(f"unknown tool {name!r}")
