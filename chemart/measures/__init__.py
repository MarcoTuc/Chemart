"""Measures on reaction networks, populations and trajectories.

    import chemart
    net = chemart.generate_network("brusselator")
    chemart.measure(net)                           # every cheap measure that applies
    chemart.measure(net, ["deficiency", "bow_tie"])
    chemart.measures.applicable(net)               # {name: None or why not}

    traj = chemart.evolve("alchemy", seed=0)
    chemart.measures.over(traj, ["richness", "shannon"])   # measures over time

    rows = chemart.measures.sweep("random-catalytic-networks", {"n": [10, 20, 40]}, seeds=range(5))
    chemart.measures.scaling(rows, "n_reactions")          # exponent of n_reactions ~ n_species^a

Every measure is registered here with the section of docs/guide/measures.md it
belongs to, what it needs from a network (the codes of that page: T topology,
S stoichiometry, C catalysts, F a food set, K rate constants, D a trajectory,
str molecular structure), what it takes as input (a network, the population
state of one frame, or a whole trajectory) and its cost. A measure that does
not apply is left out of the result, never filled in; `applicable` says why.
"""

from __future__ import annotations

import math
from collections import Counter, deque
from dataclasses import dataclass
from functools import cached_property
from itertools import product
from typing import Any, Callable, Iterable

import numpy as np

from chemart.network import Network, Reaction, Species

SECTIONS = {
    "A": "size and composition", "B": "stoichiometric structure", "C": "graph topology",
    "D": "catalysis, autocatalysis and organisation", "E": "constructiveness and growth",
    "F": "kinetics", "G": "dynamics and trajectories", "H": "robustness and redundancy",
    "I": "information and complexity", "J": "scaling",
}
COSTS = ("cheap", "moderate", "exponential")
INPUTS = ("network", "state", "trajectory")
NEEDS = {"T", "S", "C", "F", "K", "D", "str"}


@dataclass(frozen=True)
class Measure:
    name: str
    section: str
    fn: Callable
    input: str
    needs: frozenset
    cost: str
    limit: int | None
    doc: str


REGISTRY: dict[str, Measure] = {}


def register(name: str, section: str, *, input: str = "network", needs: str = "T",
             cost: str = "cheap", limit: int | None = None):
    """Add a measure to the registry. `limit` caps species + reactions for a
    network measure (above it the measure does not run unless forced)."""
    tags = frozenset(needs.split())
    if section not in SECTIONS or input not in INPUTS or cost not in COSTS or not tags <= NEEDS:
        raise ValueError(f"bad registration of measure {name!r}")

    def wrap(fn):
        doc = " ".join((fn.__doc__ or "").strip().split("\n\n")[0].split())
        REGISTRY[name] = Measure(name, section, fn, input, tags, cost, limit, doc)
        return fn

    return wrap


# --------------------------------------------------------------------------
# What a network gives the measures, computed once
# --------------------------------------------------------------------------

class Context:
    """A network and the derived objects the measures share, built on first use."""

    def __init__(self, net: Network, *, food: Iterable[str] | None = None, seed: int = 0):
        self.net = net
        self._food = None if food is None else list(food)
        self.seed = seed

    @cached_property
    def ids(self) -> list[str]:
        return [s.id for s in self.net.species]

    @cached_property
    def index(self) -> dict[str, int]:
        return {s: i for i, s in enumerate(self.ids)}

    @cached_property
    def S(self) -> np.ndarray:
        """The stoichiometric matrix, species x reactions (dense)."""
        _, R, P = self.net.matrices()
        return (P - R).toarray()

    @cached_property
    def rank(self) -> int:
        """rank(S), exactly: stoichiometries are integers, so the rank is found by
        sparse elimination modulo a large prime (it can only differ from the
        rational rank if the prime divides every relevant minor, which small
        integer coefficients never allow). Linear in practice, where a dense SVD
        of a thousands-by-thousands S takes many seconds."""
        return sparse_rank(self.net)

    @cached_property
    def graph(self):
        """The bipartite species-reaction digraph: s -> r for reactants, r -> s for products."""
        import networkx as nx

        g = nx.DiGraph()
        g.add_nodes_from((("s", s) for s in self.ids), kind="species")
        for j, r in enumerate(self.net.reactions):
            g.add_node(("r", j), kind="reaction")
            g.add_edges_from((("s", s), ("r", j)) for s in r.reactants)
            g.add_edges_from((("r", j), ("s", s)) for s in r.products)
        return g

    @cached_property
    def species_graph(self):
        """Substrate -> product digraph over species (self-loops left out)."""
        import networkx as nx

        g = nx.DiGraph()
        g.add_nodes_from(self.ids)
        for r in self.net.reactions:
            g.add_edges_from((a, b) for a in r.reactants for b in r.products if a != b)
        return g

    @cached_property
    def food(self) -> list[str]:
        return self._food if self._food is not None else food_set(self.net)

    def rng(self) -> np.random.Generator:
        return np.random.default_rng(self.seed)


_PRIME = 2_147_483_647


def sparse_rank(net: Network) -> int:
    """Rank of the stoichiometric matrix of `net`, reaction by reaction."""
    index = {s.id: i for i, s in enumerate(net.species)}
    basis: dict[int, dict[int, int]] = {}      # leading species -> vector with leading entry 1
    for r in net.reactions:
        v: dict[int, int] = {}
        for s, n in r.reactants.items():
            v[index[s]] = v.get(index[s], 0) - n
        for s, n in r.products.items():
            v[index[s]] = v.get(index[s], 0) + n
        v = {i: x % _PRIME for i, x in v.items() if x % _PRIME}
        while v:
            lead = min(v)
            b = basis.get(lead)
            if b is None:
                inv = pow(v[lead], _PRIME - 2, _PRIME)
                basis[lead] = {i: x * inv % _PRIME for i, x in v.items()}
                break
            f = v[lead]
            for i, x in b.items():
                y = (v.get(i, 0) - f * x) % _PRIME
                if y:
                    v[i] = y
                else:
                    v.pop(i, None)
    return len(basis)


def food_set(net: Network) -> list[str]:
    """What a network is fed with: extras["food"], else its inflow, else its
    buffered species, else its initial state. Empty when it says nothing."""
    known = {s.id for s in net.species}
    for source in (net.extras.get("food"), net.inflow, net.extras.get("buffered"), net.initial_state):
        if source:
            return [s for s in source if s in known]
    return []


def _has(net: Network, tag: str, ctx: Context | None) -> bool:
    if tag in ("T", "S"):
        return bool(net.reactions)
    if tag == "C":
        return any(r.catalysts for r in net.reactions)
    if tag == "K":
        return bool(net.reactions) and all(r.rate for r in net.reactions)
    if tag == "F":
        return bool((ctx.food if ctx else food_set(net)))
    if tag == "str":
        return any(s.structure for s in net.species)
    return False


_WHY = {
    "T": "the network has no reactions", "S": "the network has no reactions",
    "C": "no reaction has a catalyst", "K": "not every reaction has a rate",
    "F": "no food set (extras['food'], inflow, buffered species or initial state)",
    "str": "species carry no structure", "D": "needs a trajectory",
}


def _why_not(m: Measure, obj: Any, ctx: Context | None, force: bool) -> str | None:
    if m.input == "state":
        return None if isinstance(obj, dict) else "measures a population state (a frame's state)"
    if m.input == "trajectory":
        from chemart.trajectory import Trajectory

        return None if isinstance(obj, Trajectory) else _WHY["D"]
    net = obj.network if hasattr(obj, "frames") else obj
    if not isinstance(net, Network):
        return "measures a network"
    for tag in sorted(m.needs - {"D"}):
        if not _has(net, tag, ctx):
            return _WHY[tag]
    size = len(net.species) + len(net.reactions)
    if m.limit is not None and size > m.limit and not force:
        return f"{size} species and reactions exceed its limit of {m.limit} (pass force=True)"
    return None


def _select(names, cost: str, inputs: tuple[str, ...]) -> list[Measure]:
    if cost not in COSTS:
        raise ValueError(f"cost must be one of {COSTS}, got {cost!r}")
    if names is None:
        top = COSTS.index(cost)
        return [m for m in REGISTRY.values() if m.input in inputs and COSTS.index(m.cost) <= top]
    if isinstance(names, str):
        names = [names]
    unknown = [n for n in names if n not in REGISTRY]
    if unknown:
        import difflib

        hints = {n: difflib.get_close_matches(n, REGISTRY, n=1) for n in unknown}
        raise ValueError("unknown measure(s): " + ", ".join(
            f"{n!r}" + (f" (did you mean {h[0]!r}?)" if h else "") for n, h in hints.items()))
    return [REGISTRY[n] for n in names]


def _inputs(obj) -> tuple[str, ...]:
    if isinstance(obj, dict):
        return ("state",)
    if hasattr(obj, "frames"):
        return ("trajectory", "network")
    return ("network",)


# --------------------------------------------------------------------------
# The API
# --------------------------------------------------------------------------

def measure(obj, names=None, *, cost: str = "cheap", food: Iterable[str] | None = None,
            seed: int = 0, force: bool = False) -> dict[str, Any]:
    """Compute measures of a network, a trajectory or a population state.

    `names` picks measures by name; by default every registered measure up to
    `cost` ("cheap", "moderate" or "exponential") that takes this kind of
    input. A trajectory gets its trajectory measures plus the network measures
    of its (cumulative) network. Measures that do not apply are left out;
    `applicable` says why. `food` overrides the food set; `seed` drives the
    measures that sample.
    """
    net = obj.network if hasattr(obj, "frames") else obj
    ctx = Context(net, food=food, seed=seed) if isinstance(net, Network) else None
    out: dict[str, Any] = {}
    for m in _select(names, cost, _inputs(obj)):
        if _why_not(m, obj, ctx, force) is not None:
            continue
        if m.input == "network":
            out[m.name] = m.fn(ctx)
        else:
            out[m.name] = m.fn(obj)
    return out


def applicable(obj, names=None, *, cost: str = "exponential", food=None, force: bool = False) -> dict[str, str | None]:
    """For each measure: None if it applies to `obj`, else the reason it does not."""
    net = obj.network if hasattr(obj, "frames") else obj
    ctx = Context(net, food=food) if isinstance(net, Network) else None
    return {m.name: _why_not(m, obj, ctx, force) for m in _select(names, cost, INPUTS)}


def describe() -> list[dict[str, Any]]:
    """Every registered measure: name, section, input, needs, cost, limit, meaning."""
    return [{"name": m.name, "section": m.section, "input": m.input, "needs": sorted(m.needs),
             "cost": m.cost, "limit": m.limit, "meaning": m.doc} for m in REGISTRY.values()]


def fired_network(fired: Iterable[list], like: Network | None = None) -> Network:
    """The observed network of a list of [[reactants], [products], count] entries."""
    reactions: dict[tuple, list] = {}
    for lhs, rhs, n in fired:
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        reactions.setdefault(key, [lhs, rhs, 0])[2] += n
    names = dict.fromkeys(s for lhs, rhs, _ in reactions.values() for s in (*lhs, *rhs))
    known = {s.id: s for s in like.species} if like is not None else {}
    return Network(species=[known.get(s, Species(s)) for s in names],
                   reactions=[Reaction.of(lhs, rhs, count=n) for lhs, rhs, n in reactions.values()],
                   status="observed")


def track(frames: Iterable, names=None, *, window: int | None = 1, cost: str = "cheap",
          food=None, seed: int = 0) -> Iterable[dict[str, Any]]:
    """Measures frame by frame, for a trajectory or a live run (chemart.evolve_frames).

    Population measures read each frame's state; network measures read the
    network of the reactions fired in the last `window` frames (None: every
    frame so far). Yields {"t": ..., name: value, ...} per frame.
    """
    chosen = _select(names, cost, ("state", "network"))
    whole = [m.name for m in chosen if m.input == "trajectory"]
    if whole:
        raise ValueError(f"{', '.join(whole)}: measured on a whole run, not frame by frame; "
                         "pass the trajectory to measure() instead")
    state_measures = [m for m in chosen if m.input == "state"]
    net_measures = [m for m in chosen if m.input == "network"]
    recent: deque = deque(maxlen=window)
    total: list = []
    for frame in frames:
        row: dict[str, Any] = {"t": frame.t}
        for m in state_measures:
            row[m.name] = m.fn(frame.state)
        if net_measures:
            if window is None:
                total.extend(frame.fired)
                fired = total
            else:
                recent.append(frame.fired)
                fired = [entry for f in recent for entry in f]
            if fired:
                net = fired_network(fired)
                ctx = Context(net, food=food, seed=seed)
                for m in net_measures:
                    if _why_not(m, net, ctx, False) is None:
                        row[m.name] = m.fn(ctx)
        yield row


def over(traj, names=None, *, window: int | None = 1, cost: str = "cheap", food=None,
         seed: int = 0) -> dict[str, list]:
    """`track` over a finished trajectory, as {"t": [...], name: [...]} (None where
    a measure did not apply to a frame)."""
    rows = list(track(traj.frames, names, window=window, cost=cost, food=food, seed=seed))
    keys = list(dict.fromkeys(k for row in rows for k in row))
    return {k: [row.get(k) for row in rows] for k in keys}


def randomize(net: Network, rng: np.random.Generator, swaps_per_edge: int = 10) -> Network:
    """A null-model copy: reactant and product slots swapped between reactions,
    keeping every species' number of reactant and product slots, every
    reaction's arity and multiplicities (Maslov-Sneppen swaps on each side)."""
    sides = {"reactants": [dict(r.reactants) for r in net.reactions],
             "products": [dict(r.products) for r in net.reactions]}
    for side in sides.values():
        slots = [(j, s) for j, d in enumerate(side) for s in d]
        for _ in range(swaps_per_edge * len(slots)):
            if len(slots) < 2:
                break
            a, b = rng.integers(len(slots), size=2)
            (j1, s1), (j2, s2) = slots[a], slots[b]
            if j1 == j2 or s1 == s2 or s2 in side[j1] or s1 in side[j2] or side[j1][s1] != side[j2][s2]:
                continue
            n = side[j1].pop(s1)
            side[j2].pop(s2)
            side[j1][s2] = n
            side[j2][s1] = n
            slots[a], slots[b] = (j1, s2), (j2, s1)
    reactions = [Reaction(reactants=lhs, products=rhs) for lhs, rhs in zip(sides["reactants"], sides["products"])
                 if lhs or rhs]
    return Network(species=list(net.species), reactions=reactions, extras={
        k: v for k, v in net.extras.items() if k in ("food", "buffered")},
        initial_state=net.initial_state, inflow=net.inflow)


def zscores(net: Network, names=None, *, samples: int = 20, seed: int = 0, cost: str = "cheap") -> dict[str, float | None]:
    """How far each scalar measure lies from its null model, in standard deviations.

    None when the null model does not vary (the measure is fixed by the degrees).
    """
    rng = np.random.default_rng(seed)
    real = {k: v for k, v in measure(net, names, cost=cost, seed=seed).items() if _scalar(v)}
    null = [measure(randomize(net, rng), list(real), seed=seed) for _ in range(samples)]
    out = {}
    for name, x in real.items():
        values = np.array([float(n[name]) for n in null if _scalar(n.get(name))])
        sd = values.std(ddof=1) if len(values) > 1 else 0.0
        out[name] = None if sd == 0 else float((float(x) - values.mean()) / sd)
    return out


def _scalar(v) -> bool:
    return isinstance(v, (int, float, bool, np.integer, np.floating)) and not (
        isinstance(v, float) and math.isnan(v))


def sweep(chemistry: str, grid: dict[str, list], *, seeds: Iterable[int] = range(3), names=None,
          cost: str = "cheap", **fixed) -> list[dict[str, Any]]:
    """Measure a chemistry's networks over a grid of its arguments and seeds.

    One row per combination and seed: the arguments, the seed, and every
    measure (dict-valued measures flattened as name.key).
    """
    from chemart.api import generate_network

    keys = list(grid)
    rows = []
    for combo in product(*(grid[k] for k in keys)):
        args = dict(zip(keys, combo))
        for s in seeds:
            net = generate_network(chemistry, seed=s, **fixed, **args)
            row = {**args, "seed": s}
            for name, value in measure(net, names, cost=cost, seed=s).items():
                if isinstance(value, dict):
                    row.update({f"{name}.{k}": v for k, v in value.items()})
                else:
                    row[name] = value
            rows.append(row)
    return rows


def scaling(rows: list[dict[str, Any]], y: str, x: str = "n_species") -> dict[str, float | None]:
    """The exponent a of y ~ x^a, fitted on log-log axes over the rows where
    both are positive: {"exponent", "intercept", "r2", "n"}."""
    pairs = [(float(r[x]), float(r[y])) for r in rows
             if _scalar(r.get(x)) and _scalar(r.get(y)) and float(r[x]) > 0 and float(r[y]) > 0]
    if len(pairs) < 2 or len({p[0] for p in pairs}) < 2:
        return {"exponent": None, "intercept": None, "r2": None, "n": len(pairs)}
    lx, ly = np.log([p[0] for p in pairs]), np.log([p[1] for p in pairs])
    slope, intercept = np.polyfit(lx, ly, 1)
    fit = slope * lx + intercept
    ss = ((ly - ly.mean()) ** 2).sum()
    r2 = 1.0 - ((ly - fit) ** 2).sum() / ss if ss else 1.0
    return {"exponent": float(slope), "intercept": float(intercept), "r2": float(r2), "n": len(pairs)}


# The measures register themselves on import.
from chemart.measures import (  # noqa: E402,F401
    dynamics, graph, growth, information, kinetics, organisation, robustness, size, stoichiometry,
)

# Section order, whatever the import order: measure() and describe() list them so.
_ordered = sorted(REGISTRY.items(), key=lambda item: item[1].section)
REGISTRY.clear()
REGISTRY.update(_ordered)
