"""Molecular traveling salesman: Banzhaf (1990), book 17.2.1, PyCellChemistry MolecularTSP.py.

Catalog id: molecular-tsp.

Data strings are candidate tours s = (s0, s1..sN) with s0 = l(s), the closed
euclidean tour length (book eq. 17.1). Four fixed machines float in the soup
and act as catalysts:

    E-machine + t      -> E-machine + t'        exchange two cities
    C-machine + t      -> C-machine + t'        cut a segment, reinsert it behind a third city
    I-machine + t      -> I-machine + t'        as C, with the segment inverted
    R-machine + t1 + t2 -> R-machine + (the 2 best of t1, t2, child)

A machine releases the n_op best strings (n_op = 1 for E, C, I; 2 for R): a
single-string machine keeps its trial tour only if it is strictly shorter.

The run follows the authors' re-implementation (PyCellChemistry
MolecularTSP.py on HighOrderChem): each operation cycle picks a machine with
probability proportional to its time scale t_j (m_j = 1 machine of each sort),
draws n_op distinct strings uniformly from the M strings, and puts the released
strings back. A generation is c = ceil(M / sum_j t_j) cycles. The chemistry
is a gas with one face, ``evolve``: a frame per generation, observing the best
and mean tour length and the population overlap. The observed network holds
the reactions that changed the soup, with firing counts; species ids are
canonical tours (rotation and direction removed), while the soup keeps the
oriented strings the operators act on.
"""

from __future__ import annotations

import math
from collections import Counter

from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame

MACHINES = ("E", "C", "I", "R")
N_OP = {"E": 1, "C": 1, "I": 1, "R": 2}
DESCRIPTION = {
    "E": "exchange machine: swaps two random cities (n_op = 1)",
    "C": "cutting machine: moves the segment between two random cities behind a third (n_op = 1)",
    "I": "inversion machine: as C, with the moved segment inverted (n_op = 1)",
    "R": "recombination machine: grafts a random segment of the first tour into the second (n_op = 2)",
}


# ---------------------------------------------------------------------------
# Problem instance
# ---------------------------------------------------------------------------
def ring_cities(n: int) -> list[list[float]]:
    """TSPgraph.createRingTopo: n cities on a circle of radius n centred at (n, n)."""
    radius = n  # gridsize 2n, radius gridsize / 2
    return [[math.cos(2 * math.pi * i / n) * radius + radius,
             math.sin(2 * math.pi * i / n) * radius + radius] for i in range(n)]


def random_cities(n: int, rng) -> list[list[float]]:
    """TSPgraph.createRndTopo: integer points on a 2n x 2n grid, at least 2 apart."""
    grid, mindist = 2 * n, 2.0
    out: list[list[float]] = []
    while len(out) < n:
        x, y = int(rng.integers(0, grid)), int(rng.integers(0, grid))
        if all(math.hypot(x - a, y - b) >= mindist for a, b in out):
            out.append([float(x), float(y)])
    return out


def ring_optimum(n: int) -> float:
    """Perimeter of the regular n-gon of radius n: the optimal tour of the ring layout."""
    return n * 2 * n * math.sin(math.pi / n)


def distance_matrix(cities) -> list[list[float]]:
    return [[math.hypot(xa - xb, ya - yb) for xb, yb in cities] for xa, ya in cities]


def tour_length(tour, dist) -> float:
    """Book eq. 17.1 / paper eq. 2: closed tour, s_{N+1} = s_1."""
    return sum(dist[a][b] for a, b in zip(tour, tour[1:] + tour[:1]))


def canonical(tour) -> tuple[int, ...]:
    """The tour as a cycle: start at the smallest city, walk towards the smaller neighbour."""
    i = tour.index(min(tour))
    t = list(tour[i:]) + list(tour[:i])
    if len(t) > 2 and t[-1] < t[1]:
        t = t[:1] + t[:0:-1]
    return tuple(t)


def tour_id(tour) -> str:
    return "t" + ".".join(str(c) for c in canonical(tour))


def machine_id(kind: str) -> str:
    return f"{kind}-machine"


# ---------------------------------------------------------------------------
# Machine operators (deterministic given positions; ports of MolecularTSP.py)
# ---------------------------------------------------------------------------
def split_tour(tour, p1: int, p2: int):
    """MolecularTSP.splitTour: [p1, p2) and [p2, p1), circularly; equal positions give (tour, [])."""
    n = len(tour)
    if p1 % n == p2 % n:
        return list(tour), []
    seg1, i = [], p1
    while i != p2:
        seg1.append(tour[i])
        i = (i + 1) % n
    seg2, i = [], p2
    while i != p1:
        seg2.append(tour[i])
        i = (i + 1) % n
    return seg1, seg2


def exchange(tour, p1: int, p2: int) -> list:
    """E-machine operator: swap the cities at positions p1 and p2."""
    new = list(tour)
    new[p1], new[p2] = tour[p2], tour[p1]
    return new


def cut(tour, p1: int, p2: int, p3: int, invert: bool = False) -> list:
    """C- (invert=False) and I-machine (invert=True) operator.

    The segment [p1, p2) is removed; the rest, read from p2, is split after its
    first p3 cities (1 <= p3 < len(rest)) and the segment, inverted for the
    I-machine, is inserted there.
    """
    seg1, seg2 = split_tour(tour, p1, p2)
    if len(seg2) > 1:
        seg3, seg4 = seg2[:p3], seg2[p3:]
    else:
        seg3, seg4 = seg2, []
    if invert:
        seg1 = seg1[::-1]
    return seg3 + seg1 + seg4


def recombine(tour1, tour2, p1: int, p2: int) -> list:
    """R-machine operator: graft the segment [p1, p2) of tour1 into tour2.

    The segment's first city is the overlapping city; its other cities are
    deleted from tour2 and the segment is inserted right after the overlapping
    city of tour2.
    """
    seg1, _ = split_tour(tour1, p1, p2)
    city, seg1 = seg1[0], seg1[1:]
    removed = set(seg1)
    new = [c for c in tour2 if c not in removed]
    p3 = new.index(city)
    return new[:p3 + 1] + seg1 + new[p3 + 1:]


def _random_pair(rng, n: int) -> tuple[int, int]:
    """MolecularTSP.randomPair: two different positions."""
    p1 = int(rng.integers(n))
    p2 = p1
    while p1 == p2:
        p2 = int(rng.integers(n))
    return p1, p2


def random_exchange(tour, rng) -> list:
    return exchange(tour, *_random_pair(rng, len(tour)))


def random_cut(tour, rng, invert: bool = False) -> list:
    """MolecularTSP.cutOrInvertOperator: redraw until the string differs."""
    if len(tour) < 3:
        return list(tour)
    new = list(tour)
    while new == list(tour):
        p1, p2 = _random_pair(rng, len(tour))
        rest = len(tour) - ((p2 - p1) % len(tour))
        p3 = int(rng.integers(rest - 1)) + 1 if rest > 1 else 0
        new = cut(tour, p1, p2, p3, invert)
    return new


def random_recombine(tour1, tour2, rng) -> list:
    return recombine(tour1, tour2, *_random_pair(rng, len(tour1)))


# ---------------------------------------------------------------------------
# Reactor
# ---------------------------------------------------------------------------
def overlap(tours, n: int) -> float:
    """Paper eq. 3-4: O = sum_ij P_ij^2 / (M^2 N), P the summed upper-half adjacency matrices.

    O = 1 when all M tours share every edge and 1/M when no edge is shared.
    """
    P: Counter = Counter()
    for t in tours:
        for a, b in zip(t, t[1:] + t[:1]):
            P[(a, b) if a < b else (b, a)] += 1
    return sum(v * v for v in P.values()) / (len(tours) ** 2 * n)


def _cities(p, rng) -> list[list[float]]:
    if p.cities:
        pts = p.cities
        if (len(pts) < 3 or not all(isinstance(c, list) and len(c) == 2 for c in pts)
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for c in pts for v in c)):
            raise ValueError(f"cities must be a list of at least 3 [x, y] number pairs, got {pts!r}")
        return [[float(x), float(y)] for x, y in pts]
    if p.layout == "ring":
        return ring_cities(p.N)
    return random_cities(p.N, rng)


def evolve(p, rng):
    """The machine-string soup: a frame per generation (frame 0 is the initial soup)."""
    rates = {"E": p.t_E, "C": p.t_C, "I": p.t_I, "R": p.t_R}
    total = sum(rates.values())
    if total <= 0:
        raise ValueError("at least one machine time scale t_E, t_C, t_I, t_R must be > 0")
    if rates["R"] > 0 and p.M < 2:
        raise ValueError(f"the R-machine needs M >= 2 strings, got M={p.M}; set t_R=0 or raise M")
    cities = _cities(p, rng)
    n = len(cities)
    dist = distance_matrix(cities)
    active = [k for k in MACHINES if rates[k] > 0]
    weights = [rates[k] / total for k in active]

    # initial population of random tours (MolecularTSP.randomTour on the full mesh)
    pop = [[int(c) for c in rng.permutation(n)] for _ in range(p.M)]
    lengths = [tour_length(t, dist) for t in pop]
    start = [tour_id(t) for t in pop]
    seen: dict[str, tuple[tuple[int, ...], float]] = {}
    for t, l in zip(pop, lengths):
        seen.setdefault(tour_id(t), (canonical(t), l))

    ops = math.ceil(p.M / total)
    tally = Tally()
    successes = dict.fromkeys(active, 0)
    machines = {machine_id(k): 1.0 for k in active}
    best = [min(lengths)]
    mean = [sum(lengths) / p.M]

    def frame(g):
        state = dict(machines)
        state.update((s, float(c)) for s, c in Counter(tour_id(t) for t in pop).items())
        return Frame(t=float(g), state=state, fired=tally.flush(),
                     observables={"best_length": best[-1], "mean_length": mean[-1],
                                  "overlap": overlap(pop, n)})

    yield frame(0)
    for generation in range(1, p.generations + 1):
        for _ in range(ops):
            kind = active[int(rng.choice(len(active), p=weights))] if len(active) > 1 else active[0]
            if N_OP[kind] == 1:
                i = int(rng.integers(p.M))
                idx = [i]
                t = pop[i]
                if kind == "E":
                    new = random_exchange(t, rng)
                else:
                    new = random_cut(t, rng, invert=(kind == "I"))
                l = tour_length(new, dist)
                if not l < lengths[i]:
                    continue
                released = [(new, l)]
            else:
                i = int(rng.integers(p.M))
                j = i
                while j == i:
                    j = int(rng.integers(p.M))
                idx = [i, j]
                child = random_recombine(pop[i], pop[j], rng)
                cand = [(pop[i], lengths[i]), (pop[j], lengths[j]), (child, tour_length(child, dist))]
                released = sorted(cand, key=lambda x: x[1])[:2]   # stable: parents win ties
            lhs = [tour_id(pop[k]) for k in idx]
            rhs = [tour_id(t) for t, _ in released]
            if Counter(lhs) == Counter(rhs):
                continue
            for k, (t, l) in zip(idx, released):
                pop[k], lengths[k] = t, l
                seen.setdefault(tour_id(t), (canonical(t), l))
            successes[kind] += 1
            tally.add([machine_id(kind), *lhs], [machine_id(kind), *rhs])
        best.append(min(lengths))
        mean.append(sum(lengths) / p.M)
        yield frame(generation)

    reactions = [Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in tally.reactions()]
    species = [Species(machine_id(k), structure=DESCRIPTION[k]) for k in active]
    species += [Species(sid, structure="(" + ", ".join([f"{l:.6f}", *map(str, t)]) + ")")
                for sid, (t, l) in seen.items()]

    optimum = ring_optimum(n) if not p.cities and p.layout == "ring" else None
    analysis = {
        "generation_size": ops,
        "machine_successes": successes,
        "optimum": optimum,
    }
    if optimum is not None:
        tol = optimum * (1 + 1e-9)
        analysis["generation_optimum_found"] = next((g for g, b in enumerate(best) if b <= tol), None)
        analysis["generation_mean_within_10_percent"] = next(
            (g for g, m in enumerate(mean) if m <= 1.1 * optimum), None)
    k = min(range(p.M), key=lambda i: lengths[i])
    initial = {machine_id(m): 1 for m in active}
    initial.update(Counter(start))
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={s: float(c) for s, c in initial.items()},
        extras={
            "analysis": analysis,
            "cities": cities,
            "tour_length": {sid: l for sid, (_, l) in seen.items()},
            "final_state": dict(Counter(tour_id(t) for t in pop).most_common()),
            "best_tour": {"id": tour_id(pop[k]), "length": lengths[k]},
            "time_scales": rates,
        },
    )
