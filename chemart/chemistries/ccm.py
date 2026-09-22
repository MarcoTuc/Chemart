"""Chemical Casting Model (CCM): Kanada (1992-1996), book 17.2.2.

Catalog id: ccm.

The working memory is ONE point of the search space: a set of atoms, each with
a type and a state, possibly linked. A reaction rule LHS -> RHS rewrites the
states of a few matched atoms; it is applied only if the instance order degree
(IOD), the sum of the local order degrees (LODs) of the matched atoms, does not
decrease. Atoms that a rule matches but does not change are catalysts: more
catalysts make the rule less local (Kanada & Hirokawa 1994, HICSS-27).

Two problems, with the rules and LODs of Kanada's papers and his own demo
programs (Queens_Sort and USAmap_color, 1996; JavaScript ports on kanadas.com):

n-queens (HICSS-27 sec. 4; Queens_Sort)
    atoms q<row> with state column; LOD o(x, y) = 0 if queens x and y attack
    each other (same column or same diagonal), else 1.
    swap rule with k catalysts (k = 0..3):
        q<a>=ca + q<b>=cb + k catalysts -> q<a>=cb + q<b>=ca + the same catalysts
    variable-catalyst moving rule: q<a>=c -> q<a>=c' with all other queens as
    catalysts.

graph-coloring (FUZZ-IEEE'95 sec. 3; Kanada 1996 sec. 3.2; USAmap_color)
    atoms v<i> with state colour c<k>, linked by the graph's edges; LOD
    o(v1, v2) = 1 if not connected or of different colours, else 0.
    recolouring rule with k linked neighbours of the vertex as catalysts
    (k = 0..3), or all its neighbours (variable-catalyst rule):
        v<i>=c1 + neighbours -> v<i>=c2 + the same neighbours

Reaction test (the demos' `reaction`/`swap`): let O and O' be the IOD over the
matched pairs involving a changed atom before and after, M its maximum and F
the summed frustration of the changed atoms. A rule without catalysts always
reacts. Otherwise the instance does not react if O = M (nothing to repair);
it reacts if O' >= O - F (acceptance non-decreasing) or O' > O - F
(increasing). Frustration accumulation method (FAM, Kanada 1995-1996): a
tested instance that does not react multiplies the frustration of each changed
atom by c; a reaction resets it to f0.

Scheduling is the mathematical-random strategy: every test draws the atoms at
random. The run ends with the demos' termination test (a run of failed tests
of adaptive length) or after max_tests tests. The chemistry is a gas with one
face, ``evolve``: a frame after every accepted reaction, timed in tests, with
the global order degree (GOD) as its observable. The observed network holds
the distinct reactions that fired, with counts.
"""

from __future__ import annotations

from chemart.helpers import params as check
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame

RULES = ("no-catalyst", "single-catalyst", "double-catalyst", "triple-catalyst", "variable-catalyst")
FIXED_CATALYSTS = {"no-catalyst": 0, "single-catalyst": 1, "double-catalyst": 2, "triple-catalyst": 3}

# USA mainland map (48 states) as a graph: USAmap_color.init_wm (Kanada 1996,
# map data by Y. Sato). 106 borders, mean degree 4.42 (FUZZ-IEEE'95 sec. 3).
_USA_SOURCES = [0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11,
                12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 19, 20,
                20, 20, 21, 22, 22, 23, 23, 24, 24, 24, 24, 24, 24, 25, 26, 26, 27, 27, 28, 28, 28, 29, 29, 30,
                30, 31, 31, 31, 32, 32, 33, 33, 34, 34, 37, 37, 37, 37, 37, 38, 38, 38, 39, 40, 41, 41, 42, 42,
                42, 43, 46]
_USA_DESTINATIONS = [1, 4, 2, 3, 4, 3, 6, 4, 5, 6, 5, 7, 8, 6, 8, 9, 10, 8, 11, 12, 9, 12, 13, 10, 13, 14, 15, 15,
                     16, 12, 17, 13, 17, 18, 14, 18, 19, 15, 18, 19, 16, 19, 20, 20, 21, 18, 22, 19, 22, 23, 20,
                     23, 24, 28, 21, 24, 25, 25, 23, 26, 27, 28, 25, 28, 29, 32, 33, 34, 29, 27, 30, 28, 30, 30,
                     31, 32, 34, 35, 31, 38, 32, 38, 39, 33, 39, 34, 36, 35, 36, 38, 41, 42, 43, 44, 39, 40, 44,
                     40, 44, 42, 46, 43, 45, 46, 45, 47]
USA_STATES = 48


def usa_mainland_edges() -> list[tuple[int, int]]:
    return list(zip(_USA_SOURCES, _USA_DESTINATIONS))


# ---------------------------------------------------------------------------
# LODs and global order degree (GOD, HICSS-27 sec. 5.1: the sum of all LODs)
# ---------------------------------------------------------------------------
def queens_lod(row1: int, col1: int, row2: int, col2: int) -> int:
    """Queens_Sort.order: 0 if the queens can take each other, else 1 (HICSS-27 fig. 4)."""
    if row1 == row2 or col1 == col2 or row1 - row2 == col1 - col2 or row1 - row2 == col2 - col1:
        return 0
    return 1


def queens_god(columns) -> int:
    """Number of non-attacking pairs; the maximum N(N-1)/2 is reached exactly at solutions."""
    n = len(columns)
    return sum(queens_lod(i, columns[i], j, columns[j]) for i in range(n) for j in range(i))


def coloring_lod(connected: bool, color1: int, color2: int) -> int:
    """o(v1, v2) = 1 if not connected(v1, v2) or v1.color != v2.color, else 0."""
    return 1 if (not connected or color1 != color2) else 0


def coloring_god(colors, edges) -> int:
    """Satisfied edges (unconnected pairs always have LOD 1 and are not counted)."""
    return sum(colors[a] != colors[b] for a, b in edges)


def queen_id(row: int, col: int) -> str:
    return f"q{row}={col}"


def vertex_id(v: int, color: int) -> str:
    return f"v{v}=c{color}"


# ---------------------------------------------------------------------------
# Problem instances
# ---------------------------------------------------------------------------
def _graph(p, rng) -> tuple[int, list[tuple[int, int]]]:
    if p.graph == "usa-mainland":
        return USA_STATES, usa_mainland_edges()
    if p.graph == "random":
        n = p.V
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < p.edge_probability]
        return n, edges
    if not p.edges:
        raise ValueError("graph='custom' needs edges, a list of [u, v] vertex pairs")
    raw = check.edges("edges", p.edges)
    if any(u < 0 or v < 0 for u, v in raw):
        raise ValueError(f"edges: vertex ids must be >= 0, got {p.edges!r}")
    seen, edges = set(), []
    for u, v in raw:
        key = (min(u, v), max(u, v))
        if key not in seen:
            seen.add(key)
            edges.append(key)
    return max(max(e) for e in edges) + 1, edges


def _distinct(rng, n: int, k: int) -> list[int]:
    """k distinct indices in range(n), uniformly (rejection, as the demos do)."""
    out: list[int] = []
    while len(out) < k:
        i = int(rng.integers(n))
        if i not in out:
            out.append(i)
    return out


def _other(rng, n: int, old: int) -> int:
    """A uniform value in range(n) different from old (USAmap_color.select_color)."""
    v = int(rng.integers(n - 1))
    return v + 1 if v >= old else v


# ---------------------------------------------------------------------------
# Reactor
# ---------------------------------------------------------------------------
class _Run:
    def __init__(self, p, n_atoms: int, patience: int):
        self.p = p
        self.f0 = p.f0 if p.frustration else 0.0
        self.frustration = [self.f0] * n_atoms
        self.tests = 0
        self.attempts = 0
        self.failed = 0
        self.reactions = 0
        self.uphill = 0
        self.tally = Tally()
        self.shown = 0                          # tests at the last frame
        self.threshold = patience
        self.threshold_20 = 20 * patience
        self.terminated = False

    def decide(self, changed, before: int, after: int, pairs: int, unconditional: bool) -> bool:
        """The demos' reaction test and FAM update; True if the instance reacts."""
        if unconditional:
            return True
        if before >= pairs:
            self.failed += 1
            return False
        F = sum(self.frustration[a] for a in changed)
        ok = after >= before - F if self.p.acceptance == "non-decreasing" else after > before - F
        if not ok:
            self.failed += 1
            if self.p.frustration:
                for a in changed:
                    self.frustration[a] *= self.p.c
            return False
        return True

    def react(self, changed, lhs: list[str], rhs: list[str], before: int, after: int) -> None:
        self.reactions += 1
        self.failed = 0
        if after < before:
            self.uphill += 1
        for a in changed:
            self.frustration[a] = self.f0
        self.tally.add(lhs, rhs)

    def frame(self, state: dict[str, float], god: int) -> Frame:
        """A frame at the current test: the working memory, the reactions since the last frame, GOD."""
        self.shown = self.tests
        return Frame(t=float(self.tests), state=state, fired=self.tally.flush(), observables={"god": god})

    def done(self) -> bool:
        """Queens_Sort/USAmap_color loopTick termination test, plus the max_tests budget."""
        if self.failed >= self.threshold:
            if self.tests < self.threshold_20:
                self.terminated = True
                return True
            self.threshold_20 *= 2
            self.threshold = self.threshold_20 // 20
        return self.attempts >= self.p.max_tests


def _queens(p, rng):
    n = p.N
    rule = p.rule
    k = FIXED_CATALYSTS.get(rule)
    patterns = n if rule == "variable-catalyst" else 2 + k
    if patterns > n:
        raise ValueError(f"rule {rule!r} matches {patterns} queens but N={n}; a rule with more patterns "
                         f"than queens cannot produce a solution (HICSS-27 sec. 5.3): lower the catalysts or raise N")
    cols = [int(c) for c in rng.permutation(n)] if p.initial == "random" else list(range(n))
    initial = list(cols)
    run = _Run(p, n, 1000 + n * n)
    god = queens_god(cols)
    trace = [god]

    def state():
        return {queen_id(r, c): 1.0 for r, c in enumerate(cols)}

    yield run.frame(state(), god)

    while not run.done():
        run.attempts += 1
        run.tests += 1
        if rule == "variable-catalyst":
            a = int(rng.integers(n))
            new = _other(rng, n, cols[a])
            before = sum(queens_lod(a, cols[a], m, cols[m]) for m in range(n) if m != a)
            after = sum(queens_lod(a, new, m, cols[m]) for m in range(n) if m != a)
            if not run.decide([a], before, after, n - 1, False):
                continue
            cats = [m for m in range(n) if m != a]
            lhs = [queen_id(a, cols[a])] + [queen_id(m, cols[m]) for m in cats]
            rhs = [queen_id(a, new)] + [queen_id(m, cols[m]) for m in cats]
            god += after - before
            cols[a] = new
        else:
            a, b, *cats = _distinct(rng, n, 2 + k)
            ca, cb = cols[a], cols[b]
            before = queens_lod(a, ca, b, cb) + sum(queens_lod(a, ca, m, cols[m]) + queens_lod(b, cb, m, cols[m])
                                                    for m in cats)
            after = queens_lod(a, cb, b, ca) + sum(queens_lod(a, cb, m, cols[m]) + queens_lod(b, ca, m, cols[m])
                                                   for m in cats)
            if not run.decide([a, b], before, after, 1 + 2 * k, k == 0):
                continue
            lhs = [queen_id(a, ca), queen_id(b, cb)] + [queen_id(m, cols[m]) for m in cats]
            rhs = [queen_id(a, cb), queen_id(b, ca)] + [queen_id(m, cols[m]) for m in cats]
            god += (sum(queens_lod(a, cb, m, cols[m]) + queens_lod(b, ca, m, cols[m])
                        - queens_lod(a, ca, m, cols[m]) - queens_lod(b, cb, m, cols[m])
                        for m in range(n) if m != a and m != b)
                    + queens_lod(a, cb, b, ca) - queens_lod(a, ca, b, cb))
            cols[a], cols[b] = cb, ca
        run.react([a] if rule == "variable-catalyst" else [a, b], lhs, rhs, before, after)
        trace.append(god)
        yield run.frame(state(), god)
    if run.tests > run.shown:
        yield run.frame(state(), god)

    start = [queen_id(r, c) for r, c in enumerate(initial)]
    return run, start, cols, trace, n * (n - 1) // 2, initial, {"N": n}, species_of_queen


def species_of_queen(sid: str) -> str:
    row, col = sid[1:].split("=")
    return f"(queen, row={row}, column={col})"


def species_of_vertex(sid: str) -> str:
    v, color = sid[1:].split("=c")
    return f"(vertex, id={v}, color={color})"


def _coloring(p, rng):
    n, edges = _graph(p, rng)
    q = p.colors
    neighbours: list[list[int]] = [[] for _ in range(n)]
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    colors = [int(c) for c in rng.integers(q, size=n)] if p.initial == "random" else [0] * n
    initial = list(colors)
    rule = p.rule
    k = FIXED_CATALYSTS.get(rule)
    run = _Run(p, n, 1000)
    god = coloring_god(colors, edges)
    trace = [god]

    def state():
        return {vertex_id(i, c): 1.0 for i, c in enumerate(colors)}

    yield run.frame(state(), god)

    while not run.done():
        run.attempts += 1
        v = int(rng.integers(n))
        nb = neighbours[v]
        if rule == "variable-catalyst":
            cats = list(nb)
        elif k == 0:
            cats = []
        else:
            if len(nb) < k:                         # USAmap_color.change_color1/2
                run.failed += 1
                continue
            cats = [nb[i] for i in _distinct(rng, len(nb), k)]
        run.tests += 1
        old = colors[v]
        new = _other(rng, q, old)
        before = sum(colors[m] != old for m in cats)
        after = sum(colors[m] != new for m in cats)
        if not run.decide([v], before, after, len(cats), k == 0):
            continue
        lhs = [vertex_id(v, old)] + [vertex_id(m, colors[m]) for m in cats]
        rhs = [vertex_id(v, new)] + [vertex_id(m, colors[m]) for m in cats]
        god += sum(colors[m] != new for m in nb) - sum(colors[m] != old for m in nb)
        colors[v] = new
        run.react([v], lhs, rhs, before, after)
        trace.append(god)
        yield run.frame(state(), god)
    if run.tests > run.shown:
        yield run.frame(state(), god)

    instance = {"vertices": n, "edges": [list(e) for e in edges], "colors": q,
                "mean_degree": 2 * len(edges) / n if n else 0.0}
    start = [vertex_id(i, c) for i, c in enumerate(initial)]
    return run, start, colors, trace, len(edges), initial, instance, species_of_vertex


def evolve(p, rng):
    """Kanada's reactor: a frame after every accepted reaction, timed in tests, plus one
    at the last test if the run ended on failed tests; the observable is GOD."""
    if p.rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}, got {p.rule!r}")
    if p.frustration and p.f0 <= 0:
        raise ValueError("frustration=True needs f0 > 0 (a zero frustration never grows); set frustration=False instead")
    solve = _queens if p.problem == "n-queens" else _coloring
    run, start, final, trace, god_max, initial, instance, describe = yield from solve(p, rng)

    order: dict[str, None] = dict.fromkeys(start)
    reactions = []
    for lhs, rhs, count in run.tally.reactions():
        order.update(dict.fromkeys(rhs))
        reactions.append(Reaction.of(lhs, rhs, count=count))
    species = [Species(sid, structure=describe(sid)) for sid in order]

    god = trace[-1]
    first = next((i for i, g in enumerate(trace) if g == god_max), None)
    analysis = {
        "tests": run.tests,
        "reactions": run.reactions,
        "uphill_reactions": run.uphill,
        "terminated": run.terminated,
        "god_initial": trace[0],
        "god_final": god,
        "god_max": god_max,
        "mod_final": god / god_max if god_max else 1.0,
        "solved": god == god_max,
        "first_solution_reaction": first,
        "mean_frustration_final": sum(run.frustration) / len(run.frustration),
    }
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={s: 1.0 for s in start},
        extras={
            "analysis": analysis,
            "instance": instance,
            "initial_assignment": initial,
            "final_assignment": final,
        },
    )
