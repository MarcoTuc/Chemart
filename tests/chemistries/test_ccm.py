"""Chemical Casting Model: book 17.2.2, Kanada & Hirokawa (1994) [439], Kanada (1995, 1996), demo programs."""

import statistics
from itertools import combinations

import pytest

from chemart import generate_network
from chemart.chemistries.ccm import coloring_god, coloring_lod, queens_god, queens_lod, usa_mainland_edges

ID = "ccm"


# ---------------------------------------------------------------------------
# helpers: recompute the instance order degree of a recorded reaction from its species
# ---------------------------------------------------------------------------
def atoms(side):
    """{'q3=5': 1} -> {3: 5}; {'v2=c1': 1} -> {2: 1}."""
    out = {}
    for sid in side:
        name, state = sid.split("=")
        out[int(name[1:])] = int(state.lstrip("c"))
    return out


def iod(reaction, problem):
    """(IOD before, IOD after, number of pairs) over the matched pairs involving a changed atom."""
    before, after = atoms(reaction.reactants), atoms(reaction.products)
    assert before.keys() == after.keys()
    changed = {a for a in before if before[a] != after[a]}
    pairs = [(x, y) for x, y in combinations(sorted(before), 2) if x in changed or y in changed]

    def lod(state, x, y):
        if problem == "n-queens":
            return queens_lod(x, state[x], y, state[y])
        return coloring_lod(True, state[x], state[y])   # catalysts are linked neighbours

    return sum(lod(before, *xy) for xy in pairs), sum(lod(after, *xy) for xy in pairs), len(pairs)


def valid_queens(cols):
    return queens_god(cols) == len(cols) * (len(cols) - 1) // 2


def assert_observed_run(net, problem, strict_monotone=False, frustration=False):
    a = net.extras["analysis"]
    assert net.status == "observed"
    assert sum(r.count for r in net.reactions) == a["reactions"]
    assert len(a["god"]) == a["reactions"] + 1
    assert a["god"][-1] == a["god_final"]
    ids = {s.id for s in net.species}
    assert set(net.initial_state) <= ids
    for r in net.reactions:
        assert r.count >= 1 and r.rate is None
        o, o2, pairs = iod(r, problem)
        if pairs:
            assert o < pairs                   # an instance with nothing to repair never reacts
        if not frustration:
            assert o2 >= o                      # Kanada: the IOD is not decreased
            if strict_monotone:
                assert o2 > o                   # the book's wording: LHS sum smaller than RHS sum
    if not frustration:
        assert a["uphill_reactions"] == 0


# ---------------------------------------------------------------------------
# rules, LODs and instances
# ---------------------------------------------------------------------------
def test_lods_and_the_usa_map():
    """HICSS-27 fig. 4 (queens LOD), FUZZ-IEEE'95 sec. 3 (colouring LOD, USA map of 48 states, mean degree 4.42)."""
    assert queens_lod(0, 0, 3, 3) == queens_lod(0, 3, 3, 0) == 0      # diagonals
    assert queens_lod(0, 0, 1, 2) == 1
    assert queens_god([0, 1, 2, 3, 4, 5, 6, 7]) == 0                  # all on a diagonal: the GOD minimum
    assert valid_queens([0, 4, 7, 5, 2, 6, 1, 3])                      # a known 8-queens solution
    assert coloring_lod(False, 1, 1) == 1 and coloring_lod(True, 1, 2) == 1 and coloring_lod(True, 2, 2) == 0
    edges = usa_mainland_edges()
    assert len({v for e in edges for v in e}) == 48 and len(set(edges)) == 106
    assert round(2 * len(edges) / 48, 2) == 4.42
    assert coloring_god([0] * 48, edges) == 0


def test_six_queens_local_maximum_of_fig_12():
    """HICSS-27 fig. 12 and sec. 5.4: a 6-queens state with GOD 14 (max 15) that no swap improves,
    escaped by the one-catalyst rule matching the queens of columns 6, 5 and 4 to Queen1, 2, 3."""
    rows = [6, 1, 4, 2, 5, 3]                          # column (1-based) of the queen in each row, read off fig. 12
    s = [c - 1 for c in rows]
    assert queens_god(s) == 14
    swaps = []
    for i, j in combinations(range(6), 2):
        t = list(s)
        t[i], t[j] = t[j], t[i]
        swaps.append(queens_god(t))
    assert max(swaps) <= 14                            # GOD(S) >= GOD(S') for every swap: a local maximum
    q1, q2, q3 = rows.index(6), rows.index(5), rows.index(4)
    before = queens_lod(q2, s[q2], q3, s[q3]) + queens_lod(q1, s[q1], q2, s[q2]) + queens_lod(q1, s[q1], q3, s[q3])
    after = queens_lod(q2, s[q3], q3, s[q2]) + queens_lod(q1, s[q1], q2, s[q3]) + queens_lod(q1, s[q1], q3, s[q2])
    assert (before, after) == (2, 3)                   # the local IOD increases, so the instance reacts ...
    t = list(s)
    t[q2], t[q3] = t[q3], t[q2]
    assert queens_god(t) == 13                         # ... while the global order degree drops: a conflicting system


# ---------------------------------------------------------------------------
# runs
# ---------------------------------------------------------------------------
def test_default_run_solves_eight_queens():
    """HICSS-27 sec. 4: the single-catalyst swap rule from queens on the diagonal reaches a solution and stops."""
    net = generate_network(ID, seed=1)
    a = net.extras["analysis"]
    assert a["solved"] and a["terminated"]
    assert a["god_initial"] == 0 and a["god_max"] == 28
    assert valid_queens(net.extras["final_assignment"])
    assert net.initial_state == {f"q{i}={i}": 1.0 for i in range(8)}
    assert all(len(r.catalysts) == 1 and len(r.reactants) == 3 for r in net.reactions)
    assert_observed_run(net, "n-queens", frustration=True)


def test_eight_queens_never_fail():
    """HICSS-27 sec. 4.2: 'the problems never fail to be solved' (original CCM, random and diagonal layouts)."""
    for seed in range(4):
        for initial in ("ordered", "random"):
            net = generate_network(ID, seed=seed, frustration=False, initial=initial)
            assert net.extras["analysis"]["solved"], (seed, initial)
            assert valid_queens(net.extras["final_assignment"])
            assert_observed_run(net, "n-queens")


def test_conflicting_system_god_is_not_monotone():
    """HICSS-27 fig. 8 and sec. 5.2: N queens is a conflicting system; GOD falls although no reaction lowers its IOD."""
    drops = 0
    for seed in range(4):
        a = generate_network(ID, seed=seed, frustration=False).extras["analysis"]
        assert a["uphill_reactions"] == 0
        drops += sum(y < x for x, y in zip(a["god"], a["god"][1:]))
    assert drops > 0


def test_strict_acceptance_is_the_books_criterion():
    """Book 17.2.2: a rule applies if the LHS LOD sum is smaller than the RHS sum (HICSS-27 footnote 3 variant)."""
    net = generate_network(ID, seed=2, frustration=False, acceptance="increasing")
    assert net.extras["analysis"]["solved"]
    assert_observed_run(net, "n-queens", strict_monotone=True)


def test_rule_without_catalyst_never_stops():
    """HICSS-27 sec. 5.3: without the catalyst 'the system does not stop even when a solution is found'."""
    a = generate_network(ID, seed=0, rule="no-catalyst", max_tests=3000).extras["analysis"]
    assert not a["terminated"]
    assert a["tests"] == a["reactions"] == 3000
    net = generate_network(ID, seed=0, problem="graph-coloring", rule="no-catalyst", max_tests=2000)
    assert not net.extras["analysis"]["terminated"] and net.extras["analysis"]["reactions"] == 2000
    assert not any(r.catalysts for r in net.reactions)


def test_more_catalysts_fewer_reactions():
    """HICSS-27 figs. 9-10: the number of reactions to a solution falls as catalysts are added (8 queens)."""
    def mean_reactions(rule):
        runs = [generate_network(ID, seed=s, rule=rule, frustration=False).extras["analysis"] for s in range(10)]
        assert all(r["solved"] for r in runs)
        return statistics.mean(r["reactions"] for r in runs)
    one, three = mean_reactions("single-catalyst"), mean_reactions("triple-catalyst")
    assert three < 0.6 * one


def test_usa_map_is_colored_in_every_run():
    """FUZZ-IEEE'95 sec. 3: four colours from a one-colour map; 'a correct solution was found in every run'."""
    edges = usa_mainland_edges()
    for rule in ("single-catalyst", "variable-catalyst"):
        for seed in range(3):
            net = generate_network(ID, seed=seed, problem="graph-coloring", rule=rule)
            colors = net.extras["final_assignment"]
            assert net.extras["analysis"]["solved"] and all(colors[a] != colors[b] for a, b in edges)
            assert set(colors) <= {0, 1, 2, 3}
            assert net.initial_state == {f"v{i}=c0": 1.0 for i in range(48)}
            assert_observed_run(net, "graph-coloring", frustration=True)
            if rule == "variable-catalyst":
                v_links = {v: {b for a, b in edges if a == v} | {a for a, b in edges if b == v} for v in range(48)}
                for r in net.reactions:
                    (changed,) = [int(s[1:].split("=")[0]) for s in r.reactants if s not in r.products]
                    assert {int(s[1:].split("=")[0]) for s in r.catalysts} == v_links[changed]


def test_fam_rescues_the_variable_catalyst_rule():
    """FUZZ-IEEE'95 sec. 3, Kanada 1996 sec. 3.2-4 and the demo pages: without FAM the variable-catalyst rule
    sometimes stops in or flips between non-solutions; with FAM (f0 = 1e-5, c = 2) it almost always solves.
    With FAM the USA map needs about as many reactions as FUZZ-IEEE'95 table 2 reports (112 on average)."""
    kw = dict(problem="graph-coloring", rule="variable-catalyst", max_tests=30000)
    without = [generate_network(ID, seed=s, frustration=False, **kw).extras["analysis"] for s in range(10)]
    with_fam = [generate_network(ID, seed=s, **kw).extras["analysis"] for s in range(10)]
    assert sum(not r["solved"] for r in without) >= 2
    assert all(r["solved"] for r in with_fam)
    assert 60 < statistics.mean(r["reactions"] for r in with_fam) < 200
    queens = [generate_network(ID, seed=s, rule="variable-catalyst").extras["analysis"] for s in range(10)]
    assert all(r["solved"] for r in queens)
    assert sum(r["uphill_reactions"] for r in queens) > 0      # frustration lets the IOD fall


def test_random_and_custom_graphs():
    net = generate_network(ID, seed=3, problem="graph-coloring", graph="random", V=20, edge_probability=0.2,
                           colors=4, initial="random")
    inst = net.extras["instance"]
    assert inst["vertices"] == 20 and all(0 <= u < v < 20 for u, v in inst["edges"])
    assert net.to_dict() == generate_network(ID, seed=3, problem="graph-coloring", graph="random", V=20,
                                             edge_probability=0.2, colors=4, initial="random").to_dict()
    triangle = [[0, 1], [1, 2], [2, 0], [1, 0]]
    ok = generate_network(ID, seed=0, problem="graph-coloring", graph="custom", edges=triangle, colors=3)
    assert ok.extras["instance"]["edges"] == [[0, 1], [1, 2], [0, 2]]
    assert ok.extras["analysis"]["solved"] and sorted(ok.extras["final_assignment"]) == [0, 1, 2]
    impossible = generate_network(ID, seed=0, problem="graph-coloring", graph="custom", edges=triangle, colors=2,
                                  frustration=False, max_tests=5000)
    assert not impossible.extras["analysis"]["solved"]


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(N=4, rule="triple-catalyst"), "more patterns than queens"),
        (dict(f0=0.0), "f0 > 0"),
        (dict(problem="graph-coloring", graph="custom"), "needs edges"),
        (dict(problem="graph-coloring", graph="custom", edges=[[1, 1]]), "self-loop"),
    ],
)
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
