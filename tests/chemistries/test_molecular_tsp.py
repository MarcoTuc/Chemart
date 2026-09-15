"""Molecular TSP: book 17.2.1, Banzhaf (1990) [62], PyCellChemistry MolecularTSP.py."""

import math
from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.molecular_tsp import (
    canonical, cut, distance_matrix, exchange, overlap, random_cut, random_exchange,
    random_recombine, recombine, ring_cities, ring_optimum, tour_id, tour_length,
)

ID = "molecular-tsp"
FIG2A = [1, 2, 3, 4, 5, 6, 7]


def same_cycle(a, b):
    return canonical(list(a)) == canonical(list(b))


def test_operators_reproduce_paper_fig2():
    """Paper fig. 2 (b-e) on the sample tour 1..7, and the examples in MolecularTSP.py."""
    # b: exchange cities 4 and 6 -> 1 2 3 6 5 4 7
    assert exchange(FIG2A, 3, 5) == [1, 2, 3, 6, 5, 4, 7]
    # c: cut the part between 2 and 3 and move it behind 5 -> 1 4 5 2 3 6 7
    c = cut(FIG2A, 1, 3, 2)
    assert same_cycle(c, [1, 4, 5, 2, 3, 6, 7])
    # d: cut-inverse, the part 2 3 inserted inverted behind 6 -> 1 4 5 6 3 2 7
    assert same_cycle(cut(FIG2A, 1, 3, 3, invert=True), [1, 4, 5, 6, 3, 2, 7])
    # MolecularTSP.py invertMachine docstring: 1 2 3 4 5 6 7 ==> 1 4 5 3 2 6 7
    assert same_cycle(cut(FIG2A, 1, 3, 2, invert=True), [1, 4, 5, 3, 2, 6, 7])
    # e: the segment 1 4 5 2 grafted at the overlapping city 1 -> 1 4 5 2 7 6 3
    assert c == [4, 5, 2, 3, 6, 7, 1]                 # the cut string, read from p2
    child = recombine(c, [1, 7, 6, 5, 4, 3, 2], 6, 3)  # circular segment 1 4 5 2 of c
    assert child == [1, 4, 5, 2, 7, 6, 3]
    # MolecularTSP.py recombinationOperator docstring: [3 6 1 4 5 2 7] grafted into [1..7]
    assert recombine([3, 6, 1, 4, 5, 2, 7], FIG2A, 2, 6) == [1, 4, 5, 2, 3, 6, 7]


def test_tour_length_is_book_eq_17_1():
    """Book 17.2.1: s = (l, 3, 1, 4, 2) is a closed tour; l sums euclidean distances."""
    square = {1: (0, 0), 2: (1, 0), 3: (1, 1), 4: (0, 1)}
    cities = [square[i] for i in (1, 2, 3, 4)]
    d = distance_matrix(cities)
    as_index = [3 - 1, 1 - 1, 4 - 1, 2 - 1]          # 3 -> 1 -> 4 -> 2 -> 3
    assert tour_length(as_index, d) == pytest.approx(2 + 2 * math.sqrt(2))
    assert tour_length([0, 1, 2, 3], d) == pytest.approx(4.0)
    # rotations and reversal are the same data string species
    assert tour_id([2, 3, 0, 1]) == tour_id([1, 0, 3, 2]) == tour_id([0, 1, 2, 3]) == "t0.1.2.3"
    # the ring optimum is the regular polygon
    ring = distance_matrix(ring_cities(12))
    assert tour_length(list(range(12)), ring) == pytest.approx(ring_optimum(12))


def test_random_operators_keep_valid_tours():
    """Book: 'these operations are performed with safeguards to ensure valid tours'."""
    rng = np.random.default_rng(0)
    for _ in range(500):
        n = int(rng.integers(3, 15))
        t1, t2 = [int(c) for c in rng.permutation(n)], [int(c) for c in rng.permutation(n)]
        e = random_exchange(t1, rng)
        assert sorted(e) == sorted(t1) and sum(a != b for a, b in zip(e, t1)) == 2
        for invert in (False, True):
            c = random_cut(t1, rng, invert)
            assert sorted(c) == sorted(t1) and c != t1
        r = random_recombine(t1, t2, rng)
        assert sorted(r) == sorted(t1)


def test_overlap_is_paper_eq_4():
    pentagon, pentagram = [0, 1, 2, 3, 4], [0, 2, 4, 1, 3]     # edge-disjoint on K5
    assert overlap([pentagon] * 9, 5) == pytest.approx(1.0)
    assert overlap([pentagon, pentagram], 5) == pytest.approx(1 / 2)


def assert_machine_reaction(net, r):
    """Machine + n_op strings -> machine + the n_op best (book 17.2.1)."""
    length = net.extras["tour_length"]
    (machine,) = r.catalysts.keys() & {"E-machine", "C-machine", "I-machine", "R-machine"}
    ins = [s for s, n in r.reactants.items() if s != machine for _ in range(n)]
    outs = [s for s, n in r.products.items() if s != machine for _ in range(n)]
    n_op = 2 if machine == "R-machine" else 1
    assert len(ins) == len(outs) == n_op
    if n_op == 1:
        assert length[outs[0]] < length[ins[0]]
    else:
        (child,) = (Counter(outs) - Counter(ins)).elements()
        (lost,) = (Counter(ins) - Counter(outs)).elements()
        assert length[child] < length[lost]
        assert sorted(length[s] for s in outs) == sorted(length[s] for s in [*ins, child])[:2]


def test_default_ring_run_reaches_the_polygon():
    """Book fig. 17.2 / paper fig. 3: random tours on a ring converge to the ring in 1000 generations."""
    net = generate_network(ID, seed=3)
    a = net.extras["analysis"]
    assert net.status == "observed"
    assert a["generation_size"] == math.ceil(9 / 3.01) == 3
    assert len(a["best_length"]) == len(a["mean_length"]) == 1001
    assert a["best_length"][0] > 1.3 * a["optimum"]
    assert a["best_length"][-1] == pytest.approx(a["optimum"])
    assert net.extras["best_tour"]["id"] == "t" + ".".join(map(str, range(10)))
    # local selection: neither the best nor the mean length ever increases
    for series in (a["best_length"], a["mean_length"]):
        assert all(b <= x + 1e-9 for x, b in zip(series, series[1:]))
    assert a["overlap"][-1] > a["overlap"][0]
    for r in net.reactions:
        assert r.count >= 1
        assert_machine_reaction(net, r)
    machines = {s: n for s, n in net.initial_state.items() if s.endswith("-machine")}
    assert machines == {"E-machine": 1, "C-machine": 1, "I-machine": 1, "R-machine": 1}
    assert sum(net.initial_state.values()) - 4 == sum(net.extras["final_state"].values()) == 9
    assert sum(a["machine_successes"].values()) == sum(r.count for r in net.reactions)


def gens_to_quality(seeds, **kw):
    return [generate_network(ID, seed=s, N=20, **kw).extras["analysis"]["generation_mean_within_10_percent"]
            for s in seeds]


def test_recombination_accelerates_the_search():
    """Paper table 1a (ring, M = 9, criterion: mean 10% above optimum): E + R is far faster than E alone."""
    alone = gens_to_quality(range(3), t_C=0, t_I=0, t_R=0, generations=1500)
    with_r = gens_to_quality(range(3), t_C=0, t_I=0, t_R=0.01, generations=1500)
    assert alone == [None, None, None]
    assert all(g is not None for g in with_r)


def test_higher_recombination_frequency_is_faster():
    """Paper table 1b: t_R = 1 converges in far fewer generations than t_R = 1/1000."""
    rare = gens_to_quality(range(3), t_R=0.001, generations=1000)
    often = gens_to_quality(range(3), t_R=1.0, generations=1000)
    assert rare == [None, None, None]
    assert all(g is not None for g in often)


def test_recombination_collapses_variance_on_random_cities():
    """Paper table 2b / conclusion: frequent recombination drives the population overlap (eq. 4) up."""
    def final(t_R, s):
        a = generate_network(ID, seed=s, N=20, layout="random", t_R=t_R, generations=400).extras["analysis"]
        return a["overlap"][-1], a["best_length"][-1]
    rare = [final(0.001, s) for s in range(3)]
    often = [final(1.0, s) for s in range(3)]
    assert all(o > 0.8 for o, _ in often)
    assert all(o < 0.5 for o, _ in rare)
    # and it helps escape local minima (paper table 2a)
    assert sum(b for _, b in often) < sum(b for _, b in rare)


def test_layouts_and_explicit_cities():
    net = generate_network(ID, seed=2, layout="random", N=12, generations=5)
    pts = net.extras["cities"]
    assert len(pts) == 12 and all(float(v).is_integer() and 0 <= v < 24 for p in pts for v in p)
    assert min(math.dist(p, q) for i, p in enumerate(pts) for q in pts[i + 1:]) >= 2
    assert net.extras["analysis"]["optimum"] is None
    square = generate_network(ID, seed=0, cities=[[0, 0], [0, 1], [1, 1], [1, 0], [0.5, 2]], generations=200)
    assert square.extras["best_tour"]["length"] == pytest.approx(3 + 2 * math.hypot(0.5, 1))
    only_e = generate_network(ID, seed=0, t_C=0, t_I=0, t_R=0, M=1, generations=50)
    assert {r_c for r in only_e.reactions for r_c in r.catalysts} == {"E-machine"}


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(t_E=0, t_C=0, t_I=0, t_R=0), "at least one machine"),
        (dict(M=1), "M >= 2"),
        (dict(cities=[[0, 0], [1, 1]]), "at least 3"),
        (dict(cities=[[0, 0], [1, 1], ["a", 2]]), "at least 3"),
    ],
)
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
