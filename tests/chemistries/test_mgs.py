"""MGS transformations reproduce the published programs.

Book: Banzhaf & Yamamoto 2015, section 9.6 (sort, EcoRI).
2002: Giavitto & Michel, Data structure as topological spaces, UMC 2002, LNCS 2509
(restriction enzymes sec. 3.1, Eden sec. 3.2, Turn fig. 2).
2003: Cohen, Typing rule-based transformations over topological collections,
ENTCS 86(2) (bubble sort and Eratosthenes' sieve sec. 2.1, bead-sort sec. 5.2 fig. 4).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.mgs import PROGRAMS, Machine, bead_rows

ID = "mgs"


def primes_upto(n):
    return [k for k in range(2, n + 1) if all(k % d for d in range(2, k))]


def inert(net, present):
    """No reaction can fire from this set of species."""
    return not any(set(r.reactants) <= set(present) for r in net.reactions)


# --- sort (book 9.6; Cohen 2003 sec. 2.1) --------------------------------------------------
def test_book_sort_network_and_fixpoint():
    net = generate_network(ID)
    assert net.status == "complete"
    assert net.extras["rules"] == ["r1: x, y / x > y => y, x"]
    run = net.extras["run"]
    assert run["fixpoint"] and run["final"] == [1, 2, 3, 4, 5]
    for rx in net.reactions:                              # position-specific swaps x@i + y@i+1
        (i, x), (j, y) = sorted((int(s.split("@")[1]), int(s.split("@")[0])) for s in rx.reactants)
        assert j == i + 1 and x > y
        assert rx.products == {f"{y}@{i}": 1, f"{x}@{j}": 1}
    assert inert(net, [f"{v}@{i}" for i, v in enumerate(run["final"])])
    assert not inert(net, [f"{v}@{i}" for i, v in enumerate([4, 2, 5, 1, 3])])
    assert len(net.extras["conservation"]) == 5
    assert net.extras["space"]["collection"] == "seq"


@pytest.mark.parametrize("strategy", ["maximal-parallel", "asynchronous"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_sort_reaches_the_sorted_sequence(strategy, seed):
    values = [9, 3, 7, 3, 1, 8, 2]
    net = generate_network(ID, initial=values, strategy=strategy, seed=seed)
    assert net.extras["run"]["final"] == sorted(values)


def test_cohen_bubble_sort_notation_gives_the_same_network():
    book = generate_network(ID)
    cohen = generate_network(ID, program="custom", collection="seq", initial=[4, 2, 5, 1, 3],
                             rules=["x, y/(y<x) => y :: x :: empty_seq", "x => [x]"])
    assert [s.id for s in cohen.species] == [s.id for s in book.species]
    assert cohen.reactions == book.reactions


# --- Eratosthenes' sieve on a set (Cohen 2003 sec. 2.1) -----------------------------------
def test_sieve_leaves_the_primes():
    net = generate_network(ID, program="sieve")
    assert net.status == "complete"
    assert sorted(net.extras["run"]["final"]) == primes_upto(30)
    for rx in net.reactions:                              # x + y -> x with x | y, x a catalyst
        (x,) = rx.products
        (y,) = set(rx.reactants) - {x}
        assert int(y) % int(x) == 0 and rx.catalysts == {x: 1}
    consumed = {s for rx in net.reactions for s in rx.reactants if s not in rx.products}
    assert sorted(int(s) for s in {s.id for s in net.species} - consumed) == primes_upto(30)
    assert inert(net, [str(k) for k in primes_upto(30)])


def test_sieve_asynchronous_and_other_range():
    net = generate_network(ID, program="sieve", initial=list(range(2, 51)), strategy="asynchronous", seed=4)
    assert sorted(net.extras["run"]["final"]) == primes_upto(50)


# --- restriction enzymes (book 9.6; Giavitto & Michel 2002 sec. 3.1) -----------------------
def test_restriction_enzymes_published_tube():
    net = generate_network(ID, program="restriction-enzymes")
    assert net.status == "complete"
    assert net.extras["run"]["fixpoint"]
    assert Counter(net.extras["run"]["final"]) == Counter({"AATTCAA": 1, "TTG": 1, "CCCG": 1, "AATTCGGG": 1})
    pairs = {(tuple(rx.reactants), tuple(sorted(rx.products))) for rx in net.reactions}
    assert (("CCCGAATTCAA",), ("AATTCAA", "CCCG")) in pairs         # the book's example
    assert (("TTGAATTCGGG",), ("AATTCGGG", "TTG")) in pairs
    assert set(net.extras["reaction_rules"]) == {"EcoRI"}              # Void is elastic
    assert {s.id: s.structure for s in net.species}["TTG"] == "T,T,G"


def test_ecori_needs_flanks_and_cuts_every_site():
    m = Machine("tube", PROGRAMS["restriction-enzymes"]["rules"], ["GAATTCAA", "AGAATTCTGAATTCA"])
    species, reactions, _ = m.closure()
    assert not any(lhs == ("GAATTCAA",) for lhs, _, _ in reactions)    # X+ is non-empty
    cuts = sorted(sorted(rhs) for lhs, rhs, _ in reactions if lhs == ("AGAATTCTGAATTCA",))
    assert cuts == [["AATTCA", "AGAATTCTG"], ["AATTCTGAATTCA", "AG"]]
    final, _, fixpoint = m.run("maximal-parallel", np.random.default_rng(0))
    assert fixpoint and Counter(final) == Counter({"GAATTCAA": 1, "AG": 1, "AATTCTG": 1, "AATTCA": 1})


# --- Eden growth (Giavitto & Michel 2002 sec. 3.2, fig. 1) ---------------------------------
@pytest.mark.parametrize("collection, degree", [("grid", 4), ("hexagon", 6)])
def test_eden_on_grid_and_hexagon(collection, degree):
    net = generate_network(ID, program="eden", collection=collection)
    assert net.status == "complete"
    edges = net.extras["space"]["edges"]
    # every site can be occupied; every site but the initial cell starts empty
    assert len(net.species) == 2 * 49 - 1 and "undef@3_3" not in {s.id for s in net.species}
    assert len(net.reactions) == 2 * len(edges) - degree      # growth along every edge, never into the seed
    for rx in net.reactions:                              # true@p + undef@q -> true@p + true@q
        (occupied,) = rx.catalysts
        assert occupied.startswith("true@")
        (empty,) = set(rx.reactants) - {occupied}
        assert empty.startswith("undef@") and rx.products[empty.replace("undef", "true")] == 1
    centre = [s for s, c in net.initial_state.items() if s.startswith("true@")]
    assert centre == ["true@3_3"]
    assert sum(1 for rx in net.reactions if "true@3_3" in rx.reactants) == degree
    run = net.extras["run"]
    assert run["fixpoint"] and all(v is True for row in run["final"] for v in row)


def test_eden_maximal_parallel_growth():
    m = Machine("grid", PROGRAMS["eden"]["rules"], PROGRAMS["eden"]["initial"])
    rng, state = np.random.default_rng(3), m.start
    for t in range(1, 5):
        state = m.step(state, "maximal-parallel", rng)
        cells = [m.topology.positions[i] for i, v in enumerate(state) if v == "true"]
        assert len(cells) <= 2 ** t and all(abs(r - 3) + abs(c - 3) <= t for r, c in cells)
        if t == 1:
            assert len(cells) == 2                       # the single cell claims one neighbour


# --- Turn (Giavitto & Michel 2002 fig. 2) --------------------------------------------------
def test_turn_figure_2():
    net = generate_network(ID, program="turn")
    m = Machine("grid", PROGRAMS["turn"]["rules"], PROGRAMS["turn"]["initial"])
    rng, state = np.random.default_rng(0), m.start
    state = m.step(state, "maximal-parallel", rng)
    assert m.to_json(state) == [[None, 1, None], [2, 0, 4], [None, 3, None]]
    state = m.step(state, "maximal-parallel", rng)
    assert m.to_json(state) == [[None, 4, None], [1, 0, 3], [None, 2, None]]
    for _ in range(2):
        state = m.step(state, "maximal-parallel", rng)
    assert state == m.start                              # period 4, never a fixpoint
    assert not net.extras["run"]["fixpoint"]
    assert net.status == "complete" and all(len(rx.reactants) == 5 for rx in net.reactions)


# --- bead-sort (Cohen 2003 sec. 5.2, fig. 4) -----------------------------------------------
def test_bead_sort_figure_4():
    net = generate_network(ID, program="bead-sort")
    final = net.extras["run"]["final"]
    assert net.extras["run"]["fixpoint"]
    assert final == bead_rows([2, 2, 3, 4], 4)
    for rx in net.reactions:                              # a bead moves one row south
        moved = {s.split("@")[1] for s in rx.reactants}
        (r1, c1), (r2, c2) = sorted(tuple(map(int, s.split("_"))) for s in moved)
        assert c1 == c2 and r2 == r1 + 1
        assert rx.reactants == {f"true@{r1}_{c1}": 1, f"false@{r2}_{c2}": 1}
    net2 = generate_network(ID, program="bead-sort", initial=bead_rows([1, 5, 2, 4, 3], 5),
                            strategy="asynchronous", seed=2)
    assert net2.extras["run"]["final"] == bead_rows([1, 2, 3, 4, 5], 5)


# --- language and parameters -----------------------------------------------------------------
def test_bag_multiset_rewriting_and_torus():
    net = generate_network(ID, program="custom", collection="bag", initial=["a", "a", "b"],
                           rules=["x, y / x = y => x"])
    # species-level closure: two copies of any known value can meet
    assert [(rx.reactants, rx.products) for rx in net.reactions] == [({"a": 2}, {"a": 1}), ({"b": 2}, {"b": 1})]
    assert net.extras["run"]["final"] == {"a": 1, "b": 1}
    assert net.initial_state == {"a": 2.0, "b": 1.0}
    torus = generate_network(ID, program="eden", initial=[[True, None, None]], torus=True)
    assert len(torus.extras["space"]["edges"]) == 3


def test_bad_parameters():
    with pytest.raises(ValueError, match="only used with program='custom'"):
        generate_network(ID, rules=["x => x"])
    with pytest.raises(ValueError, match="needs a collection"):
        generate_network(ID, program="custom", rules=["x => x"], initial=[1])
    with pytest.raises(ValueError, match="as many elements"):
        generate_network(ID, program="custom", collection="seq", initial=[1, 2], rules=["x, y => x"])
    with pytest.raises(ValueError, match="unknown direction"):
        generate_network(ID, program="custom", collection="grid", initial=[[1, 2]], rules=["x |up> y => y, x"])
    with pytest.raises(ValueError, match="unexpected|cannot read"):
        generate_network(ID, program="custom", collection="bag", initial=[1],
                         rules=["x / __import__('os') => x"])
    with pytest.raises(ValueError, match="cannot read"):
        generate_network(ID, program="custom", collection="bag", initial=[1], rules=["x / os.system => x"])
    with pytest.raises(ValueError, match="unbound"):
        generate_network(ID, program="custom", collection="bag", initial=[1], rules=["x / z > 1 => x"])
    with pytest.raises(ValueError, match="only supported in a tube"):
        generate_network(ID, program="custom", collection="seq", initial=[1], rules=["X+ => X"])
    with pytest.raises(ValueError, match="torus"):
        generate_network(ID, torus=True)
    with pytest.raises(ValueError, match="integers"):
        generate_network(ID, program="sieve", initial=[2, "a"], strategy="maximal-parallel")
