"""Self-replicating loops: Langton (1984), Byl (1989), Chou-Reggia (1993), Sayama (1998-2004).

Published results reproduced here:

- the size of each automaton and of each ancestor loop (book 10.4: "8 states per
  cell where the self-reproducing configuration is of size 5 cells" for the
  Chou-Reggia loop; Langton's 8 states and 219 transitions; the table of loop
  variants in Wikipedia's "Langton's loops", sourced from the original papers);
- the replication periods: Langton 151, Byl 25, Chou-Reggia 15, evoloop 363,
  and the SDSR loop's 151 (it is Langton's loop plus structural dissolution);
- the daughter is an exact copy of the ancestor (Langton 1984);
- in bounded space SDSR loops dissolve and the population reaches a steady
  state with continuous turn-over (Sayama 1998, 1999; book 10.7.2 fig. 10.12);
- evoloops evolve towards smaller loops (Sayama 1999; Salzberg, Antony &
  Sayama 2004; book 8.2.3 and fig. 10.13).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import sr_loops as S

ID = "sr-loops"

#: rule -> (states, transitions, cells of the ancestor, replication period)
PUBLISHED = {
    "langton": (8, 219, 86, 151),
    "byl": (6, 144, 12, 25),
    "reggia-1": (8, 174, 6, 13),
    "reggia-2": (8, 65, 5, 15),
    "sdsr": (9, 142, 86, 151),
    "evoloop": (9, 132, 149, 363),
}


def lattice(rule, side):
    """The ancestor of `rule` centred on a side x side lattice."""
    pattern = S.parse_rle(S.ANCESTORS[rule])
    h, w = pattern.shape
    grid = np.zeros((side, side), dtype=np.int8)
    y, x = (side - h) // 2, (side - w) // 2
    grid[y:y + h, x:x + w] = pattern
    return grid, pattern, (y, x)


def at_step(rule, side, steps):
    """The lattice after `steps` updates, and the ancestor pattern and its position."""
    grid, pattern, home = lattice(rule, side)
    n_states, table, _ = S.rule(rule)
    for grid in S.evolve(grid, table, n_states, steps):
        pass
    return grid, pattern, home


# --- the automata and the ancestors -------------------------------------------------
@pytest.mark.parametrize("name", sorted(PUBLISHED))
def test_published_size_of_each_automaton_and_ancestor(name):
    states, transitions, cells, period = PUBLISHED[name]
    n_states, table, n_transitions = S.rule(name)
    assert (n_states, n_transitions) == (states, transitions)
    pattern = S.parse_rle(S.ANCESTORS[name])
    assert int((pattern > 0).sum()) == cells
    assert pattern.max() < n_states
    assert S.PERIODS[name] == period
    # a table is a partial function: unspecified neighbourhoods leave the cell alone
    assert table.shape == (n_states ** 5,)
    # an empty neighbourhood: quiescent, either by a transition or by not being covered
    assert table[0] in (0, -1)


def test_the_book_smallest_replicator_has_eight_states_and_five_cells():
    # book 10.4: "the smallest currently known example ... is a cellular automaton with 8
    # states per cell where the self-reproducing configuration is of size 5 cells [705]"
    n_states, _, _ = S.rule("reggia-2")
    assert n_states == 8
    assert int((S.parse_rle(S.ANCESTORS["reggia-2"]) > 0).sum()) == 5


# --- replication periods -------------------------------------------------------------
@pytest.mark.parametrize("name, side, period", [
    ("langton", 80, 147), ("byl", 40, 25), ("reggia-1", 40, 13),
    ("reggia-2", 40, 15), ("sdsr", 80, 147), ("evoloop", 140, 363),
])
def test_an_exact_copy_of_the_ancestor_appears_at_the_published_step(name, side, period):
    # the published periods are 151 (Langton, SDSR), 25 (Byl), 15 (Chou-Reggia) and
    # 363 (evoloop); for Langton and the SDSR loop the first complete copy stands on the
    # lattice at 147 and the newly born loop is a copy at 151 (see the catalog decisions)
    before, pattern, home = at_step(name, side, period - 1)
    after, _, _ = at_step(name, side, period)
    assert S.copies(before, pattern, home=home) == []
    assert S.copies(after, pattern, home=home), "a complete copy of the ancestor"


def test_byl_loop_replicates_every_25_steps():
    # Byl (1989): a 12-cell loop with a period of 25
    grid, pattern, home = lattice("byl", 60)
    n_states, table, _ = S.rule("byl")
    found = [t for t, g in enumerate(S.evolve(grid, table, n_states, 80), start=1)
             if S.copies(g, pattern, home=home)]
    assert found == [25, 50, 75]


def test_sdsr_is_langtons_loop_with_dissolution_added():
    # Sayama (1998): "a dissolving state 8 was introduced into the set of states"
    langton, table_l, _ = S.rule("langton")
    sdsr, table_s, _ = S.rule("sdsr")
    assert (langton, sdsr) == (8, 9)
    a, _, _ = lattice("langton", 80)
    b, _, _ = lattice("sdsr", 80)
    for _ in range(200):
        a = S.step(a, table_l, langton)
        b = S.step(b, table_s, sdsr)
        assert np.array_equal(a, b), "an undisturbed SDSR loop behaves like Langton's"
    assert not (b == 8).any()


# --- the observed network ------------------------------------------------------------
def test_default_network_is_langtons_loop_replicating_itself():
    net = generate_network(ID, seed=1)
    assert net.status == "observed"
    assert [s.id for s in net.species] == ["L086aaa"]
    assert [r.to_text() for r in net.reactions] == ["L086aaa -> 2 L086aaa  (x3)"]
    analysis = net.extras["analysis"]
    # Langton 1984: the loop makes an identical copy of itself in 151 steps
    assert analysis["ancestor_copy_steps"][0] == 151
    first = analysis["replications"][0]
    assert (first["step"], first["copy"], first["identified"]) == (128, 151, 151)
    assert (first["mother"], first["daughter"]) == ("L086aaa", "L086aaa")
    assert analysis["deaths"] == 0, "Langton's loops never dissolve"
    # species structure is the ancestor configuration; initial + events = final population
    rows = "/".join("".join(str(v) for v in row) for row in S.parse_rle(S.ANCESTORS["langton"]))
    assert net.species[0].structure in {rows, *(
        "/".join("".join(str(v) for v in row) for row in np.rot90(S.parse_rle(S.ANCESTORS["langton"]), k))
        for k in range(4))}
    state = Counter({k: int(v) for k, v in net.initial_state.items()})
    for r in net.reactions:
        for s, n in r.reactants.items():
            state[s] -= n * r.count
        for s, n in r.products.items():
            state[s] += n * r.count
    assert {k: v for k, v in state.items() if v} == analysis["final_population"]
    space = net.extras["space"]
    assert (space["shape"], space["states"], space["transitions"]) == ([60, 60], 8, 219)
    assert len(space["final"]) == 60 and set("".join(space["final"])) <= set("01234567")


@pytest.mark.parametrize("name, kw, events", [
    ("byl", dict(grid=40, steps=80, min_loop_cells=4), "L012aaa -> 2 L012aaa"),
    ("reggia-2", dict(grid=40, steps=80, min_loop_cells=4), "L005aaa -> 2 L005aaa"),
])
def test_the_small_loops_replicate_as_themselves(name, kw, events):
    net = generate_network(ID, rule=name, **kw)
    text = [r.to_text() for r in net.reactions]
    assert any(t.startswith(events) for t in text)
    cells = PUBLISHED[name][2]
    assert list(net.initial_state) == [f"L{cells:03d}aaa"]
    assert sum(net.extras["analysis"]["final_population"].values()) > 5


# --- bounded space: dissolution and evolution ----------------------------------------
@pytest.mark.slow
def test_sdsr_loops_dissolve_and_the_population_reaches_a_steady_state():
    # Sayama (1998, 1999): with structural dissolution the colony keeps turning over
    # instead of freezing, and the population settles in a bounded space
    net = generate_network(ID, rule="sdsr", grid=100, steps=4000, min_loop_cells=40,
                           track_every=5)
    analysis = net.extras["analysis"]
    assert analysis["births"] > 100 and analysis["deaths"] > 100
    text = [r.to_text() for r in net.reactions]
    assert any(t.startswith("L086aaa -> 2 L086aaa") for t in text)
    assert any(t.startswith("L086aaa -> ∅") for t in text)
    counts = [p["loops"] for p in analysis["population"]]
    assert counts[0] == 1 and max(counts) > 10
    second_half = counts[len(counts) // 2:]
    assert 5 <= min(second_half) and max(second_half) <= 3 * (sum(second_half) / len(second_half))
    assert analysis["final_population"].get("L086aaa", 0) >= 5, "the ancestor still dominates"


@pytest.mark.slow
def test_evoloop_colony_evolves_towards_smaller_loops():
    # Sayama (1999), Salzberg et al. (2004), book 8.2.3: "as the various loops compete for
    # space, smaller loops that have a reproductive advantage emerge and dominate"
    net = generate_network(ID, rule="evoloop", grid=200, steps=30000, min_loop_cells=20,
                           track_every=25)
    analysis = net.extras["analysis"]
    assert analysis["births"] > 500 and analysis["deaths"] > 500
    assert analysis["species_seen"] > 20, "variation on collision creates new species"
    text = [r.to_text() for r in net.reactions]
    assert any(t.startswith("L149aaa -> 2 L149aaa") for t in text), "the ancestor replicates"
    assert any(" + L" in t and t.startswith("L149aaa -> L149aaa") for t in text), "L -> L + L'"
    start, end = analysis["population"][0], analysis["population"][-1]
    assert end["loops"] > 20 > start["loops"]
    assert start["cells"] == [149]
    # the typical loop left is far smaller than the size-8 ancestor (leftover sheath
    # structures of dissolved loops are still counted, so the mean stays higher)
    cells = sorted(end["cells"])
    assert cells[0] <= 30 and cells[len(cells) // 2] < 100
    dominant = max(end["by_species"], key=lambda s: end["by_species"][s])
    assert int(dominant[1:4]) < 60, "the dominant species is much smaller than the ancestor"
    assert end["by_species"].get("L149aaa", 0) < end["by_species"][dominant]


# --- the micro reading ---------------------------------------------------------------
def test_micro_reading_is_one_reaction_per_state_transition():
    # book 10.7.2: a cell state is a molecule, the quiet state is the absence of one,
    # and state transitions are movement and reaction
    net = generate_network(ID, mode="micro", steps=200)
    assert net.status == "observed" and net.reactions
    assert [s.id for s in net.species] == [f"s{k}" for k in range(1, 8)]
    assert net.initial_state == {"s1": 17.0, "s2": 61.0, "s4": 2.0, "s7": 6.0}
    for r in net.reactions:
        gained = Counter(r.products) - Counter(r.reactants)
        lost = Counter(r.reactants) - Counter(r.products)
        assert len(gained) <= 1 and len(lost) <= 1, "exactly one cell changes state"
        assert r.count >= 1
        catalysts = Counter(r.reactants) - lost
        assert sum(catalysts.values()) <= 4, "the four neighbours are the catalysts"
    firings = net.extras["analysis"]["firings"]
    assert firings == sum(r.count for r in net.reactions)
    # the same run in macro mode covers the same lattice
    assert net.extras["space"]["rule"] == "langton"


# --- parameters ----------------------------------------------------------------------
def test_placement_uses_the_seed_and_the_run_is_otherwise_deterministic():
    a = generate_network(ID, seed=4, grid=80, steps=160, ancestors=3)
    b = generate_network(ID, seed=4, grid=80, steps=160, ancestors=3)
    c = generate_network(ID, seed=5, grid=80, steps=160, ancestors=3)
    assert a.to_dict() == b.to_dict()
    assert a.extras["space"]["final"] != c.extras["space"]["final"]
    assert sum(a.initial_state.values()) == 3
    # one ancestor: the CA is deterministic, the seed changes nothing
    assert (generate_network(ID, seed=0, steps=160).extras["space"]["final"]
            == generate_network(ID, seed=99, steps=160).extras["space"]["final"])


def test_a_custom_ancestor_can_be_given_as_rle():
    net = generate_network(ID, rule="byl", grid=40, steps=60, min_loop_cells=4,
                           ancestor_pattern=S.ANCESTORS["byl"])
    assert list(net.initial_state) == ["L012aaa"]
    with pytest.raises(ValueError, match="RLE pattern"):
        generate_network(ID, ancestor_pattern="not an rle")


def test_bad_parameters():
    with pytest.raises(ValueError, match="grid"):
        generate_network(ID, grid=7)
    with pytest.raises(ValueError, match="smaller than"):
        generate_network(ID, rule="evoloop", grid=16)
    with pytest.raises(ValueError, match="min_loop_cells"):
        generate_network(ID, rule="reggia-2", grid=40, steps=10)
    with pytest.raises(ValueError, match="cannot place"):
        generate_network(ID, grid=20, ancestors=50)
