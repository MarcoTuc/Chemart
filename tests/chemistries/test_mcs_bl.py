"""MCS.bl reproduces book 11.1.1 (tables 11.1, 11.2) and Decraene's thesis (2009) / ACS 2011.

Strings are written in the book's ASCII alphabet (# diamond, $ reversed triangle,
% triangle, ' prime); thesis glyphs are translated with to_ascii.
"""

from itertools import product

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import evolve, generate_network
from chemart.chemistries.mcs_bl import ALPHABET, C0, Rules, candidates, mutate, species_id, to_ascii, units


def reactions_by_structure(net):
    """{(enzyme-substrate multiset, product)} with species given by their strings."""
    structure = {s.id: s.structure for s in net.species}
    out = set()
    for r in net.reactions:
        lhs = tuple(sorted(structure[s] for s, n in r.reactants.items() for _ in range(n)))
        new = {s: n - r.reactants.get(s, 0) for s, n in r.products.items() if n > r.reactants.get(s, 0)}
        (prod,) = new
        out.add((lhs, structure[prod]))
    return out


# --- the broadcast language ------------------------------------------------------
# Book table 11.1 = thesis table 4.2 = CEC 2009 table I (quotes restored from the thesis).
TABLE_11_1 = [
    ("*$1:$0", "1:0", None),            # no reaction
    ("*'*0$:0$", "*0:1", "0:1"),         # inhibition
    ("*$:$", "*00:11", "*00:11"),        # universal replication
    ("*$0:$0", "*$0:$0", "*$0:$0"),      # self-replication
    ("*$1:$10", "*0:1", "*0:10"),        # concatenation
    ("*$1:$", "*0:1", "*0:"),            # cleavage
]


@pytest.mark.parametrize("enzyme, substrate, prod", TABLE_11_1)
def test_book_table_11_1(enzyme, substrate, prod):
    assert candidates(enzyme, substrate) == (() if prod is None else (prod,))


def test_table_11_1_activation_erratum():
    # The printed enzyme *$1:'*$ gives *0: (the same rule as cleavage: $ stops before the
    # matched 1); the published product *0:1 comes from the unit *$:'*$.
    assert candidates("*$1:'*$", "0:1") == ("*0:",)
    assert candidates("*$:'*$", "0:1") == ("*0:1",)
    assert units("*0:1") == [("0", "1")]


# Thesis 4.3.3, reactions r1-r20 (None = elastic).
THESIS_R = [
    ("*11'%0:1'1", "11%0", "11"), ("*11'%0:11", "1100", None), ("*11%0:1'%", "1100", "1%"),
    ("*10111:00", "10111", "00"), ("10111:00", "10110", None), ("*1011100", "10110", None),
    ("*:1011100", "10110", None), ("*10111:100", "10111", "100"), ("*10111:100", "10110", None),
    ("*10#11:00", "10011", "00"), ("*10#11:00", "10:11", "00"), ("*10#11:00", "11:11", None),
    ("*10#:00", "10", None), ("*1011#:100", "1011000", "100"), ("*1011#:100", "1011010101010", "100"),
    ("*10$:$", "10011", "011"), ("*10$:$$", "100100101", "01001010100101"),
    ("*10$:$$", "110100101", None), ("*11%0:1%", "1100", "10"), ("*11%0:1%", "1110", "11"),
]


@pytest.mark.parametrize("k, enzyme, substrate, prod", [(i + 1, *r) for i, r in enumerate(THESIS_R)])
def test_thesis_reactions_r1_r20(k, enzyme, substrate, prod):
    assert candidates(enzyme, substrate) == (() if prod is None else (prod,)), f"r{k}"


def test_thesis_unit_parsing_example():
    # thesis 4.3.2: m1 = 10*11'*%0:1% *%00*11$:11%:1 has units *11'*%0:1% and *11$:11%
    m1 = to_ascii("10∗11′∗△0 : 1△ ∗ △00∗11▽ : 11△ : 1")
    assert units(m1) == [("11'*%0", "1%"), ("11$", "11%")]
    assert units(to_ascii("11′∗△0 : 1△′∗ : 11△ : 0001♦")) == []   # a null device (fig. 4.4 m3)


# --- published string sets --------------------------------------------------------
TABLE_11_2 = {
    "c0": ["∗▽0 : ▽1", "∗▽0 : ▽0", "∗▽1 : ▽0", "∗▽1 : ▽1"],
    "c1": ["∗▽0 : ▽1", "∗▽0 : ▽0", "∗▽♦ : ▽0", "∗▽♦ : ▽1"],
    "c2": ["∗▽0 : ▽1", "∗▽0 : ▽0", "∗▽♦ : ▽△0", "∗▽♦ : ▽△1"],
    "c3": ["∗▽0 : ▽1", "∗▽0 : ▽0", "∗▽△ : ▽0", "∗▽△ : ▽1"],
}


def test_book_table_11_2_cell_types_are_closed_with_6_and_9_reactions():
    # thesis 7.2.3 / ACS 2011: c0 has 6 reactions, c1 three more; self-replication disabled.
    counts = {}
    for name, strings in TABLE_11_2.items():
        net = generate_network("mcs-bl", strings=strings)
        assert net.status == "complete"
        assert sorted(s.structure for s in net.species) == sorted(map(to_ascii, strings))  # closed
        counts[name] = len(net.reactions)
    assert counts == {"c0": 6, "c1": 9, "c2": 9, "c3": 9}
    assert [to_ascii(s) for s in TABLE_11_2["c0"]] == list(C0)
    assert generate_network("mcs-bl").to_dict()["reactions"] == generate_network("mcs-bl", strings=list(C0)).to_dict()["reactions"]


def test_c1_c2_c3_are_phenotypically_equivalent():
    # thesis 7.2.3: # and % act alike in the condition, the % of c2's actions is ignored.
    nets = {n: reactions_by_structure(generate_network("mcs-bl", strings=s)) for n, s in TABLE_11_2.items()}
    rename = {"*$#:$%0": "*$#:$0", "*$#:$%1": "*$#:$1", "*$%:$0": "*$#:$0", "*$%:$1": "*$#:$1"}

    def mapped(rs):
        return {(tuple(sorted(rename.get(s, s) for s in lhs)), rename.get(p, p)) for lhs, p in rs}

    assert mapped(nets["c2"]) == nets["c1"] == mapped(nets["c3"])
    assert nets["c0"] != nets["c1"]


def test_c0_reactions():
    s1, s2, s3, s4 = C0   # thesis table 7.1 numbering
    expected = {
        ((s1, s2), s1), ((s1, s3), s4), ((s2, s3), s3),
        ((s1, s3), s2), ((s3, s4), s3), ((s1, s4), s1),
    }
    got = reactions_by_structure(generate_network("mcs-bl"))
    assert got == {(tuple(sorted(lhs)), p) for lhs, p in expected}


def test_thesis_appendix_d2_reaction_table():
    s1, s2, s3, s4 = "*$0:$1", "*$1:$0", "*$1:$1", "*$0:$0"
    table = {  # enzyme, substrate -> product (None = elastic), thesis appendix D.2
        (s1, s1): None, (s1, s2): s3, (s1, s3): None, (s1, s4): s1,
        (s2, s1): s4, (s2, s2): None, (s2, s3): s2, (s2, s4): None,
        (s3, s1): s1, (s3, s2): None, (s3, s3): None, (s3, s4): None,
        (s4, s1): None, (s4, s3): None, (s4, s4): None,
    }
    rules = Rules(max_length=500, self_replication=False)
    for (e, s), p in table.items():
        assert rules.products(e, s) == ([] if p is None else [p]), (e, s)
    # erratum: D.2 prints s4 + s2 -> s4; the rule replicates the substrate s2
    assert rules.products(s4, s2) == [s2]
    # "No further molecular species have been produced": the closure of {s1, s2} is c0
    net = generate_network("mcs-bl", strings=[s1, s2])
    assert net.status == "complete" and sorted(s.structure for s in net.species) == sorted(C0)
    assert len(net.reactions) == 6


# Thesis table 5.1 (ALife 2008 table 2): the entries whose action has a single $.
SELF_REPLICASES_5_1 = [
    "1▽0 ∗ ▽ : ▽", "1△▽0 ∗ ▽ : ▽", ": 1 ∗ ▽ : ▽♦ : 1 ∗ ▽ : ▽♦", ": 0▽▽ ∗ ▽ : △△▽",
    "1♦▽ : ∗ ▽▽ : ▽ :", "∗ ▽ : ▽", "∗ ▽▽ : ▽", "∗ ▽0▽▽ : ▽0", "△1 ∗ ▽ : ▽♦",
]
TABLE_5_2 = ["∗▽ : ▽", "∗▽1 : ▽1", "∗▽01 : ▽01", "∗▽101 : ▽101", "∗▽0101 : ▽0101"]


@pytest.mark.parametrize("glyphs", SELF_REPLICASES_5_1 + TABLE_5_2)
def test_thesis_self_replicases(glyphs):
    s = to_ascii(glyphs)
    assert candidates(s, s) and set(candidates(s, s)) == {s}


def test_only_one_self_replicase_of_length_4():
    # thesis 5.2: the only autocatalytic molecule of minimal length is s_R0 = *$:$
    found = [s for s in map("".join, product(ALPHABET, repeat=4)) if candidates(s, s) == (s,)]
    assert found == ["*$:$"]


def test_universal_replicase_copies_anything():
    for s in ["0101", "*$1:$0", ":'#%", "*$:$"]:
        assert candidates("*$:$", s) == (s,)


def test_elongation_catastrophe_section_5_5():
    r4, r4p, r4pp, r4ppp = "*$0101:$0101", "*$0101:$00101", "*$0101:$000101", "*$0101:$00000101"
    assert candidates(r4, r4) == (r4,)
    assert candidates(r4, r4p) == (r4p,)
    assert candidates(r4p, r4) == (r4p,)
    assert candidates(r4p, r4p) == (r4pp,)
    assert candidates(r4pp, r4pp) == (r4ppp,)
    # BDL_max makes longer products elastic, so the closure stops at max_length ...
    net = generate_network("mcs-bl", strings=[r4, r4p], self_replication=True, max_length=16)
    assert net.status == "complete"
    assert [len(s.structure) for s in net.species] == [12, 13, 14, 15, 16]
    # ... and without it the species budget cuts the runaway
    net = generate_network("mcs-bl", strings=[r4, r4p], self_replication=True, max_species=20)
    assert net.status == "truncated" and len(net.species) == 20


def test_several_satisfied_units_share_the_collision():
    # CEC 2009 II-A: one satisfied unit is chosen at random
    device = "*$:$0*$:$1"
    rules = Rules(max_length=10, self_replication=True)
    assert rules.outcomes(device, "0") == (2, {"00": 1, "01": 1})
    assert rules.k(device, "0", "00") == rules.k("0", device, "01") == 0.5
    # the device on itself would give 11 symbols: both units elastic under max_length = 10
    assert rules.outcomes(device, device) == (2, {})
    net = generate_network("mcs-bl", strings=[device, "0"], max_length=10, max_species=4)
    structure = {s.id: s.structure for s in net.species}
    ks = {structure[next(iter(set(r.products) - set(r.reactants)))]: r.rate["k"]
          for r in net.reactions if set(r.reactants) == {species_id(device), species_id("0")}}
    assert ks == {"00": 0.5, "01": 0.5}


# --- dynamics ---------------------------------------------------------------------
def test_universal_replicase_has_zero_expected_growth_eqs_5_3_to_5_6():
    net = generate_network("mcs-bl", strings=["*$:$", "0101"], self_replication=True)
    assert len(net.reactions) == 2
    ids, f = rhs(net)
    for x1 in (0.1, 0.5, 0.9):
        x = np.array([x1, 1 - x1]) if ids[0] == species_id("*$:$") else np.array([1 - x1, x1])
        assert np.allclose(f(0.0, x), 0.0)


def test_network_is_the_catalytic_network_equation_5_1():
    strings = [to_ascii(s) for s in TABLE_11_2["c1"]] + ["*$:$", "0:1"]
    net = generate_network("mcs-bl", strings=strings)
    assert net.status == "complete" and len(net.species) > 6   # the whole system, not a truncation
    ids, f = rhs(net)
    structure = {s.id: s.structure for s in net.species}
    sp = [structure[i] for i in ids]
    rules = Rules(max_length=500, self_replication=False)
    rng = np.random.default_rng(3)
    x = rng.random(len(sp))
    x /= x.sum()
    prod = np.zeros(len(sp))
    for i, a in enumerate(sp):
        for j, b in enumerate(sp):
            for p in rules.products(a, b):   # alpha^k_ij = 1 (eq. 5.2)
                prod[sp.index(p)] += x[i] * x[j]
    assert np.allclose(f(0.0, x), prod - x * prod.sum())


def test_single_reactor_soup():
    traj = evolve("mcs-bl", seed=4, steps=3000, n_max=200, p_s=0.0)
    net = traj.network
    assert net.status == "observed"
    assert traj.times()[:2] == [0.0, 40.0] and traj.times()[-1] == 3000.0   # a frame per 40 molecules
    closure = reactions_by_structure(generate_network("mcs-bl"))
    assert reactions_by_structure(net) <= closure
    assert sum(r.count for r in net.reactions) == net.extras["analysis"]["productive"] > 160
    assert sum(net.extras["final_state"].values()) == 200          # grew to capacity, then displaced
    assert net.extras["analysis"]["mutant_products"] == 0


def test_soup_with_mutation_explores_new_strings():
    net = evolve("mcs-bl", seed=2, steps=3000, n_max=300, p_s=0.05).network
    assert net.extras["analysis"]["mutant_products"] > 0
    assert len(net.species) > 4
    assert any(r.rate is None for r in net.reactions)   # mutant products have no structural k
    assert all(set(s.structure) <= set(ALPHABET) for s in net.species)


def test_mutation_operators():
    rng = np.random.default_rng(0)
    assert mutate("*$0:$1", 0.0, rng) == "*$0:$1"
    lengths = {len(mutate("0000000000", 1.0, rng)) for _ in range(200)}
    assert min(lengths) < 10 < max(lengths)   # deletions and insertions both occur
    flipped = [mutate("0", 1.0, rng) for _ in range(300)]
    assert all(set(m) <= set(ALPHABET) for m in flipped)
    assert "" in flipped and any(len(m) == 2 for m in flipped) and any(m not in ("", "0") and len(m) == 1 for m in flipped)


def test_random_seed_strings_and_bad_strings():
    net = generate_network("mcs-bl", strings=[], n_random=30, random_length=6, seed=5, max_species=50)
    assert len(net.extras["seed"]) == 30 and all(len(s.structure) >= 1 for s in net.species)
    with pytest.raises(ValueError, match="strings"):
        generate_network("mcs-bl", strings=["*$0:$x"])
    with pytest.raises(ValueError, match="max_length"):
        generate_network("mcs-bl", strings=["*$0101:$0101"], max_length=5)
    with pytest.raises(ValueError, match="n_random"):
        generate_network("mcs-bl", strings=[])
