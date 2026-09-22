"""BNC cell model, Hintze & Adami (2008): molecule and reaction counts, examples, affinity and eq. 3."""

from collections import Counter
from itertools import product

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import generate_network
from chemart.chemistries import bnc_cell as bnc


def brute_molecules(L):
    """Independent enumeration: assign bond orders left to right and check every valence."""
    out = []
    for n in range(2, L + 1):
        for atoms in product((1, 2, 3), repeat=n):
            for orders in product((1, 2, 3), repeat=n - 1):
                full = (0, *orders, 0)
                if all(full[i] + full[i + 1] == atoms[i] for i in range(n)):
                    out.append("".join(map(str, atoms)))
                    break
    return out


def cut_sites(mols):
    sites = []
    for m in mols:
        orders = bnc.bonds(m)
        sites += [(m[:k], m[k:], orders[k - 1]) for k in range(1, len(m))]
    return sites


def molecule_vector(net, ids, S, atom):
    vec = np.array([net.extras["conservation"][int(atom) - 1]["vector"].get(s, 0) for s in ids])
    return vec @ S


# --- Results: 608 molecules, precursors, Methods numbering ------------------------------
def test_608_molecules_and_the_papers_numbering():
    mols = bnc.molecules(12)
    assert len(mols) == 608
    assert brute_molecules(6) == list(bnc.molecules(6))                         # independent enumeration, same order
    names = [bnc.name(m) for m in mols]
    assert names[:6] == ["1-1", "2=2", "3#3", "1-2-1", "1-3=2", "2=3-1"]    # M0..M5
    assert names[607] == "2=3-3=3-3=3-3=3-3=3-3=2"                         # M607
    assert len(mols[:53]) == 53 and {len(m) for m in mols[:53]} == set(range(2, 8))
    assert len([m for m in mols if len(m) <= 7]) == 53                     # the 53 smallest are all molecules up to 7 atoms


def test_valid_and_invalid_molecules_of_the_text():
    for text in ["1-2-2-1", "2=2", "1-2-3=3-2-1", "1-2-1", "1-2-3=3-1", "1-3=2"]:   # paper and book examples
        assert bnc.name(bnc.parse(text)) == text
    for text in ["1-3=1", "1-2-3=2-1"]:          # paper's invalid example; the book's import example (see decisions)
        with pytest.raises(ValueError, match="not a valid molecule"):
            bnc.parse(text)
    assert bnc.code(bnc.parse("1-2-3=3-2-1")) == "123321000000"         # Methods: import specificity 123321000000


# --- Results: 5,020,279 legal reactions -------------------------------------------------
def test_reaction_count_matches_the_paper_up_to_one():
    assert bnc.count_cleavage_reactions(12) == 5_020_279 + 1


@pytest.mark.parametrize("L", [4, 6, 8])
def test_reaction_count_formula_against_brute_force(L):
    sites = cut_sites(bnc.molecules(L))
    count = 0
    for (ah, at, o1), (bh, bt, o2) in product(sites, repeat=2):
        c, d = ah + bt, bh + at
        if o1 == o2 and len(c) <= L and len(d) <= L and sorted((c, d)) != sorted((ah + at, bh + bt)):
            assert bnc.bonds(c) is not None and bnc.bonds(d) is not None    # equal cut orders give valid products
            count += 1
        elif o1 != o2 and len(c) <= L:
            assert bnc.bonds(c) is None                                     # unequal orders never do
    assert bnc.count_cleavage_reactions(L) == count


def test_chemistry_network_is_every_distinct_recombination():
    L = 7
    net = generate_network("bnc-cell", max_length=L)
    mols = bnc.molecules(L)
    expected = set()
    for (ah, at, o1), (bh, bt, o2) in product(cut_sites(mols), repeat=2):
        c, d = ah + bt, bh + at
        if o1 == o2 and len(c) <= L and len(d) <= L:
            lhs, rhs_ = Counter([bnc.name(ah + at), bnc.name(bh + bt)]), Counter([bnc.name(c), bnc.name(d)])
            if lhs != rhs_:
                expected.add((frozenset(lhs.items()), frozenset(rhs_.items())))
    ours = {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in net.reactions}
    assert ours == expected and len(net.reactions) == len(expected) == 3192
    assert len(net.species) == 53 and all(r.rate is None for r in net.reactions)
    assert all((b, a) in ours for a, b in ours)                             # every reaction is reversible


def test_default_network_examples_conservation_and_fitness():
    net = generate_network("bnc-cell")
    assert len(net.species) == 87 and len(net.reactions) == 11814
    texts = {r.to_text() for r in net.reactions}
    assert "1-2-2-1 + 2=3-3=2 -> 1-2-3=2 + 2=3-2-1" in texts                # book eq. 18.21 / paper example
    assert "1-2-1 + 1-2-2-1 -> 1-1 + 1-2-2-2-1" in texts                    # Figure 1, protein C
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    for atom in "123":
        assert not np.any(molecule_vector(net, ids, S, atom))
    phi = net.extras["fitness"]["phi"]
    mols = bnc.molecules(8)
    assert phi[bnc.name(mols[52])] == 0.0 and phi[bnc.name(mols[53])] == pytest.approx(53 ** 2 / 608 ** 2)   # eq. 5
    assert net.extras["reaction_space"]["cleavage_reactions"] == 2 * 19368


# --- Figure 11: the ancestral genome's proteins ----------------------------------------
def test_figure_11_proteins():
    net = generate_network("bnc-cell", mode="cell", max_length=12, proteins=[
        "enzyme 1|2-2-2-2-1 + 1-2-2-2-2|1",
        "enzyme 1|3=3-3=3-1 + 1-3=3-3=3|1",
        "import 1-2-2-2-2-1",
        "import 1-2-2-3=3-2-2-2-2-2-3=2",
    ])
    enzymes = [r for r in net.reactions if any(s.startswith("enzyme") for s in r.reactants)]
    assert Counter(enzymes[0].products) == Counter({"1-1": 1, "1-2-2-2-2-2-2-2-2-1": 1, "enzyme_0": 1})
    assert Counter(enzymes[0].reactants) == Counter({"1-2-2-2-2-1": 2, "enzyme_0": 1})
    assert Counter(enzymes[1].products) == Counter({"1-1": 1, "1-3=3-3=3-3=3-3=3-1": 1, "enzyme_1": 1})
    # The third reaction of Figure 11 is a typo: its product breaks the bond rule and an atom is lost.
    with pytest.raises(ValueError):
        bnc.parse("1-2-3=3-2-3=2-3=2")


# --- Methods: affinity score and eqs. 2-4 ------------------------------------------------
def paper_D(molecule, domain):
    """D(M, P) = 1 - sqrt(sum_i f(m_i EQUAL p_i)^2 / 108) on 2-bit sites."""
    bits = lambda v: format(v, "02b")   # noqa: E731
    m = [int(c) for c in molecule.ljust(12, "0")]
    total = 0
    for mi, pi in zip(m, domain):
        equal = "".join("1" if x == y else "0" for x, y in zip(bits(mi), bits(pi)))
        total += int(equal, 2) ** 2
    return 1 - (total / 108) ** 0.5


def test_affinity_score_of_the_methods():
    perfect = [int(b, 2) for b in "10-01-00-00-01-10-11-11-11-11-11-11".split("-")]   # the paper's example
    assert bnc.affinity("123321", perfect) == pytest.approx(1.0)
    assert bnc.affinity("123321", [1, 2, 3, 3, 2, 1] + [0] * 6) == pytest.approx(0.0)
    rng = np.random.default_rng(0)
    for _ in range(50):
        m = bnc.molecules(12)[int(rng.integers(608))]
        dom = [int(v) for v in rng.integers(0, 4, 12)]
        assert 0.0 <= bnc.affinity(m, dom) <= 1.0
        assert bnc.affinity(m, dom) == pytest.approx(paper_D(m, dom))


def test_enzyme_kinetics_are_eqs_2_to_4():
    net = generate_network("bnc-cell", seed=5, mode="cell", proteins=[
        "enzyme 1-2|1 + 1|2-2-1",          # Figure 1, protein C
        "enzyme 1-2|2-1 + 2=3|3=2",        # book eq. 18.21
    ])
    prot = {p["id"]: p for p in net.extras["proteins"]}
    A = {}
    # domains match A, B, A' (A's head + B's tail), B' (B's head + A's tail)
    for pid, molecules_ in (("enzyme_0", ["121", "1221", "12221", "11"]), ("enzyme_1", ["1221", "2332", "1232", "2321"])):
        A[pid] = sum(paper_D(m, [int(c) for c in d]) for m, d in zip(molecules_, prot[pid]["domains"])) / 4   # eq. 4
        assert prot[pid]["affinity"] == pytest.approx(A[pid])
    assert set(net.extras["buffered"]) == {"enzyme_0", "enzyme_1"}
    rng = np.random.default_rng(1)
    x = {s.id: float(v) for s, v in zip(net.species, rng.uniform(0.2, 1.5, len(net.species)))}
    x.update(net.initial_state)
    P = net.initial_state
    kout = {"1-2-1": 1, "1-2-2-1": 2, "2=3-3=2": 1}                          # 1-2-2-1 enters both reactions
    v0 = x["1-2-1"] / kout["1-2-1"] * x["1-2-2-1"] / kout["1-2-2-1"] * A["enzyme_0"] * P["enzyme_0"]   # eq. 3
    v1 = x["1-2-2-1"] / kout["1-2-2-1"] * x["2=3-3=2"] / kout["2=3-3=2"] * A["enzyme_1"] * P["enzyme_1"]
    paper = {"1-2-1": -v0, "1-2-2-1": -v0 - v1, "1-1": v0, "1-2-2-2-1": v0,          # eq. 2
             "2=3-3=2": -v1, "1-2-3=2": v1, "2=3-2-1": v1, "enzyme_0": 0.0, "enzyme_1": 0.0}
    ids, f = rhs(net)
    ours = dict(zip(ids, f(0.0, np.array([x[s] for s in ids]))))
    assert set(ours) == set(paper)
    for s in paper:
        assert ours[s] == pytest.approx(paper[s], rel=1e-10, abs=1e-14)


# --- cell mode: Figure 1 pathway, random cells, validation --------------------------------
def test_figure_1_pathway():
    net = generate_network("bnc-cell", mode="cell", proteins=["import 1-2-1", "import 1-2-2-1", "enzyme 1-2|1 + 1|2-2-1", "export 1-1"])
    texts = [r.to_text().split("  [")[0] for r in net.reactions]
    assert texts == [
        "1-2-1_out + import_0 -> 1-2-1 + import_0",
        "1-2-2-1_out + import_1 -> 1-2-2-1 + import_1",
        "1-2-1 + 1-2-2-1 + enzyme_2 -> 1-2-2-2-1 + 1-1 + enzyme_2",
        "1-1 + export_3 -> 1-1_out + export_3",
    ]
    assert net.reactions[2].rate["kout_product"] == 1 and net.reactions[0].rate is None
    assert net.extras["removed_every_update"] == ["1-1_out"]
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    for atom in "123":
        assert not np.any(molecule_vector(net, ids, S, atom))
    assert "compartments" in net.provides and "space" in net.provides


def test_random_cell():
    net = generate_network("bnc-cell", seed=3, mode="cell", enzymes=40, importers=10, exporters=5, max_length=12)
    kinds = Counter(p["type"] for p in net.extras["proteins"])
    assert kinds == {"enzyme": 40, "import": 10, "export": 5}
    precursors = set(net.extras["precursors"])
    assert all(p["target"] in precursors for p in net.extras["proteins"] if p["type"] == "import")
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    for atom in "123":
        assert not np.any(molecule_vector(net, ids, S, atom))
    for r in net.reactions:
        if r.rate:
            assert r.rate["k"] == pytest.approx(r.rate["affinity"] / r.rate["kout_product"])
            assert r.catalysts


def test_invalid_parameters():
    with pytest.raises(ValueError, match="n_precursors"):
        generate_network("bnc-cell", max_length=6, n_precursors=53)
    with pytest.raises(ValueError, match="different order"):
        generate_network("bnc-cell", mode="cell", proteins=["enzyme 1-2|1 + 2|2"])
    with pytest.raises(ValueError, match="own reactants"):
        generate_network("bnc-cell", mode="cell", proteins=["enzyme 1-2|1 + 1-2|1"])
    with pytest.raises(ValueError, match="not a valid molecule"):
        generate_network("bnc-cell", mode="cell", proteins=["import 1-3=1"])
    with pytest.raises(ValueError, match="max_length"):
        generate_network("bnc-cell", mode="cell", proteins=["enzyme 1|2-2-2-2-1 + 1-2-2-2-2|1"])
