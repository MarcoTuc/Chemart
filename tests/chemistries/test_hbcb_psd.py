"""HBCB / PSD (book 18.3.2, Oohashi et al. 2009, ref [643]).

The 2009 paper could not be obtained (see the catalog decisions), so the
published content tested here is the book's paragraph and the paper's
abstract: the hierarchy ordered by covalent bond energy, hydrolysis of a BP
into BMs, matter that is finite and recycled, and the superiority of the
species whose self-decomposition stops at reusable biomonomers.
"""

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.hbcb_psd import class_names, rebuild_cost

ID = "hbcb-psd"


def texts(net):
    return set(net.to_text().splitlines())


# --- the HBCB hierarchy -----------------------------------------------------
def test_hierarchy_is_the_closure_of_cleavage_from_the_biopolymers():
    """|S| = n_families x n_levels, named BP, BM and the lower classes."""
    net = generate_network(ID)
    assert net.status == "complete"
    assert len(net.species) == 4 * 4
    assert class_names(4) == ["C0", "C1", "BM", "BP"]
    assert [s.id for s in net.species][:4] == ["BP_0", "BM_0", "C1_0", "C0_0"]
    # every class is units_per_class molecules of the class below it
    assert "4 x BM_0" in dict((s.id, s.structure) for s in net.species)["BP_0"]
    # a BP holds 4^3 = 64 base units
    assert "64 base units" in dict((s.id, s.structure) for s in net.species)["BP_0"]
    assert net.initial_state == {f"BP_{f}": 8.0 for f in range(4)}


def test_hydrolysis_and_synthesis_connect_neighbouring_classes():
    """Book 18.3.2 / abstract: hydrolysis decomposes a BP into BMs; re-use is its reverse."""
    net = generate_network(ID, n_families=1)
    assert texts(net) == {
        "BP_0 -> 4 BM_0", "4 BM_0 -> BP_0",
        "BM_0 -> 4 C1_0", "4 C1_0 -> BM_0",
        "C1_0 -> 4 C0_0", "4 C0_0 -> C1_0",
    }
    assert all(r.rate is None for r in net.reactions)   # no published kinetics


def test_bond_energies_are_ordered_down_the_hierarchy():
    """HBCB: the BM-BM links of a BP are the weakest bonds, deeper bonds are stronger."""
    net = generate_network(ID, n_families=1)
    per_link = net.extras["hbcb"]["bond_energy_per_link"]
    assert per_link == {"C1_0": 4.0, "BM_0": 2.0, "BP_0": 1.0}
    # a molecule's total bond energy accumulates down the hierarchy
    energies = net.extras["energies"]
    assert energies["C0_0"] == 0.0
    assert energies["C1_0"] == 3 * 4.0
    assert energies["BM_0"] == 4 * (3 * 4.0) + 3 * 2.0
    assert energies["BP_0"] == 4 * energies["BM_0"] + 3 * 1.0
    assert energies["BP_0"] > energies["BM_0"] > energies["C1_0"] > energies["C0_0"]


# --- finite, recycled matter ------------------------------------------------
def test_base_units_are_conserved_exactly():
    """The ecosystem is closed: every reaction conserves each family's base units."""
    net = generate_network(ID)
    assert net.inflow is None and net.outflow is None
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    laws = net.extras["conservation"]
    assert len(laws) == 4
    for f, law in enumerate(laws):
        assert law["vector"] == {f"BP_{f}": 64, f"BM_{f}": 16, f"C1_{f}": 4, f"C0_{f}": 1}
        m = np.array([law["vector"].get(i, 0) for i in ids], dtype=float)
        assert np.allclose(S.T @ m, 0.0)


# --- programmed self-decomposition -----------------------------------------
@pytest.mark.parametrize(
    "target, reaction, recyclable",
    [
        ("BM", "BP_0 -> 4 BM_0", True),          # hydrolysis to reusable monomers
        ("sub-BM", "BP_0 -> 16 C1_0", False),
        ("base", "BP_0 -> 64 C0_0", False),
    ],
)
def test_psd_cleaves_a_biopolymer_down_to_its_target_class(target, reaction, recyclable):
    net = generate_network(ID, psd_target=target, n_families=1)
    assert net.extras["psd"]["reactions"] == [reaction]
    assert net.extras["analysis"]["recyclable"] is recyclable
    assert reaction in texts(net)
    # PSD to BM *is* the hydrolysis step, not a second copy of it
    assert len(net.reactions) == (6 if target == "BM" else 7)


def test_a_species_without_self_decomposition_has_no_psd_reaction():
    net = generate_network(ID, psd_target="none", n_families=1)
    assert net.extras["psd"]["reactions"] == []
    assert net.extras["analysis"]["psd_target_class"] is None
    assert net.extras["analysis"]["recyclable"] is False
    assert len(net.reactions) == 6      # only the chemistry itself


def test_decomposition_to_bm_is_the_cheapest_route_back_to_a_biopolymer():
    """The paper's result: cleaving to BMs beats cleaving to classes lower than BM.

    Abstract: "the virtual species in which the self-decomposition process mainly
    involved covalent bond cleavage from a BP to BMs showed evolutionary superiority
    over other species in which the self-decomposition process involved cleavage from
    BP to classes lower than BM". Here that superiority is the re-use cost.
    """
    cost = generate_network(ID).extras["analysis"]["rebuild_cost_by_class"]
    assert set(cost) == {"BM", "C1", "C0"}
    bm = cost["BM"]
    assert bm == {"syntheses_per_bp": 1, "bond_energy_per_bp": 3.0}
    for lower in ("C1", "C0"):
        assert cost[lower]["syntheses_per_bp"] > bm["syntheses_per_bp"]
        assert cost[lower]["bond_energy_per_bp"] > bm["bond_energy_per_bp"]
    # and strictly worse the deeper the cleavage goes
    assert cost["C0"]["syntheses_per_bp"] > cost["C1"]["syntheses_per_bp"]
    assert cost["C0"]["bond_energy_per_bp"] > cost["C1"]["bond_energy_per_bp"]
    # 4 C1 -> BM sixteen times, then 4 BM -> BP: 21 syntheses per biopolymer
    assert cost["C0"]["syntheses_per_bp"] == 16 + 4 + 1


@pytest.mark.parametrize("n_levels", [2, 3, 4, 5, 6])
def test_bm_is_optimal_for_every_hierarchy(n_levels):
    energies = [float(n_levels - i) for i in range(n_levels - 1)]
    costs = [rebuild_cost(level, n_levels, 3, energies) for level in range(n_levels - 1)]
    assert costs[-1] == (1, 2 * energies[-1])          # BM: one synthesis
    assert costs == sorted(costs, reverse=True)        # strictly worse further down
    assert len({c for c in costs}) == len(costs)


def test_deeper_decomposition_is_never_cheaper_in_the_generated_network():
    for n_levels, energies in [(3, [2.0, 1.0]), (5, [8.0, 4.0, 2.0, 1.0])]:
        cost = generate_network(ID, n_levels=n_levels, bond_energies=energies,
                                n_families=1).extras["analysis"]["rebuild_cost_by_class"]
        best = min(cost.values(), key=lambda c: c["syntheses_per_bp"])
        assert cost["BM"] == best


# --- parameters -------------------------------------------------------------
def test_the_network_does_not_depend_on_the_seed():
    assert generate_network(ID, seed=1).to_dict()["reactions"] == generate_network(ID, seed=2).to_dict()["reactions"]


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(n_levels=3), "bond_energies must have 2 entries"),
        (dict(bond_energies=[1.0, 2.0, 3.0]), "decrease strictly"),
        (dict(bond_energies=[4.0, 2.0, 2.0]), "decrease strictly"),
        (dict(bond_energies=[4.0, 2.0, 0.0]), "positive"),
        (dict(n_levels=2, bond_energies=[1.0], psd_target="sub-BM"), "n_levels >= 3"),
        (dict(n_levels=1), "n_levels"),
    ],
)
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)


def test_a_two_class_hierarchy_is_just_bp_and_bm():
    net = generate_network(ID, n_levels=2, bond_energies=[1.0], n_families=1)
    assert [s.id for s in net.species] == ["BP_0", "BM_0"]
    assert texts(net) == {"BP_0 -> 4 BM_0", "4 BM_0 -> BP_0"}
    assert net.extras["analysis"]["rebuild_cost_by_class"] == {
        "BM": {"syntheses_per_bp": 1, "bond_energy_per_bp": 3.0}
    }
