"""W1 written-down networks: structure published in the book (reaction lists, invariants)."""

import numpy as np
import pytest

from chemart import generate_network


def conservation_holds(net):
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    for law in net.extras["conservation"]:
        m = np.array([law["vector"].get(s, 0) for s in ids])
        residual = m @ S
        if "modulus" in law:
            residual = residual % law["modulus"]
        if np.any(residual != 0):
            return False
    return True


def test_chameleon_counts_and_invariants():
    net = generate_network("chameleon")
    assert sum(net.initial_state.values()) == 5400
    assert net.initial_state == {"r": 2700, "g": 0, "b": 2700}
    assert conservation_holds(net)


@pytest.mark.parametrize("network", ["cno-cycle", "pp-chain", "big-bang-nucleosynthesis"])
def test_nuclear_networks_conserve_baryons_charge_and_leptons(network):
    net = generate_network("nuclear-reaction-networks", network=network)
    assert [law["name"] for law in net.extras["conservation"]][:2] == ["baryon number", "charge"]
    assert conservation_holds(net)


def test_cno_cycle_is_book_eq_20_4():
    net = generate_network("nuclear-reaction-networks")
    assert [r.to_text() for r in net.reactions] == [
        "12C + 1H -> 13N + gamma",
        "13N -> 13C + e+ + nu_e",
        "13C + 1H -> 14N + gamma",
        "14N + 1H -> 15O + gamma",
        "15O -> 15N + e+ + nu_e",
        "15N + 1H -> 12C + 4He",
    ]


def test_mechanical_self_assembly_conserves_monomers():
    net = generate_network("mechanical-self-assembly")
    assert len(net.reactions) == 9 and conservation_holds(net)


def test_n_economy_default_is_book_eq_20_6():
    net = generate_network("n-economy")
    assert [r.to_text() for r in net.reactions] == [
        "2 2 -> 2 5",
        "2 2 + 2 5 + 2 13 -> 3 130",
        "3 2 + 3 130 -> 6 260",
        "6 2 + 6 260 -> 6 104",
        "6 104 -> 13 2",
    ]
    assert {s.id: s.structure for s in net.species}["260"] == "2^2*5*13"


def test_xor_gate_follows_truth_table():
    net = generate_network("organization-computing", a=True, b=False)
    lines = net.to_text().splitlines()
    for a in (0, 1):
        for b in (0, 1):
            assert f"a{a} + b{b} -> c{a ^ b}" in lines
    assert "c0 + c1 -> ∅" in lines
    assert net.initial_state == {"a1": 1, "b0": 1}


def test_maximal_independent_set_rules():
    net = generate_network("organization-computing", problem="maximal-independent-set", graph=[[0, 1], [1, 2]])
    lines = set(net.to_text().splitlines())
    assert "s0_0 + s0_2 -> 2 s1_1" in lines          # eq. 17.6 for the middle vertex
    assert {"s1_1 -> s0_0", "s1_0 -> s0_1"} <= lines  # eq. 17.7, both directions
    assert "s0_1 + s1_1 -> ∅" in lines               # eq. 17.5


def test_quasispecies_mutation_rows_sum_to_fitness():
    net = generate_network("quasispecies", L=4, m=0.5)
    out = {}
    for r in net.reactions:
        (source,) = [s for s in r.reactants]
        out[source] = out.get(source, 0.0) + r.rate["k"]
    for genome, total in out.items():
        assert total == pytest.approx(genome.count("1") / 4)
    no_mutation = generate_network("quasispecies", L=4, m=0.0)
    assert all(len(r.products) == 1 for r in no_mutation.reactions)


def test_oregonator_real_f_keeps_rate_equations():
    net = generate_network("oregonator", f=1.5, k5=2.0)
    regen = [r for r in net.reactions if r.reactants == {"Z": 1}]
    assert sum(r.rate["k"] for r in regen) == pytest.approx(2.0)
    assert sum(r.rate["k"] * r.products.get("Y", 0) for r in regen) == pytest.approx(2.0 * 1.5)


def test_repressilator_is_cyclic():
    lines = generate_network("repressilator").to_text().splitlines()
    assert any(line.startswith("G1 + 2 P3 -> C1") for line in lines)
    assert any(line.startswith("G2 + 2 P1 -> C2") for line in lines)


def test_french_flag_three_segments():
    net = generate_network("french-flag", n_cells=9)
    fates = [next(iter(r.products)).split("_")[0] for r in net.reactions]
    assert fates == ["blue"] * 3 + ["white"] * 3 + ["red"] * 3
