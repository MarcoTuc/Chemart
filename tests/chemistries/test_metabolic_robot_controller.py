"""Reaction graphs of Ziegler & Banzhaf (2001), with catalysis as in Ziegler, Dittrich & Banzhaf (1998)."""

import numpy as np
import pytest
from odes import rhs

from chemart import generate_network
from chemart.chemistries.metabolic_robot_controller import TYPES, material_balance


def test_paper_example_violates_material_balance():
    # Sec. 3.1.2: A -> B + C, C -> A. Rows are reactions, columns A, B, C.
    M = np.array([[-1, 1, 1], [1, 0, -1]])
    assert np.linalg.matrix_rank(M) == 2 and M.shape[1] - np.linalg.matrix_rank(M) == 1   # eq. 23
    assert material_balance(M) is None                  # the only solution has m_B = 0


def test_elementary_reaction_is_balanced():
    # Eqs. 18-20: A + B -> C + D needs m_A + m_B = m_C + m_D, not m_A = m_C.
    m = material_balance(np.array([[-1, -1, 1, 1]]))
    assert m is not None and np.all(m > 0)
    assert m[0] + m[1] == pytest.approx(m[2] + m[3])


@pytest.mark.parametrize("seed", range(5))
def test_generated_graphs_fulfil_material_balance(seed):
    net = generate_network("metabolic-robot-controller", seed=seed, n_substances=12, n_reactions=15)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    (law,) = net.extras["conservation"]
    m = np.array([law["vector"][s] for s in ids], dtype=float)
    assert np.all(m > 0)
    assert np.allclose(m @ S, 0)
    graph = net.extras["reaction_graph"]
    assert len(graph) == 15
    for node in graph:                                   # every node is one of eqs. 4-7
        assert (len(node["reactants"]), len(node["products"])) == TYPES[node["type"]]
        assert sorted(node["reactants"]) != sorted(node["products"])


def test_input_and_output_sets():
    net = generate_network("metabolic-robot-controller", seed=3)
    # Table 1 terminal set {a, b, c, d, e}; sec. 5.1-5.2
    assert net.extras["input_set"] == ["b", "c", "d", "e"]
    assert net.extras["output_set"] == ["a"]
    assert {"a", "b", "c", "d", "e"} <= {s.id for s in net.species}
    assert net.outflow == {"a": 0.9}                     # Table 2: alpha = 0.9
    assert net.extras["actuators"]["a"]["threshold"] == 0.1   # Table 2: a_min = 0.1


def test_sensor_inflow_follows_eq_30():
    # Surrounded by obstacles, the inflow sums to 0.4 of the cell volume (sec. 5.1).
    net = generate_network("metabolic-robot-controller", seed=0)
    assert sum(net.inflow.values()) == pytest.approx(0.4)
    # f(s1, s2) = maxInflow * max(s1, s2) / 1023, linear in the larger reading
    net = generate_network("metabolic-robot-controller", seed=0,
                           sensor_readings={"front1": 200, "front2": 1023 / 2, "left2": 1023})
    assert net.inflow == pytest.approx({"b": 0.1, "c": 0.05, "d": 0.0, "e": 0.0})


def test_catalysed_rate_matches_the_1998_appendix():
    # w1 + w5 = k [A][B] + k kappa [A][B][C]: the catalyst multiplies the rate by (1 + kappa [C]).
    net = generate_network("metabolic-robot-controller", seed=11, n_substances=6, n_reactions=1,
                           reaction_types=[5], p_modifier=1.0, p_inhibitor=0.0, k_0=1.5, k_min=2.0, k_max=8.0)
    (node,) = net.extras["reaction_graph"]
    catalyst, kappa = node["modifier"]["species"], node["modifier"]["k"]
    assert node["modifier"]["effect"] == "catalytic" and 2.0 <= kappa <= 8.0
    assert len(net.reactions) == 2 and net.reactions[1].catalysts.get(catalyst)

    ids, f = rhs(net)
    rng = np.random.default_rng(0)
    x = rng.uniform(0.2, 2.0, len(ids))
    c = dict(zip(ids, x))
    flux = 1.5 * np.prod([c[s] for s in node["reactants"]]) * (1 + kappa * c[catalyst])
    expected = dict.fromkeys(ids, 0.0)
    for s in node["reactants"]:
        expected[s] -= flux
    for s in node["products"]:
        expected[s] += flux
    # odes.rhs includes the network's constant sensor inflow and first-order outflow
    for s, value in (net.inflow or {}).items():
        expected[s] += value
    for s in ids:
        rate = net.outflow.get(s, 0.0) if isinstance(net.outflow, dict) else (net.outflow or 0.0)
        expected[s] -= rate * c[s]
    assert dict(zip(ids, f(0.0, x))) == pytest.approx(expected, rel=1e-9, abs=1e-12)


def test_inhibition_is_recorded_on_the_rate():
    net = generate_network("metabolic-robot-controller", seed=2, n_reactions=4, p_modifier=1.0, p_inhibitor=1.0)
    assert len(net.reactions) == 4                       # no extra channel for inhibitors
    for reaction, node in zip(net.reactions, net.extras["reaction_graph"]):
        assert reaction.rate["inhibitor"] == node["modifier"]["species"]
        assert 10.0 <= reaction.rate["k_inhibition"] <= 10.0


def test_invalid_combinations_are_rejected():
    with pytest.raises(ValueError, match="n_substances"):
        generate_network("metabolic-robot-controller", n_substances=4)
    with pytest.raises(ValueError, match="disjoint"):
        generate_network("metabolic-robot-controller", actuator_map={"b": "rotate"})
    with pytest.raises(ValueError, match="k_min"):
        generate_network("metabolic-robot-controller", k_min=5.0, k_max=1.0)
    with pytest.raises(ValueError, match="sensor_readings"):
        generate_network("metabolic-robot-controller", sensor_readings={"nose": 3})
    with pytest.raises(ValueError, match="reaction_types"):
        generate_network("metabolic-robot-controller", reaction_types=[3])
