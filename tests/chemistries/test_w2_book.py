"""W2 random-topology chemistries specified by the book: published equations and model rules."""

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import generate_network
from chemart.chemistries.farmer_immune import match
from chemart.chemistries.jain_krishna import attractor


def derivatives(net, state):
    ids, f = rhs(net)
    return dict(zip(ids, f(0.0, np.array([state[s] for s in ids]))))


# --- random catalytic networks (book 7.2.9) ---------------------------------
def test_random_catalytic_network_equation_and_no_direct_replication():
    net = generate_network("random-catalytic-networks", seed=4, n=6, density=0.3)
    for r in net.reactions:
        (made,) = [s for s in r.products if r.products[s] > r.reactants.get(s, 0)]
        assert made not in r.reactants
    x = dict(zip([s.id for s in net.species], np.random.default_rng(0).dirichlet(np.ones(6))))
    production = {s: 0.0 for s in x}
    for r in net.reactions:
        flux = r.rate["k"] * np.prod([x[s] ** n for s, n in r.reactants.items()])
        for s in r.products:
            production[s] += flux * (r.products[s] - r.reactants.get(s, 0))
    phi = sum(production.values())
    ours = derivatives(net, x)
    for s in x:  # eq. 7.32 with sum x = 1
        assert ours[s] == pytest.approx(production[s] - x[s] * phi, abs=1e-12)


# --- Jain-Krishna (book 15.2.2) ----------------------------------------------
def test_jain_krishna_attractor_is_a_fixed_point_of_eq_15_7():
    C = np.zeros((4, 4), dtype=int)
    C[0, 1] = C[1, 0] = C[2, 0] = 1       # 2-cycle X1 <-> X2 feeding X3; X4 isolated
    lam, x = attractor(C)
    assert lam == pytest.approx(1.0)
    growth = C @ x
    assert np.allclose(growth - x * growth.sum(), 0.0, atol=1e-12)
    assert x[3] == pytest.approx(0.0, abs=1e-12)


def test_jain_krishna_links_follow_eq_15_7_and_acs_emerges():
    net = generate_network("jain-krishna", seed=2, m=30, p=0.02)
    assert all(len(r.reactants) == 1 and set(r.reactants) <= set(r.products) for r in net.reactions)
    evolved = generate_network("jain-krishna", seed=2, m=30, p=0.02, graph_updates=1500)
    assert evolved.extras["analysis"]["perron_frobenius_eigenvalue"] >= 1.0 - 1e-9


# --- Kauffman binary polymer model (book 6.3.1) ------------------------------
def test_kauffman_polymer_counts_and_conservation():
    net = generate_network("kauffman-autocatalytic-sets", P=0.0)
    assert len(net.species) == 2 + 4 + 8 + 16 + 32
    assert len(net.reactions) == 2 * sum((n - 1) * 2 ** n for n in range(2, 6))
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    for law in net.extras["conservation"]:
        assert not np.any(np.array([law["vector"].get(s, 0) for s in ids]) @ S)


def test_kauffman_catalysis_probability():
    everything = generate_network("kauffman-autocatalytic-sets", P=1.0, max_length=3)
    none = generate_network("kauffman-autocatalytic-sets", P=0.0, max_length=3)
    pairs = len(none.reactions) // 2
    assert len(everything.reactions) == len(none.reactions) + 2 * pairs * len(none.species)


# --- Farmer immune network (book 11.2.2) -------------------------------------
def test_match_counts_complementary_alignments():
    assert match("0000", "1111", s=4) == 1            # only the full overlap has 4 complementary bits
    assert match("0000", "1111", s=1) == 1 + 2 + 3 + 4 + 3 + 2 + 1
    assert match("1111", "1111", s=1) == 0


def test_immune_network_reproduces_eqs_11_5_and_11_6():
    p = dict(N=5, M=2, l_e=6, l_p=6, s=4, c=0.7, k1=1.3, k2=0.4, k3=2.0)
    net = generate_network("farmer-immune", seed=5, **p)
    struct = {s.id: dict(part.split("=") for part in s.structure.split()) for s in net.species}
    X = [f"X{i + 1}" for i in range(p["N"])]
    Y = [f"Y{j + 1}" for j in range(p["M"])]
    m = [[match(struct[a]["e"], struct[b]["p"], p["s"]) for b in X] for a in X]
    my = [[match(struct[a]["e"], struct[b]["p"], p["s"]) for b in X] for a in Y]
    rng = np.random.default_rng(1)
    state = {s: v for s, v in zip(X + Y, rng.uniform(0.1, 2.0, len(X) + len(Y)))}
    ours = derivatives(net, state)
    x = [state[s] for s in X]
    y = [state[s] for s in Y]
    for i in range(p["N"]):
        paper = p["c"] * (
            sum(m[j][i] * x[i] * x[j] for j in range(p["N"]))
            - p["k1"] * sum(m[i][j] * x[i] * x[j] for j in range(p["N"]))
            + sum(my[j][i] * x[i] * y[j] for j in range(p["M"]))
        ) - p["k2"] * x[i]
        assert ours[X[i]] == pytest.approx(paper, rel=1e-9, abs=1e-12)
    for j in range(p["M"]):
        paper = -p["k3"] * sum(my[j][i] * x[i] for i in range(p["N"])) * y[j]
        assert ours[Y[j]] == pytest.approx(paper, rel=1e-9, abs=1e-12)


# --- ACGP (book 16.6) ---------------------------------------------------------
def test_acgp_program_is_the_reaction_multiset():
    net = generate_network("acgp", seed=3, n_instructions=40)
    assert len(net.reactions) == len(net.extras["program"]) == 40
    inputs = set(net.extras["input_registers"])
    for r, instruction in zip(net.reactions, net.extras["program"]):
        (dst,) = r.products
        assert dst not in inputs and instruction.startswith(f"{dst} = ")
    with pytest.raises(ValueError, match="read-only"):
        generate_network("acgp", output_registers=[0])
