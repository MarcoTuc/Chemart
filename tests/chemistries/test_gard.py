"""GARD: rate equation and beta distribution of the 2018 PLoS ONE parameterisation."""

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import generate_network


def beta_matrix(net, G):
    beta = np.zeros((G, G))
    for r in net.reactions:
        if len(r.reactants) == 2 and all(s.startswith(("L", "A")) for s in r.reactants):
            (li,) = [s for s in r.reactants if s.startswith("L")]
            (aj,) = [s for s in r.reactants if s.startswith("A")]
            beta[int(li[1:]) - 1, int(aj[1:]) - 1] = r.rate["k"]
    return beta


def test_log_beta_is_normal_minus_four_four():
    beta = beta_matrix(generate_network("gard", seed=3), 100)
    log_beta = np.log(beta).ravel()
    assert log_beta.mean() == pytest.approx(-4.0, abs=0.2)
    assert log_beta.std() == pytest.approx(4.0, abs=0.2)


def test_rate_equation_without_crowding():
    G = 5
    net = generate_network("gard", seed=1, N_G=G)
    beta = beta_matrix(net, G)
    ids, f = rhs(net)
    n = np.random.default_rng(2).uniform(0, 20, G)
    state = {**{f"L{i + 1}": 0.01 for i in range(G)}, **{f"A{i + 1}": n[i] for i in range(G)}}
    ours = dict(zip(ids, f(0.0, np.array([state[s] for s in ids]))))
    for i in range(G):
        paper = 0.01 * 0.01 + (beta[i] @ n) * 0.01 - 1e-4 * n[i]
        assert ours[f"A{i + 1}"] == pytest.approx(paper, rel=1e-9)
    assert all(r.rate.get("crowding_capacity") == 100.0 for r in net.reactions if "L1" in r.reactants)


def test_catalysed_leaving_adds_reverse_catalysis():
    base = generate_network("gard", seed=1, N_G=4)
    both = generate_network("gard", seed=1, N_G=4, catalysed_leaving=True)
    assert len(both.reactions) == len(base.reactions) + 16
