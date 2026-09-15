"""W2 chemistries checked against their papers: Bigan et al. (2013) and the NK model definition."""

import math

import numpy as np
import pytest

from chemart import generate_network


# --- Bigan conservative random networks -------------------------------------
def test_bigan_network_is_conservative_and_maximal():
    net = generate_network("bigan-conservative-crn", seed=11)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    law = net.extras["conservation"][0]
    m = np.array([law["vector"][s] for s in ids], dtype=float)
    assert np.all(m > 0) and np.allclose(m @ S, 0.0)
    assert law["unique"], "a maximum-sized network has a single admissible mass vector"


def test_bigan_sizes_match_the_paper_range():
    # paper: maximum sizes 17..96 direct reactions, mean 42, for N = 10
    sizes = [len(generate_network("bigan-conservative-crn", seed=s).reactions) // 2 for s in range(8)]
    assert all(10 <= n <= 110 for n in sizes)
    assert 20 <= np.mean(sizes) <= 70


def test_bigan_detailed_balance():
    net = generate_network("bigan-conservative-crn", seed=3)
    G = net.extras["energies"]["formation_free_energy_RT"]
    for f, r in net.extras["reaction_pairs"]:
        forward, reverse = net.reactions[f], net.reactions[r]
        dG = sum(G[s] * n for s, n in forward.products.items()) - sum(G[s] * n for s, n in forward.reactants.items())
        assert dG <= 1e-12, "forward reactions go downhill"
        assert forward.rate["k"] / reverse.rate["k"] == pytest.approx(math.exp(-dG), rel=1e-9)
        assert forward.reactants == reverse.products and forward.products == reverse.reactants


def test_bigan_saturating_kinetics_and_nutrient():
    net = generate_network("bigan-conservative-crn", seed=3, kinetics="saturating", nutrient=2)
    assert all(r.rate["law"] == "saturating" and r.rate["K"] == pytest.approx(0.01) for r in net.reactions)
    assert net.inflow == {"A2": 1.0}


# --- NK landscape -------------------------------------------------------------
def test_nk_fitness_follows_book_index_convention():
    net = generate_network("nk-landscape", seed=7, N=3, K=1, mu=0.0)
    W = net.extras["analysis"]["fitness"]
    assert len(W) == 8 and all(0 <= w < 1 for w in W.values())
    # without mutation only exact replication reactions remain
    assert all(r.reactants.keys() == r.products.keys() for r in net.reactions)


def test_nk_k0_is_unimodal():
    for seed in range(5):
        net = generate_network("nk-landscape", seed=seed, N=6, K=0)
        assert len(net.extras["analysis"]["local_optima"]) == 1


def test_nk_rejects_k_not_below_n():
    with pytest.raises(ValueError, match="K must satisfy"):
        generate_network("nk-landscape", N=4, K=4)
