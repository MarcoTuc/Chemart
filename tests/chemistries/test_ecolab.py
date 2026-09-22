"""Ecolab (Standish): eq. 1 of arXiv nlin/0404011 and the mutation operator of the EcoLab
technical report / ecolab_model.cc (boundedness, caps, speciation rule, rounding)."""

import numpy as np
import pytest
from chemart.simulate import integrate, rhs

from chemart import generate_network
from chemart.chemistries import ecolab


def beta_of(net):
    eco = net.extras["ecosystem"]
    index = {s: i for i, s in enumerate(eco["species"])}
    B = np.diag(eco["beta_diag"])
    for a, b, v in eco["beta_offdiag"]:
        B[index[a], index[b]] = v
    return eco["species"], np.array(eco["r"]), B, np.array(eco["mu"])


def test_network_rate_equations_are_eq_1():
    # nlin/0404011 eq. 1: n_i' = r_i n_i - n_i sum_j beta_ij n_j
    net = generate_network("ecolab", seed=3)
    X, r, B, _ = beta_of(net)
    assert len(X) > 3 and net.status == "complete"
    ids, f = rhs(net)
    assert ids == X
    rng = np.random.default_rng(0)
    for _ in range(3):
        n = rng.uniform(0, 200, len(X))
        assert np.allclose(f(0.0, n), r * n - n * (B @ n), rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_mutation_preserves_boundedness_and_caps(seed):
    # technical report, mutation steps 4-6: beta_ii > 0 with |beta_ii| >= |r|/(0.1 INT_MAX),
    # beta_ij + beta_ji >= 0 (report's sufficient criterion), mu <= mut_max
    net = generate_network("ecolab", seed=seed, mut_max=0.02, cycles=60, gen_bias=0.5)
    X, r, B, mu = beta_of(net)
    assert net.extras["analysis"]["speciations"] > 20
    assert np.all(np.diag(B) > 0)
    assert np.all(np.diag(B) >= np.abs(r) / (0.1 * ecolab.INT_MAX))
    assert np.all(B + B.T >= -1e-15)
    assert np.all(mu <= 0.02) and np.all(mu > 0)
    # hence sum_i n_i' = r.n - n.beta.n < 0 for every large enough positive n
    rng = np.random.default_rng(seed)
    for _ in range(50):
        n = rng.exponential(1.0, len(X)) * 1e6
        assert n @ B @ n > 0 and r @ n - n @ B @ n < 0


def test_populations_stay_bounded_under_the_ode():
    # nlin/0404011: boundedness ensures populations never exceed a limit
    net = generate_network("ecolab", seed=5, nsp=8, cycles=30)
    X, r, B, _ = beta_of(net)
    start = {s: 1e4 for s in X}
    x, _ = integrate(net, 2000, x0=start)
    totals = np.sum([x[s] for s in X], axis=0)
    assert totals[-1] < totals[0] and totals.max() <= totals[0] * (1 + 1e-6)


def test_only_growing_species_speciate_and_offspring_resemble_parents():
    # nlin/0404011: mutants proportional to n_i r_i mu_i, so r_i <= 0 never speciates;
    # "offspring resemble their parents" (book 8.2.3)
    net = generate_network("ecolab", seed=11, cycles=60)
    phylo = net.extras["analysis"]["phylogeny"]
    children = [v for v in phylo.values() if v["parent"] is not None]
    assert len(children) > 20
    assert all(phylo[c["parent"]]["r"] > 0 for c in children)
    jumps = [abs(c["r"] - phylo[c["parent"]]["r"]) for c in children]
    assert np.median(jumps) < 0.05 * 0.2   # well below the r range repro_max - repro_min


def test_no_mutation_no_new_species_and_diversity_bookkeeping():
    # Standish 2000: with mutation switched off diversity can only fall to a stable ecology
    off = generate_network("ecolab", seed=4, mut_max=0.0, cycles=50)
    a = off.extras["analysis"]
    assert a["speciations"] == 0 and all(v["parent"] is None for v in a["phylogeny"].values())
    assert np.all(np.diff([20] + a["history"]["diversity"]) <= 0)
    # nlin/0404011 eq. 3: diversity = speciations - extinctions (plus the founders)
    on = generate_network("ecolab", seed=4, cycles=40).extras["analysis"]
    h = on["history"]
    assert h["diversity"] == list(20 + np.cumsum(h["speciations"]) - np.cumsum(h["extinctions"]))
    assert on["diversity"] == h["diversity"][-1] == len(generate_network("ecolab", seed=4, cycles=40).species)


def test_probabilistic_rounding():
    # technical report: "0.25 is rounded up to 1 25% of the time and down to 0 75% of the time"
    rng = np.random.default_rng(0)
    out = ecolab._round(rng, np.full(200_000, 0.25))
    assert set(np.unique(out)) == {0, 1}
    assert out.mean() == pytest.approx(0.25, abs=0.005)
    assert list(ecolab._round(rng, np.array([-3.0, 0.0, 7.0]))) == [0, 0, 7]


def test_zero_cycles_is_a_bounded_random_ecosystem():
    net = generate_network("ecolab", seed=2, cycles=0, nsp=30)
    X, r, B, mu = beta_of(net)
    assert len(X) == 30 and np.allclose(np.diag(B), 0.001) and np.all(mu == 0.01)
    assert np.all(B + B.T >= 0) and np.all((r >= -0.1) & (r <= 0.1))
    assert all(v == 100.0 for v in net.initial_state.values())


def test_invalid_combinations_are_rejected():
    with pytest.raises(ValueError, match="conn"):
        generate_network("ecolab", nsp=3, conn=3)
    with pytest.raises(ValueError, match="repro_min"):
        generate_network("ecolab", repro_min=0.2)
