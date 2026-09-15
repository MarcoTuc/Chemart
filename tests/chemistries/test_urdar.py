"""Urdar: Gerlee & Lundh, Evolution 64:2716 (2010) (book 8.2.3, ref [314]).

Model: eq. 1, "Implementation", appendices A-B. Results: figs. 4 and 6 and
the metabolic-depth numbers (flow dependence), fig. 7 (129/145 coexistence),
and the invasion-experiment text (antidiagonal neutrality, 105 over 109,
c(150, 109) = 1, the intransitive triple 126/134/141).
"""

from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import urdar as U

ID = "urdar"


def bits(s):
    return np.array([int(c) for c in s], dtype=np.uint8)


def run(rules, abundances, *, n=256, flow=0.03, mu=0.0, updates=400, seed=0):
    p = SimpleNamespace(rules=rules, abundances=abundances, n_organisms=n)
    rng = np.random.default_rng(seed)
    return U.simulate(rng, U.initial_population(p), 5 * n, 100, flow, mu, 0.1, 0.95,
                      updates, record=False, gain_only=True)


def c(i, j, seeds=(0,)):
    """Invasion matrix entry: final frequency of rule i started at 9:1 against j (mu = 0, gamma = 0.03)."""
    return float(np.mean([np.mean(run([i, j], [9, 1], seed=s)["population"] == i) for s in seeds]))


def test_elementary_ca_rule_30_of_figure_a1():
    assert U.rule_table(30) == "00011110"
    rows = ["00000100000", "00001110000", "00011001000", "00110111100"]
    x = bits(rows[0])
    for expected in rows[1:]:
        x = U.ca_step(x, 30)
        assert "".join(map(str, x)) == expected


def test_entropy_of_appendix_b_and_inflow_energy():
    # s = 0 for ...010101..., s -> 1 bit for random strings
    assert U.entropy(bits("01" * 50))[0] == pytest.approx(0.0, abs=1e-12)
    assert U.entropy(bits("0" * 100))[0] == 0.0
    rng = np.random.default_rng(0)
    assert U.entropy(rng.integers(0, 2, (5, 4000))).min() > 0.99
    # eq. 2 with p0 = 0.95: E0 = 1 - s0 = 0.714; the paper quotes E0 of about 0.8
    e0 = U.energy(U.fresh_strings(rng, 2000, 100, 0.95)).mean()
    assert 0.7 < e0 < 0.8
    # the inflow is balanced between strings dominated by ones and by zeros
    ones = U.fresh_strings(rng, 2000, 100, 0.95).mean(axis=1) > 0.5
    assert 0.45 < ones.mean() < 0.55


def test_reproduction_probability_eq_1():
    beta = 0.1
    dE = np.array([-0.2, 0.0, 0.001, 0.005])
    expected = [0.0, 0.0, *((1 - np.exp(-dE[2:] / beta)) / (1 - np.exp(-beta)))]
    assert np.allclose(U.reproduction_probability(dE, beta), expected)
    # beta -> 0: any gain gives (at least) certain reproduction; larger beta: the gain matters
    assert U.reproduction_probability(0.01, 1e-3) >= 1
    assert U.reproduction_probability(0.01, 2.0) < U.reproduction_probability(0.1, 2.0) < 1


def test_antidiagonal_rules_make_complementary_strings_of_equal_entropy():
    # "rule 145 and its antidiagonal partner 255 - 145 = 110 ... output strings are inverses"
    rng = np.random.default_rng(3)
    x = U.fresh_strings(rng, 50, 100, 0.8)
    for k in (145, 30, 86, 150):
        a, b = U.ca_step(x, k), U.ca_step(x, 255 - k)
        assert np.array_equal(a, 1 - b)
        assert np.allclose(U.entropy(a), U.entropy(b))


def test_observed_network_events_follow_the_model():
    net = generate_network(ID, seed=2)
    assert net.status == "observed"
    struct = {s.id: s.structure for s in net.species}
    energies = net.extras["energies"]
    births = 0
    for r in net.reactions:
        orgs_in = Counter({s: n for s, n in r.reactants.items() if s.startswith("R")})
        orgs_out = Counter({s: n for s, n in r.products.items() if s.startswith("R")})
        mets_in = [s for s in r.reactants if s.startswith("m")]
        mets_out = [s for s in r.products if s.startswith("m")]
        if not orgs_in:                                   # flow: m -> ∅ or ∅ -> m
            assert len(mets_in) + len(mets_out) == 1
            continue
        (old,), (new,) = mets_in, mets_out
        catalysts = [s for s in orgs_in if orgs_out[s] >= 1]
        assert any(np.array_equal(U.ca_step(bits(struct[old]), int(a[1:])), bits(struct[new]))
                   for a in catalysts)
        # organisms only metabolise strings they extract energy from (entropy increases)
        assert energies[old] - energies[new] > U.TOL
        assert sum(orgs_in.values()) in (1, 2) and sum(orgs_in.values()) == sum(orgs_out.values())
        births += r.count * (sum(orgs_in.values()) == 2)
    a = net.extras["analysis"]
    assert births == sum(a["births_per_update"]) > 0
    assert sum(a["final_population"].values()) == 128
    assert sum(n for s, n in net.initial_state.items() if s.startswith("R")) == 128
    assert sum(n for s, n in net.initial_state.items() if s.startswith("m")) == 640
    sample = [s for s in struct if s.startswith("m")][:200]
    assert np.allclose([energies[s] for s in sample], U.energy(np.array([bits(struct[s]) for s in sample])))
    assert a["shannon_index"][0] == pytest.approx(np.log(128))


def test_diversity_efficiency_and_depth_fall_with_flow_productivity_rises():
    # figs. 4 and 6; strings pass through about 1 metabolic step at gamma = 0.3 and more at 0.003
    def stats(flow):
        h = run(list(range(256)), [], flow=flow, mu=0.01, updates=300, seed=4)["history"]
        mean = lambda xs: float(np.mean([x for x in xs[20:] if x is not None]))
        return (mean(h["shannon"]), mean(h["births"]), mean(h["efficiency"]),
                mean(h["removed_depth"][150:]))
    H_hi, rho_hi, eta_hi, depth_hi = stats(0.3)
    H_lo, rho_lo, eta_lo, depth_lo = stats(0.003)
    assert H_lo > H_hi
    assert rho_hi > 2 * rho_lo
    assert eta_lo > 5 * eta_hi
    assert depth_hi < 1.0 < 3.0 < depth_lo


def test_no_flow_freezes_the_dynamics():
    # "Without this flow, no organismal divisions can occur and the dynamics are frozen"
    # (strings only gain entropy; rare mutants that can still digest them leave a trickle of births)
    h = run(list(range(256)), [], n=64, flow=0.0, mu=0.01, updates=300, seed=1)["history"]
    early, late = np.mean(h["births"][:10]), np.mean(h["births"][-50:])
    assert early > 10 and late < early / 20
    # no inflow and only energy-extracting steps: the pool's mean energy never rises
    assert np.all(np.diff(h["pool_energy"]) <= 1e-12) and h["pool_energy"][-1] < h["pool_energy"][0]


def test_rules_129_and_145_coexist_at_1_to_4():
    # fig. 7 / fig. 8: c(129, 145) = 0.2 from a 9:1 start; one interior fixed point, so c_ij + c_ji ~ 1
    c12, c21 = c(129, 145, seeds=(0, 1)), c(145, 129, seeds=(0, 1))
    assert 0.08 < c12 < 0.35 and 0.65 < c21 < 0.92
    assert abs(c12 + c21 - 1) < 0.15


def test_published_dominance_relations():
    # "species 105 outcompetes 109" in the pairwise experiments; "c(150, 109) = 1"
    assert c(105, 109) == 1.0 and c(109, 105) < 0.5
    assert c(150, 109) == 1.0


def test_rock_paper_scissors_triple_126_134_141():
    # 141 > 126, 126 > 134, 134 > 141: i > j when c_ij > 0.5 and c_ji < 0.5
    def dominates(i, j):
        return c(i, j) > 0.5 and c(j, i) < 0.5
    assert dominates(141, 126) and dominates(126, 134) and dominates(134, 141)


def test_parameters_are_checked():
    with pytest.raises(ValueError, match="0..255"):
        generate_network(ID, rules=[300])
    with pytest.raises(ValueError, match="one weight per rule"):
        generate_network(ID, rules=[1, 2], abundances=[1])
    with pytest.raises(ValueError, match="beta"):
        generate_network(ID, beta=0.0)
    net = generate_network(ID, seed=1, rules=[129, 145], abundances=[9, 1], mutation_rate=0.0,
                           updates=5, transform="always")
    assert set(s.id for s in net.species if s.id.startswith("R")) == {"R129", "R145"}
