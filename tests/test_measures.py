"""chemart.measures: reference values on small networks worked out by hand."""

import numpy as np
import pytest

import chemart
from chemart import measures
from chemart.helpers.explicit import network
from chemart.trajectory import Frame, Trajectory


def m(net, *names, **kw):
    out = chemart.measure(net, list(names), **kw)
    return out[names[0]] if len(names) == 1 else out


# --- B. stoichiometry and chemical reaction network theory ----------------------------
def test_michaelis_menten_textbook_values():
    # E + S <-> ES -> E + P: total enzyme and total substrate are conserved; deficiency zero
    net = network([("E + S -> ES", 1.0), ("ES -> E + S", 1.0), ("ES -> E + P", 1.0)])
    assert m(net, "conservation_laws") == 2
    assert m(net, "stoichiometric_rank") == 2
    assert m(net, "deficiency") == 0
    assert m(net, "weakly_reversible") is False       # E + P never goes back
    assert m(net, "conservative") is True
    assert m(net, "flux_dimension") == 1


@pytest.mark.parametrize("reactions, delta, wr", [
    (["A -> 2 B", "2 B -> A"], 0, True),                         # 2 complexes, 1 class, rank 1
    (["A -> B", "B -> A", "C -> D", "D -> C"], 0, True),         # 4 complexes, 2 classes, rank 2
    (["A + B -> C", "C -> A + B", "C -> A"], 0, False),          # 3 complexes, 1 class, rank 2
    (["A -> 2 A", "A -> "], 1, False),                           # 3 complexes, 1 class, rank 1
])
def test_deficiency_by_hand(reactions, delta, wr):
    net = network([(r, 1.0) for r in reactions])
    assert m(net, "deficiency") == delta
    assert m(net, "weakly_reversible") is wr


def test_growth_is_not_conservative():
    assert m(network([("A -> 2 A", 1.0)]), "conservative") is False


# --- C. graph topology ------------------------------------------------------------------------
def test_bow_tie_and_flow_hierarchy_on_a_toy_graph():
    # A feeds the loop B <-> C, which drains into D
    net = network([("A -> B", 1.0), ("B -> C", 1.0), ("C -> B", 1.0), ("C -> D", 1.0)])
    assert m(net, "bow_tie") == {"core": 0.5, "in": 0.25, "out": 0.25, "other": 0.0}
    chain = network([("A -> B", 1.0), ("B -> C", 1.0)])
    assert m(chain, "flow_hierarchy") == 1.0                     # nothing feeds back
    assert m(net, "flow_hierarchy") == pytest.approx(4 / 8)      # the 4 edges of the loop are cyclic


def test_nodf_of_a_perfectly_nested_incidence():
    # species A in r1, r2, r3; B in r1, r2; C in r1: every row and column pair nested
    net = network([("A + B -> C", 1.0), ("A -> B", 1.0), ("A -> ", 1.0)])
    assert m(net, "nodf") == pytest.approx(100.0)


def test_cycle_rank_and_degrees():
    net = network([("A -> B", 1.0), ("B -> A", 1.0)])
    assert m(net, "cycle_rank") == 1                             # A-r1-B-r2-A
    assert m(net, "degree_cv") == 0.0 and m(net, "degree_gini") == 0.0
    assert m(net, "reciprocity") == 1.0


def test_catalytic_spectral_radius_finds_the_jain_krishna_cycle():
    # the Jain-Krishna page's network: a 3-cycle plus one parasite, lambda1 = 1
    net = chemart.generate_network("jain-krishna", seed=1, graph_updates=768)
    assert m(net, "catalytic_spectral_radius") == pytest.approx(1.0)
    assert net.extras["analysis"]["perron_frobenius_eigenvalue"] == pytest.approx(1.0)
    # a catalysed chain has no cycle
    chain = network([("A + X -> A + B", 1.0), ("B + Y -> B + C", 1.0)])
    assert m(chain, "catalytic_spectral_radius") == 0.0


def test_spectral_gap_is_zero_for_a_disconnected_network():
    net = network([("A -> B", 1.0), ("C -> D", 1.0)])
    assert m(net, "spectral_gap") == 0.0
    assert m(network([("A -> B", 1.0), ("B -> A", 1.0)]), "spectral_gap") > 0


# --- D. organisation --------------------------------------------------------------------------
def test_max_raf_by_hand():
    # food {a, b}; ab catalyses its own ligation (a RAF); a second ligation needs
    # catalyst c, which nothing makes, so it is left out
    net = network([("a + b + ab -> 2 ab", None), ("ab + a + c -> aba + c", None)],
                  extras={"food": ["a", "b"]})
    assert m(net, "max_raf_fraction") == 0.5
    # the scope ignores catalysts, as the RAF closure does: ab, then aba; never c
    assert m(net, "scope_fraction") == pytest.approx(4 / 5)
    assert m(net, "expansion_depth") == 2


def test_food_falls_back_to_inflow_buffered_and_initial_state():
    net = network([("A -> B", 1.0)], initial_state={"A": 1.0})
    assert measures.food_set(net) == ["A"]
    buffered = network([("A -> B", 1.0)], extras={"buffered": ["B"]}, initial_state={"A": 1.0})
    assert measures.food_set(buffered) == ["B"]
    assert m(net, "scope_fraction", food=["B"]) == 0.5


# --- F. kinetics ------------------------------------------------------------------------------
def test_wegscheider_on_a_triangle():
    balanced = network([("A -> B", 2.0), ("B -> A", 1.0), ("B -> C", 3.0), ("C -> B", 1.0),
                        ("C -> A", 1.0), ("A -> C", 6.0)])        # 2 * 3 / 6 = 1 around the cycle
    assert m(balanced, "wegscheider_residual") == pytest.approx(0.0, abs=1e-12)
    unbalanced = network([("A -> B", 2.0), ("B -> A", 1.0), ("B -> C", 1.0), ("C -> B", 1.0),
                          ("C -> A", 1.0), ("A -> C", 1.0)])
    assert m(unbalanced, "wegscheider_residual") > 0.1
    assert m(unbalanced, "rate_spread") == pytest.approx(np.log10(2.0))


# --- G. dynamics ------------------------------------------------------------------------------
def test_population_measures():
    state = {"A": 2.0, "B": 1.0, "C": 1.0}
    assert m(state, "richness") == 3
    assert m(state, "dominance") == 0.5
    assert m(state, "shannon") == pytest.approx(-(0.5 * np.log(0.5) + 2 * 0.25 * np.log(0.25)))


def test_trajectory_measures_and_over():
    frames = [Frame(0.0, {"A": 5, "B": 5, "C": 5, "D": 5}), Frame(1.0, {"A": 10, "B": 10}, [[["C"], ["A"], 5]]),
              Frame(2.0, {"A": 20}, [[["B"], ["A"], 10]])]
    traj = Trajectory(network=network([("C -> A", None), ("B -> A", None)]), frames=frames, method="evolve")
    assert m(traj, "turnover") == pytest.approx((0.5 + 0.5) / 2)
    assert m(traj, "collapse_time") is None                     # 4 -> 1 species is not below 10% of 4
    assert m(traj, "final_richness_ratio") == 0.25
    series = measures.over(traj, ["richness", "n_reactions"])
    assert series["richness"] == [4, 2, 1]
    assert series["n_reactions"] == [None, 1, 1]


def test_over_a_gas():
    traj = chemart.evolve("alchemy", seed=1)
    series = measures.over(traj, ["richness", "shannon", "n_species"], window=5)
    assert len(series["richness"]) == len(traj.frames)
    assert series["richness"][0] == 100 and series["richness"][-1] == 32
    with pytest.raises(ValueError, match="whole run"):
        measures.over(traj, ["turnover"])


# --- the registry and its API -----------------------------------------------------------------
def test_every_measure_is_documented():
    for info in measures.describe():
        assert info["meaning"] and info["section"] in measures.SECTIONS, info["name"]


def test_inapplicable_measures_are_left_out_with_a_reason():
    net = network([("A -> B", None)])
    out = chemart.measure(net)
    assert "rate_spread" not in out and "max_raf_fraction" not in out
    why = measures.applicable(net)
    assert why["rate_spread"] == "not every reaction has a rate"
    assert why["max_raf_fraction"] == "no reaction has a catalyst"
    assert why["n_species"] is None
    with pytest.raises(ValueError, match="did you mean 'deficiency'"):
        chemart.measure(net, ["deficency"])


def test_null_model_keeps_degrees_and_arities():
    net = chemart.generate_network("random-catalytic-networks", seed=1)
    rnd = measures.randomize(net, np.random.default_rng(0))
    def slots(n, side):
        return sorted(sum(getattr(r, side).get(s.id, 0) for r in n.reactions) for s in n.species)
    for side in ("reactants", "products"):
        assert slots(rnd, side) == slots(net, side)
    assert sorted(sum(r.reactants.values()) for r in rnd.reactions) == sorted(sum(r.reactants.values()) for r in net.reactions)
    z = measures.zscores(net, ["nodf", "n_species"], samples=5)
    assert z["n_species"] is None and isinstance(z["nodf"], float)


def test_sweep_and_scaling():
    rows = measures.sweep("random-catalytic-networks", {"n": [5, 10, 20]}, seeds=[1, 2],
                          names=["n_species", "n_reactions"])
    assert len(rows) == 6 and rows[0]["n"] == 5 and rows[0]["n_species"] == 5
    assert rows == measures.sweep("random-catalytic-networks", {"n": [5, 10, 20]}, seeds=[1, 2],
                                  names=["n_species", "n_reactions"])
    fit = measures.scaling([{"x": x, "y": 3 * x ** 2} for x in (1, 2, 4, 8)], "y", "x")
    assert fit["exponent"] == pytest.approx(2.0) and fit["r2"] == pytest.approx(1.0)


# --- moderate and exponential ---------------------------------------------------------------
def mx(net, *names, **kw):
    return m(net, *names, cost="exponential", **kw)


def test_p_invariants_and_flux_modes():
    mm = network([("E + S -> ES", 1.0), ("ES -> E + S", 1.0), ("ES -> E + P", 1.0)])
    assert mx(mm, "p_invariants") == 2                          # E + ES, and S + ES + P
    chain = network([("A -> B", None), ("B -> C", None)], extras={"food": ["A"]})
    assert mx(chain, "elementary_flux_modes") == {"count": 1, "mean_length": 2.0}
    branch = network([("A -> B", None), ("A -> C", None), ("B -> D", None), ("C -> D", None)],
                     extras={"food": ["A"]})
    assert mx(branch, "elementary_flux_modes") == {"count": 2, "mean_length": 2.0}
    assert mx(network([("A -> B", None), ("C -> D", None)], extras={"food": ["A"]}), "blocked_fraction") == 0.5


def test_modularity_and_motifs():
    two = network([("A -> B", None), ("B -> A", None), ("C -> D", None), ("D -> C", None)])
    assert mx(two, "modularity") > 0.3
    net = chemart.generate_network("random-catalytic-networks", seed=1)
    profile = mx(net, "motif_profile")
    assert len(profile) == 13 and sum(v * v for v in profile.values()) == pytest.approx(1.0)


def test_irreducible_rafs_and_autocatalytic_cores():
    two = network([("a + b + ab -> 2 ab", None), ("c + d + cd -> 2 cd", None)],
                  extras={"food": ["a", "b", "c", "d"]})
    assert mx(two, "irreducible_rafs") == 2
    assert mx(network([("X + A -> 2 X", None)]), "autocatalytic_cores") == 1
    assert mx(network([("A -> B", None), ("B -> C", None)]), "autocatalytic_cores") == 0
    # C makes A and B back; A + B -> 2 C closes either cycle
    cycle = network([("A + B -> 2 C", None), ("C -> A", None), ("C -> B", None)])
    assert mx(cycle, "autocatalytic_cores") == 2


def test_organisations_by_hand():
    assert mx(network([("A -> B", None), ("B -> A", None)]), "organisations") == {"count": 2, "largest": 1.0}
    # A makes B, which only decays: {A, B} is closed but not self-maintaining without a source of A
    assert mx(network([("A -> B", None), ("B -> ", None)]), "organisations") == {"count": 1, "largest": 0.0}


def test_kinetics_of_the_brusselator():
    # a = 1, b = 3: the fixed point (1, 3) has trace b - 1 - a^2 = 1 and determinant a^2 = 1
    net = chemart.generate_network("brusselator", seed=1)
    out = mx(net, "stability", "steady_states", "oscillation", "sloppiness")
    assert out["stability"]["max_real"] == pytest.approx(0.5, abs=1e-6)
    assert out["steady_states"] == 1
    assert out["oscillation"]["oscillates"] and 6 < out["oscillation"]["period"] < 9
    assert out["sloppiness"] > 0
    stable = chemart.generate_network("brusselator", seed=1, b=1.5)
    assert mx(stable, "stability")["max_real"] < 0 and not mx(stable, "oscillation")["oscillates"]


def test_entropy_production_is_zero_at_detailed_balance():
    balanced = network([("A -> B", 1.0), ("B -> A", 2.0)], initial_state={"A": 1.0})
    assert mx(balanced, "entropy_production") == pytest.approx(0.0, abs=1e-9)
    driven = network([("A -> B", 2.0), ("B -> A", 1.0), ("B -> C", 2.0), ("C -> B", 1.0),
                      ("C -> A", 2.0), ("A -> C", 1.0)], initial_state={"A": 1.0})
    assert mx(driven, "entropy_production") > 0.1
    assert mx(driven, "flux_concentration") >= 0


def test_attractor_types():
    from chemart import simulate

    cycle = simulate.ode(chemart.generate_network("brusselator"), 60, points=800)
    assert mx(cycle, "attractor_type") == "cycle"
    settled = simulate.ode(chemart.generate_network("brusselator", b=1.5), 60, points=800)
    assert mx(settled, "attractor_type") == "fixed point"


def test_knockouts_and_degeneracy_by_hand():
    # food A; B is made directly (r1) or through C (r2, r3)
    net = network([("A -> B", None), ("A -> C", None), ("C -> B", None)], extras={"food": ["A"]})
    assert mx(net, "knockout_tolerance") == pytest.approx(2 / 3)     # only A -> C is essential
    assert mx(net, "synthetic_lethal_pairs") == 1.0                  # r1 and r3 back each other up
    assert mx(net, "degeneracy") == pytest.approx((2 + 1) / 2)       # two routes to B, one to C


def test_structure_function_mutual_information():
    net = chemart.generate_network("alchemy", seed=1)
    assert mx(net, "structure_function_mi") >= 0
    small = chemart.measure(network([("A -> B", None)]), ["structure_function_mi"])
    assert small == {}                                   # species carry no structure


def test_limits_hold_back_expensive_measures():
    net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
    why = measures.applicable(net, ["organisations"])
    assert "exceed its limit" in why["organisations"]
    assert chemart.measure(net, ["organisations"]) == {}


# --- every cheap measure on every chemistry ------------------------------------------------
from chemart import catalog  # noqa: E402

IMPLEMENTED = sorted(c.id for c in catalog.load() if c.implemented)


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_cheap_measures_run_on_every_chemistry(cid):
    """Each cheap measure returns a value or is left out with a reason, in seconds."""
    import time

    net = chemart.generate_network(cid, seed=1)
    start = time.perf_counter()
    values = chemart.measure(net)
    elapsed = time.perf_counter() - start
    reasons = measures.applicable(net, cost="cheap")
    for name, why in reasons.items():
        if measures.REGISTRY[name].input == "network":
            assert (name in values) == (why is None), (name, why)
    assert elapsed < 10, f"cheap measures took {elapsed:.1f} s"


def test_measures_page_is_current():
    import subprocess
    import sys

    from chemart.catalog import ROOT

    done = subprocess.run([sys.executable, str(ROOT / "tools" / "gen_measures_page.py"), "--check"],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stdout
