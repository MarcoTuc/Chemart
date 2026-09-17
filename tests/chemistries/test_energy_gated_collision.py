"""The Arrhenius-gated collision algorithm (book 18.3.3), validated against the
equations it is built from. There are no published numbers for this algorithm -
the book gives it as three lines of prose - so the validation is that the
simulation reproduces the closed-form kinetics of chapter 2:

- eq. 2.32, k = A exp(-Ea/RT), with R = 8.31451 J/K/mol (sec. 2.2.6): the
  measured fraction of effective collisions is the Arrhenius factor, over a grid
  of barriers and temperatures, and an Arrhenius plot of ln k against 1/T
  recovers the barrier and the prefactor that were put in;
- the law of mass action (sec. 2.2.3, eq. 2.10-2.12), which the algorithm
  degenerates to when the barrier is zero: a second-order decay then follows
  1/n = 1/n0 + 2 P s / (M (M - 1)) over s collisions of a vessel of M molecules;
- eq. 2.33, K = exp(-dG/RT), reached by a reversible pair whose barriers differ
  by dG as fig. 2.4 prescribes;
- eq. 2.28, E = Ek + Ep, conserved exactly by the `conserved` energy model.

The `conserved` model also measures the discrepancy the book's literal wording
carries: gating on the *sum* of two Maxwell-Boltzmann kinetic energies accepts
with the Gamma-2 tail (1 + Ea/RT) exp(-Ea/RT), not with the Arrhenius factor.
"""

from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.energy_gated_collision import (
    R_GAS, acceptance_probability, arrhenius_factor, default_system, run,
    validated_system)
from chemart.kinetics import RATE_LAWS

FORWARD = "X1+X2->Y1+Y2"
REVERSE = "Y1+Y2->X1+X2"


def gated_run(Ea, T, seed=11, steps=8000, delta_G=-3.0, population=400, model="bath"):
    """Run the algorithm on the default reversible system; return (system, result)."""
    system = validated_system(default_system(Ea, delta_G))
    result = run(system, np.random.default_rng(seed), temperature=T,
                 population=population, steps=steps, energy_model=model)
    return system, result


def acceptance(result, rid):
    tally = result["stats"][rid]
    return tally["effective"] / tally["attempts"], tally["attempts"]


def sigmas(observed, predicted, n):
    """How many binomial standard errors the measurement is from the prediction."""
    return abs(observed - predicted) / np.sqrt(predicted * (1 - predicted) / n)


def one_way(Ea):
    """X1 + X2 -> Y1 + Y2 only, so the reactants decay instead of equilibrating."""
    return {"energies": {"X1": 0.0, "X2": 0.0, "Y1": -3.0, "Y2": 0.0},
            "reactions": [{"id": "fwd", "reactants": ["X1", "X2"],
                           "products": ["Y1", "Y2"], "Ea": Ea}],
            "initial": {"X1": 1.0, "X2": 1.0}}


def decayed_acceptance(Ea, seed, steps=4000, population=400):
    """Invert 1/n = 1/n0 + 2 P s / (M (M-1)) for the acceptance probability P.

    The mass-action rate of X1 + X2 -> ... in the collision algorithm is the
    chance of drawing that pair, 2 n1 n2 / (M (M-1)), times the chance the
    collision clears the barrier; with n1 = n2 = n this integrates to the
    second-order decay law above.
    """
    system = validated_system(one_way(Ea))
    result = run(system, np.random.default_rng(seed), temperature=300.0,
                 population=population, steps=steps)
    n0, nf = population // 2, Counter(result["final"])["X1"]
    return (1.0 / nf - 1.0 / n0) * population * (population - 1) / (2.0 * steps)


# --- the gate itself --------------------------------------------------------------
def test_gas_constant_and_gate_match_the_book():
    """R = 8.31451 J/K/mol (sec. 2.2.6, under eq. 2.32) and the two acceptance laws."""
    assert R_GAS * 1000 == pytest.approx(8.31451)
    assert RATE_LAWS["arrhenius"] == ("A", "Ea")          # the law's exact parameter names
    assert R_GAS * 300 == pytest.approx(2.494353)          # RT at room temperature, kJ/mol

    for Ea, T in [(2.0, 300.0), (8.0, 400.0), (0.0, 300.0)]:
        x = Ea / (R_GAS * T)
        assert arrhenius_factor(Ea, T) == pytest.approx(np.exp(-x))
        assert acceptance_probability(Ea, T, "bath") == pytest.approx(np.exp(-x))
        assert acceptance_probability(Ea, T, "conserved") == pytest.approx((1 + x) * np.exp(-x))
    # A barrierless collision always reacts, whichever model books the energy.
    assert acceptance_probability(0.0, 300.0, "bath") == 1.0
    assert acceptance_probability(0.0, 300.0, "conserved") == 1.0


@pytest.mark.parametrize("Ea, T", [(1.0, 300.0), (2.0, 300.0), (5.0, 300.0),
                                   (2.0, 200.0), (2.0, 500.0), (8.0, 400.0)])
def test_measured_frequency_is_the_arrhenius_factor(Ea, T):
    """Eq. 2.32: the fraction of collisions of a reactant pair that are effective
    is exp(-Ea/RT). Both channels of the reversible system are measured, so each
    run checks two different barriers (Ea and Ea - dG = Ea + 3)."""
    system, result = gated_run(Ea, T)
    checked = 0
    for r in system["reactions"]:
        observed, attempts = acceptance(result, r["id"])
        if attempts < 200:
            continue
        predicted = arrhenius_factor(r["Ea"], T)
        assert sigmas(observed, predicted, attempts) < 4, (
            f"{r['id']}: Ea={r['Ea']} T={T} measured {observed:.4f}, "
            f"Arrhenius {predicted:.4f} over {attempts} collisions")
        checked += 1
    assert checked == 2


def test_arrhenius_plot_recovers_the_barrier_and_the_prefactor():
    """The published way of validating eq. 2.32: ln k against 1/T is a straight
    line of slope -Ea/R and intercept ln A. Here k = A x (effective collisions /
    collisions of the reactants), measured at five temperatures."""
    Ea, A = 4.0, 2.5
    inverse_T, ln_k = [], []
    for T in (250.0, 300.0, 350.0, 420.0, 500.0):
        _, result = gated_run(Ea, T, seed=23, steps=12000)
        observed, attempts = acceptance(result, FORWARD)
        assert attempts > 400
        inverse_T.append(1.0 / T)
        ln_k.append(np.log(A * observed))
    slope, intercept = np.polyfit(inverse_T, ln_k, 1)
    assert -slope * R_GAS == pytest.approx(Ea, rel=0.15)      # slope = -Ea/R
    assert np.exp(intercept) == pytest.approx(A, rel=0.20)    # intercept = ln A
    residuals = ln_k - (slope * np.array(inverse_T) + intercept)
    assert np.max(np.abs(residuals)) < 0.12, "the Arrhenius plot must be straight"


def test_raising_the_barrier_or_lowering_the_temperature_suppresses_the_reaction():
    """Ea and T enter only through Ea/RT, and the gate is exponential in it."""
    by_barrier = [acceptance(gated_run(Ea, 300.0, seed=31)[1], FORWARD)[0]
                  for Ea in (0.0, 1.0, 2.0, 4.0, 8.0)]
    assert by_barrier == sorted(by_barrier, reverse=True), by_barrier
    assert by_barrier[0] == 1.0                     # no barrier: every collision reacts
    assert by_barrier[-1] < 0.1 * by_barrier[1]     # 8 kJ/mol is ~3.2 RT

    by_temperature = [acceptance(gated_run(5.0, T, seed=31)[1], FORWARD)[0]
                      for T in (150.0, 250.0, 400.0, 700.0)]
    assert by_temperature == sorted(by_temperature), by_temperature
    assert by_temperature[0] < 0.1 * by_temperature[-1]


# --- the mass-action limit --------------------------------------------------------
def test_zero_barrier_reduces_to_the_naive_collision_algorithm():
    """With Ea = 0 the gate always opens and the algorithm is exactly the naive
    bimolecular collision algorithm of sec. 2.3.3: every collision of a reactive
    pair is effective, and the only elastic collisions are the pairs with no rule."""
    system, result = gated_run(0.0, 300.0, seed=41, delta_G=0.0, steps=5000)
    attempted = 0
    for r in system["reactions"]:
        observed, attempts = acceptance(result, r["id"])
        assert observed == 1.0, "a barrierless collision can never be rejected"
        attempted += attempts
    effective = sum(t["effective"] for t in result["stats"].values())
    assert effective == attempted < 5000            # the rest drew a non-reactive pair


@pytest.mark.parametrize("Ea", [0.0, 2.0, 5.0])
def test_second_order_decay_measures_the_arrhenius_rate_constant(Ea):
    """Mass action (eq. 2.12) times the Arrhenius factor (eq. 2.32): the decay of
    a one-way second-order reaction over s collisions gives back exp(-Ea/RT).
    The comparison is against the mean-field solution, which is biased high by a
    few percent at this population size, hence the loose tolerance."""
    estimates = [decayed_acceptance(Ea, 100 + s) for s in range(5)]
    assert np.mean(estimates) == pytest.approx(arrhenius_factor(Ea, 300.0), rel=0.25)


@pytest.mark.slow
def test_equilibrium_constant_is_the_exponential_of_minus_delta_G():
    """Eq. 2.33: forward and reverse barriers differing by dG (fig. 2.4) drive the
    vessel to K = kf/kr = exp(-dG/RT). Measured as [Y1][Y2]/[X1][X2] at the end of
    a long run, averaged over independent vessels."""
    for delta_G, Ea in [(-3.0, 2.0), (2.0, 3.0)]:
        measured = []
        for s in range(6):
            _, result = gated_run(Ea, 300.0, seed=200 + s, delta_G=delta_G, steps=20000)
            n = Counter(result["final"])
            measured.append(n["Y1"] * n["Y2"] / (n["X1"] * n["X2"]))
        predicted = np.exp(-delta_G / (R_GAS * 300.0))
        assert np.mean(measured) == pytest.approx(predicted, rel=0.25), (delta_G, measured)


# --- energy bookkeeping -----------------------------------------------------------
def test_conserved_model_conserves_total_energy_exactly():
    """Eq. 2.28, E = Ek + Ep: what the barrier does not consume goes to the
    products, so the closed vessel keeps its total energy to floating-point
    exactness, while the exothermic reactions move it from potential to kinetic."""
    for delta_G in (0.0, -3.0, 2.0):
        _, result = gated_run(3.0, 300.0, seed=5, delta_G=delta_G, steps=8000,
                              model="conserved")
        budget = result["budget"]
        assert abs(budget["drift"]) < 1e-9 * abs(budget["total_initial"])
        assert budget["total_final"] == pytest.approx(budget["total_initial"])
        assert budget["kinetic_initial"] == pytest.approx(
            400 * R_GAS * 300.0, rel=0.15)        # M molecules of mean energy RT
        if delta_G < 0:                            # exergonic: the vessel heats up
            assert budget["kinetic_final"] > budget["kinetic_initial"]
            assert budget["potential_final"] < budget["potential_initial"]

    # The bath model instead books the reaction heat against the thermostat.
    _, bath = gated_run(3.0, 300.0, seed=5, delta_G=-3.0, steps=8000)
    assert bath["budget"]["kinetic_tracked"] is False
    assert bath["budget"]["supplied_by_bath"] == pytest.approx(
        bath["budget"]["potential_final"] - bath["budget"]["potential_initial"])


@pytest.mark.parametrize("Ea", [1.0, 2.0, 4.0])
def test_conserved_model_follows_the_pair_sum_tail_not_the_arrhenius_factor(Ea):
    """The book's literal wording in sec. 2.2.6 gates on the *sum* of the two
    molecules' kinetic energies. Two Exponential(RT) energies sum to a Gamma(2,
    RT), so the acceptance is (1 + Ea/RT) exp(-Ea/RT) - well above the Arrhenius
    factor. Measured on a thermoneutral system, where the kinetic pool keeps its
    Maxwell-Boltzmann shape."""
    system, result = gated_run(Ea, 300.0, seed=5, delta_G=0.0, steps=8000,
                               model="conserved")
    predicted = acceptance_probability(Ea, 300.0, "conserved")
    factor = arrhenius_factor(Ea, 300.0)
    assert predicted > 1.3 * factor, "the two readings must be distinguishable"
    for r in system["reactions"]:
        observed, attempts = acceptance(result, r["id"])
        assert attempts > 200
        assert sigmas(observed, predicted, attempts) < 4, (observed, predicted)
        assert abs(observed - factor) > abs(observed - predicted)


# --- the generated network --------------------------------------------------------
def test_default_network_is_the_observed_gated_run():
    net = generate_network("energy-gated-collision", seed=1)
    analysis = net.extras["analysis"]
    assert net.status == "observed"
    assert [s.id for s in net.species] == ["X1", "X2", "Y1", "Y2"]
    assert net.initial_state == {"X1": 200.0, "X2": 200.0}
    assert len(net.reactions) == 2 and all(r.count > 0 for r in net.reactions)

    for r in net.reactions:
        assert r.rate["law"] == "arrhenius"
        assert r.rate["A"] == 1.0 and r.rate["T"] == 300.0 and r.rate["R"] == R_GAS
    forward = next(r for r in net.reactions if r.reactants == {"X1": 1, "X2": 1})
    assert forward.products == {"Y1": 1, "Y2": 1} and forward.rate["Ea"] == 2.0

    energies = net.extras["energies"]
    assert energies["potential"] == {"X1": 0.0, "X2": 0.0, "Y1": -3.0, "Y2": 0.0}
    assert energies["activation"] == {FORWARD: 2.0, REVERSE: 5.0}   # Ea - dG = 2 + 3
    assert energies["units"] == "kJ/mol"

    # Counts, tallies and the book's accounting of collisions all agree.
    effective = sum(r.count for r in net.reactions)
    assert analysis["effective_collisions"] == effective
    assert analysis["elastic_collisions"] == 4000 - effective
    assert analysis["generations"] == 10.0                          # 4000 / M, sec. 2.6.1
    assert analysis["elastic_fraction"] > 0.5, "the book's caveat: mostly elastic"
    for rid, measured in analysis["reactions"].items():
        assert measured["acceptance"] == pytest.approx(
            measured["effective"] / measured["attempts"])
        assert measured["k_observed"] == pytest.approx(1.0 * measured["acceptance"])
        assert measured["k_arrhenius"] == pytest.approx(
            arrhenius_factor(measured["Ea"], 300.0))
        assert sigmas(measured["acceptance"], measured["arrhenius_factor"],
                      measured["attempts"]) < 4
    assert set(net.provides) <= {"topology", "stoichiometry", "rate-constants",
                                 "rate-law", "energies", "initial-state"}


def test_same_seed_same_run_and_parameters_reach_the_gate():
    a = generate_network("energy-gated-collision", seed=5)
    assert a.to_dict() == generate_network("energy-gated-collision", seed=5).to_dict()
    assert generate_network("energy-gated-collision", seed=6).extras["analysis"] \
        != a.extras["analysis"]

    # A colder vessel with a higher barrier fires far less often.
    cold = generate_network("energy-gated-collision", seed=5, temperature=200.0, Ea=6.0)
    assert cold.extras["analysis"]["effective_collisions"] < \
        0.2 * a.extras["analysis"]["effective_collisions"]
    # The prefactor scales every rate constant and nothing else.
    scaled = generate_network("energy-gated-collision", seed=5, A=10.0)
    assert scaled.extras["analysis"]["effective_collisions"] == \
        a.extras["analysis"]["effective_collisions"]
    assert all(r.rate["A"] == 10.0 for r in scaled.reactions)


def test_user_supplied_system_and_invalid_ones():
    system = one_way(1.0)
    net = generate_network("energy-gated-collision", seed=3, system=system, steps=2000)
    assert len(net.reactions) == 1 and net.reactions[0].rate["Ea"] == 1.0
    assert net.extras["energies"]["potential"]["Y1"] == -3.0

    with pytest.raises(ValueError, match="below its energy change"):
        generate_network("energy-gated-collision", Ea=1.0, delta_G=3.0)
    with pytest.raises(ValueError, match="must be bimolecular"):
        generate_network("energy-gated-collision", system={
            "energies": {"a": 0.0, "b": 0.0},
            "reactions": [{"reactants": ["a"], "products": ["b"], "Ea": 1.0}]})
    with pytest.raises(ValueError, match="no entry in system\\['energies'\\]"):
        generate_network("energy-gated-collision", system={
            "energies": {"a": 0.0, "b": 0.0},
            "reactions": [{"reactants": ["a", "b"], "products": ["c"], "Ea": 1.0}]})
    with pytest.raises(ValueError, match="thermodynamically consistent"):
        generate_network("energy-gated-collision", system={
            "energies": {"a": 0.0, "b": 0.0, "c": -4.0, "d": 0.0},
            "reactions": [
                # consistency (fig. 2.4) would need Ea(r) = Ea(f) - dG = 2 + 4 = 6
                {"id": "f", "reactants": ["a", "b"], "products": ["c", "d"], "Ea": 2.0},
                {"id": "r", "reactants": ["c", "d"], "products": ["a", "b"], "Ea": 5.0}]})
    with pytest.raises(ValueError, match="does not change the multiset"):
        generate_network("energy-gated-collision", system={
            "energies": {"a": 0.0, "b": 0.0},
            "reactions": [{"reactants": ["a", "b"], "products": ["b", "a"], "Ea": 1.0}]})
    with pytest.raises(ValueError, match="temperature"):
        generate_network("energy-gated-collision", temperature=0.0)
