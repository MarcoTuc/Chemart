"""W1 written-down networks: published dynamics, reproduced by integrating the generated network.

Each test names the book statement or figure it checks.
"""

import math

import numpy as np
import pytest
from odes import integrate

from chemart import generate_network


def test_dimerization_reaches_equilibrium_constant():
    # book eq. 2.27: [A.B] / ([A][B]) = k_f / k_r
    x, _ = integrate(generate_network("dimerization", k_f=2.0, k_r=0.2), 50)
    assert x["C"][-1] / (x["A"][-1] * x["B"][-1]) == pytest.approx(10.0, rel=1e-4)


def test_brusselator_oscillates_only_above_hopf():
    # book 17.4.2: sustained oscillations for b > 1 + a^2 (all k = 1, a = 1)
    x, _ = integrate(generate_network("brusselator", b=3.0), 60, t_eval=np.linspace(40, 60, 400))
    assert np.ptp(x["X"]) > 1.0
    x, _ = integrate(generate_network("brusselator", b=1.7), 300)
    assert x["X"][-1] == pytest.approx(1.0, rel=1e-3)
    assert x["Y"][-1] == pytest.approx(1.7, rel=1e-3)


def test_logistic_plateaus_at_carrying_capacity():
    # book eq. 7.11 / figure 7.2
    for K in (1.0, 5.0):
        x, _ = integrate(generate_network("logistic-chemistry", K=K), 40)
        assert x["X"][-1] == pytest.approx(K, rel=1e-4)


def test_replication_death_is_exponential():
    # book eq. 7.4
    x, _ = integrate(generate_network("replication-death", b=1.0, d=0.5, x0=2.0), 3)
    assert x["X"][-1] == pytest.approx(2.0 * math.exp(1.5), rel=1e-5)


def test_selection_equation_survival_of_the_fittest():
    # book 7.2.5: competitive exclusion under constant total concentration
    x, _ = integrate(generate_network("selection-equation"), 40)
    assert x["X3"][-1] > 0.99
    assert x["X1"][-1] + x["X2"][-1] + x["X3"][-1] == pytest.approx(1.0, rel=1e-6)


def test_selection_equation_subexponential_growth_keeps_everybody():
    # book 7.2.5: c < 1 gives "survival of everybody"
    x, _ = integrate(generate_network("selection-equation", c=0.5), 200)
    assert min(x["X1"][-1], x["X2"][-1], x["X3"][-1]) > 0.05


def test_replicator_rock_paper_scissors_conserves_product():
    # zero-sum RPS replicator dynamics conserves x1 x2 x3
    start = {"X1": 0.5, "X2": 0.3, "X3": 0.2}
    x, _ = integrate(generate_network("replicator-equation"), 20, x0=start)
    assert x["X1"][-1] * x["X2"][-1] * x["X3"][-1] == pytest.approx(0.5 * 0.3 * 0.2, rel=1e-5)


def test_lotka_volterra_matches_book_reactions_and_conserves_invariant():
    net = generate_network("lotka-volterra")
    assert net.to_text().splitlines() == [
        "X1 -> 2 X1  [mass-action k=1.0]",       # eq. 7.25 with G folded in
        "X2 -> ∅  [mass-action k=1.0]",          # eq. 7.27
        "X2 + X1 -> 2 X2  [mass-action k=1.0]",  # eq. 7.26
    ]
    x, _ = integrate(net, 20)
    def V(i):
        return x["X1"][i] - math.log(x["X1"][i]) + x["X2"][i] - math.log(x["X2"][i])
    assert V(-1) == pytest.approx(V(0), rel=1e-5)


def test_repressilator_needs_cooperativity_to_oscillate():
    # book 19.3.2: oscillations for n = 2 (figure 19.18), none for n = 1
    x, _ = integrate(generate_network("repressilator"), 400, t_eval=np.linspace(300, 400, 500))
    assert np.ptp(x["P1"]) > 1.0
    x, _ = integrate(generate_network("repressilator", n=1), 2000, t_eval=np.linspace(1800, 2000, 500))
    assert np.ptp(x["P1"]) < 0.01 * x["P1"].mean()


def test_michaelis_menten_abridged_matches_elementary_when_enzyme_is_scarce():
    # book eq. 18.4: quasi-steady state holds for E0 << S0
    kwargs = dict(E0=0.01, S0=10.0, ka=10.0, ka_rev=1.0, kb=1.0)
    elementary, _ = integrate(generate_network("michaelis-menten", **kwargs), 200)
    abridged, _ = integrate(generate_network("michaelis-menten", form="abridged", **kwargs), 200)
    assert abridged["P"][-1] == pytest.approx(elementary["P"][-1], rel=0.02)


def test_hill_elementary_equilibrium_matches_hill_function():
    # book eq. 18.8: C / G0 = P^n / (K^n + P^n) at equilibrium, with P held at its level
    net = generate_network("hill-kinetics", n=4, K=1.2, P0=1000.0)
    x, _ = integrate(net, 50)
    P = x["P"][-1]
    assert x["C"][-1] == pytest.approx(P**4 / (1.2**4 + P**4), rel=1e-4)


def flip_time(**params):
    t = np.linspace(0, 120, 1201)
    x, _ = integrate(generate_network("okamoto-switch", **params), 120, t_eval=t, method="Radau")
    crossing = t[np.argmax(x["I1"] <= x["I2"])]
    return crossing, t[np.argmax(x["A"] < 0.5)], x


def test_okamoto_switch_flips_after_inputs_cross():
    # book figure 17.5: A holds while I1 > I2 and flips to B after the inputs cross;
    # slowing the I1 -> I2 conversion delays the switch (bottom panel).
    crossing, flip, x = flip_time()
    assert crossing == pytest.approx(17.5, abs=0.5)
    assert x["A"][100] > 0.99 and x["B"][100] < 0.01   # t = 10 s
    assert crossing < flip and x["A"][-1] < 0.01 and x["B"][-1] > 0.99
    slow_crossing, slow_flip, _ = flip_time(k_conv=0.003)
    assert slow_crossing > crossing and slow_flip > flip


@pytest.mark.parametrize("function, x1, x2", [
    ("sqrt", 9.0, 0.0), ("square", 3.0, 0.0), ("add", 2.0, 5.0), ("multiply", 2.0, 5.0), ("divide", 6.0, 4.0),
])
def test_analog_functions_compute_at_steady_state(function, x1, x2):
    net = generate_network("analog-function-crn", function=function, x1=x1, x2=x2)
    x, _ = integrate(net, 60)
    assert x["Y"][-1] == pytest.approx(net.extras["readout"]["steady_state"], rel=1e-4)


def test_disperser_converges_to_network_average():
    # book 17.3.1: 1000 jobs injected at node 4 end at ~250 per node
    x, _ = integrate(generate_network("disperser"), 60)
    for node in "1234":
        assert x[f"X{node}"][-1] == pytest.approx(250.0, rel=1e-4)
