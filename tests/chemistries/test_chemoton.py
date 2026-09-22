"""The generated Chemoton reproduces model I of Fernando & Di Paolo (2004), eqs. 1-14."""

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import generate_network
from chemart.chemistries.chemoton import RATES


def paper_derivatives(c, N, k):
    """Eqs. 1-14 as printed (X and Y are held constant, so they are omitted)."""
    A1, A2, A3, A4, A5 = (c[f"A{i}"] for i in range(1, 6))
    X, Y, V, R, T0, Ts, T, S = (c[s] for s in ("X", "Y", "V", "R", "T0", "T*", "T", "S"))
    pV = [c[f"pV{r}"] for r in range(N)]
    propagation = sum(k["k7"] * pV[r] * V for r in range(1, N))
    d = {
        "A1": 2 * (k["k5"] * A5 - k["k5r"] * A1 * A1) - k["k1"] * A1 * X + k["k1r"] * A2,
        "A2": k["k1"] * A1 * X - k["k1r"] * A2 - k["k2"] * A2 + k["k2r"] * A3 * Y,
        "A3": k["k2"] * A2 - k["k2r"] * A3 * Y - k["k3"] * A3 + k["k3r"] * A4 * V,
        "A4": k["k3"] * A3 - k["k3r"] * A4 * V - k["k4"] * A4 + k["k4r"] * A5 * T0,
        "A5": k["k4"] * A4 - k["k4r"] * A5 * T0 - k["k5"] * A5 + k["k5r"] * A1 * A1,
        "pV0": 2 * k["k7"] * pV[N - 1] * V + k["k6r"] * pV[1] * R - k["k6"] * pV[0] * V,
        "pV1": k["k6"] * pV[0] * V - k["k6r"] * pV[1] * R - k["k7"] * pV[1] * V,
        "V": k["k3"] * A3 - k["k3r"] * A4 * V + k["k6r"] * pV[1] * R - k["k6"] * pV[0] * V - propagation,
        "R": k["k6"] * pV[0] * V - k["k6r"] * pV[1] * R + k["k9r"] * T - k["k9"] * Ts * R + propagation,
        "T0": k["k4"] * A4 - k["k4r"] * A5 * T0 - k["k8"] * T0,
        "T*": k["k8"] * T0 - k["k9"] * Ts * R + k["k9r"] * T,
        "T": k["k9"] * Ts * R - k["k9r"] * T - k["k10"] * T * S,
        "S": k["k10"] * T * S,
    }
    for r in range(2, N):
        d[f"pV{r}"] = k["k7"] * pV[r - 1] * V - k["k7"] * pV[r] * V
    return d


@pytest.mark.parametrize("N", [2, 3, 25])
def test_network_rate_equations_match_the_paper(N):
    net = generate_network("chemoton", N=N)
    ids, f = rhs(net)
    rng = np.random.default_rng(0)
    state = dict(zip(ids, rng.uniform(0.5, 5.0, len(ids))))
    state["V"] = 40.0  # above the threshold, so initiation is active in the paper's model
    x = np.array([state[s] for s in ids])
    ours = dict(zip(ids, f(0.0, x)))
    expected = paper_derivatives(state, N, RATES)
    for species, value in expected.items():
        assert ours[species] == pytest.approx(value, rel=1e-9, abs=1e-9), species


def test_initiation_carries_its_threshold():
    net = generate_network("chemoton", V_threshold=30.0)
    (initiation,) = [r for r in net.reactions if r.reactants == {"pV0": 1, "V": 1}]
    assert initiation.rate["threshold_species"] == "V" and initiation.rate["threshold"] == 30.0


def test_rate_overrides_are_validated():
    assert generate_network("chemoton", rates={"k1": 5.0}).reactions[0].rate["k"] == 5.0
    with pytest.raises(ValueError, match="unknown constant"):
        generate_network("chemoton", rates={"k99": 1.0})
