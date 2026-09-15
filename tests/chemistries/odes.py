"""Minimal ODE integration of generated networks, for testing published dynamics.

Supports the rate laws Chemart emits for explicit chemistries, the
`constant-total` dilution flux, and species held constant via
extras['buffered']. This is test scaffolding, not a simulator.
"""

import numpy as np
from scipy.integrate import solve_ivp

from chemart.network import CONSTANT_TOTAL


def rhs(net):
    ids, R, P = net.matrices()
    index = {s: i for i, s in enumerate(ids)}
    S = (P - R).toarray().astype(float)
    reactants = [[(index[s], n) for s, n in r.reactants.items()] for r in net.reactions]
    rates = [r.rate for r in net.reactions]
    buffered = [index[s] for s in net.extras.get("buffered", [])]

    def propensity(rate, reac, x):
        law = rate["law"]
        if law == "mass-action":
            return rate["k"] * np.prod([x[i] ** n for i, n in reac])
        (i, _), = reac
        if law == "power":
            return rate["k"] * max(x[i], 0.0) ** rate["order"]
        if law == "michaelis-menten":
            return rate["vmax"] * x[i] / (rate["km"] + x[i])
        raise NotImplementedError(law)

    def hill(rate, x):
        p = max(x[index[rate["regulator"]]], 0.0)
        h = p ** rate["n"] / (rate["K"] ** rate["n"] + p ** rate["n"])
        return rate["vmax"] * (h if rate["mode"] == "activation" else 1.0 - h)

    def f(t, x):
        v = np.array([
            hill(rate, x) if rate["law"] == "hill" else propensity(rate, reac, x)
            for rate, reac in zip(rates, reactants)
        ])
        dx = S @ v
        if net.outflow == CONSTANT_TOTAL:
            dx = dx - x * dx.sum() / x.sum()
        dx[buffered] = 0.0
        return dx

    return ids, f


def integrate(net, t_end, x0=None, t_eval=None, method="LSODA"):
    ids, f = rhs(net)
    start = x0 if x0 is not None else (net.initial_state or {})
    x = np.array([float(start.get(s, 0.0)) for s in ids])
    sol = solve_ivp(f, (0.0, t_end), x, method=method, t_eval=t_eval, rtol=1e-8, atol=1e-10)
    assert sol.success, sol.message
    return {s: sol.y[i] for i, s in enumerate(ids)}, sol.t
