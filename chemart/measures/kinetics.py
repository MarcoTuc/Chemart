"""F. Kinetics: measures that need rate constants."""

from __future__ import annotations

from typing import Any

import numpy as np

from chemart.measures import register


def _mass_action_k(ctx) -> list[float]:
    return [float(r.rate["k"]) for r in ctx.net.reactions
            if r.rate["law"] == "mass-action" and float(r.rate["k"]) > 0]


@register("rate_spread", "F", needs="K")
def rate_spread(ctx) -> float | None:
    """Orders of magnitude spanned by the mass-action rate constants,
    log10(max k / min k). None without mass-action rates."""
    k = _mass_action_k(ctx)
    return float(np.log10(max(k) / min(k))) if k else None


def reversible_pairs(net) -> list[tuple[int, int]]:
    """Index pairs (forward, backward) of reactions that are each other's reverse."""
    seen: dict[tuple, int] = {}
    pairs = []
    for j, r in enumerate(net.reactions):
        key = (frozenset(r.reactants.items()), frozenset(r.products.items()))
        back = (key[1], key[0])
        if back in seen:
            pairs.append((seen.pop(back), j))
        else:
            seen[key] = j
    return pairs


@register("wegscheider_residual", "F", needs="K S")
def wegscheider_residual(ctx) -> float | None:
    """How far the reversible mass-action pairs are from allowing detailed
    balance (Wegscheider's conditions): the least-squares residual of
    ln(k+/k-) against the reactions' stoichiometry, 0 when some chemical
    potentials make every pair balance. None without reversible
    mass-action pairs."""
    rows, target = [], []
    for f, b in reversible_pairs(ctx.net):
        rf, rb = ctx.net.reactions[f].rate, ctx.net.reactions[b].rate
        if rf["law"] == rb["law"] == "mass-action" and rf["k"] > 0 and rb["k"] > 0:
            rows.append(ctx.S[:, f])
            target.append(np.log(rf["k"] / rb["k"]))
    if not rows:
        return None
    A, y = np.array(rows, dtype=float), np.array(target)
    # detailed balance at amounts x*: k+ x*^R = k- x*^P, i.e. ln(k+/k-) = S_jᵀ ln x*
    mu, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.linalg.norm(A @ mu - y))


# --------------------------------------------------------------------------
# Moderate: steady states and dynamics, from the rate equations
# --------------------------------------------------------------------------
#
# These measures run the network's rate equations (chemart.simulate), so they
# need rates and an initial state. They look at the *active* species: those a
# reaction consumes and that are not buffered. Pure products (waste such as
# the Brusselator's D and E) never act back on the dynamics and grow without
# bound, so a steady state is sought for the active species only.

class Dynamics:
    """The rate equations of a network, restricted to its active species."""

    def __init__(self, ctx):
        from scipy.linalg import null_space

        from chemart.simulate import rhs

        net = ctx.net
        self.net = net
        self.ids, f = rhs(net)
        self.f = f
        buffered = set(net.extras.get("buffered", []))
        consumed = {s for r in net.reactions for s, n in r.reactants.items() if n > r.products.get(s, 0)}
        self.active = [i for i, s in enumerate(self.ids) if s in consumed and s not in buffered]
        self.x0 = np.array([float((net.initial_state or {}).get(s, 0.0)) for s in self.ids])
        Sa = ctx.S[self.active].astype(float)
        self.L = null_space(Sa.T).T if Sa.size else np.zeros((0, len(self.active)))
        self._fixed = None

    def full(self, xa):
        x = self.x0.copy()
        x[self.active] = xa
        return x

    def fa(self, xa):
        return self.f(0.0, self.full(xa))[self.active]

    def jacobian(self, xa, h=1e-7):
        n = len(xa)
        J = np.zeros((n, n))
        for k in range(n):
            step = h * max(1.0, abs(xa[k]))
            up, dn = xa.copy(), xa.copy()
            up[k] += step
            dn[k] = max(dn[k] - step, 0.0)
            J[:, k] = (self.fa(up) - self.fa(dn)) / (up[k] - dn[k])
        return J

    def solve(self, start, target):
        """A fixed point of the active species in the compatibility class of `target`."""
        from scipy.optimize import least_squares

        L = self.L

        def residual(xa):
            return np.concatenate([self.fa(xa), L @ (xa - target)])

        res = least_squares(residual, np.maximum(start, 0.0), bounds=(0.0, np.inf), xtol=1e-12, ftol=1e-12,
                            gtol=1e-12, max_nfev=2000)
        scale = max(1.0, float(np.abs(target).max()))
        return res.x if np.abs(residual(res.x)).max() < 1e-7 * scale else None

    def fixed_point(self, t_end: float = 200.0):
        """The fixed point reached, or approached, from the initial state (cached)."""
        if self._fixed is None:
            from chemart.simulate import integrate

            traj, _ = integrate(self.net, t_end)
            end = np.array([traj[s][-1] for s in self.ids])[self.active]
            self._fixed = self.solve(end, self.x0[self.active])
        return self._fixed


def _dynamics(ctx) -> Dynamics:
    d = getattr(ctx, "_dynamics", None)
    if d is None:
        d = ctx._dynamics = Dynamics(ctx)
    return d


def _needs_state(ctx) -> bool:
    return bool(ctx.net.initial_state)


def _spectrum(d, x) -> np.ndarray:
    """Jacobian eigenvalues at x, without the zero modes of the conservation laws."""
    ev = np.linalg.eigvals(d.jacobian(x))
    keep = np.argsort(np.abs(ev))[len(d.L):]
    return ev[keep]


@register("stability", "F", needs="K", cost="moderate", limit=400)
def stability(ctx) -> dict[str, float] | None:
    """At the fixed point of the active species reached from the initial state:
    the largest real part of the Jacobian's eigenvalues (negative: stable) and
    the stiffness ratio (fastest over slowest relaxation rate). None without an
    initial state or when no fixed point is found."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    d = _dynamics(ctx)
    x = d.fixed_point()
    if x is None:
        return None
    ev = _spectrum(d, x)
    if not len(ev):
        return {"max_real": 0.0, "stiffness": 1.0}
    rates = np.abs(ev.real[np.abs(ev.real) > 1e-12])
    return {"max_real": float(ev.real.max()),
            "stiffness": float(rates.max() / rates.min()) if len(rates) else 1.0}


@register("steady_states", "F", needs="K", cost="moderate", limit=200)
def steady_states(ctx, starts: int = 12) -> int | None:
    """Distinct non-negative fixed points of the active species, found from
    `starts` (12) random starts that keep the initial state's conserved totals."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    d = _dynamics(ctx)
    target = d.x0[d.active]
    rng = ctx.rng()
    Sa = ctx.S[d.active].astype(float)
    scale = max(float(target.mean()), 1.0)
    found: list[np.ndarray] = []
    first = d.fixed_point()
    candidates = [first] if first is not None else []
    for _ in range(starts):
        for _ in range(100):                      # a random start in the same compatibility class
            start = target + Sa @ rng.normal(0.0, scale, Sa.shape[1])
            if (start >= 0).all():
                break
        else:
            start = target * rng.uniform(0.5, 1.5, len(target))
        x = d.solve(start, target)
        if x is not None:
            candidates.append(x)
    for x in candidates:
        if not any(np.allclose(x, y, rtol=1e-4, atol=1e-6 * scale) for y in found):
            found.append(x)
    return len(found)


def _run_length(d) -> float:
    """Long enough to see the dynamics near the fixed point: 20 periods of its
    fastest rotation and 30 growth times of its instability, at least 200 and
    at most 20,000 time units."""
    x = d.fixed_point()
    if x is None:
        return 200.0
    ev = _spectrum(d, x)
    omega = float(np.abs(ev.imag).max()) if len(ev) else 0.0
    growth = float(ev.real.max()) if len(ev) else 0.0
    t = max(200.0, 20 * 2 * np.pi / omega if omega > 1e-9 else 0.0, 30 / growth if growth > 1e-9 else 0.0)
    return float(min(t, 20000.0))


@register("oscillation", "F", needs="K", cost="moderate", limit=400)
def oscillation(ctx) -> dict[str, Any] | None:
    """Whether the rate equations settle into sustained oscillation from the
    initial state, and the period: at least three peaks of an active species in
    the second half of a run long enough for the fixed point's own time scales,
    with a swing above 0.1% of its level."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    from chemart.simulate import integrate

    d = _dynamics(ctx)
    t_end = _run_length(d)
    t = np.linspace(0.0, t_end, 4001)
    traj, times = integrate(ctx.net, t_end, t_eval=t)
    late = times > t_end / 2
    best = None
    for i in d.active:
        x = traj[d.ids[i]][late]
        swing = x.max() - x.min()
        if swing <= 1e-3 * max(abs(x.mean()), 1e-12):
            continue
        peaks = np.flatnonzero((x[1:-1] > x[:-2]) & (x[1:-1] >= x[2:]) & (x[1:-1] > x.mean())) + 1
        if len(peaks) >= 3 and (best is None or swing > best[0]):
            best = (swing, float(np.diff(times[late][peaks]).mean()))
    return {"oscillates": best is not None, "period": best[1] if best else None}


def _rates_at(ctx, d, x):
    from chemart.simulate import _Rates

    return _Rates(ctx.net)(d.full(x))


@register("flux_concentration", "F", needs="K", cost="moderate", limit=2000)
def flux_concentration(ctx) -> float | None:
    """Gini coefficient of the reaction rates at the fixed point: whether a few
    reactions carry most of the flux."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    from chemart.measures.graph import gini

    d = _dynamics(ctx)
    x = d.fixed_point()
    return None if x is None else gini(np.abs(_rates_at(ctx, d, x)))


@register("entropy_production", "F", needs="K", cost="moderate", limit=2000)
def entropy_production(ctx) -> float | None:
    """Σ (J+ - J-) ln(J+/J-) over the reversible pairs at the fixed point (in
    units of the gas constant times temperature): 0 at detailed balance,
    positive when the network runs driven, away from equilibrium."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    pairs = reversible_pairs(ctx.net)
    d = _dynamics(ctx)
    x = d.fixed_point() if pairs else None
    if x is None:
        return None
    J = _rates_at(ctx, d, x)
    total = 0.0
    for f, b in pairs:
        if J[f] > 0 and J[b] > 0:
            total += (J[f] - J[b]) * np.log(J[f] / J[b])
    return float(total)


@register("sloppiness", "F", needs="K", cost="moderate", limit=60)
def sloppiness(ctx, t_end: float = 20.0, points: int = 21) -> float | None:
    """Orders of magnitude spanned by the eigenvalues of the Fisher information
    of the log mass-action constants, from the sensitivity of the active
    species' trajectories (Gutenkunst et al. 2007): large means a few parameter
    combinations matter and most barely do."""
    if not _needs_state(ctx) or not _dynamics(ctx).active:
        return None
    from chemart.network import Network
    from chemart.simulate import integrate

    d = _dynamics(ctx)
    t = np.linspace(0.0, t_end, points)
    idx = [j for j, r in enumerate(ctx.net.reactions) if r.rate["law"] == "mass-action" and r.rate["k"] > 0]
    if not idx:
        return None

    def run(net):
        traj, _ = integrate(net, t_end, t_eval=t)
        return np.concatenate([traj[d.ids[i]] for i in d.active])

    base = run(ctx.net)
    J = []
    for j in idx:
        data = ctx.net.to_dict()
        data["reactions"][j]["rate"] = dict(data["reactions"][j]["rate"], k=data["reactions"][j]["rate"]["k"] * np.exp(1e-4))
        J.append((run(Network.from_dict(data)) - base) / 1e-4)
    ev = np.linalg.eigvalsh(np.array(J) @ np.array(J).T)
    ev = ev[ev > 1e-12 * ev.max()] if ev.max() > 0 else ev[:0]
    return float(np.log10(ev.max() / ev.min())) if len(ev) else None
