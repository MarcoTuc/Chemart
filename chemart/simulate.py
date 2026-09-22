"""Simulate a reaction network: rate equations (ODE) or Gillespie's SSA.

    from chemart import simulate
    net = chemart.generate_network("brusselator", seed=1)
    traj = simulate.ode(net, t_end=40)                    # deterministic
    traj = simulate.ssa(net, t_end=40, volume=100, seed=0)  # one stochastic path

Both return a `chemart.Trajectory`. Networks without rates, or with rates to
change, get them from `assign` (also reachable as the `rates=` and `x0=`
arguments of `ode` and `ssa`): a number, a distribution such as
``{"dist": "lognormal", "mean": 0, "sigma": 1}``, a table keyed by reaction text
or index (or by species for x0), a .json or .csv file, or a function
``(reaction, rng) -> rate``.

Supported rate laws: mass-action, power, michaelis-menten, hill, saturating,
and arrhenius when a temperature and gas constant are known (on the rate dict
as "T" and "R", or passed as `temperature=` and `gas_constant=`). Flows
(`inflow`, `outflow`, the "constant-total" dilution) and `extras["buffered"]`
species are handled by both simulators.

Units: amounts in `initial_state` and in trajectory frames are per unit volume.
The SSA works on molecule counts, count = amount * volume * avogadro, and
reports amounts again; with the defaults volume = avogadro = 1 they are counts.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from chemart import kinetics
from chemart.network import CONSTANT_TOTAL, Network, Reaction
from chemart.trajectory import Frame, Trajectory

DISTRIBUTIONS = {
    "constant": ("value",),
    "uniform": ("low", "high"),
    "loguniform": ("low", "high"),
    "lognormal": ("mean", "sigma"),      # of the underlying normal (numpy's convention)
    "exponential": ("scale",),
    "gamma": ("shape", "scale"),
}


class NotSimulable(Exception):
    """The network cannot be simulated as it stands; the message says why."""


# --------------------------------------------------------------------------
# Checks and rate laws
# --------------------------------------------------------------------------

def check(net: Network, *, temperature: float | None = None, gas_constant: float | None = None) -> None:
    """Fail early and specifically, rather than silently producing a flat line."""
    name = net.chemistry or "network"
    if not net.reactions:
        raise NotSimulable(f"{name} has no reactions; it defines an interaction law instead "
                           "(see extras['interaction_law'])")
    unrated = [i for i, r in enumerate(net.reactions) if r.rate is None]
    if len(unrated) == len(net.reactions):
        raise NotSimulable(f"{name} carries no rate constants; give it some with "
                           "rates=... (a number, a distribution, a table or a file)")
    if unrated:
        raise NotSimulable(f"{len(unrated)} of {len(net.reactions)} reactions have no rate "
                           f"(indices {unrated[:8]}{'...' if len(unrated) > 8 else ''}); rate "
                           "them with assign(..., fill_only=True) or pass rates=")
    for r in net.reactions:
        if r.rate["law"] == "arrhenius":
            _arrhenius_k(r.rate, temperature, gas_constant)


def _arrhenius_k(rate: dict, temperature, gas_constant) -> float:
    T = rate.get("T", temperature)
    R = rate.get("R", gas_constant)
    if T is None or R is None:
        raise NotSimulable("arrhenius rates need a temperature and a gas constant in the units "
                           "of Ea: put 'T' and 'R' on the rate, or pass temperature= and gas_constant=")
    return float(rate["A"]) * float(np.exp(-float(rate["Ea"]) / (float(R) * float(T))))


class _Rates:
    """Reaction rates v(x) for concentrations x, vectorised for mass action."""

    def __init__(self, net: Network, temperature=None, gas_constant=None):
        ids, R, _ = net.matrices()
        self.index = {s: i for i, s in enumerate(ids)}
        self.n = len(ids)
        self.reactants = [[(self.index[s], n) for s, n in r.reactants.items()] for r in net.reactions]
        self.rates = [r.rate for r in net.reactions]
        mass = [j for j, r in enumerate(self.rates) if r["law"] in ("mass-action", "arrhenius")]
        self.mass = np.array(mass, dtype=int)
        self.k = np.array([
            _arrhenius_k(self.rates[j], temperature, gas_constant) if self.rates[j]["law"] == "arrhenius"
            else float(self.rates[j]["k"]) for j in mass
        ])
        self.R_mass = R[:, mass].T.tocsr() if mass else None
        self.other = [j for j in range(len(self.rates)) if j not in set(mass)]

    def __call__(self, x: np.ndarray) -> np.ndarray:
        v = np.empty(len(self.rates))
        xp = np.maximum(x, 0.0)
        if self.R_mass is not None:
            with np.errstate(divide="ignore"):
                v[self.mass] = self.k * np.exp(self.R_mass @ np.log(xp))
        for j in self.other:
            v[j] = self._law(self.rates[j], self.reactants[j], xp)
        return v

    def _law(self, rate, reac, x) -> float:
        law = rate["law"]
        if law == "saturating":
            return rate["k"] * np.prod([(x[i] / (1.0 + x[i] / rate["K"])) ** n for i, n in reac])
        if law == "hill":
            # The regulator is usually not a reactant, so it is named on the rate dict.
            p = x[self.index[rate["regulator"]]]
            h = p ** rate["n"] / (rate["K"] ** rate["n"] + p ** rate["n"])
            return rate["vmax"] * (h if rate.get("mode", "activation") == "activation" else 1.0 - h)
        (i, _), = reac                              # power and michaelis-menten: one substrate
        if law == "power":
            return rate["k"] * x[i] ** rate["order"]
        if law == "michaelis-menten":
            return rate["vmax"] * x[i] / (rate["km"] + x[i])
        raise NotSimulable(f"rate law {law!r} is not supported")

    def jacobian(self, x: np.ndarray) -> np.ndarray | None:
        """dv/dx for all-mass-action networks; None when another law is present."""
        if self.other:
            return None
        xp = np.maximum(x, 0.0)
        D = np.zeros((len(self.rates), self.n))
        for row, j in enumerate(self.mass):
            reac = self.reactants[j]
            for i, n in reac:
                rest = np.prod([xp[l] ** m for l, m in reac if l != i])
                D[j, i] = self.k[row] * n * xp[i] ** (n - 1) * rest
        return D


def _flows(net: Network, index: dict[str, int]) -> tuple[np.ndarray, np.ndarray]:
    influx = np.zeros(len(index))
    for s, value in (net.inflow or {}).items():
        influx[index[s]] = value
    efflux = np.zeros(len(index))
    if isinstance(net.outflow, dict):
        for s, value in net.outflow.items():
            efflux[index[s]] = value
    elif isinstance(net.outflow, (int, float)) and not isinstance(net.outflow, bool):
        efflux[:] = float(net.outflow)
    return influx, efflux


def rhs(net: Network, *, temperature=None, gas_constant=None):
    """Return (species_ids, f(t, x)) for the network's mass-balance ODEs."""
    ids, R, P = net.matrices()
    S = (P - R).tocsr().astype(float)
    v = _Rates(net, temperature, gas_constant)
    influx, efflux = _flows(net, v.index)
    buffered = [v.index[s] for s in net.extras.get("buffered", []) if s in v.index]
    dilute = net.outflow == CONSTANT_TOTAL

    def f(_t, x):
        dx = S @ v(x) + influx - efflux * x
        if dilute:
            total = x.sum()
            if total > 0:
                dx = dx - x * dx.sum() / total
        dx[buffered] = 0.0
        return dx

    return ids, f


def jacobian(net: Network, *, temperature=None, gas_constant=None):
    """The analytic Jacobian J(t, x) of `rhs` for all-mass-action networks without
    the constant-total dilution, or None (the solver then estimates it)."""
    if net.outflow == CONSTANT_TOTAL:
        return None
    ids, R, P = net.matrices()
    S = (P - R).toarray().astype(float)
    v = _Rates(net, temperature, gas_constant)
    if v.other:
        return None
    _, efflux = _flows(net, v.index)
    buffered = [v.index[s] for s in net.extras.get("buffered", []) if s in v.index]

    def J(_t, x):
        out = S @ v.jacobian(x) - np.diag(efflux)
        out[buffered, :] = 0.0
        return out

    return J


# --------------------------------------------------------------------------
# Assigning rates and initial states
# --------------------------------------------------------------------------

def reaction_key(r: Reaction) -> str:
    """The text a rate table uses for a reaction: 'A + B -> 2 C'."""
    side = lambda d: " + ".join(s if n == 1 else f"{n} {s}" for s, n in d.items()) or "∅"  # noqa: E731
    return f"{side(r.reactants)} -> {side(r.products)}"


def draw(spec: dict, rng: np.random.Generator) -> float:
    """One value from a distribution spec such as {"dist": "uniform", "low": 0, "high": 1}."""
    name = spec.get("dist")
    if name not in DISTRIBUTIONS:
        raise ValueError(f"unknown distribution {name!r}; choose one of {sorted(DISTRIBUTIONS)}")
    missing = [k for k in DISTRIBUTIONS[name] if k not in spec]
    if missing:
        raise ValueError(f"{name} distribution needs {missing}")
    if name == "constant":
        return float(spec["value"])
    if name == "uniform":
        return float(rng.uniform(spec["low"], spec["high"]))
    if name == "loguniform":
        return float(np.exp(rng.uniform(np.log(spec["low"]), np.log(spec["high"]))))
    if name == "lognormal":
        return float(rng.lognormal(spec["mean"], spec["sigma"]))
    if name == "exponential":
        return float(rng.exponential(spec["scale"]))
    return float(rng.gamma(spec["shape"], spec["scale"]))


def _load(spec: Any) -> Any:
    """A spec given as a path is read from a .json or a two-column .csv file."""
    if not isinstance(spec, (str, Path)):
        return spec
    path = Path(spec)
    if path.suffix == ".json":
        return json.loads(path.read_text())
    if path.suffix == ".csv":
        table = {}
        with path.open(newline="") as fh:
            for row in csv.reader(fh):
                if len(row) < 2 or not row[0].strip():
                    continue
                try:
                    table[row[0].strip()] = float(row[1])
                except ValueError:
                    continue                     # a header line
        return table
    raise ValueError(f"{path}: rate and state files must be .json or .csv")


def _value(spec: Any, keys: list, subject: Any, rng) -> Any:
    """Resolve a spec for one reaction or species; None keeps what is there."""
    if callable(spec):
        return spec(subject, rng)
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        return float(spec)
    if isinstance(spec, list):
        i = keys[0]
        return _value(spec[i], keys, subject, rng) if i < len(spec) else None
    if isinstance(spec, dict):
        if "dist" in spec:
            return draw(spec, rng)
        if "law" in spec:
            return dict(spec)
        for key in [*keys, *(str(k) for k in keys), "*"]:
            if key in spec:
                return _value(spec[key], keys, subject, rng)
        return None
    raise ValueError(f"cannot read a rate or state from {spec!r}")


def assign(net: Network, rates: Any = None, x0: Any = None, *,
           rng: np.random.Generator | None = None, fill_only: bool = False) -> Network:
    """A copy of `net` with rates and/or an initial state assigned.

    `rates` and `x0` take: a number (a mass-action k, or an amount for every
    species), a distribution spec, a dict keyed by reaction index or text
    ('A + B -> C') for rates or by species id for x0 ("*" for the rest), a list
    by index, a .json or .csv file, or a function (reaction or species id, rng).
    A rate value may be a number (mass action) or a full rate dict. With
    `fill_only`, only reactions without a rate are rated.
    """
    rng = rng if rng is not None else np.random.default_rng()
    out = Network.from_dict(json.loads(json.dumps(net.to_dict())))
    if rates is not None:
        spec = _load(rates)
        for j, r in enumerate(out.reactions):
            if fill_only and r.rate is not None:
                continue
            value = _value(spec, [j, reaction_key(r)], r, rng)
            if value is None:
                continue
            rate = {"law": "mass-action", "k": float(value)} if not isinstance(value, dict) else value
            problem = kinetics.rate_problem(rate)
            if problem:
                raise ValueError(f"reaction {j} ({reaction_key(r)}): {problem}")
            r.rate = rate
    if x0 is not None:
        spec = _load(x0)
        state = dict(out.initial_state or {})
        for s in out.species:
            value = _value(spec, [s.id], s.id, rng)
            if value is not None:
                state[s.id] = float(value)
        out.initial_state = {s: v for s, v in state.items() if v}
    Network.__post_init__(out)
    return out


def _prepare(net, rates, x0, rng, fill_only=False):
    if rates is not None or x0 is not None:
        net = assign(net, rates, x0, rng=rng, fill_only=fill_only)
    if net.initial_state is None:
        raise NotSimulable("no initial state: the chemistry prescribes none, so pass x0= "
                           "(a number for every species, a table, a distribution or a file)")
    return net


# --------------------------------------------------------------------------
# Deterministic: rate equations
# --------------------------------------------------------------------------

def integrate(net: Network, t_end: float, x0: dict | None = None, t_eval=None,
              method: str = "LSODA", *, temperature=None, gas_constant=None):
    """Integrate the rate equations to t_end; return ({species: array}, times).

    The low-level form used by the tests; `ode` wraps it in a Trajectory.
    """
    from scipy.integrate import solve_ivp

    check(net, temperature=temperature, gas_constant=gas_constant)
    ids, f = rhs(net, temperature=temperature, gas_constant=gas_constant)
    start = x0 if x0 is not None else (net.initial_state or {})
    x = np.array([float(start.get(s, 0.0)) for s in ids])
    jac = jacobian(net, temperature=temperature, gas_constant=gas_constant) if method != "RK45" else None
    sol = solve_ivp(f, (0.0, t_end), x, method=method, t_eval=t_eval, jac=jac, rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise NotSimulable(f"integration failed: {sol.message}; stiff systems may need "
                           "solver='Radau' or 'BDF'")
    return {s: sol.y[i] for i, s in enumerate(ids)}, sol.t


def ode(net: Network, t_end: float, *, rates: Any = None, x0: Any = None, points: int = 200,
        solver: str = "LSODA", seed: int | None = None, temperature=None, gas_constant=None,
        fill_only: bool = False) -> Trajectory:
    """Integrate the network's rate equations from its initial state to t_end.

    `rates` and `x0` are assigned first (see `assign`); `seed` drives any
    distribution in them. Frames are `points` evenly spaced times.
    """
    rng = np.random.default_rng(seed)
    net = _prepare(net, rates, x0, rng, fill_only)
    t_eval = np.linspace(0.0, t_end, max(points, 2))
    traj, t = integrate(net, t_end, t_eval=t_eval, method=solver,
                        temperature=temperature, gas_constant=gas_constant)
    ids = list(traj)
    frames = [Frame(t=float(t[k]), state={s: float(traj[s][k]) for s in ids if traj[s][k] != 0.0})
              for k in range(len(t))]
    return Trajectory(network=net, frames=frames, method="ode", clock="time", settings={
        "t_end": t_end, "points": points, "solver": solver, "seed": seed,
        "rates": _describe(rates), "x0": _describe(x0),
    })


def _describe(spec: Any) -> Any:
    """A spec as it can be stored in the trajectory's settings (JSON)."""
    if isinstance(spec, Path):
        return str(spec)
    if spec is None or isinstance(spec, (int, float, str, list, dict)):
        return spec
    return repr(spec)


# --------------------------------------------------------------------------
# Stochastic: Gillespie's direct method
# --------------------------------------------------------------------------

def ssa(net: Network, t_end: float, *, rates: Any = None, x0: Any = None, volume: float = 1.0,
        avogadro: float = 1.0, points: int = 200, seed: int | None = None,
        max_events: int = 1_000_000, temperature=None, gas_constant=None,
        fill_only: bool = False) -> Trajectory:
    """One stochastic path of the network, sampled at `points` evenly spaced times.

    Mass-action propensities use c = kinetics.k_to_c(k, reactants, volume,
    avogadro) and a = c * prod_i C(x_i, n_i); other laws use a = Ω f(x/Ω), with
    Ω = volume * avogadro. Inflow is a zeroth-order source, outflow a
    first-order sink; with the constant-total dilution a random molecule is
    removed whenever an event raises the total above its starting value.
    Stops after `max_events` events (settings["stopped"] then says so).
    """
    from scipy.special import comb

    rng = np.random.default_rng(seed)
    net = _prepare(net, rates, x0, rng, fill_only)
    check(net, temperature=temperature, gas_constant=gas_constant)
    omega = volume * avogadro
    v = _Rates(net, temperature, gas_constant)
    ids = [s.id for s in net.species]
    n_species = len(ids)
    buffered = np.zeros(n_species, dtype=bool)
    for s in net.extras.get("buffered", []):
        if s in v.index:
            buffered[v.index[s]] = True
    x = np.array([np.rint(float((net.initial_state or {}).get(s, 0.0)) * omega) for s in ids])

    # Every event channel: the reactions, then inflow sources and outflow sinks.
    changes, labels = [], []
    for r in net.reactions:
        delta = {}
        for s, n in r.reactants.items():
            delta[v.index[s]] = delta.get(v.index[s], 0) - n
        for s, n in r.products.items():
            delta[v.index[s]] = delta.get(v.index[s], 0) + n
        changes.append([(i, d) for i, d in delta.items() if d and not buffered[i]])
        labels.append([[s for s, n in r.reactants.items() for _ in range(n)],
                       [s for s, n in r.products.items() for _ in range(n)]])
    influx, efflux = _flows(net, v.index)
    sources = [i for i in range(n_species) if influx[i] > 0]
    sinks = [i for i in range(n_species) if efflux[i] > 0]
    for i in sources:
        changes.append([] if buffered[i] else [(i, 1)])
        labels.append([[], [ids[i]]])
    for i in sinks:
        changes.append([] if buffered[i] else [(i, -1)])
        labels.append([[ids[i]], []])

    mass = v.mass
    c = np.array([kinetics.k_to_c(v.k[row], net.reactions[j].reactants, volume, avogadro)
                  for row, j in enumerate(mass)])
    entries = [(row, i, n) for row, j in enumerate(mass) for i, n in v.reactants[j]]
    e_row = np.array([e[0] for e in entries], dtype=int)
    e_species = np.array([e[1] for e in entries], dtype=int)
    e_n = np.array([e[2] for e in entries], dtype=float)
    n_channels = len(changes)
    n_reactions = len(net.reactions)
    dilute = net.outflow == CONSTANT_TOTAL
    free = ~buffered
    target = x[free].sum()

    def propensities(x):
        a = np.zeros(n_channels)
        if len(mass):
            with np.errstate(divide="ignore"):
                logs = np.log(comb(x[e_species], e_n))
            a[mass] = c * np.exp(np.bincount(e_row, weights=logs, minlength=len(mass)))
        if v.other:
            conc = v(x / omega)
            a[v.other] = omega * conc[v.other]
        if sources:
            a[n_reactions:n_reactions + len(sources)] = influx[sources] * omega
        if sinks:
            a[n_reactions + len(sources):] = efflux[sinks] * x[sinks]
        return a

    samples = np.linspace(0.0, t_end, max(points, 2))
    frames: list[Frame] = []
    window = np.zeros(n_channels, dtype=np.int64)
    k, t, events, stopped = 0, 0.0, 0, None

    def record(time):
        fired = [[*labels[ch], int(window[ch])] for ch in np.flatnonzero(window)]
        frames.append(Frame(t=float(time), state={ids[i]: float(x[i] / omega) for i in np.flatnonzero(x)},
                            fired=fired))
        window[:] = 0

    while True:
        a = propensities(x)
        a0 = a.sum()
        t_next = t + rng.exponential(1.0 / a0) if a0 > 0 else np.inf
        while k < len(samples) and samples[k] < t_next:
            record(samples[k])
            k += 1
        if k == len(samples):
            break
        if events >= max_events:
            stopped = f"max_events={max_events} reached at t={t:.6g}"
            break
        ch = int(np.searchsorted(np.cumsum(a), rng.random() * a0, side="right"))
        ch = min(ch, n_channels - 1)
        for i, d in changes[ch]:
            x[i] += d
        window[ch] += 1
        events += 1
        if dilute:
            while x[free].sum() > target:
                pool = np.where(free, x, 0.0)
                i = int(np.searchsorted(np.cumsum(pool), rng.random() * pool.sum(), side="right"))
                x[min(i, n_species - 1)] -= 1
        t = t_next

    if stopped:
        record(t)
    return Trajectory(network=net, frames=frames, method="ssa", clock="time", settings={
        "t_end": t_end, "points": points, "volume": volume, "avogadro": avogadro,
        "seed": seed, "events": events, "stopped": stopped,
        "rates": _describe(rates), "x0": _describe(x0),
    })
