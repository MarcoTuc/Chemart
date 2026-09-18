#!/usr/bin/env python3
"""Integrate a Chemart network as ODEs.

Chemart ships no simulator on purpose: its job is to hand you a correct
network. This is the supported way to get a quick trajectory out of one.
(The repository's own ODE helper lives under tests/ and is not importable as
part of the package, so this script stands alone.)

    uv run python .claude/skills/chemart/scripts/simulate.py brusselator --t-end 40
    uv run python .claude/skills/chemart/scripts/simulate.py repressilator \
        --seed 1 --t-end 300 --plot repr.png

As a library:

    from simulate import integrate
    traj, t = integrate(net, t_end=40.0, points=400)

Handles mass-action, power, michaelis-menten, hill and saturating rate laws,
inflow/outflow, the "constant-total" dilution flux, and extras["buffered"].
Raises a clear error for arrhenius, which needs a temperature the rate dict
does not carry.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp

from chemart.network import CONSTANT_TOTAL, Network


class NotIntegrable(Exception):
    """The network cannot be integrated as written, with a reason that names names."""


def check(net: Network) -> None:
    """Fail early and specifically, rather than silently producing a flat line."""
    if not net.reactions:
        raise NotIntegrable(
            f"{net.chemistry or 'network'} has no reactions; "
            "it defines an interaction law instead — see extras['interaction_law']"
        )
    unrated = [i for i, r in enumerate(net.reactions) if r.rate is None]
    if len(unrated) == len(net.reactions):
        raise NotIntegrable(
            f"{net.chemistry or 'network'} carries no rate constants "
            f"(provides: {', '.join(net.provides)}); it is a topology-only chemistry. "
            "Filter on 'rate-constants' to find integrable ones."
        )
    if unrated:
        raise NotIntegrable(
            f"{len(unrated)} of {len(net.reactions)} reactions have no rate "
            f"(indices {unrated[:8]}{'…' if len(unrated) > 8 else ''}). "
            "Partially rated networks need a modelling decision, not a default."
        )
    laws = {r.rate["law"] for r in net.reactions}
    if "arrhenius" in laws:
        raise NotIntegrable(
            "arrhenius rates need a temperature and gas constant that are not part "
            "of the law (chemart.kinetics is unit-agnostic). Entries using it "
            "validate their kinetics in closed form; see to_decide.md."
        )


def rhs(net: Network):
    """Return (species_ids, f(t, x)) for the network's mass-balance ODEs."""
    ids, R, P = net.matrices()
    index = {s: i for i, s in enumerate(ids)}
    S = (P - R).toarray().astype(float)
    reactants = [[(index[s], n) for s, n in r.reactants.items()] for r in net.reactions]
    rates = [r.rate for r in net.reactions]
    buffered = [index[s] for s in net.extras.get("buffered", []) if s in index]

    influx = np.zeros(len(ids))
    for s, value in (net.inflow or {}).items():
        influx[index[s]] = value

    efflux = np.zeros(len(ids))
    if isinstance(net.outflow, dict):
        for s, value in net.outflow.items():
            efflux[index[s]] = value
    elif isinstance(net.outflow, (int, float)) and not isinstance(net.outflow, bool):
        efflux[:] = float(net.outflow)

    def propensity(rate, reac, x):
        law = rate["law"]
        if law == "mass-action":
            return rate["k"] * np.prod([max(x[i], 0.0) ** n for i, n in reac])
        if law == "saturating":
            return rate["k"] * np.prod(
                [(max(x[i], 0.0) / (1.0 + max(x[i], 0.0) / rate["K"])) ** n for i, n in reac]
            )
        (i, _), = reac                      # the remaining laws are single-substrate
        xi = max(x[i], 0.0)
        if law == "power":
            return rate["k"] * xi ** rate["order"]
        if law == "michaelis-menten":
            return rate["vmax"] * xi / (rate["km"] + xi)
        raise NotIntegrable(f"rate law {law!r} is not supported by this script")

    def hill(rate, x):
        # The regulator is usually not a reactant, so it is named on the rate dict.
        p = max(x[index[rate["regulator"]]], 0.0)
        h = p ** rate["n"] / (rate["K"] ** rate["n"] + p ** rate["n"])
        return rate["vmax"] * (h if rate.get("mode", "activation") == "activation" else 1.0 - h)

    def f(_t, x):
        v = np.array([
            hill(rate, x) if rate["law"] == "hill" else propensity(rate, reac, x)
            for rate, reac in zip(rates, reactants)
        ])
        dx = S @ v + influx - efflux * x
        if net.outflow == CONSTANT_TOTAL:
            total = x.sum()
            if total > 0:
                dx = dx - x * dx.sum() / total
        dx[buffered] = 0.0
        return dx

    return ids, f


def integrate(net: Network, t_end: float, x0=None, points: int | None = None,
              method: str = "LSODA"):
    """Integrate to t_end. Returns ({species_id: array}, times)."""
    check(net)
    ids, f = rhs(net)
    start = x0 if x0 is not None else (net.initial_state or {})
    if not start:
        raise NotIntegrable(
            "no initial state: the chemistry prescribes none, so pass x0="
            "{'X': 1.0, ...} explicitly."
        )
    x = np.array([float(start.get(s, 0.0)) for s in ids])
    t_eval = np.linspace(0.0, t_end, points) if points else None
    sol = solve_ivp(f, (0.0, t_end), x, method=method, t_eval=t_eval, rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise NotIntegrable(
            f"integration failed: {sol.message}. These systems are often stiff — "
            "try method='Radau' or 'BDF'."
        )
    return {s: sol.y[i] for i, s in enumerate(ids)}, sol.t


def _param(text: str):
    name, sep, raw = text.partition("=")
    if not sep:
        raise argparse.ArgumentTypeError(f"expected NAME=VALUE, got {text!r}")
    try:
        return name, json.loads(raw)
    except json.JSONDecodeError:
        return name, raw


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("chemistry")
    ap.add_argument("--t-end", type=float, default=40.0)
    ap.add_argument("--points", type=int, default=200, help="output points (0 = solver's own)")
    ap.add_argument("--seed", type=int)
    ap.add_argument("-p", "--param", action="append", type=_param, default=[],
                    metavar="NAME=VALUE", help="chemistry parameter")
    ap.add_argument("--x0", action="append", type=_param, default=[], metavar="SPECIES=VALUE",
                    help="override an initial concentration")
    ap.add_argument("--method", default="LSODA")
    ap.add_argument("--csv", help="write the trajectory here")
    ap.add_argument("--plot", help="write a PNG here (needs matplotlib)")
    ap.add_argument("--species", nargs="*", help="restrict printing/plotting to these")
    args = ap.parse_args(argv)

    import chemart

    net = chemart.generate_network(args.chemistry, args.seed, **dict(args.param))
    if net.status == "observed":
        print(f"warning: {net.chemistry} is an 'observed' network — these are the reactions "
              f"that fired in one run, not a rate system. Integrating it is rarely meaningful.",
              file=sys.stderr)

    x0 = dict(net.initial_state or {})
    x0.update({k: float(v) for k, v in args.x0})

    try:
        traj, t = integrate(net, args.t_end, x0=x0 or None,
                            points=args.points or None, method=args.method)
    except NotIntegrable as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    names = args.species or list(traj)
    buffered = set(net.extras.get("buffered", []))

    print(net.summary())
    if buffered:
        print(f"buffered (held constant by the model): {', '.join(sorted(buffered))}")
    print()
    width = max(len(n) for n in names)
    print(f"{'t':>10}  " + "  ".join(f"{n:>{max(width, 10)}}" for n in names))
    step = max(1, len(t) // 10)
    for i in range(0, len(t), step):
        row = "  ".join(f"{traj[n][i]:>{max(width, 10)}.4g}" for n in names)
        print(f"{t[i]:>10.4g}  {row}")

    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("t," + ",".join(names) + "\n")
            for i in range(len(t)):
                fh.write(f"{t[i]}," + ",".join(str(traj[n][i]) for n in names) + "\n")
        print(f"\nwrote {args.csv}")

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("error: matplotlib is not installed; use --csv instead", file=sys.stderr)
            return 2
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for n in names:
            ax.plot(t, traj[n], label=n + (" (buffered)" if n in buffered else ""))
        ax.set_xlabel("time")
        ax.set_ylabel("concentration")
        ax.set_title(f"{net.chemistry} (seed={net.seed})")
        ax.legend(loc="best", fontsize="small")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"wrote {args.plot}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
