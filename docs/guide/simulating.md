# Simulating dynamics

!!! info "Chemart ships no simulator, on purpose"
    Its job is to hand you a correct, sourced network. What you integrate it
    with is your choice — and for anything substantial you want a real solver,
    not something bundled here.

## Which networks can be simulated at all

You need rate constants. **44** of the 98 entries carry them or a rate law; the
rest are topology or structure only, by the nature of the model. Tierra has no
rate constants because Tierra is not that kind of object.

```python
import chemart

net = chemart.generate_network("brusselator", seed=1)
"rate-constants" in net.provides         # True -> integrable
all(r.rate for r in net.reactions)       # check every reaction, not just some
```

A network can be *partially* rated. Reactions with `rate=None` have no
propensity, so an integrator cannot fire them; treat a mixed network as a
modelling question rather than a bug.

**29 entries** integrate cleanly at default parameters — fully rated, fixed
species set, no spatial or compartment structure to lose:

`analog-function-crn`, `arn`, `bagley-farmer`, `bigan-conservative-crn`,
`brusselator`, `chameleon`, `dimerization`, `ecolab`, `farmer-immune`,
`hill-kinetics`, `jain-krishna`, `kappa-calculus`, `logistic-chemistry`,
`lotka-volterra`, `matrix-chemistry`, `mcs-bl`, `mechanical-self-assembly`,
`metabolic-robot-controller`, `michaelis-menten`, `naming-game-ac`,
`nk-landscape`, `okamoto-switch`, `quasispecies`, `random-catalytic-networks`,
`replication-death`, `replicator-equation`, `repressilator`,
`selection-equation`, `synthon`.

!!! warning "Six more are integrable only as a mean field"
    `chemoton`, `disperser`, `gard` and `isologous-diversification` define
    **compartments**; `flow-ac` and `oregonator` define **space**. Their
    networks carry rate constants, so an integrator will happily run them — and
    silently discard structure that is part of the published model. Treat a
    well-mixed run of these as an approximation you chose, not as the model.

!!! note "`michaelis-menten` and `hill-kinetics` emit elementary mass action"
    Both give the underlying mechanism (`E + S ⇌ ES → E + P`) rather than the
    abridged saturating law, which is the more faithful choice — the abridged
    form has fewer species and is an approximation. In fact **no** entry emits
    `michaelis-menten`, `hill`, `saturating` or `power` at default parameters:
    the only laws that appear are mass-action (39 entries) and arrhenius (2).

## Quick trajectories

The Claude Code skill in this repository bundles a small ODE script for exactly
this, at `.claude/skills/chemart/scripts/simulate.py`:

```bash
uv run python .claude/skills/chemart/scripts/simulate.py brusselator --t-end 40
uv run python .claude/skills/chemart/scripts/simulate.py brusselator \
    --t-end 40 --plot bruss.png --csv bruss.csv
```

It handles mass-action — which is what every rated network actually emits —
plus `power`, `michaelis-menten`, `hill` and `saturating` for completeness,
along with `inflow`/`outflow`, the `constant-total` dilution flux and buffered
species. It refuses clearly, naming the problem, when a chemistry has no usable
kinetics, and it cannot integrate `arrhenius` (see below).

## Writing your own integrator

The record gives you everything directly. The mass-balance ODEs are
`dx/dt = S·v(x)`, where `S = P - R`:

```python
import numpy as np
from scipy.integrate import solve_ivp

ids, R, P = net.matrices()
S = (P - R).toarray().astype(float)
index = {s: i for i, s in enumerate(ids)}

def propensity(r, x):
    # mass action: k * prod(x_i ** n_i)
    return r.rate["k"] * np.prod([x[index[s]] ** n for s, n in r.reactants.items()])

def f(t, x):
    return S @ np.array([propensity(r, x) for r in net.reactions])

x0 = np.array([net.initial_state.get(s, 0.0) for s in ids])
sol = solve_ivp(f, (0, 40), x0, method="LSODA", rtol=1e-8)
```

That covers mass action. For the other laws see
[the rate-law vocabulary](../reference/record.md#rate-laws), and mind the flows
and buffering below.

## Flows and buffered species

Three mechanisms appear in the record, and all three must be applied or your
trajectories will be wrong:

- **`inflow`** — a constant influx per species, added to `dx/dt`.
- **`outflow`** — a first-order removal rate: a per-species dict, a single
  number for all species, or the string `"constant-total"`.
- **`extras["buffered"]`** — species held at constant concentration. Set their
  derivative to zero.

`"constant-total"` is the book's standard flow reactor: a non-selective dilution
flux Φ that holds total concentration constant,

```
dx_i/dt = production_i - x_i * Phi / sum_k(x_k)
```

```python
from chemart.network import CONSTANT_TOTAL
if net.outflow == CONSTANT_TOTAL:
    dx = dx - x * dx.sum() / x.sum()
```

## Pitfalls

!!! warning "Buffered species look broken but aren't"
    In the Brusselator, `A` and `B` are buffered — their trajectories are flat
    lines because they are reservoirs. Check `net.extras.get("buffered")`
    before reporting that a concentration never changes.

!!! warning "`observed` networks are samples, not systems"
    Integrating an `observed` network treats one run's firing record as if it
    were a rate system. It will run; the result is rarely meaningful. Check
    `net.status` first.

**Stiffness.** Several of these systems are stiff. `LSODA` copes; if a run fails
or crawls, try `Radau` or `BDF` before concluding the network is wrong.

**Units are the model's own.** Chemart imposes none. Two chemistries' rate
constants are not comparable just because both are called `k`.

## Stochastic simulation

For Gillespie-style simulation the record has what you need — `matrices()` for
stoichiometry, `initial_state` for counts — but remember the conversion from a
macroscopic rate constant to a mesoscopic one:

```python
from chemart.kinetics import k_to_c

k_to_c(1.0, {"A": 1, "B": 1}, volume=1e-15)   # heterodimer
k_to_c(1.0, {"A": 2},        volume=1e-15)    # homodimer: gains a factor 2
```

The combinatorial factor for identical reactants is the part that is easy to get
wrong, so it is implemented once here.

## Exporting to a real tool

```python
d = net.to_dict()    # plain JSON: species, reactions, rates, flows, extras
```

From there it is a short hop to SBML, Antimony, Tellurium, COPASI or a
hand-written SSA. The record is deliberately plain so this conversion is
boring — reactant and product multiplicities are already separate, which is the
part that usually goes wrong in translation.
