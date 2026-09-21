# Simulating dynamics

Chemart ships **no simulator**. It hands you a correct, sourced network; what
you integrate it with is your choice. This file covers the bundled ODE script
for quick work, and what to do when you outgrow it.

## Contents
- [Which networks can be simulated at all](#which-networks-can-be-simulated-at-all)
- [The bundled script](#the-bundled-script)
- [Using it as a library](#using-it-as-a-library)
- [What it handles](#what-it-handles)
- [Pitfalls](#pitfalls)
- [Going further: SSA, and real solvers](#going-further-ssa-and-real-solvers)

## Which networks can be simulated at all

You need rate constants. Many entries provide them; the rest are topology
or structure only, by nature of the model — the Chemical Casting Model has no
rate constants because it is not that kind of object.

```python
import chemart
net = chemart.generate_network("brusselator", seed=1)
"rate-constants" in net.provides          # True -> integrable
all(r.rate for r in net.reactions)        # check every reaction, not just some
```

A network can be partially rated. Reactions with `rate=None` have no
propensity and the integrator cannot fire them; treat a mixed network as a
modelling question, not a bug.

**These integrate cleanly** at default parameters (fully rated, fixed
species set, no structure to lose): `analog-function-crn`, `bagley-farmer`, `bigan-conservative-crn`, `brusselator`,
`chameleon`, `farmer-immune`, `hill-kinetics`, `jain-krishna`, `kappa-calculus`,
`matrix-chemistry`, `mcs-bl`, `mechanical-self-assembly`,
`metabolic-robot-controller`, `michaelis-menten`, `okamoto-switch`,
`random-catalytic-networks`, `repressilator`, `synthon`.

Some more are rated but carry structure a well-mixed integrator throws away:
`chemoton`, `disperser`, `gard` (compartments) and `flow-ac`, `oregonator`
(space). They will run; the result is a mean-field
approximation you chose, not the published model. Say so if you report it.

Note that `michaelis-menten` and `hill-kinetics` emit **elementary mass
action** (`E + S ⇌ ES → E + P`), not the abridged laws. No entry emits
`michaelis-menten`, `hill`, `saturating` or `power` at default parameters —
only mass-action and arrhenius.

## The bundled script

```bash
# print a small table of the trajectory
uv run python .claude/skills/chemart/scripts/simulate.py brusselator --t-end 40

# pass chemistry parameters through, and a seed
uv run python .claude/skills/chemart/scripts/simulate.py repressilator \
    --seed 1 --t-end 300 --points 600

# save a PNG (needs matplotlib) and/or CSV
uv run python .claude/skills/chemart/scripts/simulate.py brusselator \
    --t-end 40 --plot bruss.png --csv bruss.csv

# override the initial state
uv run python .claude/skills/chemart/scripts/simulate.py brusselator \
    --t-end 40 --x0 X=1.0 --x0 Y=1.0
```

If a chemistry has no usable kinetics the script says so and names the
offending reactions instead of silently producing a flat line.

## Using it as a library

```python
import sys
sys.path.insert(0, "<repo>/.claude/skills/chemart/scripts")   # absolute path
from simulate import integrate

import chemart
net = chemart.generate_network("brusselator", seed=1)
traj, t = integrate(net, t_end=40.0, points=400)
traj["X"]        # numpy array aligned with t
```

Use an **absolute** path, or run from the repository root where the relative
`.claude/skills/chemart/scripts` resolves. A relative path from any other
working directory fails with `ModuleNotFoundError: No module named 'simulate'`,
which is a confusing way to learn this. Inside a script, derive it from the
file instead — `examples/dynamics.py` shows the robust form:

```python
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
```

`integrate(net, t_end, x0=None, points=None, method="LSODA")` returns
`(dict_of_species_arrays, times)`. `x0` defaults to `net.initial_state`.

## What it handles

- the rate laws `mass-action`, `power`, `michaelis-menten`, `hill`,
  `saturating` (see `network-record.md` for their exact parameter names)
- `inflow` as a constant influx term
- `outflow` as a first-order removal, per-species dict or single number
- `outflow == "constant-total"`, the book's non-selective dilution flux
- `extras["buffered"]`, species pinned at their initial value

It does **not** handle `arrhenius` — it raises with a clear message. That law
needs a temperature and gas constant that are not part of the rate dict, which
is an open question in the repository's `to_decide.md`.

Because flows and buffering are applied for you, do not add dilution or
clamping terms yourself on top; you will double-count them.

## Pitfalls

**Buffered species look broken but aren't.** In the Brusselator, `A` and `B`
are buffered, so their trajectories are flat lines. That is the model: they are
reservoirs. Check `net.extras.get("buffered")` before reporting "the
concentration never changes".

**`observed` networks are samples, not systems.** Integrating an `observed`
network treats one run's firing record as if it were the full rate system. It
will run, and the result is usually meaningless. Check `net.status` first.

**Defaults are small.** A default network may be too small to show the
published behaviour. `describe_chemistry(id)["params"]` gives each parameter's
description including a `range:` note with the paper's scale.

**Stiffness.** Several of these systems are stiff; the default `LSODA` copes.
If a run fails or crawls, try `method="Radau"` or `"BDF"` before concluding the
network is wrong.

**Units are the model's own.** Chemart does not impose units. Two chemistries'
rate constants are not comparable just because both are called `k`.

## Going further: SSA, and real solvers

For stochastic simulation, the record gives you everything Gillespie needs —
`matrices()` for stoichiometry, `r.rate["k"]` for propensities, and
`initial_state` for counts. Remember the mass-action → stochastic conversion:
`chemart.kinetics.k_to_c(k, reactants, volume, avogadro)` implements the
standard correction, including the combinatorial factor for homodimers
(`2A → …` gains a factor 2).

For anything substantial, export and use a real tool:

```python
d = net.to_dict()        # plain JSON: species, reactions, rates, flows, extras
```

From there it is a short hop to SBML, Antimony, Tellurium, COPASI, or a
hand-written SSA. The record is deliberately simple so this conversion is
boring — reactant/product multiplicities are already separate, which is the
part that usually goes wrong.
