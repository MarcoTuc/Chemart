# Simulating dynamics

A given network, or one built by a generator, is simulated by `chemart.simulate`:
`ode` integrates its rate equations, `ssa` samples one stochastic path with
Gillespie's direct method. Both return a [`Trajectory`](../reference/record.md#the-trajectory-record):
a list of frames, each with the time and the amount of every species.

```python
import chemart
from chemart import simulate

net = chemart.generate_network("brusselator", seed=1)
traj = simulate.ode(net, t_end=40)                          # rate equations
path = simulate.ssa(net, t_end=40, volume=100, seed=0)      # one stochastic path

ids, t, X = traj.array(["X", "Y"])          # dense arrays for plotting
traj.frames[-1].state                       # {species: amount} at the last time
```

From the command line:

```bash
uv run chemart simulate brusselator --t-end 40
uv run chemart simulate brusselator --method ssa --volume 100 --seed 0 --format csv > bruss.csv
```

Turing gases are not simulated this way: their network is not known before the
run. They are *evolved*, and the evolution is itself the simulation.

## Which networks can be simulated

A network needs a rate for every reaction and an initial state. Many entries
carry both; the [Catalog](../catalog/index.md) counts those with their own
kinetics. The rest are topologies, by the nature of the model: Kauffman's
autocatalytic sets or the RAF model say which reactions exist, not how fast
they run.

`simulate` refuses clearly, naming what is missing, rather than returning a
flat line:

```python
net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
simulate.ode(net, t_end=1)
# NotSimulable: no initial state: the chemistry prescribes none, so pass x0= ...
```

These integrate cleanly at default parameters: fully rated, fixed species set,
no spatial or compartment structure to lose.

`analog-function-crn`, `bagley-farmer`, `bigan-conservative-crn`, `brusselator`,
`chameleon`, `farmer-immune`, `hill-kinetics`, `jain-krishna`, `kappa-calculus`,
`matrix-chemistry`, `mcs-bl`, `mechanical-self-assembly`,
`metabolic-robot-controller`, `michaelis-menten`, `okamoto-switch`,
`random-catalytic-networks`, `repressilator`, `synthon`.

!!! warning "Some more are integrable only as a mean field"
    `chemoton`, `disperser` and `gard` define **compartments**; `flow-ac`
    and `oregonator` define **space**. Their networks carry rate constants, so
    the simulator will run them, and silently discard structure that is part
    of the published model. Treat a well-mixed run of these as an
    approximation you chose, not as the model.

!!! note "`michaelis-menten` and `hill-kinetics` emit elementary mass action"
    Both give the underlying mechanism (`E + S ⇌ ES → E + P`) rather than the
    abridged saturating law, which is the more faithful choice. The simulator
    also supports the `power`, `michaelis-menten`, `hill` and `saturating`
    laws for networks you build yourself; `arrhenius` needs a temperature and a
    gas constant in the units of `Ea`, given on the rate dict as `T` and `R`, or
    passed as `temperature=` and `gas_constant=`.

## Rates and initial states: from the chemistry, a distribution, a table or a file

`rates=` and `x0=` give a network the rates and starting amounts it lacks, or
replace the ones it has. Both accept the same kinds of value:

| value | meaning |
|---|---|
| `2.0` | the same mass-action constant for every reaction (or amount for every species) |
| `{"dist": "lognormal", "mean": 0, "sigma": 1}` | one draw per reaction (or species); also `uniform`, `loguniform`, `exponential`, `gamma`, `constant` |
| `{"A + B -> C": 3.0, "*": 0.5}` | a table by reaction text or index (by species id for `x0`); `"*"` covers the rest |
| `[0.1, 0.2, ...]` | a list by reaction index |
| `"rates.csv"`, `"x0.json"` | a file: a two-column CSV (key, value) or a JSON table |
| `lambda reaction, rng: ...` | a function returning a number or a full rate dict |

A value may be a number (a mass-action constant) or a full rate dict such as
`{"law": "hill", ...}`. Draws use the `seed` of the call, so a run is
reproducible. `fill_only=True` rates only the reactions without a rate, the
explicit decision a partially rated network needs.

```python
net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
traj = simulate.ode(net, t_end=10, x0=1.0,
                    rates={"dist": "lognormal", "mean": 0, "sigma": 1}, seed=0)
```

`simulate.assign(net, rates=..., x0=...)` does the same and returns the rated
network, without simulating it.

A generator is a source of networks: vary its arguments, then simulate each
network it builds. The chemistry's own kinetic arguments, such as the
Brusselator's rate constants, are the first place to set rates.

## Stochastic simulation

`ssa` works on molecule counts. The amounts in `initial_state` and in the
frames are per unit volume, so a count is `amount × volume × avogadro`. With
the defaults (`volume=1`, `avogadro=1`) amounts *are* counts; a larger volume
makes the path smoother and closer to the rate equations.

Mass-action rate constants become stochastic ones with the conversion that is
easy to get wrong, implemented once in `chemart.kinetics.k_to_c`: for `m`
reactant molecules, `c = k / (N_A V)^(m-1) × Π l_i!`, where `l_i` is the
multiplicity of each reactant species (a homodimerisation `2 A -> B` gains a
factor 2). The propensity is then `c × Π C(x_i, l_i)`. Other rate laws are
scaled as `Ω f(x/Ω)`, with `Ω = volume × avogadro`.

Each frame of a stochastic path also records the reactions that fired since the
previous frame (`frame.fired`), so fluxes are available per time window.
`max_events` (a million by default) stops a run that would otherwise take too
long; `traj.settings["stopped"]` then says so.

## Flows and buffered species

Three mechanisms appear in the record, and both simulators apply all three:

- **`inflow`**: a constant influx per species. In the SSA it is a source that
  adds one molecule per event.
- **`outflow`**: a first-order removal rate, as a per-species dict, a single
  number for all species, or the string `"constant-total"`.
- **`extras["buffered"]`**: species held at constant amount (their
  derivative is zero; in the SSA their count never changes).

`"constant-total"` is the book's standard flow reactor: a non-selective
dilution flux Φ that holds the total concentration constant,

```
dx_i/dt = production_i - x_i * Phi / sum_k(x_k)
```

In the SSA, a molecule chosen at random is removed whenever an event raises the
total above its starting value.

## Pitfalls

!!! warning "Buffered species look broken but aren't"
    In the Brusselator, `A` and `B` are buffered: their trajectories are flat
    lines because they are reservoirs. Check `net.extras.get("buffered")`
    before reporting that a concentration never changes.

!!! warning "`observed` networks are samples, not systems"
    A gas's network is the record of what fired in one run. Simulating it as a
    rate system will run; the result is rarely meaningful. Check `net.status`
    first.

**Stiffness.** Several of these systems are stiff. The default `LSODA` copes;
if a run fails or crawls, pass `solver="Radau"` or `"BDF"`. For networks that
use only mass action, without the constant-total dilution, the simulator
supplies the exact Jacobian to the solver.

**Units are the model's own.** Chemart imposes none. Two chemistries' rate
constants are not comparable just because both are called `k`.

## Exporting to another tool

```python
d = net.to_dict()        # plain JSON: species, reactions, rates, flows, extras
d = traj.to_dict()       # the trajectory, with the network and the settings
```

From there it is a short hop to SBML, Antimony, Tellurium or COPASI. The record
is deliberately plain so this conversion is boring: reactant and product
multiplicities are already separate, which is the part that usually goes wrong
in translation.
