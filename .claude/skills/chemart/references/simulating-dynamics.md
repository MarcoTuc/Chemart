# Simulating and evolving dynamics

`chemart.simulate` runs a given network, or one a generator built: `ode`
integrates the rate equations, `ssa` samples one stochastic path (Gillespie's
direct method). Turing gases are not simulated this way: their network only
exists as the record of a run, so they are *evolved* with `chemart.evolve`.
All three return a `chemart.trajectory.Trajectory`.

## Contents
- [Quick start](#quick-start)
- [Which networks can be simulated](#which-networks-can-be-simulated)
- [Giving a network rates and an initial state](#giving-a-network-rates-and-an-initial-state)
- [The trajectory](#the-trajectory)
- [What the simulators handle](#what-the-simulators-handle)
- [Evolving a gas](#evolving-a-gas)
- [Pitfalls](#pitfalls)

## Quick start

```python
import chemart
from chemart import simulate

net = chemart.generate_network("brusselator", seed=1)
traj = simulate.ode(net, t_end=40, points=400)           # rate equations
path = simulate.ssa(net, t_end=40, volume=100, seed=0)   # one stochastic path
ids, t, X = traj.array(["X", "Y"])                        # dense arrays for plotting
```

```bash
uv run chemart simulate brusselator --t-end 40
uv run chemart simulate repressilator --seed 1 --t-end 300 --points 600
uv run chemart simulate brusselator --method ssa --volume 100 --seed 0 --format csv > bruss.csv
uv run chemart simulate brusselator --t-end 40 --x0 '{"X": 1.0, "Y": 1.0}'
```

`examples/dynamics.py` is a runnable tour.

## Which networks can be simulated

A network needs a rate for every reaction and an initial state. Many entries
carry both; the rest are topologies by the nature of the model (Kauffman sets,
the RAF model, the Chemical Casting Model). `simulate` raises `NotSimulable`
naming what is missing instead of producing a flat line.

**These integrate cleanly** at default parameters (fully rated, fixed species
set, no structure to lose): `analog-function-crn`, `bagley-farmer`,
`bigan-conservative-crn`, `brusselator`, `chameleon`, `farmer-immune`,
`hill-kinetics`, `jain-krishna`, `kappa-calculus`, `matrix-chemistry`,
`mcs-bl`, `mechanical-self-assembly`, `metabolic-robot-controller`,
`michaelis-menten`, `okamoto-switch`, `random-catalytic-networks`,
`repressilator`, `synthon`.

Some more are rated but carry structure a well-mixed simulator throws away:
`chemoton`, `disperser`, `gard` (compartments) and `flow-ac`, `oregonator`
(space). They will run; the result is a mean-field approximation you chose,
not the published model. Say so if you report it.

`michaelis-menten` and `hill-kinetics` emit **elementary mass action**
(`E + S ⇌ ES → E + P`), not the abridged laws.

## Giving a network rates and an initial state

`rates=` and `x0=` (on `ode`, `ssa` and `simulate.assign`) accept:

- a number: one mass-action constant for every reaction, or one amount for
  every species;
- a distribution: `{"dist": "lognormal", "mean": 0, "sigma": 1}`, also
  `uniform`, `loguniform`, `exponential`, `gamma`, `constant`; one draw per
  reaction or species, from the call's `seed`;
- a table by reaction text or index, or by species id: `{"A + B -> C": 3.0, "*": 0.5}`;
- a list by reaction index;
- a `.json` or two-column `.csv` file;
- a function `(reaction or species id, rng) -> value`.

A rate value may be a number (mass action) or a full rate dict. `fill_only=True`
rates only the reactions without a rate: the explicit decision a partially
rated network needs, never a silent default.

```python
net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
traj = simulate.ode(net, t_end=10, x0=1.0, seed=0,
                    rates={"dist": "lognormal", "mean": 0, "sigma": 1})
```

## The trajectory

`traj.frames` is a list of `Frame(t, state, fired, observables)`: `state` maps
species to amounts (zeros left out), `fired` lists the reactions that fired
since the previous frame (SSA only), `observables` holds quantities a chemistry
reports itself. `traj.array(species)` gives `(ids, t, X)`, `traj.window(i, w)`
the network of reactions fired in a range of frames, `traj.to_dict()` plain
JSON. `traj.settings` records the seed, solver and specs used.

## What the simulators handle

- the rate laws `mass-action`, `power`, `michaelis-menten`, `hill`,
  `saturating` (see `network-record.md` for their parameter names), and
  `arrhenius` when `T` and `R` are on the rate dict or `temperature=` and
  `gas_constant=` are passed;
- `inflow` as a constant influx (a source in the SSA);
- `outflow` as first-order removal, per-species dict or single number;
- `outflow == "constant-total"`, the book's non-selective dilution flux (in
  the SSA a random molecule is removed whenever an event raises the total);
- `extras["buffered"]`, species pinned at their initial amount.

Because flows and buffering are applied for you, do not add dilution or
clamping terms yourself on top; you will double-count them.

For the SSA, amounts are per unit volume and counts are `amount × volume ×
avogadro` (both 1 by default, so amounts are counts). Mass-action constants
are converted with `chemart.kinetics.k_to_c`, which includes the combinatorial
factor for identical reactants (`2A → …` gains a factor 2). `max_events`
(10^6 by default) stops long runs; `traj.settings["stopped"]` says so.

## Evolving a gas

A chemistry with an evolve face (every gas, plus a few lattices and one SSA
generator; `describe_chemistry(id)["faces"]`) runs its own process:

```python
run = chemart.evolve("alchemy", seed=1, collisions=4000)   # evolve-face params allowed
run.clock                      # "collisions": t is counted in the chemistry's own unit
run.frames[-1].state           # the soup at the end, {species: count}
run.series("prime_fraction")   # a chemistry's own observable (here: prime-number-chemistry)
run.network                    # what fired over the whole run, status "observed"

for frame in chemart.evolve_frames("bff", seed=1, every=4):   # live; merges 4 frames into one
    ...
```

Every frame carries `fired`, the reactions since the previous one, so
`traj.window(i, w)` gives the network of a stretch of the run, and
`chemart.measures.over(run, names, window=w)` follows measures frame by frame:
population measures (`richness`, `shannon`, `dominance`) read each frame's
state; network measures read the reactions fired in the last `w` frames
(`window=None`: everything so far).

```bash
uv run chemart evolve alchemy --seed 1 --track shannon --track n_species --window 5
uv run chemart evolve prime-number-chemistry --seed 1 --format csv --species n2 n3 n5
```

Clocks differ between gases (collisions, epochs, generations), so to compare
two gases in time, use `traj.turnover()` (cumulative reactions fired per
molecule) as the shared clock.

## Pitfalls

**Buffered species look broken but aren't.** In the Brusselator, `A` and `B`
are buffered, so their trajectories are flat lines. Check
`net.extras.get("buffered")` before reporting "the concentration never
changes".

**`observed` networks are samples, not systems.** A gas's network is one run's
firing record. Simulating it as a rate system will run, and the result is
usually meaningless. Check `net.status` first.

**Defaults are small.** A default network may be too small to show the
published behaviour. `describe_chemistry(id)["params"]` gives each parameter's
description including a `range:` note with the paper's scale.

**Stiffness.** The default solver `LSODA` copes with most systems; if a run
fails or crawls, pass `solver="Radau"` or `"BDF"`.

**Units are the model's own.** Two chemistries' rate constants are not
comparable just because both are called `k`.
