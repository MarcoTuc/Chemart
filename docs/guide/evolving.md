# Evolving a chemistry

A Turing gas has no network to begin with. Its molecules are programs, strings
or terms, and a procedure decides what two of them make when they meet, so the
reactions exist only as they happen. What you study is the process itself: the
population changes, new species appear, and some organisations take over.
Chemart runs that process with `chemart.evolve`. It records the run as a
[`Trajectory`](../reference/record.md#the-trajectory-record), a list of frames
in the chemistry's own time.

The chemistries you can evolve are the ones with an **evolve face**:

- **The gases.** Some have only this face: bff, combinatory-chemistry, ccm,
  music-ac, rbn-world, sr-loops, srsim, dorin-korb-ecosystem and
  molecular-tsp. The rest also have a generate face that returns the closure
  of their rule, among them alchemy, matrix-chemistry, nac, rna-folding-ac,
  squirm3, stringmol and typogenetics.
- **Chemistries with a reactor of their own.** For example, the lattices of
  autopoiesis-vmu and ono-ikegami-protocell, and the stochastic run of synthon.

`describe_chemistry(id)["faces"]` says which faces a chemistry has, its
`clock` the unit its time is counted in and its `duration` the parameter that
says how long it runs. The [catalog](../catalog/index.md) groups entries by
type.

## Running one

```python
import chemart

traj = chemart.evolve("alchemy", seed=1)            # the whole run
traj.clock                                          # 'collisions'
traj.frames[-1].state                               # {species: amount} at the end
traj.network                                        # every reaction that fired, with counts

traj = chemart.evolve("bff", seed=1, epochs=200, every=10)    # keep one frame in ten
```

Parameters work as in `generate_network`. A parameter that belongs to one face
only is marked in the catalog (*evolve only*), and passing it to the other face
is an error rather than a silent no-op. `every` merges frames: the reactions
fired in between add up, and the last state stands.

For a live run, `chemart.evolve_frames` yields the frames as they come:

```python
run = chemart.evolve_frames("combinatory-chemistry", seed=2)
for frame in run:
    print(frame.t, len(frame.state), frame.observables.get("reductions"))
```

From the command line:

```bash
uv run chemart evolve alchemy --seed 1 --track shannon --track n_species
uv run chemart evolve bff --seed 1 -p epochs=200 --every 10 --format csv > bff.csv
uv run chemart evolve bff --seed 1 --format json > run.json    # the whole trajectory
```

To watch a run instead, use the [simulation pit](../hub.md#the-simulation-pit)
(`uv run chemart-hub pit`). It plots the population, the chemistry's
observables and the measures you pick while the run is still going.

## Running until you stop

Every process has a parameter that says how long it runs — alchemy's
`collisions`, bff's `epochs`, the lattices' `steps`. The catalog records which
one it is (`duration`), and **0 means it runs until you stop reading its
frames**:

```python
run = chemart.evolve_frames("alchemy", seed=1, collisions=0)
for frame in run:
    print(frame.t, len(frame.state))
    if bored(frame):
        break                       # closing the generator ends the process
```

In the [pit](../hub.md#the-simulation-pit) that is the **run until I stop**
box, with the Stop button ending it; what ran is kept, so it can still be
downloaded. `chemart.evolve` refuses an endless run, because the whole
trajectory would never arrive — iterate `evolve_frames` instead.

A process may still end on its own before you stop it: a fraglets program goes
inert, proof-ac finds its proof, a population dies out. That is the chemistry
finishing, not the budget running out.

## What a frame holds

| field | contents |
|---|---|
| `t` | the time of the frame, in the chemistry's clock |
| `state` | `{species: amount}`: the population at that time |
| `fired` | the reactions fired since the previous frame, as `[[reactants], [products], count]` |
| `observables` | quantities the chemistry itself reports (bff's `high_order_entropy`, combinatory-chemistry's `reductions`, the membranes of autopoiesis-vmu) |

The first frame is the initial population and has fired nothing. Frame by
frame, the `fired` counts add up to the counts in `traj.network`, and the same
seed gives the same frames. The contract checks both for every chemistry.

`traj.series(name)` is one observable over time, and `traj.array(species)` is
the amounts as dense arrays for plotting. `traj.window(i, width)` is the
network of the reactions fired in frames `i-width+1` to `i`.

## Time

Each chemistry keeps its native clock. For example, AlChemy counts
`collisions`, bff counts `epochs`, combinatory-chemistry counts `iterations`
and the lattices count `steps` or `sweeps`. The catalog records the clock of
each entry. Clocks do not compare across chemistries: a thousand collisions
of λ-terms and a thousand bff epochs are different amounts of chemistry.

`traj.turnover()` is a clock every gas shares: the cumulative number of
reactions fired per molecule of the population. Plotted against turnover
rather than native time, two gases are compared by how many times their
population has been rewritten.

## Measures in evolutionary time

A gas's measures change as it runs, so you follow them:

```python
from chemart import measures

series = measures.over(traj, ["shannon", "dominance", "n_species", "cycle_rank"], window=5)
series["t"], series["shannon"]          # one value per frame (None where it did not apply)

measures.measure(traj)                  # trajectory measures plus those of the cumulative network
```

Two kinds of measure are involved:

- **Population measures** (richness, Shannon diversity, dominance, population)
  read each frame's `state`. They describe what the gas *is* at that time,
  and they are the first axes to compare gases on.
- **Network measures** read the network of the reactions fired in the last
  `window` frames (`window=None` for everything so far). They describe what
  the gas *does*: how it is organised, and whether its reactions are
  autocatalytic. A measure that needs a food set applies to a window only if
  you pass `food=`.

`measures.track(chemart.evolve_frames(...), names)` does the same on a live
run, one row per frame. Trajectory measures such as `turnover`,
`novelty_rate`, `complexity_drift`, `collapse_time` and `attractor_type`
summarise a whole run. See the [measures page](measures.md) for all of them.

!!! warning "An observed network is a sample"
    The network of a window is what happened to fire in it, so it mixes the
    chemistry with how much of it the run happened to see. A short window
    on a rich gas looks sparse because it is small, not because the
    chemistry is. Compare against a null model built in the same window
    (`measures.zscores(traj.window(i, w))`), keep the window size fixed when
    comparing gases, and read population measures first.

## Two faces of one chemistry

Many gases also have a generate face. For them, `generate_network` returns
something different from a run: the closure of the rule from a seed set, cut
off by `max_species`. That is everything the rule *can* make, while a run shows
what a finite population *does* make. The two answer different questions, and
their networks differ in status (`complete` or `truncated` against
`observed`). For a gas with only an evolve face, `generate_network` runs the
process to the end and returns the observed network.

## Writing an evolve face

A chemistry module can define, besides or instead of `generate(p, rng)`:

```python
def evolve(p, rng):
    population = ...
    tally = Tally()                      # chemart.soup.Tally counts what fired
    yield Frame(t=0, state=counts(population))
    for step in range(1, p.steps + 1):
        ...                              # react, and tally.add(lhs, rhs, 1) per reaction
        if step % p.every_n == 0:
            yield Frame(t=step, state=counts(population), fired=tally.flush())
    return network_of(tally.reactions())  # the observed Network
```

The catalog entry then needs a `clock`. For a well-stirred soup,
`chemart.soup.stir` already does the loop: it yields the population and the
tally at every frame. [Contributing a chemistry](../contributing.md) has the
rest of the contract.
