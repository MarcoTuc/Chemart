## Introduction

Most chemistries in this catalog are written down: someone lists the molecules
and the rules by which they react. This one is *read off* instead. You start
from a cellular automaton, a row of cells that are each 0 or 1 and that all
update at once, each cell looking only at its close neighbours. You run it, and
you draw the successive rows under one another to make a *space-time diagram*.
In some automata that picture falls apart into large regular patches, all 0s,
all 1s, or a checkerboard `0101…`, and between the patches there are thin
boundaries that move at constant speed, meet, and disappear or turn into other
boundaries. Treat the boundaries as molecules and their meetings as reactions,
and you have a reaction network such as `gamma + delta -> ∅` or
`beta + gamma -> eta`. Nobody programmed it: it is what the automaton does,
described at a coarser level.

The description comes from a group at the Santa Fe Institute and Berkeley. James
Crutchfield and James Hanson developed *computational mechanics* for cellular
automata in the early 1990s: a method for finding an automaton's regular patches
(*domains*), filtering them out of the space-time diagram, and studying what is
left over (*particles*). Melanie Mitchell, Crutchfield and colleagues (Peter
Hraber, then Rajarshi Das) used a genetic algorithm to evolve automata that
perform a global task, *density classification*: decide whether the starting
row holds more 1s or more 0s, and say so by turning every cell into 1 or into
0. That work began as a re-examination of earlier claims about "computation at
the edge of chaos" (Mitchell, Hraber and Crutchfield 1993). The best evolved
automata did it with particles, and the question the two lines of work answer
together is how a system in which no cell sees more than seven cells can
compute a property of the whole row. Crutchfield, Mitchell and Das (1998) answer it by
listing each automaton's domains, particles and particle interactions in a
*particle catalog*; Hordijk, Crutchfield and Mitchell (1996, 1998) then showed
that the catalog alone, run as a simple model of colliding particles, predicts
how well the automaton performs.

Banzhaf and Yamamoto describe the idea in a short passage, "Embedded Particles
in Cellular Automata as Molecules", inside book §10.7.2 on cellular automata
seen as chemistries, with one filtered space-time diagram (figure 10.14, from
Hordijk, Crutchfield and Mitchell 1998). The contrast they draw is with the
[self-replicating loops](sr-loops.md) of the same section, where a molecule is
a connected group of cells in particular states: here a molecule is a
*boundary* between two domains, and it becomes visible only as a structure in
space-time. Of the loops, the book remarks that once molecules are emergent
patterns, their definition "could be dependent on an observer or on an
interpretation"; here too a particle exists only relative to the domains the
observer chose to filter out. The neighbouring entries of §10.7 use automata
differently again: in [Avida](avida.md) each cell holds a self-replicating
program, and in [bondable cellular automata](bondable-ca.md) each atom is a
whole small automaton.

So this entry is an **analysis method**, not a simulation model that someone
proposed as a chemistry. Chemart offers it both ways: it hands back the
published particle catalog of an evolved automaton as a complete reaction
network, and it runs the automaton and reports the collisions a domain filter
actually sees, with their counts.

## How it works

### The automaton and its task

The automata here are one-dimensional rows of `N` cells joined into a ring (the
last cell is next to the first). At every time step every cell looks at its
*neighbourhood*: itself and the three cells on each side, seven cells in all
(the *radius* is 3). Seven binary cells can be in 2⁷ = 128 patterns, and the
automaton's *rule* is a look-up table giving the new state of the centre cell
for each of them. The papers print these 128-bit tables as 32 hexadecimal
digits; Chemart embeds them verbatim.

The task is *density classification*. Write ρ₀ (rho-zero) for the fraction of
1s in the initial configuration (IC). If ρ₀ > 1/2, the automaton should turn
the whole ring into 1s; otherwise into 0s; and it must do so within
`T_max = 2N` steps. The answer counts only if the ring has reached all-1s or
all-0s; there is no partial credit. A rule's *performance* `P_N` is the
fraction of random ICs on `N` cells it classifies correctly. Counting is
trivial for an ordinary computer; the difficulty is that no cell can count, and
a cell's local majority is often wrong about the global one.

Three rules are included. `phi-par-a` (φ_par^a in the papers) and `phi-par-b`
are two of the particle-based automata evolved by the genetic algorithm of
Crutchfield, Mitchell and Das (1998). `gkl` is the rule of Gács, Kurdyumov and
Levin, designed by hand for a different purpose (reliable computation) but
behaving much like the evolved ones.

### Domains: the solvent

Watch φ_par^a from a random IC and within a few steps the ring is covered by
three kinds of patch: all 0s, all 1s, and the checkerboard `0101…`. These are
its *regular domains*, written `Λ0 = 0+`, `Λ1 = 1+` and `Λ2 = (01)+` (the `+`
means one or more repetitions). Formally a regular domain is a pattern that the
rule maps onto itself, so a patch of it stays a patch of it, and that can occur
anywhere on the ring. Crutchfield, Mitchell and Das (1998) point out that a
domain does very little computing: it just keeps producing its own pattern. In
the chemical reading the domains are the solvent, not molecules. (φ_par^b has
the same `Λ0` and `Λ1`, but a striped `(011)+` in place of the checkerboard.)

### Particles: the molecules

A *domain filter* marks every cell that belongs to one of the domains and leaves
the rest. What is left are the *walls* between domains. When a wall stays
narrow and moves at a steady speed it is a *particle*. A particle is named by
the two domains it separates, left then right, written `p ~ ΛiΛj`, and its
velocity is in cells per step (negative is leftwards). The catalog of φ_par^a
(Crutchfield, Mitchell and Das 1998, Table 3) has six:

| particle | wall | velocity | what it is |
|---|---|---|---|
| `alpha` | `Λ0Λ1` | unstable | 0s on the left, 1s on the right; decays at once |
| `beta` | `Λ1Λ0` | 0 | 1s on the left, 0s on the right; stands still |
| `gamma` | `Λ0Λ2` | −1 | 0s, then checkerboard |
| `delta` | `Λ2Λ0` | −3 | checkerboard, then 0s |
| `eta` | `Λ1Λ2` | +3 | 1s, then checkerboard |
| `mu` | `Λ2Λ1` | +1 | checkerboard, then 1s |

A wall between two patches of the *same* domain (for example two checkerboards
out of step with each other) is a *dislocation*. The published catalogs leave
these out; Chemart keeps them as species named `wL2L2`, `wL0L0` and so on.

### Reactions: collisions

When two particles meet, the result is fixed by the rule. φ_par^a's catalog
lists six interactions:

```
alpha -> gamma + mu        decay
beta + gamma -> eta        reaction
mu + beta -> delta         reaction
eta + delta -> beta        reaction
eta + mu -> ∅              annihilation, leaving Λ1
gamma + delta -> ∅         annihilation, leaving Λ0
```

The walls must add up. `beta` is `Λ1Λ0` and `gamma` is `Λ0Λ2`; when they meet,
the `Λ0` patch between them is squeezed out and what remains is `Λ1Λ2`, which is
`eta`. `gamma + delta` is `Λ0Λ2` followed by `Λ2Λ0`: the checkerboard between
them vanishes, and so do both walls, leaving plain `Λ0`. The catalog says it
leaves out possible three-particle interactions.

This is how the density computation works (Crutchfield, Mitchell and Das 1998,
section 6.3). A `beta` or an `alpha` sits where a region of 1s meets a region
of 0s, a place where the density is locally ambiguous. The `alpha` decays into
a `gamma` moving left and a `mu` moving right, at the same speed, so the growing
checkerboard between them is centred on where the `alpha` was. Whichever of the
two reaches a `beta` first has crossed the smaller region: if `mu` wins, the
black (1s) region was the smaller, and `mu + beta` makes a `delta` that races
left at three times `gamma`'s speed, catches the `gamma` and annihilates with
it, leaving 0s. If `gamma` wins, `beta + gamma` makes an `eta` that catches the
`mu`, leaving 1s. Repeated over larger and larger regions, these contests end
with one domain covering the ring.

### A worked example

Here is φ_par^a on 75 cells (the size of the book's figure 10.14) from an IC
with 36 ones, ρ₀ = 0.48. On the left is the raw configuration (`1` for a 1,
`.` for a 0); on the right the same row after domain filtering: `0`, `1` and
`c` mark cells in `Λ0`, `Λ1` and the checkerboard `Λ2`, and `#` marks cells in
a wall. The code is under *Using it* below.

```
 0 111..1.1.1.111.1..1..111.1.....1.1.1....1..1.1..1...1.1111.11.1111..1...1.1   #######cc########################c#########################################
 1 1.11.11.1.1.11.1..1.1111.1....1.1.........11.........1.11111.111.111.....1.   #########c##########################000########000#########################
 2 ....11.111....11.111.111.1...1.............1........1.1.1111111111..1...1.1   #################################0000000#######00##########1111############
 3 ....1.11..11..1.11111111.11........................1.1.1.11111111111.......   0##################11#########000000000000000000#####cc#####11111######0000
 4 ...1.....1....11.1111111111.......................1.1.1.1.111111111.1......   ####################1111######00000000000000000#####cccc#####111########000
 5 ..............1.1111111111.1.....................1.1.1.1.1.11111111.1......   00000000000########1111########000000000000000#####cccccc#####11########000
 6 .............1.1.111111111.1....................1.1.1.1.1.1.1111111.1......   0000000000##########111########00000000000000#####cccccccc#####1########000
 7 ............1.1.1.11111111.1...................1.1.1.1.1.1.1.111111.1......   000000000#####cc#####11########0000000000000#####cccccccccc#############000
 8 ...........1.1.1.1.1111111.1..................1.1.1.1.1.1.1.1.11111.1......   00000000#####cccc#####1########000000000000#####cccccccccccc############000
 9 ..........1.1.1.1.1.111111.1.................1.1.1.1.1.1.1.1.1.1111.1......   0000000#####cccccc#############00000000000#####cccccccccccccc###########000
10 .........1.1.1.1.1.1.11111.1................1.1.1.1.1.1.1.1.1.1.111.1......   000000#####cccccccc############0000000000#####cccccccccccccccc##########000
11 ........1.1.1.1.1.1.1.1111.1...............1.1.1.1.1.1.1.1.1.1.1.11.1......   00000#####cccccccccc###########000000000#####cccccccccccccccccc#########000
12 .......1.1.1.1.1.1.1.1.111.1..............1.1.1.1.1.1.1.1.1.1.1....11......   0000#####cccccccccccc##########00000000#####ccccccccccccccccc###########000
13 ......1.1.1.1.1.1.1.1.1.11.1.............1.1.1.1.1.1.1.1.1.1........1......   000#####cccccccccccccc#########0000000#####ccccccccccccccc#####00#######000
14 .....1.1.1.1.1.1.1.1.1....11............1.1.1.1.1.1.1.1.1..................   00#####ccccccccccccc###########000000#####ccccccccccccc#####000000000000000
15 ....1.1.1.1.1.1.1.1........1...........1.1.1.1.1.1.1.1.....................   0#####ccccccccccc#####00#######00000#####ccccccccccc#####000000000000000000
16 ...1.1.1.1.1.1.1......................1.1.1.1.1.1.1........................   #####ccccccccc#####0000000000000000#####ccccccccc#####000000000000000000000
17 ..1.1.1.1.1.1........................1.1.1.1.1.1...........................   ####ccccccc#####000000000000000000#####ccccccc#####00000000000000000000000#
18 .1.1.1.1.1..........................1.1.1.1.1..............................   ##cccccc#####00000000000000000000#####ccccc#####0000000000000000000000000##
19 1.1.1.1............................1.1.1.1.................................   ###cc#####0000000000000000000000#####ccc#####000000000000000000000000000###
20 .1.1..............................1.1.1...................................1   #######000000000000000000000000#####c#####00000000000000000000000000000####
21 1................................1.1.....................................1.   ####00000000000000000000000000#########0000000000000000000000000000000#####
22 ................................1.......................................1..   #0000000000000000000000000000#######000000000000000000000000000000000######
23 ...........................................................................   000000000000000000000000000000000000000000000000000000000000000000000000000
```

Step by step:

- **Steps 0 to 4.** The random IC is not yet made of domains; most cells are
  wall. The first step at which every cell is either in a domain or in a narrow
  wall between two domains is the *condensation time* `t_c`. Here `t_c = 5`.
- **Step 5.** Reading left to right, the ring holds an `alpha` (0s then 1s,
  around cell 14), a `beta` (1s then 0s, around cell 26), a `gamma` (0s then
  checkerboard, around cell 48), a `mu` (checkerboard then 1s, around cell 59)
  and another `beta` (around cell 67). The right-hand checkerboard is already
  growing: its `gamma` edge moves left one cell per step and its `mu` edge right
  one cell per step. It is the product of an earlier `alpha` decay, at step 2.
- **Step 7.** The left `alpha` has decayed: a second checkerboard opens at cell
  14, with a `gamma` on its left and a `mu` on its right (`alpha -> gamma + mu`).
  On the right, the `mu` has reached the `beta` at cell 67 and the black patch
  between them is gone: `mu + beta -> delta`.
- **Steps 8 to 14.** The left `mu` reaches its `beta` (cell 26) and makes a
  second `delta`, which takes a few steps to form. Both `delta`s then move left
  three cells per step, eating the checkerboards from the right, while the
  `gamma`s ahead of them retreat at one cell per step.
- **Steps 19 to 23.** Each `delta` catches its `gamma` (the left pair meets
  across the wrap-around at cell 0), and both pairs annihilate:
  `gamma + delta -> ∅`, leaving `Λ0`. At step 23 the ring is all 0s, the right
  answer for ρ₀ = 0.48.

The 1s lost because the black patches were small: each checkerboard stood for a
contest between a patch of 0s and a patch of 1s, and in both contests the `mu`
reached its `beta` first.

### How Chemart reads reactions from a run

`chemart.evolve` runs the published look-up table and, at every step, does
what the example above did by eye.

1. **Filter.** A cell is assigned to a domain when the seven cells around it
   (the window is set by `filter_window`, default 3 on each side) match one
   phase of that domain's pattern. The papers use a finite-state *transducer*
   for this (a machine that reads the row and outputs domain or wall for each
   cell); Chemart's window test is simpler, as the decisions below explain.
2. **Name the walls.** Each maximal run of unassigned cells is a wall, named by
   the domains on either side: a catalogued particle if the pair is in the
   catalog, `wLiLj` otherwise.
3. **Follow them.** A particle is matched to the nearest particle of the same
   type within 4 cells at the next step. This uses no published velocity, so
   the velocities it yields are a measurement.
4. **Collect events.** A collision takes several steps in the automaton, during
   which the row cannot be split cleanly into domains and particles. Whatever
   disappears and appears within 14 cells and 4 steps of an event is collected
   into that one event, and a particle that disappears and comes back within
   that window is treated as having survived. Each finished event is a
   reaction: reactants are what vanished, products what appeared.
5. **Discard the transient.** Events that begin before `t_c` are dropped as
   filter debris.

The last two rules have a side effect, visible in the example. An event that
opens before `t_c` stays open for as long as something keeps changing within 14
cells of it, and when it is finally dropped it takes later collisions with it.
In the run above Chemart reports only `gamma + delta -> ∅` (twice) and
`beta + mu -> delta` (once): the left `alpha -> gamma + mu` at step 7 and the
second `mu + beta -> delta` at step 9 were absorbed into an event that started
at step 1 and were discarded with it. Observed counts are therefore a lower
bound on what happened after condensation.

In the formal specification below, the reactor is labelled `lattice-2d`, the
closest type the catalog has: the lattice itself is the one-dimensional ring
described here, and it is its space-time diagram that is two-dimensional.

## Using it

The chemistry has two faces. The call printed above,
`chemart.generate_network("ca-embedded-particles")`, runs nothing: it returns
φ_par^a's published interaction table as a complete network
(`status="complete"`), and its only parameter is `rule`.
`chemart.evolve("ca-embedded-particles")` runs the automaton and returns a
trajectory with a frame per CA iteration; `lattice`, `steps`, `density` and
`filter_window` belong to it.

#### The published catalog as a network

```python
net = chemart.generate_network("ca-embedded-particles")
for s in net.species:
    print(s.id, s.structure)
for r in net.reactions:
    print(r.to_text())
```

```
alpha Λ0Λ1 (v=unstable)
beta Λ1Λ0 (v=0.0)
gamma Λ0Λ2 (v=-1.0)
delta Λ2Λ0 (v=-3.0)
eta Λ1Λ2 (v=3.0)
mu Λ2Λ1 (v=1.0)
alpha -> gamma + mu
beta + gamma -> eta
mu + beta -> delta
eta + delta -> beta
eta + mu -> ∅
gamma + delta -> ∅
```

With `rule="phi-par-b"` the same call returns φ_par^b's catalog (Crutchfield,
Mitchell and Das 1998, Table 4): the same six reactions over differently placed
walls, with `alpha` at velocity 0, `beta` +1, `gamma` 0, `delta` −3, `eta` +3
and `mu` +3/2. That is the only way to use φ_par^b: its printed look-up table
does not classify density when run (see the decisions), so
`chemart.evolve(..., rule="phi-par-b")` raises an error. The seed plays no
part here.

#### A run of the automaton

The default run is φ_par^a on 149 cells for 298 = 2N steps (the task's answer
time) from an IC with 72 ones (ρ₀ = 0.483); the published space-time figures
use ρ₀ = 0.48 and 0.51. The automaton classifies this IC correctly.

```python
traj = chemart.evolve("ca-embedded-particles", seed=1)
net = traj.network
for r in net.reactions:
    print(r.to_text())
```

```
∅ -> wL2L2  (x1)
beta + mu -> delta  (x1)
delta -> ∅  (x1)
delta + gamma -> ∅  (x1)
eta -> beta  (x1)
mu -> delta  (x1)
wL2L2 -> ∅  (x1)
```

These are the events the filter saw after `t_c = 5`. On a ring this size only
a handful of collisions happen, and several of the events are one-sided
(`delta -> ∅`, `mu -> delta`): a particle vanishing or appearing with its
partner assigned to a different event, or a dislocation (`wL2L2`) forming and
healing.

The frames follow the run step by step:

```python
t_c = net.extras["analysis"]["condensation_time"]          # 5
traj.frames[t_c].state      # {'beta': 1.0, 'delta': 1.0, 'gamma': 2.0, 'mu': 1.0, 'wL0L0': 1.0}
[(f.t, f.fired) for f in traj.frames if f.fired][:2]
# [(11.0, [[['mu'], ['delta'], 1]]), (18.0, [[['delta'], [], 1], [['eta'], ['beta'], 1]])]
[round(d, 3) for d in traj.series("density")[::50]]
# [0.483, 0.356, 0.02, 0.0, 0.0, 0.0]
```

A frame's `state` counts the particles present at that iteration, by name,
and its `fired` lists the events recorded in that iteration. An event is
recorded once it has settled, up to 4 iterations after the collision (events
still open when the run ends are in the last frame), and never before `t_c`.
The observable `density` is the fraction of cells in state 1; here the ring
reaches all 0s by iteration 102 and stays there.

What to read in the network:

- `net.species`: the six catalogued particles plus any other wall seen; each
  `structure` gives the wall and its published velocity, e.g.
  `Λ0Λ2 (v=-1.0)`.
- `net.reactions`: the observed events, each with its count, most frequent
  first. `net.status` is `observed`, since this is what one run showed, not a
  complete network.
- `net.initial_state`: the particles present at `t_c`, the same as
  `traj.frames[t_c].state`; the first frame is the particles at iteration 0.
- `net.extras["analysis"]`: `condensation_time`; `classification` (`correct`,
  `incorrect` or `undecided`); `initial_density`; `particles_seen` (how many
  times each particle type was seen, summed over steps); `measured_velocities`
  (mean displacement per step of every particle followed for at least 5
  steps); `interactions_observed`; and the rule's `catalog`, `lookup_hex` and
  `published_performance`.
- `net.extras["space"]`: the `initial` and `final` rows as strings of 0 and 1,
  and `filtered`, the filtered space-time diagram as rows of `.` (domain) and
  `#` (wall), thinned to about 200 rows on long runs.

The velocities measured in a whole run are rough, averaged over the one to
three tracks per particle type that last long enough: in the default run
`gamma` and `mu` come out at exactly −1 and +1, but `delta` at −1.9 (three
tracks) and `eta` at +0.33 (one track). Use `probe`, below, for clean values.

#### The worked example

The run of *How it works*, with the replay that prints each configuration and
its filtered version (takes about a second):

```python
import numpy as np
import chemart
from chemart.chemistries import ca_embedded_particles as C

net = chemart.evolve("ca-embedded-particles", seed=7, lattice=75, steps=150).network
for r in net.reactions:
    print(r.to_text())
print(net.extras["analysis"]["condensation_time"], net.extras["analysis"]["classification"])

# replay the same run, printing each configuration next to its filtered version
spec = C.RULES["phi-par-a"]
table = C.lookup(spec["hex"])
x = np.array([int(c) for c in net.extras["space"]["initial"]], dtype=np.int8)
for t in range(24):
    labels = C.domain_labels(x, spec["domains"], 3)
    raw = "".join("1" if v else "." for v in x)
    filtered = "".join("#01c"[v + 1] for v in labels)
    print(f"{t:2d} {raw}   {filtered}")
    x = C.step(x, table)
```

```
delta + gamma -> ∅  (x2)
beta + mu -> delta  (x1)
5 correct
```

followed by the 24 rows shown above.

#### Measuring one particle: `probe`

`C.probe(rule, particle)` builds a ring that is one domain on the left half and
the other on the right, so that it holds exactly the wall of interest, runs it
for 20 steps and follows that wall:

```python
for p in C.ORDER:
    print(C.probe("phi-par-a", p))
```

```
{'particle': 'alpha', 'wall': 'Λ0Λ1', 'velocity': 0.5, 'steps_followed': 4, 'decays_into': ['gamma', 'mu']}
{'particle': 'beta', 'wall': 'Λ1Λ0', 'velocity': 0.0, 'steps_followed': 20, 'decays_into': []}
{'particle': 'gamma', 'wall': 'Λ0Λ2', 'velocity': -0.975, 'steps_followed': 20, 'decays_into': []}
{'particle': 'delta', 'wall': 'Λ2Λ0', 'velocity': -3.025, 'steps_followed': 20, 'decays_into': []}
{'particle': 'eta', 'wall': 'Λ1Λ2', 'velocity': 2.975, 'steps_followed': 20, 'decays_into': []}
{'particle': 'mu', 'wall': 'Λ2Λ1', 'velocity': 1.025, 'steps_followed': 20, 'decays_into': []}
```

The five stable particles move at their published velocities, within 0.025,
and `alpha` survives 4 steps before it is replaced by a `gamma` and a `mu`.

#### Performance on the task

`C.performance(rule, n, trials, rng)` measures `P_N` the way the papers do:
unbiased random ICs (each cell 0 or 1 with probability 1/2), `T_max = 2N`, no
partial credit. It does not build a network.

```python
import numpy as np
for rule in ["phi-par-a", "gkl"]:
    print(rule, C.performance(rule, n=149, trials=1000, rng=np.random.default_rng(1)))
print(C.performance("phi-par-a", n=599, trials=200, rng=np.random.default_rng(1)))
```

```
phi-par-a 0.783
gkl 0.795
0.795
```

The two runs at N = 149 take about 19 seconds together, the 200 ICs at
N = 599 about 8.5 seconds. With 1,000 ICs the sampling error is about ±0.013;
with 200 it is about ±0.03.

#### Collisions over several runs

The six-seed run of the slow test, on 599 cells for 1,198 steps (about 6
seconds):

```python
from collections import Counter
fired = Counter()
for seed in range(6):
    net = chemart.evolve("ca-embedded-particles", seed=seed, lattice=599,
                         steps=1198, density=0.48 if seed % 2 else 0.52).network
    for r in net.reactions:
        fired[r.to_text().split("  (")[0]] += r.count
```

The most frequent events were `delta + gamma -> ∅` (28), `beta + mu -> delta`
(24), `beta + gamma -> eta` (17), `eta + mu -> ∅` (16) and
`delta + eta -> beta` (8), then dislocations appearing and disappearing (8 and
7) and a tail of one-sided and multi-particle events. `alpha -> gamma + mu`
was recorded once: `alpha` decays within a few steps, usually before `t_c`.
Five of the six ICs were classified correctly; `t_c` ranged from 6 to 9.

#### The GKL rule

`chemart.evolve(..., rule="gkl")` runs the Gács-Kurdyumov-Levin rule. It has the same three domains,
but no particle catalog is published for it in these sources, and its walls do
not move like φ_par^a's, so no Greek names are borrowed: every wall is a
`wLiLj` species with a measured velocity, and `generate_network` raises an
error for it. With `seed=1` it reports `wL1L2 + wL2L1 -> ∅` twice,
`wL0L1 + wL2L0 -> wL2L1` once, and a dislocation forming and healing; this IC
it classifies incorrectly.

The other parameters are in the table below: `lattice` and `steps` set the
ring and the run length (the default run takes a fraction of a second, and
runs grow with `lattice × steps`), `density` the IC, and `filter_window` how wide a patch
must be to count as domain (larger values report wider walls and different
particle counts).

## Results

### Evolving automata that classify density

Crutchfield, Mitchell and Das (1998) report 300 runs of their genetic algorithm
on radius-3 rules, with populations of 100 rules scored on 100 ICs per
generation for 100 generations at N = 149. The runs ended with three kinds of
strategy. In 11 runs the best rule was a *default* strategy: go to all 0s (or
all 1s) whatever the IC, which scores about 0.5. In 280 runs it was *block
expanding*: go to one answer unless the IC contains a block of equal cells
about as long as the neighbourhood (7 or 8 cells), and then expand that block.
These score about 0.6 at N = 149 and fall towards 0.5 on longer rings, because long blocks become common at any
density. In 9 runs the genetic algorithm found *embedded-particle* strategies,
and these generalise to longer rings. Their Table 1 gives, over 10⁴ unbiased
ICs:

| rule | P_149 | P_599 | P_999 |
|---|---|---|---|
| φ_par^a | 0.775 | 0.740 | 0.728 |
| φ_par^b | 0.766 | 0.687 | 0.641 |
| block expanding (φ_exp^a) | 0.656 | 0.523 | 0.504 |
| default (φ_def^a) | 0.500 | 0.500 | 0.500 |

The plain radius-3 majority vote scores zero: patches of 1s and 0s form, and
nothing decides between them. Crutchfield and Mitchell (1995) give GKL's
performance as 0.816, 0.766 and 0.757.

### The particle catalog explains the strategy

The catalog turns the space-time diagram into an account of the computation
(section 6 of the 1998 paper). The mechanism of *How it works* is theirs, and
they derive quantities from it. In one recurring pattern, a white region of
length `W` and a black region of length `B` between two `beta`s go through
`beta alpha beta → beta gamma mu beta → eta delta → beta`, and the new `beta`
ends up `2(B − W)` cells to the right of the original `alpha`: the black
region has gained `B − W` cells, which is the information later collisions
use.

The same analysis says why φ_par^b is worse. Its striped domain plays the part
of the checkerboard, but its velocities are asymmetric: `gamma` closes on
`beta` at one cell per step, `mu` at only half a cell per step. On one IC
with ρ₀ = 0.52 that condensed into 85 black cells next to 64 white ones, the
`gamma` reached its `beta` after about 85 steps and the `mu` after about 103,
so the white region won and the IC was misclassified. In φ_par^b's version of
the `beta gamma mu beta` pattern the white region gains `2W − B` cells, so it
grows only if it is at least half the size of the black one. φ_par^a makes two
kinds of mistakes: either the configuration at `t_c` already has the wrong
majority, or small islands get cut off by the geometry, which becomes more
common as N grows and is why `P_N` falls. On a set of 10⁴ ICs φ_par^a
classified 81% of the low-density and 74% of the high-density ones, taking 81
steps on average and 227 at most.

### The catalog predicts performance

Hordijk, Crutchfield and Mitchell turned the catalog into an *embedded-particle
model*: start from a random configuration of particles at `t_c`, move them at
their catalogued velocities, and whenever two meet replace them at once by the
catalogued product, until all have annihilated or time runs out; the domain
left over is the answer. The model rests on five stated simplifications: the
dynamics before `t_c` matter only through the particle distribution they
produce; particles have zero width; only two particles interact at a time;
interactions are instantaneous; and where the outcome depends on the
particles' phases, it is drawn at random with measured probabilities.

In the 1996 version, with particles placed at random, the model matched a
synchronisation rule φ_sync2 and the density rule φ_dens1 within 1% and 3%, but
underestimated the better density rule φ_dens2 by about 23%, because the
distances between particles, which the model ignored, encode the sizes of the
regions. In the 1998 version the initial particle configuration is generated
by the domain-particle transducer with measured transition probabilities. Over
five density rules and five synchronisation rules (their Table 2, 10 × 10⁴ ICs
each at N = 149) the model is within 5% of the automaton for nine, and 7.6% off
for φ_dens4; for the best density rule, φ_dens5, the automaton scores 0.7702
and the model 0.7689. Its average condensation time at N = 149 is
about 12 steps. (φ_dens5's printed look-up table is the φ_100 of Crutchfield and
Mitchell 1995, which differs from φ_par^a in two hexadecimal digits; the
decisions above judge the difference a transcription slip, since both versions
reach the published performance.) Crutchfield, Mitchell and Das
(1998, Table 6) compare six rules along the evolutionary line that led to
φ_par^a: model and automaton agree within a few per cent except for one
intermediate rule (0.691 against 0.747), whose long-lived transient domain is
missing from its catalog. They conclude that the particle-level description
captures how the automata compute.

Later, Hordijk, Shalizi and Crutchfield (2001) proved an upper bound on the
number of distinct products that particle interactions in a one-dimensional
automaton can generate, controlled by a measure of how much space-time
information a particle stores, and found it tight on several automata.

### Limits of the task

No two-state automaton of finite radius classifies density perfectly for all
N (Land and Belew 1995). Crutchfield, Mitchell and Das (1998) give the best
known radius-3 rules as reaching P_149 ≈ 0.85, and in their related-work
section report a genetic-programming rule by Andre, Bennett and Koza at 0.828
and a coevolved rule by Juillé and Pollack at 0.863 (whose performance falls
faster with N than GKL's).

### What Chemart reproduces

- **The catalogs.** φ_par^a's and φ_par^b's domains, particles, velocities and
  interaction tables are embedded exactly as printed, and
  `generate_network` returns them (tests
  `test_published_particle_catalog_is_reproduced_exactly` and
  `test_published_interaction_table_is_reproduced_exactly`). A test also
  checks that every interaction composes its reactants' walls, as in the
  `beta + gamma -> eta` example.
- **The look-up tables.** The GKL table printed by Crutchfield and Mitchell
  (1995) expands bit for bit to GKL's definition, which checks the convention
  used to read all the tables.
- **The velocities.** `probe` measures all five published velocities within
  0.05 and sees `alpha` decay into `gamma + mu` (tested).
- **Condensation.** Runs condense into domains and particles within a few
  steps: `t_c` is 5 in the default run and 6 to 9 on 599 cells, against the
  published average of about 12 for φ_dens5 at N = 149. Chemart's `t_c` uses
  its own width criterion (every wall at most 14 cells wide), not the
  published transducer, so the two are not directly comparable; the test only
  checks that `t_c` is positive and at most 30 in the default run.
- **The interactions fire.** A slow test runs six ICs on 599 cells and checks
  that all six published interactions are observed and that published
  interactions make up more than 80% of the observed two-particle events; in
  the run above every two-particle collision between catalogued particles was
  a published one. A faster test checks that no two-particle collision between
  catalogued particles ever contradicts the table.
- **Performance.** Slow tests check P_149 against 0.775 (φ_par^a) and 0.816
  (GKL) within 0.08 on 250 ICs, and P_599 against 0.740 within 0.10. The runs
  above gave 0.783 and 0.795 on 1,000 ICs, and 0.795 at N = 599 on 200 ICs,
  which is above the published 0.740 by about two sampling errors.

Not reproduced: the genetic algorithm itself and the evolutionary history of
φ_par^a; the automatic discovery of domains (the *ε-machine reconstruction* of
computational mechanics; here the domains are given); the embedded-particle
model and its predicted performances; interaction-result probabilities and
three-particle interactions; running φ_par^b, whose printed table does not
work; and the synchronisation rules, whose printed tables did not synchronise
in Chemart's implementation.

## Further reading

- Crutchfield, J. P. & Hanson, J. E. (1993). Turbulent pattern bases for
  cellular automata. *Physica D* 69:279–301. The general construction of the
  domain filter.
- Das, R., Mitchell, M. & Crutchfield, J. P. (1994). A genetic algorithm
  discovers particle-based computation in cellular automata. In *Parallel
  Problem Solving from Nature — PPSN III*, LNCS 866, 344–353. Springer.
- Gács, P. (1985). Nonergodic one-dimensional media and reliable computation.
  *Contemporary Mathematics* 41:125. The reference Crutchfield, Mitchell and
  Das give for the purpose of the GKL rule.
- Land, M. & Belew, R. K. (1995). No perfect two-state cellular automata for
  density classification exists. *Physical Review Letters* 74(25):5148.
- Andre, D., Bennett, F. H. III & Koza, J. R. (1996). Evolution of intricate
  long-distance communication signals in cellular automata using genetic
  programming. In *Artificial Life V*. MIT Press.
- Juillé, H. & Pollack, J. B. (1998). Coevolutionary learning: a case study.
  In *Proceedings of the Fifteenth International Conference on Machine
  Learning*.
