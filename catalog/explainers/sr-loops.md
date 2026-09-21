## Introduction

A *cellular automaton* (CA) is a grid of cells, each in one of a few states,
that all update at once by the same table: a cell's next state depends only on
its own state and on those of its immediate neighbours. Self-replicating loops
are patterns in such a grid that copy themselves. A closed ring of cells
carries a circulating string of signals; the signals push out an arm and bend
it round until it meets its own root, and the closed-off arm becomes a second
ring, identical to the first, which starts doing the same. Nothing in the table
mentions loops. The loop exists only as a configuration of cell states, yet it
behaves like a thing that reproduces.

The line of work starts with John von Neumann, who, following a suggestion of
Stanislaw Ulam, invented cellular automata to show that a machine can build a
copy of itself. His automaton had 29 states per cell, its core design would
have occupied about 200,000 cells (book §10.7.2, after Kemeny), and it was a
*universal constructor*, able to build any machine it was given a description
of. Christopher Langton (1984) dropped universality: his 86-cell loop, in an
8-state automaton, builds only itself. Byl (1989) found a 12-cell loop in 6
states, and Reggia, Armentrout, Chou and Peng (1993) loops of 6 and 5 cells in
8 states, the smallest the book knows of (§10.4). The book notes the price:
without a general constructor, the potential for evolutionary growth of
complexity "is correspondingly less general".

Hiroki Sayama (1998, 1999) then gave Langton's loop a way to die. His
*structurally dissolvable* (SDSR) loop falls apart when it cannot complete a
copy, so a population must keep reproducing to hold its numbers. His
*evoloop* adds variation: loops that collide while replicating can produce
different offspring, so a colony shows natural selection although the
automaton is fully deterministic.

Banzhaf and Yamamoto use these loops (§10.7.2) to show that *what counts as a
molecule depends on the observer*. Seen cell by cell, a cell state is a
molecule and a state change is a reaction. Seen pattern by pattern, a whole
loop is a molecule and replication and dissolution are its reactions, but
someone has to decide how to recognise a loop. Chemart generates both
readings. The nearest catalog neighbour, [embedded particles in cellular
automata](ca-embedded-particles.md), reads molecules out of a one-dimensional
CA as moving domain boundaries; [the machine-tape chemistry](ikegami-hashimoto.md)
keeps von Neumann's split between constructor and description; and
[Squirm3](squirm3.md) gets replication from atoms that move and bond rather
than from a fixed grid of states.

## How it works

### The cellular automaton

Each cell looks at its four nearest neighbours, north, east, south and west
(the *von Neumann neighbourhood*), and all cells update together from a
*transition table*: entries of the form "a cell in state C with neighbours N,
E, S, W becomes C′". An entry also covers its three rotations, so a loop
behaves the same whichever way it faces, and a neighbourhood no entry covers
leaves the cell unchanged. State 0 is the *quiescent* background.

In Langton's family of loops, state 2 is the *sheath*, a wall on both sides of
a channel of *core* cells in state 1, along which signals travel: a 7 paired
with a 0 means "grow the arm one cell straight on", a 4 paired with a 0 "turn
left". This is the default ancestor, as Chemart stores it, one row per line:

```
022222222000000
217014014200000
202222220200000
272000021200000
212000021200000
202000021200000
272000021200000
212222221222220
207107107111112
022222222222220
```

The square on the left is the loop; the tail at the bottom right is its arm.
Signals circulate round the channel, and each one passing the arm's root is
copied into the arm and extends or turns it. After three left turns the tip
meets the root, the new loop separates, and both loops start again.

### Two readings of "molecule"

**Micro.** A cell state is a molecule (`s1` to `s7`, or `s8` for Sayama's
rules) and the quiescent state is its absence, so it is left out of formulas.
Each transition that fires is a reaction in which one cell changes and its
four neighbours act as catalysts, unchanged on both sides. The three most
frequent reactions in 151 steps of the default loop are:

```
s7 + 2 s2 + s1 -> 2 s2 + s1  (x1460)
s7 + 2 s2 + s1 -> 2 s1 + s7 + 2 s2  (x1448)
s1 + s7 + 2 s2 -> 2 s7 + 2 s2  (x1179)
```

In the first, a cell in state 7 with two sheath neighbours, one core
neighbour and one empty neighbour becomes empty (`s7` is gone on the right).
In the second, an empty cell with neighbours 7, 2, 2, 1 becomes a 1 (the
right side gains an `s1`). In the third, a core cell next to a 7 becomes a 7.
Together they move a signal: as it passes, a channel cell goes 1 → 7 → 0 → 1.
The book reads exactly this, a state vanishing in one cell and appearing in
the next, as a molecule moving.

**Macro.** A loop is a molecule. After each step (or every `track_every`
steps) Chemart looks for *loops*: connected groups of non-quiescent cells,
diagonal contact included, of at least `min_loop_cells` cells. Each loop is
followed from look to look by the cells it shares with its earlier self. A new
group that split off a loop is a birth; a loop that is gone is a death.

A *species* is registered from the configuration a loop has at birth (and
from the ancestor's starting configuration); any loop that later shows a
registered configuration, in any rotation, belongs to that species. So an
exact copy is its mother's species and a mutant is a new one. Names give the
cell count at registration and a letter label: `L086aaa` is the first species,
registered with 86 cells. The reactions are:

```
L -> 2 L        replication: the daughter is the mother's species
L -> L + L'     replication with variation: the daughter is a new species
L -> ∅          death: the loop has disappeared
∅ -> L          a loop appeared without splitting off a tracked loop
```

The last kind needs a warning. The specification and decisions below say that
a loop appearing detached from every tracked loop is given the nearest living
loop as its mother; the current code instead records it as `∅ -> L`. Such
events are common in crowded runs: the SDSR recipe below has 33 of
`∅ -> L086aaa`.

The default run is Langton's loop on a 60×60 lattice for 320 steps. Its first
recorded birth is

```
{'step': 128, 'identified': 151, 'copy': 151, 'mother': 'L086aaa', 'daughter': 'L086aaa'}
```

The daughter separates at step 128, unfinished, and is provisionally a species
of its own. At step 151 it shows the ancestor's exact configuration, so it is
merged into `L086aaa` and counted as a copy: Langton's period of 151 steps.
(The mother itself shows its starting configuration again at step 147, a
matter of the phase in which the ancestor pattern is distributed.) Mother and
daughter reproduce again at steps 275 and 279, and the run condenses into one
reaction, `L086aaa -> 2 L086aaa`, fired three times.

### Space, dissolution and variation

The lattice is finite: by default it wraps round into a torus, and with
`boundary="quiescent"` everything outside it is background. That is what makes
loops compete, since each needs free room to build a daughter. When room runs
out, Langton's loops stop and become static. Sayama's SDSR table adds a ninth
state, 8, the *dissolver*, which spreads through a connected loop structure
and erases it; in Salzberg, Antony and Sayama's (2004) account of the
evoloop, it is triggered by configurations outside the normal replication
cycle and "typically arises from shortage of space due to overcrowding". An
undisturbed SDSR loop behaves exactly like Langton's. In the evoloop,
collisions during replication can also change the signal sequence an
offspring inherits, and since that sequence describes the offspring's shape,
the offspring can come out a different size. That is the only source of
variation. The run is deterministic; the seed only places several ancestors
(`ancestors`) at random positions.

## Using it

The default call above is Langton's loop replicating three times. The species'
`structure` is the loop's configuration, and the run is in `net.extras`:

```python
a = net.extras["analysis"]
a["births"], a["deaths"], a["loops_final"]    # (3, 0, 4)
a["ancestor_copy_steps"]                      # [151, 298, 302]
a["population"]                               # snapshots: loops, their sizes, count per species
net.extras["space"]["final"]                  # the final lattice, one string per row
```

With `mode="micro"` the same run gives the cell-state network instead: 7
species and 89 distinct transitions (6,347 firings) in 151 steps.

The loop detector cannot tell two touching loops from one larger loop. When a
colony packs the lattice, a loop that merges into a neighbour's group is
recorded as `L -> ∅`, and fragments can register as new species. In Langton's
and Byl's loops, which never dissolve, every recorded death is of this kind.
Raising `min_loop_cells` hides debris; the parameter table gives values per
rule.

**The small loops.** Byl's and the Chou-Reggia loops need `min_loop_cells=4`:

```python
net = chemart.generate_network("sr-loops", rule="byl", grid=40, steps=80,
                               min_loop_cells=4)
[r.to_text() for r in net.reactions]    # ['L012aaa -> 2 L012aaa  (x7)']
```

Eight Byl loops stand on the lattice at the end. With `rule="reggia-2"` the
5-cell loops fill the same torus within 80 steps, and 35 firings of
`L005aaa -> 2 L005aaa` come with touching artefacts such as `L005aaa -> ∅`.

**Freezing against turn-over.** Langton's loop and the SDSR loop on a 100×100
torus for 4,000 steps take 3 to 5 seconds each:

```python
net = chemart.generate_network("sr-loops", rule="sdsr", grid=100, steps=4000,
                               min_loop_cells=40, track_every=5)
```

The SDSR run records 274 `L086aaa -> 2 L086aaa` and 296 `L086aaa -> ∅`. In the
recorded snapshots the number of loops first reaches 28 near step 900, then
fluctuates between 8 and 28; 13 loops, all `L086aaa`, remain at the end.
Short-lived species such as `L049aac` are pieces of dissolving loops. With
`rule="langton"` the count settles at 13 by step 1,600, and between steps
3,000 and 3,100 only 20 of 4,431 occupied cells change, against 1,255 of 1,507
in the SDSR lattice.

**Evolution in the evoloop.** This run takes about 30 seconds:

```python
net = chemart.generate_network("sr-loops", rule="evoloop", grid=200, steps=30000,
                               min_loop_cells=20, track_every=25)
```

It records 2,973 births, 3,281 deaths and 923 species. At the end there are
76 loops, 39 of them species `L033atf`; species registered at 20 and 33 cells
carry most of the population. One evoloop replication takes 363 steps, so runs
of a few hundred steps show only replication. Sayama's runs, 10⁴ to 10⁷ steps
on lattices up to 1000×1000, are possible but slow.

## Results

**Replication without a universal constructor.** Langton's 86-cell loop makes
an exact copy of itself in 151 steps, Byl's 12-cell loop in 25, the 5-cell
Chou-Reggia loop in 15 (the 6-cell one in 13), and the evoloop ancestor in
363. Chemart's tests reproduce every period by running the tables, and the
published sizes: 8 states and 219 transitions for Langton, 6 states for Byl,
8 states and 5 cells for the book's smallest replicator. They also find Byl
copies at steps 25, 50 and 75, and check, for Langton's loop, that a complete
copy first stands on the lattice at step 147 and the new loop is exact at 151.

**Dissolution: from freezing to a steady state.** In a bounded space a colony
of Langton loops grows until the space is full and then stops. SDSR loops
that can no longer replicate dissolve, so, in the book's words, they "will
have to permanently keep reproducing in order to hold their numbers at a
certain concentration level" (§10.7.2). Sayama showed that in a bounded space
this makes loops compete, which leads to evolutionary effects. Chemart's tests
check that an undisturbed SDSR loop is Langton's loop, step for step, for 200
steps, and that on a 100×100 torus over 4,000 steps the SDSR colony has more
than 100 births and more than 100 deaths while the loop count stays bounded in
the second half; the recipe above shows Langton's colony freezing instead.

**Selection for smaller loops.** As evoloops compete for space, "smaller loops
that have a reproductive advantage emerge and dominate" (book §8.2.3).
Salzberg, Antony and Sayama (2004) read an evoloop's signals as a genome, G
for a straight-growth gene (`0 7 1`), T for a left-turn gene (`0 4 1`), C for a
plain core cell, and define a loop's *size* as its number of G's. A viable
loop of size n has n G's and a pair of T's with no G between them, followed by
a G; loops of size 3 or less cannot replicate. Their Figure 7 starts from a
size-8 ancestor on a 1000×1000 lattice. Chemart's evoloop ancestor, Golly's
pattern, is larger: 13 straight-growth signals around an 11×11 empty middle, a
size-13 loop of 149 cells. After the 30,000-step run above, counting the
square holes enclosed in the final lattice (done for this page; Chemart does
not report it) finds 35 with a 2×2 middle and 35 with a 3×3 middle, the
middles of size-4 and size-5 loops. Chemart's test checks coarser things: more
than 500 births and deaths, more than 20 species, variation
(`L149aaa -> L149aaa + L…`), and a dominant species and median loop much
smaller than the ancestor.

**Genetic diversity and microevolution.** Salzberg et al. (2004) estimate the
number of viable species of size n as the binomial coefficient
C(2n−2, n−2) (their Table 1): 15 of size 4, 3,003 of size 8, over two billion
by size 18. In a round robin of the 15 size-4 species, each pair started at
opposite ends of a 1000×1000 space for 100,000 updates, survival differed
significantly even between loops of equal size. Colonies grow parabolically,
"due to the geometric constraint of the 2D space", and the fitted growth
coefficient correlates with survival (0.674 for updates 0–2000); both rise
during evolution, so selection acts on more than replication time. Any
subsequence of the form G{C}T{C}TG (`{C}` being any number of C's) survives
mutation and fixes a minimum size; a run seeded with one forcing size 15 kept
changing its dominant species for over six million steps on a 401×401
lattice, with 7,106 species observed, 58 of them self-replicating. Chemart
reproduces none of this, because it names species by configuration and does
not read the G/T/C genome.

**Islands, concentrations and robustness.** The book's Figure 10.13, from
Sayama (1999), plots changing concentrations of loop species in a bounded
space, and §8.2.3 says that in large spaces the species cluster into islands.
Chemart records per-species counts in `extras["analysis"]["population"]` and
the final lattice in `extras["space"]["final"]`, from which both can be drawn,
but no test checks them. The book also credits dissolution with giving
evoloops "some degree of robustness to underlying hardware errors"; Chemart
has no way to inject such errors and does not test it. Not implemented: the
study of evoloops under hostile pathogens and mass killings (book ref [739]),
von Neumann's 29-state automaton, and the later loops of Tempesti, Perrier
and others.

## Further reading

- von Neumann, J. (1966). *Theory of Self-Reproducing Automata*. Urbana, IL:
  University of Illinois Press.
- Sipper, M. (1998). Fifty years of research on self-replication: an
  overview. *Artificial Life* 4, 237–257.
- Golly, the cellular automaton simulator whose rule and pattern files supply
  Chemart's tables and ancestors: <https://golly.sourceforge.net>
