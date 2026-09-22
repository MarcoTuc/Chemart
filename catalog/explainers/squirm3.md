## Introduction

Squirm3 is a two-dimensional world of "atoms" that wander at random, bump into
each other and, when a rule says so, bond together or break apart. Tim Hutton
built it in the early 2000s to show that a molecule can copy itself in such a
world by nothing more than local collisions: no program reads the molecule, no
machine assembles the copy. Place a short chain of atoms in a soup of loose
atoms and, a few thousand steps later, there are several identical chains.

Each atom has a *type*, one of the letters `a` to `f`, which never changes, and
a *state*, a small number that reactions do change (the book suggests reading
it as something like energy content). A bonded group of atoms is a *molecule*, and moves as one piece. Eight
reaction rules, each looking at just two neighbouring atoms, are enough for any
chain that starts with an `e` atom, ends with an `f` atom and carries any
sequence of `a`, `b`, `c` and `d` atoms in between to make a copy of itself
from the loose atoms around it. The middle letters play no part in the
copying, so they can carry information, like the bases of DNA.

Hutton (2002) took up Bedau and colleagues' open challenge to "achieve the
transition to life in an artificial chemistry in silico", with the long-term
goal of seeing complexity grow over evolutionary time. He chose an artificial
chemistry because traditional cellular automata, he wrote, "do not easily
permit creatures to interact with each other". Later papers added enzymes whose function is written in the chain's sequence
(2003) and a membrane around the chain, giving cells that copy their gene,
grow and divide (2004, 2007).

Squirm3 is a simulation model, and an example of a *constructive* chemistry:
molecules of any size can form, so the set of species is not fixed in advance.
Banzhaf and Yamamoto describe it in chapter 11, among the bio-inspired
chemistries with spatial structure (book §11.4.1). Its nearest
neighbours in the catalog are the [self-replicating loops](sr-loops.md), where
the replicator is a pattern of cell states rather than a bonded object, and
the lattice membrane models of [Varela, Maturana and Uribe](autopoiesis-vmu.md)
and [Ono and Ikegami](ono-ikegami-protocell.md). In Squirm3 the central
mechanism is template replication of an arbitrary sequence. And since atoms
are never created or destroyed, the number of atoms of each type is exactly
conserved.

## How it works

### Atoms, molecules and the world

The world is a square grid. Each grid point holds at most one atom. At every
time step each atom tries to move to one of its eight neighbouring points (the
*Moore neighbourhood*); the move is allowed only if the point is empty and the
atom stays next to every atom it is bonded to. Reactions happen between an atom and
the atoms directly above, below, left or right of it (the four-point *von
Neumann neighbourhood*). The 2007 paper also uses a Moore reaction
neighbourhood with longer bonds, or continuous space with round atoms.

### Reading a rule

The rules use a notation close to ordinary chemistry, except that the number
after a letter is the atom's *state*, not a count. Writing two atoms side by
side means they are bonded; a `+` between them means they are not. The first
rule of the replicator is

```
R1: e8 + e0 -> e4e3
```

"When an `e` atom in state 8 meets a loose `e` atom in state 0, they bond and
take states 4 and 3." A rule that makes a bond is an *association*, one that
breaks a bond is a *dissociation*, and one that only changes the states of two
bonded atoms is a *transformation*. The letters `x` and `y` are *type
variables*: `x4y1 -> x2y5` applies to any two bonded atoms in states 4 and 1,
whatever their types. The same letter twice means the same type, so
`x5 + x0 -> x7x6` bonds an atom in state 5 only to a loose atom of *its own*
type. A rule ignores whatever else the two atoms are bonded to.

The eight rules of the replicator (book table 11.5, Hutton 2002 table 1) are:

```
R1: e8 + e0 -> e4e3     specific association
R2: x4y1 -> x2y5        general transformation
R3: x5 + x0 -> x7x6     general homo-association
R4: x3 + y6 -> x2y3     general hetero-association
R5: x7y3 -> x4y3        general transformation
R6: f4f3 -> f8 + f8     specific dissociation
R7: x2y8 -> x9y1        general transformation
R8: x9y9 -> x8 + y8     general dissociation
```

Loose atoms in state 0 are the food. State 1 marks a base resting in a chain,
and the other states are signals that travel along the molecule.

### Worked example: one full copy

This is the complete replication of `e8-a1-b1-f1`, taken from a run of the
default world (`chemart.evolve("squirm3", seed=6)`). A species name lists the
atoms of a molecule. For a simple chain it is the chain itself, read from the
end that sorts first. For any
other shape, such as the ladder a half-copied chain forms, the atoms are listed
in a fixed order and followed by the bonds, as pairs of positions counted from
0. Each line gives the rule that fired and the molecule-level reaction.

```
R1  e8-a1-b1-f1 + e0 -> e3-e4-a1-b1-f1
R2  e3-e4-a1-b1-f1 -> e3-e2-a5-b1-f1
R3  e3-e2-a5-b1-f1 + a0 -> a6-a7-b1-e2-e3-f1/0.1,1.2,1.3,2.5,3.4
R4  ... -> a3-a7-b1-e2-e2-f1/...
R5  ... -> a3-a4-b1-e2-e2-f1/...
R2  ... -> a2-a3-b5-e2-e2-f1/...
R3  ... + b0 -> a2-a3-b6-b7-e2-e2-f1/...
      (R4, R5, R2 again, then R3 with an f0, R4, R5)
R6  ... -> a2-a2-b2-b2-e2-e2-f8-f8/...
R7, R7, R8, R7, R7, R8 -> f1-b1-a8-e2-e2-a8-b1-f1
R7  f1-b1-a8-e2-e2-a8-b1-f1 -> f1-b1-a1-e9-e2-a8-b1-f1
R7  f1-b1-a1-e9-e2-a8-b1-f1 -> f1-b1-a1-e9-e9-a1-b1-f1
R8  f1-b1-a1-e9-e9-a1-b1-f1 -> 2 e8-a1-b1-f1
```

Step by step:

1. **Start.** R1 bonds a loose `e0` to the seed's `e8`, which becomes `e4`;
   the newcomer is `e3`. This is the first rung of a ladder.
2. **Signal.** R2 turns the bonded pair `e4`–`a1` into `e2`–`a5`: state 5
   marks the next base to copy.
3. **Pair.** R3 bonds `a5` to a loose atom of its own type, `a0`, giving `a7`
   and `a6`. In the third line above, numbering the listed atoms 0 to 5,
   bond `0.1` is this new rung `a6`–`a7` and `3.4` the rung `e2`–`e3`.
4. **Link.** R4 bonds the copy's `e3` to its `a6`, giving the copy a backbone
   of its own; R5 then turns `a7` into `a4`, and R2 passes the signal on to
   `b1`. Steps 3 and 4 repeat for `b` and `f`. Because each base pairs with its
   own type, the copy has the same sequence as the template.
5. **Unzip.** R6 breaks the rung between the two `f` atoms
   (`f4f3 -> f8 + f8`). Back up the ladder, R7 turns each bonded pair in
   states 2 and 8 into 9 and 1, and R8 breaks each rung whose atoms are both
   in state 9. The last R8 separates the two `e` atoms and leaves two chains
   `e8-a1-b1-f1`, ready to go again.

A four-atom chain takes 23 reactions and one loose atom per atom of the
chain. For a lone molecule, Hutton notes, only one reaction can apply at each
step. In a crowded world R4 can also join atoms in states 3 and 6 that belong to *different* molecules. Most
such tangles are dead ends, but Hutton found that some separate again into
working replicators of new lengths, a form of crossover.

### Floods, cosmic rays and enzymes

Two outside influences drive evolution. A **flood** every `T` steps dissolves
everything in one half (or quarter) of the world back into loose atoms in
state 0, clearing out stuck molecules and supplying fresh food. A **cosmic
ray** sets an atom's state to a random value with a small probability per atom
per step, leaving its type and bonds alone.

The later papers add **enzymes** and, with them, rules on three atoms. An
enzyme is an atom `z` whose state `i` encodes a whole reaction between two
other atoms, of types `x` and `y`: namely their states
before (`g`, `h`) and after (`j`, `k`), and whether they are bonded before
(`b1`) and after (`b2`). The formula (book eq. 11.15) packs these into one
number like the digits of a mixed-base number:

```
i = 2(2(|T|(|T|(|S|(|S|(|S| g + h) + j) + k) + x) + y) + b1) + b2 + |S|
```

Here `|T|` = 6 is the number of types (`a` = 0 to `f` = 5) and `|S|` the
number of ordinary states; the final `+ |S|` puts enzyme states above them.
Enzymes are read off the chain, whose bases are digits in base 4 (`a` = 0 to
`d` = 3). With the 2007 value
`|S|` = 38, the gene `bdca` reads `1320` in base 4, which is 120, and gives the
enzyme state 120 + 38 = 158. (In the book's version of this example,
"13204 + |S|", the final 4 is the subscript marking base 4.) Because the gene decides which
reactions its enzymes catalyse, the reactions a molecule can use become
heritable. In the 2007 cells, 41 rules copy the gene inside a loop of membrane
atoms (loose atoms can pass through the loop, larger molecules cannot), pinch
the membrane in two, read the gene out into enzymes and, rarely,
insert or delete a base (mutation).

## Using it

Squirm3 has two faces. `chemart.generate_network`, called above, returns the
*closure*: it ignores the grid, lets any molecules meet, and collects every
molecule-level reaction reachable from the seed molecule and one loose atom of
each type. `chemart.evolve` runs the world itself and records what actually
happens in it.

**The closure.** The default network above starts along the replication path:
its first four reactions are the first four steps of the worked example (R1 to
R4), and the next three already join two half-copied molecules into one. Such
tangles multiply, so the closure never completes and stops at
`max_species` molecules (default 60, status `truncated`). It takes about a
second, and the time grows faster than the budget: 100 molecules take about
five seconds. Its reactions carry no counts, since nothing is run. The closure
takes a molecule as its seed, not a `cell:`, and the world's parameters (the
grid, the food, the steps, floods and cosmic rays) belong to `chemart.evolve`
only.

**The world.** The default world is Hutton's experiment 1 set-up with the
book's example molecule: `e8-a1-b1-f1` on a 20×20 grid with 75 loose atoms of
random type, for 3,000 steps. It takes under a second.

```python
traj = chemart.evolve("squirm3", seed=1)
net = traj.network
print(net.summary())
print(net.extras["rule_counts"])
```

```
squirm3: 74 species, 67 reactions, status=observed
provides: initial-state, mass-conservation, space, stoichiometry, topology
seed: 1
extras: conservation, dissolved_by_flood, final_state, floods, reaction_rules, rule_counts, rules, space
{'R7': 34, 'R2': 25, 'R3': 25, 'R4': 23, 'R5': 23, 'R8': 14, 'R1': 9, 'R6': 7}
```

Species are named as in the worked example. The reactions are the
molecule-level events that happened, with how often (`count`) but no rate
constants, since rates come from collisions in the world. In `net.extras`:

- `rule_counts`: how often each rule fired.
- `reaction_rules`: each reaction with the rule or rules behind it.
- `final_state` and `net.initial_state`: the molecules at the end and at the
  start, with counts.
- `conservation`: one exact conservation law per atom type, plus the total
  atom count.
- `space`: the geometry of the world; `floods` and `dissolved_by_flood` count
  the floods and the molecules they dissolved.

The run is also a sequence of frames, `traj.frames`, with time in steps: the
start, after the first step, every `steps // 200` steps from there (15 in the
default run) and the end. Each frame holds the molecules present, the
reactions that fired since the previous frame and, as
`traj.series("molecules")`, the number of molecules of more than one atom.

```python
print(len(traj.frames), traj.times()[:4], traj.times()[-1])
print(traj.frames[2].fired)
```

```
202 [0.0, 1.0, 16.0, 31.0] 3000.0
[[['e8-a1-b1-f1', 'e0'], ['e3-e4-a1-b1-f1'], 1], [['e3-e4-a1-b1-f1'], ['e3-e2-a5-b1-f1'], 1], [['e3-e2-a5-b1-f1', 'a0'], ['a6-a7-b1-e2-e3-f1/0.1,1.2,1.3,2.5,3.4'], 1]]
```

Between steps 1 and 16 the seed took the first three steps of the worked
example: R1, R2 and R3. Over the whole run it splits into two copies four
times, but at step 3,000 every copy is busy copying itself again or tangled
with another.

**The 2002 experiment 1 molecule.** Hutton started from `e8-a1-b1-c1-f1` and
showed the world after 3,544 steps. The same settings, over five seeds:

```python
g = "e8-a1-b1-c1-f1"
for s in range(1, 6):
    net = chemart.evolve("squirm3", seed=s, seed_molecule=g, steps=3544).network
    splits = sum(r["count"] for r in net.extras["reaction_rules"]
                 if r["products"].get(g) == 2)
    print(s, splits, net.extras["final_state"].get(g, 0))
# seed, times the molecule split into two copies, finished copies at the end
# 1 6 0
# 2 1 1
# 3 1 0
# 4 2 0
# 5 5 0
```

Each run takes about a second. Floods (`flood_period`, `flood_sectors`),
cosmic rays (`cosmic_ray`) and a 100×100 world set up experiments 2 and 3,
with the published values in the parameter table; at tens or hundreds of
thousands of steps, those runs take minutes to hours.

**The 2007 cell.** `rules="membrane"` loads the 41 rules of Hutton (2007),
which need `n_states` of at least 38. The seed prefixed with `cell:` becomes
the gene of the paper's starting cell: a loop of 18 membrane atoms `a36` with
two anchors `a37`, one bonded to each end of the gene.

```python
net = chemart.evolve("squirm3", seed=5, rules="membrane", n_states=38,
                     seed_molecule="cell:e1-b1-c1-a1-f1",
                     space="lattice-moore", width=16, height=16,
                     food=90, steps=400).network
print(net.extras["rule_counts"])
# {'R35': 11, 'R6': 3, 'R4': 3, 'R5': 2, 'R7': 2, 'R8': 2, 'R9': 2, 'R1': 1, 'R2': 1, 'R3': 1}
```

This takes about a second. R1, R6, R2 and R3 are the opening moves the paper
describes (a second `e` atom joins the membrane), R4 to R9 copy bases, and R35
lets the membrane gain and lose atoms. 2,500 steps on a 20×20 grid with 150
loose atoms take about ten seconds. Finally, `rules="custom"` with `rule_text`
runs your own pair rules, written in the notation above, in either face.

## Results

**Replication of any sequence (2002).** Hutton's central result is that the
eight rules copy any chain of the form `e8 {x1}* f1` (an `e8`, any number of
bases in state 1, an `f1`) in a soup of atoms in state 0, which he wrote as
`e8 {x1}* f1 + {x0}* -> 2 e8 {x1}* f1`. Chemart's tests
check that `e8-a1-f1`, `e8-a1-b1-f1` and `e8-c1-a1-d1-f1` each split from a
double strand into two copies of themselves, and that R1, R3 and R4 build
loose atoms into molecules.

**188 reactions out of 7,200.** Written out without type variables, the eight
rules become 188 explicit reactions (different variable letters may take the
same type). With six types and ten states there are 60 × 60 × 2 = 7,200
possible inputs (two atoms, bonded or not), so the other 7,012 do nothing. For
`n` bases the paper prints "5n² + 21n + 4" needed reactions, but expanding the
rules gives 5n² + 21n + 24, which is the paper's own 188 for n = 4: the
printed constant is a typo. The tests check 188 and the corrected formula for
n = 1 to 4.

**Experiment 1: copies and crossover.** On a 20×20 grid with 75 loose atoms,
Hutton's `e8-a1-b1-c1-f1` had produced 11 copies after 3,544 steps: eight
finished and unable to copy further because no `e0` was left, three stuck
halfway for lack of `b0`. In another run two molecules tangled after 435 steps and
later separated into four molecules that were not all identical: two
`e-a-b-c-f`, one `e-a-b-f` and one `e-c-a-b-c-f`, a base `c` having moved from
one to the other. Chemart does not reproduce the 11 copies. In the runs shown
under *Using it*, the molecule splits one to six times, but by step 3,544 most
copies are stuck mid-copy or tangled. Crossover is not tested.

**Experiment 2: the shortest replicators win.** With a flood every T = 2,000
steps, alternating between the two halves of the world, shorter replicators
took over, because they spread faster before the next flood. By 34,500 steps
the world was dominated by the shortest possible replicators. A 100×100 run
with T = 20,000 shows the length of a seven-atom ancestor's descendants
falling over the generations. Then evolution stopped: the bases had no
function, and Hutton concluded that the chemistry needed more reactions. Chemart's tests check only that floods happen on
schedule and dissolve molecules; the selection experiment is not run.

**Experiment 3: replicators from nothing.** A 100×100 world of loose atoms in
state 0 does nothing, since no rule applies. With cosmic rays at
p = 0.00001 per atom per step, and floods every 10,000 steps to clear inactive
clumps, the first replicators appeared at around 400,000 steps and the
reaction rate jumped. Chemart's test checks only the start: a world of loose
atoms does nothing in 400 steps, while a high cosmic-ray rate (0.01) makes
reactions happen. The spontaneous replicator, which would need hundreds of
thousands of steps, is not reproduced.

**Enzymes and cells (2003–2007).** Free-floating enzyme producers, Hutton
(2007) explains, also help parasites, shorter molecules that replicate faster,
and his earlier system typically ended in global extinction; a membrane keeps
the enzymes for the cell that made them. In the 2007 paper's first experiment, rule R3
(`e6e3 -> e2e3`) was removed, so cells could divide only with an enzyme coded
in their gene. On a 70×70 world, with a one-in-a-million mutation chance per
opportunity and floods every 8,000 steps, an overnight run logged 15,772
divisions and 145 distinct genomes. The necessary gene was kept against
mutation, and the ancestor was finally replaced by a neutral variant, with an
extra `a` in front, coding for the same enzyme. In the second, cells with an
extra gene that let them eat a second kind of food outcompeted the ancestor
despite reproducing more slowly. Hutton judged the system not evolvable
enough: a new enzyme needs about 14 bases, and unused bases tend to be lost.

Chemart implements the 41 rules and the enzyme encoding. The tests check the
formula against the book, that it inverts exactly, the two enzymes worked out
in the 2004 paper (`d5731` codes for `e0 + a0 -> e2a3` and `d6740227` for
`a8 + a0 -> a8a7`, with 18 states) and the 2007 read-outs (`bca` gives 62,
`bdca` 158). The cell test checks the starting cell's shape and that division
starts (R1, R6, R2, R3 and R35 fire). Longer runs on a 20×20 grid with 150
loose atoms (seed 5 for 2,500 steps, seeds 1 and 2 for 3,000) get as far as
the last membrane rule, R34, but each ends as one bonded structure of 38 to 46
atoms rather than two cells, and no enzyme is ever read out (R37 never fires).
Chemart does not forbid crossed bonds, which the 2007 lattice does, and runs
none of the evolution experiments.

**Conservation.** Every reaction conserves the atoms of each type, and the
tests check that each vector in `extras["conservation"]` is an exact
conservation law of the network, which makes Squirm3 a positive control for
conservation-law detection.

**Later work.** Hutton also evolved biosynthetic pathways of increasing
length in a related chemistry (2003) and built an "enzyme artificial
chemistry" in which the atoms of a replicating string carry reactions (2005).
The book recommends Lucht (2012) for an account of the Squirm3 models.

## Further reading

- Hutton, T. J. (2003). Simulating evolution's first steps. In W. Banzhaf,
  T. Christaller, P. Dittrich, J. T. Kim and J. Ziegler (eds.), *Proc. Seventh
  European Conference on Artificial Life*, 51–58. Dortmund, Germany.
- Hutton, T. J. (2003). Information-replicating molecules with programmable
  enzymes. In *Proc. Sixth International Conference on Humans and Computers*,
  170–175. University of Aizu, Aizu-Wakamatsu, Japan.
- Lucht, M. W. (2012). Size selection and adaptive evolution in an artificial
  chemistry. *Artificial Life* 18(2), 143–163.
- The squirm3 program, running in the browser: <https://timhutton.github.io/squirm3>
