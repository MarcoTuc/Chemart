## Introduction

Bondable Cellular Automata (BCA) is an artificial chemistry in which every atom
is a small running dynamical system, and whether two atoms react is decided by
watching that system rather than by looking up a table. It was proposed by
Mark Hatcher, Wolfgang Banzhaf and Tina Yu at Memorial University of
Newfoundland, in a paper at the European Conference on Artificial Life
(ECAL 2011). It is a simulation model, and at the stage the book describes it
is an early one: the paper presents the model and some first measurements, and
leaves full reaction runs to future work.

The atom is a **one-dimensional cellular automaton** (CA): a ring of, say, 12
cells, each 0 or 1, all updated together at every tick by one rule that looks
at a cell and its two neighbours. There are exactly 256 such rules, numbered
0 to 255 in Stephen Wolfram's scheme, and in BCA the rule *is* the atom's
type, so there are 256 elements. From its cells an atom gets a **polarity**,
the number of 1s minus the number of 0s, which runs from −12 (all zeros) to
+12 (all ones). As the CA runs, the polarity changes, and its average over time
settles to a value that is a property of the rule. Atoms whose average
polarities have opposite signs may bond, like opposite charges. Once bonded,
cells of one ring read their neighbours from the other ring, so the two CAs
disturb each other's dynamics. That can change their polarities, and a bond
can later break without anything having hit the molecule.

The paper's motivation is stated in its introduction. In most artificial
chemistries the reaction rules are written by hand and the atoms have "little
or no internal dynamic". Giving each atom and molecule an internal dynamical
system creates a feedback loop: reactions change the systems, and the systems
decide the reactions. Because polarity can be measured for a single atom or a
molecule of any size (the authors call it "resolution independent"), any two
bodies can in principle be tested against each other with the same rule.
Banzhaf and Yamamoto call this kind of model a *subsymbolic* artificial
chemistry, an idea first discussed for random Boolean networks in RBN-World
(Faulconbridge et al.). Sorting the 256 rules by their settled polarity gives
the book's figure 10.16, a "periodic table" of CA atoms.

BCA sits in the book's chapter on automata and machines, §10.7.3, next to the
other chemistries built on cellular automata. Its closest neighbour in the
catalog is [RBN World](rbn-world.md): the same
subsymbolic idea, with random Boolean networks as atoms and properties such as
cycle length as the bonding criterion. [Embedded particles in cellular
automata](ca-embedded-particles.md) (§10.7.2) also uses CAs, but the other way
round: there the whole world is one CA and the "molecules" are moving
structures inside it, while in BCA every molecule contains its own CAs.
Unlike [Avida](avida.md) or [Tierra](tierra.md), BCA's atoms are not programs
and nothing replicates.

## How it works

### An atom and its polarity

An atom is a ring of `w` cells (12 by default) run by one Wolfram rule. The
rule number, written in binary, is a lookup table: for each of the eight
possible neighbourhoods (left, self, right) it gives the cell's next state.
Rule 0 turns every cell to 0; rule 255 turns every cell to 1; rule 204 leaves
every cell as it is; rule 90 sets each cell to the exclusive-or of its two
neighbours.

The **polarity** of a configuration is `ones − zeros`, on the scale `−w … +w`.
A ring of 12 cells has only 4,096 configurations, so from any start the CA
eventually revisits a state and cycles for ever. Chemart runs each atom from a
fixed start, half the ring set to 1 in one contiguous block (`111111000000`,
polarity 0), until the cycle closes. The **settled mean polarity** is the
average polarity over that cycle, and the **settling iteration** is how many
ticks it took to enter the cycle. Rule 90, for instance, settles after two
ticks into a two-state cycle whose states both have eight 1s, so its settled
mean polarity is 8 − 4 = +4. Rules that preserve the number of 1s, such as the
identity rule 204 and the "traffic" rule 184, keep the balanced start's
polarity and settle at 0.

That start is not stated in the book. It was recovered from figure 10.16: of
all 4,096 possible starts, only the half-ring block (up to rotation) puts all
256 rules within half a column of where the figure draws them. The original
paper used a different start, a single live cell; see Results.

### Bonding strength

Two atoms that meet are compared cell by cell. The book's figure 10.15 marks
the **bonding strength** as the "largest contiguous sequence of complementary
binary numbers": the longest unbroken run of positions where one ring has a 1
and the other a 0. The book does not say how the rings are lined up, so
Chemart tries all `w` rotations of one ring against the other and keeps the
best. The result is symmetric and at most `w`. A bond forms only if this run
is at least `bond_threshold` cells long (10 by default).

### Coupling: what a bond does

The book says only that once bonded, "each cell in an atom is influenced by the
neighborhood cells of the other (bonded) atom". Chemart's reading: inside the
bonded run, each cell takes the neighbour on the partner's side from the
partner's facing cell instead of from its own ring. The neighbourhood stays
three cells wide, so the 256 rules still apply unchanged, but information now
flows between the rings. The run is recomputed at every tick from the current
cells, so a bond can lengthen or shorten as the coupled system runs.

### Molecules, reactions and the reactor

A molecule is a chain of atoms bonded end to end. Two kinds of reaction exist:

- **Association**, `A + B → AB`. Two molecules join when their polarities pass
  the polarity test (by default, opposite signs; zero counts as neither, so a
  polarity-0 molecule is inert) and the facing atoms at the joined ends reach
  the bonding strength. The new molecule is then run, coupled, until it
  settles.
- **Dissociation**, `AB → A + B`. If, after settling, a bond's run has fallen
  below `bond_threshold`, the molecule breaks there and each fragment is run
  again on its own.

A molecule's species name lists its atom types, such as `r165-r90`. Because
the cells matter, the same skeleton can occur in several cell states; these
are numbered `r165-r90#2`, `r165-r90#3`, and so on. A single atom can also
reappear in a new state after a molecule breaks, as `r90#2`.

The model's reactor is a well-stirred pot in which molecules collide at
random. Chemart does not sample collisions. It starts from one atom of each
seeded type and applies every possible association and dissociation to every
molecule found so far, until nothing new appears: the reaction *closure* of
the seed set. There are no rates, since the book gives none.

### A worked example

The default run seeds eight atom types: 0, 18, 90, 110, 165, 184, 204 and 255.
Its first reaction is the simplest possible bond:

```
r0 + r255 -> r0-r255
```

Rule 0 settles to all zeros (polarity −12) and rule 255 to all ones (+12). The
signs are opposite, and every one of the 12 facing cells is complementary, so
the bonding strength is 12. Neither rule reads its neighbours, so the coupling
changes nothing and the bond stays at 12 for ever. The molecule's polarity is
the average over its atoms, (−12 + 12) / 2 = 0, so under the default test it
is inert: it never bonds again.

The next reactions show a bond that changes itself:

```
r90 + r165 -> r165-r90
r165-r90 -> r165#2 + r90#2
r90#2 + r165#3 -> r165-r90#3
```

Rule 90 settles at `110011110011` (+4) and rule 165 at `001100001100` (−4).
These two rings are exact complements, so they bond with strength 12. Once
coupled, both rings start to change. After Chemart's default budget of 256
ticks the pair has not yet settled, and at that moment the bond's run is only
3 cells, below the threshold of 10, so the molecule breaks into two atoms in
new states, `r165#2` and `r90#2`. Two such fragments can then bond again,
giving `r165-r90#3`, whose run of 10 survives and which settles at polarity 0.

That decay depends on the budget. Given 2,048 ticks, the same coupled pair
settles after 1,830 ticks with a run of 10, and does not break (see *Using
it*). The book says bonds may "increase" or "decrease"; in the default run the
decreases all come from molecules that are still in their transient when the
budget runs out.

### Two details of Chemart's code

The polarity test at association compares each molecule's polarity in the
single configuration Chemart stores for it (the state at which its cycle
closed), not its mean over the cycle. For most seeded atoms the two agree. For
rule 18 the cycle mean is −5.14 and the stored state reads −8; for rule 110 the
mean is +0.67 but the stored state reads −2, so rule 110 is treated as negative
although the periodic table places it on the positive side.

Species identity treats a chain and its reverse as the same molecule, but the
coupled dynamics do not: bonding atom A on the right of B is not the mirror
image of bonding it on the left, because the rings are not reflected. That is
why the default network contains both `r90 + r165 -> r165-r90` and
`r165 + r90 -> r165-r90#2`: the same pair, met in the other order, settles to
a different state.

## Using it

The default call above builds the closure of the eight seed types described
in the worked example: 23 species and 11 reactions (7 associations, 4
dissociations), in about two seconds. It reproduces no particular published
reaction experiment, because none was published; the paper only tested atoms
in isolation and in pairs. What it does reproduce is the periodic table of
figure 10.16, which every run computes for all 256 rules and stores in
`net.extras["analysis"]`:

```python
import chemart
from collections import Counter

net = chemart.generate_network("bondable-ca", seed=1)
table = net.extras["analysis"]["periodic_table"]      # 256 rows, sorted by polarity
cols = Counter(round(row["mean_polarity"]) for row in table)
print(cols[-12], cols[0], cols[12])                   # 20 54 20
print(max(row["settling_iteration"] for row in table))  # 28
```

Each row gives `rule`, `mean_polarity`, `settling_iteration`, `cycle_length`
and `settled_config` (cell 0 first). For example rule 110 has mean polarity
0.667, settling iteration 28 and a cycle of 9 states.

Other fields of `net.extras` worth knowing:

- `analysis["molecule_polarity"]`: the polarity of each species, averaged over
  its atoms (`r0-r255`: 0.0, `r165-r90#3`: 0.0).
- `analysis["bond_strengths"]`: the current run length of each bond in each
  molecule (`r165-r90`: `[3]`, `r165-r90#3`: `[10]`).
- `analysis["unsettled_molecules"]`: molecules whose coupled dynamics did not
  settle within `settle_iterations` (default run: `r165-r90`, `r165-r90#2`).
- `conservation`: one law per seeded atom type. Reactions only rearrange
  atoms, so the number of atoms of each rule is conserved.
- `interaction_law`: the settings of this run in words.
- Each species' `structure` holds the whole molecule, e.g.
  `r0:000000000000~+0~r255:111111111111`: rule and cells of each atom, and the
  rotation offset of each bond.

The generator is deterministic; the seed changes nothing.

### Recipes

**Let molecules settle.** With a larger budget the default seeds give only the
three associations and no decay:

```python
net = chemart.generate_network("bondable-ca", seed=1, settle_iterations=2048)
# 11 species, 3 reactions:
#   r0 + r255 -> r0-r255
#   r90 + r165 -> r165-r90
#   r165 + r90 -> r165-r90#2
```

**Switch the interaction off.** `coupling="none"` keeps bonding but stops the
rings from influencing each other, so nothing ever decays: 10 species and two
reactions (`r0 + r255 -> r0-r255`, `r90 + r165 -> r165-r90`).
`coupling="xor"`, the other reading of the book (the partner's cell is
exclusive-ored into the neighbour rather than replacing it), gives 18 species
and 10 reactions and is the only one of the three that builds a three-atom
molecule from the default seeds (`r0-r165-r90`, `r165-r90-r90`).

**Loosen the rules.** Each of these hits the species budget
(`max_species=200`, status `truncated`), because far more pairs can bond:

| setting | time | reactions | associations |
|---|---|---|---|
| `polarity_rule="any"` | 12 s | 175 | 103 |
| `polarity_rule="sum-zero"` | 4 s | 164 | 98 |
| `bond_threshold=9` | 5 s | 157 | 115 |
| `bond_threshold=8` | 9 s | 177 | 116 |
| `max_atoms=4` | 10 s | 139 | 67 |
| `settle_iterations=2048, bond_threshold=8` | 56 s | 203 | 147 |

Wide seed sets behave the same way: every eighth rule
(`atom_types=list(range(0, 256, 8))`) truncates at 200 species in 2 s.
Larger `settle_iterations` makes loose settings much slower, since each coupled
molecule may run the whole budget. In the truncated runs,
`unsettled_molecules` also lists molecules that were built but did not make it
into the 200 species.

**Reproduce the paper's isolated-atom tally.** The helper `settle` runs one atom
from any start. From the paper's start, a single live cell, the rules split as
in the paper's Table 3:

```python
from collections import Counter
from chemart.chemistries.bondable_ca import settle

signs = Counter()
for rule in range(256):
    m = settle(rule, 12, config=1)["mean_polarity"]
    signs["positive" if round(m) > 0 else "negative" if round(m) < 0 else "neutral"] += 1
print(signs)    # Counter({'negative': 115, 'positive': 75, 'neutral': 66})
```

The parameter table below lists every setting and its range.

## Results

### What the paper showed

Hatcher, Banzhaf and Yu (2011) built the model and then asked one question
about its observable: would mean polarity be too smooth, "essentially reducing
in most cases to a static value", which would freeze the chemistry? They first
tried **instant polarity**, the polarity at the moment of collision, and found
it "too stochastic for some transition rules": for rule 30 it jumps erratically
from tick to tick, so whether two bodies bonded depended on the random moment
of collision. The running mean smooths this out while staying cheap to
compute, so they adopted it. They also rejected cycle length, the observable
used in RBN-World, because it only changes when a body reacts and costs a
full simulation to measure.

**Isolated atoms.** All 256 rules were run on 12-cell rings for 4,096 ticks
from a single live cell (polarity −10), with the mean taken from the first
tick. For 141 of the 256 rules (56%) the sign of the mean changed during the
run. At the end 75 rules (29%) were positive, 115 (45%) negative and 66 (26%)
neutral, and so inert. The authors noted that most of the neutral rules took
more than 256 ticks to get there, leaving time for them to bond first, and
that the further a rule's final polarity lay from its start, the longer it
tended to take to settle. The paper's appendix draws these results as a "periodic
table" of the rules.

**Bonded pairs.** They then tried to bond every pair of different rules, 32,640
pairs, and ran each bonded pair on. 8,625 pairs (25%) bonded. In 81% of the
bonded atoms (14,049) bonding changed the atom's mean polarity, and in 34%
(5,799) it changed its sign; the pair's own mean polarity changed in 91% of
pairs and its sign in 45%. Bonding also moved pairs away from neutrality,
leaving fewer inert bodies than isolated atoms would. From this they
concluded that bonding "could be a self-sustaining process, keeping the system
active", rather than locking atoms into inert structures. (The paper gives
the length of these pair runs as 10,000 ticks in the text and 100,000 in its
Table 2.)

The paper's model differs from Chemart's in several ways: bonds join the
longest run of 1s in the positive atom to the longest run of 0s in the
negative one, with no threshold; a linked cell sees its partner's whole
neighbourhood; bonds break when the bonded atoms no longer have opposite
polarity; and molecules are binary trees in which two molecules may bond
through up to two pairs of atoms. The paper leaves "full simulation runs",
with many molecules reacting to form larger bodies, to further work. The book
cites no later publication on BCA; its figure 10.16 comes from unpublished
work by Hatcher (2013).

### The book's periodic table, and what Chemart reproduces

Figure 10.16 of the book orders all 256 rules of width 12 by the value at which
their mean polarity settles, from −12 to +12, and shades each box by how long
it took to settle, in buckets up to 2,048 ticks. Twenty rules settle at −12
(the ring dies to all zeros), twenty at +12 (all ones) and 54 at 0; the 54 are
printed as a block below the table because they do not fit in one column.

- **The periodic table: reproduced.** The tests compare Chemart's value for
  every rule with the figure as read off the book's page: all 256 lie within
  0.5 of the figure's column and 226 match exactly. The other 30 have a
  fractional cycle mean, drawn at the nearest column; eight of them sit
  exactly halfway, and the figure rounds four toward zero and four away, so no
  single rounding reproduces it. The tests also check the 20/54/20 split and
  that the half-ring start is what makes the figure come out (random balanced
  starts and the alternating start fall well short).
- **Density-preserving rules at 0: reproduced.** Rules 204, 170, 240, 184, 232
  and 51 keep the start's number of 1s (or, for 51, swap 1s and 0s every tick)
  and so settle at 0; the tests check each.
- **Settling within 28 ticks: reproduced, with a caveat.** Under Chemart's
  definition (the ticks before the cycle starts) no rule takes more than 28
  ticks, which the tests check and which fits inside the figure's scale. The
  figure's shading, however, spreads rules across buckets up to 2,048, which
  looks like the time for a running mean to settle, as in the paper. Chemart
  does not reproduce the shading.
- **Reactivity read off the dynamics: reproduced by construction.** Bonds are
  decided by polarity and complementary cells, never by a table.
- **Bonding strength symmetric, at most `w`: reproduced.** The tests check both
  on 400 random pairs, and that a half-ring block is fully complementary to
  itself rotated by half a ring.
- **Coupling changes the dynamics, bonds strengthen or weaken, molecules
  decay: reproduced.** The tests check that coupled pairs settle to different
  states than uncoupled ones for many pairs, and that switching coupling off
  changes the network. As the worked example shows, in the default network the
  decays come from molecules cut off by the 256-tick budget.
- **Paper's Table 3: reproduced, not tested.** From the paper's single-cell
  start, Chemart's CA gives the paper's 75 / 115 / 66 split (recipe above).
  The 141 sign changes are not reproduced (Chemart's CA gives 123 counting
  exact signs, 149 counting rounded ones).
- **Not reproduced:** the paper's bonded-pairs experiment and its bond
  mechanism, the running mean, tree-shaped molecules, collision-driven runs
  with reaction success rates and lifetimes, and any kinetics. Chemart's model
  follows the book; the paper was found only after it was written, and the
  choices it forced are listed under *Implementation decisions*.

## Further reading

- Faulconbridge, A., Stepney, S., Miller, J. & Caves, L. (2010).
  RBN-World: the hunt for a rich artificial chemistry. In *Proceedings of the Twelfth
  International Conference on the Synthesis and Simulation of Living Systems
  (ALife XII)*, 285–292. The follow-up to RBN-World that BCA cites alongside
  the original.
