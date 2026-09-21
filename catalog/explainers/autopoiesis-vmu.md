## Introduction

This is the first computer model of *autopoiesis*, published by Francisco
Varela, Humberto Maturana and Ricardo Uribe in 1974. Autopoiesis means
"self-production". Maturana and Varela, biologists at the University of Chile,
used the word for what they took to be the defining organisation of a living
cell: a network of reactions that keeps producing its own components,
including the boundary that holds the network together and separates it from
its surroundings. On this view a system is alive if it maintains itself. It
need not reproduce or evolve (book §6.1.5).

The model was built to show that this abstract idea can be realised in
something concrete and very simple. It is a square grid of sites, each holding
one particle. Particles wander at random and come in three kinds: *substrate*,
*catalyst* and *link*. A catalyst turns two neighbouring substrates into a
link. Links bond to one another into chains, and a chain that closes on itself
is a *membrane*. Links also fall apart at random. The picture the authors
wanted is this: a closed chain of links surrounds a catalyst. Substrate can
pass through the chain, links cannot, so links made by the catalyst pile up
inside. When a link of the membrane falls apart, one of those spare links
drifts into the gap and closes it. The "cell" keeps its shape while its parts
are replaced. The 1974 paper reported both this repair and the spontaneous
formation of such a cell from a lone catalyst in substrate.

The model has a second history. In the 1990s Barry McMullin could not make a
re-implementation work: the links inside the cell bonded to each other,
stopped moving, and could not repair anything. A printout of an early version
of the original program, found among Varela's papers, showed why. The program
contained a rule that no published description mentioned, *chain-based bond
inhibition*, which stops a free link from bonding while it sits next to an
existing chain. McMullin and Varela (1997) restored it in a new program, SCL
(Substrate-Catalyst-Link), and the self-repair came back. Chemart implements
SCL's version of the chemistry.

It is a simulation model. Banzhaf and Yamamoto present it in chapter 6, "The
Essence of Life", next to [the chemoton](chemoton.md), Gánti's minimal cell
(§6.1.4). The chemoton couples a metabolism, a replicating template and a
growing membrane through exact reaction equations, with no space. This model
has no genetic part and no reproduction, and everything happens in space: its
reaction list is three lines, and the boundary is a physical ring of
particles. The closest relative in the catalog is the
[Ono–Ikegami protocell](ono-ikegami-protocell.md), a later lattice model in
which membranes form from repelling particle types and the cells can grow and
divide. Unlike reaction-network models such as
[Kauffman's autocatalytic sets](kauffman-autocatalytic-sets.md) or
[AlChemy](alchemy.md), which have self-producing sets of molecules but no
space, this one has a boundary built by the reactions it encloses. That
boundary is what the model was built to show.

## How it works

### The particles

The world is a square lattice that wraps around at the edges (a *torus*), so
it has no border. Every site holds exactly one thing:

- a **hole**, an empty site;
- a **substrate** particle (written `O` in the book, `S` in Chemart), the raw
  material;
- a **catalyst** (`*`, `C`), which is never made or destroyed;
- a **link** (`[O]`), the membrane material.

A link can hold up to two **bonds** to neighbouring links, so links form
chains. Chemart names a link by its number of bonds: `L0` is a free link, `L1`
the end of a chain, `L2` a *chain link* in the middle of a chain. A link can
also carry one substrate particle inside it (the next section says why), and
then gets an `S` suffix: `L0S`, `L1S`, `L2S`. A closed chain of six or more
links is a membrane; a smaller closed chain is a *cluster*, which is too small
to enclose anything. A membrane around a catalyst is a *cell*.

### What happens in one time step

Each time step every particle gets one turn, in random order. On its turn a
particle first may **move**: it picks one of its four side neighbours (up,
down, left, right, the *von Neumann* neighbourhood) and swaps places with it
with probability `sqrt(m_i × m_j)`, where `m` is a mobility factor per kind of
particle. Substrate and holes move easily (factor 0.5), catalysts and free
links slowly (0.1). A bonded link never moves. Then the particle acts,
according to its kind. Actions look at all eight surrounding sites (the
*Moore* neighbourhood, which includes diagonals):

- **Production.** A catalyst with two substrates among its eight neighbours
  turns one of them into a link and the other into a hole:
  `* + 2 O -> * + [O]`. The catalyst is unchanged.
- **Disintegration.** A link falls apart with a small probability per step
  (0.01 by default). It becomes substrate again, and a second substrate goes
  into a neighbouring hole if there is one. Its bonds disappear with it.
- **Bonding.** A link with fewer than two bonds picks one neighbour at random.
  If that neighbour is also a link with a free bond site, the two bond, unless
  the inhibition rule below forbids it.
- **Absorption and emission.** A link can take in a neighbouring substrate,
  leaving a hole, and later release it into a neighbouring hole on any side.

The last pair is how a membrane lets substrate through. A bonded link cannot
move, so nothing can squeeze past it; but a substrate on the outside can be
absorbed by a membrane link and emitted on the inside, where the catalyst can
use it. Because motion only goes sideways, while a chain can run diagonally,
a closed chain whose consecutive links touch at least at a corner seals off
its interior for every moving particle. So a membrane is permeable to
substrate and impermeable to catalysts and links.

### Chain-based bond inhibition

The rule missing from the 1974 paper is a restriction on bonding: **a free
link may not bond while a chain link (`L2`) is among its eight neighbours.**
Only free links are affected. A link at the end of a chain (`L1`) can still
bond.

The consequence is the whole point. The links produced inside a cell are
always next to the membrane, whose links are chain links, so they stay free and
mobile. When a membrane link disintegrates, the gap is flanked by two chain
ends, which are `L1` and do not inhibit anything. A free link that wanders into
the gap bonds to both ends and closes the membrane. Without the rule, the free
links inside bond to each other and to the membrane as soon as they meet,
become immobile, fill up the cell, and none is left to close a gap.

Here is the same ready-made cell in Chemart after 60 steps, with disintegration
switched off, with and without the rule (seed 1, a 15 × 15 world). In the grid
`o` is substrate, `.` a hole, `*` the catalyst, `L` a free link, `l` a link
with one bond, `=` a chain link, and `M`, `m`, `#` the same three holding an
absorbed substrate:

```
with the rule           without it
o.oooo=##M.o.oo         oooooo###mooooo
ooooo#*ML#oo.oo         .o.oo##*l#ooooo
ooo.o#L.M#ooooo         ooooo###.#ooooo
oo.o.#MLM#ooooo         ooooo#...#.oooo
o.oo.o###oooooo         oo.ooo##=oooooo
```

These are rows 5–9 of the grid. The membrane started as a ring of twelve `=`
around a 3 × 3 interior (rows 6–8, columns 6–8). With the rule, seven links sit free inside the
cell (`L`, `M`) next to the catalyst. Without it, the four links made inside
have bonded to the membrane and to each other, and no free link is left.

One side effect of the square lattice shows here too. Where the ring turns a
corner, the site diagonally outside the corner is one of the catalyst's eight
neighbours when the catalyst sits in the inner corner, so a catalyst inside a
cell can make a link *outside* it. The `M` at the top right of the left grid
and the `m` at the top right of the right grid were made that way.

### A worked example: reading the reactions

A run records every reaction that fired, with a count. These are from the
default run below:

```
C + 2 S -> C + L0  (x39)
L0 + S -> L0S  (x45)
L1 + L0 -> L2 + L1  (x5)
2 L2 + L1 -> L0 + L1 + 2 S  (x1)
```

The first is production: 39 links were made. The second is a free link
absorbing a substrate. The third is bonding: a free link bonds to the end of a
chain, which becomes a chain link (`L1 -> L2`) while the newcomer becomes the
new end (`L0 -> L1`). The fourth is a disintegration inside a chain. A chain
link (`L2`) falls apart; its two partners, a chain link and a chain end, each
lose a bond and become an `L1` and an `L0`; and the dead link leaves two
substrates. Species on both sides of a reaction are the bond states of the
links involved, which is why one physical event can appear as several
different reactions.

One recording error affects this list: when two links of the *same* species
bond, the reaction is stored with one of each instead of two. The default run
contains `L1S -> L2S (x7)`, which is really `2 L1S -> 2 L2S`, and likewise
`L1 -> L2`, `L0 -> L1` and `L0S -> L1S`. The lattice itself is not affected.

### What the formal specification below means

The specification uses the book's symbols (`O`, `*`, `[O]`) for substrate,
catalyst and link, and `[O]O` for a link holding an absorbed substrate.
"Molecules per reaction: 1, 2, 3" counts the particles on the left of a
reaction: one for a decay or an emission, two for a bond or an absorption,
three for production. The probabilities in the parameter table are chances per particle
per time step in the lattice algorithm, not rate constants: no rate law has
been published for this model, so the reactions carry no rates.

## Using it

The default run is the book's figure 6.3 at the start: one catalyst in the
centre of a 30 × 30 world full of substrate. It runs 120 steps, enough to see
the elementary reactions but not a cell; with seed 1 the only closed chain at
the end is a three-link cluster. Everything about the run is in `net.extras`:

```python
a = net.extras["analysis"]
a["events"]        # {'absorption': 146, 'bond': 36, 'disintegration': 21, 'emission': 116,
                   #  'motion': 2851, 'production': 39, 'substrates_released': 41}
a["membrane"]      # closed_chains, chain_lengths, membranes (6+ links), clusters,
                   # enclosed_catalysts, first_enclosure_step, steps_enclosed, ruptures, repairs
a["permeability"]  # substrate / link / catalyst crossings of the cell boundary
a["per_step"]      # series: substrates, links, free_links, chain_links, closed_chains, ...
net.extras["space"]["grid"]    # the final lattice as strings, legend in space["legend"]
```

A catalyst counts as *enclosed* when some closed chain cuts it off from the
rest of the torus. `ruptures` counts the steps at which enclosure is lost,
`repairs` the steps at which it is regained after the first enclosure.

**Start from a cell.** `initial="cell"` places a ready-made membrane of twelve
links around the catalyst, enclosing a 3 × 3 interior (the cell of Von Kamp's
figure 2). With a 15 × 15 world and `disintegration_probability=0.001` this is
close to the set-up of McMullin and Varela (1997) described under Results:

```python
net = chemart.generate_network("autopoiesis-vmu", seed=1, initial="cell",
                               width=15, height=15, steps=2000,
                               disintegration_probability=0.001)
net.extras["analysis"]["membrane"]["steps_enclosed"]   # 111
net.extras["analysis"]["permeability"]
# {'substrate_crossings': 32, 'link_crossings': 0, 'catalyst_crossings': 0}
```

**Switch off the missing rule.** Add `bond_inhibition=False` to reproduce the
failure of the published algorithm; the grids above show what it does. The
original FORTRAN program inhibited a free link only when two or more chain
links were next to it (McMullin 1997); `chain_inhibit_count=2` gives that
threshold, while the default 1 is SCL's.

**Watch a cell form.** Spontaneous formation takes hundreds to thousands of
steps. On a 14 × 14 world with `steps=2000`, every one of seeds 0–5 closed a
chain around the catalyst at some point, first at steps 90 to 1,443. On the
default 30 × 30 world, seed 1 with `steps=3000` first encloses the catalyst at
step 374.

**The reaction list alone.** `mode="reactions"` returns the 16 reactions the
chemistry defines, without running the lattice: production, three bonding
reactions between unloaded links, and disintegration, absorption and emission
for each link species.

All of these are fast: about a second per thousand steps on the 30 × 30
lattice, less on smaller ones. Disintegration does not always conserve
matter: when a decaying link has no hole next to it, the second substrate
(and an absorbed one) is lost, so the total slowly drops.

## Results

### The 1974 paper: formation and repair

According to McMullin's (1997) analysis of the paper, Varela, Maturana and
Uribe (1974) showed two sequences of lattice snapshots from one run. The first, instants 0 to 6, is the book's figure 6.3: a lone
catalyst in substrate produces links, visible at t = 2, which bond into a
closed boundary around it by t = 6. The second, instants 44 to 47, shows the
boundary rupturing and being repaired. The caption of that figure, as quoted
by Von Kamp (2002), reads "Ongoing production of links re-establishes the unity
under changes of form and turnover of components." The paper also gave a
six-point test for whether an entity is autopoietic (book §6.1.5) and a
natural-language algorithm of the simulation.

McMullin's 2004 review places the paper historically. The word autopoiesis
was coined by mid-1971, and the simulation was written that year with Uribe's
help. The paper was rejected by several journals before *BioSystems* accepted
it in 1974. It was the first English-language publication of the theory, and
McMullin calls it the first application of agent-based computer modelling to
separating the living from the non-living.

**In Chemart.** Spontaneous formation is reproduced qualitatively: a test runs
a 14 × 14 world (seed 3) for 2,000 steps and checks that a membrane of six or
more links closes around the catalyst. The timing is not reproduced, and
should not be expected to be. McMullin (1997) found that the published
snapshots match the rediscovered program rather than the published
algorithm, and that program let a catalyst make up to eight links per step
from a 25-site neighbourhood. In Chemart a catalyst makes at most one link per
step, and closure takes hundreds of steps. Enclosures formed this way are
short-lived: 5 to 83 steps in total in the six 14 × 14 runs above.

### The flaws, and the missing rule (McMullin 1997; McMullin and Varela 1997)

McMullin's 1997 working paper compared the published algorithm with the
published snapshots and found them inconsistent. The number of links should
always equal the number of holes, yet the snapshots show about six links per
hole; six links appear in the first step although the algorithm allows one per
catalyst per step; and free links sit next to each other without bonding,
although the algorithm makes bonding practically certain. A rekeyed copy of a
FORTRAN IV listing of an early version (`EXP29.FOR`) explained all of these,
and contained the unpublished bonding restriction.

McMullin and Varela (1997) then tested the chemistry itself with SCL, a new
program written on the Swarm simulation system. Each experiment was five runs
from a 12-link membrane around one catalyst in a 15 × 15 torus, with
disintegration probability 0.001, which gives an expected 84 steps before the
first rupture.

- **Experiment 1, the published chemistry.** All five runs failed without a
  single repair. The interior filled with bonded links until no room was left
  for production. In run 1 only two open sites remained inside by step 110.
  In run 2 the only two links made inside bonded to each other at step 69 and
  boxed in the catalyst. Runs 3 to 5 were clogged by steps 282, 126 and 165. The
  authors add that two other independent re-implementations had failed the
  same way.
- **Experiment 2, with chain-based bond inhibition.** In three of five runs a
  self-repairing cell was established. In run 4 one membrane shape lasted from
  step 199 to 1,437 through 12 ruptures and repairs, and in run 1 the
  membrane lasted until it fragmented at step 1,746. In the other two runs the cell failed early.

Since the two experiments differed only in this rule, the authors concluded
that the phenomenon "relies critically on the presence of this interaction".
They stressed that this corrected the historical record without changing the
concept, and argued that computational models should be published with their
code.

**In Chemart.** The mechanism is reproduced; the long-lived cells are not.

- Keeping links free: a test runs a ready-made cell for 80 steps without decay
  on six seeds. With the rule at least five links stay free and at most two
  bonds form; without it at most three links stay free and at least four bonds
  form. The grids under *How it works* show the same thing.
- Repair: a test removes one membrane link from a cell after 20 steps and
  checks for re-closure within 60 steps: at least 6 of 10 seeds must repair
  with the rule, at most 3 without. Over seeds 0–19 I measured 14 of 20 with
  the rule and 4 of 20 without.
- Lifetime: in the McMullin–Varela set-up (15 × 15, decay 0.001, 3,000 steps,
  seeds 0–9) Chemart's cells do not last a thousand steps. The catalyst is
  enclosed for a median of 122 steps in total with the rule, and 376 without
  it; the first opening comes at a median of step 72 and 117. So the rule
  makes repair after a single break far more likely, but it does not make
  cells live longer in Chemart. Without the rule, the catalyst is often
  re-enclosed by chains that have grown from the interior, which the
  enclosure count cannot tell apart from a repair. Why the cells with the rule
  die sooner has not been established. Two differences from SCL are known: the
  catalyst can make links outside the membrane across its corners (see *How
  it works*), and the neighbourhoods themselves are Chemart's choice, since the
  sources do not state SCL's. Chemart's default decay rate, 0.01, is also
  different: it is the standard value of Von Kamp's extended system, ten times
  the 0.001 McMullin and Varela used.

### Permeability

The membrane passes substrate by absorption and emission and holds catalysts
and links, because bonded links cannot move. This is the book's "closed
membranes that are nevertheless permeable to substrate molecules". A Chemart
test runs a ready-made cell without decay for 120 steps on five seeds and
checks that the twelve-link membrane stays closed throughout, that substrate
crosses it, and that no link or catalyst does. Reproduced.

### Is it really autopoietic?

Later work doubted that the model's cells meet the definition they were made
to illustrate.

- **Individuation (McMullin 2000).** McMullin proposed that autopoiesis is
  roughly "collective autocatalysis plus spatial individuation", and tested
  whether two identical cells side by side keep separate identities. In
  unpublished SCL experiments, which he calls "at best, inconclusive", the
  cells were unstable even alone, and chain-based bond inhibition hindered the
  upkeep of two adjacent membranes, so that "adjacent agents tend positively to
  merge rather than to maintain their individuality". He concluded that SCL's
  cells do not pass his test for full autopoiesis. Chemart can place several
  catalysts (`n_catalysts`) but has no test or measurement of this.
- **The catalyst and the broken membrane (Von Kamp 2002).** Von Kamp, in
  McMullin's group, argued that original SCL has no autopoietic entities for
  two reasons: the catalyst is never produced by the cell, and a cell whose
  membrane breaks has no boundary, so the cell that re-forms is a new one
  rather than a repaired one.
- **Size (McMullin 2004).** The inhibition rule keeps the interior links free
  only while the cell is small, so this form of cell cannot easily grow or
  divide.

### What came after

McMullin's review traces the later models. Milan Zeleny re-implemented the
model in APL between 1975 and 1978 and reported growth, oscillations and
self-reproduction, but his code is lost and his changes may have broken the
locality of the rules. Breyer, Ackermann and McCaskill (1999) let bonded links
move and rearrange bonds, which removes the need for bond inhibition and lets
cells grow. McMullin and Groß (2001) disabled bonding between two free links
instead, obtaining longer-lived and growing cells at the cost of spontaneous
formation. Von Kamp's SCL-DIV (2002) added membrane growth, catalyst
production and fission, and obtained self-reproducing cells. Ono and Ikegami's
lattice chemistry, a coarser model that allows many particles per site, gave
cells that maintain and reproduce themselves; see
[the Ono–Ikegami protocell](ono-ikegami-protocell.md). The book lists more than a
dozen works that tried to reproduce and improve on the original experiments (§6.1.5)
and notes that the ideas were taken into the wet-lab construction of minimal
cells. None of these extensions are in this entry.

## Further reading

The Results above draw on the full texts of McMullin (1997) and McMullin and
Varela (1997). The implementation had only the abstract of the second and
none of the first, as the sources list says.

- McMullin, B. (2004). 30 years of computational autopoiesis: a review.
  *Artificial Life* 10(3), 277–296 (book ref [569]). The history of the model,
  its flaws and its descendants.
- Breyer, J., Ackermann, J. & McCaskill, J. S. (1999). Evolving
  reaction-diffusion ecosystems with self-assembling structures in thin films.
  *Artificial Life* 4(1), 25–40 (book ref [136]).
- McMullin, B. & Groß, D. (2001). Towards the implementation of evolving
  autopoietic artificial agents. In J. Kelemen & P. Sosík (eds.), *Advances
  in Artificial Life, ECAL 2001*, 440–443. Springer.
