## Introduction

This chemistry is a proposal by Alan Dorin and Kevin Korb (Monash University,
2007) for building a virtual ecosystem from the bottom up, starting at the level
of atoms. In most artificial-life ecosystems the organisms are agents with
senses and behaviours, and the environment they live in is a different kind of
thing: a space, perhaps with "abiotic furniture" such as grains of sugar that
the agents pick up. Dorin and Korb point out that this split leaves out what
ecologists consider the core of an ecosystem, the *biogeochemical cycles* by
which elements pass from the environment into organisms and back again.
Simulated organisms, they write, "often die, but they seldom rot to return
materials to their environment". Their answer is to make organisms and
environment out of the same stuff: square virtual atoms on a grid that bond,
unbond, and store or release energy as they do.

The picture is a closed world of four kinds of atom, `A`, `B`, `C` and `O`,
plus a few catalyst atoms. An `A`–`B` bond plays the role of sugar: it takes a
lot of energy to make and gives that energy back when it breaks. `C`–`C` bonds
play the role of biomass. A chlorophyll-like catalyst uses sunlight to turn
`A`–`O` and `B`–`O` into sugar; an enzyme breaks sugar and releases its
energy; the released energy lets carbon chains grow; other enzymes break carbon
chains and inorganic bonds, returning the atoms to the pool. An "organism" is
nothing more than a large molecule that happens to carry catalysts in the right
places: the authors hand-design a photosynthetic plant-like body and a
decomposer body. Whether a molecule counts as organic, as part of an organism or
as environment depends only on its bonds and surroundings.

It is best described as a simulation model offered as a feasibility argument.
The paper specifies the chemistry and the organism bodies, argues that producers
(autotrophs, which make their own food), consumers and decomposers (heterotrophs,
which live on what others made) all fit in it, and says plainly that a
simulation with all organism types together "remains to be performed". It
reports no runs and no numbers. Chemart's version is therefore a
reconstruction: the rules are the paper's, the numbers are Chemart's choices.

Banzhaf and Yamamoto describe it at the end of book §8.2.3, "Ecological
Modeling with Artificial Chemistries", as "a more general approach to ecosystem
modeling". Its neighbours there model ecology at a coarser level:
[EVOLVE](evolve-series.md) also conserves matter and energy but gives
organisms genomes whose genes are looked up in a function table,
[Ecolab](ecolab.md) evolves the coefficients of Lotka–Volterra population
equations, and [Urdar](urdar.md) has cellular-automaton organisms feeding on
bitstrings. Dorin and Korb's system is the only one in that group built from
atoms and bonds on a grid, and the only one with no evolution: its footnote
says the model "is not concerned with reproduction, evolution or the
self-assembly of organisms". Mechanically, the closest relative in the catalog
is [Squirm3](squirm3.md), another grid of atoms that bond to their neighbours
and move as whole molecules, whose aim is self-replication rather than energy flow.

## How it works

### Atoms, bonds and the energy they carry

Each atom fills one square of a grid (in Chemart, a torus: the edges wrap
around). Two atoms can bond when they share an edge. How many bonds an atom
can hold, its *valence*, comes from a small electron-shell model in the
paper's appendix: `A`, `B` and `O` take one bond each, `C` takes four, and
every catalyst takes two. So `A`–`O` is a finished molecule, while carbon can
form chains and branches.

Only five kinds of bond react, and the paper's table 1 gives each of them a
probability of forming and of breaking per time step, with and without a
catalyst, and a bond energy:

| bond | stands for | forms | breaks | catalysed | energy |
|---|---|---|---|---|---|
| `A-B` | sugar | low | low | chlorophyll makes it, sugar enzyme breaks it | − high |
| `C-C` | biomass | moderate | low | organic enzyme breaks it | − low |
| `A-O` | inorganic | high | low | chlorophyll or inorganic enzyme breaks it | + low |
| `B-O` | inorganic | high | low | chlorophyll or inorganic enzyme breaks it | + low |
| `C-O` | carbon–oxygen (the decomposer's spacer) | low | moderate | organic enzyme makes it | + low |

The sign convention is the table's own: a positive bond gives out its energy
when it forms and must be paid to break; a negative bond is the opposite, it
costs energy to form and gives the energy back when it breaks. The negative
bonds are the ones that *store* energy, which is why `A-B` is sugar and `C-C`
biomass. `A` and `B` bind readily to `O` ("high" probability, and it releases
energy), so free sugar ingredients are scarce unless something pries them off
oxygen. The paper only grades the entries as low, moderate and high; Chemart
turns them into numbers (see *Using it*).

### Catalysts

There are four catalyst atoms, named in Chemart after their jobs: `K`
(chlorophyll), `EAB` (the sugar-breaking enzyme), `ECC` (the organic
decomposer) and `EO` (the inorganic decomposer). A catalyst acts when it
touches either of the two atoms of a bond; it then raises that bond's
probability to "high" and is not changed by the reaction.

### Energy is local and cannot be saved

Energy has to come from somewhere nearby. When a bond releases energy, that
energy goes into a pool shared by every atom that touches the reaction site,
directly or through a chain of touching neighbours (in Chemart, the *contact
cluster*). A bond that costs energy can form or break only if its cluster's
pool holds enough. Whatever is not spent in the same time step is lost for
good, so the only way to keep energy is to lock it into a negative bond.
Sunlight is the outside source: every atom receives an amount that rises and
falls with a sine wave over time, but only chlorophyll-catalysed reactions may
use it.

### One time step

1. **Movement.** Each molecule, with some probability, tries to slide one
   square in a random direction. Bonded atoms move together, keeping their
   shape, and a move that would bump into anything is cancelled.
2. **Reactions.** Every pair of touching atoms whose bond is in the table is
   considered once, in random order. A pair that is bonded may break; a pair
   that is not, and has free valence on both sides, may bond. The reaction
   happens with the table's probability (the catalysed one if a catalyst is
   touching) and, if it costs energy, only when the cluster can pay.

### The paper's three reactions, as bond events

The paper writes three reactions:

```
AO + BO  -(chlorophyll & sunlight)->  AB + 2 O     photosynthesis
O + AB   -(enzyme)->                  A + BO + energy   respiration
C + C    -(energy)->                  C2           biosynthesis
```

In the model these are not extra rules but sequences of single bond events
from the table. Photosynthesis is chlorophyll breaking `A-O`, breaking `B-O`,
then joining the freed `A` and `B`, all three paid for by sunlight (or, when
there is not enough, by energy released nearby).
Respiration is the enzyme breaking `A-B`, which releases the stored energy,
and `B` then bonding to a free `O`. Biosynthesis is a `C-C` bond forming with
energy some other reaction has just released next to it.

Here are the corresponding lines from the default run (`x13` means the
reaction fired 13 times; the number after each is the energy it released in
Chemart's units, negative for a cost):

```
AO + K -> A + O + K          (x47)   -1   chlorophyll pries A off O
BO + K -> B + O + K          (x52)   -1
A + B + K -> AB + K          (x13)   -8   chlorophyll makes sugar
AB + EAB -> A + B + EAB      (x7)    +8   the enzyme releases the stored energy
2 C -> C2                    (x3)    -1   biosynthesis, paid by energy released nearby
C2 + C4OECC_5c53 -> 2 C + C4OECC_5c53  (x1)  +1   the decomposer breaks biomass
```

`C4OECC_5c53` is the seeded decomposer body: four carbons, an oxygen and an
`ECC` enzyme (the suffix is a hash that tells molecules of the same formula
apart). It appears on both sides because it acts as the catalyst.

### The organism bodies

The paper draws two bodies, which Chemart can place on the grid:

- **The photosynthetic autotroph** (figure 5): a closed square wall of carbon
  around a cavity (a *vacuole*). Chlorophyll is fixed to one inner wall and the
  sugar enzyme to the opposite one. `A-O` and `B-O` in the vacuole are turned
  into sugar at the chlorophyll wall and stay trapped; when a sugar molecule
  drifts to the enzyme wall it is respired, and the energy spreads through the
  wall where it can pay for new `C-C` bonds. In Chemart the wall is 16 carbons
  and the vacuole starts with one `A-O` and one `B-O`.
- **The decomposer** (figure 6): a short carbon chain carrying the organic
  enzyme on an oxygen spacer, so that the enzyme can break other structures'
  `C-C` bonds without touching its own.

The bonds that hold these bodies together are *anchors* in Chemart: they never
react, because the paper gives them no table entry. The bodies can still gain
and lose atoms through ordinary table bonds.

## Using it

The default run is not a published experiment, since there is none. It is an
18×18 grid seeded with both organism bodies, 116 free atoms and 16 free
catalysts, run for 150 steps (about 2 seconds). The network it returns is the
set of reactions that actually fired, each with its count. Species names are
the atom formula (`AO`, `C2`) for small molecules, and formula plus a hash for
larger ones: `C16KEAB_d33a` is the autotroph body at the start, and
`C19KEAB_7460` is the same body at the end of the run after it had bonded
three more carbons.

The trophic summary counts each step of the cycle over the whole run. To
follow the world step by step, `chemart.evolve` returns a trajectory with one
frame per time step, the first being the seeded grid. Each frame holds the
molecules present (`state`), the bond events of that step (`fired`) and four
`observables`: the `sugar_bonds` (`A-B`), `biomass_bonds` (`C-C`) and
`inorganic_bonds` (`A-O` plus `B-O`) present, and the `free_atoms`, atoms
with no bond. Its network is the one `generate_network` returns:

```python
traj = chemart.evolve("dorin-korb-ecosystem", seed=1)
net = traj.network
a = net.extras["analysis"]
a["trophic"]["sugar_made"], a["trophic"]["sugar_respired"]         # (13, 7)
a["trophic"]["biomass_built"], a["trophic"]["biomass_decomposed"]  # (12, 5)
a["trophic"]["inorganic_split"]                                    # 146
max(traj.series("sugar_bonds"))  # 5: A-B bonds present at once, at most
traj.frames[-1].observables
# {'sugar_bonds': 0, 'biomass_bonds': 23, 'inorganic_bonds': 33, 'free_atoms': 63}
net.extras["energies"]["ledger"]
# {'consumed': 271, 'dissipated': 246, 'light_incident': 62014, 'light_lost': 61797,
#  'light_spent': 217, 'released': 300, 'balanced': True}
```

The ledger is the energy account: energy `released` by bonds plus sunlight
`light_spent` equals energy `consumed` by bonds plus energy `dissipated`
unused at the end of a step. Nearly all the incident light is lost; light is
not the limiting factor at the default settings. `a["event_counts"]` breaks
every bond event down by bond, action and catalyst, and `extras["events"]`
gives the same for each reaction in the network. Note that in this run all 13
sugar molecules were made by free-floating chlorophyll atoms in the soup, not by
the autotroph's own chlorophyll; the body grew on energy released by reactions
in its contact cluster.

The numbers behind "low", "moderate" and "high" are parameters (`p_low`,
`p_moderate`, `p_high`, `energy_low`, `energy_high`), as are the light
(`light_amplitude`, `light_period`) and the starting inventory (`atoms`,
`catalysts`, `structures`); see the table below.

**Switching the light off.** With `light_amplitude=0.0` the chlorophyll has
nothing to work with. Over seeds 0–9 at default settings, no sugar bond ever
formed, against 5 to 15 made per run with the light on. Biosynthesis still
happens a little (0 to 4 `C-C` bonds per run), paid by energy that other bonds
release nearby.

**The autotroph on its own.** Remove everything but the figure 5 body:

```python
net = chemart.generate_network("dorin-korb-ecosystem", seed=5,
                               atoms={}, catalysts={}, structures="photoautotroph")
```
```
AO + C16KEAB_d33a -> A + O + C16KEAB_d33a  (x19)
BO + C16KEAB_d33a -> B + O + C16KEAB_d33a  (x19)
O + B -> BO  (x18)
A + O -> AO  (x19)
A + B + C16KEAB_d33a -> AB + C16KEAB_d33a  (x1)
AB + C16KEAB_d33a -> A + B + C16KEAB_d33a  (x1)
```

Chlorophyll on the wall split the trapped `A-O` and `B-O` 19 times each, and
each time the freed atoms bonded back to oxygen before `A` and `B` could meet.
Once they did meet at the chlorophyll, sugar formed; later it reached the
enzyme wall and was respired. The 8 units it released were lost at the end of
the step: with no free carbon there is nothing to build. This outcome is not typical:
over seeds 0–19, sugar was made in only 4 runs (once each with seeds 5, 11
and 17, three times with seed 13).

**The soup without bodies.** `structures="none"` gives a lively inorganic soup
(6 to 12 sugar bonds made per run over seeds 0–4) but builds almost no biomass
(1 to 3 `C-C` bonds), because a free carbon can only bond where a neighbouring
reaction has just released energy.

**Bigger worlds.** Scale the grid and the inventory together. A 32×32 grid for
400 steps with about three times the atoms takes about 2 seconds; a 64×64 grid
for 1,000 steps with about 1,650 atoms took 22 seconds and made 551 sugar bonds and
151 `C-C` bonds. Much larger grids, as the paper envisaged, will be slow.

## Results

**What the paper claims.** Dorin and Korb (2007) present no simulation
results. Their contribution is the design: a bond table, an energy rule and
hand-built organism bodies from which, they argue, a complete ecosystem can
be assembled, with organisms that "naturally fall into trophic levels,
generate energy from chemical bonds and transform material elements in the
process". Besides the photosynthetic autotroph and the decomposer, they
describe a *chemosynthetic autotroph*, which gets its energy without sunlight
by splitting `A-O` and `B-O` with a surface catalyst and letting the atoms
rejoin against its body, and a heterotroph that lives on sugar made by others.
Their conclusion is that "at least in principle there is nothing preventing"
a run with all organism types from succeeding. They hoped to adjust the rules
later so that structures could self-assemble, replicate and evolve. The book
repeats the design and the authors' view that it "does more justice to
ecosystem mechanism" than other artificial-life ecosystem models; it cites no
later work built on it.

What Chemart can show is whether the rules as written do what the paper says
they should. Its tests check the following.

**A trophic cycle emerges from the bookkeeping.** In the default run,
chlorophyll fixes sunlight into sugar, the enzyme breaks sugar, carbon chains
grow, the decomposer breaks them, and inorganic bonds are split and reformed.
The tests require each of these to happen (at least 5 sugar bonds made, 2
respired, 5 `C-C` bonds built, 1 decomposed, 20 inorganic bonds split), that
each catalyst performs exactly the reactions table 1 allows it and no others,
and that removing every organic-decomposer atom stops all catalysed `C-C`
breaking. A slow test repeats the checks on a 32×32 grid over 400 steps.
Whether this counts as the "emergence" of trophic levels is a matter of
reading: the catalysts' roles are written into the table, and what the run
shows is that the roles connect into a cycle of matter and energy.

**Sugar needs an energy source.** With the light off, the test run (seed 1)
makes no sugar at all, and none formed in the ten default-density runs above.
This is not a hard rule of the model, though. Any reaction may be paid from
the local pool, so in a denser soup (a 20×20 grid with 320 building-block
atoms, seed 0) energy released by other bonds paid for five sugar bonds in the
dark. Sunlight is the main source, not the only one.

**The autotroph works alone, sometimes.** The tests reproduce the figure 5
story with seed 5: the body photosynthesises the atoms in its vacuole and then
respires the sugar, and its wall stays intact. As shown above, this happens in
a minority of seeds, because photosynthesis is three separate bond events and
the freed `A` and `B` usually rebond to oxygen first. The paper gives no
probabilities, so it neither supports nor contradicts that ratio. Biosynthesis
does not happen in this setting.

**Exact conservation.** The model's premise is a closed ledger. The tests
check that every reaction conserves every atom type, that the final inventory
equals the seeded one, that each reaction makes or breaks exactly one bond
with the energy table 1 assigns, and that the energy account balances to the
unit (Chemart uses integer energies for this reason).

**Not reproduced.** Chemart does not build the chemosynthetic autotroph or a
sugar-eating heterotroph as bodies, although free chlorophyll, inorganic enzyme
and sugar enzyme atoms play those roles in the soup. It does not attempt what
the paper leaves for future work: all organisms together over long runs,
self-assembly, reproduction or evolution. The seeded bodies do not maintain
themselves against decay either, since their structural bonds are anchors that
never break; the paper's point that an organic structure "must produce
sufficient energy to sustain itself against natural decay" applies in Chemart
only to the bonds a body adds during the run. The paper's bond probabilities
and energies are qualitative, so any number Chemart produces depends on its
chosen values.
