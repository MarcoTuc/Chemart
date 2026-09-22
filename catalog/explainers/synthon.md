## Introduction

The Synthon artificial chemistry is a way of letting a computer invent a
reaction network by itself, starting from a handful of atoms and a short list
of *kinds* of reaction, while keeping the bookkeeping of real chemistry:
atoms, bonds, lone pairs, unpaired electrons and charge. It was built by Tom
Lenaerts and Hugues Bersini at the IRIDIA laboratory of the Université Libre
de Bruxelles. Their 2009 paper in *Artificial Life* is the one the book cites;
an earlier report
(2005) describes the same framework and applies it to the chemistry of
interstellar clouds, and that report is what Chemart is built from.

The central idea comes from mathematical chemistry: the synthon model of
Kvasnička and Koča, a formal model of the logical structure of organic
chemistry, and the Dugundji-Ugi model it extends. A molecule is a graph:
atoms are vertices, covalent bonds are edges, and each atom also carries its
lone pairs (electron pairs that take no part in a bond) and its unpaired
electrons. Now put the reactants of a reaction side by side and treat them as
one big graph, an *ensemble*. The products contain exactly the same atoms and
the same electrons; only the edges differ. So a reaction is an
*isomerisation* of the whole ensemble, in the same sense that two isomers are
molecules with the same atoms arranged differently. Take H₂⁺ + H₂ → H₃⁺ + H:
four hydrogen atoms and three electrons go in, four hydrogen atoms and three
electrons come out, and one hydrogen has moved from one partner to the other.
Conservation of mass and charge is not a check performed afterwards; a
reaction that broke it could not be written at all.

On top of this representation the authors put two network generators. One
applies every reaction class to every molecule and pair of molecules, again
and again, until no new molecule appears. The other runs a stochastic
simulation and only lets molecules that are actually present react, so the
network grows only where the kinetics goes. The authors' motivation, in the
2009 abstract, is that experiments meant to give insight into chemical
problems "should be expressed in a syntax that remains as close as possible to
real chemistry", and that studying the structure of reaction networks
"requires growing models that incorporate features realistic from a
biochemical perspective". In the 2005 report the application is interstellar
chemistry, where laboratory conditions are hard to reproduce, many reactions
are unknown and rate constants are "often merely estimates".

Banzhaf and Yamamoto mention Synthon in one paragraph of §18.3.2, among
artificial chemistries for modelling biochemical pathways. They set it next to
[ToyChem](toychem.md), the other graph-based chemistry in the catalog that
aims to be close to real chemistry. ToyChem computes the energy of each
molecule with an approximate quantum-chemical method and derives its rates
from it; Synthon takes its rates from a published table and instead
focuses, in the book's words, on "the dynamical aspects of the system, including stochastic chemical kinetics and the
coevolutionary dynamics of the system": the network and the concentrations
change together. It is a generator of concrete chemistry, with real elements
and published rate constants, which sets it apart from the abstract chemistries
of most of the catalog.

## How it works

### Molecules: atoms, electrons and a canonical name

Each atom records three things: its element, its number of lone pairs, and
its number of unpaired electrons (radicals). Bonds are covalent, with an
integer order (single, double, triple), or *ionic*. The ionic link is one of
two extensions Lenaerts and Bersini made to the original synthon model: a
molecule like H₂⁺ has only one electron shared between two protons, which an
ordinary two-electron bond cannot express, so they draw it as two disjoint
parts of one molecule, a proton and a hydrogen atom (H⁺ · H), and H₃⁺ as
H⁺ · H₂. Their other extension makes unpaired electrons explicit, so that
radicals and ions are different graphs. The free electron e⁻ is a species of
its own.

Charge is never stored. It is computed from the graph by the rule of the
Dugundji-Ugi model: an atom's formal charge is its
number of valence electrons, minus two per lone pair, minus one per unpaired
electron, minus the orders of its bonds. Oxygen has six valence electrons, so
an O with two lone pairs and one bond to H is OH⁺ (6 − 4 − 0 − 1 = +1), and
the same with one extra unpaired electron is the neutral radical OH.

Every molecule is written as a canonical code, so that the same graph always
gets the same name however it was built (the papers use Weininger's CANGEN
algorithm, the one behind unique SMILES strings). An atom is its element
followed by one `:` per lone pair and one `*` per unpaired electron; bonded
atoms are written next to each other, `=` and `#` mark double and triple
bonds, `.` marks an ionic link and parentheses mark branches. In Chemart's
H/O network:

| code | molecule |
|---|---|
| `H*` | hydrogen atom H (one unpaired electron) |
| `H` | proton H⁺ (no electrons at all) |
| `HH` | H₂ |
| `H.H*` | H₂⁺, a proton ionically linked to an H atom |
| `H.HH` | H₃⁺, a proton ionically linked to H₂ |
| `O::**` | oxygen atom (two lone pairs, two unpaired electrons) |
| `HO::*` | OH radical |
| `HO::` | OH⁺ |
| `HO::H` | water |
| `e-` | free electron |

Molecules with the same formula but a different split of electrons between
lone pairs and unpaired electrons are different species, called *electronic
isomers* below. Water `HO::H` and `HO:**H`, a water whose oxygen has one lone
pair and two unpaired electrons, are both H₂O. This is the Dugundji-Ugi notion
of isomerism, and it is why the network has more species than formulas.

### Reaction classes

A *reaction class* is itself a small graph: the atoms and electrons involved,
and for each edge whether it is formed (+1) or deleted (−1). The 2005 report
draws three such graphs (radiative association, charge transfer, dissociative
recombination) and lists eleven classes, G1 to G11, taken with their rate
constants from Duley and Williams's textbook *Interstellar Chemistry* (1984):

- **G1–G3, cosmic-ray ionisation**: a molecule loses an electron (`A → A⁺ + e⁻`),
  possibly falling apart (`AB → A + B⁺ + e⁻`);
- **G4, ion-molecule exchange**: `A⁺ + BC → AB⁺ + C`;
- **G5, charge rearrangement**: `A⁺ + B → A + B⁺`;
- **G6, G7, recombination**: an ion captures an electron, with or without
  breaking up (`AB⁺ + e⁻ → A + B`, `A⁺ + e⁻ → A`);
- **G8, neutral reactions**: `A + BC → AB + C`;
- **G9, G10, photodissociation**: a photon splits a molecule (G10 is H₂
  specifically, with its own rate);
- **G11, grain-surface reaction**: two H atoms meet on a dust grain and
  leave as H₂.

Cosmic rays, photons and grains carry no atoms, so they are not species: the
ionisation and photodissociation classes act on one molecule. A *reaction
object* is one concrete instance of a class, such as `H* → e- + H`, cosmic-ray
ionisation of a hydrogen atom, which is object 1 of the report's table 2.

### A worked example

Here is object 32 of table 2, as the default Chemart network produces it:

```
HH + H.H* -> H.HH + H*  [mass-action k=1e-09 classes=G4 units=cm^3 s^-1]
```

It reads H₂ + H₂⁺ → H₃⁺ + H, an ion-molecule exchange. Step by step:

1. The two reactants are merged into one ensemble: four hydrogen atoms, the
   H–H bond of `HH`, the ionic link of `H.H*`, and three electrons (two in the
   bond, one unpaired on `H*`).
2. The H–H bond of `HH` breaks *homolytically*, one electron to each side,
   leaving two hydrogen atoms with one unpaired electron each.
3. One of them bonds to the `H*` end of `H.H*`. Both partners have an unpaired
   electron, so the new bond is covalent: the proton now hangs, by its ionic
   link, on an H₂. That is `H.HH`, H₃⁺.
4. The other hydrogen atom leaves as `H*`.

Before: charges +1 and 0, electrons 1 and 2. After: charges +1 and 0,
electrons 2 and 1. The ensemble has the same atoms and electrons on both sides;
only edges changed. Its rate is table 1's constant for class G4.

### The two generators

The **deterministic network generator** (DNG) starts from the initial species,
applies every class to every molecule and every pair, keeps the products that
satisfy the size limits, and repeats with the new molecules until nothing new
appears; the result is the *closure* of the initial species. The limits
are those of the report: at most N atoms per molecule,
at most `es` lone pairs or unpaired electrons per molecule and `ep` per atom.
The report adds an *observational constraint*, a list of the molecules allowed
in the system, to keep its figure drawable.

The **MC-sampling network generator** (MCNG, MC for Monte Carlo), after
Faulon and Sault (2001),
simulates a finite population of molecules with Gillespie's stochastic
simulation algorithm (SSA). At each step it picks the next reaction event in
proportion to rate constant × number of possible reactant combinations, then
updates the counts. Reactions are only ever sought among species present at
that moment, so the network recorded is the part of the closure the kinetics
actually visits.

In Chemart the two generators are the two ways of running the chemistry:
`chemart.generate_network` runs the DNG and returns the closure, and
`chemart.evolve` runs the MCNG and returns its trajectory, one frame per
reaction event.

## Using it

The default call above is the report's figure 4 experiment: the DNG, the
eleven interstellar classes, starting from one hydrogen atom and one oxygen
atom, with the figure's thirteen molecules as the observational constraint.
It closes at 31 species and 161 reactions: the thirteen molecules plus their
electronic isomers. `net.extras` translates the codes:

```python
f = net.extras["formulas"]                   # {"H.HH": "H3^+", "HO::H": "H2O", ...}
net.extras["charges"]["HO::"]                # 1
net.extras["electrons"]["HO::H"]             # 8
net.extras["reaction_classes_used"][5]       # ['G4']  the class(es) of reactions[5]
net.extras["conservation"]                   # vectors: atoms of H, atoms of O, electrons, charge
```

When several classes produce the same stoichiometric reaction, their rate
constants are summed and `reaction_classes_used` lists them all. Rates are
mass-action constants in the table's units: s⁻¹ for one-molecule classes,
cm³ s⁻¹ for two-molecule ones. The initial state holds the densities of Duley
and Williams's setup, n(H) = 1000 cm⁻³ and n(O) = 0.44 cm⁻³.

**Letting kinetics prune the network (figures 5 and 6).** Run the MCNG with
`chemart.evolve`. The report does not publish its number of molecules or of
steps, so they are parameters of this run only (`molecules`, default 400;
`steps`, default 2000):

```python
traj = chemart.evolve("synthon", seed=3)
net = traj.network
net.summary()                  # synthon: 10 species, 7 reactions, status=observed
net.extras["final_state"]      # {'H*': 377, 'HH': 11, 'HO::*': 1}
```

Only 7 of the 161 reactions fired. Each reaction's `count` says how often:

```
982 HH + O::** -> H* + HO::*
980 H* + HO::* -> HH + O::**
 16 2 H* -> HH
  9 HH + HO::* -> H* + HO::H
  9 H* + HO::H -> HH + HO::*
  3 HH -> 2 H*
  1 HO::* -> H* + O::**
```

Out of 400 molecules the oxygen density gives oxygen a single atom (at least
one is always kept), and almost every event moves that atom between O and OH.
`extras["elapsed_time"]` is the simulated time in seconds, here about
3.2 × 10¹³ s. The run takes under two seconds.

The trajectory has a frame for the initial molecules and one after each of
the 2000 events: `frame.t` is the simulated time in seconds (the trajectory's
`clock` is `"s"`), `frame.state` the molecule counts and `frame.fired` the
reaction that fired. `every=k` keeps one frame in k, which turns the run
into a short time series of the counts, the kind of curve figure 6 plots:

```python
traj = chemart.evolve("synthon", seed=3, every=500)
for f in traj.frames:
    print(f"{f.t:9.3g}", f.state)
```

```
        0 {'H*': 400.0, 'O::**': 1.0}
 1.49e+13 {'H*': 397.0, 'HH': 1.0, 'HO::*': 1.0}
 2.42e+13 {'H*': 394.0, 'HH': 3.0, 'O::**': 1.0}
 2.93e+13 {'H*': 386.0, 'HH': 7.0, 'O::**': 1.0}
 3.22e+13 {'H*': 377.0, 'HH': 11.0, 'HO::*': 1.0}
```

**The combinatorial explosion.** `allowed=[]` removes the observational
constraint. Keep `max_atoms` small:

| call | status | species | reactions | time |
|---|---|---|---|---|
| `max_atoms=3` | complete | 22 | 95 | 1.6 s |
| `allowed=[], max_atoms=3, max_species=400` | complete | 110 | 3,060 | 2.8 s |
| `allowed=[], max_atoms=4, max_species=400` | truncated | 400 | 23,796 | 61 s |

Without the constraint, three-atom molecules already bring in O₃, HO₂, H⁻,
O⁻ and their ions (24 formulas). At four atoms the closure hits the species
budget.

**The figure 3 classes.** `templates="fig3"` uses only radiative association,
charge transfer and dissociative recombination. From H and O it builds four
reactions (`2 H* -> HH`, `H* + O::** -> HO::*`, `2 O::** -> O::*O::*`,
`H* + HO::* -> HO::H`). The report gives no rates for these graphs, so the
reactions have no rate and `chemart.evolve` refuses this set.

Other elements can be added with `initial` (H, He, C, N, O, F, Ne, S, Cl, Ar),
but then `allowed` must list their molecules, or be emptied.

## Results

**The Artificial Life paper (2009).** The paper the book cites could not be
obtained for Chemart, only its abstract. The abstract describes "a
coevolutionary model" joining "the logical structure of constitutional
chemistry and its kinetics" with "the topological evolution of the chemical
reaction network", and two illustrative examples showing "that one needs to be
careful in making general claims concerning the structure of chemical reaction
networks". The book adds that its experiments show "the partition of the
molecule sets into two categories of chemical reaction networks", a phenomenon
that could help explain the origin of homochirality: the fact that living
things use only one of the two mirror-image forms of many molecules. Chemart
does not reproduce this result, or anything else specific to the 2009 paper,
because its content is unknown here.

Everything below comes from the 2005 report.

**An H/O network grown from two atoms.** From H and O alone, with the eleven
classes and the constraints N = 6, es = 10, ep = 7, the DNG grew the network of
the report's figure 4, restricted to thirteen molecules: H, H⁺, H₂, H₂⁺, H₃⁺,
O, OH⁺, OH, H₂O, O₂, H₂O⁺, H₃O⁺ and e⁻. The authors note that even with
this restriction "the reaction network can become very large", and that the
network alone does not say which reactions matter. Chemart's tests check that
exactly these thirteen formulas appear, that each of the eleven classes
contributes, and that N = 6 gives the same species as the default N = 4. The
figure's reaction objects cannot be counted from the available text, so the
161 reactions of Chemart's network are not compared with it.

**Reactions conserve everything.** The report states that its isomeric
reactions "conserve mass, charge and radicals". Chemart's tests check, for
every reaction, that the atoms of each element, the electrons and the charge
balance, and that each conservation vector annihilates the stoichiometric
matrix. They also check that the canonical code does not depend on the order
in which a molecule is built, and that electronic isomers are separate species.

**Ionic bonds and H₃⁺.** Table 2 lists sample reaction objects. The tests
find object 1 (H + c.r. → H⁺ + e⁻, c.r. standing for a cosmic ray) and object 32 (H₂⁺ + H₂ → H₃⁺ + H) in the
network, with H₂⁺ and H₃⁺ written as `H.H*` and `H.HH`, the report's own
H⁺ · H and H⁺ · H₂. Object 8, H⁺ + H → H + H⁺, is an electron hopping from
the atom to the proton. It changes nothing in the counts of species, so
Chemart leaves it out of the network; the tests check that class G5 does
produce it at the level of the graphs. Object 56 is garbled in the text that
was available and is not tested.

**Kinetics prunes the network.** With MCNG, starting from n(H) = 1000 cm⁻³
and n(O) = 0.44 cm⁻³, the report found "a drastic reduction in the size of the
reaction network keeping only the most relevant reactions" (figure 5), and
simulated concentrations in which "H and H2 are the dominant species"
(figure 6). For that run the authors also removed reactions like object 8, to
stay close to Duley and Williams's setup. Chemart's test runs 2000 events and
checks that the observed network has fewer than a quarter of the closure's
reactions, that H and H₂ are the two most abundant species at the end, and
that the trajectory has one frame per event.

**The combinatorial explosion.** The report opens with the drawback of formal
network generators: the number of generated molecules grows exponentially
with the number of atoms, which is why the observational constraint and the
kinetic generator are needed. Chemart's (slow) test reproduces this at three
atoms: without the constraint the closure has more than four times the
species and more than ten times the reactions, including O₃, HO₂ and H⁻.

**What Chemart adds.** The report does not say when an ionic link may form,
whether bonds break into radicals or into ions, or which molecules are
physically possible. Chemart's answers, and the guards it added because the
unguarded closure filled with species like O⁶⁺, are listed in the
implementation decisions above.

## Further reading

- Duley, W. W. & Williams, D. A. (1984). *Interstellar Chemistry*. Academic
  Press. The source of the reaction classes and rate constants of table 1.
- Faulon, J.-L. & Sault, A. G. (2001). Stochastic generator of chemical
  structure. 3. Reaction network generator. *J. Chem. Inf. Comput. Sci.* 41,
  894–908. The DNG and MCNG generators.
- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical
  reactions. *J. Phys. Chem.* 81(25), 2340–2361.
