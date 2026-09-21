## Introduction

ToyChem is an artificial chemistry that tries to look and behave like the
organic chemistry taught in an introductory course. Its molecules are
structural formulas: graphs whose vertices are atoms (hydrogen, carbon,
nitrogen, oxygen) and whose edges are single, double or triple bonds. Its
reactions are reaction mechanisms written as graph rewriting rules: "wherever
you find this pattern of atoms and bonds, rewire it like this". And, unlike
almost every other chemistry in the catalog, every molecule has an energy that
is computed from its structural formula alone, by a stripped-down version of a
quantum-chemical method. Nobody assigns energies or rate constants by hand;
they follow from the graphs.

Gil Benkö, Christoph Flamm and Peter Stadler, in Vienna and Leipzig,
described it in a series of papers from 2003 to 2005. Their question was
about reaction networks, not about single molecules. Some real networks, such
as a cell's metabolism, look "small-world"
(most species are a few reactions apart, and neighbours of a species tend to be
neighbours of each other), but the few networks known in detail are sampled
with a bias, because chemists only record the reactions they care about. The
2003 paper therefore asks for "a computational model for generating them": a
chemistry realistic enough that its networks mean something, cheap enough to
build large ones, and one whose chemistry can itself be varied. Their answer
keeps what they call the "look-and-feel" of chemistry: atoms are conserved,
reactions follow known mechanisms, and energy decides which reactions matter.

Concretely, you give ToyChem a few starting molecules and mechanisms, and it
builds every molecule and reaction they lead to: a network in which each
species has a structure and an energy and each reaction has a barrier.

Banzhaf and Yamamoto give it one paragraph and one figure in their chapter on
bio-inspired chemistries, as the first of the network models (book §11.3.1).
Its catalog neighbours make the contrast clear.
[Conservative random reaction networks](bigan-conservative-crn.md) (§11.3.2)
also care about thermodynamics, but draw species and energies at random, with
no structure. [Suzuki's network chemistry](nac.md) (§11.3.3) uses graphs too,
but there one graph is both the space and the molecules. The
[synthon chemistry](synthon.md) is the closest relative: it also treats
molecules as atom graphs that conserve atoms and electrons; ToyChem adds an
energy computed from each graph.
ToyChem's rewriting scheme was later reused by the
[RNA-folding ribozyme chemistry](rna-folding-ac.md) of book §18.1.1, in which
folded RNAs act as catalysts for ToyChem-style reactions.

## How it works

### Molecules and their energy

A molecule is written as a SMILES string, the standard one-line notation for a
structural formula: `C=O` is formaldehyde, `OCC=O` is glycolaldehyde
(HO-CH2-CHO), `C=CC=C` is butadiene. Hydrogens are implied. Two molecules are
the same species if their graphs are the same up to renumbering, which is
decided by comparing their *canonical* SMILES, a unique spelling of each
graph. That is how the 2003 paper tests for identity too.

To get an energy, ToyChem first turns the structural formula into an **orbital
graph**. Each atom gets its outer orbitals, the places its bonding electrons
can sit: one `1s` orbital for hydrogen, and four for carbon, nitrogen and
oxygen. Which four depends on how many neighbours and lone pairs the atom has,
following the VSEPR rules of textbook chemistry: a carbon with four neighbours
has four `sp3` hybrid orbitals, a carbon in a double bond has three `sp2`
hybrids and one leftover `p` orbital. The orbitals become the vertices of the
new graph, and an edge joins two orbitals on neighbouring atoms that overlap: a
strong *sigma* overlap for the two hybrids pointing along a bond, a *pi*
overlap between the `p` orbitals of a double bond, and weaker ones: a
*semi-direct* overlap when only one of the two orbitals points along the bond,
*hyperconjugation* between a `p` orbital and an `sp3` orbital next to it, and
the bent *banana bonds* of strained three- and four-membered rings. The 2003
paper credits this kind of graph to the chemist O. E. Polansky, hence "the
orbital graph of Polansky" in the specification below.

Every overlap gets a number taken from a table (Table 1 of the 2003 paper),
rather than being computed from atomic positions. From these numbers the model
builds the matrices of Extended Hückel Theory, a simple molecular-orbital
method, and solves them. The result is a list of molecular orbital energies.
Filling the lowest ones with two electrons each gives the molecule's energy.
Subtracting the energy of the separated atoms gives its **total atomization
energy** (TAE), the energy released when the molecule is assembled from free
atoms. It is negative for every bound molecule, and more negative means more
stable. Two orbitals matter for reactions: the **HOMO**, the highest orbital
that holds electrons, and the **LUMO**, the lowest empty one.

Because nothing is embedded in space, the model cannot tell cis from trans
isomers or a molecule from its mirror image, and it handles only neutral
molecules with paired electrons.

### Reactions as graph rewriting

A rule has a *left graph* (the bonds as they are before), a *right graph* (the
bonds after) and a *context* (atoms that must be present but keep their
bonds). Every rule must obey two conservation laws: it may not create, destroy
or change atoms, and it must keep the total bond order, which stands for the
number of valence electrons. A reaction between two molecules is written as two
half rules, one per partner, plus a join rule that makes the new bond between
them.

Chemart implements three rules:

- **Diels-Alder**, copied from the rule printed in the 2003 paper: a chain
  `C=C-C=C` (a *diene*) and a `C=C` (a *dienophile*) close into a six-membered
  ring. Butadiene plus ethene gives cyclohexene.
- **Keto-enol tautomerism**: a hydrogen moves so that `H-C-C=O` becomes
  `C=C-O-H` (an *enol*), or back.
- **Aldol condensation**: the end carbon of an enol attacks the carbon of a
  `C=O` group in another molecule, joining the two by a new carbon-carbon bond.

### How reactive is a reaction?

For each pair of molecules a rule applies to, ToyChem needs to say how easily
they react. The papers use Frontier Molecular Orbital theory: a reaction goes
more easily when the HOMO of one partner lies close in energy to the LUMO of
the other. The papers write the reactivity as `ξ / (E_LUMO − E_HOMO)`, with a
constant `ξ` per mechanism, and say it can be turned into rate constants with
the Arrhenius law. Chemart takes the gap `E_LUMO − E_HOMO` itself as the
**activation energy** `Ea`, the energy barrier in the Arrhenius rate
`A exp(−Ea/RT)`. For a reaction of one molecule the gap is that molecule's own
HOMO-LUMO gap.

### Growing a network

Starting from the seed molecules, the generator tries every rule on every
molecule and every pair of molecules, adds any new products, and repeats with
the new molecules until nothing new appears or a species budget is reached.
The 2003 paper calls this "orderly generation".

### A worked example

The default network is the formose reaction, which turns formaldehyde into
sugars. Its first reactions are printed further down. Here is what the first
steps mean, with numbers from the default run.

The aldol rule needs an enol, and formaldehyde (`C=O`) cannot become one: it
has no carbon neighbour to take a hydrogen from. So glycolaldehyde
tautomerises first:

```
O=CCO -> OC=CO        reaction energy +25.9 kcal/mol, Ea 131.1 kcal/mol
```

Its product, the *ene-diol* HO-CH=CH-OH, is less stable (TAE −3428.8 against
−3454.7 kcal/mol), so the step is uphill. Its barrier is glycolaldehyde's own
frontier gap: LUMO −9.94 eV minus HOMO −15.63 eV is 5.69 eV, or 131.1 kcal/mol.
Now the ene-diol can attack formaldehyde:

```
C=O + OC=CO -> O=CC(O)CO        reaction energy −23.2 kcal/mol, Ea 107.0 kcal/mol
```

The product is glyceraldehyde, a three-carbon sugar, and the step is downhill.
Its barrier is the smaller of the two cross gaps: formaldehyde's LUMO
(−9.13 eV) minus the ene-diol's HOMO (−13.77 eV), 4.64 eV or 107.0 kcal/mol.
Both molecules keep their atoms: one carbon, two hydrogens and one oxygen plus
two carbons, four hydrogens and two oxygens make C3H6O3. Glyceraldehyde can
enolise in turn, attack or be attacked again, and so on. Following the 2003
paper's figure 8, no product with more than four carbons is formed, so the
network closes with 11 species and 16 reactions.

## Using it

The default call builds the formose network of the 2003 paper's figure 8 from
formaldehyde and glycolaldehyde, with keto-enol and aldol rules. It is
complete, not truncated. Species names are canonical SMILES; `O=C(CO)CO` is
dihydroxyacetone and `O=CC(O)C(O)CO` a four-carbon sugar (a tetrose, such as
erythrose; the model has no stereochemistry, so it cannot tell which).

Everything else is in `net.extras`:

```python
e = net.extras["energies"]
e["total_atomization_energy"]["C=O"]    # -1722.5   kcal/mol, one per species
e["homo_ev"]["C=O"], e["lumo_ev"]["C=O"] # about (-15.64, -9.13)  frontier levels in eV
e["reaction_energy"][2]                 # -23.2     one per reaction, in order
net.extras["reaction_rules"][2]         # ['aldol']
net.extras["conservation"]              # atoms of type C, H and O: one vector each
net.extras["analysis"]["substrate_graph"]  # network statistics, see Results
```

Each reaction carries its rate law, `A exp(−Ea/RT)` with `Ea` in kcal/mol, so
the network can go to an ODE solver or a stochastic simulator. `A` is 1 and no
initial concentrations are set, because the papers give neither.

**The Diels-Alder network.** `network="diels-alder"` starts from the five
molecules of the 2003 paper's figure 7 (cyclobutadiene, ethenol, phthalic
anhydride, methylbutadiene, cyclohexa-1,3-diene). Every Diels-Alder product can
react again, so this network never closes, and `max_species` decides where it
stops. With 40 species, the size of the published network, it takes about
5 seconds:

```python
net = chemart.generate_network("toychem", seed=1, network="diels-alder", max_species=40)
len(net.species), len(net.reactions), net.status    # (40, 41, 'truncated')
```

The default budget of 60 species takes about 10 seconds. Larger budgets
grow the time quickly.

**Your own chemistry.** `network="custom"` takes any seed molecules over H, C,
N and O and any of the three rules. Butadiene and ethene with the Diels-Alder
rule:

```python
net = chemart.generate_network("toychem", seed=1, network="custom",
                               seed_molecules=["C=CC=C", "C=C"], rules=["diels-alder"],
                               max_species=8)
```

```
2 C=CC=C -> C=CC1CC=CCC1               reaction energy -92.4, Ea 175.1
C=CC=C + C=C -> C1=CCCCC1              reaction energy -112.1, Ea 175.1
C=CC=C + C=CC1CC=CCC1 -> C1=CCC(C2CC=CCC2)CC1    -97.4, Ea 120.9
...
```

The first two are the classic adducts, 4-vinylcyclohexene and cyclohexene.
After that, butadiene keeps adding to any leftover double bond, so the network
grows without end, and the budget stops it (`status` is `truncated`).

**Other knobs.** `barrier_cutoff` drops every reaction whose barrier exceeds
it. `rewrite_mode="random"` keeps one product channel per reacting pair, drawn
with the seed; `"priority"` keeps the most exothermic one (10 species and 12
reactions on the default network). `kappa` and `semi_direct_scale` change the
energy model.

## Results

**Networks that look like chemistry.** The 2003 paper's main demonstration is
that the model produces sensible chemistry from two textbook systems: a
network of repeated Diels-Alder reactions (its figure 7, built in three rounds
of orderly generation) and the formose network (figure 8), in which
formaldehyde and glycolaldehyde condense into sugars. Chemart's tests check
that the Diels-Alder rule gives exactly cyclohexene from butadiene and ethene
and 4-vinylcyclohexene from two butadienes, that the formose network contains
glyceraldehyde, dihydroxyacetone and a tetrose and nothing longer than four
carbons, and that every reaction conserves atoms and total bond order. The
networks are smaller than the published ones: Chemart's formose network has 11
species where the ECAL paper counts 48. Chemart does not implement the
"dismutations" (Cannizzaro reactions) that the paper's formose network also
uses.

**Energies from structure.** The 2003 paper checks its energy model against
experiment in figure 3: computed and experimental TAEs of the alkanes from
methane to hexane lie on a straight line, and twelve isomers of C6H10 fall in
a narrow band (the caption says "C4H10", but the molecules it lists are C6H10).
Chemart has no experimental values, so it checks what it can. Each added CH2
group changes the TAE by the same amount, about 577 kcal/mol. The twelve isomers span 97 kcal/mol around a
mean of −3203, and 1-methylcyclopentene, the most stable in the caption's
list, is also the most stable in Chemart.

**Absolute energies are not reproduced.** The 2005 paper tabulates TAEs in
kcal/mol, for example −415.95 for ethene and −789.62 for butadiene. Those were
computed with a later parameter set for C, H, N, O, P and S that was only
distributed with a technical report that is no longer online. With the
published 2003 parameters, Chemart gives −1095.6 and −2023.0. For the four
molecules in that table the ratio of published to computed energy stays
between 0.378 and 0.390, and this is what the tests check. Reaction energies
fare worse. The 2005 paper reports that butadiene plus ethene releases 16.76
kcal/mol and two butadienes 8.15; Chemart gives 112.1 and 92.4. Both are
exothermic and ranked in the published order, which the tests check, but the
ratio between them is 1.2 rather than 2.1.

**The size of a network is set by its barriers.** The 2003 ECAL paper grows
the Diels-Alder network by lowering a reactivity threshold (ln k = 108, 91.2
and 55.4 in its figure 3) and follows the network measures up to 246 species
(figure 4). Chemart's `barrier_cutoff` plays that role, and the test checks
that the formose network never grows as the cutoff is tightened. The effect is
abrupt, not gradual: at 131.5 kcal/mol or more the full network forms; below
131.1, the barrier of the first step, only the two seed molecules remain. The
ECAL series is not reproduced.

**Small worlds, or not.** The ECAL paper studies the *substrate graph*: species
are vertices, and two species are joined if they take part in the same
reaction. It compares the mean degree `<k>`, mean path length `<L>` and
clustering `<C>` (how often two neighbours of a species are neighbours of each
other) with a random graph of the same size. Its Table 1:

| network | n | `<k>` | `<L>` | `<L_rand>` | `<C>` | `<C_rand>` |
|---|---|---|---|---|---|---|
| Formose | 48 | 3.25 | 3.55 | 3.28 | 0.15 | 0.068 |
| Diels-Alder | 40 | 4.65 | 2.15 | 2.40 | 0.72 | 0.110 |
| E. coli metabolism | 282 | 7.35 | 2.9 | 3.04 | 0.32 | 0.026 |

Diels-Alder and E. coli are small-world (much higher clustering than random,
and paths no longer); formose is not, because many of its species cannot react
until a keto-enol step makes them reactive, which lengthens paths and lowers
clustering. The paper concludes that chemical networks do not all fall into one
class of small-world networks. It
also finds the degree distributions of both networks close to a power law,
with a slope of −1.19 in its figure 5.

Chemart's 40-species Diels-Alder network gives `<k>` 4.45, `<L>` 2.08 against
2.47 for the random graph, and `<C>` 0.81 against 0.11: small-world, as
published. A slow test checks the same criterion on the default 60-species
budget. The formose result is not
reproduced. Chemart's 11-species formose network also comes out small-world by
the same criterion (`<L>` 2.07 against 2.25, `<C>` 0.47 against 0.29), because
it is much smaller and denser than the published one. The degree distribution
is not tested.

**Collisions instead of rules.** The 2005 paper replaces prescribed mechanisms
by bond-by-bond "collisions". For butadiene and ethene it finds the
Diels-Alder adducts, plus hydrogen shifts it calls artifacts of the
parameters. Chemart does not implement this.

**What the reconstruction simplifies.** Two limits matter if you use the
energies. In Chemart's orbital graph the pi overlaps sit only on double and
triple bonds, so ethene and butadiene come out with the same frontier levels
(HOMO −13.75 eV, LUMO −6.16 eV); conjugation across a single bond does not
lower the gap. And the barriers are not thermodynamically consistent: forward
and reverse barriers of a step do not differ by its reaction energy. The
catalog's implementation decisions below list the other choices, including
how rates are derived from the frontier gap.
