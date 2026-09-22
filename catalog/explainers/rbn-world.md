## Introduction

RBN World is an artificial chemistry whose atoms are small dynamical systems.
Each atom is a *random Boolean network* (RBN): a handful of on/off switches,
or *nodes*, wired together at random, each reading a few others through a
truth table drawn at random. Started from any pattern, such a network is
deterministic and has finitely many patterns, so it soon falls into a cycle
it repeats for ever, its *attractor* (the [random Boolean networks](rbn.md)
page explains these networks on their own). RBN World gives each network two
*bonding sites*. Two atoms bond when a number measured on their attractors
meets a criterion, and the bond breaks again when their dynamics change.
Molecules are thus trees of Boolean networks.

RBN World was built by Adam Faulconbridge, Susan Stepney, Julian Miller and
Leo Caves at York between 2009 and 2011; Faulconbridge's PhD thesis (2011) is
its fullest description. The book mentions it twice. Its section on random
Boolean networks (§18.4.2) states the aim: "to construct ACs able to exhibit
rich and complex dynamics, with the emergence of higher-level structure able
to transition to higher orders of organization". Its section on bondable
cellular automata (§10.7.3) credits it with the idea of a *subsymbolic*
artificial chemistry: whether two atoms react is not written in a table but
emerges from the atoms' behaviour.

Its closest relative in the catalog is [Bondable Cellular
Automata](bondable-ca.md), which uses one-dimensional cellular automata as
atoms in the same way. Unlike the [classic RBN](rbn.md), whose network is
drawn once and then written down, RBN World is a *Turing gas*: molecules carry
structure, a procedure computes what a collision produces, and the soup and
its network change as a run goes on.

## How it works

An atom of RBN World is a *bonding RBN* (bRBN): a random RBN of N nodes with K
inputs each, in which two inputs, chosen at random, are rewired to two
*bonding sites*. An empty site reads 0 and a filled one reads 1, so bonding
changes what the network computes. When two atoms are bonded inside a
molecule, the site wiring becomes reciprocal inputs between the two networks,
and the pair runs as one larger network. A molecule is a tree: the atoms are
the leaves, and each internal node is a composite bRBN made of the bRBNs below
it.

To decide whether two bRBNs bond, each is run from its stored state until a
state repeats, and a number is computed over its attractor cycle. The default
*bonding rule*, `proportion-sum-one`, uses the fraction of nodes that are on,
averaged over the cycle, and requires the two fractions to add up to 1.

A collision picks a free site on each molecule. It tests the two atoms that
hold the sites, climbing to larger enclosing composites if the test fails.
If some pair passes, the sites are filled and the test is repeated, because
filling a site changes the dynamics. If it still passes, the molecules join;
if not, the sites are emptied again, but the atoms keep the states they moved
to. Finally every bond in the molecule is rechecked, and those that no longer
pass break, which can split a molecule. The reactor mixes random pairs of
molecules for a fixed number of collisions and records each distinct outcome
as a reaction. Nothing enters or leaves the reactor, so every atom present at
the start is still there at the end, bonded or not.

The default run (seed 1) has five elements, `A` to `E`, 20 copies of each.
Their settled proportions are 0.6, 0.5, 0.55, 0.35 and 0.6. Only `B` pairs
with anything, and only with itself, since 0.5 + 0.5 = 1. So every reaction in
this run involves `B`. The most frequent ones are:

```
2 B.1 -> B.1 + B.2  (x15)
B.2 + E.1 -> B.1 + E.1  (x10)
2 B.1 -> (B-B).1  (x4)
```

The number after the dot numbers the different states of the same structure
in the order they were first seen: `B.1` is `B` settled on its attractor,
`B.2` a second state. In the first reaction two `B` atoms tried to bond, the
test failed once the sites were filled, and one of them was left in a new
state. This is the pattern of the thesis's example `A1 + B → A2 + B`, which it
calls a simple form of catalysis, "not designed into any part of the system":
`B` is unchanged, yet the collision turns `A1` into `A2`. In the second, a
collision with `E` makes `B.2` run its dynamics again, which returns it to
`B.1`. The third is synthesis: two `B` atoms bonded, and the bond survived.

## Using it

RBN World has one face, a run of its reactor. The call printed above,
`chemart.generate_network("rbn-world", seed=1)`, runs the default reactor to
the end and returns the network of every reaction that fired, with `count`
for how often; its `initial_state` is the starting soup, 20 copies of each
settled atom. `chemart.evolve` returns the whole run as a trajectory, with a
frame every 100 collisions (one per molecule at the start):

```python
traj = chemart.evolve("rbn-world", seed=1)
[f.t for f in traj.frames][:3]                 # [0.0, 100.0, 200.0]
traj.frames[1].fired
# [[['B.1', 'B.1'], ['(B-B).1'], 1], [['B.1', 'B.1'], ['B.1', 'B.2'], 1]]
traj.frames[-1].state
# {'A.1': 20.0, 'E.1': 20.0, 'D.1': 20.0, 'B.1': 12.0, 'C.1': 20.0, '(B-B).1': 4.0}
traj.series("largest_molecule_atoms")          # [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
```

A frame's `state` is the soup at that moment, `fired` the reactions since the
previous frame (reactants, products, count), and the observable
`largest_molecule_atoms` the number of atoms in the largest molecule present.
In this run the first bond forms in the first 100 collisions, and the soup
ends with four `B` dimers.

The species ids are readable structures (`(B-B).1`), and each species'
`structure` field holds the full code: the tree, which sites are bonded, and
the node states at every level. `traj.network.extras` also has `elements`
(the wiring and tables of each atom, with `site0` and `site1` marking the
bonding sites), `conservation` (one law per element: atoms are never created
or destroyed) and, in `analysis`, the `final_population`,
`largest_molecule_atoms` (over every species seen) and
`distinct_collision_outcomes`.

**Decomposition.** A larger soup, 100 copies of each element for 5,000
collisions, shows molecules breaking apart. With seed 2 the most frequent
decomposition is

```python
net = chemart.evolve("rbn-world", seed=2, copies=100, collisions=5000).network
```

```
(B-C).1 + B.1 -> C.1 + 2 B.1  (x5)
```

a `B` atom colliding with the dimer `(B-C)` and the bond breaking.

**A chemistry that bonds.** Because the atoms are random, few of them can bond
under the default rule: seeds 1 to 4 give only one to four distinct synthesis
reactions each. `bonding_rule="cycle-length-equal"`, the original 2009
chemistry, bonds far more readily. With the default settings and seed 1, 98 of
the 100 atoms end up in a single molecule, and the frames show it growing:

```python
traj = chemart.evolve("rbn-world", seed=1, bonding_rule="cycle-length-equal")
traj.series("largest_molecule_atoms")   # [1, 4, 21, 34, 98, 98, 98, 98, 98, 98, 98]
```

Seeds 2 and 3 end with molecules of 62 and 39 atoms.

**Speed.** The default run and the 5,000-collision run each take a fraction
of a second. Large molecules make each collision slower, since every bRBN on
the way up the tree is run to its attractor: the cycle-length run with seed 3
takes about 3 seconds.

## Results

**Reactions from dynamics.** Faulconbridge's thesis (2011, ch. 6)
walks through example reactions of increasing complexity: synthesis
`A + B → (A–B)`, bonding onto an existing molecule at an atom, `(A–B–C)`, or at
the whole molecule, `((A–B)–C)`, decomposition `(A–B) + C → A + B + C`, and
`A1 + B → A2 + B`. In the last, an attempted bond fails but leaves `A` in a
different state; the thesis calls multiple states "a key feature of RBN-World
that is often missing from other Artificial Chemistries". Chemart's tests
rebuild the thesis's structure-tree examples, including the nine-atom example
of Faulkner et al. (2018, fig. 12), check that every surviving bond still meets the
rule, and find both synthesis and state-changing catalysis in the default
network on six seeds, with the atom counts conserved, both in the network and
in every frame of a run.

**Linking is not associative.** Because the bonding rule is tested on the
bRBNs at every level of the tree, the order in which atoms join matters:
`((A–B–C)–D)` is in general not `((A–B)–(C–D))` (Faulkner et al. 2018). The
tests check the structural side of this, that bonding `C` to the atom `B` of
`(A–B)` gives `(A–B–C)` while bonding it to the composite gives `((A–B)–C)`.

**Which chemistries are worth studying.** The ALife XII paper (Faulconbridge et
al. 2010) treats the design choices as a search space: six bonding properties
(cycle length, flashing, flashes, total, magnitude, proportion), several
comparison criteria, and atoms of 5 to 25 nodes with K = 2 or 3. It tested
200 alternative chemistries, 10,000 random samples each, for five low-level
behaviours: synthesis, self-synthesis, decomposition, substitution and
catalysis. A chemistry passed a test if its samples showed variation, some
passing and some failing. Only 19
passed all five (the paper's table 5), and they all use either proportion with
sum one or total with sum zero; the original cycle-length/equal choice did not
pass. Chemart offers these two and the original, and its tests reproduce the
paper's worked example of the six properties on a four-node cycle (table 2:
cycle length 6, flashing 3, flashes 8, total −4, magnitude 14, proportion
0.417). It does not rerun the 200-chemistry screen.

**Evolved atoms and long loops.** The aim of RBN World, in the book's words, is
chemistries "able to produce autocatalytic sets, hypercycles and
heteropolymers". Rather than sampling atoms at random, the thesis (ch. 8) uses a
genetic algorithm, 100 vessels for 300 generations, each seeded with 1,000
atoms of each of five types, to find atom sets whose reaction networks contain
long reaction loops. One network analysed in detail had 1,286 reactions and
645 molecular species, and its longest loop had 8 reactions (thesis §8.5.1.1,
repeated in Faulkner et al. 2018). Chemart reproduces neither this search nor
the thesis's Gillespie-style timing; its atoms are random, as in the 2009 and
2010 papers. Nor does it implement the temperature analogue with which
Faulkner et al. report preliminary experiments.
