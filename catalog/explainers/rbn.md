## Introduction

A random Boolean network (RBN) is a set of switches wired together at random.
Each switch, or *node*, is either on or off. It reads the current values of a
few other nodes, its *inputs*, and a fixed rule, a truth table drawn at
random, says whether it will be on or off at the next tick. Stuart Kauffman
proposed these networks in 1969 as a "formal genetic net": a caricature of a
gene regulatory network, in which a node is a gene that is either expressed or
silent, and the wiring says which genes regulate which. The book presents it in
its chapter on modelling genetic regulatory networks (§18.4.2).

Because the number of on/off patterns of N nodes is finite (2^N) and the
update is deterministic, every run eventually repeats a pattern and then cycles
for ever. These repeating cycles are the network's *attractors*. The book reports
two results of Kauffman's. First, how much the network moves depends
mostly on K, the number of inputs per node: with K = 1 it freezes, with large K
it wanders erratically, and around K = 2 it sits at what he called the "edge of
chaos". Second, he conjectured that the few attractors of such a network could
be the distinct *cell types* of an organism: every cell carries the same
genes, but each settles into a different stable pattern of gene activity (book
§18.5).

Chemart offers the RBN as a chemistry in two different ways, under one entry:

- **The RBN written as reactions.** Each node's two states become two
  species, and each row of a truth table becomes a reaction in which the input
  nodes act as catalysts that switch the node. This is the classic model, in a
  form that reaction-network tools can read.
- **RBN World**, built by Adam Faulconbridge, Susan Stepney, Julian Miller and
  Leo Caves at York between 2009 and 2011. Here a whole RBN is one *atom*.
  Two atoms bond when a property of their dynamics, measured on their
  attractors, meets a criterion, and bonds break again when the dynamics change.
  Molecules are thus trees of Boolean networks. The book calls this a
  *subsymbolic* artificial chemistry (§10.7.3): whether two atoms react is not
  written in a table but emerges from the atoms' behaviour.

The first form differs from random-network neighbours such as [random
catalytic networks](random-catalytic-networks.md) in being purely logical,
with no rates or concentrations. The second is the sibling of [Bondable Cellular Automata](bondable-ca.md), which uses
one-dimensional cellular automata as atoms in the same way. The archived entry
organization-computing uses the same two-species-per-variable encoding.

## How it works

### The classic network

An RBN has N nodes. Node i reads K distinct inputs (it may read itself) and
has a truth table with 2^K rows, one for each combination of input values.
Each output in the table is 1 with probability `function_bias`, which is ½ for
the book's "uniformly random" functions. At each tick every node computes its
table at the same time; this *synchronous* update is the book's model.

Chemart turns node i into two species, `x{i}_0` (off) and `x{i}_1` (on), with
exactly one of them present at any time. A table row becomes a reaction: the
input species of that row appear on both sides, so they are catalysts, and the
node's species changes to the row's output. Here is node 1 of the default
network (seed 1). It reads nodes 0 and 8, and its table is `1011`: row r is read
as the binary number `x0 x8`, so the rows are 00 → 1, 01 → 0, 10 → 1,
11 → 1. That gives four reactions:

```
x0_0 + x8_0 + x1_0 -> x0_0 + x8_0 + x1_1      row 00: switch node 1 on
x0_0 + x8_1 + x1_1 -> x0_0 + x8_1 + x1_0      row 01: switch it off
x0_1 + x8_0 + x1_0 -> x0_1 + x8_0 + x1_1      row 10: switch it on
x0_1 + x8_1 + x1_0 -> x0_1 + x8_1 + x1_1      row 11: switch it on
```

Each row appears only in the direction that changes something: if node 1 is
already on, row 00 does nothing, so there is no reaction for it. When a node
reads itself, the rows in which it already shows the output are dropped for
the same reason.

These reactions fire one at a time, so on their own they describe an
*asynchronous* network, one node updated at a time. That is not the book's
model, and asynchronous RBNs generally have different attractors. Chemart
therefore also computes the synchronous map on all 2^N states and reports its
attractors, which is why N is capped at 16. Seed 1 is a quiet network: from its
random initial state it reaches a fixed point after four ticks, and every one
of the 1,024 states ends there.

```
t=0  0001111100
t=1  1101111001
t=2  1101110000
t=3  1101100100
t=4  1101100101      and it stays here
```

Each string lists node 0 first. Seed 3, with the same N and K, is livelier:
it has four attractors, two of length 2 and two of length 4, whose *basins* (the
sets of states that end on each one) hold 490, 294, 220 and 20 of the 1,024
states. In Kauffman's reading, that network has four cell types.

### Frozen, critical, chaotic

The *sensitivity* λ measures how far a small disturbance spreads: it is the
average number of nodes whose output changes when one input is flipped
(Drossel 2008). If λ < 1, disturbances die out and the network freezes; if
λ > 1, they grow and it becomes chaotic. For random tables with bias p,
averaging over the ensemble gives the *annealed* value λ = 2 K p (1 − p), so
K = 2 with p = ½ gives λ = 1 exactly, the critical point. Chemart reports both
this annealed value and the actual average sensitivity of the network drawn.
Lowering p is a second way to make a network orderly, besides lowering K.

The book also suggests replacing the fixed K by a distribution with mean K;
Chemart draws in-degrees from a Poisson or a power-law distribution
(`K_distribution`).

### RBN World

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
as a reaction.

The default run (seed 1) has five elements, `A` to `E`, 20 copies of each.
Their settled proportions are 0.6, 0.5, 0.55, 0.35 and 0.6. Only `B` pairs
with anything, and only with itself, since 0.5 + 0.5 = 1. So every reaction in
this run involves `B`. The most frequent ones are:

```
B.1 + B.1 -> B.1 + B.2       15 times
B.2 + E.1 -> B.1 + E.1       10 times
B.1 + B.1 -> (B-B).1          4 times
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

The default call above is the classic network: N = 10 nodes, K = 2 inputs,
unbiased tables, so it sits at the critical point λ = 1. The species `x{i}_b`
say that node i has value b, and `net.initial_state` holds one of each pair,
the random starting pattern. The dynamics are in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["inputs"][1], a["functions"][1]        # ([0, 8], '1011')
a["attractors"]                          # [{'length': 1, 'basin_size': 1024, 'cycle': ['1101100101']}]
a["initial_attractor"], a["transient_length"]            # (0, 4)
a["average_sensitivity"], a["annealed_sensitivity"]      # (0.8, 1.0)
```

Attractors are listed largest basin first; `initial_attractor` is the index of
the one the initial state falls into. Cycles longer than 64 states are cut to
their first 64.

**From order to chaos.** Sweep K on 12-node networks, 20 seeds each, and look at
the median attractor length (averaged over states, so weighted by basin size),
the median number of attractors and the mean sensitivity:

```python
import statistics, chemart
for K in (1, 2, 3, 5, 12):
    runs = [chemart.generate_network("rbn", seed=s, N=12, K=K).extras["analysis"]
            for s in range(20)]
    L = [sum(x["length"] * x["basin_size"] for x in a["attractors"]) / 2**12 for a in runs]
    print(K, round(statistics.median(L), 1),
          statistics.median(len(a["attractors"]) for a in runs),
          round(statistics.mean(a["average_sensitivity"] for a in runs), 2))
```

```
1 1.0 1.0 0.55
2 2.4 2.5 1.05
3 5.7 4.0 1.56
5 21.1 4.0 2.52
12 41.5 5.0 6.0
```

For K = N = 12, Drossel's estimate of the longest attractors is of order
2^(N/2) = 64. The run takes a few seconds.

**Other ensembles.** Same setup, 12 nodes and 20 seeds. Biased tables move the
critical point: with K = 4, `function_bias=0.15` gives λ = 2 × 4 × 0.15 × 0.85
= 1.02, and the median attractor length falls from 8.1 (at p = ½, λ = 2) to
2.0. With K = 3, the median length is 5.7 for constant in-degree, 4.8 for
`K_distribution="poisson"` and 3.1 for `"power-law"`, although the average
sensitivity is about 1.5 in all three. A power law on 1..N needs
K < (N + 1)/2; larger K raises an error.

**RBN World.** Set `model="rbn-world"`. The species ids are readable
structures (`(B-B).1`), and each species' `structure` field holds the full code: the tree,
which sites are bonded, and the node states at every level. `net.extras` also
has `elements` (the wiring and tables of each atom, with `site0` and `site1`
marking the bonding sites), `conservation` (one law per element: atoms are
never created or destroyed) and, in `analysis`, the `final_population` and
`largest_molecule_atoms`. The default run takes about a second; 5,000
collisions on 100 copies of each element take about two. With seed 2 and those
settings, decomposition appears among the observed reactions:

```
(B-C).1 + B.1 -> C.1 + B.1 + B.1       5 times
```

Because the atoms are random, few of them can bond under the default rule:
seeds 1 to 4 give only one to four distinct synthesis reactions each. `bonding_rule="cycle-length-equal"`,
the original 2009 chemistry, bonds far more readily: with the default settings
and seed 1, 98 of the 100 atoms end up in a single molecule (seeds 2 and 3: 62
and 39).

## Results

**Attractors as cell types.** Kauffman (1969) found cycle attractors in RBNs
near the edge of chaos and conjectured that they correspond to cell types, a
cell's differentiation being its convergence onto one of them; in 1969 there
was no way to test this (book §18.5). Chemart cannot test a biological
conjecture; it supplies the attractors and basins the conjecture is about, and
its tests check that they are exactly the cycles of the synchronous map, with
basins covering all 2^N states.

**Order and chaos.** The book summarises Kauffman's finding that as K grows
from 1 to N, RBNs pass from frozen to complex to chaotic dynamics, with
interesting behaviour at K = 2 or 3. Drossel's review (2008) makes the chaotic
end precise: a K = N network behaves like a random map on its 2^N states, with
attractor lengths of order 2^(N/2). Chemart's tests reproduce the trend in the
median attractor length of 12-node networks over 15 seeds, rising from K = 1
to K = 2 to K = N, with the K = N value above 2^(N/2)/8.

**The sensitivity.** Drossel gives the annealed sensitivity λ = 2 K p (1 − p)
for random tables with bias p, with the frozen and chaotic phases separated by
λ = 1. The tests check that the average sensitivity of drawn networks matches
this formula within 0.08, over 30 seeds, for K = 2, p = ½ (λ = 1) and
K = 3, p = 0.2 (λ = 0.96).

**Scale-free networks.** The book reports, citing Fox and Hill, that power-law
in-degree distributions give more ordered dynamics and fit biological data
better. Chemart's tests only check that its Poisson and power-law in-degrees
have mean K; the recipe above shows the more ordered dynamics in a small run,
but that is not a test, and only the abstract of Fox and Hill's paper was
available when the generator was written.

**RBN World: reactions from dynamics.** Faulconbridge's thesis (2011, ch. 6)
walks through example reactions of increasing complexity: synthesis
`A + B → (A–B)`, bonding onto an existing molecule at an atom, `(A–B–C)`, or at
the whole molecule, `((A–B)–C)`, decomposition `(A–B) + C → A + B + C`, and
`A1 + B → A2 + B`. In the last, an attempted bond fails but leaves `A` in a
different state; the thesis calls multiple states "a key feature of RBN-World
that is often missing from other Artificial Chemistries". Chemart's tests
rebuild the thesis's structure-tree examples, including the nine-atom example
of Faulkner et al. (2018, fig. 12), check that every surviving bond still meets the
rule, and find both synthesis and state-changing catalysis in the default
network on six seeds, with the atom counts conserved.

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

## Further reading

- Fox, J. J. & Hill, C. C. (2001). From topology to dynamics in biochemical
  networks. *Chaos* 11(4), 809–815. Book reference [291], the source of the
  scale-free result.
- Kauffman, S. A. (1969). Metabolic stability and epigenesis in randomly
  constructed genetic nets. *Journal of Theoretical Biology* 22(3), 437–467.
  Book reference [445].
