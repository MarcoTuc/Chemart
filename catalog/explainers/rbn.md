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

Chemart offers the RBN written as reactions. Each node's two states become
two species, and each row of a truth table becomes a reaction in which the
input nodes act as catalysts that switch the node. This is the classic model,
in a form that reaction-network tools can read. It differs from random-network
neighbours such as [random catalytic networks](random-catalytic-networks.md)
in being purely logical, with no rates or concentrations. The archived entry
organization-computing uses the same two-species-per-variable encoding.

The same section of the book also describes a second use of RBNs in
chemistry, [RBN World](rbn-world.md), in which a whole RBN is one *atom* and
atoms bond when their dynamics are compatible. It has its own entry.

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

## Using it

The default call above draws a network of N = 10 nodes with K = 2 inputs and
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

## Further reading

- Fox, J. J. & Hill, C. C. (2001). From topology to dynamics in biochemical
  networks. *Chaos* 11(4), 809–815. Book reference [291], the source of the
  scale-free result.
- Kauffman, S. A. (1969). Metabolic stability and epigenesis in randomly
  constructed genetic nets. *Journal of Theoretical Biology* 22(3), 437–467.
  Book reference [445].
