## Introduction

The Network Artificial Chemistry (NAC) is Hideaki Suzuki's proposal, first
presented in 2004, for representing *space* in an artificial chemistry with a
graph. Every node of the graph is a molecule, and an edge between two nodes
means that the two molecules are in contact. Molecules move by having their
edges rewired, so there are no coordinates and no lattice: who touches whom is
all the space there is.

The book introduces NAC through water (book §11.3.3). Much of what makes cells
work, from membrane self-assembly to protein folding, comes from hydrophobic
and hydrophilic interactions: water-loving molecules gather together, and
water-avoiding ones are pushed into clusters of their own. Simulating that
atom by atom is expensive. NAC keeps only the outcome: each node is
*hydrophilic* or *hydrophobic*, and a node may only be wired to a node of its
own kind by the passive rewiring. The book summarises the purpose as modelling
such interactions "in an abstract way using a graph rewiring approach".

Suzuki's 2004 slides give the argument from the other side. They list ten
conditions that he held an artificial-life system must satisfy, six of which
concern spatial structure (a cell wall, cells that can change size, signal
transmission, randomisation, selective transport, compartments that mingle),
and conclude that "very few systems satisfy all of the Conditions (5) to
(10)". He borrowed the rewiring rule from a model of how people get to know
each other through a common friend, which is known to produce a *small-world*
network, and argued that the way acquaintances form resembles the way
molecules collide.

NAC then grew a second, *active* layer: nodes holding small programs that
rewire strong (covalent and hydrogen) edges on purpose. With it Suzuki folded
node chains into working "enzymes" (2007) and made a hydrophilic "network
cell" divide on signals from centrosome-like agents (2008). Chemart implements
only the first, passive layer, because the papers that specify the programs
could not be obtained. What Chemart offers is a simulation model of a graph
that rewires itself under a hydrophilic/hydrophobic constraint.

NAC sits in the book's chapter on bio-inspired chemistries, in the section on
networks. Its neighbours there are quite different: [ToyChem](toychem.md) uses
graphs as structural formulas of single molecules and computes their energies,
and [the conservative random network](bigan-conservative-crn.md) is a network
of *reactions*, not of molecules. In NAC the graph is the reaction vessel.
Spatial chemistries in the next section, such as [Squirm3](squirm3.md), place
their molecules on a grid instead. Suzuki's earlier work with Ono,
[SAC](sac.md), is a string chemistry with cells.

## How it works

### Molecules, contacts and clusters

The state is one undirected graph. A node is a molecule with a single
property, its polarity: `i` for hydrophilic, `o` for hydrophobic. An edge is a
weak bond, the water-mediated contact between two molecules. Its kind follows
from its ends, so no edge type is stored: an edge between two `i` nodes is a
hydrophilic contact and one between two `o` nodes a hydrophobic contact. A
random starting graph also has *mixed* `i`–`o` edges; the rule below removes
them.

Chemistry needs species and reactions, and the book never writes NAC that
way. Chemart's reading is that a *species* is a connected cluster of nodes, up
to isomorphism: two clusters with the same polarities wired the same way are
the same species, whatever the node numbers. A species id lists the
polarities, then the edges, both in a canonical order. `iiioo:0-1.1-2.2-3.3-4`
is five nodes, three hydrophilic then two hydrophobic, wired in a path
0–1–2–3–4. A lone hydrophilic node is just `i`.

### The local rewiring rule

At every step (Suzuki's slides 9 to 13):

1. a *starting node* A is chosen at random among all the nodes;
2. a *stopping node* B is chosen at random among A's neighbours;
3. a *new stopping node* C is chosen at random among the nodes at distance
   two from A, that is, neighbours of neighbours that are not already
   neighbours;
4. edge A–B is deleted and edge A–C is created.

A's contact slides one step further along the graph. NAC adds one condition
(slides 18 and 19): "a hydrophilic node and a hydrophobic node cannot be
wired by the passive rewiring", so when C has the other polarity from A, the
whole step "is canceled" and A–B stays.

The step does nothing when A has no neighbour or nothing at distance two from
it. Two things follow from the rule, and they shape everything a run does:

- **A step never joins two clusters.** A, B and C are always in the same
  cluster, so the number of nodes and edges in a cluster can only be split,
  never pooled. If A–B was the only link between two parts (a *bridge*),
  deleting it splits the cluster in two. A node that loses its last edge stays
  alone for ever, because no rule can reach it.
- **A mixed edge can be deleted but never created.** Edges between `i` and `o`
  nodes that exist at the start disappear one by one as their nodes are
  rewired, and none replace them. This is what makes the two kinds demix.

As a reaction, a step is therefore either an *isomerisation* (one cluster
becomes a differently wired cluster of the same nodes) or a *fragmentation*
(one cluster becomes two).

### A worked example

Take the five-node path `i–i–i–o–o`, nodes 0 to 4, and run 40 steps:

```python
net = chemart.generate_network("nac", seed=1, polarities="iiioo",
                               edges=[[0, 1], [1, 2], [2, 3], [3, 4]], steps=40)
```

Three reactions fired:

```
iiioo:0-1.1-2.2-3.3-4 -> iiioo:0-2.1-2.2-3.3-4            (x5)
iiioo:0-2.1-2.2-3.3-4 -> iiioo:0-1.1-2.2-3.3-4            (x5)
iiioo:0-1.1-2.2-3.3-4 -> iii:0-1.0-2.1-2 + oo:0-1         (x1)
```

The first is A = 0, B = 1, C = 2: edge 0–1 goes, edge 0–2 comes, and node 2
becomes a hub joined to 0, 1 and 3. The second undoes it: from the hub
shape, A = 1 moves its edge from 2 to 0, which gives a path again. The third
is the one that ends the run. A = 2 picks its hydrophobic neighbour B = 3 and
the hydrophilic node C = 0 at distance two. Edge 2–3, the only mixed edge and a
bridge, is deleted and 0–2 is created. The result is a hydrophilic triangle and
a hydrophobic pair. In a triangle every node is next to every other, so nothing
is at distance two, and a pair has no distance two either: both are frozen.

The run's tally, in `net.extras["analysis"]["attempts"]`, shows how the other
steps went: `{'cancelled-polarity': 15, 'moved': 12, 'no-distance-two': 13}`.
The cancelled ones tried to wire an `i` node to an `o` node, for instance A = 1
reaching for C = 3. Twelve steps moved an edge but only eleven appear as
reactions: a move that returns the same cluster up to isomorphism (A = 2
moving its edge from 1 to 0 turns the path 0–1–2–3–4 into the path 1–0–2–3–4)
changes nothing a species id can see, so it is not recorded.

### The reactor

There are no rates. No source gives a time scale for a rewiring, so Chemart
counts steps. Because A is a random *node*, a cluster is rewired in proportion
to its size. Method `rewiring` runs the rule on one graph and returns the
reactions that fired with their counts; method `closure` instead lists every
cluster that can be reached from the starting clusters by single rewirings,
without simulating anything.

## Using it

The default run starts from a random graph of 12 nodes, six of each polarity,
with 18 edges (mean degree 3), and makes 600 rewiring attempts. The first
reaction printed above already sheds a lone `o` node from the one big starting
cluster. What matters is in `net.extras`:

```python
a = net.extras["analysis"]
a["attempts"]
# {'cancelled-polarity': 73, 'moved': 175, 'no-neighbour': 87, 'no-distance-two': 265}
a["trace"][0]
# {'mixed_edges': 8, 'clusters': 1, 'largest_hydrophilic_cluster': 0, 'clustering': 0.15, 'path_length': 2.136364}
a["trace"][-1]
# {'mixed_edges': 0, 'clusters': 4, 'largest_hydrophilic_cluster': 5, 'clustering': 0.777778, 'path_length': 1.1}
net.extras["final_state"]
# {'iiiii:0-1.0-2.0-3.0-4.1-2.1-3.1-4.2-3.2-4.3-4': 1, 'ooooo:0-3.0-4.1-2.1-3.1-4.2-3.2-4.3-4': 1, 'i': 1, 'o': 1}
```

`trace` samples the graph every `n_nodes` steps (`a["sampled_every"]`, here
12). Its fields are the number of mixed `i`–`o` edges, the number of clusters,
the size of the largest all-hydrophilic cluster, the clustering coefficient and
the mean path length (both defined under Results). Here the eight mixed edges
are gone by step 156, and the graph ends as a five-node hydrophilic
cluster with every pair linked, a five-node hydrophobic cluster missing two of
its ten possible links, and one stray node of each kind. The 18 edges are all
still there (10 + 8). Once a cluster is fully linked nothing is at distance
two, which is why most late attempts are `no-distance-two`: the run has
frozen. `extras["space"]` holds the node polarities and the initial and final
edge lists, and `extras["conservation"]` the three conserved quantities
(hydrophilic nodes, hydrophobic nodes, edges) per species.

**Switching the constraint off.** `polarity_constraint=False` is the plain
acquaintance-network rule, a Chemart option rather than a published variant.
On the same seed nine mixed edges remain after 600 steps, and the largest
cluster holds five `i` and three `o` nodes:

```python
loose = chemart.generate_network("nac", seed=1, polarity_constraint=False)
loose.extras["analysis"]["trace"][-1]["mixed_edges"]      # 9
```

**A larger graph.** With 40 nodes, mean degree 4 and 4,000 steps (about 5 s),
mixed edges fall from 37 to 2 by step 400 and to 0 by step 1,000. The run ends
with a 13-node hydrophilic cluster, a 15-node hydrophobic one and 12 isolated
nodes (7 hydrophilic, 5 hydrophobic). Without the constraint (about 9 s) 42
mixed edges remain at the end:

```python
net = chemart.generate_network("nac", seed=1, n_nodes=40, mean_degree=4, steps=4000)
```

Almost every event produces a new cluster shape here, so this run lists 2,918
species; for bigger graphs the species list stops being readable and the
module's graph functions (next recipe) are the better tool.

**The small-world measurement at N = 200.** A network of 200-node clusters is
not useful, so the published measurement is run with the module's helpers.
With one polarity only the constraint never applies, and this is Suzuki's
setting of slide 14 (200 nodes, mean degree K = 10). It takes about 0.3 s:

```python
import numpy as np
from chemart.chemistries.nac import random_graph, rewire, clustering, mean_path_length
rng = np.random.default_rng(0)
labels, adj = random_graph(200, 1.0, 10, rng)
for step in range(2000):
    rewire(labels, adj, rng)
```

```
rewirings  clustering  path length
0          0.053       2.55
500        0.119       2.59
1000       0.141       2.60
2000       0.151       2.59
5000       0.149       2.58
```

**Explicit graphs and the closure.** `polarities` and `edges` (see the
parameter table) set the starting graph by hand, as in the worked example.
`method="closure"` on that same path finds the same three reactions and stops,
`status=complete`. On a random graph of 8 nodes and 6 edges
(`n_nodes=8, mean_degree=1.5`) the closure is complete at 16 species and 24
reactions; at 8 edges (`mean_degree=2.0`) it hits the `max_species` budget of
300 in under a second and returns `status=truncated`.

## Results

### What Chemart reproduces

**Local rewiring makes a small world.** A network is called *small-world* when
it is highly clustered, like a regular lattice, yet any two nodes are only a
few links apart, like a random graph. The two standard measures are the
*clustering coefficient* C (for each node, the fraction of pairs of its
neighbours that are themselves linked, averaged over nodes) and the *mean path
length* L (the average number of links on the shortest path between two
nodes). A random graph with mean degree K and N nodes has C close to K/N,
0.05 for Suzuki's N = 200, K = 10. Slide 14 plots C and L against the number
of rewirings and states: "The local rewiring increases C while keeping L
constant", so that "repetition of the local rewiring can make a small-world
network". This was Suzuki's reason for adopting the rule: it builds a
clustered, local structure out of random contacts, which is what he wanted of a
space. The archived slide's axis values are not legible, so Chemart checks the
claim rather than a curve. In the run above C triples from 0.053 to 0.151
within 2,000 rewirings while L stays near 2.6. Run longer, C drifts slowly back
down (0.131 after 20,000 rewirings and 0.119 after 50,000, per the catalog
decisions), a regime the slide does not show. The test
`test_local_rewiring_makes_a_small_world` asserts that after 2,000 rewirings C
is more than 2.5 times its random value, L is within 15% of its start, and the
edge count is unchanged.

**Hydrophilic and hydrophobic nodes demix.** Because a mixed edge can be
deleted but never created, the number of `i`–`o` edges only falls, and it
reaches zero. The two kinds end up in separate clusters, the qualitative
picture of the book's Figure 11.14(a)–(b), where randomly mixed hydrophilic
and hydrophobic nodes separate and the hydrophilic ones gather. That figure,
however, comes from Suzuki (2008), whose model reaches it with programmed
agents; the passive-rule version here is Chemart's own consequence of slides
18–19, not a published experiment. It also differs in outcome: the rewiring
strands some nodes on their own, so the hydrophilic nodes do not always end in
one cluster (13 of 20 in the 40-node run above). The test
`test_hydrophilic_and_hydrophobic_nodes_demix` checks, on the default seed,
that mixed edges start above zero, never increase and end at zero, that all
but at most one hydrophilic node share a cluster, that clustering rises, and
that with the constraint off mixed edges remain.

**Rewiring conserves nodes and edges.** Every reaction acts on one cluster and
turns it into one or two clusters with the same nodes of each polarity and the
same number of edges. Tests check this on the rule itself, on the stoichiometry
of generated networks against the three conservation vectors, and on the
closure, where a product never has more nodes than its reactant.

### What Chemart does not reproduce

Everything that uses NAC's active layer is left out, because the papers that
specify it (the book's [824], [825] and [827]) are closed access and no open
copy could be found. What is known of them comes from the book and their
abstracts:

- **Folding node chains into enzymes** (Suzuki 2007, book [824]). A genotype
  written as a chain of nodes carrying symbols is "agglomerated, tangled, and
  finally folded into a node cluster that works as a control-flow machine"
  (abstract). The book compares such clusters both to parallel computers and to
  enzymes with a catalytic site. Two were built: a *splitase*, which splits a
  cluster in two (in §18.1.1, splitting a set of molecules into hydrophilic and
  hydrophobic clusters), and a *replicase*, which copies a node chain. The book
  returns to this in §18.1.2, naming NAC as one of the closest artificial
  chemistries to protein folding.
- **Replication and partitioning** (slides 20–21). The 2004 slides sketch
  copying a chain of covalently linked nodes with two active nodes, a
  *polymerase* and a *helicase*, and splitting a hydrophilic cluster walled by
  hydrophobic nodes with a splitase.
- **A dividing network cell** (Suzuki 2008, book [825]). Agents with
  assembler programs of three types (centrosome, hydrogen and van der Waals)
  move through the network. The abstract reports a pseudo-lattice among
  hydrophilic nodes, repulsion between hydrophilic and hydrophobic nodes, and
  "the division of a network cell actualized as a hydrophilic cluster". In the
  book's Figure 11.14 two centrosomes pull the hydrophilic cluster apart into
  two daughter clusters.
- **A network energy** (Suzuki 2006, Australian Journal of Chemistry). To make
  the graph imitate molecules moving in three-dimensional space, Suzuki
  formulated a "network energy", rewired the graph to minimise it and compared
  the result with a hard-sphere random walk (abstract). No formula could be
  obtained, so Chemart has no energy.

### Standing

In the slides Suzuki judged that NAC "almost satisfies" his spatial conditions
(5) to (10), with a question mark on compartments mingling, and listed as
future problems a translation from template to program, a network cell that
synchronises replication and partitioning, making evolution happen, and a
theory linking the rewiring rule to spatial dimension. The book describes the
2007 and 2008 papers and does not mention later uses of NAC by other authors.

## Further reading

- Suzuki, H. (2011). Artificial chemistry and molecular network. In
  *Biological Functions for Information and Communication Technologies*,
  Studies in Computational Intelligence, Springer, pp. 87–161.
  <https://doi.org/10.1007/978-3-642-15102-6_3>. A later chapter by Suzuki
  on the same programme (not consulted for this page).
