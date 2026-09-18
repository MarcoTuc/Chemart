# Network Artificial Chemistry (NAC)

`nac` · *Suzuki, 2004-2009*

*Also known as:* *network artificial chemistry*, *Suzuki's network AC*, *graph rewiring artificial chemistry*

A chemistry in which the graph is both the space and the molecules. Every node is a molecule - hydrophilic or hydrophobic - and every edge a weak bond standing in for water-mediated interaction. A step picks a node, an adjacent node and a node two steps away, and rewires one edge between them. Repeating that local move turns a random graph into a small-world one, and forbidding bonds between unlike nodes makes the two kinds demix into clusters. Hydrophobic effects without atomistic cost.

| | |
|---|---|
| **family** | network |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.3.3; folding 18.1.1 |
| **refs** | [824], [825], [827], doi:10.1016/j.biosystems.2006.09.005, doi:10.1016/j.biosystems.2008.05.017, doi:10.1071/CH06182 |
| **provides** | `topology`, `stoichiometry`, `initial-state`, `space`, `mass-conservation`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): a graph stands in for the solvent: every node is a molecule, classified as hydrophilic (i) or hydrophobic (o), and every edge is a weak bond (the hydrophilic or hydrophobic interaction that water mediates). A species is a connected cluster of nodes taken up to isomorphism; its id is the canonical labels followed by the canonical edges, e.g. iio:0-1.1-2 (a lone hydrophilic node is i), and its structure is nodes=iio edges=0-1,1-2.

**R — reactions** (implicit, arity 1): one local rewiring event acting on one cluster: (1) a starting node A is chosen at random, (2) a stopping node B at random among the nodes adjacent to A, (3) a new stopping node C at random among the nodes at distance two from A, (4) edge A-B is eliminated and edge A-C is created. A hydrophilic and a hydrophobic node cannot be wired by this passive rewiring, so the step is cancelled when C has the other polarity. As a reaction this is an isomerisation (cluster -> the same cluster rewired) or a fragmentation (cluster -> cluster + cluster).

**A — reactor**: graph-rewrite
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("nac", seed=1)
```

```
nac: 61 species, 64 reactions, status=observed
provides: initial-state, mass-conservation, space, stoichiometry, topology
seed: 1
extras: analysis, conservation, final_state, space
```

First reactions:

```
iiiiiioooooo:0-1.1-2.1-3.1-10.2-3.2-6.2-7.3-6.3-8.4-6.5-7.5-9.6-7.6-10.7-10.8-9.8-11.9-10 -> iiiiiiooooo:0-1.1-2.1-3.1-10.2-3.2-6.2-7.3-6.3-8.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10 + o  (x1)
iiiiiiooooo:0-1.1-2.1-3.1-10.2-3.2-6.2-7.3-6.3-8.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10 -> iiiiiiooooo:0-1.1-2.1-3.1-6.1-8.2-3.2-10.3-6.3-7.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10  (x1)
iiiiiiooooo:0-1.1-2.1-3.1-6.1-8.2-3.2-10.3-6.3-7.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10 -> iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-8.2-6.2-7.3-10.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10  (x1)
iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-8.2-6.2-7.3-10.4-6.5-7.5-9.6-7.6-8.6-10.7-10.8-9.9-10 -> iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-9.2-6.2-7.3-10.4-6.5-7.5-8.6-7.6-9.6-10.8-9.8-10.9-10  (x1)
iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-9.2-6.2-7.3-10.4-6.5-7.5-8.6-7.6-9.6-10.8-9.8-10.9-10 -> iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-9.2-6.2-7.3-8.4-6.5-7.5-10.6-7.6-9.6-10.8-9.8-10.9-10  (x1)
iiiiiiooooo:0-1.0-2.1-2.1-3.1-6.1-9.2-6.2-7.3-8.4-6.5-7.5-10.6-7.6-9.6-10.8-9.8-10.9-10 -> iiiiiiooooo:0-2.1-2.1-3.2-3.2-4.2-7.2-9.3-6.3-7.4-8.5-6.5-10.6-7.7-9.7-10.8-9.8-10.9-10  (x1)
iiiiiiooooo:0-2.1-2.1-3.2-3.2-4.2-7.2-9.3-6.3-7.4-8.5-6.5-10.6-7.7-9.7-10.8-9.8-10.9-10 -> iiiiiiooooo:0-2.1-2.1-3.2-3.2-4.2-7.2-10.3-6.4-8.5-6.5-9.6-10.7-8.7-9.7-10.8-9.8-10.9-10  (x1)
iiiiiiooooo:0-2.1-2.1-3.2-3.2-4.2-7.2-10.3-6.4-8.5-6.5-9.6-10.7-8.7-9.7-10.8-9.8-10.9-10 -> iiiiiiooooo:0-2.1-2.1-3.2-3.2-4.2-7.2-8.3-6.4-10.5-6.5-9.6-10.7-8.7-9.7-10.8-9.8-10.9-10  (x1)
… and 56 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `rewiring` | structural | rewiring: run the local rewiring rule on one graph and return the reactions that fired with their counts; closure: every cluster reachable from the initial ones by single rewirings (chemart.expand, one alternative per rewiring) <br>one of `rewiring`, `closure` |
| `n_nodes` | `int` | `12` | population | number of molecular nodes in the graph <br>`2` … `400` · *range:* the published small-world measurement uses N = 200 with K = 10 (module helpers random_graph, rewire, clustering, mean_path_length); closures are readable up to about 8 nodes |
| `hydrophilic_fraction` | `float` | `0.5` | structural | fraction of the nodes that are hydrophilic (round(n_nodes * fraction) of them, at random positions); the rest are hydrophobic <br>`0.0` … `1.0` |
| `mean_degree` | `float` | `3.0` | spatial | initial random wiring: round(n_nodes * mean_degree / 2) distinct undirected weak edges drawn uniformly <br>`0.0` … `100.0` · *range:* K = 10 in the published small-world measurement |
| `polarities` | `str` | `` | structural | explicit node polarities over i (hydrophilic) and o (hydrophobic); overrides n_nodes and hydrophilic_fraction <br>*range:* e.g. iiiooo |
| `edges` | `list` | `` | spatial | explicit initial wiring as [[u, v], ...] over node ids 0..n-1; overrides mean_degree <br>*range:* e.g. [[0, 1], [1, 2], [2, 3]] |
| `steps` | `int` | `600` | population | rewiring only: number of attempted rewiring steps, cancelled ones included <br>`0` … `10000000` · *range:* the small-world measurement needs a few thousand rewirings at N = 200; the default graph has demixed after about 600 attempts |
| `polarity_constraint` | `bool` | `True` | structural | true: a rewiring whose new stopping node C has the other polarity is cancelled (the NAC hydrophilic/hydrophobic interaction); false: rewire regardless, which is the plain acquaintance-network rule and does not demix |
| `max_species` | `int` | `300` | structural | closure only: species budget; reactions whose new clusters would exceed it are dropped and the status becomes truncated <br>≥ `1` · *range:* 8 nodes and 6 edges close at 16 clusters; 8 nodes and 8 edges pass 300 in a fraction of a second |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- repeated local rewiring turns a random network into a small-world one: the clustering coefficient rises while the average path length stays about the same (slide 14; at N = 200, K = 10, C goes from 0.053 - the random value K/N - to 0.151 in 2000 rewirings while L stays at 2.55-2.59)
- hydrophilic/hydrophobic demixing: with the constraint on, edges between the two kinds of node are destroyed and never created, so their number falls to zero and the hydrophilic nodes end up in clusters of their own - the 'network cell' of Suzuki (2008) and figure 11.14(b) of the book; with the constraint off, mixed edges persist
- a rewiring conserves the nodes and the edge count of the cluster it acts on: every reaction is an isomerisation or a fragmentation, never a merge
- not reproduced (active layer, sources not obtainable): folding of a node chain into a control-flow cluster, the designed splitase and replicase, centrosome-driven division of a network cell

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 11.3.3 (molecules occupy the nodes of a graph, nodes hydrophilic or hydrophobic, covalent and hydrogen bonds as directed edges and van der Waals interactions as undirected edges, node programs, splitase and replicase, figure 11.14 cell division) and section 18.1.1 (the mathematical folding of [824]: strings folding into a graph of covalent, ionic, hydrogen and van der Waals bonds, the folded enzyme as a control-flow computer, a splitase splitting molecules into hydrophilic and hydrophobic clusters).
- Suzuki, H. (2004). Artificial chemistry on small-world networks. Proceedings of the 18th Annual Conference of JSAI, 2H4-03. Author's presentation slides, the only NAC source obtained in full: the local rewiring rule step by step (slides 9-13), its effect on the clustering coefficient and the average path length at N = 200, K = 10 (slide 14), the passive/active node structure with the neighbour lists cv[], hy[], il[], ob[] (slide 16), the four edge types and which are rewired actively or passively (slide 17), the hydrophilic/hydrophobic cancellation of a rewiring (slides 18-19), replication of a covalent node chain by polymerase and helicase (slide 20) and partitioning of a hydrophilic cluster by splitase (slide 21). http://web.archive.org/web/0id_/http://www.nis.atr.jp/~hsuzuki/papers/2004.06_JSAI.presentation.pdf
- Davidsen, J., Ebel, H. & Bornholdt, S. (2002). Emergence of a small world from local interactions: modeling acquaintance networks. Physical Review Letters 88(12):128701. The transitive local rewiring that Suzuki's rule is taken from (cited on slide 8).
- Suzuki, H. (2007). Mathematical folding of node chains in a molecular network. Biosystems 87(2-3):125-135 (book [824]). Abstract only (paywalled): a genotype expressed as a node chain with symbol data is agglomerated, tangled and folded into a node cluster that works as a control-flow machine; splitase and replicase are the two examples. https://pubmed.ncbi.nlm.nih.gov/17045733/
- Suzuki, H. (2008). A network cell with molecular agents that divides from centrosome signals. Biosystems 94(1-2):118-125 (book [825]). Abstract only (paywalled): molecular agents with assembler programs move through the network and rewire edges according to node type; with the agent types centrosome, hydrogen and van der Waals it organises a pseudo-lattice between hydrophilic nodes, repulsion between hydrophilic and hydrophobic nodes, and the division of a network cell actualised as a hydrophilic cluster. https://pubmed.ncbi.nlm.nih.gov/18616973/
- Suzuki, H. (2006). An approach toward emulating molecular interaction with a graph. Australian Journal of Chemistry 59(12):869-873 (the journal version of [827]). Abstract only (paywalled): NAC's rewiring rules are designed so that edges imitate the relations between molecules in three-dimensional space, a 'network energy' expressing that constraint is formulated, the graph is rewired to minimise it and the result is compared with a hard-sphere random walk. https://doi.org/10.1071/CH06182

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- None of the three references of the book section could be obtained: [824] and [825] (Biosystems, Elsevier) and [827] (ECAL 2005 workshop CD-ROM, never published online; its journal version is the paywalled Australian Journal of Chemistry paper above) are all closed access, no open version exists in Unpaywall, OpenAlex, OpenAIRE or Semantic Scholar, the author's ATR page (www.nis.atr.jp/~hsuzuki) stops in February 2006 and its 2004-2005 paper PDFs are not in the Internet Archive, whose CDX search was offline during this work. Everything below therefore rests on the book, the three abstracts and Suzuki's own JSAI 2004 slides; nothing here is quoted from [824], [825] or [827].
- Scope: Chemart implements NAC's passive layer only - the solvent graph, the local rewiring rule and the hydrophilic/hydrophobic constraint, all of which the slides state in full. The active layer (nodes holding von Neumann programs that rewire covalent and hydrogen edges; the designed polymerase, helicase, splitase, replicase and centrosome agents; the folding of a node chain into a control-flow cluster of [824]) is specified only in the papers that could not be obtained, and is not implemented rather than invented.
- Molecules and reactions: the book calls the nodes molecules and the edges bonds but never writes NAC as a reaction network. Chemart makes a species out of a connected cluster of weakly bonded nodes, identified up to isomorphism, and a reaction out of one rewiring event: isomerisation when the cluster stays whole, fragmentation when the deleted edge was a bridge. This is a Chemart reading of the model, not a published formulation.
- Edge types: the slides give four (covalent and hydrogen, rewired by programs; hydrophilic and hydrophobic, rewired passively) and the book calls the first two directed and the weak ones undirected. Without the active layer only the passive edges remain, and a weak edge's type is fixed by its endpoints (both hydrophilic or both hydrophobic), since the constraint forbids a mixed one from ever being created. Chemart therefore stores one undirected weak edge set and no direction.
- The rewiring rule is taken verbatim from slides 9-13 (A at random; B at random among A's neighbours; C at random among the nodes at distance two from A; A-B deleted, A-C created). The distance-two set is computed before the deletion, as in the slides' order. The step is inert when A has no neighbour or nothing at distance two; extras.analysis.attempts counts moved, cancelled-polarity, no-neighbour and no-distance-two.
- Small-world check of slide 14: read literally (C drawn among all the nodes at distance two from A), the rule takes a random network of N = 200, K = 10 from C = 0.053, L = 2.55 to C = 0.151, L = 2.59 in about 2000 rewirings - the rise of the slide with the path length unchanged - after which C drifts slowly back down (0.148 at 5x10^3, 0.131 at 2x10^4, 0.119 at 5x10^4 rewirings). The slide's plot range is not legible in the archived copy, so the test asserts the published claim in the regime where it holds rather than a plateau. A reading in which C is instead drawn among the neighbours of B (the 'common friend' of the acquaintance network) raises C several times more slowly (0.094 after 10^4 rewirings), so the literal reading is the one implemented.
- The polarity constraint follows slides 18-19 ('a hydrophilic node and a hydrophobic node cannot be wired by the passive rewiring ... the rewiring is canceled'): the whole step is cancelled, so A-B survives. The slides show the cancellation for C's polarity only; Chemart does not test B, which is being unwired. Consequence, and the reason the network demixes: a mixed edge can be deleted but never created.
- Selection is uniform over the nodes for A, over A's neighbours for B and over the distance-two set for C, as the slides say ('randomly chosen'); no weights are given anywhere. chemart.soup.soup is not used because the rule picks a node, not a molecule, so larger clusters react proportionally more often.
- Initial state: a uniformly random simple graph with round(n_nodes * mean_degree / 2) edges and round(n_nodes * hydrophilic_fraction) hydrophilic nodes at random positions, matching the book's figure 11.14(a) ('hydrophilic and hydrophobic nodes are randomly mixed in a network with randomly generated edges'). The published run of slide 14 is N = 200, K = 10; the defaults are much smaller so that cluster ids stay readable and the run is fast.
- No kinetics: no source gives a rate or a time scale for a rewiring, so no reaction carries a rate. Counts in the rewiring method are firing counts of the recorded events.
- Conservation: a rewiring moves one edge inside one cluster, so hydrophilic nodes, hydrophobic nodes and weak edges are conserved by every reaction; the three vectors are in extras.conservation. Clusters never merge, so the closure terminates.
- The v1 parameter initial_edges (matrix) became edges (a list of [u, v] pairs) plus mean_degree; programs (callable) is dropped with the active layer; method, steps, polarity_constraint and max_species are added. polarity_constraint false is a Chemart option, not a published variant: it gives the plain acquaintance-network rewiring of Davidsen et al. and is what shows that the demixing comes from the constraint.
- The 'network energy' of the Australian Journal of Chemistry paper, which NAC minimises to make the graph imitate three-dimensional space, is only described in its abstract; no formula could be obtained, so none is implemented and no energy is reported.

## Notes

An abstract stand-in for water-mediated hydrophobic/hydrophilic interaction without atomistic cost: the graph is both the space and the set of weak bonds. The module also exports canonical, cluster_id, random_graph, rewire, clustering, mean_path_length and mixed_edges, so the published small-world measurement can be run at N = 200 without building a network of 200-node clusters.

---

*Specification: `catalog/chemistries/nac.yaml` · generator: `chemart/chemistries/nac.py` · tests: `tests/chemistries/test_nac.py`*
