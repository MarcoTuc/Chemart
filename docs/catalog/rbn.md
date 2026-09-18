# Random Boolean Networks (RBN) and RBN World

`rbn` · *Kauffman, 1969; RBN World: Faulconbridge, Stepney, Miller & Caves, 2009-2011*

Kauffman's random Boolean networks, read as chemistry two different ways. In the first, each node's on and off states are species and a truth-table row becomes a reaction in which the inputs catalyse a switch of the output; attractors of the network are then interpreted as cell types, and tuning the number of inputs moves the system between frozen, critical and chaotic. In RBN World the network itself is an atom: two such atoms bond when their dynamics are compatible, so molecules are trees of Boolean networks.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 18.4.2; subsymbolic AC 10.7.3 |
| **refs** | [445], [243], [348], [447], [591], [263], [264], [291] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `mass-conservation`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): RBN: two species per node, x{i}_0 and x{i}_1 (node i off / on). RBN World: an atom is a bonding RBN (bRBN: an RBN plus 2 bonding sites); a molecule is a tree of bRBNs joined at bonding sites, and each tree level carries its own node states

**R — reactions** (implicit, arity 2): RBN: s_i(t+1) = f_i(inputs of i). Truth-table row r of node i is the reaction x{j1}_{b1} + ... + x{jk}_{bk} + x{i}_{1-f_i(r)} -> x{j1}_{b1} + ... + x{i}_{f_i(r)} (the inputs catalyse the switch). RBN World: two molecules collide at a random empty site each; a bond forms between the smallest bRBNs whose emergent attractor property satisfies the bonding rule, before and after the sites are filled; then every bond that no longer satisfies it breaks, e.g. A + B -> (A-B), A1 + B -> A2 + B, (A-B) + C -> A + B + C

**A — reactor**: maximally-parallel, well-stirred-multiset
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("rbn", seed=1)
```

```
rbn: 20 species, 40 reactions, status=complete
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: analysis
```

First reactions:

```
x4_0 + x5_0 + x0_0 -> x4_0 + x5_0 + x0_1
x4_0 + x5_1 + x0_0 -> x4_0 + x5_1 + x0_1
x4_1 + x5_0 + x0_0 -> x4_1 + x5_0 + x0_1
x4_1 + x5_1 + x0_0 -> x4_1 + x5_1 + x0_1
x0_0 + x8_0 + x1_0 -> x0_0 + x8_0 + x1_1
x0_0 + x8_1 + x1_1 -> x0_0 + x8_1 + x1_0
x0_1 + x8_0 + x1_0 -> x0_1 + x8_0 + x1_1
x0_1 + x8_1 + x1_0 -> x0_1 + x8_1 + x1_1
… and 32 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `model` | `enum` | `rbn` | structural | rbn: Kauffman's random Boolean network written as reactions (book 18.4.2); rbn-world: the sub-symbolic chemistry of bonding RBNs (observed reactions of a collision run) <br>one of `rbn`, `rbn-world` |
| `N` | `int` | `10` | structural | rbn: number of nodes (attractors are enumerated over all 2^N states); rbn-world: Boolean nodes per atom <br>`1` … `16` · *range:* RBN World papers use N in {5, 10, 15, 20, 25}; the thesis's evolved atoms use N = 10 |
| `K` | `int` | `2` | structural | inputs per node, or the mean in-degree when K_distribution is not constant <br>`1` … `12` · *range:* 1 <= K <= N; K = 1 frozen, K = 2 critical ('edge of chaos') at function_bias 0.5, K = N chaotic; RBN World uses K in {2, 3} |
| `K_distribution` | `enum` | `constant` | structural | rbn only: in-degrees fixed at K, 1 + Poisson(K - 1) capped at N, or a power law P(k) ~ k^-gamma on 1..N with mean K (scale-free GRNs, Fox & Hill) <br>one of `constant`, `poisson`, `power-law` |
| `function_bias` | `float` | `0.5` | structural | probability that a truth-table output is 1; the annealed sensitivity is 2 K p (1 - p) <br>`0` … `1` |
| `bonding_rule` | `enum` | `proportion-sum-one` | structural | rbn-world only: emergent property and criterion that decide bonds. proportion-sum-one and total-sum-zero passed all low-level tests (ALife XII); cycle-length-equal is the original ECAL 2009 chemistry <br>one of `proportion-sum-one`, `total-sum-zero`, `cycle-length-equal` |
| `elements` | `int` | `5` | structural | rbn-world only: number of random atom types (bRBNs), named A, B, C, ... <br>`1` … `26` · *range:* Faulkner et al. (2018) seed the reactor with five atom types |
| `copies` | `int` | `20` | population | rbn-world only: initial molecules of each element <br>`1` … `1000` · *range:* the thesis uses 1000 copies of each element |
| `collisions` | `int` | `1000` | population | rbn-world only: random pairwise collisions simulated; the observed network grows with it <br>`0` … `100000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- cycle attractors interpreted as cell types (Kauffman's cell-differentiation conjecture)
- order/chaos transition with K: K = 1 frozen, K = 2 critical, K = N chaotic with attractor lengths of order 2^(N/2) (Drossel)
- annealed sensitivity lambda = 2 K p (1 - p), critical at lambda = 1
- scale-free in-degree distributions give more ordered dynamics (Fox & Hill)
- RBN World: synthesis, decomposition, and emergent catalysis with state change, A1 + B -> A2 + B (thesis 6.2.6)
- RBN World linking is not associative: ((A-B-C)-D) != ((A-B)-(C-D)) in general
- only proportion/sum-one and total/sum-zero chemistries show variation in all five low-level tests (synthesis, self-synthesis, decomposition, substitution, catalysis; ALife XII)
- RBN World aims at autocatalytic sets, hypercycles and heteropolymers built from bRBNs; GA-evolved atoms gave a network of 1286 reactions and 645 particles with a longest loop of 8 reactions (Faulkner et al. 2018)

## Sources

- Banzhaf & Yamamoto (2015), section 18.4.2: the classic RBN (N nodes, K inputs with 1 <= K <= N, random functions, synchronous update, distributions of K with mean K)
- Faulconbridge, A. (2011). RBN-World: sub-symbolic artificial chemistry for artificial life. PhD thesis, University of York. Chapters 5-8 (bRBN construction, reaction pseudocode figs. 6.1.1-6.1.4, example reactions 6.2.1-6.2.8, bonding properties and criteria, reactor). https://etheses.whiterose.ac.uk/2061/
- Faulconbridge, A., Stepney, S., Miller, J. F. & Caves, L. (2010). RBN-World: the hunt for a rich AChem. ALife XII, 261-268 (book ref [263]). Tables 1-5. https://www-users.york.ac.uk/~ss44/bib/ss/nonstd/alife12_rbn.pdf
- Faulkner, P., Krastev, M., Sebald, A. & Stepney, S. (2018). Sub-symbolic artificial chemistries. In Inspired by Nature, Springer. Section 5 (bRBN-World). https://www-users.york.ac.uk/~ss44/bib/ss/nonstd/ssachem2018.pdf
- Drossel, B. (2008). Random Boolean networks (book ref [243]), eq. 10: sensitivity lambda = 2 K p (1 - p); K = N networks are chaotic. https://arxiv.org/abs/0706.3351
- Aldana, M. (2003). Boolean dynamics of networks with scale-free topology. Physica D 185:45-66: power-law in-degrees P(k) = k^-gamma / zeta(gamma). https://arxiv.org/abs/cond-mat/0209571

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Two models under one id: model=rbn is the book's RBN; model=rbn-world is RBN World, which the book names in one paragraph and is reconstructed from the thesis, the ALife XII paper and Faulkner et al. (2018).
- RBN network form: node i has species x{i}_0 and x{i}_1 (presence/absence, as in organization-computing), and each truth-table row whose output differs from node i's state in that row is one reaction in which the row's input species catalyse the switch. Rows where node i already shows the output would be null reactions and are omitted. No rates: the model is Boolean.
- A synchronous update of all nodes is not a reaction network. With one molecule per node the reactions give node-by-node (asynchronous) dynamics, so the book's synchronous map is analysed in extras.analysis: attractors (length, basin size, cycle states, listing capped at 64 states), the random initial state's attractor and transient, inputs, functions, in-degrees, average and annealed sensitivity. N is capped at 16 because attractors are enumerated over all 2^N states.
- Inputs are K distinct nodes and may include the node itself: the book allows 1 <= K <= N, and Faulkner et al. (2018) use 'no duplicates, possibly itself'. Kauffman (1969) excluded self-inputs and the thesis allowed duplicates. Truth-table row r reads the inputs as bits with the first input highest (thesis table 5.2.1).
- function_bias is not in the book's text. It is the biased-function ensemble of Drossel [243], where the order/chaos boundary is 2 K p (1 - p) = 1. The book's 'uniformly random function' is function_bias = 0.5.
- K distributions: the book asks for distributions with mean K, and Fox & Hill [291] compare delta, Poisson and power-law in-degrees (only their abstract was accessible). poisson: k = 1 + Poisson(K - 1), capped at N, which keeps k >= 1 as the book requires. power-law: P(k) ~ k^-gamma truncated to 1..N (Aldana's form), with gamma solved so the mean is K. That needs K < (N + 1) / 2, otherwise ValueError; K = 1 gives all k = 1. Only in-degrees vary; inputs are still drawn uniformly.
- v1 parameters dropped: update (the book specifies synchronous update, and the reaction network is the same for any update scheme) and bonding_sites (b = 2 in every RBN World source).
- bRBN atoms (thesis 5.3.1; Faulkner et al. 5.1.2): a random N-node, K-input RBN plus two bonding nodes, each replacing one input of a random node. Not stated: the two sites take distinct input slots. Node initial states are random with probability 0.5 (thesis 5.2.1) and shared by every copy of an element. Atoms are settled onto their attractor before entering the reactor, so the first species of each element is its settled state A.1.
- Bonding rules: only the three published chemistries are offered. proportion/sum one (thesis ch. 8 and Faulkner et al. 2018, the default) and total/sum zero are the only property/criterion pairs that passed all low-level tests (ALife XII table 5); cycle length/equal is the original ECAL 2009 choice, taken as exact equality ('bonds only exist between bRBNs that have the same cycle length', thesis ch. 6). Sum tolerance is 0.001. The other combinations of ALife XII table 3 (equal/similar/different on flashing, flashes, magnitude, etc.) are not offered: they normalise by a property range p_min..p_max whose value for composites, which differ in size and cycle length, is not specified. cycle_properties() still computes all six properties.
- Property evaluation (thesis 6.2.1, 6.2.6): a bRBN is run from its stored state until a state repeats; that first repeated state becomes the stored state, and the property is taken over the cycle. Empty sites read 0, filled sites read 1, and inside a composite a bond becomes reciprocal inputs between the two nodes fed by the bonded sites (thesis 5.3.3). Every tree level stores its own node states. Proportion and total count the N nodes of each atom in the bRBN, bonding nodes excluded. If no attractor is found within 65536 steps the property is undefined and the rule fails.
- Reaction (thesis figs. 6.1.1-6.1.4): pick a random empty site on each molecule. Test the atoms holding them, then climb to the larger enclosing composite (both when sizes are equal) until the rule holds or no larger composite is left. Fill the two sites and retest; if the rule fails, empty them again, keeping the changed states (A1 + B -> A2 + B). Otherwise zip the two trees level by level (Faulkner et al. fig. 12). Finally break every bond whose two child bRBNs fail the rule, repeating until none breaks; split disconnected composites and remove single-child ones. A new composite takes its states from the nodes it replaces, and a split one inherits them (not stated).
- Faulkner et al. (2018)'s temperature analogue is described as preliminary and is not implemented. Linking and decomposition succeed with probability 1, as in all the earlier RBN World work.
- Reactor: iterative mixing (thesis fig. 8.2.1) of copies molecules per element for the given number of collisions, recording each distinct non-elastic reaction with its count (status observed). The thesis's Gillespie-like timing and its genetic-algorithm search for rich atom sets (ch. 8) are not implemented; atoms are random, as in ECAL 2009 and ALife XII. RBN World defines no rate constants. extras.conservation gives one atom-count law per element.
- Species identity: a canonical code of the tree (children sorted), with site-level bonds and the node states at every level. Two molecules with the same element tree but different bonding sites or states are different species. Ids are the readable structure plus an index in order of first observation (A.1, (A-B).2), extending the thesis's state numbering to whole molecules. Isomorphic molecules with internal symmetries could in rare cases get two ids.
- The ECAL 2009 paper [264] was not accessible (Springer). Its model is the one described in chapters 5-6 of the thesis by the same author.

---

*Specification: `catalog/chemistries/rbn.yaml` · generator: `chemart/chemistries/rbn.py` · tests: `tests/chemistries/test_rbn.py`*
