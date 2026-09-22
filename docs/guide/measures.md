# Measures for comparing chemistries

This page collects the measures you can compute on a reaction network to place
it in a **morphospace**: a space whose axes are measurable properties, where
every chemistry becomes a point, or rather a cloud of points, one per seed.

## Why a morphospace

Chemart gives every artificial chemistry the same output, a reaction network.
That makes it possible to ask a question the field has never answered
systematically: **how many of these chemistries are really different?** Two
chemistries with different names, molecules and reactors may still produce
networks that are indistinguishable by every measure below. In a morphospace
that redundancy is visible: clouds that overlap. Empty regions are visible
too: kinds of chemistry nobody has built.

There are two morphospaces worth building, and comparing them is the point:

- a **design morphospace**, whose axes are how a chemistry is built. They are
  read from the catalog: `type`, `family`, `kind`, `constructive`,
  `provides`, the molecule type in `S`;
- a **behavioural morphospace**, whose axes are the measures on this page,
  computed on generated networks.

Chemistries far apart in design but close in behaviour are the redundant ones.
Chemistries close in design but far apart in behaviour show which design choice
matters.

Real chemical networks (a metabolic model, a combustion mechanism, the formose
reaction, an atmospheric network) can be imported into the same record and
placed in the behavioural morphospace as reference points. They also test the
measures themselves: if a set of measures cannot tell combustion from
metabolism, it is too coarse to judge two artificial chemistries.

## How to use this page

Every measure on this page is a function of `chemart.measures`, listed in the
tables below by the name you call it with. The tables are generated from the
code, so they cannot drift from it.

```python
import chemart
from chemart import measures

net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
chemart.measure(net)                                  # every cheap measure that applies
chemart.measure(net, ["deficiency", "max_raf_fraction"])
measures.applicable(net)                              # {name: None, or why it does not apply}
measures.zscores(net, ["nodf", "clustering"])         # against the null model of rule 2

traj = chemart.evolve("alchemy", seed=1)
measures.over(traj, ["richness", "shannon", "n_reactions"], window=5)   # over time

rows = measures.sweep("random-catalytic-networks", {"n": [10, 20, 40]}, seeds=range(5))
measures.scaling(rows, "n_reactions")                 # section J
```

A measure takes a network, the population state of one frame (a dict of
amounts), or a whole trajectory: the **input** column. Its **cost** is
*cheap* (run by default), *moderate* or *exponential*; pass `cost=` to include
the dearer ones. Measures with a node limit do not run above it unless you pass
`force=True`.

### What each measure needs

Every measure needs some part of the [network record](../reference/record.md).
The **needs** column of each table uses these codes:

| code | what | where it is in Chemart |
|---|---|---|
| **T** | topology: which species take part in which reactions | `net.reactions` |
| **S** | stoichiometry, S = P − R | `ids, R, P = net.matrices()` |
| **C** | catalysts | `reaction.catalysts` |
| **F** | a food set: what is supplied from outside | `net.extras["food"]`, else `net.inflow`, the buffered species or `net.initial_state` (`measures.food_set(net)`) |
| **K** | rate constants | `reaction.rate` |
| **D** | a trajectory | `chemart.simulate.ode`/`ssa` (see [Simulating dynamics](simulating.md)) or `chemart.evolve` |
| **str** | the molecules' internal structure | `species.structure` |

`net.summary()` lists what a network `provides`. Many artificial chemistries
have no rates and no catalysts, so the kinetic and catalytic dimensions will be
missing for them. That is a finding, not a gap to fill: leave the value
missing rather than imputing it.

### Five rules

1. **Choose one graph and keep it.** Use the bipartite species–reaction graph:
   one node per species, one per reaction, an edge from each reactant to its
   reaction and from each reaction to its products. Projecting to a
   species-only graph turns every reaction into a clique that is not in the
   chemistry, and networks with different reaction arities stop being
   comparable.
2. **Report measures relative to a null model.** Raw values mostly measure
   size. Compare each measure with networks randomised so that each species
   keeps its degree and each reaction keeps its numbers of reactants and
   products, and report the z-score.
3. **Clouds, not points.** Compute everything over at least 20 seeds. The
   spread across seeds is itself a feature: some chemistries are far less
   reproducible than others.
4. **Mind the scale.** Networks here range from four species (the
   Brusselator) to thousands. Prefer size-normalised measures, or use how a
   measure *scales* with size as the feature ([section J](#j-scaling-relationships)).
5. **Freeze the set before you look.** Choose the measures before seeing where
   the points land, or the morphospace will just echo what you expected.

## A. Size and composition

These are for normalising and sanity checks. On their own they say how big a
network is, not what kind it is, so they are reference values rather than axes.

<!-- measures A -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `n_species` | Number of species. | T | network | cheap |
| `n_reactions` | Number of reactions. | T | network | cheap |
| `reaction_density` | Reactions per species, r / n. | T | network | cheap |
| `arity` | Share of reactions by molecularity, "reactants->products" (e.g. "2->1": 0.4), counting molecules with their multiplicity. | S | network | cheap |
| `mean_reactants` | Mean number of reactant molecules per reaction. | S | network | cheap |
| `mean_products` | Mean number of product molecules per reaction. | S | network | cheap |
| `reversible_fraction` | Share of reactions whose exact reverse is also in the network. | S | network | cheap |
| `catalysed_fraction` | Share of reactions with a catalyst: a species on both sides. | T | network | cheap |
<!-- /measures -->

## B. Stoichiometric structure

These come from linear algebra on the stoichiometric matrix S = P − R. They are
standard in the theory of real reaction networks, which makes them the most
reliable bridge between artificial and real chemistry.

<!-- measures B -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `stoichiometric_rank` | Rank of the stoichiometric matrix: the dimension of the space the concentrations can move in. | S | network | cheap |
| `rank_ratio` | rank(S) / n: the share of independent directions of change. | S | network | cheap |
| `conservation_laws` | n - rank(S): the number of independent conserved linear combinations of amounts (dimension of the left null space of S). | S | network | cheap |
| `conservative` | Whether a strictly positive mass is conserved: some m > 0 with mᵀS = 0 (a linear programme). | S | network | cheap |
| `deficiency` | Feinberg's deficiency δ = n_c - ℓ - rank(S): complexes minus linkage classes minus the rank. With δ = 0 and weak reversibility, mass action has exactly one positive steady state in each stoichiometric compatibility class, and it is stable, whatever the rates. | S | network | cheap |
| `weakly_reversible` | Whether every reaction lies on a cycle of the complex graph (each linkage class is strongly connected). | S | network | cheap |
| `flux_dimension` | r - rank(S): the number of independent steady-state flux patterns of the closed network (dimension of the right null space of S). | S | network | cheap |
| *p_invariants* (planned) | the semi-positive conservation laws (conserved moieties): extreme rays of {m ≥ 0, mᵀS = 0} | S | network | exponential |
| *elementary_flux_modes* (planned) | number and mean length of the minimal pathways through the network with its food boundary | S F | network | exponential |
| *blocked_fraction* (planned) | share of reactions that can never carry flux at steady state (flux variability analysis) | S F | network | moderate |
<!-- /measures -->

## C. Graph topology

These treat the network as a graph (rule 1). Most of them only mean something
relative to a null model (rule 2).

<!-- measures C -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `degree_cv` | Coefficient of variation (sd / mean) of species degrees: how unevenly species take part in reactions. | T | network | cheap |
| `degree_gini` | Gini coefficient of species degrees. | T | network | cheap |
| `assortativity` | Degree assortativity of the undirected bipartite graph: whether high-degree species meet high-degree reactions. None when undefined (all degrees equal). | T | network | cheap |
| `clustering` | Mean bipartite clustering of the species (Latapy et al. 2008): how much the reaction neighbourhoods of species that share a reaction overlap. | T | network | cheap |
| `reciprocity` | Share of substrate -> product links between species that also run back. None when there are no such links. | T | network | cheap |
| `cycle_rank` | Independent loops of the undirected bipartite graph: edges - nodes + components. | T | network | cheap |
| `flow_hierarchy` | Share of edges on no cycle (Luo & Magee 2011): 1 is purely feed-forward, 0 means everything feeds back. | T | network | cheap |
| `bow_tie` | Shares of species in the core (the largest strongly connected set of the substrate -> product graph), upstream of it (in), downstream of it (out), and elsewhere (other). | T | network | cheap |
| `efficiency` | Global efficiency of the undirected bipartite graph: the mean inverse shortest-path length over all pairs of nodes. | T | network | moderate, up to 5000 nodes |
| `spectral_gap` | Second-smallest eigenvalue of the normalised Laplacian of the undirected bipartite graph (its algebraic connectivity): 0 when it falls apart, larger when it is well connected. | T | network | cheap, up to 20000 nodes |
| `spectral_entropy` | Shannon entropy (nats) of the normalised Laplacian eigenvalues, read as a distribution. | T | network | moderate, up to 3000 nodes |
| `catalytic_spectral_radius` | Largest eigenvalue modulus of the catalytic graph, an arrow from each catalyst to each net product of the reactions it catalyses. At least 1 exactly when the graph has a cycle: an autocatalytic set (Jain & Krishna 1998). | C | network | cheap, up to 3000 nodes |
| `nodf` | Nestedness (NODF, Almeida-Neto et al. 2008) of the species x reaction incidence matrix, 0 to 100: whether rarer species take part only in reactions that commoner species also take part in. | T | network | cheap |
| `mean_hyperedge_size` | Mean number of distinct species per reaction, reactions read as hyperedges. | S | network | cheap |
| *modularity* (planned) | Louvain modularity Q of the bipartite graph, relative to the null model | T | network | moderate |
| *motif_profile* (planned) | significance profile of the triads of the substrate -> product graph (Milo et al. 2004) | T | network | moderate |
<!-- /measures -->

## D. Catalysis, autocatalysis and organisation

These are the central ideas of artificial chemistry: sets of molecules that
make each other and so sustain themselves. They are where artificial
chemistries should differ most, and where comparing them with origin-of-life
chemistry means most.

<!-- measures D -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `max_raf_fraction` | Share of the reactions that belong to the maximal RAF set. | C F | network | cheap |
| `scope_fraction` | Share of species the network can make from its food set (its scope), catalysts not required. | F S | network | cheap |
| `expansion_depth` | Generations network expansion takes to reach the scope of the food set. | F S | network | cheap |
| *irreducible_rafs* (planned) | a lower bound on the number of irreducible RAFs, by sampling | C F | network | exponential |
| *autocatalytic_cores* (planned) | minimal stoichiometric autocatalytic cores (Blokhuis, Lacoste & Nghe 2020) | S | network | exponential |
| *organisations* (planned) | number of chemical organisations and the size of the largest (Dittrich & Speroni di Fenizio 2007) | S | network | exponential |
<!-- /measures -->

## E. Constructiveness and growth

A constructive chemistry keeps producing species that did not exist before. A
single static network cannot show this; a sweep or a trajectory can.

<!-- measures E -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `mean_structure_length` | Mean length of the species' structure strings (a λ-term, a sequence, a fold): a crude size of the molecules in their own representation. | str | network | cheap |
| `novelty_rate` | Species never seen before, per unit of the trajectory's clock. | D | trajectory | cheap |
| `complexity_drift` | Slope over time of the abundance-weighted mean length of species ids: do molecules get bigger as the run goes on? | D | trajectory | cheap |
<!-- /measures -->

## F. Kinetic measures

Only for chemistries with rate constants.

<!-- measures F -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `rate_spread` | Orders of magnitude spanned by the mass-action rate constants, log10(max k / min k). None without mass-action rates. | K | network | cheap |
| `wegscheider_residual` | How far the reversible mass-action pairs are from allowing detailed balance (Wegscheider's conditions): the least-squares residual of ln(k+/k-) against the reactions' stoichiometry, 0 when some chemical potentials make every pair balance. None without reversible mass-action pairs. | K S | network | cheap |
| *steady_states* (planned) | number of steady states found by multi-start root finding | K | network | moderate |
| *stability* (planned) | largest real part of the Jacobian's eigenvalues at steady state, and the stiffness ratio | K | network | moderate |
| *oscillation* (planned) | whether a simulation settles into sustained oscillation, and its period | K | network | moderate |
| *flux_concentration* (planned) | Gini coefficient of the steady-state fluxes | K | network | moderate |
| *sloppiness* (planned) | eigenvalue spread of the Fisher information of the rate constants | K | network | moderate |
| *entropy_production* (planned) | entropy production at steady state over the reversible pairs | K | network | moderate |
<!-- /measures -->

## G. Dynamics and trajectories

For algorithmic chemistries such as AlChemy or Tierra the network is only a
record of what happened, and the interesting part is how the population
changed. These measures capture that.

<!-- measures G -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `richness` | Number of species present. | T | state | cheap |
| `shannon` | Shannon diversity (nats) of the population's composition. | T | state | cheap |
| `dominance` | Share of the most abundant species (Berger-Parker index). | T | state | cheap |
| `population` | Total amount of all species. | T | state | cheap |
| `turnover` | Mean Jaccard distance between the species sets of consecutive frames: how fast the population's composition changes. None with one frame. | D | trajectory | cheap |
| `collapse_time` | First time, after its peak, at which richness falls to `fraction` (10%) of the peak. None if it never does. | D | trajectory | cheap |
| `final_richness_ratio` | Richness at the end over the peak richness of the run. | D | trajectory | cheap |
| *attractor_type* (planned) | fixed point, cycle or chaos, from the recurrence of the trajectory | D | trajectory | moderate |
<!-- /measures -->

## H. Robustness and redundancy

These measure redundancy *inside* a network: how many ways it has of doing the
same thing. They are the measures closest to the redundancy question.

<!-- measures H -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `production_multiplicity` | Mean number of reactions with a net production of each species. | S | network | cheap |
| `percolation` | Area under the curve of the share of species in the largest connected piece as species are removed, at random and highest-degree first (Albert, Jeong & Barabási 2000). A robust network keeps a large area under both. | T | network | moderate, up to 3000 nodes |
| *degeneracy* (planned) | structurally different pathways to each species from the food set | S F | network | moderate |
| *knockout_tolerance* (planned) | share of reactions whose removal leaves the scope and the maxRAF unchanged | S F | network | moderate |
| *synthetic_lethal_pairs* (planned) | pairs of reactions that back each other up (sampled on large networks) | S F | network | moderate |
<!-- /measures -->

## I. Information and algorithmic complexity

<!-- measures I -->
| measure | meaning | needs | input | cost |
|---|---|---|---|---|
| `compressibility` | Compressed size over raw size of the canonical reaction list (zlib, level 9): lower means more regular. | T | network | cheap |
| `degree_entropy` | Shannon entropy (nats) of the species degree distribution. | T | network | cheap |
| *structure_function_mi* (planned) | mutual information between reactant and product structure features | str | network | moderate |
<!-- /measures -->

## J. Scaling relationships

For any measure above, compute it across many instances of the same chemistry
(seeds, or a sweep of its size parameter) and fit how it scales with network
size. Use the exponent as the feature. That removes the size problem of rule 4,
and it is how real biochemistry has been compared across levels of biological
organisation, which makes these the best features for placing artificial
chemistries next to real ones.

`measures.sweep(chemistry, grid, seeds=...)` measures a chemistry's networks
over a grid of its arguments and seeds, one row each;
`measures.scaling(rows, y, x="n_species")` fits the exponent of `y ~ x^a` on
log-log axes and reports it with its R².

## Where to start

About a dozen measures are cheap, apply to real networks as well, and cover
most sections:

- `rank_ratio`, `conservation_laws`, `deficiency`, `arity` (A–B);
- `bow_tie`, `flow_hierarchy`, `clustering`, `nodf`, each against its null model with `zscores` (C);
- `max_raf_fraction` where catalysts exist, `scope_fraction` and `expansion_depth` (D);
- `production_multiplicity` (H);
- the growth exponent: `scaling` over a `sweep` of the chemistry's size argument (E, J).

Add the moderate and exponential measures afterwards, and only for chemistries
small enough to compute them on.

## References

!!! warning "Not yet checked"
    These pointers were written from memory. Check each one before citing it or
    building a test on it.

- Raup, D. M. (1966). Geometric analysis of shell coiling. *Journal of
  Paleontology*. The original theoretical morphospace.
- Feinberg, M. (2019). *Foundations of Chemical Reaction Network Theory*.
  Springer. Deficiency, weak reversibility, the deficiency-zero theorem.
- Schuster, S., Fell, D. A. & Dandekar, T. (2000). A general definition of
  metabolic pathways useful for systematic organization and analysis of
  complex metabolic networks. *Nature Biotechnology*. Elementary flux modes.
- Milo, R. et al. (2004). Superfamilies of evolved and designed networks.
  *Science*. Motif significance profiles.
- Luo, J. & Magee, C. L. (2011). Detecting evolving patterns of
  self-organizing networks by flow hierarchy measurement. *Complexity*.
- Jain, S. & Krishna, S. (1998). Autocatalytic sets and the growth of
  complexity in an evolutionary model. *Physical Review Letters*. Spectral
  radius and autocatalytic sets.
- Hordijk, W. & Steel, M. (2004). Detecting autocatalytic, self-sustaining sets
  in chemical reaction systems. *Journal of Theoretical Biology*. RAF sets.
- Blokhuis, A., Lacoste, D. & Nghe, P. (2020). Universal motifs and the
  diversity of autocatalytic systems. *PNAS*. Stoichiometric autocatalytic
  cores.
- Dittrich, P. & Speroni di Fenizio, P. (2007). Chemical organisation theory.
  *Bulletin of Mathematical Biology*.
- Handorf, T., Ebenhöh, O. & Heinrich, R. (2005). Expanding metabolic networks:
  scopes of compounds, robustness, and evolution. *Journal of Molecular
  Evolution*.
- Albert, R., Jeong, H. & Barabási, A.-L. (2000). Error and attack tolerance of
  complex networks. *Nature*.
- Edelman, G. M. & Gally, J. A. (2001). Degeneracy and complexity in biological
  systems. *PNAS*.
- Almaas, E. et al. (2004). Global organization of metabolic fluxes in the
  bacterium *Escherichia coli*. *Nature*.
- Kim, H. et al. (2019). Universal scaling across biochemical networks on
  Earth. *Science Advances*; Gagler, D. C. et al. (2022). Scaling laws in
  enzyme function reveal a new kind of biochemical universality. *PNAS*.
- Marshall, S. M. et al. (2021). Identifying molecules as biosignatures with
  assembly theory and mass spectrometry. *Nature Communications*. Assembly
  index.
