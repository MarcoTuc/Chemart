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

### What each measure needs

Every measure needs some part of the [network record](../reference/record.md).
The **needs** column of each table uses these codes:

| code | what | where it is in Chemart |
|---|---|---|
| **T** | topology: which species take part in which reactions | `net.reactions` |
| **S** | stoichiometry, S = P − R | `ids, R, P = net.matrices()` |
| **C** | catalysts | `reaction.catalysts` |
| **F** | a food set: what is supplied from outside | `net.inflow`, or `net.initial_state` |
| **K** | rate constants | `reaction.rate` |
| **D** | a trajectory | see [Simulating dynamics](simulating.md), or the observed runs of soup chemistries |
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

| measure | meaning | how to compute | needs |
|---|---|---|---|
| species *n*, reactions *r* | size | counts | T |
| *r / n* | reaction density | ratio | T |
| reaction arity distribution | how many molecules go in and come out: 1→1, 2→1, 2→2, … | histogram of reactant and product counts per reaction | S |
| fraction of reversible pairs | how much of the network runs both ways | match each reaction with its reverse | S |
| fraction of catalysed reactions | how central catalysis is | reactions with a catalyst, or with a species on both sides | C |

## B. Stoichiometric structure

These come from linear algebra on the stoichiometric matrix S = P − R. They are
standard in the theory of real reaction networks, which makes them the most
reliable bridge between artificial and real chemistry.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| rank(S) / *n* | how many independent directions the concentrations can move in | `numpy.linalg.matrix_rank`; exact rank with sympy for small integer matrices | S |
| number of conservation laws, *n* − rank(S) | quantities that never change: mass, atoms, moieties | dimension of the left null space of S | S |
| conservative network (yes/no) | whether some strictly positive "mass" is conserved | LP: find *m* > 0 with *mᵀS* = 0 (`scipy.optimize.linprog`) | S |
| semi-positive conservation laws (P-invariants) | the actual conserved moieties | extreme rays of {*m* ≥ 0, *mᵀS* = 0}, e.g. with `pycddlib` | S |
| deficiency δ = *n_c* − ℓ − rank(S) | Feinberg's measure from chemical reaction network theory (CRNT). δ = 0 together with weak reversibility guarantees a unique, stable steady state under mass-action kinetics, whatever the rates | *n_c* is the number of distinct complexes (the reactant and product sides of reactions), ℓ the number of connected components of the complex graph | S |
| weak reversibility | every reaction lies on a cycle of the complex graph | the strongly connected components of the complex graph equal its connected components | S |
| flux cone dimension | how many independent steady-state flux patterns exist | nullity of S with boundary reactions for the food set added | S, F |
| number and length of elementary flux modes (EFMs) | how many minimal pathways run through the network; a direct measure of redundancy | `efmtool`, or cobrapy for small models. The count explodes combinatorially: small networks only | S, F |
| fraction of blocked reactions | reactions that can never carry flux at steady state | flux variability analysis in cobrapy, with bounded fluxes | S, F |

## C. Graph topology

These treat the network as a graph (rule 1). Most of them only mean something
relative to a null model (rule 2).

| measure | meaning | how to compute | needs |
|---|---|---|---|
| degree distributions (in/out, per species and per reaction) | how unevenly species take part | histograms, summarised by the coefficient of variation or the Gini coefficient. Avoid power-law exponents: small networks cannot support a fit | T |
| degree assortativity | whether hubs connect to hubs | `nx.degree_assortativity_coefficient` | T |
| clustering | local redundancy of connections | `networkx.algorithms.bipartite` clustering | T |
| mean shortest path, global efficiency | how many reaction steps separate species | `nx.global_efficiency` on the largest component | T |
| bow-tie fractions | shares of species upstream of (IN), inside (core) and downstream of (OUT) the largest strongly connected component; the classic shape of metabolism | strongly connected components plus reachability in networkx | T |
| flow hierarchy | fraction of edges on no cycle: 1 is purely feed-forward, 0 means everything feeds back | edges between different strongly connected components ÷ all edges | T |
| cycle rank, *m* − *n* + *c* | number of independent loops | on the undirected graph: edges − nodes + components | T |
| reciprocity | share of direct back-and-forth links | `nx.reciprocity` | T |
| modularity Q | whether the network splits into semi-independent subsystems | `nx.community.louvain_communities`, then `nx.community.modularity`; report Q relative to the null model | T |
| motif significance profile | which small subgraphs are over- or under-represented: a structural fingerprint | `nx.triadic_census`, as z-scores against randomised graphs, normalised into a vector | T |
| spectral radius λ₁ of the catalytic graph | λ₁ ≥ 1 exactly when the catalysis graph contains a cycle, i.e. an autocatalytic set | largest eigenvalue of the catalyst → product adjacency matrix | C |
| Laplacian spectral gap, spectral entropy | connectivity and how fast things spread; spectra compare across representations more easily than most measures | `scipy.sparse.linalg.eigsh` | T |
| nestedness (NODF) | whether rare species' partners are a subset of common species' partners | NODF on the species × reaction incidence matrix | T |
| hyperedge size distribution | reactions as hyperedges, without the information lost by projecting to a graph | XGI or HyperNetX | S |

## D. Catalysis, autocatalysis and organisation

These are the central ideas of artificial chemistry: sets of molecules that
make each other and so sustain themselves. They are where artificial
chemistries should differ most, and where comparing them with origin-of-life
chemistry means most.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| maxRAF size / *n* | the share of the network that forms a RAF set: reactions that are all catalysed from within the set and fed from the food set | Hordijk–Steel algorithm, polynomial time, about 50 lines | C, F |
| number of irreducible RAFs | how many *different* ways the network can sustain itself | sampling or enumeration; exponential in the worst case | C, F |
| stoichiometric autocatalysis (exists / number of minimal cores) | autocatalysis defined by stoichiometry, without labelled catalysts | existence: LP for a flux *v* ≥ 0 with S_M·*v* > 0 on a subset M of species. Enumerating the minimal cores is expensive | S |
| scope size from the food set | everything reachable from the food set (network expansion) | `from chemart.expand import expand` | S, F |
| expansion depth | how many generations the expansion takes to close | number of expansion steps | S, F |
| number of chemical organisations, height of their lattice | sets of species that are closed and self-maintaining; many organisations means many alternative stable states | closure plus an LP for self-maintenance per candidate set. Exponential: small or reduced networks only | S |
| size of the largest organisation / *n* | how much of the network can persist | as above | S |

## E. Constructiveness and growth

A constructive chemistry keeps producing species that did not exist before. A
single static network cannot show this; a sweep or a trajectory can.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| growth exponent of species count | *n* ~ size^α or *n* ~ time^α | line fit on log–log axes over a sweep of the chemistry's size parameter or run length | T, over a sweep |
| novelty rate | new species per step or per collision | from a trajectory | D |
| structural complexity of molecules | how elaborate molecules get: string length, Lempel–Ziv complexity, λ-term depth, assembly index | computed in each representation's own terms, then rank-normalised so representations compare | str |
| complexity drift | whether molecules become more complex over time | slope of mean complexity along the trajectory | D, str |

## F. Kinetic measures

Only for chemistries with rate constants.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| spread of rate constants | how heterogeneous the kinetics are | range of log₁₀ *k*, or entropy of the *k* distribution | K |
| thermodynamic consistency | whether the rates allow detailed balance (Wegscheider conditions) | for each cycle of reversible reactions, the product of *k⁺/k⁻* around it must be 1 | K, S |
| number of steady states | multistability | multi-start root finding; δ = 0 gives uniqueness (section B) | K |
| stability, stiffness | largest real part of the Jacobian's eigenvalues at steady state; ratio of fastest to slowest timescale | Jacobian eigenvalues | K |
| oscillation (yes/no, period) | dynamical complexity | eigenvalues crossing into instability (a Hopf bifurcation), or peaks in a simulation | K, D |
| flux concentration | whether a few reactions carry most of the flux | Gini coefficient of the steady-state fluxes | K |
| sloppiness | how many parameter combinations actually matter | eigenvalue spread of the Fisher information matrix | K, D |
| entropy production at steady state | how far from equilibrium the system runs | Σ (*J⁺* − *J⁻*) ln(*J⁺*/*J⁻*) over reversible pairs | K |

## G. Dynamics and trajectories

For algorithmic chemistries such as AlChemy or Tierra the network is only a
record of what happened, and the interesting part is how the population
changed. These measures capture that.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| diversity over time | Shannon entropy and richness of the population | per time window | D |
| turnover | how fast the population's composition changes | Jaccard distance between successive windows | D |
| dominance | whether a few species take over | Berger–Parker index: the largest species' share | D |
| time to collapse or fixation | how long diversity survives | first time diversity drops below a threshold | D |
| attractor type | fixed point, cycle or chaos | recurrence analysis, largest Lyapunov exponent | D |
| variance across seeds | how reproducible the chemistry is | spread of any measure on this page across seeds | D |

## H. Robustness and redundancy

These measure redundancy *inside* a network: how many ways it has of doing the
same thing. They are the measures closest to the redundancy question.

| measure | meaning | how to compute | needs |
|---|---|---|---|
| production multiplicity | average number of reactions producing each species | column sums of P (as a 0/1 matrix) | S |
| degeneracy | structurally different pathways doing the same job | EFMs per target species, or node-disjoint paths from the food set to each species | S, F |
| single-knockout tolerance | share of reactions whose removal leaves the scope, maxRAF or organisation unchanged | remove each reaction in turn and recompute | S, C, F |
| synthetic-lethal pairs | pairs of reactions that back each other up | pairwise knockouts; quadratic cost, so sample on large networks | S, F |
| random vs targeted percolation | how the network falls apart when nodes are removed at random vs hubs first | size of the giant component against the fraction removed; report the area under each curve | T |

## I. Information and algorithmic complexity

| measure | meaning | how to compute | needs |
|---|---|---|---|
| compressibility of the reaction list | how regular the network is | compressed size (zlib, LZMA) of a canonical encoding ÷ raw size. Relabel species canonically first, or the naming will bias it | T |
| graph entropies | disorder of the degree distribution or of the spectrum | Shannon entropy of the degree or eigenvalue distribution | T |
| structure–function mutual information | whether a reaction's outcome is predictable from its reactants' structure | mutual information between reactant and product features | str |

## J. Scaling relationships

For any measure above, compute it across many instances of the same chemistry
(seeds, or a sweep of its size parameter) and fit how it scales with network
size. Use the exponent as the feature. That removes the size problem of rule 4,
and it is how real biochemistry has been compared across levels of biological
organisation, which makes these the best features for placing artificial
chemistries next to real ones.

## Where to start

About a dozen measures are cheap, apply to real networks as well, and cover
most sections:

- rank(S) / *n*, number of conservation laws, deficiency, arity distribution (A–B);
- bow-tie fractions, flow hierarchy, null-normalised modularity, motif profile (C);
- maxRAF fraction where catalysts exist, scope size and expansion depth (D);
- production multiplicity and single-knockout tolerance (H);
- the growth exponent (E).

Add EFMs, organisations and the kinetic measures afterwards, and only for
chemistries small enough to compute them on.

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
