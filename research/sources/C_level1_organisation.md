# C. Level 1 and the step to level 2: what an organisation is, what makes it a unit, and why organisations fail to compose

*Literature report for Chemart, covering the redundancy hypothesis and measures. Compiled 2026-09-23.*

**How this was checked.** Every reference in §6 was checked in this session. The checks were one or more of: a Crossref DOI lookup, the arXiv/PMC/Europe PMC abstract page, or the full text. Full texts read: Fontana & Buss 1994a, 1994b and 1996 (the authors' own copies); Mathis et al. 2024; Vimal et al. 2025; Agüera y Arcas et al. 2024; Kruszewski & Mikolov 2021; Hordijk, Steel & Kauffman 2012; Steel 2015; Dittrich & Speroni di Fenizio (2005 preprint); Montévil & Mossio 2015; McMullin 2004; Speroni di Fenizio, Dittrich & Banzhaf 2001 (abstract and introduction); Vasas et al. 2012; and Banzhaf & Yamamoto chapters 9, 12 and 13.

Text in quotation marks is verbatim. Anything I add that does not come from a source is marked **[synthesis]**.

A note on terms. In Fontana & Buss, *Level 0* is not "molecules". It is "self-copying objects or simple ensembles of copying objects" (hypercycles). *Level 1* is a constructive, self-maintaining organisation that arises when copying is suppressed. *Level 2* is a "metaorganization" of Level-1 organisations. Marco's hierarchy (molecules → networks → networks of networks) therefore puts molecules below F&B's Level 0. His "networks" are F&B's L1, and his "networks of networks" are F&B's L2.

---

## 1. Landscape in brief

1. **"Organisation" at level 1 has about ten formal definitions, and almost all are maintenance criteria.** Each asks whether every member of a set is regenerated from inside the set (closure of production), sometimes with an added flux-feasibility condition. The definitions differ in four ways:
   - what counts as "production": catalysis (RAF), stoichiometry (COT and autocatalytic cores), constraints (closure of constraints), or efficient causes ((M,R) systems);
   - whether dynamics enter;
   - whether a boundary is required;
   - whether reproduction is required.

   The formal links that have been proved:
   - a closed RAF is a COT organisation, but not the other way round [HSD18];
   - under mild conditions any RAF is stoichiometrically autocatalytic [Gol26, preprint];
   - every ODE fixed point is an instance of an organisation [DS07];
   - an attractor that does not touch the boundary lies inside one organisation [PD11].
2. **Being maintained is not the same as being a unit.** The properties that let a level-1 entity act as a building block come from other traditions:
   - a boundary and individuation (autopoiesis, chemoton, closure of constraints);
   - reproduction of the whole (dynamic kinetic stability, chemoton, major transitions);
   - heredity with many possible states, and selectability (Szathmáry; Vasas et al.).

   None of the network measures Chemart already has (maxRAF, irrRAFs, cores, organisations, scope, knockout tolerance) tests these properties.
3. **In the primary sources, the AlChemy ceiling is mainly a problem of missing units, not missing networks.** Fontana & Buss's L1 organisations "are self-maintaining, but not reproducing", and "in no sense can one identify multiple instances of the same L1 organization in our flow reactor". L2 was built by the experimenter, who merged two separately grown L1 organisations. It helped to give the glue different collision rules, and spontaneous L2 was rare. Mathis et al. (2024) repeated the merger on 455 pairs and found coexistence rare.
4. **Every robust higher-level organisation in an artificial chemistry that I found had space, compartments, or group reproduction built into the design.** Examples are spatial hypercycles, Stringmol in space, a combinator chemistry on a planar graph, the Ono–Ikegami cells, and Avida colonies. I found no well-mixed AC that robustly produced level 2 from scratch (§4.7).
5. **A strict "closed and self-maintaining" test is the wrong detector for a finite constructive soup.** Fontana & Buss say this themselves. L1 organisations are infinite invariant subspaces, and only a finite part is ever present, "hence the set of object species in the reactor is no longer closed under interaction". Mathis et al. repeat the point. Graded, time-windowed, flux-weighted measures and perturbation assays are needed instead (§5).
6. **There is a direct link to the redundancy hypothesis.** Fontana & Buss's necessary "condition for organizing" is many-to-one convergence: the number of distinct products must grow more slowly than the number of collision histories. If every collision produces a new object, "no network can ever be formed". This convergence can be measured, and it is the first place where tunable redundancy should show an effect (§4.9, §5.6).

---

## 2. Organisation concepts compared

| Concept (source) | Criterion (exact where quoted) | What it needs as input | What it captures | What it misses (relevant to Chemart) |
|---|---|---|---|---|
| **λ-organisation, Levels 0/1/2** (Fontana & Buss 1994a, 1994b, 1996) | "an organization is a kinetically self-maintaining algebraic structure". This has two parts: (i) "algebraic, a network of mutual production pathways that is a fixed-point under applicative interaction", and (ii) "kinetic, the concentrations of the expressions in the network core are maintained positive" (1996). The organisation is characterised by "(i) boundaries established by the invariances, (ii) strong self-repair capabilities responsible for a robustness to perturbation, and (iii) a center, defined as the smallest set kinetically persistent and self-maintaining generator set of the algebra" (1994a). L2 is "self-maintaining metaorganizations composed of Level 1 organizations". | The reaction rule itself (so the product of any pair can be computed); flow-reactor trajectories; grammar and "laws" worked out by hand | Organisations that are infinite and only partly present ("kinetic confinement"); self-repair; a syntactic boundary; generator and seeding sets; hierarchical nesting | Reproduction (L1 and L2 are "self-maintaining, but not self-reproducing"); mass conservation; space; different reaction rates. There is no algorithmic test. |
| **Chemical organisation (COT)** (Dittrich & Speroni di Fenizio 2007) | *Closed*: every reaction whose reactants are all in O produces only members of O. *Self-maintaining*: there is a flux v with v_r > 0 for every reaction whose reactants are in O, v_r = 0 for every other reaction, and (Sv)_i ≥ 0 for every i in O. The 2005 preprint calls this condition "mass-maintaining" and uses "self-maintaining" for a weaker version. | A finite, explicit stoichiometric network, including inflow and outflow; an ODE for the dynamical theorems | The lattice of organisations; flux feasibility; fixed points and interior attractors lie inside organisations | Catalysis gets no special role (a catalyst is just a reactant and a product). It asks whether some v *exists*, not whether it is realised. It misses stability, small copy numbers, infinite constructive networks and units. The join of two organisations is not their union (§4.8). |
| **RAF** (Hordijk & Steel 2004; wording from Steel 2015) | R′ is non-empty and, for every r in R′, "each reactant of r and at least one catalyst of r is either present in F or able to be constructed from F by using just reactions from within the set R′". Variants: CAF (the catalyst must already be present), pseudo-RAF (can persist but cannot arise from F), closed RAF. | Reactions, which molecule catalyses which reaction, and a food set F | Collective catalytic closure; decomposition into irrRAFs and subRAFs; detection in polynomial time | Stoichiometry, dynamics, inhibition and outflow. A RAF is not closed, so products can leak. "Since the union of any collection of RAFs is also an RAF", combining RAFs always succeeds and so says nothing about integration. In an infinite system an infinite RAF can contain no finite RAF (Steel 2015). It does not distinguish maintenance from reproduction (Kruszewski & Mikolov). |
| **Closed RAF ⇒ organisation** (Hordijk, Steel & Dittrich 2018) | A closed RAF "contains each reaction r for which each reactant and at least one catalyst are either generated by another reaction or are part of the food set F". Result: "closed RAFs are chemical organizations, but … the converse is not necessarily true." | A RAF plus stoichiometry | A bridge between the two theories: COT tools can list all closed RAFs | The same dynamical blind spots as COT |
| **Stoichiometric autocatalysis and cores** (Blokhuis, Lacoste & Nghe 2020; Andersen, Flamm, Merkle & Stadler 2021; Golnik et al. 2026) | The submatrix is "autonomous" (each reaction has a reactant and a product in the set), and some flux makes the net production of every member strictly positive. "An autocatalytic core is an autocatalytic motif which is minimal because it does not contain any smaller autocatalytic motif." There are five core types. | A stoichiometric network with catalysis written out as explicit steps | Self-amplification (growth rather than just persistence); minimal motifs; kinetic viability; autocatalysis that only appears across coupled compartments | "Stoichiometric conditions do not guarantee that autocatalysts within motifs amplify." Definitions differ in subtle ways between authors: the difficulty is to "distinguish autocatalysis e.g. from the superposition of a cycle that consumes and produces equal amounts of X and a pathway that produces X" (Andersen et al.). No boundary and no unit. |
| **Viable core plus periphery** (Vasas et al. 2012) | A core is one or more linked autocatalytic loops that are strongly connected. It is viable if its autocatalysts use reactants not produced only by the loop itself. The periphery is the set of species the core catalyses. | A catalytic network, food, and compartment dynamics | A candidate heritable "genotype" of a network; several cores give several heritable states | It is only selectable with compartments that divide, and it gives only limited heredity |
| **ACS core and periphery** (Jain & Krishna 2001, 2002) | An autocatalytic set is a subgraph in which every node has a catalytic link coming in from inside the subgraph. It has an "irreducible 'core' surrounded by a parasitic 'periphery'". Its size is tracked by λ1, the largest (Perron–Frobenius) eigenvalue of the catalytic graph. | A catalytic graph plus population dynamics with extinction and replacement | Growth by accretion; fragility through "keystone" species; "core shifts"; λ1 as a signal of fragility | Stoichiometry. Composition only happens as one set grows. |
| **Autopoiesis** (Varela, Maturana & Uribe 1974; McMullin 2004) | A network of productions of components that "(i) participate recursively in the same network of productions of components which produced these components, and (ii) realize the network of productions as a unity in the space in which the components exist". This wording is quoted from secondary sources; the original PDF was not read. The paper also gives a six-point key. | A spatial model in which components and a boundary can be identified | A boundary the system makes itself, and so individuation. The unity is "alive regardless of whether it reproduces or not" (1974 abstract). | Reproduction and heredity, deliberately left out. McMullin found that the 1974 model depended on an undocumented "chain-based bond inhibition", and that such a cell "will not easily support growth, or, more particularly, cell reproduction by fission". |
| **Closure of constraints** (Montévil & Mossio 2015; Moreno & Mossio 2015) | A constraint is something that acts on a process and is conserved at the relevant timescale. A set C of constraints realises closure if every C_i "depends directly on at least one other constraint belonging to C" and "there is at least one other constraint C_j belonging to C which depends on C_i". *Strict* closure adds: "C cannot be split into two closed sets." | A split into processes and constraints, with timescales. In an AC, catalysts are the natural constraints. | Irreducibility (one system, not two); a graded "tendency to closure" for drawing boundaries; levels defined as closures that are "both separated and hierarchically nested" | No general algorithm; stoichiometry and kinetics are left implicit |
| **(M,R) systems, closure to efficient causation** (Rosen; Letelier et al. 2006; Jaramillo et al. 2010) | Metabolism f, repair Φ and replication β produce one another. f is the solution of a fixed-point functional equation. | Maps between sets. Small networks can be written in RAF terms (Jaramillo et al.). | Catalysts that the system itself produces ("organizational invariance") | A solution exists only if the set of admissible metabolisms is "drastically smaller than Rosen's own analysis suggested". Uniqueness of the replication map "becomes harder and harder to obtain as the network grows in size". No dynamics. |
| **Chemoton** (Gánti 1975; B&Y §6.1; Szathmáry 2006) | Three autocatalytic subsystems coupled by stoichiometry: a metabolic cycle, template polymerisation, and a membrane. Together they grow and divide. | A stoichiometric model with a compartment | A unit of reproduction with a boundary and a hereditary subsystem. Szathmáry: minimal living systems "must comprise at least a metabolic subsystem, a hereditary subsystem and a boundary". | It is a design, not a detection criterion. It presupposes its parts. |
| **Dynamic kinetic stability (DKS)** (Pross; Pascal, Pross & Sutherland 2013) | "a particular kind of stability … which is not usually observed in regular chemistry, and which is reflected in the persistence of entities capable of self-reproduction", maintained far from equilibrium | The population dynamics of reproducing entities | Persistence through reproduction: the stability at population level that selection needs | No structural criterion and no unit boundary |

**How the concepts relate.**
- **Fontana & Buss vs autopoiesis.** Fontana & Buss (1996) call their organisation "arguably indistinguishable" from autopoiesis, except for one thing: "Our organization[s] are indeed bounded, but bounded syntactically (i.e., λ-organizations are special invariant subspaces of λ-space). A bounding is indeed a necessary feature of organizations, but the space need not be 3-space."
- **COT vs Fontana & Buss.** COT formalised Fontana & Buss's closure and self-maintenance for explicit, finite sets of species. It dropped two of their ingredients: kinetic confinement and the syntactic boundary. This loss is behind the detection problem in §5.4.
- **The maintenance/reproduction split.** Fontana & Buss, COT, RAF and autopoiesis are maintenance criteria. DKS, the chemoton, and Kruszewski & Mikolov's "self-reproducing" metabolisms add reproduction. Kruszewski & Mikolov note that "while the concept of autocatalytic set captures both patterns that perpetuate themselves in time and patterns that also multiply their numbers, it does not explicitly differentiate between them."

---

## 3. What makes a level-1 entity a unit at level 2

### 3.1 Criteria from the literature

| # | Criterion | Primary statement (source) | Why level 2 needs it | Chemart status and possible measure |
|---|---|---|---|---|
| U1 | **Closure of production** (self-maintenance) | F&B; COT; RAF (§2) | The minimum: the entity regenerates itself | Exists (organisations, RAF, cores) but only in strict form. Add graded versions (M1, M2 in §5.6). |
| U2 | **Kinetic persistence under flow or dilution** | F&B's centre is the smallest "kinetically persistent" generator set. COT's flux condition. DKS. | A unit must outlast the turnover of its parts | New: how long members stay; similarity of the soup to itself after a lag (Mathis) (M3) |
| U3 | **Self-repair, regrowth from a small generator** | "Organizations often repair themselves following removal of even large portions of their component expressions" (F&B 1996). F&B describe "seeding sets" and "generators of the algebra". | Regrowing from a few copies is exactly what splitting into offspring needs | Partly present: Chemart's knockout tolerance is a network-level measure (I did not check whether it runs dynamics). Add a dynamic delete-and-regrow test and measure seeding-set size (M4). |
| U4 | **A membership boundary** | F&B: invariances act as "abstract boundaries of the organization. They determine membership"; a foreign expression is either "stably integrated" or "diluted out". Autopoiesis requires a spatial boundary. Speroni di Fenizio et al. 2001 found membranes made of "elastic" (non-reacting) molecules. | If foreign objects are integrated indiscriminately, cross-products dissolve the unit | New: rate at which injected molecules are incorporated; a "grammar" signature of members (M4) |
| U5 | **Individuation: several separable instances** | F&B: "In no sense can one identify multiple instances of the same L1 organization in our flow reactor." McMullin: "autopoiesis requires self-generated 'individuation'", tested by asking "whether two putatively individual cells, in direct contact with each other, can reliably maintain their separate identities". | Selection needs a population of instances | New. Needs space or compartments. Count instances; run a contact test (M9). |
| U6 | **Reproduction of the whole** | F&B: "Reproduction at the new object level requires a means for separating two instances of a L1 object." Szathmáry 2015: "new units of reproduction emerge, and establishment of such units requires high fidelity of reproduction (as opposed to mere replication)". Ono & Ikegami: cells that divide. | DKS; growth of the population of units | New: growth in the number of instances; division events |
| U7 | **Heredity with more than a few states** | Szathmáry 2000: "Systems with limited heredity have only a limited evolutionary potential because the number of available types is too low." Vasas et al. 2010 on compositional genomes: "replication of compositional information is so inaccurate that fitter compositional genomes cannot be maintained by selection". Vasas et al. 2012 on cores. | Adaptations must be able to accumulate | New: number of alternative viable organisations or cores; how faithfully composition passes from parent to offspring |
| U8 | **Selectability** | Vasas et al. 2012: "only when a chemical reaction network consists of many such viable cores, can it be evolvable". Loops inside one core "cannot be independent targets for natural selection". Takeuchi & Hogeweg 2009 on multilevel selection. | Level-2 structure has to be favoured, not just possible | New: variance in persistence or growth across instances that is inherited |
| U9 | **Ability to compose through glue** | F&B's glue. Hordijk, Steel & Kauffman's "meta-RAF": "one set enabling (catalyzing) the existence of another, in mutually beneficial ways". Peng et al. 2020: "pairs of autocatalytic cycles can exhibit competitive, predator-prey, or mutualistic associations just like biological species". | Level 2 means components linked by the products of their interaction | New: merger assay, glue measures, and a test for mutual enablement (M7, M8) |
| U10 | **Separated levels and loss of autonomy** (strong level 2 only) | Montévil & Mossio: "Two closed regimes constitute two different levels of organisation if they are both separated and hierarchically nested." Goldsby et al. 2012: "individuals cease to be able to perform tasks in isolation". Szathmáry 2015. | Distinguishes a real level 2 from coexistence | New: after a merger, put each component back in isolation and check whether it still persists (M7) |
| U11 | **Individuality in information terms** | Krakauer et al. 2020: individuals are "aggregates that preserve a measure of temporal integrity, i.e., 'propagate' information from their past into their futures". This allows individuals that "do not necessarily have physical boundaries". | Finds units in well-mixed or distributed systems where no membrane exists | New: information-theoretic measure (M9) |

**Criteria the literature says are hard to meet in well-mixed networks.** U5–U8 in particular:
- The Vasas et al. 2012 model needed "many compartments" that grow and divide.
- Hordijk, Steel & Kauffman describe the recombination of subRAFs as happening in "(presumed to be compartmentalized replicating entities)".
- Szathmáry 2000 names the "paradox of specificity": "while their abstract feasibility seems to require a high number of molecular types, the harmful effect of side reactions calls for a small system size. No satisfactory solution to this problem is known."

### 3.2 Two strengths of "level 2" (keep them apart in Chemart)

- **Weak level 2: an ecological composite.** This covers F&B's L2, the meta-RAF, and Peng/Baum's mutualistic cycles.
  - Two organisations A and B remain organisations inside a larger one, and are joined by a glue that is "not self-maintaining, but which acts to knit the self-maintaining Level 1 sets into a higher-order self-maintaining entity" (F&B 1994b).
  - F&B stress that the parts "retain their capacity for self-maintenance" (1994a).
  - The glue "is not an organization. It does not preserve its structure when the supporting L1-organizations are removed" (1994a).
- **Strong level 2: a transition in individuality.**
  - A new unit of reproduction appears at the higher level (Szathmáry 2015).
  - The components lose autonomy (Goldsby et al. 2012).
  - The two closures are separated and nested (Montévil & Mossio 2015).
  - No AlChemy L2 reported so far is strong in this sense. **[synthesis]** Marco's "networks of networks" should say which of the two it means. The measures differ: U9 for weak level 2; U5, U6 and U10 for strong level 2.

---

## 4. The ceiling: what the primary sources observed

### 4.1 Fontana & Buss (1994a Bull Math Biol; 1994b PNAS; 1996 "Barrier of objects")

**Set-up.**
- A flow reactor of 1000–3000 λ-expressions. A collision (A)B is reduced to its normal form, the product is checked against filters, and a random object is removed to keep the size constant (this removal acts as dilution).
- The model "does not consider spatial constraints, conservation laws, or unequal reaction rates" (1994b).

**Level 0.** Copy actions produce hypercycles, which are unstable when perturbed. F&B note that the spatial version of Boerlijst & Hogeweg (1991) counters parasitism.

**Level 1.**
- "Level 1 experiments are identical to Level 0 experiments except that copying functions (i.e., Level 0 entities) are barred from action."
- Copying pre-empts level 1: "If identity functions are either present upon initialization of the system or allowed to arise early during the course of an experiment, the system typically does not reach Level 1." Level 1 can still arise if copying is restricted or if copiers "also support constructive interactions".
- Properties: self-repair, seeding sets, emergent "laws" (a grammar plus an algebra), and robustness to injected random objects.
- "Variation in λ-organizations is highly constrained. Specifically, perturbations rarely alter the algebraic laws that characterize a system" (1994a).

**Level 2: protocol and outcome** (1994b, verbatim).
- "Level 2 experiments are initiated with the products of two different Level 1 experiments. The procedure is otherwise identical to the Level 1 protocol, except that the system is increased to a constant size of 3000 objects. Such experiments have one of two outcomes: either a single Level 1 organization comes to dominate the system or a new self-maintaining metaorganization arises (hereafter referred to as Level 2)."
- The glue "contains objects that result from the communication (cross-interaction) between the Level 1 organizations and that do not belong to either organization."
- L2 "laws" are new, and "the organizational description … is not a superposition of the descriptions of the Level 1 organizations".
- The shared metabolism changes the diversity inside each component.
- L2 is "resistant to small perturbations". Where perturbation did change something, it *simplified* the centre of one component.

**How L2 was actually obtained** (1994a, verbatim).
- "Our experience is that L2-organizations are more readily attained by first generating L1-organizations separately, and subsequently combining them into the same reactor. L2-organizations, however, do not always arise under such conditions. It greatly facilitates the construction of L2-organizations, if the glue, or even the constituent L1-organizations operate under different boundary conditions, in particular different collision rules."
- When L2 fails: "one of the two organizations displaces the other, depending on the relative magnitude of the transformation flows between them (if any), as well as on their internal growth rate."
- One spontaneous case is shown: two L1 organisations arose in the same reactor of 1000 objects and formed an L2. No rate is given.
- Banzhaf & Yamamoto (§9.1) summarise: "The spontaneous emergence of a L2-organization from a randomly initialized population is extremely rare. However, it can be synthetically generated by merging two independently emerging L1-organizations."
- A compact statement from 1994a: "Centers compose, organizations not."

**What the authors say is missing.**
- "While L1 organizations maintain themselves kinetically and constructively, they do not reproduce. In no sense can one identify multiple instances of the same L1 organization in our flow reactor … Reproduction at the new object level requires a means for separating two instances of a L1 object" (1994a).
- The L1 and L2 results depend on two conventions: reactants are not consumed ("food"), and products that cannot react are banned ("waste"). For both, "If this assumption is relaxed, self-maintaining organizations fail to emerge" (1994b).
- The 1996 paper lists further limits of the model: no shape, no symmetry, no mass action, no reaction classes, no rate constants. Its linear-logic successor "MC2 has yet to be implemented". Kruszewski & Mikolov (2021) knew of no empirical work on it.

**A condition for organisation that bears on redundancy** (1994a §7).
- An extreme case: if every collision product is its own equivalence class, "no network can ever be formed. The system's diversity explodes like an ever branching tree."
- The condition: "for organization to occur it is necessary that the total number of normal forms … grows with lower order than the total number of collision sequences".
- The number of collision sequences always grows exponentially. For organisations maintained in the flow reactor, the number of distinct normal forms "increases not faster than polynomially".
- The chemistry therefore has to be many-to-one. F&B's abstract lists "chemistry's diversity of equivalence classes, that many different reactants can yield the same stable product" as one of the two abstractions the whole theory rests on.

### 4.2 Mathis, Patel, Weimer & Forrest 2024 (Chaos; arXiv 2408.12137)

**Exact claim** (abstract): "complex, stable organizations emerge more frequently than previously expected, that these organizations are robust against collapse into trivial fixed-points, but that these stable organizations cannot be easily combined into higher order entities."

**What they ran and measured.**
- The original code, recompiled. Soups of 1000 expressions, run for about 10^6 collisions. Then 100 random expressions were added and the run continued; this was repeated five times (6×10^6 collisions in total).
- *L0 runs (no filters).* Many runs ended with 10 to hundreds of distinct expressions (one steady state had about 380) and were unaffected by perturbation. The distribution spans orders of magnitude.
- *Robustness.* A random p% of the soup was replaced with the identity λx.x. The organisations "required very large perturbations to be destroyed"; even 90% replacement destroyed "only some, but not most".
- *L1 runs.* Similarity was measured with the Jaccard index on expression sets. They note that it "does not account for the relative abundance".
  - "In general, even when the number of unique expressions in the simulation remains stable over time, the internal structure of the organization is not as stable—the expressions are continually changing."
  - Changing only the random seed made runs "slowly drift apart into distinct states", which were then stable.
  - Bistability was rare.
- *L2 runs.* Pairs of end-state L1 soups (L and R) were combined in a reactor of up to 2500 objects and run for 10^6 collisions, for 455 pairs. Outcomes were classified as:
  - Dominance: similarity to one input only;
  - Coexistence: mean similarity above 0.1 to both inputs;
  - Mutual Destruction: below 0.1 to both.
  - Conclusion: "the organizations produced by AlChemy can possibly coexist, but it rarely occurs for organizations evolved in different simulations."
  - The exact frequencies are in a table inside the image of Fig. 4B and could not be extracted as text.
- *Their detectors* are not closure tests: number of distinct expressions, survival (anything left besides the identity), Jaccard similarity to self over time and to the inputs, population entropy, and mean expression length.
- *On closure:* "L1 organizations are not necessarily closed under interaction; at any given time in an L1 organization the composition of two expressions could yield a new expression not currently in the system. Over time, however, these new expressions will be diluted out."
- *On detecting L2:* "L2 organizations may be difficult to identify when they emerge spontaneously because deciding whether an organization can be split into two distinct stable organizations is a difficult problem to solve without resorting to trial and error."

### 4.3 Vimal, Mathis, Weimer & Forrest 2025 (arXiv 2509.03534, "Prebiotic functional programs")

- They steer AlChemy with "amplifier" functions that act like unit tests, synthesising Church addition and successor from S/K/I/P combinators.
- "Trickster" functions arise that pass the tests without computing the target.
- "We propose an amplification solution, not a generation solution. Indeed, a key difficulty in amplifying the 'add' function was the fact that it rarely emerged spontaneously."
- Their diagnosis of the wider ceiling: "prebiotic chemical experiments often exhibit diverse forms of self-organization, but open-ended evolution of those forms is limited because they lack genetic systems (Szathmáry, 2006)."
- Selection here acts on molecules, not on organisations.

### 4.4 Agüera y Arcas et al. 2024 (arXiv 2406.19108, "Computational Life": BFF, Forth, Z80, 8080, SUBLEQ)

**What happens after replicators appear.**
- A "state transition": the number of distinct tokens drops sharply and "high-order entropy" (a compression-based measure) rises. At 0.024% mutation, this happened in 40% of runs within 16k epochs.
- In the traced BFF run:
  - the first replicator appears;
  - a "zero-poisoning" phase follows, in which "Replication stagnates and complexity degrades";
  - then a more robust family of replicators takes over;
  - "The soup is a very busy place now: it's full of different replicator versions, constantly overwriting each other."
- On a 2D Z80 grid: "Some of these self-replicators form a sort of symbiotic ecosystems, while other compete for domination … We often observe a series of state transition-like events when more and more capable self-replicators or replicator collectives overtake the soup multiple times."

**The authors' own summary.**
- "We also showed anecdotal evidence that this is the beginning of more complex dynamics."
- "The behavior of such systems is markedly different from auto-catalytic networks and biologically-inspired systems."
- They do no network-level (level-1) analysis.

**[synthesis]** In F&B's terms, BFF is a world dominated by Level 0. F&B found that copiers pre-empt Level 1 unless copying is suppressed or the copiers also build other things. So BFF's first barrier to *networks* may be copying pre-emption, before any brittleness of molecules comes into play. Mathis et al. partly weaken this: in their AlChemy runs without filters, complex organisations often survived even when flooded with the identity.

### 4.5 Kruszewski & Mikolov 2021 (Artificial Life 27(3–4), published 2022; combinatory chemistry)

**System.** SKI combinators with conservation of mass. Weak reduction allows recursive structures that never finish reducing.

**Structures found.** Simple autopoietic, recursively growing, and self-reproducing "metabolisms", for example (AA)+3A ⇒* 2(AA)+φ(A).

**Their critique of AlChemy.**
- "each level of organization was only reached after external interventions."
- Without the food and waste conventions, "complex organizations fail to emerge."

**Blind spots of RAF** that they identify.
- Recursively growing structures cannot be detected "because the resulting expression is not exactly equal to the original one, it still involves a structure that preserves in time its functionality".
- RAF does not separate maintenance from reproduction.

**Their detector.** Rates at which reactants are consumed, scored as pointwise information against a chance model of random collisions and cleavages.

**Result.** They report no level 2.

### 4.6 Well-mixed network models of evolvability

- **Jain & Krishna 2001, 2002.**
  - "A small autocatalytic set, appearing by chance, provides the seed for the spontaneous growth of connectivity and cooperation"; the set then spreads through the whole graph.
  - Crashes follow "the chance elimination of 'keystone' species".
  - "The largest eigenvalue of the adjacency matrix of the graph is an important signal of network fragility or robustness."
  - Composition happens by the growth of one set, not by combining units.
- **Segré, Ben-Eli & Lancet 2000 (GARD).** Assemblies with "compositional genomes".
- **Vasas, Szathmáry & Santos 2010.** Such systems lack evolvability: the system "cannot substantially depart from the asymptotic steady-state solution already built-in in the dynamical equations."
- **Vasas et al. 2012.**
  - "The concept of an autocatalytic or RAF set, although important for questions of self-organization, does not directly address heredity or selectability."
  - Inhibition "produced multiple attractors in an autocatalytic set that cannot be selected for".
  - Evolvability appears only with many viable cores inside compartments that divide: "Acquisition of cores by rare chemical events, and loss of cores at division, allows macromutation, limited heredity and selectability."

### 4.7 Where higher-level organisation did happen, and what enabled it

| System | What the higher level looked like | What enabled it | Emergent or imposed |
|---|---|---|---|
| AlChemy (Fontana & Buss 1994a,b) | L1 + L1 + glue metaorganisation | The experimenter merged them; a larger reactor; different collision rules for the glue; one spontaneous case | Mostly imposed |
| AlChemy re-run (Mathis et al. 2024) | Coexistence in some merged pairs | Merger | Imposed, and rare |
| Combinator chemistry on a planar graph (Speroni di Fenizio, Dittrich & Banzhaf 2001) | "effectively separated components … kept separated by elastic reactions from molecules generated inside the component itself", read as proto-cells. "Special conditions may arise under which different molecules that would have destroyed themselves in a well-stirred reactor can instead coexist … This gives rise to higher organizational levels." | Locality: molecules only interact along the edges of a graph | Emergent |
| Hypercycles on a cellular automaton (Boerlijst & Hogeweg 1991) | Spiral waves that resist parasites | Space | Emergent |
| RNA-like replicators (Takeuchi & Hogeweg 2009) | Travelling waves (implicit higher level) compared with protocells (explicit). Both "achieve the macroscopic stability of a replicator system through the evolutionary dynamics on mesoscopic entities that counteract that of microscopic entities". | Space, or compartments | Emergent (space) or imposed (vesicles) |
| Stringmol in space (Hickinbotham, Stepney & Hogeweg 2021) | Complex replication strategies, discrimination of self from non-self, "complex ecosystems". Parasitism "would lead to extinction unless prevented by compartmentalization or spatial patterning". | Space, starting from a "hand-designed replicator" | Emergent from a designed seed |
| Lattice cell model (Ono & Ikegami 2000) | "A metabolic cycle will produce a self-assembling membrane that will enclose the metabolic cycle"; the cells can reproduce | Space plus membrane chemistry | Emergent. B&Y note "no explicit genetic material and no evolution mechanism". |
| Viable-core model (Vasas et al. 2012) | Compartments carrying several cores, selected against each other | Compartments that divide | Imposed |
| Avida (Goldsby et al. 2012) | Colonies with division of labour, whose individuals lose autonomy | "400 competing 'colonies'"; "A colony that collects a designated number of units of resources … divides into two colonies" | Imposed |
| Z80 soup on a 2D grid (Agüera y Arcas et al. 2024) | "replicator collectives", "symbiotic ecosystems" | Locality in 2D | Emergent, described as anecdotal |

**Answer to "did any AC robustly produce level 2, and what enabled it?"**
- I found no well-mixed artificial chemistry that robustly produced level-2 organisation from random initial conditions.
- Robust cases always involved (a) space, (b) compartments or boundaries, whether imposed or emergent, or (c) group-level reproduction imposed by the design.
- The only well-mixed cases are the merger experiments in AlChemy, and they need the experimenter.
- Banzhaf & Yamamoto (§13.2.2) sketch the same route. The outcome of competition inside a reactor is contingent, so "Now suppose this is part of a larger system, maybe a system of spatially distributed reactors, that are allowed to interact with each other, perhaps on a slower time-scale. Their interaction might lead to different outcomes, depending on which of the two self-replicators has won."

### 4.8 Why organisations fail to compose: the mechanisms named in the sources

1. **Competitive exclusion when two organisations share one well-mixed capacity.** F&B: "one of the two organizations displaces the other, depending on the relative magnitude of the transformation flows between them (if any), as well as on their internal growth rate." Mathis et al. found dominance to be common. Peng et al. (2020) show that autocatalytic cycles behave like species, including "ecological precedence, which makes a system's trajectory historically contingent on the order in which cycles are seeded".
2. **Cross-reactions reopen closure.**
   - In COT, the join of two organisations is generated from their union: close it, then contract to the largest self-maintaining subset. The result can gain species (novel cross-products) and lose them (B&Y §12.1.3: "In the most general case A ⊄ D nor D ⊄ A").
   - The static lattice only gives an upper bound. The dynamics can "move down" to A or to B alone.
   - F&B: "Centers compose, organizations not."
3. **No individuation, so no selection at L1.** A soup holds one instance of each organisation (F&B), so there is no population of L1 units for L2 to be built on by evolution. It can only be assembled by chance or by hand.
4. **Limited heredity even when units exist.** Holistic or ensemble replicators have few heritable types (Szathmáry 2000), and compositional heredity is too inaccurate (Vasas et al. 2010).
5. **Parasites and side reactions.**
   - F&B's L0 ensembles fail through parasitism.
   - Ensemble replicators face the "paradox of specificity" (Szathmáry 2000).
   - Spatial patterning or compartments are the known fixes (Boerlijst & Hogeweg; Hickinbotham et al.).
6. **No separation of scales.** Levels need closures that are "separated and hierarchically nested" (Montévil & Mossio). A single well-mixed reactor offers no second spatial or temporal scale on which a higher closure could sit.

**[synthesis]** Mechanisms 1, 3 and 6 are all properties of well-mixed reactors, not of the molecules.

### 4.9 What this implies for the redundancy experiment [synthesis]

- **Split the hypothesis in two.**
  - (H1) Redundancy increases how often level-1 organisations form and how robust they are. This is plausible, through F&B's need for convergence onto equivalence classes and through resistance to parasites and copiers.
  - (H2) Redundancy lifts the L1→L2 ceiling. The literature above puts that ceiling mainly in missing individuation, reproduction and heredity at L1 (U5–U8). Every case where these were overcome had space or compartments.
- **Suggested design.** A factorial experiment: redundancy × {well-mixed, lattice or graph, dividing compartments}. Score each run on U1–U11, not on counts of organisations.
- **Before interpreting "0.00–0.02 organisations per run", measure convergence.** Track how the number of distinct products grows against the number of reaction events, using the equivalence relation the neural chemistry actually uses (M5). If products hardly ever coincide, F&B predict no network at all, whatever the redundancy of the molecules.
- **Check two conventions.** (a) Do reactants survive reactions? This is F&B's "food". (b) Are inert or unreactive products removed? This is F&B's "waste". F&B found that "self-maintaining organizations fail to emerge" without both.
- **Check copying.** If the soup contains copiers, they may pre-empt level 1.

---

## 5. Detecting organisations in dynamics

### 5.1 What the theory guarantees, and what it does not

- **Fixed points lie in organisations** (Dittrich & Speroni di Fenizio). For an ODE ẋ = Sv(x), every fixed point is an instance of an organisation. The set of species is read off with a threshold Θ, and "for practical reasons, it makes often sense to apply a positive threshold greater zero, e.g., when we take into consideration that the number of molecules in a reaction vessel is finite."
- **The lattice partitions state space.** Every state generates exactly one organisation, G(φ(x)), so "a lattice of organizations partitions the state space X". Dynamics can then be read as movement through the lattice. Speroni di Fenizio & Dittrich 2002 is the origin of this idea (bibliographic entry verified; text not read); B&Y §12.4 develops it.
- **The converse does not hold.** "Even if each fixed point is an instance of an organization, an organization does not necessarily possess a fixed point."
- **Interior attractors** (Peter & Dittrich 2011, from the abstract as indexed; full text not read): any attractor that does not touch the boundary of state space, in particular any periodic attractor, lies inside one organisation.
- **Static organisations are a superset of the stable ones.** B&Y §12.7: "the organizations that can be found by the structural calculations and tests discussed here are a superset of what is really interesting, namely, dynamically stable organizations."
- **Finite simulations can end somewhere else.** B&Y §13.2.1: ODEs "always produce the smallest closed superset of the initial distribution of types, before falling back to the largest self-maintaining subset of this superset … That is not necessarily the case for an explicit molecular simulation where the finiteness of the reaction vessel might lead to the disappearance first of certain sorts of molecules and thus ultimately to other results."

### 5.2 B&Y §13.3 "Observing Organizations"

**The method.**
- Give each known organisation an ID and a concentration vector o_i, normalised so that o_i·o_i = 1.
- Normalise the current state x(t) the same way.
- Report the "strength" of each organisation as the dot product, o_i(t) = x_N(t)·o_i, which is a cosine similarity.

**Limits the authors state.**
- The organisations and their concentration signatures have to be known in advance.
- "Some organizations have a near-infinite number of different expressions in terms of concentrations, because they are dependent on initial conditions … it is not possible to define organization vectors in a unique way."

**[synthesis]** In a constructive chemistry, the organisations are not known in advance. The method still works for tracking organisations found in earlier runs, for example across a merger assay.

### 5.3 B&Y §13.4 "Probabilistic Notions of Closure and Self-Maintenance"

**The definitions.**
- Degree of closure: α_C = m_C / r, where m_C is the number of the set's r reactions whose products stay inside the set.
- Degree of self-maintenance: α_S = m_S / r.

**Their use in constructive systems.** Closure cannot be defined when the set of objects may be infinite, "But can we define a degree of closure and a degree of self-maintenance? Based on empirical evidence, perhaps we can!" This means observing the multiset over time and computing the degrees at each step. The authors also note that keeping these degrees high "naturally leads to the necessity of filters (and membranes)."

**Caution [synthesis].** As printed, α_S counts "reactions … that produce an element of S_B". That makes it the same quantity as α_C. A real degree of self-maintenance has to measure regeneration of the set's members: for example, the fraction of members (or of lost mass) that is re-produced from inside the set within a time window. See M2.

### 5.4 Why a noisy finite soup can have zero strict organisations and still be organised

1. **Kinetic confinement.** The organisation is an infinite invariant subspace, and only a changing finite core of it is present. F&B: the set of species "is no longer closed under interaction". Mathis et al.: new expressions appear all the time and are "diluted out … because they are not being generated consistently". A strict closure test fails at every snapshot, even though the organisation persists.
2. **Members turn over while the organisation persists.** Mathis et al.: "the expressions are continually changing". F&B: the algebra "persists through a fluctuating, yet stably sustained, finite set of expressions". A fixed set of species is the wrong invariant. The invariant is a grammar or equivalence class, or a set of laws.
3. **Finite samples of reactions.** A network built from observed events holds only the reactions that happened to occur in the window. Rare products break closure, and rarely produced members break self-maintenance. Filisetti et al. (2011) handle this with "a temporal threshold that defines how long a specific reaction is kept in the reaction graph", so that cycles can be defined in an asynchronous stochastic simulation.
4. **Small copy numbers.** The fixed-point theorems are ODE results. Kreyssig et al. (2014) had to extend the theory to a "discrete chemical organization theory" to handle single molecules. B&Y §13.2.1 makes the same point about finite reactors (quoted above).
5. **Identity of molecules in continuous chemistries [synthesis].** If molecules are small neural networks, two products are "the same species" only under a chosen equivalence (quantisation, clustering or behaviour). Under exact identity, nearly every reaction produces a new species, closure is impossible, and F&B's convergence condition fails by construction. The number of organisations then depends on how fine the equivalence is. Report it as a curve over that tolerance, not as one number.
6. **How outflow is booked.** If decay or outflow is written as reactions, the flux condition asks for net production of every member, and rare members fail. If it is left out, trivial organisations appear. B&Y §13.1 shows that under the logical definition the empty set becomes an organisation, and that a decaying pair {A, B} passes the logical test but not the flux test.
7. **Infinite networks.** In an open-ended system, an infinite RAF can contain no finite RAF, and some infinite RAFs contain no irrRAF (Steel 2015). Any truncation of a constructive chemistry, such as a maximum length or a finite window, can hide the organisation.
8. **Structures that do not repeat exactly.** A metabolism that grows recursively "preserves in time its functionality" without any member recurring exactly (Kruszewski & Mikolov). Any criterion based on sets of species misses it.
9. **The count can also mislead the other way.** Strict organisations include trivial ones: single self-replicators, and the empty set if decay is modelled. A non-zero count therefore does not mean the soup is organised.

### 5.5 Detectors used in the artificial-chemistry literature

| Detector | Source | What it detects |
|---|---|---|
| Emergent laws (grammar plus algebra, via Knuth–Bendix completion), seeding sets, perturbation schedules (inject random objects, delete parts) | F&B 1994a,b, 1996 | Organisations as invariant subspaces; self-repair; resistance to foreign objects |
| Number of distinct expressions; survival when flooded with the identity; Jaccard similarity to self after a lag and to merged inputs; entropy; mean length | Mathis et al. 2024 | Stable organisations; robustness; outcome of mergers |
| Reactant-consumption rates as pointwise information against a chance model | Kruszewski & Mikolov 2021 | Metabolisms, including ones that grow recursively |
| High-order entropy (compression-based); number of distinct tokens | Agüera y Arcas et al. 2024 | State transitions when replicators appear |
| λ1 of the catalytic graph | Jain & Krishna 2002 | Whether an autocatalytic set exists, and how fragile it is |
| Temporal threshold on the reaction graph | Filisetti et al. 2011 | Cycles in stochastic asynchronous runs |
| Discrete chemical organisation theory | Kreyssig et al. 2014 | Stabilisation by single molecules |
| Degrees of closure and self-maintenance; cosine match to known organisations | B&Y §13.3–13.4 | How organised a soup is over time |
| Tendency to closure δK(V, l) | Montévil & Mossio 2015 | Boundaries between interacting closed systems |
| Information-theoretic individuality | Krakauer et al. 2020 | Units without a physical boundary |

### 5.6 A set of measures for Chemart [synthesis, each tied to a source]

All of these should be computed on equivalence classes of molecules and over sliding windows W.

- **M1. Empirical closure in a window.**
  - Definition: α_C(S, W) = (reaction events in W with all reactants in S and all products in S) ÷ (reaction events in W with all reactants in S).
  - Source: B&Y §13.4, weighted by how often each reaction actually occurred.
- **M2. Empirical self-maintenance in a window.**
  - Definition: the fraction of members s in S whose production in W is at least their loss in W (consumption plus dilution or removal), with a version weighted by mass.
  - Source: COT's flux condition, applied to realised fluxes. This replaces the printed α_S.
- **M3. Persistence.**
  - Definition: Jaccard similarity of the set to itself after a lag n (Mathis), together with an abundance-weighted similarity (because Jaccard ignores abundance) and a half-life of membership.
- **M4. Perturbation assays.**
  - (a) Survival when p% of the soup is replaced by a trivial molecule (Mathis).
  - (b) Rate at which injected random molecules become permanent members (F&B's "resistance to the addition").
  - (c) Delete a fraction and measure how fully the set regrows, and from how small a seed (F&B's self-repair and seeding sets).
  - These carry Chemart's network-level knockout tolerance over into the running dynamics.
- **M5. Convergence, which is where redundancy should show first.**
  - Definition: how the number of distinct products grows against the number of reaction events. Equivalently, how |S_n| grows under Chemart's existing scope expansion from the seed, compared with the number of collision histories.
  - Source: F&B's condition that normal forms grow polynomially while collision sequences grow exponentially.
  - For neural chemistries, compute it over the equivalence tolerance ε.
- **M6. λ1 of the observed catalytic graph in each window** (Jain & Krishna).
- **M7. Merger assay** (F&B; Mathis). Evolve two organisations separately, merge them, and record:
  - the outcome class (dominance, coexistence, mutual destruction, or fusion with glue);
  - the size of the glue, |G|/|O|, where G = O \ (A ∪ B);
  - mutual enablement: whether glue products feed both A and B;
  - an isolation test: return each component to a reactor on its own and check whether it still persists. Loss of autonomy is the signature of strong level 2 (U10).
- **M8. Static decomposition test.**
  - COT form: an organisation O is weak level 2 if O = A ⊔ B for two proper sub-organisations A and B that are not contained in one another, and the glue G is non-empty and not self-maintaining on its own.
  - RAF form: R = R1 ∪ R2 ∪ G, where R1 is a RAF only given products of R2, and the other way round. This is the meta-RAF of Hordijk, Steel & Kauffman.
  - Both can be computed on the lattice or subRAF poset that Chemart already builds. They are necessary conditions, not sufficient ones (§4.8, item 2).
- **M9. Individuation, in spatial variants only.**
  - The number of separated instances that share an organisation signature (U5).
  - McMullin's contact test.
  - Krakauer's information measure of individuality.
  - Montévil & Mossio's tendency to closure across spatial scales.

---

## 6. References

All entries below were checked in this session. "Bibliographic only" means the metadata (and usually the abstract) was confirmed but the full text was not read.

1. **Fontana W, Buss LW (1994a).** "The arrival of the fittest": Toward a theory of biological organization. *Bulletin of Mathematical Biology* 56(1):1–64. doi:10.1007/BF02458289. **VERIFIED** (Crossref; SFI working paper 93-09-055 abstract; full text read from the authors' copy at https://sites.santafe.edu/~walter/Papers/arrival.US.ps.gz).
2. **Fontana W, Buss LW (1994b).** What would be conserved if "the tape were played twice"? *PNAS* 91(2):757–761. doi:10.1073/pnas.91.2.757. **VERIFIED** (Crossref; PMC43028; full text read from https://sites.santafe.edu/~walter/Papers/tape.US.ps.gz).
3. **Fontana W, Buss LW (1996).** The barrier of objects: From dynamical systems to bounded organizations. In Casti J, Karlqvist A (eds), *Boundaries and Barriers*, pp. 56–116. Addison-Wesley. Also SFI working paper 96-05-035 and IIASA WP-96-027. https://sfi-edu.s3.amazonaws.com/sfi-edu/production/uploads/sfi-com/dev/uploads/filer/e7/12/e71200e6-546d-449c-adfa-2415e6f73220/96-05-035.pdf. **VERIFIED** (full text of the working paper read; book details from B&Y's bibliography and RePEc).
4. **Banzhaf W, Yamamoto L (2015).** *Artificial Chemistries*. MIT Press. §6.1, §6.3.2, §9.1, chapter 12, §13.1–13.5. **VERIFIED** (local text of the book).
5. **Dittrich P, Speroni di Fenizio P (2007).** Chemical organisation theory. *Bulletin of Mathematical Biology* 69:1199–1231. doi:10.1007/s11538-006-9130-8. Preprint: arXiv:q-bio/0501016 (2005). **VERIFIED** (Crossref; preprint full text read).
6. **Speroni di Fenizio P, Dittrich P (2002).** Artificial chemistry's global dynamics. Movement in the lattice of organisation. *The Journal of Three Dimensional Images* 16(4):160–163. **VERIFIED, bibliographic only** (Dittrich's publication list); content taken from B&Y §12.4.
7. **Peter S, Dittrich P (2011).** On the relation between organizations and limit sets in chemical reaction systems. *Advances in Complex Systems* 14(1):77–96. doi:10.1142/S0219525911002895. **VERIFIED, bibliographic only** (Crossref); the claim about attractors comes from the abstract as indexed by the search engine.
8. **Kreyssig P, Wozar C, Peter S, Veloz T, Ibrahim B, Dittrich P (2014).** Effects of small particle numbers on long-term behaviour in discrete biochemical systems. *Bioinformatics* 30(17):i475–i481. doi:10.1093/bioinformatics/btu453. **VERIFIED** (Crossref; abstract).
9. **Filisetti A, Graudenzi A, Serra R, Villani M, De Lucrezia D, Füchslin RM, Kauffman SA, Packard N, Poli I (2011).** A stochastic model of the emergence of autocatalytic cycles. *Journal of Systems Chemistry* 2:2. doi:10.1186/1759-2208-2-2. **VERIFIED** (Crossref; abstract).
10. **Hordijk W, Steel M (2004).** Detecting autocatalytic, self-sustaining sets in chemical reaction systems. *Journal of Theoretical Biology* 227(4):451–461. doi:10.1016/j.jtbi.2003.11.020. **VERIFIED** (Crossref; definition checked against the formal restatement in Steel 2015).
11. **Steel M (2015).** Self-sustaining autocatalytic networks within open-ended reaction systems. *Journal of Mathematical Chemistry* 53:1687–1701. doi:10.1007/s10910-015-0512-8; arXiv:1501.05731. **VERIFIED** (Crossref; full text read).
12. **Hordijk W, Steel M, Dittrich P (2018).** Autocatalytic sets and chemical organizations: modeling self-sustaining reaction networks at the origin of life. *New Journal of Physics* 20:015011. doi:10.1088/1367-2630/aa9fcd. **VERIFIED** (Crossref; abstract; publisher page).
13. **Hordijk W, Steel M, Kauffman S (2012).** The structure of autocatalytic sets: Evolvability, enablement, and emergence. *Acta Biotheoretica* 60:379–392. doi:10.1007/s10441-012-9165-1; arXiv:1205.0584. **VERIFIED** (Crossref; full text read).
14. **Blokhuis A, Lacoste D, Nghe P (2020).** Universal motifs and the diversity of autocatalytic systems. *PNAS* 117:25230–25236. doi:10.1073/pnas.2013527117. **VERIFIED** (Crossref; PMC7568248).
15. **Andersen JL, Flamm C, Merkle D, Stadler PF (2021).** Defining autocatalysis in chemical reaction networks. arXiv:2107.03086. **VERIFIED** (preprint; no journal version found on Crossref).
16. **Golnik R, Gatter T, Hordijk W, Stadler PF, Vassena N (2026).** Bridging two theoretical frameworks of autocatalysis: RAF sets and stoichiometric autocatalysis. arXiv:2605.25523. **VERIFIED** (preprint from May 2026, not peer-reviewed).
17. **Vasas V, Fernando C, Santos M, Kauffman S, Szathmáry E (2012).** Evolution before genes. *Biology Direct* 7:1. doi:10.1186/1745-6150-7-1. **VERIFIED** (Crossref; full text via Europe PMC).
18. **Vasas V, Szathmáry E, Santos M (2010).** Lack of evolvability in self-sustaining autocatalytic networks constrains metabolism-first scenarios for the origin of life. *PNAS* 107(4):1470–1475. doi:10.1073/pnas.0912628107. **VERIFIED** (Crossref, where the title carries the typo "constraints"; abstract).
19. **Segré D, Ben-Eli D, Lancet D (2000).** Compositional genomes: Prebiotic information transfer in mutually catalytic noncovalent assemblies. *PNAS* 97(8):4112–4117. doi:10.1073/pnas.97.8.4112. **VERIFIED**.
20. **Jain S, Krishna S (2001).** A model for the emergence of cooperation, interdependence, and structure in evolving networks. *PNAS* 98(2):543–547. doi:10.1073/pnas.98.2.543. Also **Jain S, Krishna S (2002).** Crashes, recoveries, and "core shifts" in a model of evolving networks. *Physical Review E* 65:026103. doi:10.1103/PhysRevE.65.026103. **VERIFIED** (both via Crossref and abstracts).
21. **Szathmáry E (2000).** The evolution of replicators. *Philosophical Transactions of the Royal Society B* 355:1669–1676. doi:10.1098/rstb.2000.0730. **VERIFIED**.
22. **Szathmáry E (2006).** The origin of replicators and reproducers. *Philosophical Transactions of the Royal Society B* 361:1761–1776. doi:10.1098/rstb.2006.1912. **VERIFIED** (Crossref; abstract).
23. **Szathmáry E (2015).** Toward major evolutionary transitions theory 2.0. *PNAS* 112(33):10104–10111. doi:10.1073/pnas.1421398112. **VERIFIED**.
24. **Varela FG, Maturana HR, Uribe R (1974).** Autopoiesis: the organization of living systems, its characterization and a model. *BioSystems* 5(4):187–196. doi:10.1016/0303-2647(74)90031-8. **VERIFIED** (Crossref; abstract). The definition in §2 is quoted from secondary sources because the original PDF was not accessible.
25. **McMullin B (2004).** Thirty years of computational autopoiesis: A review. *Artificial Life* 10(3):277–295. doi:10.1162/1064546041255548. **VERIFIED** (Crossref; full text of the author's copy at https://www.eeng.dcu.ie/~alife/bmcm-alj-2004/html-single/). The individuation test comes from McMullin B (2000), Remarks on autocatalysis and autopoiesis, *Annals of the New York Academy of Sciences* 901:163–174, doi:10.1111/j.1749-6632.2000.tb06276.x (Crossref-verified, read through the 2004 review).
26. **Montévil M, Mossio M (2015).** Biological organisation as closure of constraints. *Journal of Theoretical Biology* 372:179–191. doi:10.1016/j.jtbi.2015.02.029. **VERIFIED** (full text at montevil.org). Related book: **Moreno A, Mossio M (2015).** *Biological Autonomy: A Philosophical and Theoretical Enquiry*. Springer. doi:10.1007/978-94-017-9837-2. **VERIFIED, bibliographic and table of contents only** (includes a chapter "Organisms and levels of autonomy").
27. **Letelier JC, Soto-Andrade J, Guíñez Abarzúa F, Cornish-Bowden A, Cárdenas ML (2006).** Organizational invariance and metabolic closure: Analysis in terms of (M,R) systems. *Journal of Theoretical Biology* 238(4):949–961. doi:10.1016/j.jtbi.2005.07.007. **VERIFIED**.
28. **Jaramillo S, Honorato-Zimmer R, Pereira U, Contreras D, Reynaert B, Hernández V, Soto-Andrade J, Cárdenas ML, Cornish-Bowden A, Letelier JC (2010).** (M,R) systems and RAF sets: Common ideas, tools and projections. *Proceedings of ALIFE XII* (Odense), MIT Press, pp. 94–100. https://www.research.ed.ac.uk/en/publications/mr-systems-and-raf-sets-common-ideas-tools-and-projections/. **VERIFIED** (Edinburgh repository record and abstract; that record lists the authors in a different order).
29. **Gánti T (1975).** Organization of chemical reactions into dividing and metabolizing units: The chemotons. *BioSystems* 7(1):15–21. doi:10.1016/0303-2647(75)90038-6. **VERIFIED, bibliographic only**; the description comes from B&Y §6.1 and Szathmáry 2006.
30. **Pascal R, Pross A, Sutherland JD (2013).** Towards an evolutionary theory of the origin of life based on kinetics and thermodynamics. *Open Biology* 3:130156. doi:10.1098/rsob.130156. **VERIFIED**.
31. **Krakauer D, Bertschinger N, Olbrich E, Flack JC, Ay N (2020).** The information theory of individuality. *Theory in Biosciences* 139:209–223. doi:10.1007/s12064-020-00313-7. **VERIFIED**.
32. **Peng Z, Plum AM, Gagrani P, Baum DA (2020).** An ecological framework for the analysis of prebiotic chemical reaction networks. *Journal of Theoretical Biology* 507:110451. doi:10.1016/j.jtbi.2020.110451. **VERIFIED**.
33. **Baum DA, Peng Z, Dolson E, Smith E, Plum AM, Gagrani P (2023).** The ecology–evolution continuum and the origin of life. *Journal of the Royal Society Interface* 20:20230346. doi:10.1098/rsif.2023.0346. **VERIFIED**. This paper argues that meta-ecosystems of autocatalytic chemical ecosystems can evolve before bounded individuals exist, and that "adaptive evolution can explain the emergence of self-bounded units".
34. **Mathis C, Patel D, Weimer W, Forrest S (2024).** Self-organization in computation and chemistry: Return to AlChemy. *Chaos* 34(9):093142. doi:10.1063/5.0207358; arXiv:2408.12137. **VERIFIED** (Crossref; full text of arXiv v2 read).
35. **Vimal D, Mathis C, Weimer W, Forrest S (2025).** Prebiotic functional programs: Endogenous selection in an artificial chemistry. arXiv:2509.03534. **VERIFIED** (preprint; full text read).
36. **Agüera y Arcas B, Alakuijala J, Evans J, Laurie B, Mordvintsev A, Niklasson E, Randazzo E, Versari L (2024).** Computational life: How well-formed, self-replicating programs emerge from simple interaction. arXiv:2406.19108. **VERIFIED** (full text read).
37. **Kruszewski G, Mikolov T (2021).** Emergence of self-reproducing metabolisms as recursive algorithms in an artificial chemistry. *Artificial Life* 27(3–4):277–299 (2021 volume, published online 16 March 2022). doi:10.1162/artl_a_00355; arXiv:2103.08245. **VERIFIED** (Crossref; full text read).
38. **Hickinbotham SJ, Stepney S, Hogeweg P (2021).** Nothing in evolution makes sense except in the light of parasitism: Evolution of complex replication strategies. *Royal Society Open Science* 8:210441. doi:10.1098/rsos.210441. **VERIFIED** (abstract via Europe PMC).
39. **Boerlijst MC, Hogeweg P (1991).** Spiral wave structure in pre-biotic evolution: Hypercycles stable against parasites. *Physica D* 48:17–28. doi:10.1016/0167-2789(91)90049-F. **VERIFIED, bibliographic only**.
40. **Takeuchi N, Hogeweg P (2009).** Multilevel selection in models of prebiotic evolution II: A direct comparison of compartmentalization and spatial self-organization. *PLoS Computational Biology* 5(10):e1000542. doi:10.1371/journal.pcbi.1000542. **VERIFIED**.
41. **Speroni di Fenizio P, Dittrich P, Banzhaf W (2001).** Spontaneous formation of proto-cells in an universal artificial chemistry on a planar graph. In Kelemen J, Sosík P (eds), *ECAL 2001*, LNAI 2159, pp. 206–215. Springer. doi:10.1007/3-540-44811-X_22. **VERIFIED** (Crossref; abstract and introduction read from the proceedings text).
42. **Ono N, Ikegami T (2000).** Self-maintenance and self-reproduction in an abstract cell model. *Journal of Theoretical Biology* 206(2):243–253. doi:10.1006/jtbi.2000.2121. **VERIFIED**.
43. **Goldsby HJ, Dornhaus A, Kerr B, Ofria C (2012).** Task-switching costs promote the evolution of division of labor and shifts in individuality. *PNAS* 109(34):13686–13691. doi:10.1073/pnas.1202233109. **VERIFIED** (abstract and methods via PMC3427090).

*Not included:* Tierra and Coreworld (outside this subtopic's needs), Hutton 2007 (cells built by design, not emergent), and Kreyssig et al. 2012 (a marginal link between cycles and transitions between organisations).
