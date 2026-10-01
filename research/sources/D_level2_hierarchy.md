# Level 2: what an "organisational hierarchy" or "higher-order entity" is, and how to detect and measure one

*Literature report D for Chemart. Written 2026-09-23. All 72 references were checked in this session (§8 says how each was checked). Quotes are verbatim from the source unless marked as a paraphrase. My own proposals, which are not claims from the literature, are marked **[proposal]**.*

---

## Why this report exists

Chemart has two jobs: to be a morphospace of artificial chemistries (ACs), and to test the premise behind the redundancy hypothesis. That premise is that algorithmic chemistries such as AlChemy and BFF hit a ceiling, and that more redundant molecules might lift it. Marco frames the ceiling as a three-level picture taken from Fontana & Buss: molecules (level 0), then networks of molecules (level 1), then networks of networks (level 2).

The original definitions are short. Fontana & Buss's abstract (SFI working paper 93-09-055, spelling normalised) says: *"Level 0 is defined by self-copying objects or simple ensembles of copying objects. Level 1 denotes a new object class, whose objects are self-maintaining organizations made of Level 0 objects, and Level 2 is defined by self-maintaining metaorganizations composed of Level 1 organizations."* Their PNAS companion lists as a generic outcome that *"self-maintaining organizations, once established, can combine into higher-order self-maintaining organizations"* (Fontana & Buss 1994b).

Thirty years later, Mathis et al. (2024) re-ran AlChemy. They found that level-1 organizations are common, *"but that these stable organizations cannot be easily combined into higher order entities."* They merged pairs of level-1 soups (455 pairs) and scored each outcome as dominance, coexistence (mean Jaccard similarity > 0.1 to both inputs) or mutual destruction. Coexistence *"rarely occurs for organizations evolved in different simulations."* They also give the most operational description of level 2 found anywhere: *"L2 organizations are characterized by the existence of two or more L1 organizations and additional expressions (called 'glue' in the original work). The 'glue' expressions could not exist without at least one of the L1 organizations and are produced by composing functions from different organizations."* They add that *"deciding whether an organization can be split into two distinct stable organizations is a difficult problem to solve without resorting to trial and error."*

So level 2 is the target, but nothing in Chemart can see one. Chemart's hierarchy-flavoured measures are flow hierarchy, modularity, NODF nestedness, bow-tie, and the list of chemical organisations. All of them describe one reaction graph. This report covers four things: what "hierarchy" and "higher-order entity" mean in the literatures that have thought hardest about them, what each literature counted as a new level, what it measured, and which of those measures could run on Chemart data.

---

## 0. The short version (eight points)

1. **"Hierarchy" names at least four different things, and only two of them together make Marco's level 2.** The four are: order/control/flow (who drives whom), compositional nesting (parts within wholes), individuality or levels of selection (wholes made of former individuals), and Salthe's specification hierarchy. Level 2 is **compositional + individuated**. Every flow-hierarchy measure is blind to it, including Chemart's flow hierarchy and, from the literature, global reaching centrality, trophic coherence and treeness/feedforwardness/orderability. Nested-module methods (Ravasz, the hierarchical random graph, Sales-Pardo, the nested SBM) see nested *graph structure*, not self-maintaining wholes (§1, §6).
2. **The AC-native definition of level 2 is Fontana & Buss's "self-maintaining metaorganization composed of Level 1 organizations", and Mathis et al.'s "glue" makes it operational.** Two formal tools can generate candidates. In chemical organisation theory (COT), organisations form a lattice for catalytic flow systems, and AlChemy is one. RAF theory has the subRAF poset and a polynomial-time test for whether a RAF is the union of two proper subRAFs (Hordijk, Steel & Kauffman 2012). Hordijk et al. also name the target, *"an autocatalytic set of autocatalytic sets"*, but they do not give a criterion for detecting one (§6.3).
3. **McShea (2001) gives the best structural scoring rubric.** It has levels of nestedness × degree of individuation, with sublevels a (monomorphic aggregate), b (differentiated) and c (with intermediate-level parts), each with operational criteria. It also shows that hierarchy is best tracked as a *trend in the maximum* (§2.6).
4. **Evolutionary individuality can only be measured when collectives reproduce.** This covers Maynard Smith & Szathmáry, Michod and Okasha's MLS2. The measurable quantities are a Price-equation partition, Okasha's stages of fitness decoupling, and Michod's covariance effect. A single well-stirred reactor has no collective reproduction. You need compartments or space (the stochastic corrector, Hogeweg's waves and vesicles, P-system division) or ecological scaffolding: patches plus dispersal (Black, Bourrat & Rainey 2020) (§3, §7).
5. **The most promising general detector is closure of a candidate coarse-graining.** The route runs from informational closure (Bertschinger et al. 2006; Pfante et al. 2014) to informational, causal and computational closure tied to Markov-chain lumpability (Rosas et al. 2024), plus Krakauer et al.'s (2020) individuality measures. Each needs a trajectory and a *candidate* partition. They are feasible only after the soup has been compressed to a few macro-variables. The space of coarse-grainings grows super-exponentially, so candidates must come from structure (COT/RAF). Krakauer's measures never decrease as the candidate system grows, so a size penalty is needed (§5).
6. **Watch the redundancy trap.** Rosas et al.'s (2020) practical emergence criteria Ψ and Δ are whole-minus-sum quantities. In their words, *"redundancy will drive Ψ and Δ more negative"*, and the criteria *"double-count redundancy up to n times"*. They are therefore biased against exactly the redundant chemistries the redundancy hypothesis wants to build. The effective-information approach behaves the other way round: it finds *more* causal emergence in degenerate (convergent) biological networks (Klein & Hoel 2020) (§5.6–5.7).
7. **The ALife "dynamical hierarchies" programme learned early that loose definitions are trivially satisfiable.** Critics built "higher-order" structures from simpler parts (Gross & McMullin 2001) and built an infinite self-assembling hierarchy whose new properties were trivial (Dorin & McCormack 2002). Rasmussen et al.'s Ansatz, the conjecture that more levels need more complex primitives, is close to Marco's premise, and it was contested (§4.6).
8. **Wimsatt's robustness gives the acceptance rule.** Call something a level only if independent detectors agree: structural, dynamical-informational, perturbational and, where possible, selective. Report each as its own morphospace axis rather than folding them into one number (§2.2, §7).

---

## 1. The landscape: what people mean by "hierarchy"

| Sense | Question it answers | What a "level" is | Canonical sources | Typical measures | Can it see level 2? |
|---|---|---|---|---|---|
| **(a) Order / control / flow** | Who drives or constrains whom? | A rank in a directed order | Simon's "formal hierarchy"; Pattee's control hierarchy; Luo & Magee 2011; Mones et al. 2012; Johnson et al. 2014; MacKay et al. 2020; Corominas-Murtra et al. 2013 | flow hierarchy, global reaching centrality, trophic levels and coherence, treeness/feedforwardness/orderability | **No.** It measures directionality of influence, not wholes |
| **(b) Compositional / nested ("scalar")** | What is part of what? | A whole made of parts | Simon 1962; Salthe (scalar hierarchy); McShea 2001 (nestedness); Ravasz et al. 2002; Clauset et al. 2008; Sales-Pardo et al. 2007; Peixoto 2014; P-system membranes; the COT organisation lattice | nested modules, dendrograms, containment trees, lattices of organisations | **Half.** It gives nesting without individuation, except where the parts are defined as self-maintaining (COT, RAF) |
| **(c) Individuals / levels of selection** | Which units are individuals, and which units does selection act on? | An entity made of former individuals that is now itself an individual | Maynard Smith & Szathmáry 1995; Szathmáry 2015; Michod; Okasha 2006; Queller & Strassmann 2009; Godfrey-Smith 2009 | MLS2 Price partition, fitness decoupling, covariance effect, the Darwinian space for reproducers | **Yes, if collectives reproduce** |
| **(d) Specification ("integrative levels")** | Which realm is a refinement of which? | A subclass subsumed by a more general class, such as {physical {chemical {biological}}} | Salthe | none quantitative | No |
| *(x) "Networks of networks" = interdependent / multilayer networks* | How do coupled networks fail or percolate together? | A network layer | Gao, Buldyrev, Stanley & Havlin 2012 | cascades, mutual percolation | **Unrelated: a search trap.** Gao et al. review a percolation framework for networks that interact with and depend on other networks. The phrase does not mean networks whose nodes are networks. |

**How these senses relate.**

- **Simon (1962)** explicitly *widened* "hierarchy" beyond (a). By a hierarchic system he means *"a system that is composed of interrelated subsystems, each of the latter being, in turn, hierarchic in structure until we reach some lowest level of elementary subsystem"* (p. 468). He includes *"systems in which there is no relation of subordination among subsystems"* and reserves *"'formal hierarchy'"* for the boss-and-subordinates sense (p. 468).
- **Eronen & Brooks (SEP, 2018/2023)** name the classic split as *"levels of composition and levels of control (Simon 1962; Pattee 1973)"*. The former is *"the nested compositionality typically identified with levels of organization"*. The latter is *"the idea that higher levels impose constraints on the processes at lower levels, for example by limiting the degrees of freedom of the system at a lower level (Pattee 1973: 85)"*.
- **McShea (2001)** separates (b) from (c) sharply. Maynard Smith & Szathmáry's criterion (lower-level entities lose independent replication) means *"their scheme constitutes what might be called a process hierarchy ... or a control hierarchy, rather than a structural hierarchy; in other words, their concern is with a fundamentally different phenomenon"* (p. 408).
- **Fontana & Buss's level 2 is (b) plus individuation, and does not require (c).** An organization in their sense *"is self-maintaining, and is characterized by (i) boundaries established by the invariances, (ii) strong self-repair capabilities responsible for a robustness to perturbation, and (iii) a center, defined as the smallest set kinetically persistent and self-maintaining generator set of the algebra."* They stress that organizations arise *"without appeal to natural selection"*, and they draw the analogy *"from self-replication to self-maintaining procaryotic organizations to ultimately yield self-maintaining eucaryotic organizations."* Their level 2 is therefore the analogue of McShea's step from prokaryote to eukaryote (§2.6), and of Szathmáry's *egalitarian* transitions (§3.2).

**Working definition used in the rest of this report [proposal, assembled from the sources above].** A *level-2 entity* in an AC is a set of ≥ 2 level-1 organizations with four properties:
- **nested:** each part is itself self-maintaining;
- **integrated:** the whole contains glue that no part produces alone, and the parts depend on each other;
- **individuated:** the whole is closed, self-repairing and bounded as a dynamical unit;
- **more than an aggregate:** it has its own dynamics at a slower timescale.

Sense (c) is an optional fifth property, available when composites reproduce as units.

---

## 2. Conceptual foundations

### 2.1 Simon (1962): near-decomposability

**What it is.** Hierarchies are *nearly decomposable*: systems *"in which the interactions among the subsystems are weak, but not negligible"* (p. 474). Two propositions follow: *"(a) in a nearly decomposable system, the short-run behavior of each of the component subsystems is approximately independent of the short-run behavior of the other components; (b) in the long run, the behavior of any one of the components depends in only an aggregate way on the behavior of the other components"* (p. 474).

**How it works.** Order the interaction matrix so that *"all its large elements lie inside a string of square submatrices along the main diagonal ... We may take some small number, ε, as the upper bound of the extradiagonal elements. We shall call a matrix having these properties a nearly decomposable matrix"* (p. 475). Simon's own physico-chemical example sets the level by ε: *"If we select an epsilon just a little smaller than the magnitude of a covalent bond, the system will decompose into subsystems—the constituent molecules"* (p. 475). His evolutionary argument (the watchmaker parable) is that *"complex systems will evolve from simple systems much more rapidly if there are stable intermediate forms than if there are not. The resulting complex forms in the former case will be hierarchic"* (p. 473).

**What counts as a level.** A diagonal block with strong internal and weak external coupling. Each level is defined relative to a threshold ε, and each level-crossing is a separation of timescales.

**On AC data.**
- Inputs: the Jacobian of the mass-action ODE at a steady state, or an empirical species-by-species flux matrix counted from reaction events in a trajectory window.
- Measurable quantities: the ratio of off-block to on-block coupling for a candidate partition, and the gap between within-block relaxation times and between-block times.
- Cost: cheap for ODE chemistries, moderate for stochastic soups (event counting per window).

**Limits.** A nearly decomposable block is a *module*, not an individual: nothing requires it to be self-maintaining. Watson & Pollack (2005), in the 2005 dynamical-hierarchies issue, stress that *"modularity need not imply that intermodule dependences are weak or unimportant"*. Simon's long-run aggregate interaction can be strong. A level-2 composite should look Simon-like: strong coupling inside each level-1 part, weaker but essential aggregate coupling between parts.

### 2.2 Wimsatt: levels as local maxima of regularity and predictability

**What it is.** As quoted in the SEP: *"[l]evels of organization can be thought of as local maxima of regularity and predictability in the phase space of alternative modes of organization of matter"* (Wimsatt 1976a: 209; the SEP notes this is *"the closest that [Wimsatt] will come to a definition"*). Levels are *"constituted by families of entities usually of comparable size and dynamical properties, which characteristically interact primarily with one another"* (Wimsatt 1994 [2007: 204]). In the SEP's gloss: *"if we plot regularity and predictability against (size) scale, then levels of organization will appear as peaks in the plot."*

**Robustness.** Levels and their entities *"should be detectable, measurable, derivable, definable, and so on, in a variety of independent ways"* (SEP, paraphrasing Wimsatt 1981/1994).

**What this means for Chemart [proposal].**
- The "peak of predictability against scale" is exactly what the closure and causal-emergence measures in §5 formalise: predictability of a macro-description as a function of coarse-graining.
- Robustness is the acceptance rule for a level-2 claim. Structural, informational and perturbational detectors should agree on the same composite.

**Limits.** The SEP asks which regularities count and how one would count them. Wimsatt is a guide to what to measure, not a measure.

### 2.3 Eronen & Brooks, SEP "Levels of Organization in Biology"

**What it is.** A survey of three accounts: the layer-cake account (Oppenheim & Putnam), levels of mechanisms (Craver), and Wimsatt's local-maxima account. It also covers deflationary critiques. The key critique is Potochnik & McGill: levels demand *"not only the ubiquity, but also the uniformity, of part-whole composition"*, and nature does not provide that.

**Where it bears on selection.** For levels of selection, the SEP reports Okasha's requirement (following McShea and Sober & Wilson) that entities at a level have fitness-affecting interactions. It also reports his requirement that they be *"homologous with organisms in a free-living state, either extant or extinct."*

**What this means for Chemart.** Don't look for a global ladder of layers across the whole soup. Level-2 entities will be *local* wholes, which may coexist with free level-1 organizations and plain molecules in the same reactor.

### 2.4 Salthe: scalar versus specification hierarchies

**What it is.** From Salthe's own summary (2001), formalised in Salthe (2012):
- *"The scalar hierarchy is one of parts nested within wholes ... where [higher level [focal level [lower level]]]"*.
- *"The specification hierarchy is one of classes and subclasses, as e.g., {material world {biological world {social world}}}"*.

**Criteria.** *"In the scalar hierarchy components at different levels differ in size roughly by orders of magnitude"*. Levels in the specification hierarchy *"mark the qualitative differences of different realms of being"*.

**Dynamics.** Because of the order-of-magnitude differences, *"dynamics at different levels do not directly interact or exchange energy, but transact by way of mutual constraint"*. Higher levels act as slowly varying boundary conditions. Salthe calls this *"screening off"*.

**Two points that matter for level 2.**
- *"If the parts are functional in some given analysis, they are referred to as components, if not they are constituents."* In a level-2 entity, the level-1 organizations should be *components*.
- On how new levels form: *"The actual process of formation of a level would involve the cohesion of entities out of lower level units guided by higher level boundary conditions. This process is little understood since this hierarchy has largely been used for synchronic analyses."*

**On AC data [proposal].** Salthe's criterion becomes a timescale test. The macro-variables of a level-2 candidate, such as the relative abundances of its level-1 parts, should vary much more slowly than those parts' internal species abundances. The internal abundances should in turn vary more slowly than single-molecule turnover. Autocorrelation times computed from trajectories are cheap.

### 2.5 Pattee: control hierarchy

Pattee (1973, in the volume *Hierarchy Theory*) is the classic source for *levels of control*: higher levels as constraints that limit lower-level degrees of freedom. The SEP (§1) summarises him, and I did not read the chapter itself. The distinction matters because a *control* level can exist without *compositional* nesting (a regulator), and the reverse also holds.

### 2.6 McShea (2001): a scale of hierarchical structure

McShea gives the most operational structural definition anywhere, and he designed it to be scored on fossils. His abstract describes it as *"a higher-resolution scale ... in which hierarchical structure is decomposed into levels and sublevels, with levels reflecting number of layers of nestedness, and sublevels reflecting degree of individuation at the highest level."*

**Nestedness.**
- *"Nestedness refers to physical containment or inclusion: a higher-level individual contains entities from the next lower level, and includes them as parts"* (p. 408).
- Extra requirements: lower-level entities must *"be bounded, and also ... be spatially aggregated and attached to each other, an indication that they interact in some way and that their behaviors are correlated"* (p. 408).
- *"Implicit in this numbering scheme is a requirement that lower-level entities must be homologous with organisms in a free-living state, either extant or extinct"* (p. 408).
- Numbering: level 1 = prokaryotic cell; 2 = aggregates of prokaryotic cells (eukaryotic cell); 3 = aggregates of level-2 organisms (multicellular); 4 = aggregates of level-3 organisms (colonial individuals).

**Individuation.** *"the degree to which the highest-level entity in a nested sequence constitutes a unified whole"* (p. 409). It has three criteria: *"(1) connectedness (a generalization of what Beklemishev calls 'integration'), (2) differentiation ('polymorphism'), and (3) the presence of intermediate-level parts ('cormidia')"* (p. 409).

**The sublevels** (Fig. 2 caption):

| Sublevel | Definition (verbatim) |
|---|---|
| **a** | monomorphic aggregates: *"aggregates of two or more lower-level entities that remain reliably attached over some significant portion of their existence"* |
| **b** | differentiated aggregates: *"attached aggregates with two or more different morphological types"* |
| **c** | differentiated aggregates with intermediate-level parts: *"parts may consist either of a subgroup of two or more lower-level entities or a single lower-level entity that is hypertrophied or elaborated in various ways"* |

The full scale runs 1c (solitary prokaryote), 2a, 2b, 2c (solitary eukaryotic cell), 3a, 3b, 3c (solitary metazoan), 4a, 4b, 4c (individuated metazoan colony).

**Caveats McShea gives.**
- Overconnectedness: *"extreme connectedness may tend to undermine individuation and thus to undermine hierarchical structure"* (p. 410).
- The focus is structural, not the process or control hierarchy (§1).
- The trend he documents is a trend in the **maximum**: *"an increase in the degree of hierarchical structure present in the hierarchically deepest organism on Earth."*

**Mapping to ACs [proposal].**

| McShea | AC analogue |
|---|---|
| lower-level entity | a level-1 organization (a COT organisation or irreducible RAF that persists under perturbation) |
| "homologous with free-living" | each part has been observed self-maintaining on its own (the input soups in Mathis et al.'s merge experiment are exactly this) |
| "attached over a significant portion of existence" | co-persistence over a window much longer than the parts' turnover time |
| a (monomorphic) | several copies of the same organization held together. This needs compartments, because in one well-stirred reactor, copies merge |
| b (differentiated) | ≥ 2 distinct level-1 organizations plus glue: Fontana & Buss level 2 |
| c (intermediate parts) | sub-composites inside a larger composite: the start of level 3 |
| overconnectedness | every part reacts with every other, so the parts dissolve into one undifferentiated organization (back to level 1) |

The output is an ordinal (level, sublevel) score per composite. Its **maximum over a run** is a morphospace axis, following McShea's "trend in the maximum". It needs the decomposition tools in §6.3.

---

## 3. Major transitions and individuality: which theories give a number?

### 3.1 Maynard Smith & Szathmáry (1995)

**What it is.** A two-fold characterisation, as quoted by Okasha (2022):
- a *"change in the way that information is stored and transmitted"*;
- *"entities that were capable of independent replication before the transition can replicate only as part of a larger whole after it"* (p. 8).

Okasha (2022) recommends the second, hierarchical criterion: a major transition is *"the formation of a new higher-level biological unit from a group of lower-level units"*.

**On AC data [proposal].** "Replicate only as part of a larger whole" becomes a **loss-of-autonomy test**. After the composite forms, does level-1 organization A still self-maintain when B is removed? Mathis et al.'s merge experiments measure the forward direction (can A and B coexist?). The reverse test, knocking out a part of an established composite, has not been done in AlChemy as far as I found.

### 3.2 Szathmáry (2015): "Toward major evolutionary transitions theory 2.0"

**Egalitarian and fraternal transitions** (from the abstract): *"the concept of fraternal and egalitarian transitions (lower-level units like and unlike, respectively)"*. Following Queller: *"In the first, like units join or remain joined, reaping the first benefits from the economy of scale, and then evolving division of labor by differentiation. In the second, unlike units come together, complementing their functions in a higher unit ... The main control of conflicts is ensured by kinship and fairness in reproduction for the fraternal and egalitarian transitions, respectively."*

**Diagnostic of success:** *"After the transition, units show strong cooperation and very limited realized conflict."* Phases: *"origin, maintenance, and transformation (i.e., further evolution) of the higher level units."*

**Reproduction versus replication:** *"During transitions, new units of reproduction emerge, and establishment of such units requires high fidelity of reproduction (as opposed to mere replication)."* Also: *"In simple forms, reproduction is compositional (only numbers of different particle types matter)"*.

**Egalitarian transitions have no stage 3:** *"particle fitness values cannot go down to zero, but they need to be tightly controlled through the mediation of conflicts (reproductive leveling) ... There is no stage 3 for egalitarian transitions because no reproductive division of labor can exist."*

**Protocells, MLS1 versus MLS2:** *"All passive models of compartmentation are examples of multilevel selection models of type 1 (MLS1) where the focal units are still the individual replicators rather than the groups"*. By contrast, the stochastic corrector *"is a clear case of multilevel selection of the second type (MLS2) where the focal units are groups (or collectives)"*.

**A warning that applies directly to ACs:** *"there is a notoriously recurring error in the literature equating any collectively autocatalytic network with hypercycles ... Cross-catalytic peptides or anabolic ligases are collective autocatalysts but their members are not cooperators in the evolutionary sense."*

**What this means for Chemart [proposal].**
- Fontana & Buss level 2, two *unlike* level-1 organizations plus glue, is an **egalitarian** transition.
- Its known failure modes are conflict between unlike partners, which matches Mathis et al.'s "dominance" and "mutual destruction".
- The known stabiliser is fair or synchronised reproduction of the partners *as a unit*. A single well-stirred reactor has no such mechanism, which is a structural reason to expect a level-2 ceiling there regardless of how redundant the molecules are.

### 3.3 Michod: fitness decoupling and the covariance effect

**Concept** (Michod & Nedelcu 2003, abstract): the problem is *"how a group of individuals becomes a new kind of individual, possessing the property of heritable variation in fitness at the new level of organization"*. *"Only cooperation transfers fitness from lower levels (costs to group members) to higher levels (benefits to the group)."* Szathmáry (2015) quotes them (p. 66): *"as the evolutionary transition proceeds, group fitness becomes decoupled from the fitness of its lower-level components."*

**A measurable quantity** (Michod 2006, PNAS). The group's viability `V` and fecundity `B` are averages of the cells' efforts `v` and `b`, and group fitness is `W = V·B`. Then:

```
W − w̄ = −Cov(v, b)     # w̄ = mean cell fitness (mean of v·b over cells)
```

In Michod's words: *"the fitness of the cell group, W (taken as the product of V and B), is greater than the average fitness of member cells, w̄, by an amount equal to the negative covariance of the fitness components at the cell level"*. Trade-offs, where the covariance is negative, make the group more than the average of its parts.

**On AC data.** This needs per-part fitness components (for example, a level-1 organization's persistence as "viability" and its contribution to founding new compartments as "fecundity"). It is therefore feasible only where composites reproduce (§7, stage E).

### 3.4 Okasha: MLS1 versus MLS2, and the Price-equation partition

**The terms.** Terminology from Damuth & Heisler (1988), made central by Okasha (2006):
- **collective fitness₁** = the average fitness of the collective's particles (MLS1: particles are the focal units);
- **collective fitness₂** = the number of offspring *collectives* (MLS2: collectives are the focal units).

**Okasha's stages** (Okasha 2005, p. 1023, quoted by Szathmáry 2015): *"(Stage 1) Collective fitness defined as average particle fitness (cooperation spreads among particles). (Stage 2) Collective fitness not defined as average particle fitness, but still proportional to average particle fitness (collectives start to emerge as entities in their own right). (Stage 3.) Collective fitness neither defined as nor proportional to average particle fitness (collectives have fully emerged; fitnesses are decoupled)."*

**The accounting identity** (Price 1972, eq. A15, for groups of sizes `n` and offspring sizes `n′`):

```
ΔP = cov_n(s̃, p) + ave_n′(Δp)
```

Here `P` is the population-wide frequency, `p` the within-group frequency and `s̃` the relative group growth factor. The first term is change due to *between-group* selection and the second is change *within* groups. Price: this *"does not depend upon any assumptions about mechanisms of heredity or anything else of that sort, but holds because it is an identity"*. Because the second term is itself a change in a mean, the same identity can be applied inside each group, which gives a recursive multilevel partition. Okasha (2006) evaluates multilevel selection with both the Price equation and contextual analysis.

**What is measurable.**
- MLS1 partition: needs groups (compartments, patches or spatial clusters) and particle counts across one "generation".
- MLS2: needs a parent–offspring map among collectives.
- Okasha stage diagnostic [proposal]: regress collective fitness₂ on mean particle fitness across collectives. Identity means stage 1, proportional means stage 2, uncorrelated means stage 3.

**Limit.** Price partitions are statistical, not causal (Okasha's reason for considering contextual analysis). Egalitarian composites cannot reach stage 3 (Szathmáry 2015).

### 3.5 Queller & Strassmann (2009): organismality as high cooperation and low conflict

**Definition** (abstract): *"all organisms originated from groups of simpler units that now show high cooperation among the parts and are nearly free of conflicts. We suggest that this near-unanimous cooperation be taken as the defining trait of organisms."*

**What is not required.** Traits they list as *not* essential: *"physical contiguity, indivisibility, clonality or high relatedness, development from a single cell, short-term and long-term genetic cotransmission, germ-soma separation and membership in the same species."*

**On AC data [proposal].** Two axes:
- *cooperation*: the fraction of cross-part reactions that feed or catalyse the other part;
- *conflict*: whether cheaters invade. A cheater is a species that receives catalysis from the composite but supplies none, which is Boerlijst & Hogeweg's parasite. This can be tested by injecting such species.

This is attractive for Chemart because it needs no reproduction and no membrane.

### 3.6 Godfrey-Smith (2009): Darwinian spaces

**Three parameters.** For collective reproducers, as stated in his own chapter summarising the book (Godfrey-Smith 2012 preprint; the figure is from the 2009 book):
- **B**, bottleneck: *"a narrowing that marks the divide between generations"*;
- **G**, germ line: *"the degree of reproductive specialization within a collective"*;
- **I**, integration: *"a general division of labor (aside from that in G), the mutual dependence of parts, and the maintenance of a boundary between a collective and what is outside it"*.

*"Clear cases of collective reproduction are associated with 'high scores' on all these features. Marginal cases are associated with low scores."* He also notes that *"The evolution of new individuals partly 'de-Darwinizes' the old ones that make them up."*

**On AC data.** B and G need collective reproduction. **I is scoreable in a single reactor** (mutual dependence plus boundary, where the boundary is the closure boundary of an organization).

### Summary: which theories give a measurable quantity?

| Theory | Quantity | Inputs | Needs collective reproduction? |
|---|---|---|---|
| Maynard Smith & Szathmáry | loss of independent replication | knockout experiments | no (for the loss-of-autonomy test) |
| Szathmáry 2015 | cooperation/conflict; fidelity of composite reproduction | parasite injection; parent–offspring composition | partly |
| Michod | `W − w̄ = −Cov(v,b)`; decoupling | per-part fitness components | yes |
| Okasha / Price | between/within partition; stages 1–3 | groups, offspring counts | MLS1: groups only; MLS2: yes |
| Queller & Strassmann | cooperation vs conflict | interaction signs; invasion tests | no |
| Godfrey-Smith | B, G, I | life-cycle data | B and G yes; I no |

---

## 4. ALife and prebiotic models of level-making: what counted, and what was measured

| Model | What counted as the new level | What was measured to show it |
|---|---|---|
| **Hypercycle** (Eigen & Schuster 1977) | Replicators coupled cyclically by catalysis. In Szathmáry's (2015) summary, *"Each member grows due to a combination of autocatalytic effect and heterocatalytic aid ... ecologically stable, but evolutionarily unstable because of the parasite problem."* | ODE analysis of coexistence and selection among replicators |
| **Stochastic corrector** (Szathmáry & Demeter 1987) | Reproducing compartments: *"replicative templates are competing within replicative compartments, whose selective values depend on the internal template composition"* | *"An Eigen equation at the compartment level is set up and calculated"*. *"The genetic information of templates is evaluated at both levels, and the higher (compartment) level successfully constrains the lower (template) one."* This is an explicit two-level selection accounting. |
| **Spiral waves** (Boerlijst & Hogeweg 1991) | Spiral waves emerging in a cellular automaton of hypercycle members | (i) Lineage labelling: *"the molecules from the middle have taken over the complete spiral region and the molecules from the periphery have disappeared"*, so the spiral core is the source of growth. (ii) Parasite challenge: *"We infected the situation ... randomly with 100 'deadly' parasites"*. Parasites were wiped out or, at worst, remained *"present as a cyst"*, i.e. *"positive selection for a strong altruistic property."* (iii) Robustness to diffusion. |
| **Waves versus vesicles** (Hogeweg & Takeuchi 2003; Takeuchi & Hogeweg 2009, 2012) | *Implicit* levels, where *"new levels of selection arise as large scale spatial patterns with a dynamics of their own"* (2003). *Explicit* levels imposed by vesicles. | *"macroscopic stability of a replicator system through the evolutionary dynamics on mesoscopic entities that counteract that of microscopic entities"*; long-term evolutionary trends that differ between the two mechanisms; *"a sharp transition ... in the long-term evolutionary dynamics of the compartmentalized system as a function of replicator mutation rate"* (2009). |
| **AlChemy** (Fontana & Buss 1994; Mathis et al. 2024) | L2 = *"self-maintaining metaorganizations composed of Level 1 organizations"*: ≥ 2 level-1 organizations plus "glue" | Original: organizations characterised by invariant grammars and algebra, self-repair and centre; L2 found by merging level-1 soups or spontaneously. 2024: Jaccard similarity of the merged soup to each input over 10⁶ collisions (inputs of 1000 expressions each, cap 2500), classed as dominance, coexistence or destruction over 455 pairs; repeated random perturbation for level 1. |
| **Hyperstructures** (Baas 1994, as restated by Rennard) | A second-order structure `S² = R(S¹ᵢ, Obs¹, Int¹)`. A property `P` is emergent iff `P ∈ Obs²(S²)` but `P ∉ Obs²(S¹ᵢ)` for all i. Iterating gives `Sᴺ`, a *"hyperstructure"* (p. 525). | Conceptual: no measurement procedure |
| **Ansatz** (Rasmussen, Baas, Mayer & Nilsson 2001) | A 2-D molecular-dynamics lattice gas: *"water and monomers at level one, polymers and water at level two, and micelles (polymer aggregates) and water at level three"* | The emergence of robust structures with new observables. The Ansatz itself is *"a conjecture of a necessary minimal complexity within the fundamental interacting structures"*. |
| **Critiques** (Gross & McMullin 2001; Dorin & McCormack 2002) | — | Gross & McMullin: *"essentially comparable phenomena can be produced with relatively simpler primitive objects. We also question the order classification of the micellar structures."* Dorin & McCormack built an *"infinitely-levelled, self-assembling dynamical hierarchy"* from fixed-complexity elements and concluded *"these new properties are trivial ... since the definition of the problem in the literature admits such trivial possibilities, more specific definitions are required."* |
| **Dynamical-hierarchies special issue** (*Artificial Life* 11(4), eds. Lenaerts, Chu & Watson 2005) | Levels as coarse-grainings under which the dynamics close | McGregor & Fernando: *"each level is a near-state-determined system, and levels are related to one another in a partial ordering"*, with *"a state-dependence criterion enforcing predictability within a level"* and *"a distinctness criterion enforcing the idea that the higher-level description must do more than just throw information away."* Jacobi: levels are *"smooth projective maps"* where *"each level describe[s] a self-contained deterministic dynamical system"*, and *"a necessary and sufficient condition ... is that the kernel of the differential of the map is tangent to an invariant manifold with respect to the flow."* Rowe, Vose & Wright: groupings *"compatible with those dynamics"*, with *"necessary and sufficient conditions"* in the linear case. Watson & Pollack: modularity with strong aggregate coupling. The issue also contains Prokopenko et al. and Altenberg. |
| **Open problem 8** (Bedau et al. 2000) | *"Create a formal framework for synthesizing dynamical hierarchies at all scales"* | Framed as a need: *"a formal framework for consistently describing such hierarchical systems, with some coarse-graining procedure for moving between levels"*, noting that *"it is unclear how to simulate either the generation or the functioning of such multi-level systems, because their dynamics are characterized by multiple length and time scales"*. |
| **Avida division of labour** (Goldsby, Dornhaus, Kerr & Ofria 2012) | Groups whose members specialise until they lose autonomy | Division of labour as Shannon mutual information between tasks and individuals (per the Dryad dataset notes). Loss of autonomy: *"individuals cease to be able to perform tasks in isolation, instead requiring the context of other group members. The simultaneous loss of functionality at a lower level and emergence of new functionality at a higher level"* |
| **DISHTINY** (Moreno & Ofria 2019) | Kin groups at two hierarchical levels, **specified by the platform** (*"provides simple cell-like organisms with the ability and incentive to unite into new individuals"*, with strategies *"manually designed"*) | *"reproductive division of labor and close cooperation among cells, including resource-sharing, aggregation of resource endowments for propagules, and emergence of an apoptosis response to somatic mutation"*, and which level resources were directed to. The paper notes transitions *"are challenging to induce or detect, even with computational organisms"*: the platform avoids the detection problem by defining the groups. |
| **Hash Chemistry** (Sayama 2019) | A higher-order entity is *"a multiset ... of entities"* | Cumulative count of unique replicating entities (grew almost linearly); *"the number of individual entities involved in a single replication event gradually increased over time"*; controls with random-number fitness. Caveats: the model *"does not have an explicit representation of 'higher-order entities' and thus they are not protected"*, and multisets of multisets (level 2 in Marco's sense) are left open. |

**Lessons for Chemart.**

- Every convincing level-making case added one of two things: a **space or compartment structure that makes the collective a unit of reproduction** (stochastic corrector, Hogeweg), or a **designed-in unit** (DISHTINY). The well-stirred cases (AlChemy, Hash Chemistry) either failed to combine level-1 units (Mathis) or produced unprotected aggregates (Sayama).
- The strongest evidence always combined **a collective-level dynamics or selection equation** (the compartment-level Eigen equation), **lineage or constraint evidence** (lineage labelling, the upper level constraining the lower) and **a challenge test** (parasites, perturbations).
- The dynamical-hierarchies line converges on one formal idea, **a level is a coarse-graining under which the dynamics close**, which is the bridge to §5. Its critics show that "new observable at the higher level" alone is far too weak.

---

## 5. Formal and information-theoretic identification of levels and individuals

**Common setup.**
- `X_t` is the micro-state: for Chemart, the vector of species counts from a trajectory.
- A **candidate** coarse-graining is `Z_t = f(X_t)`, or a system/environment split `(S_t, E_t)`.

**Feasibility preamble for a soup of 10²–10³ molecules (rules of thumb, not from the sources).** The raw micro-state is a multiset over an open species set, so its state space is astronomically large. Mutual-information estimates need far more samples than joint states. Every method below is therefore practical only after compressing the soup to a few (≲ 5–10) discrete macro-variables, such as presence or binned abundance of candidate level-1 organizations. Samples should come from ensembles of replicate runs, or from long windows where the organization is stationary, with bias correction and shuffle nulls. Evolutionary runs are non-stationary, so estimate within windows or across replicates at matched times.

### 5.1 Informational closure and autonomy (Bertschinger, Olbrich, Ay & Jost 2006; 2008)

**What it is.** Closure measures how much new information flows from the environment into the system. In Krakauer et al.'s (2020) notation:

```
nC = I(S_{n+1} ; E_n | S_n)          # = 0  ⇔  informationally closed
```

*"Note that closure does not require causal independence, it only states that all influences from the environment are predictable by the system"* (Krakauer et al. 2020). Bertschinger et al. (2008) propose autonomy as *"the conditional mutual information between consecutive states of the system conditioned on the history of the environment"*, or, when the system controls its environment, the mutual information between consecutive system states.

**Non-trivial closure.** Low `nC` is not enough, because an isolated system is trivially closed. The fix is non-trivial informational closure (NTIC): large information about the environment together with closure. In the notation of Krakauer's preprint: `A* − A = I(S_{n+1};E_n) − nC`.

**Inputs and feasibility.** A time series of `(S, E)` for a candidate split. It is feasible for a coarse `S`, e.g. `S` = abundance classes of the candidate composite's parts and `E` = the rest of the soup.

### 5.2 Comparison of level-identification methods (Pfante, Bertschinger, Olbrich, Ay & Jost 2014)

**What it is.** *"Levels of a complex system are characterized by the fact that they admit a closed functional description in terms of concepts and quantities intrinsic to that level ... we present four of these approaches and investigate their mutual relationships"*, for discrete dynamical systems linked by coarse-graining.

**Result** (as summarised by Rosas et al. 2024): they relate informational closure to causal closure (*"which they call 'commutativity'"*) and Markovianity. *"Their results imply that informational closure is equivalent to causal closure for spatial coarse-grainings."* Their closure condition is `I(Z_{t+1}; X_t | Z_t) = 0` for Markov `X`.

**Takeaway.** For Markov chemistries, informational closure is the practical test. It coincides with the causal version for coarse-grainings of the current state.

### 5.3 The information theory of individuality (Krakauer, Bertschinger, Olbrich, Flack & Ay 2020)

**What it is.** *"individuals are aggregates that preserve a measure of temporal integrity, i.e., 'propagate' information from their past into their futures."*

**How it works.** Decompose

```
I(S_n, E_n ; S_{n+1}) = I(S_{n+1};S_n) + I(S_{n+1};E_n|S_n)
                      = I(S_{n+1};E_n) + I(S_{n+1};S_n|E_n)
```

and define:

```
Organismal individuality   A* := I(S_{n+1} ; S_n)
Colonial individuality     A  := I(S_{n+1} ; S_n | E_n)
Environmental determination nC := I(S_{n+1} ; E_n | S_n)
```

With a partial information decomposition (shared SI, unique UI, complementary CI):

```
A*  = SI(S_{n+1}; S_n, E_n) + UI(S_{n+1}; S_n \ E_n)
A   = CI(S_{n+1}; S_n, E_n) + UI(S_{n+1}; S_n \ E_n)
nC  = CI(S_{n+1}; S_n, E_n) + UI(S_{n+1}; E_n \ S_n)
NTIC = SI − CI     # "environmental coding"
```

**What counts as an individual.** A partition `S` that maximises A* (organismal) or A (colonial) relative to alternatives.

**The critical caveat, verbatim:** *"Both individuality measures can only grow or stay constant with increasing system size ... Thus, they are not sufficient to detect the precise boundaries between individuals. In order to obtain precise boundaries we would need to impose a cost function—or regularizer—on system size"*, and *"Our objectives here are not to find the optimal partition"*.

**On AC data [proposal].** Compute A*, A and nC for:
- the candidate composite as `S`;
- each of its level-1 parts alone;
- random species sets of the same size (the size control the caveat demands).

A level-2 entity should beat its parts on A* or A per unit size and beat size-matched random sets. The partial-information terms depend on which decomposition is chosen (the authors acknowledge this), so rely on A*, A and nC, which do not.

### 5.4 Software in the natural world: informational, causal and computational closure (Rosas, Geiger, Luppi, Seth, Polani, Gastpar & Mediano 2024, arXiv:2402.09090)

**Definition 4, information closure.** A coarse-graining `Z` of `X` is informationally closed if

```
I( X⃗_t ; Z⃗^L_{t+1} | Z⃗_t ) = 0     for all L ∈ ℕ
```

(`X⃗_t` and `Z⃗_t` are past trajectories.) *"knowing the corresponding micro-state does not provide additional information about its future evolution over what can be obtained from the past of the coarse-graining itself."*

**Definition 3, causal closure.** *"A coarse-graining Z is said to be causally closed if the coarse-grainings induced by its ϵ-machine and υ-machine are equivalent"*. The ε-machine predicts `Z` from `Z`'s own past. The υ-machine predicts `Z` optimally from the micro past. When they coincide, *"all the causes of Zt are within its own level"*: the macro process is *"software-like"*.

**Definition 5, computational closure.** Coarse-grainings of inputs and states are *"computationally closed if the resulting transitions between coarse-grained states ... is also a deterministic automaton"*. Equivalently, coarse-graining and computing the ε-machine commute.

**Structure.**
- *"computationally closed levels are hierarchically organised into a lattice of nested computational structures ordered by coarse-graining relationships"*.
- Links to lumpability: closed levels exist iff causal states are strongly lumpable. The authors note this *"open[s] the door for efficient algorithms"*.
- Candidate explosion: *"the lattice of all coarse-grainings grows super-exponentially with the number of values a process can take"*.

**Caveat** (on the closely related dynamical independence of Barnett & Seth 2023): *"noise processes may satisfy dynamical independence while having no underlying causal structure"*. Closure must be non-trivial.

**On AC data.**
- The chemical master equation is a Markov chain, and a candidate level is closed iff the chain is (strongly) lumpable with respect to the partition `Z` induces.
- The deterministic analogue is exact lumping of kinetics. Wei & Kuo (1969) treat monomolecular systems; Jacobi's (2005) invariant-manifold condition covers smooth flows.
- Exact tests are cheap-to-moderate for Chemart's fixed ODE chemistries (Brusselator-type), but infeasible on the raw soup of constructive chemistries.
- Practical proxy [proposal]: estimate `I(X̃_t ; Z_{t+1} | Z_t)` with `X̃` a richer but finite summary (e.g. top-k species counts) and `Z` = abundances of the candidate's level-1 parts plus glue, across replicate runs. Alternatively, fit a Markov chain on a few "which organizations are present" states and test approximate lumpability.
- The main value of this framework for Chemart is its **lattice**: nested closed levels are exactly "networks of networks" in the dynamical sense.

### 5.5 Lumping, the chemistry-native route (Wei & Kuo 1969; Jacobi 2005; Rowe, Vose & Wright 2005)

Exact lumpability of reaction kinetics, where aggregated species obey their own closed kinetics, is the oldest formal version of "a level" in chemistry. Rowe et al.'s linear-dynamics conditions and Jacobi's smooth-flow condition generalise it. They give a yes/no test and a way to generate candidates for ODE chemistries. They say nothing about individuation.

### 5.6 Effective information and causal emergence (Hoel, Albantakis & Tononi 2013; Klein & Hoel 2020)

**What it is.** Effective information (EI) is the mutual information between the states of a system when its past state is set by intervention to maximum entropy. It is *"higher the more the mechanisms constrain the system's possible past and future states"*. Causal emergence is *"the gain in EI when moving from a micro to a macro level of analysis"* (Hoel et al. 2013). It occurs when macro mechanisms are *"more deterministic and/or less degenerate"*.

**The network version** (Klein & Hoel 2020), with each node's out-weights `Wᵢᵒᵘᵗ` summing to 1:

```
EI          = H(⟨W_i^out⟩) − ⟨H(W_i^out)⟩
determinism = log2(N) − ⟨H(W_i^out)⟩
degeneracy  = log2(N) − H(⟨W_i^out⟩)
```

Macro-nodes are found by a greedy search (*"Checking all possible groupings is computationally intractable for all but the smallest networks"*).

**Relevance to redundancy.** Biological networks show *"a significantly greater propensity for causal emergence"*. Their lower EI *"supports long-standing hypotheses about the role of redundancy, degeneracy, and noise in biological systems"*. Note that "degeneracy" here is a property of the graph's convergent paths, not molecular robustness in Marco's sense.

**Inputs and feasibility.** A transition matrix, or a weighted digraph turned into a random walk. For Chemart this could be the species graph with flux-normalised weights. It is cheap (greedy grouping on ≤ 10³ nodes).

**Limits.** A macro-node is an informative coarse-graining of a random walk on the graph, *not* a self-maintaining individual. Use it only as a candidate generator or as an independent check (Wimsatt).

### 5.7 Reconciling emergences: ΦID criteria (Rosas, Mediano, Jensen, Seth, Barrett, Carhart-Harris & Bor 2020)

**What it is.** Order-1 practical criteria for a candidate supervenient feature `V`:

```
Ψ(V) := I(V_t ; V_t′) − Σ_j I(X^j_t ; V_t′)
Δ(V) := max_j [ I(V_t ; X^j_t′) − Σ_i I(X^i_t ; X^j_t′) ]
Γ(V) := max_j I(V_t ; X^j_t′)
```

*"Ψ > 0 is a sufficient condition for V_t to be causally emergent ... Δ > 0 is a sufficient condition for V_t to exhibit downward causation ... Ψ > 0 and Γ = 0 is sufficient for causal decoupling."* The criteria scale *"linearly with system size (for k = 1)"*. They were demonstrated on the Game of Life, boids and ECoG.

**The redundancy bias, verbatim:** *"if there is redundancy in the system it will be harder to detect emergence, since redundancy will drive Ψ and Δ more negative. Furthermore, by summing all marginal mutual informations ... these measures effectively double-count redundancy up to n times"*. And: *"these criteria are unable to rule out emergence, as they are sufficient but not necessary conditions."*

**What this means for Chemart.** Ψ is cheap and tempting, but in a chemistry built to be redundant a negative Ψ means nothing. Use the ΦID versions that commit to a specific decomposition, or use closure tests.

### 5.8 Dynamical independence (Barnett & Seth 2023)

This identifies emergent macro-variables by optimising a projection for information closure. It is a practical optimiser for continuous or linear-Gaussian data, and it shares the triviality caveat above (Rosas et al. 2024).

### 5.9 Coarse-graining and slow variables (Flack 2012; 2017)

**What it is.** Empirically, in macaque societies, *"coarse-grained, statistical representations of collective dynamics are more predictive of the future state of the system than the constantly in-flux behavioural patterns at the individual level ... As an interaction history accumulates the coarse-grained representations consolidate. This constrains individual behaviour and provides the foundations for new levels of organization"* (Flack 2012). Conceptually, *"components collectively compute their macroscopic worlds through coarse-graining"*. There is *"downward causation when components tune behaviour in response to estimates of collectively computed macroscopic properties"*, with a weak and a strong form tied to *"the origins of new organizational levels"* (Flack 2017).

**On AC data.** A slow-variable test: the candidate macro-variables have much longer autocorrelation times and better out-of-sample prediction than the micro-variables. It is cheap.

**A speculative point [proposal].** Molecules that are small neural networks could *sense* collective state. That would make Flack's *strong* downward causation designable in Chemart, where it is impossible in λ-calculus soups.

### 5.10 The individuality of a glider (Beer 2014; 2015)

**What it is.** Beer (2015) formulates Game-of-Life entities *"as self-constructing networks of interdependent processes that maintain their own boundaries"*. Beer (2014) classifies perturbations: *"the set of possible perturbations to it can be divided into destructive and nondestructive subsets"*, and the non-destructive responses map its *"cognitive domain"*.

**On AC data [proposal].** This becomes a perturbation-domain test for a candidate composite. Apply knock-outs, knock-ins and random-molecule injections (Fontana & Buss's and Mathis's perturbation protocol). Record which perturbations destroy the composite, which destroy only one part, and which it repairs. The **whole's** repair domain should exceed its parts'. It needs interventions and is moderate in cost.

### Summary of section 5

| Method | Micro time series? | Candidate partition? | Interventions? | Output | Cost on a soup of 10²–10³ |
|---|---|---|---|---|---|
| Informational closure (Bertschinger; Pfante; Rosas Def. 4) | yes (or a transition matrix) | yes | no | `nC` or `I(X;Z′|Z)` ≥ 0, where 0 means closed | moderate, after compression |
| Causal / computational closure (Rosas 2024) | yes, ε/υ-machines | yes | no | yes/no plus a lattice | expensive; exact only on small chains |
| Krakauer A*, A, nC, NTIC | yes, (S,E) | yes | no | bits | moderate; needs a size control |
| Lumping (Wei & Kuo; Jacobi; Rowe) | ODE or linear model | yes | no | yes/no | cheap for fixed ODE chemistries |
| EI and causal emergence (Hoel; Klein & Hoel) | transition matrix or weighted graph | greedy search | intervention implicit | bits, macro-nodes | cheap |
| Ψ, Δ, Γ (Rosas 2020) | yes (parts and V) | candidate V | no | bits (sufficient only) | cheap, **redundancy-biased** |
| Slow variables (Flack) | yes | candidate aggregates | no | timescales, prediction gain | cheap |
| Perturbation domain (Beer) | full dynamics | candidate entity | **yes** | classes of perturbations | moderate |

---

## 6. Network hierarchy measures on reaction networks: what they see

### 6.1 Order and flow ("who drives whom")

- **Flow hierarchy** (Luo & Magee 2011; already in Chemart): the fraction of edges not on cycles.
- **Global reaching centrality** (Mones, Vicsek & Vicsek 2012). *"The local reaching centrality, C_R(i), of node i is the proportion of all nodes in the graph that can be reached from node i via outgoing edges"*, and

  ```
  GRC = Σ_{i∈V} [ C_R^max − C_R(i) ] / (N − 1)
  ```

  Cheap (BFS from every node).
- **Trophic levels and coherence.** Johnson et al. (2014): *"The trophic level of a species can be defined as the average trophic level of its prey, plus one"*, and incoherence `q` is *"the standard deviation of the distribution of trophic distances"*. MacKay, Johnson & Sansom (2020) remove the need for basal nodes:

  ```
  u_n = w_n^in + w_n^out,  v_n = w_n^in − w_n^out,  Λ = diag(u) − W − Wᵀ
  levels h solve  Λ h = v
  F0 = Σ_mn w_mn (h_n − h_m − 1)² / Σ_mn w_mn     # trophic coherence = 1 − F0
  ```

  Cheap (a sparse linear solve). On a reaction graph with the food set as source, `h` is roughly synthesis depth from food.
- **Treeness T, feedforwardness F, orderability O** (Corominas-Murtra, Goñi, Solé & Rodríguez-Caso 2013). The graph is condensed into its strongly connected components (SCCs), *"the so-called node weighted condensed graph"*.
  - `O` = *"the fraction of the nodes of the graph G that does not belong to any cycle"*;
  - `F` weights *"the impact of cyclic modules on the feedforward structure ... cyclic modules closer to the top of G will introduce a larger penalty"*;
  - `T ∈ [−1,1]` weighs *"how pyramidal is the structure and how unambiguous is its chain of command"*, from path entropies top-down versus bottom-up.
  - Together they form *"a 3D morphospace of hierarchies"*.
  - Cost: moderate. They average over paths of the condensed graph, and path counts can blow up in dense DAGs.

**Verdict: graph structure only (sense a).** Note especially that T/F/O **collapse autocatalytic cores, which are SCCs, into single nodes and score cycles as violations of hierarchy.** The very thing level 1 is made of counts against "hierarchy" here. These measures are good *morphospace axes* for Chemart (T/F/O was explicitly designed as a morphospace) but they are not level-2 detectors.

### 6.2 Nested modules ("modules within modules")

- **Hierarchical modularity** (Ravasz et al. 2002). The signature is *"C(k) ~ k⁻¹, in contrast to the k-independent C(k) predicted by both the scale-free and modular networks. This provides direct evidence for an inherently hierarchical organization."* Cheap, but it is only a global signature and names no modules.
- **Hierarchical random graphs** (Clauset, Moore & Newman 2008). A dendrogram model in which *"vertices divide into groups that further subdivide into groups of groups"*, fitted by MCMC. Moderate cost.
- **Sales-Pardo, Guimerà, Moreira & Amaral (2007).** They target *"inclusion hierarchies"*, i.e. nested organisation in the graph. Node affinity comes from co-classification across local maxima of the modularity landscape, followed by *"box-clustering"* level by level, with a null model. They argue any method *"must have a null output for networks, such as Erdős-Rényi random graphs"*, and criticise methods that *"yield a tree even for networks with no internal structure."*
- **Nested stochastic block model** (Peixoto 2014). *"a nested generative model that, through a complete description of the entire network hierarchy at multiple scales"*, avoids the resolution limit. It *"is capable of separating signal from noise, and thus will not lead to the identification of spurious modules"*, and *"is not restricted to purely assortative mixing patterns, directed or undirected graphs"*. It scales to large graphs. This is the best tool in this family.

**Verdict: nested graph structure (sense b without individuation).** A block is a set of nodes with similar edge statistics, not a closed, self-maintaining set. **Use** [proposal]: run the nested SBM on the reaction graph and check whether its blocks coincide with COT organisations or irrRAFs. Agreement is Wimsatt-style robustness; disagreement is informative.

### 6.3 Hierarchy inside AC formalisms

**Chemical organisation theory** (Dittrich & Speroni di Fenizio 2007).
- An organisation is *"a closed and self-maintaining set of components."*
- For catalytic flow systems: *"the set of all (semi-) organizations of a catalytic flow system forms an algebraic lattice ..., which has already been noted by Fontana and Buss"*. Closed sets always form a lattice, with join = closure of the union.
- Fixed points of the reaction ODE are instances of organisations (their Theorem 1).
- COT names AlChemy explicitly: *"Examples of catalytic flow system are the replicator equation, the hypercycle, the more general catalytic network equation, or AlChemy"*. In such systems *"we can easily check, whether a set O is an organization by just checking whether it is closed and whether each molecule in that set is produced by that set."* (Mathis et al. 2024 describe AlChemy's collisions as `A + B → A + B + C`, with a random expression removed to keep N constant.)

**[proposal]** In COT terms, Fontana & Buss's L2 is an organisation `O` that is the **join of ≥ 2 incomparable smaller organisations**, `O = O₁ ⊔ O₂`, whose **glue** `O \ (O₁ ∪ O₂)` is non-empty. The glue species are produced only by cross-organisation reactions. This matches Mathis's "glue" and is *computable from Chemart's existing organisation list*.

What the lattice alone does not give: whether the join is dynamically realised and stable (Mathis found it usually is not), and whether the parts remain distinct rather than overconnected (McShea). Enumerating *all* organisations can be exponential in the worst case: n independent self-maintaining species already give 2ⁿ organisations. Checking one candidate set is cheap in catalytic flow systems (the quote above).

**RAF theory** (Hordijk, Steel & Kauffman 2012).
- A maxRAF *"could possibly consist of several smaller (independent or overlapping) subsets which themselves are RAF sets (subRAFs). If such a subRAF cannot be reduced any further without losing the RAF property, we refer to it as an irreducible RAF (irrRAF)."*
- SubRAFs form a poset (a Hasse diagram).
- *"there exist polynomial time ... algorithms that solve the following problems: (i) generate a list of all the maximal proper subRAFs of R′; (ii) determine whether or not R′ is the union of two proper subRAFs, and if so find all such pairs of subRAFs"*.
- On level 2: *"one could imagine a collection of mutually dependent RAF sets forming a meta-RAF set: one set enabling (catalyzing) the existence of another, in mutually beneficial ways. In other words, self-sustaining, functionally closed structures can arise at a higher level (an autocatalytic set of autocatalytic sets)"*. This is speculative in the paper, and **no detection criterion is given.**

**[proposal]** Build a *RAF-level graph*: nodes are irrRAFs or closed subRAFs, and an edge `i → j` means `i` supplies a catalyst or needed reactant for a reaction of `j`. A meta-RAF candidate is a set of RAF-nodes in which every node is supported by another node in the set, i.e. RAF-like at the node level. Test (ii) above is the cheap first filter. This needs chemistries that expose catalysts.

**P systems** (Păun 2000; Păun 2010).
- The membrane structure is *"a hierarchical arrangement of membranes ... delimiting compartments where multisets of objects are placed"*; *"compartments can contain other compartments"*.
- Compositional nesting is **given by construction**: the depth of the membrane tree.
- Individuation still has to be tested: is each compartment's content an organisation, and are the compartments' fates coupled?
- With membrane division, P systems support MLS2 accounting directly (compartment genealogies).

### 6.4 Verdict table: which measures can see "networks of networks"?

| Measure | Hierarchy sense | Sees composition? | Sees individuation? | Verdict for level 2 | Suggested cost tier |
|---|---|---|---|---|---|
| Flow hierarchy (Luo & Magee) | a | no | no | graph structure only | cheap (exists) |
| GRC (Mones) | a | no | no | graph structure only | cheap |
| Trophic levels and coherence (Johnson; MacKay) | a | no | no | graph structure only | cheap |
| T/F/O (Corominas-Murtra) | a | no (collapses SCCs) | no | graph structure only | moderate |
| Modularity, NODF (existing) | b (one level) | weakly | no | graph structure only | cheap (exist) |
| C(k) ~ k⁻¹ (Ravasz) | b (signature) | no entities | no | graph structure only | cheap |
| HRG (Clauset); Sales-Pardo | b | nested modules | no | graph structure only | moderate |
| Nested SBM (Peixoto) | b | nested blocks | no | graph structure only; good cross-check | moderate |
| EI macro-nodes (Klein & Hoel) | b/dynamics | coarse-grained nodes | no | candidate generator | cheap |
| COT lattice plus glue | b over self-maintaining sets | **yes** | **partly** (closure, self-maintenance) | **structural L2 candidates** | exponential to enumerate; cheap to check one set |
| subRAF poset / meta-RAF | b over autocatalytic sets | **yes** | **partly** | **structural L2 candidates** | polynomial (tests) |
| P-system membrane tree | b by construction | **yes** | no (must test) | nesting given | cheap |
| Closure tests (§5.1–5.4) on L2 candidates | dynamics | via the candidate | **yes** (closure, non-triviality) | **level-ness test** | moderate–expensive |
| Knockout and perturbation domain (§3.1, §5.10) | dynamics | via the candidate | **yes** (dependence, self-repair) | **individuation test** | moderate |
| MLS2 Price, Okasha stages, Michod covariance (§3.3–3.4) | c | via collectives | **yes** (evolutionary) | **individuality test, needs reproduction** | expensive |

---

## 7. Synthesis: what a level-2 detector for an artificial chemistry would need

This section is a design **[proposal]**, but every stage is tied to the sources that motivate it. Stages A–D need only what Chemart already produces: the reaction network with catalysts and food set, trajectories, and structure strings. Stage E needs compartments, space or scaffolding.

**Stage A: find the level-1 units (structural).**
Compute COT organisations and/or irrRAFs and closed subRAFs (Dittrich & Speroni di Fenizio 2007; Hordijk et al. 2012). Keep those that are dynamically realised and self-repairing under repeated random perturbation. This is Fontana & Buss's criterion as operationalised by Mathis et al. 2024.

**Stage B: generate level-2 candidates (structural).**
Candidates are organisations or RAFs that are **joins or unions of ≥ 2 incomparable level-1 units** (the polynomial RAF test; the COT join) and that have non-empty **glue**, meaning species produced only by cross-unit reactions (Mathis et al. 2024; Fontana & Buss). Score each with McShea's rubric:
- differentiation = the number of distinct unit types (sublevel b);
- intermediate parts = sub-composites (sublevel c);
- nesting depth = composites of composites (level 3).

Keep a per-run **maximum**, following McShea's trend in the maximum.

**Stage C: individuation tests (dynamical, interventional).**
1. **Mutual dependence and loss of autonomy.** Knock out each part. Does the other part still persist, and does the glue vanish? This follows Maynard Smith & Szathmáry's criterion and Goldsby et al.'s loss-of-autonomy measure. Symmetric dependence indicates integration; one-way dependence indicates a host with a passenger.
2. **Self-repair domain.** Map destructive and non-destructive perturbations for the whole versus its parts (Beer 2014; Fontana & Buss's self-repair).
3. **Conflict.** Inject parasites (species catalysed by the composite that give nothing back) and cheater-like mutants, and measure invasion. This is Boerlijst & Hogeweg's challenge and Queller & Strassmann's "low conflict".
4. **Not overconnected.** The parts should remain distinct organisations inside the whole, with glue as a minority. Otherwise the result is one big level-1 organisation (McShea's overconnectedness; Simon's near-decomposability).

**Stage D: level-ness (information and timescales).**
Let `Z` be the abundances of the parts plus glue.
- **Closure:** `I(X̃_t; Z_{t+1} | Z_t) ≈ 0`, and non-trivial (NTIC > 0) (Bertschinger et al. 2006, 2008; Pfante et al. 2014; Rosas et al. 2024).
- **Individuality:** Krakauer's A*, A and nC for the composite versus each part versus size-matched random species sets, which serves as the regulariser (Krakauer et al. 2020).
- **Timescale separation:** Simon, Salthe and Flack.
- **Do not** rely on Ψ and Δ alone, because of the redundancy bias (Rosas et al. 2020).

**Stage E: evolutionary individuality (only where collectives reproduce).**
Collective reproduction comes from compartments or P-system division, from space (Hogeweg), or, for well-stirred chemistries, from **ecological scaffolding**: many reactors or patches with periodic dispersal and bottlenecks. Black, Bourrat & Rainey (2020) show that *"a minimal ecological structure comprising patchily distributed resources and between-patch dispersal can scaffold Darwinian-like properties on collectives"*, and call it *"an ecological recipe for experimental realization of evolutionary transitions"*. The quantities to measure then are:
- the Price partition between and within collectives (Price 1972; Okasha 2006);
- Okasha's stage, from how collective fitness₂ relates to mean particle fitness;
- Michod's covariance effect;
- Godfrey-Smith's B, G and I;
- the fidelity of composite reproduction, which is Szathmáry's "reproduction versus mere replication".

For Fontana & Buss's egalitarian L2, Szathmáry's theory predicts that stabilisation needs *fair or synchronised reproduction of the partners*. Stage E is where that prediction can be tested.

**Acceptance rule (Wimsatt).**
Call a composite level 2 only if Stages B, C and D agree. Report the stage scores as **separate morphospace axes**:
- structural depth and differentiation;
- dependence and repair;
- closure and individuality;
- where available, selective individuality.

Separate axes let Chemart place chemistries in a hierarchy morphospace and test the redundancy hypothesis. The prediction to test: more redundant molecules raise the fraction of merged pairs that coexist, the glue fraction, and the closure and individuality scores of composites.

**Where Chemart stands.**
- Flow hierarchy, modularity, NODF and bow-tie sit in §6.1–6.2: graph structure only.
- The organisation list is Stage A input; the lattice and glue step (Stage B) is not yet used.
- Nothing yet does Stages C–E.
- The cheapest high-value additions are Stage B (glue plus decomposability, on existing organisations and RAFs) and the Stage C knockout test, which together directly test Mathis et al.'s ceiling claim across chemistries.

**Pitfalls, with sources.**
1. Graph hierarchy is not an entity (§6).
2. "New observable" definitions are trivially satisfiable (Gross & McMullin 2001; Dorin & McCormack 2002).
3. Redundancy pushes Ψ and Δ negative (Rosas et al. 2020), while EI-based emergence is higher in degenerate networks (Klein & Hoel 2020). The two families respond oppositely to redundancy.
4. Individuality scores never decrease with size (Krakauer et al. 2020), so control for size.
5. There are super-exponentially many coarse-grainings (Rosas et al. 2024), so generate candidates structurally.
6. Trivial closure: noise and isolated systems are "closed" (Rosas et al. 2024; Bertschinger et al. 2008), so require non-trivial closure.
7. Collectively autocatalytic is not a hypercycle (Szathmáry 2015).
8. Overconnectedness dissolves levels (McShea 2001).
9. Evolutionary runs are non-stationary, so estimate across replicates or windows.

---

## 8. References

Status key: **VERIFIED** means authors, year, title and venue were confirmed in this session through Crossref, Europe PMC or arXiv metadata, or the publisher's or author's page. The note after each entry says whether the full text, the abstract, or only the metadata was read. No entry is UNVERIFIED. Where content was taken from a secondary source, the note says so.

**Anchors: artificial chemistries**
1. Fontana, W. & Buss, L. W. (1994a). "The arrival of the fittest": Toward a theory of biological organization. *Bulletin of Mathematical Biology* 56(1): 1–64. https://doi.org/10.1007/BF02458289. **VERIFIED.** Crossref metadata; abstract read via SFI Working Paper 93-09-055 (https://econpapers.repec.org/RePEc:wop:safiwp:93-09-055).
2. Fontana, W. & Buss, L. W. (1994b). What would be conserved if "the tape were played twice"? *PNAS* 91(2): 757–761. https://doi.org/10.1073/pnas.91.2.757. **VERIFIED.** Abstract read.
3. Mathis, C., Patel, D., Weimer, W. & Forrest, S. (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos* 34(9): 093142. https://doi.org/10.1063/5.0207358; arXiv:2408.12137. **VERIFIED.** Full text read (arXiv v2).

**Disambiguation**
4. Gao, J., Buldyrev, S. V., Stanley, H. E. & Havlin, S. (2012). Networks formed from interdependent networks. *Nature Physics* 8: 40–48. https://doi.org/10.1038/nphys2180. **VERIFIED.** Metadata and abstract summary.
5. Luo, J. & Magee, C. L. (2011). Detecting evolving patterns of self-organizing networks by flow hierarchy measurement. *Complexity* 16(6): 53–61. https://doi.org/10.1002/cplx.20368. **VERIFIED.** Abstract read.

**Conceptual foundations**
6. Simon, H. A. (1962). The architecture of complexity. *Proceedings of the American Philosophical Society* 106(6): 467–482. JSTOR: http://links.jstor.org/sici?sici=0003-049X%2819621212%29106%3A6%3C467%3ATAOC%3E2.0.CO%3B2-1. **VERIFIED.** Full text read.
7. Wimsatt, W. C. (1976a). Reductionism, levels of organization, and the mind-body problem. In G. Globus, G. Maxwell & I. Savodnik (eds), *Consciousness and the Brain*, Plenum, pp. 205–267. https://doi.org/10.1007/978-1-4684-2196-5_9. **VERIFIED.** Metadata; quotes via SEP (ref. 9).
8. Wimsatt, W. C. (1994). The ontology of complex systems: levels of organization, perspectives, and causal thickets. *Canadian Journal of Philosophy*, Supplementary Vol. 20: 207–274. https://doi.org/10.1080/00455091.1994.10717400. **VERIFIED.** Metadata; quotes via SEP.
9. Eronen, M. I. & Brooks, D. S. (2018; substantive revision 2023). Levels of Organization in Biology. *Stanford Encyclopedia of Philosophy*. https://plato.stanford.edu/entries/levels-org-biology/. **VERIFIED.** Full text read.
10. Salthe, S. N. (2001). Summary of the Principles of Hierarchy Theory. Web text, http://www.nbi.dk/~natphil/salthe/Summary_of_the_Principles_o.pdf. **VERIFIED.** Full text read.
11. Salthe, S. N. (2012). Hierarchical structures. *Axiomathes* 22: 355–383. https://doi.org/10.1007/s10516-012-9185-0. **VERIFIED.** Metadata only.
12. Pattee, H. H. (1973). The physical basis and origin of hierarchical control. In H. H. Pattee (ed.), *Hierarchy Theory: The Challenge of Complex Systems*, New York: George Braziller, pp. 71–108. **VERIFIED.** Bibliographic entry via SEP bibliography and book listings; content via SEP (not read directly).
13. McShea, D. W. (2001). The hierarchical structure of organisms: a scale and documentation of a trend in the maximum. *Paleobiology* 27(2): 405–423. https://doi.org/10.1666/0094-8373(2001)027<0405:THSOOA>2.0.CO;2. **VERIFIED.** Full text pp. 405–410 read.

**Major transitions and individuality**
14. Maynard Smith, J. & Szathmáry, E. (1995). *The Major Transitions in Evolution*. Oxford: W. H. Freeman/Spektrum (reissued by OUP, https://global.oup.com/academic/product/the-major-transitions-in-evolution-9780198502944). **VERIFIED.** Bibliographic; both characterisations (one with p. 8) quoted via Okasha 2022 (ref. 19).
15. Szathmáry, E. (2015). Toward major evolutionary transitions theory 2.0. *PNAS* 112(33): 10104–10111. https://doi.org/10.1073/pnas.1421398112. **VERIFIED.** Full text read (PMC4547294).
16. Michod, R. E. & Nedelcu, A. M. (2003). On the reorganization of fitness during evolutionary transitions in individuality. *Integrative and Comparative Biology* 43(1): 64–73. https://doi.org/10.1093/icb/43.1.64. **VERIFIED.** Abstract read; p. 66 quote via Szathmáry 2015.
17. Michod, R. E. (2006). The group covariance effect and fitness trade-offs during evolutionary transitions in individuality. *PNAS* 103(24): 9113–9117. https://doi.org/10.1073/pnas.0601080103. **VERIFIED.** Full text read (PMC1482575).
18. Okasha, S. (2005). Multilevel selection and the major transitions in evolution. *Philosophy of Science* 72(5): 1013–1025. https://doi.org/10.1086/508102. **VERIFIED.** Abstract read; p. 1023 quote via Szathmáry 2015.
19. Okasha, S. (2022). The Major Transitions in Evolution—A Philosophy-of-Science Perspective. *Frontiers in Ecology and Evolution* 10: 793824. https://doi.org/10.3389/fevo.2022.793824. **VERIFIED.** Full text read.
20. Okasha, S. (2006). *Evolution and the Levels of Selection*. Oxford: Clarendon Press. https://doi.org/10.1093/acprof:oso/9780199267972.001.0001. **VERIFIED.** Metadata; MLS1/MLS2 definitions via reviews and Szathmáry 2015.
21. Damuth, J. & Heisler, I. L. (1988). Alternative formulations of multilevel selection. *Biology & Philosophy* 3(4): 407–430. https://doi.org/10.1007/BF00647962. **VERIFIED.** Metadata only.
22. Price, G. R. (1972). Extension of covariance selection mathematics. *Annals of Human Genetics* 35: 485–490. https://doi.org/10.1111/j.1469-1809.1957.tb01874.x (publisher DOI as registered). **VERIFIED.** Full text pp. 485–487 read.
23. Queller, D. C. & Strassmann, J. E. (2009). Beyond society: the evolution of organismality. *Philosophical Transactions of the Royal Society B* 364(1533): 3143–3155. https://doi.org/10.1098/rstb.2009.0095. **VERIFIED.** Abstract read.
24. Godfrey-Smith, P. (2009). *Darwinian Populations and Natural Selection*. Oxford University Press. https://doi.org/10.1093/acprof:osobl/9780199552047.001.0001. **VERIFIED.** Metadata; B/G/I content via ref. 25.
25. Godfrey-Smith, P. (2012 preprint). Darwinian Individuals. To appear in F. Bouchard & P. Huneman (eds), *From Groups to Individuals*, MIT Press. https://petergodfreysmith.com/PGS_Darwinian_Individuals.pdf. **VERIFIED.** Full text read.
26. Black, A. J., Bourrat, P. & Rainey, P. B. (2020). Ecological scaffolding and the evolution of individuality. *Nature Ecology & Evolution* 4: 426–436. https://doi.org/10.1038/s41559-019-1086-9. **VERIFIED.** Abstract read.

**ALife and prebiotic models of level-making**
27. Eigen, M. & Schuster, P. (1977). The hypercycle. A principle of natural self-organization. Part A: Emergence of the hypercycle. *Naturwissenschaften* 64(11): 541–565. https://doi.org/10.1007/BF00450633. **VERIFIED.** Metadata; content via Szathmáry 2015.
28. Szathmáry, E. & Demeter, L. (1987). Group selection of early replicators and the origin of life. *Journal of Theoretical Biology* 128(4): 463–486. https://doi.org/10.1016/S0022-5193(87)80191-1. **VERIFIED.** Abstract read.
29. Boerlijst, M. C. & Hogeweg, P. (1991). Spiral wave structure in pre-biotic evolution: hypercycles stable against parasites. *Physica D* 48(1): 17–28. https://doi.org/10.1016/0167-2789(91)90049-F. **VERIFIED.** Full text read (https://tbb.bio.uu.nl/pdf/Boerlijst.pd91-48.pdf).
30. Hogeweg, P. & Takeuchi, N. (2003). Multilevel selection in models of prebiotic evolution: compartments and spatial self-organization. *Origins of Life and Evolution of the Biosphere* 33(4–5): 375–403. https://doi.org/10.1023/A:1025754907141. **VERIFIED.** Abstract read.
31. Takeuchi, N. & Hogeweg, P. (2009). Multilevel selection in models of prebiotic evolution II: a direct comparison of compartmentalization and spatial self-organization. *PLoS Computational Biology* 5(10): e1000542. https://doi.org/10.1371/journal.pcbi.1000542. **VERIFIED.** Abstract read.
32. Takeuchi, N. & Hogeweg, P. (2012). Evolutionary dynamics of RNA-like replicator systems: a bioinformatic approach to the origin of life. *Physics of Life Reviews* 9(3): 219–263. https://doi.org/10.1016/j.plrev.2012.06.001. **VERIFIED.** Abstract read.
33. Baas, N. A. (1994). Emergence, hierarchies, and hyperstructures. In C. G. Langton (ed.), *Artificial Life III*, Reading, MA: Addison-Wesley, pp. 515–537. **VERIFIED.** Bibliographic, cited identically by refs. 34 and 37; formalism as restated in ref. 34 (original not read).
34. Rennard, J.-P. (2006; arXiv 2007). Artificiality in social sciences. arXiv:cs/0701087 (draft chapter for *Handbook of Research on Nature Inspired Computing for Economics and Management*). https://arxiv.org/abs/cs/0701087. **VERIFIED.** Full text read; used only for its restatement of Baas 1994.
35. Rasmussen, S., Baas, N. A., Mayer, B. & Nilsson, M. (2001). Ansatz for dynamical hierarchies. *Artificial Life* 7(4): 329–353. https://doi.org/10.1162/106454601317296988. **VERIFIED.** Abstract read.
36. Gross, D. & McMullin, B. (2001). Is it the right ansatz? *Artificial Life* 7(4): 355–365. https://doi.org/10.1162/106454601317296997. **VERIFIED.** Abstract read.
37. Dorin, A. & McCormack, J. (2002). Self-assembling dynamical hierarchies. In R. K. Standish, H. A. Abbass & M. A. Bedau (eds), *Artificial Life VIII*, MIT Press, pp. 423–428. http://alife8.alife.org/proceedings/sub3843.pdf. **VERIFIED.** Full text read.
38. Lenaerts, T., Chu, D. & Watson, R. (2005). Dynamical hierarchies (guest editors' introduction). *Artificial Life* 11(4): 403–405. https://doi.org/10.1162/106454605774270606. **VERIFIED.** Metadata and issue table of contents; editorial text not read.
39. McGregor, S. & Fernando, C. (2005). Levels of description: a novel approach to dynamical hierarchies. *Artificial Life* 11(4): 459–472. https://doi.org/10.1162/106454605774270615. **VERIFIED.** Abstract read.
40. Jacobi, M. N. (2005). Hierarchical organization in smooth dynamical systems. *Artificial Life* 11(4): 493–512. https://doi.org/10.1162/106454605774270598. **VERIFIED.** Abstract read.
41. Rowe, J. E., Vose, M. D. & Wright, A. H. (2005). State aggregation and population dynamics in linear systems. *Artificial Life* 11(4): 473–492. https://doi.org/10.1162/106454605774270624. **VERIFIED.** Abstract read.
42. Watson, R. A. & Pollack, J. B. (2005). Modular interdependency in complex dynamical systems. *Artificial Life* 11(4): 445–457. https://doi.org/10.1162/106454605774270589. **VERIFIED.** Abstract read.
43. Bedau, M. A., McCaskill, J. S., Packard, N. H., Rasmussen, S., Adami, C., Green, D. G., Ikegami, T., Kaneko, K. & Ray, T. S. (2000). Open problems in artificial life. *Artificial Life* 6(4): 363–376. https://doi.org/10.1162/106454600300103683. **VERIFIED.** Full text read (https://people.reed.edu/~mab/papers/ALife.6.4.pdf).
44. Goldsby, H. J., Dornhaus, A., Kerr, B. & Ofria, C. (2012). Task-switching costs promote the evolution of division of labor and shifts in individuality. *PNAS* 109(34): 13686–13691. https://doi.org/10.1073/pnas.1202233109. **VERIFIED.** Abstract read; metric via Dryad dataset https://doi.org/10.5061/dryad.f8j02.
45. Moreno, M. A. & Ofria, C. (2019). Toward open-ended fraternal transitions in individuality. *Artificial Life* 25(2): 117–133. https://doi.org/10.1162/artl_a_00284. **VERIFIED.** Abstract read.
46. Sayama, H. (2019). Cardinality leap for open-ended evolution: theoretical consideration and demonstration by Hash Chemistry. *Artificial Life* 25(2): 104–116. https://doi.org/10.1162/artl_a_00283; arXiv:1806.06628. **VERIFIED.** Abstract and arXiv full text read.

**Formal and information-theoretic identification**
47. Krakauer, D., Bertschinger, N., Olbrich, E., Flack, J. C. & Ay, N. (2020). The information theory of individuality. *Theory in Biosciences* 139(2): 209–223. https://doi.org/10.1007/s12064-020-00313-7. **VERIFIED.** Full text read (PMC7244620; arXiv:1412.2447).
48. Bertschinger, N., Olbrich, E., Ay, N. & Jost, J. (2006). Information and closure in systems theory. In *Explorations in the Complexity of Possible Life: Proceedings of the 7th German Workshop of Artificial Life*, Amsterdam: IOS Press, pp. 9–21. **VERIFIED.** Bibliographic (reference list of ref. 50 and the author's ResearchGate listing); content via refs. 47 and 50 (original not read).
49. Bertschinger, N., Olbrich, E., Ay, N. & Jost, J. (2008). Autonomy: an information theoretic perspective. *BioSystems* 91(2): 331–345. https://doi.org/10.1016/j.biosystems.2007.05.018. **VERIFIED.** Abstract read.
50. Rosas, F. E., Geiger, B. C., Luppi, A. I., Seth, A. K., Polani, D., Gastpar, M. & Mediano, P. A. M. (2024). Software in the natural world: a computational approach to hierarchical emergence. arXiv:2402.09090 (v2). https://arxiv.org/abs/2402.09090. **VERIFIED.** Full text read.
51. Pfante, O., Bertschinger, N., Olbrich, E., Ay, N. & Jost, J. (2014). Comparison between different methods of level identification. *Advances in Complex Systems* 17(2): 1450007. https://doi.org/10.1142/S0219525914500076. **VERIFIED.** Abstract read; results as summarised by ref. 50.
52. Hoel, E. P., Albantakis, L. & Tononi, G. (2013). Quantifying causal emergence shows that macro can beat micro. *PNAS* 110(49): 19790–19795. https://doi.org/10.1073/pnas.1314922110. **VERIFIED.** Abstract read.
53. Klein, B. & Hoel, E. (2020). The emergence of informative higher scales in complex networks. *Complexity* 2020: 8932526. https://doi.org/10.1155/2020/8932526; arXiv:1907.03902. **VERIFIED.** arXiv full text read.
54. Rosas, F. E., Mediano, P. A. M., Jensen, H. J., Seth, A. K., Barrett, A. B., Carhart-Harris, R. L. & Bor, D. (2020). Reconciling emergences: an information-theoretic approach to identify causal emergence in multivariate data. *PLoS Computational Biology* 16(12): e1008289. https://doi.org/10.1371/journal.pcbi.1008289. **VERIFIED.** arXiv full text read (arXiv:2004.08220).
55. Barnett, L. & Seth, A. K. (2023). Dynamical independence: discovering emergent macroscopic processes in complex dynamical systems. *Physical Review E* 108(1): 014304. https://doi.org/10.1103/PhysRevE.108.014304. **VERIFIED.** Metadata; content via ref. 50.
56. Flack, J. C. (2012). Multiple time-scales and the developmental dynamics of social systems. *Philosophical Transactions of the Royal Society B* 367(1597): 1802–1810. https://doi.org/10.1098/rstb.2011.0214. **VERIFIED.** Abstract read.
57. Flack, J. C. (2017). Coarse-graining as a downward causation mechanism. *Philosophical Transactions of the Royal Society A* 375(2109): 20160338. https://doi.org/10.1098/rsta.2016.0338. **VERIFIED.** Abstract read.
58. Beer, R. D. (2014). The cognitive domain of a glider in the Game of Life. *Artificial Life* 20(2): 183–206. https://doi.org/10.1162/ARTL_a_00125. **VERIFIED.** Abstract read.
59. Beer, R. D. (2015). Characterizing autopoiesis in the Game of Life. *Artificial Life* 21(1): 1–19. https://doi.org/10.1162/ARTL_a_00143. **VERIFIED.** Abstract read.
60. Wei, J. & Kuo, J. C. W. (1969). Lumping analysis in monomolecular reaction systems: analysis of the exactly lumpable system. *Industrial & Engineering Chemistry Fundamentals* 8(1): 114–123. https://doi.org/10.1021/i160029a019. **VERIFIED.** Metadata only.

**Network hierarchy measures and hierarchy inside AC formalisms**
61. Corominas-Murtra, B., Goñi, J., Solé, R. V. & Rodríguez-Caso, C. (2013). On the origins of hierarchy in complex networks. *PNAS* 110(33): 13316–13321. https://doi.org/10.1073/pnas.1300832110; arXiv:1303.2503. **VERIFIED.** arXiv full text read.
62. Mones, E., Vicsek, L. & Vicsek, T. (2012). Hierarchy measure for complex networks. *PLoS ONE* 7(3): e33799. https://doi.org/10.1371/journal.pone.0033799. **VERIFIED.** Full text read.
63. Johnson, S., Domínguez-García, V., Donetti, L. & Muñoz, M. A. (2014). Trophic coherence determines food-web stability. *PNAS* 111(50): 17923–17928. https://doi.org/10.1073/pnas.1409077111; arXiv:1404.7728. **VERIFIED.** arXiv full text read.
64. MacKay, R. S., Johnson, S. & Sansom, B. (2020). How directed is a directed network? *Royal Society Open Science* 7(9): 201138. https://doi.org/10.1098/rsos.201138; arXiv:2001.05173. **VERIFIED.** arXiv full text read.
65. Ravasz, E., Somera, A. L., Mongru, D. A., Oltvai, Z. N. & Barabási, A.-L. (2002). Hierarchical organization of modularity in metabolic networks. *Science* 297(5586): 1551–1555. https://doi.org/10.1126/science.1073374; arXiv:cond-mat/0209244. **VERIFIED.** arXiv full text read.
66. Clauset, A., Moore, C. & Newman, M. E. J. (2008). Hierarchical structure and the prediction of missing links in networks. *Nature* 453(7191): 98–101. https://doi.org/10.1038/nature06830; arXiv:0811.0484. **VERIFIED.** Abstract read.
67. Sales-Pardo, M., Guimerà, R., Moreira, A. A. & Amaral, L. A. N. (2007). Extracting the hierarchical organization of complex systems. *PNAS* 104(39): 15224–15229. https://doi.org/10.1073/pnas.0703740104; arXiv:0705.1679. **VERIFIED.** arXiv full text read.
68. Peixoto, T. P. (2014). Hierarchical block structures and high-resolution model selection in large networks. *Physical Review X* 4(1): 011047. https://doi.org/10.1103/PhysRevX.4.011047; arXiv:1310.4377. **VERIFIED.** Abstract read.
69. Dittrich, P. & Speroni di Fenizio, P. (2007). Chemical organisation theory. *Bulletin of Mathematical Biology* 69(4): 1199–1231. https://doi.org/10.1007/s11538-006-9130-8; arXiv:q-bio/0501016. **VERIFIED.** Abstract and arXiv full text read.
70. Hordijk, W., Steel, M. & Kauffman, S. (2012). The structure of autocatalytic sets: evolvability, enablement, and emergence. *Acta Biotheoretica* 60(4): 379–392. https://doi.org/10.1007/s10441-012-9165-1; arXiv:1205.0584. **VERIFIED.** arXiv full text read.
71. Păun, G. (2000). Computing with membranes. *Journal of Computer and System Sciences* 61(1): 108–143. https://doi.org/10.1006/jcss.1999.1693. **VERIFIED.** Metadata only.
72. Păun, G. (2010). Membrane computing. *Scholarpedia* 5(1): 9259. https://doi.org/10.4249/scholarpedia.9259. **VERIFIED.** Full text read.

*Note on count: 72 entries, all cited in the text. This is more than the 20–35 target. The brief named about 45 works. The rest either supply the verified text for a named work (the SEP for Wimsatt and Pattee; Okasha 2022 for Maynard Smith & Szathmáry; Rennard for Baas; Godfrey-Smith 2012 for 2009; Rosas 2024 for Bertschinger 2006 and Pfante 2014; Szathmáry 2015 for Okasha 2005 and Michod & Nedelcu) or carry a load-bearing caveat or recipe (Gross & McMullin; Dorin & McCormack; Klein & Hoel; Black et al.; Michod 2006; the 2005 special-issue papers).*
