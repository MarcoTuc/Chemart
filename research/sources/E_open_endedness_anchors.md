# E. Measuring open-endedness and anchoring "complex" in real chemistry

This report covers three things: how to measure complexity growth and open-ended evolution (OEE) over evolutionary time; what organisation, complexity and hierarchy measures show on real chemical and biochemical networks; and some brief notes on morphospace method. Each measure is judged by two questions. Can it run on a Chemart `Trajectory`? A trajectory is a list of frames, and each frame holds `state` = {species: amount}, `fired` = [[reactants], [products], count] and `observables`; the trajectory also provides `turnover()` and `window()`. And would the measure register a new level (networks of networks, Fontana & Buss level 2), or only more species at the same level?

**Verification.** I checked every reference in this session against Crossref, OpenAlex, Semantic Scholar or Europe PMC metadata, and against the abstract or full text. Quotes come from the sources I actually read. The few details I saw only in secondary coverage are flagged inline and in §5.

---

## 1. Landscape summary

**The OEE measures come in four families.**

1. **Component statistics over a population.** Bedau & Packard's *evolutionary activity* led to Bedau–Snyder–Packard's classes of long-term dynamics. Channon refined the normalisation, Dolson et al.'s MODES toolbox (change, novelty, complexity, ecology) built on the same idea, and Droop & Hickinbotham's QNN is a shadow-free variant. Sayama's cardinality-leap counts belong here too. All of them need a declared unit, the "component" (allele, genotype, species, family), and a way to separate adaptive persistence from drift: a neutral *shadow* run, a *persistence filter*, or the fossil record's own survivorship bias.
2. **Conceptual classifications relative to a model.** Banzhaf et al. 2016 sort novelty into *variation*, *innovation* and *emergence*: novelty within the model, novelty that changes the model, and novelty that changes the meta-model. Taylor 2019 renames these *exploratory*, *expansive* and *transformational*. The York and Tokyo workshop categories (Taylor et al. 2016; Packard et al. 2019) sit alongside. Stepney 2021 and Stepney & Hickinbotham 2024 draw the consequence: a measure defined inside a model can see only variation.
3. **Formal dynamical-systems criteria.** Adams et al. 2017 define unbounded evolution and innovation against the counterfactual trajectories of an isolated system.
4. **Complexity growth.** Standish 2003 uses information measured against an equivalence class. McShea 1994 gives tests that tell a *passive* trend from a *driven* one. McShea 2001 gives a hierarchy scale of nesting depth, the only real-data anchor for "levels over time".

**The field agrees on four points:**
- OEE is plural. Hallmarks, meaning observed behaviour, must be kept apart from mechanisms (Taylor et al. 2016).
- A rising measure is necessary but not sufficient. A trivial string system passes every proposed criterion (Hintze 2019).
- Nobody has a validated measure of major transitions. The MODES authors explicitly "would welcome a measurement of a system's potential to produce major transitions in individuality".
- The standard empirical anchor for "unbounded" activity is the Phanerozoic fossil record, which Bedau et al. 1998 classify as class 3. Echo, the artificial system they classified, fell into class 1 or 2. Channon (2006) later classified Geb as unbounded under a modified normalisation.

**For Marco's hierarchy question the key finding is negative.** No standard OEE statistic registers a new level by itself. Every one of them counts or weights types at a level the analyst has already declared. Banzhaf et al. make a new level, by definition, a change to the meta-model. Such a change can only be caught by a recognizer: one written outside the simulation, one pre-coded into it, or one the simulation evolves itself (their §6.5). So Chemart has to build the level-1 recognizer and apply the same statistics one level up. Two nulls are needed:
- **A neutral shadow.** In Sayama's biased-random control, the count of novel higher-order entities still rose without selection.
- **Random unions of level-1 units.** Kim et al. 2019 used random merges of genomes, and those merges were distinguishable from real ecosystems.

**"Hierarchy" has three meanings in this literature, and only one of them is Marco's:**
- **Flow hierarchy**, i.e. feed-forwardness. Examples are Luo–Magee's measure (Chemart's `flow_hierarchy`) and Corominas-Murtra et al.'s treeness, feedforwardness and orderability. On that axis, metabolic networks sit inside the cloud of random graphs.
- **Hierarchical modularity**: nested modules within one network, e.g. Ravasz et al.'s C(k) ~ k⁻¹ and Holme et al.'s subnetwork shells.
- **Compositional (nested) hierarchy of individuals**: McShea's levels, Banzhaf's level-N entities, Rasmussen's dynamical hierarchies, and Fontana & Buss's levels 0/1/2.

The first two are properties of a single network. They do not imply the third.

**"Complex" also has several meanings, which the measures below keep apart.** It can mean size (the number of parts), information relative to an equivalence class (Standish; MODES's informative sites), diversity or ecology (Shannon entropy of persistent types), hierarchical depth (McShea), or construction depth of a molecule (the assembly index; Marshall et al. 2021).

**Empirical anchors.** Generic topology does not cleanly separate living from non-living chemistry:
- The Earth's atmosphere, the interstellar medium, lifeless atmospheres and metabolism are all small worlds.
- Scale-free claims fail severe tests (Broido & Clauset 2019; Smith et al. 2021; Wong et al. 2023).
- Tracing carbon atoms instead of letting currency metabolites act as shortcuts makes the *E. coli* "small world" not small (Arita 2004).

What does separate them is systematic departure from size-matched null ensembles:
- Earth's atmospheric network is "the most nonrandom" (Wong et al. 2023) and hierarchically modular (Solé & Munteanu 2004).
- Biochemical scaling laws are not reproduced by random samples of the same reaction universe (Kim et al. 2019).
- At the molecule level, an assembly index of about 15 or more was proposed as a biosignature (Marshall et al. 2021), but abiotic minerals reach 21 (Hazen et al. 2024).

Across the levels individual → ecosystem → biosphere, network topology shows *no sharp transition* (Smith et al. 2021). The levels differ only in how some measures scale with size, relative to null merges (Kim et al. 2019).

---

## 2. OEE and complexity-growth measures

Each entry gives the definition (quoted where precision matters), the inputs, whether the measure needs types, fitness or persistence, whether it can run on Chemart trajectories, and whether it detects a new level.

### 2.1 Evolutionary activity statistics (Bedau & Packard 1992; Bedau, Snyder & Packard 1998)

**Idea.** The 1992 paper defines evolutionary activity "as the rate at which useful genetic innovations are absorbed into the population". Activity counters are attached to components, and "an innovation 'make[s] a difference' if it persists and continues to be used". In the original model the counters were gene *usage* counters.

**Definitions (1998, the version used for classification):**
- activity increment δᵢ(t) = 1 if component i exists at t (age counting); aᵢ(t) = Σ_{k≤t} δᵢ(k)
- diversity D(t) = #{i : aᵢ(t) > 0}
- total cumulative activity A_cum(t) = Σᵢ aᵢ(t); mean cumulative activity Ā_cum = A_cum / D
- new activity A_new(t) = (1/D(t)) Σ_{i: a₀ ≤ aᵢ(t) ≤ a₁} aᵢ(t). The window [a₀, a₁] holds "the lowest activity values that can be interpreted as evidence that a component has positive adaptive significance", and it is set where the Echo and neutral-shadow activity distributions cross.

**Classes of long-term dynamics (verbatim):**
- "Class 1. No adaptive evolutionary activity: diversity D is bounded, new activity A_new is zero, and mean activity Ā_cum is zero."
- "Class 2. Bounded adaptive evolutionary activity: diversity D is bounded, new activity A_new is positive, and mean activity Ā_cum is bounded."
- "Class 3. Unbounded adaptive evolutionary activity: diversity D is unbounded, new activity A_new is positive, and mean activity Ā_cum is bounded."

**Result.** The Phanerozoic fossil record, with taxonomic families as components, is class 3. Holland's Echo is class 1 or 2, and all neutral shadows are class 1. The authors conclude: "Echo lacks the unbounded growth in adaptive activity observed in the fossil record".

**Neutral shadow.** Births, deaths and mutations mirror the real run, but "all selection in the shadow model is random". A "shadow child inherits its parent's genotype unless a mutation gives the child a new, unique genotype". The activity data are then normalised against the shadow.

- **Inputs.** Component presence (or usage) over time, plus a shadow run or another neutrality criterion. For the fossil record the authors argue that appearing in the record is itself evidence of adaptive persistence.
- **Needs types?** Yes, a declared component. The authors note that "a system could exhibit different classes of evolutionary dynamics at different levels of analysis".
- **Needs fitness or persistence?** No fitness. It does need a neutral baseline, the shadow, to set a₀.
- **Runs on Chemart?** Yes. D and Ā_cum come straight from `state`. A shadow can be driven by `fired`: every product event is a birth and every removal a death. *My suggestion:* give the shadow product the identity of a randomly chosen current shadow molecule, unless the real product was never seen before, in which case it gets a fresh label. That keeps population size, turnover and the novelty flux, and removes any identity-dependent advantage in production. Use `turnover()` as the generation clock, so that runs are comparable across chemistries.
- **Detects a new level?** Only if the components are defined at that level. Run with species as components, a level-2 organisation shows up at most indirectly, as long-lived member species.

### 2.2 Channon's refinements (2006) and complexity scaling in Geb (2019)

**2006.** Channon applied the Bedau test to Geb, an agent world with neural-network controllers: "Geb exhibits unbounded evolutionary dynamics, making it the first autonomous artificial system to pass this test". He then criticised the test "most significantly with regard to its normalization method for artificial systems" and proposed "component activity normalisation". MODES (Table 1 and §2.1) credits his related 2001 and 2003 conference papers, not listed separately here, with two further changes. First, he added classes 4b (unbounded per-component activity with bounded diversity) and 4c (both unbounded) beside 4a, which is Bedau's original unbounded class. Second, he advocated the *median* rather than the mean of per-component activity, because a single component under stabilising selection inflates the mean indefinitely.
- **Inputs and verdict:** as in §2.1.

**2019 ("Maximum individual complexity is indefinitely scalable in Geb").** Geb's individuals are neural networks. "Maximum individual complexity is found to be asymptotically bounded when scaling either parameter alone" (world length, which bounds population size, or the maximum number of neurons). It is "indefinitely scalable … when scaling both … together". It "scale[s] logarithmically with (the lower of) maximum population size and maximum number of neurons per individual".
- **Relevance to Chemart.** This is the closest published anchor to the neural-network-molecule idea. Molecule capacity and population size both have to scale, and growth may be only logarithmic. Design sweeps over both, and fit bounded against unbounded growth models, as in §2.3.

### 2.3 MODES toolbox (Dolson, Vostinar, Wiser & Ofria 2019)

MODES has four metrics, all computed on components that pass a filter.

**Persistence filter.** Tag lineages at time A. A component is "persistent" if it has descendants at A + t. The filter length t is measured in generations; coalescence theory guides it, since median neutral coalescence is 2N. Shadow runs "are also a viable filter option". MODES also reduces each genome to its *meaningful sites* (knockout or null-substitution lowers fitness), so that functionally identical genomes count once.

**The four metrics:**
- **Change**: change = Σ_{c∈F} [c ∉ F′], where F is the set of persistent components now and F′ the set at the previous timepoint.
- **Novelty**: novelty = Σ_{c∈F} [c ∉ S], where S is every component that has ever passed the filter. It is "functionally equivalent to A_new".
- **Complexity**: "the highest observed count of informative sites across all components in the population that make it through the filter".
- **Ecology**: ecology = −Σ_{c∈F} P(c) log₂ P(c). It is "equivalent to D", the Shannon entropy over persistent types.

**Boundedness.** Do not judge from a plot that seems to plateau. Instead "use statistics to determine what mathematical model best fits the observed data … [and] classify … based on the limit of the best-fitting mathematical model".

**Results.** In NK landscapes complexity saturates at N. In Avida, complexity keeps rising in the Logic-9 environment but falls in the empty one. Ecology is above baseline only under fitness sharing.

**Stated caveats:**
- Single-site knockouts miss epistasis, and "it may cause fragile genomes to appear more complex than robust ones".
- The metrics cannot separate class 1 from class 2, or class 3 from class 4b.
- There is no metric for major transitions.

- **Needs types?** Yes.
- **Needs fitness or persistence?** Persistence, which is lineage-based. The complexity metric also needs a fitness or function oracle for the knockouts.
- **Runs on Chemart?** Change, novelty and ecology: yes, with *species-level* persistence (present across frames spanning at least t turnovers). Chemart has no individual lineages. `fired` gives reactant→product provenance, but descendant closure in a well-mixed gas quickly covers everything. Complexity: only if Chemart defines an operational knockout, for example replacing a sub-term with a null term and checking whether the molecule's reactions against the current population change.
- **Relevance to the redundancy hypothesis.** Knockout-based complexity mechanically rates redundant, robust molecules as less complex. Test redundancy at the level of the organisation, not by molecular complexity.
- **Detects a new level?** No. Units at a higher level could be substituted, but nothing in MODES finds them.

### 2.4 QNN evolutionary activity (Droop & Hickinbotham 2012)

**What it is.** A quantitative, non-neutral activity measure. It needs no shadow run, works with intrinsic or extrinsic fitness, and was applied to Tierra and **Stringmol**, which is in Chemart's catalog.

**How it is computed** (from the authors' R package, `franticspider/qnn`, function `qnn.activity`). Take each species' proportion pᵢ(t) of the population. The expected proportion is last frame's value. Keep only positive deviations, max(0, pᵢ(t) − pᵢ(t−1)), square them, and sum over species and time. The paper's own formulae are paywalled; I checked the code, not the PDF.

- **Needs types?** Yes.
- **Needs fitness or persistence?** No. The neutral expectation is implicit: proportions stay unchanged.
- **Runs on Chemart?** Yes, directly on `traj.array()`, which is X[time, species].
- **Detects a new level?** No.

### 2.5 Novelty classifications: Banzhaf et al. 2016, Taylor 2019, York/Tokyo categories

**Banzhaf et al. 2016 (Theory in Biosciences).** This is the framework closest to the book Chemart is built on. Entities are organised into levels: "level-0 entities are atomic entities and level-N > 0 entities are system entities that contain at least two lower-level entities, of which at least one is level-N−1". The three types of novelty are:
- "**Variation: novelty within the model.** Variation is a change to an instance of the model, a change to the values of a variable that conforms to the model."
- "**Innovation: novelty that changes the model.** Innovation is a change to the model: a change that adds a new type or relationship that conforms to the meta-model, or possibly eliminates an existing one." Adding a new species to an ecosystem is their example.
- "**Emergence: novelty that changes the meta-model.** Emergence is a change to the meta-model: a change that adds a new meta-type or relationship". "One particularly important phenomenon is the emergence of a new level of organisation."

Two further points follow:
- **Only the first instance of a level is emergence.** "Once the first level-N+1 entity has emerged as a type-2 novelty, similar entities appearing later will be added to a pre-existing level. Thence, the next level-N+1 entity to appear is only a type-1 novelty".
- **What counts as OEE.** "Open-ended system: a system with the ability to continually produce open-ended events", where innovation and emergence are the open-ended events and variation is not.

The authors state that component measures such as Bedau's "could then potentially be used to quantify the respective amounts of variation, innovation, or emergence activity". Which of the three is being counted depends on the model in use.

They list five ways a simulation can capture novelty:
1. recognised outside the simulation;
2. "anticipated and recognized, with pre-coded recognizers";
3. anticipated and captured with pre-coded structural rules;
4. emergently recognised by self-modifying code;
5. emergently captured by self-modifying code.

They also discuss kinds of levels, citing Craver: mereology, aggregativity, spatial containment and mechanism. For "networks of networks", the relevant sense is *levels of mechanism*, where the components' behaviours are relevant to the behaviour of the whole. Mere aggregation does not count.

**Taylor 2019 (Artificial Life).** Taylor renames the three types:
- "Exploratory Novelty: A novelty that can be described using the current model";
- "Expansive Novelty: A novelty that necessitates a change in the model but still using concepts present in the current meta-model";
- "Transformational Novelty: A novelty that introduces a new concept, necessitating a change in the meta-model" (for example, a major transition).

Each "…Open-Endedness" is "the ongoing production of **adaptive** … novelties". He departs from Banzhaf in two ways. Novelty is defined relative to the *initial* model, so repeated transitions stay transformational. And ongoing exploratory novelty counts as a kind of OEE. His routes to expansive and transformational novelty are "multiple domains of behaviour, transdomain bridges, and non-additive compositional systems". His claim that a standard Darwinian analysis reaches only exploratory OEE is directly relevant to asking whether neural-network molecules add non-additive composition.

**York and Tokyo categories.** Taylor et al. 2016 stressed "pluralism about OEE" and "distinguishing observable behavioral hallmarks … from hypothesized underlying mechanisms". As reproduced in Packard et al. 2019, the York list is:
1. ongoing generation of adaptive novelty: (a) new adaptations, (b) new kinds of entities, (c) emergence of major transitions, (d) evolution of evolvability;
2. ongoing growth of complexity: (a) entity complexity, (b) interaction complexity.

The Tokyo revision is: 1. interesting new kinds of entities and interactions; 2. evolution of evolvability; 3. major transitions; 4. semantic evolution. Packard et al. call emergent hierarchical structure (dynamical hierarchies) a "holy grail" for ALife models.

- **Needs types?** The whole scheme is types relative to an explicit model and meta-model.
- **Needs fitness?** Taylor requires novelties to be adaptive; Banzhaf does not.
- **Runs on Chemart?** Not as a number. It is a classification that tells you what your counts mean.
- **Detects a new level?** Yes, by definition (type 2, or transformational), but only through a recognizer.

### 2.6 Unbounded evolution and innovation (Adams, Zenil, Davies & Walker 2017)

The system is split into an "organism" subsystem o and an environment e.
- "**Definition 1** Unbounded evolution (UE): A system U that can be decomposed into two interacting subsystems o and e, exhibits unbounded evolution if there exists a recurrence time such that the state-trajectory or the rule-trajectory of o is non-repeating for t_r > t_P or t_r′ > t_P … where t_P is the Poincaré recurrence time for an equivalent isolated (non-perturbed) system o."
- "**Definition 2** Innovation (INN): … exhibits innovation if there exists a recurrence time t_r such that the state-trajectory is not contained in the set of all possible state trajectories for an equivalent isolated (non-perturbed) system."

In their cellular-automaton variants, "state-dependent dynamics … statistically out-performs other candidate mechanisms, and is the only mechanism to produce open-ended evolution in a scalable manner".

- **Inputs.** State and rule trajectories, plus the counterfactual trajectory set of isolated fixed-rule systems.
- **Needs types?** No.
- **Needs fitness?** No, but it needs the counterfactual.
- **Runs on Chemart?** Only as a coarse proxy, such as the recurrence of species-set states over windows, because t_P is astronomical for real gases. The conceptual link is useful: in AlChemy and BFF the molecules are also the rules, which is exactly their "state-dependent" case.
- **Detects a new level?** No.

### 2.7 Cardinality leap and Hash Chemistry (Sayama 2019)

**Claim.** "Facilitating formation of higher-order entities is a generalizable, effective way to cause a 'cardinality leap' in the set of possibilities". A finite possibility set becomes countably infinite, and a countably infinite one becomes uncountable.

**Measures:**
- the cumulative number of unique entity *types*, counted separately for individual entities (which saturate at 1,000) and for higher-order entities that replicated (which grow "almost linearly … without an apparent bound");
- the maximum and average **number of individual entities involved in each replication event**, which rose over time.

Growth is judged by fitting bounded and unbounded models.

**Controls.** With a fair random-fitness control, populations went extinct. With a biased random control, "the number of novel higher-order entities did increase with time … indicating that this novelty production effect of the cardinality leap would not require selection or adaptation per se". Only the quantitative difference from the main runs supports adaptation.

- **Needs types?** Yes, at both levels.
- **Needs fitness or a replication criterion?** Yes, a replication criterion, plus a random-fitness control.
- **Runs on Chemart?** Partly. Replication events have to be identified in `fired`, for example events in which a species is produced while also acting as a reactant or operator, or the collective production of an autocatalytic set within a window.
- **Detects a new level?** Closest of the quantitative measures. The size of the unit of replication is a direct proxy for higher-order entities. But higher-order novelty also rises under a biased neutral control, so it must always be compared with a shadow.

### 2.8 Critiques and necessary conditions

**Hintze 2019.** "A simple evolving computational system that satisfies all such requirements is presented. Doing so reveals a shortcoming in the definitions". The system is strings over {R, L, F} drawn as LOGO-like paths, with fitness that rewards novel paths. Its diversity (Levenshtein distance) and complexity (compression) both keep rising. The lesson: "our definitions of complexity and diversity are more important … than the fact that something runs indefinitely".

**Stepney 2021 (OEE4 workshop).** Standard measures "are measuring 'Flatland' … they are confined within a model … and cannot see outside. They are measuring type-0 variational change". She proposes looking for *discontinuities* in measures, like order parameters at phase transitions. She also proposes looking for a measure that *diverges* when the model becomes insufficient: Crutchfield's statistical complexity diverges when a stack must be added to the model. And she proposes casting models and meta-models in a formal (e.g. graph) language, so that model change itself can be measured.

**Stepney & Hickinbotham 2024 (Artificial Life).** "Attempting to quantify open-endedness misses the point: … an open-ended system will eventually move outside its current model of behavior, and hence outside any measure based on that model". They apply generic and then system-specific measures to spatial **Stringmol** runs with a parasite arms race, and recommend exactly that sequence.

**Soros & Stanley 2014 (Chromaria).** They propose four necessary conditions, quoted here as listed in Hintze 2019:
1. individuals must meet a nontrivial "minimal criterion (MC) before they can reproduce";
2. "the evolution of new individuals should create novel opportunities for satisfying the MC";
3. "decisions about how and where individuals interact with the world should be made by the individuals themselves";
4. "the potential size and complexity of the individuals' phenotypes should be (in principle) unbounded".

Chromaria stagnates when any one of them is removed. These are mechanisms, not measures, but they give Chemart a checklist of design fields.

### 2.9 Complexity growth and levels

**Standish 2003.** Standish writes: "Most people equate open-ended evolution with complexity growth, although a priori these seem to be different things". He defines complexity as C(x) = lim_{s→∞} [s·log₂N − log₂ω(s,x)], where ω(s,x) is the number of descriptions of length s that are equivalent to x under an interpreter, i.e. the neutral class. Random strings have zero complexity. In a size-neutral Tierra run, "no increase in organismal complexity was observed, although organism size did increase".
- **Relevance to Chemart.** (i) Chemart's `complexity_drift`, the slope of mean species-string length, measures size, not complexity. (ii) Under this measure, a redundant code (large ω) is *less* complex per unit length. Test the redundancy hypothesis on organisational complexity, not molecular complexity.
- **Runs on Chemart?** Only if ω can be estimated, for example by sampling variants of a molecule that react the same way against the current population.

**McShea 1994: passive versus driven trends.** "In a driven trend, the distribution mean increases on account of a force … In a passive system … the mean increases because change in one direction is blocked by a boundary". The tests: "a system is driven if the minimum increases, if increases significantly outnumber decreases among ancestor-descendant pairs [far from the lower bound], and if the mean skew of subclades is significantly positive".
- **Runs on Chemart?** Yes, and directly. Molecules have a minimum length, so a rising mean length can be purely passive. Track the minimum or low quantiles of length. Use `fired` to form reactant→product pairs as the ancestor–descendant sample, and compare the share of events whose product is longer than its reactants with the share where it is shorter, weighted by count and excluding pairs near the minimum size.

**McShea 2001: hierarchy scale.** "The degree of hierarchical structure of organisms—the number of levels of nesting of lower-level entities within higher-level individuals". He proposes "levels reflecting number of layers of nestedness, and sublevels reflecting degree of individuation at the highest level". Plotted against the fossil record, the scale shows "a long-term trend extending from the Archean through the early Phanerozoic" in the **maximum**.
- **Relevance to Chemart.** This is the template for a Chemart statistic of "maximum nesting depth over time". The sublevel idea (degree of individuation of the top level) gives graded evidence before a full new level appears.

**Rasmussen, Baas, Mayer & Nilsson 2001: dynamical hierarchies.** They build a framework to "simultaneously investigate different levels of description together with their interrelationship". In a lattice-gas molecular-dynamics model, three levels emerge (monomers, polymers, micelles). They close with "a conjecture of a necessary minimal complexity within the fundamental interacting structures". This is a published antecedent of the redundancy hypothesis: the properties of the lowest-level objects limit which higher levels can emerge.

**Krakauer, Bertschinger, Olbrich, Flack & Ay 2020: information theory of individuality.** "Individuals are aggregates that preserve a measure of temporal integrity, i.e., 'propagate' information from their past into their futures." There are three forms: organismal, colonial and driven. The approach "allows for the identification of individuals at all levels of organization". It needs no predefined types, only discretised time series of a candidate aggregate and its environment. It is the most principled detector of a new level, and also the most expensive.

**Level-1 units in chemistry.** A chemical organisation is "a closed and self-maintaining set of components", and stationary states are instances of organisations (Dittrich & Speroni di Fenizio 2007). Organisations are Chemart's natural level-1 entities. Note that Chemart's `organisations` measure has exponential cost (limit 200 species).

### 2.10 Summary table

| Measure | Inputs | Needs types? | Needs fitness / persistence / null? | Runs on Chemart trajectories? | Registers a new level? |
|---|---|---|---|---|---|
| Activity statistics D, A_new, Ā_cum; classes 1–3 (Bedau 1992/1998) | component presence or usage over time | yes | neutral shadow (or survivorship) | yes; shadow built from `fired` | only if components are level-N units |
| Channon normalisation; classes 4a–c; median activity | same | yes | shadow | yes | no |
| MODES change / novelty / ecology | persistent components | yes | persistence filter (lineage) or shadow | yes, with species persistence in turnovers | no |
| MODES complexity (informative sites) | knockouts with fitness effect | yes | fitness or function oracle | only with an operational molecular knockout | no; rates robust, redundant molecules as less complex |
| QNN | species × time abundances | yes | none (implicit neutral expectation) | yes, `traj.array()` | no |
| Constructiveness: count of species ever seen is increasing (Dittrich et al. 2001) | species sets | yes | none | yes (the integral of `novelty_rate`) | no |
| Variation / innovation / emergence (Banzhaf 2016); exploratory / expansive / transformational (Taylor 2019) | explicit model and meta-model | model types | Taylor: adaptive | as a classification only | yes, by definition, via a recognizer |
| UE and INN (Adams 2017) | state and rule trajectories; isolated-system counterfactuals | no | counterfactual | coarse proxy only | no |
| Cardinality leap; replication-unit size (Sayama 2019) | replication events and their participants | yes, two levels | replication criterion and random control | partly | **yes (closest)** |
| Information complexity (Standish 2003) | description plus equivalence classes | yes | interpreter / neutral set | if ω can be sampled | no; a larger neutral set (more redundancy) lowers C |
| Passive vs driven trend tests (McShea 1994) | trait distribution; ancestor–descendant pairs | a trait | none | **yes** (`fired` pairs) | no (applies at any level) |
| Hierarchy scale, trend in the maximum (McShea 2001) | criteria for individuality | per level | criteria | via a recognizer | **yes** |
| Information individuality (Krakauer 2020) | time series of aggregate and environment | no | none | expensive | **yes** |

### 2.11 What this means for Chemart (suggestions, not from the literature unless cited)

1. **Harden the existing measures.**
   - **`novelty_rate`.** Count a species as novel only once it persists for at least t turnovers (MODES), and report it next to its value from a neutral shadow (Bedau).
   - **`complexity_drift`.** Add the minimum and low-quantile trajectory and an ancestor–descendant bias from `fired` (McShea 1994). Relabel it as size drift (Standish 2003).
   - **Boundedness.** Fit bounded against unbounded models; do not read plateaus by eye (MODES §3.1.3).
2. **Add** activity statistics with a shadow, QNN, and MODES change and ecology. All three run on existing trajectories.
3. **A level pipeline for the hierarchy question:**
   - (a) Identify level-1 units in `traj.window(i, w)`: organisations, irreducible RAFs, or self-maintaining closed sets. Match them over time by Jaccard similarity of their species sets.
   - (b) Rerun (1–2) with these units as components. This gives innovation at level 1.
   - (c) Call a candidate level-2 entity a persistent, self-maintaining set that contains at least two level-1 units, with cross-reactions that both require. This follows Banzhaf's definition of level-N and Fontana & Buss's level 2.
   - (d) Report the maximum nesting depth over time (McShea 2001) and the number of level-1 units per self-maintaining unit (in the spirit of Sayama).
   - (e) Test against two nulls: a neutral shadow, and random unions of co-occurring level-1 units, which is the logic of Kim et al.'s random-genome networks.

   The claim in Mathis et al. 2024 that AlChemy organisations "cannot be easily combined into higher order entities" then becomes a measurable frequency: how often a union of two organisations remains an organisation, compared with the rate under the null.

---

## 3. Empirical anchors: what these measures show on real chemistry and biochemistry

| Study | System | Measure(s) | Finding | Separates biotic from abiotic? |
|---|---|---|---|---|
| Jeong et al. 2000 | metabolic networks of 43 organisms, all three domains | degree distribution; diameter against size; hub identity; attack and error tolerance | γ_in = γ_out ≈ 2.2; the diameter "is the same for all 43 organisms, irrespective of the number of substrates"; hub ranking "practically identical"; only ~4% of substrates are present in all species | No abiotic comparator; the authors stress similarity to *non*-biological networks |
| Wagner & Fell 2001 | *E. coli* metabolism | small-world indices; degree distribution; centrality | small world; power-law metabolite connectivity; tricarboxylic acid (TCA) cycle central | No comparator |
| Arita 2004 | *E. coli*, carbon-atom-traced pathways | path length with conserved structural moieties | "the average path length … is much longer than previously thought and … the metabolic world of this organism is not small" | **Pitfall:** currency metabolites create shortcuts, so representation decides the result |
| Ravasz et al. 2002 | 43 organisms | ⟨C⟩ against N; C(k) | ⟨C⟩ is about an order of magnitude above a scale-free model of the same size and independent of size; C(k) ~ k⁻¹ (hierarchical modularity) | No comparator (see the Earth's atmosphere below) |
| Holme, Huss & Jeong 2003 | 43 organisms, metabolic and whole-cell networks | recursive betweenness-based decomposition | "a few core-clusters centred around the most highly connected substances enclosed by other substances in outer shells" | No |
| Solé & Munteanu 2004 | atmospheres of Earth, Mars, Venus and Titan; giant-planet hydrocarbons; interstellar medium; *E. coli* | ⟨L⟩ and ⟨C⟩ against random; assortativity r; P(k); modularity | All are small worlds. Earth: V = 248, ⟨C⟩ = 0.31 against 0.025 random, r = −0.31, scale-free (γ ≈ 2.16), modular. Mars, Venus, Titan and the hydrocarbon network: single-scale, not modular. Interstellar medium (V = 400): broad-scale, not modular. *E. coli*: ⟨C⟩ = 0.183 against 0.008, modular | **Suggestive yes** (the living planet against dead ones), but with n = 1 biosphere and small curated mechanisms |
| Wong et al. 2023 | updated atmospheric networks of Solar System bodies | global, centrality, community and cluster analysis | Earth is *not* scale-free in the updated models, but it has "the most nonrandom topology" and in some metrics is "more similar to biological networks" | **Partly:** yes by non-randomness, no by scale-freeness |
| Kim et al. 2019 | 28,146 genomes and metagenomes; 8,658 KEGG reactions; individuals, ecosystems, biosphere | scaling of ⟨L⟩, ⟨C⟩, betweenness, assortativity, #reactions, #edges and #EC classes with size; three null ensembles; unipartite and bipartite representations | Scaling laws are shared across levels. Random reaction networks (flat samples of KEGG) do not reproduce the scaling of clustering and assortativity. Levels are distinguishable in #reactions, #edges, #EC and clustering (P < 10⁻⁵ in most cases), but not in betweenness or path length. Assortativity separates real ecosystems from random merges of genomes, so ecosystems are not mere unions. The biosphere is not predicted by ecosystem scaling. Domain of life is predicted with >80% accuracy | Separates life from *randomly assembled biochemistry* (not from abiotic chemistry) |
| Smith, Kim & Walker 2021 | 785 metagenomes and 1,082 genomes | Broido–Clauset scale-free tests; individual against ecosystem | "no more than a few biochemical networks are any more than super-weakly scale-free"; "no sharp transition in the structure of biochemical networks across these levels" | N/A; the levels do not differ in kind topologically |
| Gagler et al. 2022 | 11,955 metagenomes; 1,282 archaea; 11,759 bacteria; 200 eukaryotes | scaling of enzyme counts per EC class with total enzyme count | oxidoreductases and hydrolases superlinear; transferases and ligases sublinear across all datasets; "not explained by … components shared across known examples of life"; used to predict LUCA | Within life only |
| Goldford et al. 2017 | KEGG reaction set; network expansion from prebiotic seed compounds | scope (network expansion) | a "phosphate-independent core metabolism", about 315 reactions and 260 metabolites (Boston University coverage), enriched for iron–sulphur enzymes; thioesters overcome thermodynamic bottlenecks | N/A; a method anchor for Chemart's `scope_fraction` and `expansion_depth` |
| Xavier et al. 2020 | O₂-independent prokaryotic network (5,994 reactions, 5,723 metabolites); an acetogen; a methanogen | maxRAF size against food set | maxRAF by food set: inorganic only, 8 reactions; plus formate, methanol, acetate and pyruvate, 16; plus 8 LUCA cofactors, 914; all cofactors, 1,335 ("25% of the starting anaerobic network"). Acetogen 394, methanogen 209, intersection 172. Amino acids and bases without organic catalysts give 33. Removing NAD roughly halves the maxRAF; removing ATP has no effect | N/A; RAF size is set by small-molecule catalysts, which calibrates `max_raf_fraction` |
| Wołos et al. 2020 | Allchemy forward synthesis from water, N₂, NH₃, HCN, CH₄ and H₂S | generation-wise network; emergent catalysis, cycles, surfactants | "three forms of nontrivial chemical emergence": catalysts of downstream reaction types, self-regenerating cycles (the iminodiacetic acid cycle, validated experimentally), and surfactants | Reported: biotic molecules are more hydrophilic, more thermodynamically stable and more balanced in H-bond donors and acceptors. *This detail comes from a search-engine rendering of the Science summary; I could not read the full text.* |
| Robinson et al. 2022 | formose network in a flow reactor | product composition against environment | "the compositional complexity of the reaction products [can] be controlled as a function of factors such as feedstock and catalyst availability"; Breslow's cycle feeds C₂ building blocks | Abiotic only; no network-topology comparison found |
| Saylam, Hadj Ali & Fikri 2020 | combustion mechanisms | degree centrality | identifies "principal species" for mechanism reduction | None; I found no study comparing combustion-network topology with biological networks |
| Marshall et al. 2021 | molecules in diverse samples, measured by MS² | molecular assembly index (MA) | "only biologically produced samples produce MA above a certain threshold" (about 15); MA 15–20 has chance formation of about 1 in 10²³ | **Claimed yes** at the molecule level |
| Hazen et al. 2024 | mineral heteropolyanions | MA | MA ranges from 2 up to 21; "values of molecular assembly indices ≥15 do not represent unambiguous biosignatures" | Counterexample to the above |
| Bedau, Snyder & Packard 1998 | fossil record at family level; Echo | activity statistics | biosphere class 3; Echo class 1 or 2 | Separates the biosphere from an ALife model |
| McShea 2001 | body-fossil record | hierarchy scale (nesting plus individuation) | long-term trend in the *maximum* level from the Archean to the early Phanerozoic | The only real-data anchor for "levels over time" |
| Corominas-Murtra et al. 2013 | 125 real networks in 13 classes | treeness, feedforwardness, orderability | metabolic, neural, linguistic and some social networks sit "clearly embedded within the cloud of random graphs" | Flow hierarchy is **not** a biotic signature |

### 3.1 How to read these anchors for Chemart

- **Raw topology is a weak discriminator.** Small-world structure, heavy tails, clustering and flow hierarchy show up in abiotic chemistry (the interstellar medium, lifeless atmospheres), in engineered systems and in random graphs. Claims about power laws often fail severe tests (Broido & Clauset 2019; Smith et al. 2021; Wong et al. 2023). And conclusions flip with representation (Arita 2004). Kim et al. 2019 therefore compute both unipartite and bipartite projections.
- **Deviation from matched nulls is a better discriminator.** The discriminating studies all ask whether a network departs from a null of the same size built from the same parts: Kim's flat and frequency-weighted random reaction networks and random merges of genomes, Wong's "most nonrandom", and Solé & Munteanu's comparison against random graphs. Chemart can apply the same test to each chemistry, using as the null random networks sampled from that chemistry's own reaction universe, at matched size.
- **Levels in real biochemistry.** Individual → ecosystem → biosphere is a real-data analogue of molecules → networks → networks of networks, one step up. It shows up as scaling differences and as departures from random-merge nulls, not as a qualitative change in topology (Kim 2019; Smith 2021). Transferred to Chemart: a candidate level-2 organisation should be tested against random unions of level-1 organisations of the same number and size.
- **Calibration points Chemart can compute through its own pipeline:**
  - *E. coli* metabolism and Earth against Mars and Venus atmospheres (Solé & Munteanu's table gives targets, but recompute them in Chemart's representation);
  - KEGG ensembles (Kim et al.);
  - maxRAF fraction against food set (Xavier et al.);
  - scope from prebiotic seed compounds (Goldford et al.);
  - assembly index on string-molecule analogues, used with Hazen's caveat.
- **Gaps.** I found no network-topology study of the formose network that compares it with biology; Robinson et al. is compositional. For combustion, network measures are used only for mechanism reduction.

---

## 4. Morphospace method notes and pitfalls

**Classic idea.** A *theoretical* morphospace is generated by a model of form and then compared with where real forms fall. Raup's coiled shells (1966) are the standard example; McGhee 1999 systematised the approach. Most of the generated space is empty. Avena-Koenigsberger et al. 2015, citing Wright and McGhee: "most empirical morphospaces are mostly empty".

**For networks:**
- **Network morphospace.** Avena-Koenigsberger et al. 2015 use axes that are network traits. They separate regions that are geometrically impossible, functionally impossible, and possible but unrealised, and they explore the space by rewiring under constraints.
- **Hierarchy morphospace.** Corominas-Murtra et al. 2013 build a 3D space of treeness, feedforwardness and orderability, with random-model clouds as the reference.
- **Pareto archetypes.** Shoval et al. 2012: under trade-offs, "best-trade-off phenotypes are weighted averages of archetypes". With two tasks, phenotypes lie on a line; with three, in a triangle whose vertices are the archetypes.

**Design axes for artificial chemistries.** Dittrich, Ziegler & Banzhaf 2001: "an artificial chemistry can be defined by a triple (S, R, A), where S is the set of all possible molecules, R is a set of collision rules …, and A is an algorithm describing the reaction vessel or domain and how the rules are applied". Their characteristics give natural categorical axes: molecules, reaction laws and dynamics each explicit or implicit; analogous or abstract; constructive (weakly: new components generated randomly; strongly: generated by the action of other components); random chemistries; the measure of time; pattern matching; spatial topology. They also give an operational test for constructiveness: "the system can be considered constructive if the number of elements that appeared in the system at least once increases in time".

**Pitfalls:**
1. **Axis choice is hypothesis-laden.** "Constructing a morphospace involves a choice of structural features to analyse and a mathematical model to measure them … made with respect to the hypothesis under consideration" (Avena-Koenigsberger et al.). Many network measures are strongly correlated, especially through size and density, so decorrelate them or choose a small set on purpose.
2. **Size and density effects.** "Graph measures may be influenced by the number of nodes (N) and the average degree (k) … Direct comparisons … can therefore yield spurious results". Normalising by random surrogates "may even increase the sensitivity to differences in N and k" for clustering and the small-world index, and "none of the here-investigated methods allows for a reliable and fully unbiased comparison" (van Wijk et al. 2010). Prefer scaling fits across sweeps, as Kim et al. do and as Chemart's `measures.scaling` already supports.
3. **Null models.** Place a cloud of random or null networks in the space and interpret only departures from it. Metabolic networks sit inside the random cloud on hierarchy (Corominas-Murtra et al.). For Chemart, the nulls should be random reaction networks from the same chemistry's reaction universe, as in Kim et al.
4. **Occupied against empty regions.** Emptiness can be geometric, functional, historical (unexplored), or a sampling artefact. The catalog is the book's selection, not a random sample of artificial chemistries, and many entries are variants of one another. *My note:* treat related chemistries as non-independent points, as comparative biology does with phylogeny.
5. **Metric structure.** Raup's space "does not possess a Euclidean structure and a meaningful interpretation of the spread and spacing of taxa within it is not guaranteed" (Gerber 2017). Distances and disparity mean something only when the axes are commensurable. For mixed catalog fields, choose the distance on purpose.
6. **Representation.** Unipartite against bipartite, and with or without currency species, can reverse conclusions (Arita 2004; Kim et al. 2019). Check that each axis is robust across projections.
7. **Trade-off geometry.** If behavioural-space occupancy forms a polytope, its vertices suggest the "tasks" (Shoval et al.). The inference assumes the points were shaped by optimisation, which is doubtful for a catalog of human designs.

---

## 5. References

**Open-endedness and complexity growth**

1. Bedau, M. A., & Packard, N. H. (1992). Measurement of evolutionary activity, teleology, and life. In C. Langton, C. Taylor, D. Farmer & S. Rasmussen (Eds.), *Artificial Life II* (pp. 431–461). Addison-Wesley. https://people.reed.edu/~mab/publications/papers/alife2.pdf. **VERIFIED** (full text read).
2. Bedau, M. A., Snyder, E., & Packard, N. H. (1998). A classification of long-term evolutionary dynamics. In C. Adami et al. (Eds.), *Artificial Life VI* (pp. 228–237). MIT Press. https://people.reed.edu/~mab/publications/papers/alife6.pdf. **VERIFIED** (full text read; the preprint header also states pp. 189–198).
3. Channon, A. (2006). Unbounded evolutionary dynamics in a system of agents that actively process and transform their environment. *Genetic Programming and Evolvable Machines*, 7(3), 253–281. https://doi.org/10.1007/s10710-006-9009-3. **VERIFIED** (metadata and abstract).
4. Channon, A. (2019). Maximum individual complexity is indefinitely scalable in Geb. *Artificial Life*, 25(2), 134–144. https://doi.org/10.1162/artl_a_00285. **VERIFIED** (abstract).
5. Droop, A., & Hickinbotham, S. (2012). A quantitative measure of non-neutral evolutionary activity for systems that exhibit intrinsic fitness. In *Artificial Life 13* (pp. 45–52). MIT Press. https://doi.org/10.7551/978-0-262-31050-5-ch007. Code: https://github.com/franticspider/qnn. **VERIFIED** (metadata; computation read from the authors' code, not the paper).
6. Dolson, E. L., Vostinar, A. E., Wiser, M. J., & Ofria, C. (2019). The MODES toolbox: Measurements of open-ended dynamics in evolving systems. *Artificial Life*, 25(1), 50–73. https://doi.org/10.1162/artl_a_00280. **VERIFIED** (full text read).
7. Taylor, T., Bedau, M., Channon, A., Ackley, D., Banzhaf, W., Beslon, G., … Wiser, M. (2016). Open-ended evolution: Perspectives from the OEE workshop in York. *Artificial Life*, 22(3), 408–423. https://doi.org/10.1162/ARTL_a_00210. **VERIFIED** (abstract; York categories via ref. 8).
8. Packard, N., Bedau, M. A., Channon, A., Ikegami, T., Rasmussen, S., Stanley, K. O., & Taylor, T. (2019). An overview of open-ended evolution: Editorial introduction to the Open-Ended Evolution II special issue. *Artificial Life*, 25(2), 93–103. https://doi.org/10.1162/artl_a_00291 (arXiv:1909.04430). **VERIFIED** (full text read). A companion three-page OEE I editorial also exists: *Artificial Life* 25(1), 1–3, https://doi.org/10.1162/artl_e_00282.
9. Banzhaf, W., Baumgaertner, B., Beslon, G., Doursat, R., Foster, J. A., McMullin, B., de Melo, V. V., Miconi, T., Spector, L., Stepney, S., & White, R. (2016). Defining and simulating open-ended novelty: Requirements, guidelines, and challenges. *Theory in Biosciences*, 135(3), 131–161. https://doi.org/10.1007/s12064-016-0229-7. **VERIFIED** (author manuscript read, from the MMU repository).
10. Taylor, T. (2019). Evolutionary innovations and where to find them: Routes to open-ended evolution in natural and artificial systems. *Artificial Life*, 25(2), 207–224. https://doi.org/10.1162/artl_a_00290 (arXiv:1806.01883). **VERIFIED** (full text read).
11. Adams, A., Zenil, H., Davies, P. C. W., & Walker, S. I. (2017). Formal definitions of unbounded evolution and innovation reveal universal mechanisms for open-ended evolution in dynamical systems. *Scientific Reports*, 7, 997. https://doi.org/10.1038/s41598-017-00810-8. **VERIFIED** (full text read).
12. Sayama, H. (2019). Cardinality leap for open-ended evolution: Theoretical consideration and demonstration by Hash Chemistry. *Artificial Life*, 25(2), 104–116. https://doi.org/10.1162/artl_a_00283. **VERIFIED** (full text read).
13. Hintze, A. (2019). Open-endedness for the sake of open-endedness. *Artificial Life*, 25(2), 198–206. https://doi.org/10.1162/artl_a_00289. **VERIFIED** (full text read).
14. Stepney, S. (2021). Modelling and measuring open-endedness. OEE4: Fourth Workshop on Open-Ended Evolution, ALIFE 2021 (online). http://workshops.alife.org/oee4/papers/stepney-oee4-camera-ready.pdf. **VERIFIED** (full text read).
15. Stepney, S., & Hickinbotham, S. (2024). On the open-endedness of detecting open-endedness. *Artificial Life*, 30(3), 390–416. https://doi.org/10.1162/artl_a_00399. **VERIFIED** (abstract).
16. Soros, L. B., & Stanley, K. O. (2014). Identifying necessary conditions for open-ended evolution through the artificial life world of Chromaria. In *ALIFE 14* (pp. 793–800). MIT Press. https://doi.org/10.7551/978-0-262-32621-6-ch128. **VERIFIED** (metadata; the four conditions are quoted as listed in ref. 13).
17. Standish, R. K. (2003). Open-ended artificial evolution. *International Journal of Computational Intelligence and Applications*, 3(2), 167–175. https://doi.org/10.1142/S1469026803000914 (arXiv:nlin/0210027). **VERIFIED** (full text read).
18. McShea, D. W. (1994). Mechanisms of large-scale evolutionary trends. *Evolution*, 48(6), 1747–1763. https://doi.org/10.1111/j.1558-5646.1994.tb02211.x. **VERIFIED** (abstract).
19. McShea, D. W. (2001). The hierarchical structure of organisms: A scale and documentation of a trend in the maximum. *Paleobiology*, 27(2), 405–423. https://doi.org/10.1666/0094-8373(2001)027<0405:THSOOA>2.0.CO;2. **VERIFIED** (abstract).
20. Rasmussen, S., Baas, N. A., Mayer, B., & Nilsson, M. (2001). Ansatz for dynamical hierarchies. *Artificial Life*, 7(4), 329–353. https://doi.org/10.1162/106454601317296988. **VERIFIED** (abstract).
21. Krakauer, D., Bertschinger, N., Olbrich, E., Flack, J. C., & Ay, N. (2020). The information theory of individuality. *Theory in Biosciences*, 139(2), 209–223. https://doi.org/10.1007/s12064-020-00313-7. **VERIFIED** (abstract).
22. Dittrich, P., & Speroni di Fenizio, P. (2007). Chemical organisation theory. *Bulletin of Mathematical Biology*, 69(4), 1199–1231. https://doi.org/10.1007/s11538-006-9130-8. **VERIFIED** (abstract).
23. Mathis, C., Patel, D., Weimer, W., & Forrest, S. (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos*, 34(9), 093142. https://doi.org/10.1063/5.0207358. **VERIFIED** (abstract).

**Empirical anchors**

24. Jeong, H., Tombor, B., Albert, R., Oltvai, Z. N., & Barabási, A.-L. (2000). The large-scale organization of metabolic networks. *Nature*, 407, 651–654. https://doi.org/10.1038/35036627 (arXiv:cond-mat/0010278). **VERIFIED** (preprint read).
25. Wagner, A., & Fell, D. A. (2001). The small world inside large metabolic networks. *Proceedings of the Royal Society B*, 268(1478), 1803–1810. https://doi.org/10.1098/rspb.2001.1711. **VERIFIED** (abstract).
26. Arita, M. (2004). The metabolic world of *Escherichia coli* is not small. *PNAS*, 101(6), 1543–1547. https://doi.org/10.1073/pnas.0306458101. **VERIFIED** (abstract).
27. Ravasz, E., Somera, A. L., Mongru, D. A., Oltvai, Z. N., & Barabási, A.-L. (2002). Hierarchical organization of modularity in metabolic networks. *Science*, 297(5586), 1551–1555. https://doi.org/10.1126/science.1073374 (arXiv:cond-mat/0209244). **VERIFIED** (preprint read).
28. Holme, P., Huss, M., & Jeong, H. (2003). Subnetwork hierarchies of biochemical pathways. *Bioinformatics*, 19(4), 532–538. https://doi.org/10.1093/bioinformatics/btg033. **VERIFIED** (abstract).
29. Solé, R. V., & Munteanu, A. (2004). The large-scale organization of chemical reaction networks in astrophysics. *Europhysics Letters*, 68(2), 170–176. https://doi.org/10.1209/epl/i2004-10241-3 (arXiv:cond-mat/0406137). **VERIFIED** (preprint read).
30. Wong, M. L., Prabhu, A., Williams, J., Morrison, S. M., & Hazen, R. M. (2023). Toward network-based planetary biosignatures: Atmospheric chemistry as unipartite, unweighted, undirected networks. *JGR Planets*, 128(6), e2022JE007658. https://doi.org/10.1029/2022JE007658. **VERIFIED** (abstract).
31. Kim, H., Smith, H. B., Mathis, C., Raymond, J., & Walker, S. I. (2019). Universal scaling across biochemical networks on Earth. *Science Advances*, 5(1), eaau0149. https://doi.org/10.1126/sciadv.aau0149. **VERIFIED** (full text read, PMC6357746).
32. Smith, H. B., Kim, H., & Walker, S. I. (2021). Scarcity of scale-free topology is universal across biochemical networks. *Scientific Reports*, 11, 6542. https://doi.org/10.1038/s41598-021-85903-1. **VERIFIED** (abstract).
33. Broido, A. D., & Clauset, A. (2019). Scale-free networks are rare. *Nature Communications*, 10, 1017. https://doi.org/10.1038/s41467-019-08746-5. **VERIFIED** (abstract).
34. Gagler, D. C., Karas, B., Kempes, C. P., Malloy, J., Mierzejewski, V., Goldman, A. D., Kim, H., & Walker, S. I. (2022). Scaling laws in enzyme function reveal a new kind of biochemical universality. *PNAS*, 119(9), e2106655119. https://doi.org/10.1073/pnas.2106655119. **VERIFIED** (full text read, PMC8892295).
35. Goldford, J. E., Hartman, H., Smith, T. F., & Segrè, D. (2017). Remnants of an ancient metabolism without phosphate. *Cell*, 168(6), 1126–1134.e9. https://doi.org/10.1016/j.cell.2017.02.001. **VERIFIED** (abstract; the 315/260 figures come from Boston University press coverage).
36. Xavier, J. C., Hordijk, W., Kauffman, S., Steel, M., & Martin, W. F. (2020). Autocatalytic chemical networks at the origin of metabolism. *Proceedings of the Royal Society B*, 287(1922), 20192377. https://doi.org/10.1098/rspb.2019.2377. **VERIFIED** (full text read, PMC7126077).
37. Wołos, A., Roszak, R., Żądło-Dobrowolska, A., Beker, W., Mikulak-Klucznik, B., Spólnik, G., Dygas, M., Szymkuć, S., & Grzybowski, B. A. (2020). Synthetic connectivity, emergence, and self-regeneration in the network of prebiotic chemistry. *Science*, 369(6511), eaaw1955. https://doi.org/10.1126/science.aaw1955. **VERIFIED** (abstract). The biotic-against-abiotic property detail is **UNVERIFIED** (seen only in a search-engine rendering of the summary).
38. Robinson, W. E., Daines, E., van Duppen, P., de Jong, T., & Huck, W. T. S. (2022). Environmental conditions drive self-organization of reaction pathways in a prebiotic reaction network. *Nature Chemistry*, 14, 623–631. https://doi.org/10.1038/s41557-022-00956-7. **VERIFIED** (abstract).
39. Saylam, A., Hadj Ali, K., & Fikri, M. (2020). Degree centrality of combustion reaction networks for analysing and modelling combustion processes. *Combustion Theory and Modelling*, 24(3), 442–459. https://doi.org/10.1080/13647830.2019.1699167. **VERIFIED** (abstract).
40. Marshall, S. M., Mathis, C., Carrick, E., Keenan, G., Cooper, G. J. T., Graham, H., Craven, M., Gromski, P. S., Moore, D. G., Walker, S. I., & Cronin, L. (2021). Identifying molecules as biosignatures with assembly theory and mass spectrometry. *Nature Communications*, 12, 3033. https://doi.org/10.1038/s41467-021-23258-x. **VERIFIED** (full text read, PMC8144626).
41. Hazen, R. M., Burns, P. C., Cleaves, H. J., Downs, R. T., Krivovichev, S. V., & Wong, M. L. (2024). Molecular assembly indices of mineral heteropolyanions: Some abiotic molecules are as complex as large biomolecules. *Journal of the Royal Society Interface*, 21(211), 20230632. https://doi.org/10.1098/rsif.2023.0632. **VERIFIED** (abstract).

**Morphospace**

42. Raup, D. M. (1966). Geometric analysis of shell coiling: General problems. *Journal of Paleontology*, 40(5), 1178–1190. No DOI (JSTOR). **VERIFIED** (bibliographic record only).
43. McGhee, G. R., Jr. (1999). *Theoretical Morphology: The Concept and Its Applications*. Columbia University Press. ISBN 0-231-10616-5. **VERIFIED** (bibliographic record only).
44. Avena-Koenigsberger, A., Goñi, J., Solé, R., & Sporns, O. (2015). Network morphospace. *Journal of the Royal Society Interface*, 12(103), 20140881. https://doi.org/10.1098/rsif.2014.0881. **VERIFIED** (full text read, PMC4305402).
45. Corominas-Murtra, B., Goñi, J., Solé, R. V., & Rodríguez-Caso, C. (2013). On the origins of hierarchy in complex networks. *PNAS*, 110(33), 13316–13321. https://doi.org/10.1073/pnas.1300832110 (arXiv:1303.2503). **VERIFIED** (preprint read).
46. Shoval, O., Sheftel, H., Shinar, G., Hart, Y., Ramote, O., Mayo, A., Dekel, E., Kavanagh, K., & Alon, U. (2012). Evolutionary trade-offs, Pareto optimality, and the geometry of phenotype space. *Science*, 336(6085), 1157–1160. https://doi.org/10.1126/science.1217405. **VERIFIED** (abstract).
47. Dittrich, P., Ziegler, J., & Banzhaf, W. (2001). Artificial chemistries—a review. *Artificial Life*, 7(3), 225–275. https://doi.org/10.1162/106454601753238636. **VERIFIED** (full text read).
48. Gerber, S. (2017). The geometry of morphospaces: Lessons from the classic Raup shell coiling model. *Biological Reviews*, 92(2), 1142–1155. https://doi.org/10.1111/brv.12276. **VERIFIED** (abstract).
49. van Wijk, B. C. M., Stam, C. J., & Daffertshofer, A. (2010). Comparing brain networks of different size and connectivity density using graph theory. *PLoS ONE*, 5(10), e13701. https://doi.org/10.1371/journal.pone.0013701. **VERIFIED** (abstract).

*Why there are 49 references.* The brief named about 33 sources. I added one group that fixes Chemart's current growth measures (McShea 1994, Standish 2003 and Droop & Hickinbotham 2012, the last because QNN runs on Stringmol). I added a second group that tests the published anchors: Arita 2004, Broido & Clauset 2019, Smith et al. 2021, Wong et al. 2023 and Hazen et al. 2024. The rest concern the level question: McShea 2001, Rasmussen et al. 2001, Krakauer et al. 2020 and Stepney & Hickinbotham 2024.
