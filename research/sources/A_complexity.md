# A. What "complexity" means, and how it is measured

*Literature report for Chemart, subtopic A (written 2026-09-23). Research only: no project files were changed.*

**How to read this.** Section 1 is the map; read it first. Section 2 has one entry per concept or measure. Section 3 lists where the literature disagrees and what is still open. Section 4 lists the references, each marked VERIFIED or UNVERIFIED. I checked every reference in this session: metadata through Crossref, PubMed or the arXiv API, and definitions from the paper itself, its abstract, or a copy hosted by the author or the publisher. Where I could see only an abstract or an index summary, the entry says so. Statements tagged **[inference]** are mine, not a source's.

**Levels used throughout.** The brief lists four kinds of data. They map onto Marco's three levels like this:

| code | what is measured | Marco's level (Fontana & Buss) |
|---|---|---|
| **M** | one molecule's internal structure (string, tree, graph) | level 0: molecules |
| **N** | a reaction network (species, reactions, catalysts, rates, food set) | level 1: networks (structure) |
| **T** | a population: a snapshot of the soup, or abundance trajectories over time | level 1: networks (dynamics) |
| **S** | a system of interacting networks (several organisations and what links them) | level 2: networks of networks |

---

## 1. Landscape

### 1.1 One word, five families

Lloyd (2001) collected about forty measures of complexity and sorted them by three questions: "1) How hard is it to describe? 2) How hard is it to create? 3) What is its degree of organization?" He argued that "the many measures of complexity represent variations on a few underlying themes." Building on Lloyd, the sources below fall into five families. The grouping is mine.

| family | the question it answers | examples | what scores highest |
|---|---|---|---|
| **1. Randomness** (description length) | How much information does it take to specify this exact object? | Shannon entropy, Kolmogorov complexity, LZ/zlib/brotli compression, BDM | pure noise |
| **2. Structure** | How much pattern (memory, regularity) is there beyond noise? | excess entropy (= effective measure complexity), statistical complexity, effective complexity | mixtures of order and randomness; zero at both ends |
| **3. Scale and hierarchy** | At how many scales, or how many nested levels, is there structure? | TSE neural complexity, complexity profile, Huberman–Hogg tree diversity, McShea's hierarchy levels, Simon's near-decomposability | systems that are both segregated and integrated; many nested levels |
| **4. History** (construction) | How much work, or how many steps, did it take to make? | logical depth, thermodynamic depth, assembly index and assembly | objects that needed a long, non-trivial production history |
| **5. Function** | How much information is *about* something (an environment, a task)? | physical complexity, functional information | rare configurations that work |

The families overlap. Logical depth is both history and structure, and TSE complexity is both structure and scale. Lloyd's three questions line up roughly with them. His "description" list (entropy, algorithmic complexity, Lempel–Ziv complexity) is family 1. His "creation" list (logical depth, thermodynamic depth) is family 4. His "organization" list (effective complexity, excess entropy, effective measure complexity, ε-machine size, hierarchical complexity, tree subgraph diversity, mutual information) is families 2 and 3.

### 1.2 Why the split matters for artificial chemistries

1. **A random soup maximises family 1.** Gell-Mann & Lloyd put it bluntly: algorithmic information content "is not properly a measure of complexity, since randomness is not what we usually mean when we speak of complexity. Another name for AIC, 'algorithmic randomness,' is somewhat more apt." AlChemy and BFF start from random soups, so family-1 measures *fall* when organisation appears. Even combinatory chemistry, which starts from only S, K and I, shows the same fall once structures take over: diversity "explodes at the beginning … peaking very early on, then decreasing" (Kruszewski & Mikolov). Chemart's `compressibility` (zlib on the canonical reaction list) and `shannon` belong to family 1. They will rank a random network or soup as the most "complex".
2. **Replication makes copies, and the families disagree about copies.** Copies push Shannon diversity towards zero. They push BFF's "high-order entropy" up (its stated Property 2). Assembly theory's copy-number term requires them. Depth should ignore them: Lloyd & Pagels require that "a complex object together with a copy is not much deeper than the object alone." So a level-0 takeover by one replicator looks trivial, complex or unchanged, depending on which family you picked.
3. **Vanishing at both ends is not enough.** Feldman & Crutchfield (1998) showed that a measure built only to vanish at the ordered and disordered extremes can end up "a trivial function of the entropy density". A useful measure "must also be defined in a setting that gives a clear interpretation to what structures are quantified."
4. **No measure is general.** Mitchell (2009, pp. 109–110) closes her survey this way: each measure "captures something about our notion of complexity but all have both theoretical and practical limitations, and have so far rarely been useful for characterizing any real-world system." Choose a measure for each claim; there is no single number for "complexity".

### 1.3 Which notion of "complex" the redundancy hypothesis needs [inference]

The hypothesis makes three claims, and each needs a different measurement.

**(a) Level 0: algorithmic molecules are brittle, biological codes are redundant.** This is a claim about the map from a molecule's structure to what it does (robustness, neutrality). It is not a claim about complexity. The literature has a precise word for the genetic code's "redundancy": **degeneracy**, "the ability of elements that are structurally different to perform the same function or yield the same output" (Edelman & Gally 2001). It keeps this apart from **redundancy** in the narrow sense, which "occurs when the same function is performed by identical elements" (Tononi, Sporns & Edelman 1999). Measure level-0 redundancy directly: the share of neutral point mutations, and degeneracy. Expect family-5 measures to *drop* as redundancy rises. Physical complexity gives neutral sites zero information, and functional information falls when more configurations work. That drop is predicted, so it is not evidence against the hypothesis.

**(b) There is a ceiling on the complexity of organisation.** The measure must:

- be near zero for a random soup;
- stay low for a takeover by one replicator (copies must not count);
- grow with the number of differentiated, interdependent components and of nested levels.

Four candidates meet these conditions:

- McShea's hierarchy level with individuation sublevels. It is discrete, and his "trend in the maximum" is exactly a ceiling statistic.
- TSE complexity, or the complexity profile, computed on abundance dynamics. It is continuous.
- The length of the description of an organisation's regularities (effective complexity). Fontana & Buss's grammar-plus-algebraic-laws characterisation is an informal version.
- The functional information of reaching level *k*.

**(c) Level-1 organisations cannot be combined into level-2 ones.** This needs a measure that tells integration apart from mere coexistence:

- The complexity profile of two independent subsystems is the sum of their profiles (Allen, Stacey & Bar-Yam 2017), so any excess over the sum is integration.
- TSE complexity is zero across independent parts.
- Mathis et al. (2024) used set similarity to classify what happens when two organisations are combined.
- Fontana & Buss looked for "glue".

Two rules follow.

- **Measure cause and effect independently.** McShea (1996): "Are complex organisms more evolvable? To answer such questions, we must be able to measure the two variables involved … independently." Redundancy (the independent variable) and organisational complexity (the dependent variable) need separate operations. Otherwise the test is circular.
- **One formal result points the hypothesis' way.** Allen, Stacey & Bar-Yam: "A large-scale behavior requires redundant information among the many components engaged in that behavior." Because the total scale-weighted information is conserved, redundancy shifts information from small scales to large ones. Their "redundancy" is information shared among components, not neutral mutations. Showing that the two coincide is the hypothesis' job.

### 1.4 Summary table

M, N, T, S as defined above. "At the extremes" gives the behaviour stated in the source.

| measure | family | at the extremes | M | N | T | S | needs | computable? | used on ACs/CRNs |
|---|---|---|---|---|---|---|---|---|---|
| Size, counts (species, length, part types) | size | none; an amoeba has ~225× the base pairs of a human (Mitchell) | ✓ | ✓ | ✓ | ✓ | counts | trivial | yes: max/mean length (K&M, Mathis); Chemart `complexity_drift` |
| Shannon entropy, diversity | 1 | max for a uniform distribution, 0 for one species | symbols | degree dist. | ✓ | – | a distribution | cheap | yes: Mathis "population entropy"; K&M distinct expressions; Chemart `shannon` |
| Kolmogorov complexity via compression | 1 | max for incompressible strings | ✓ | serialised | soup | – | a string, a compressor | upper bounds only; poor on short strings | yes: inside BFF's HOE; Chemart `compressibility` |
| BDM (block decomposition) | 1 | as Kolmogorov | ✓ | adjacency | – | – | strings, arrays | approximation | none found |
| High-order entropy (BFF) | 1-derived: redundancy | →0 for i.i.d. noise; →H(D) for many copies | – | – | soup | – | soup string, compressor | cheap | BFF |
| Excess entropy / EMC | 2 | 0 for i.i.d.; log₂(period) for periodic | long strings | – | ✓ | – | long stationary symbol series | data-hungry | none found |
| Statistical complexity C_μ | 2 | 0 for i.i.d. and for constant | – | – | ✓ | – | long symbol series; ε-machine reconstruction | hard | none found |
| Effective complexity | 2 | small for random and for simple | ✓ | ✓ | – | ✓ | a model class; a judgment of what is regular | not computable | informal analogue in Fontana & Buss |
| Logical depth | 4 (and 2) | random and trivial objects are shallow | ✓ | – | ✓ | – | shortest program and its run time | not computable | none found |
| Thermodynamic depth | 4 | claimed shallow at both ends; critics: tracks entropy rate | – | – | ✓ | – | ensemble of histories; chosen macrostates | arbitrary | none found |
| LMC (entropy × disequilibrium) | "one hump" | 0 at both ends by construction | – | – | ✓ | – | a distribution | cheap | not recommended |
| TSE neural complexity C_N | 3 (and 2) | low if parts fully independent or fully dependent | – | linearised | ✓ | ✓ | joint distribution or covariance of many parts | exponential in parts; sample subsets | none found |
| Degeneracy D_N (Tononi) | 3 (and 5) | low for independent and for identical-redundant parts | – | ✓ | ✓ | ✓ | perturbations and a defined output | as TSE | none found (Chemart has a structural analogue) |
| Complexity profile C(k) | 3 | independent parts: all information at scale 1; rigid system: at large scale | – | – | ✓ | ✓ | entropies of all subsets | exponential; approximations exist | none found |
| Huberman–Hogg tree diversity | 3 | uniform and random trees minimal | parse trees | – | – | hierarchy trees | a rooted tree | cheap | none found |
| McShea hierarchy level | 3 | count of nested levels, plus individuation | – | – | – | ✓ | operational criteria for "individual" | cheap once defined | Fontana & Buss levels are the AC analogue |
| Near-decomposability (Simon) | 3 | – | – | ✓ | – | ✓ | interaction strengths, a partition | cheap | none found |
| Physical complexity (Adami) | 5 | ≈0 without selection; neutral sites count 0 | population of aligned molecules | – | – | – | aligned population near equilibrium | cheap given alignment | yes: Avida |
| Functional information | 5 | high when few configurations work | ✓ | ✓ | – | ✓ | a graded function; sampling of configurations | sampling cost | yes: Avida (Hazen 2007) |
| Bertz / Böttcher indices | chemistry | size-dominated; symmetry lowers them | molecular graphs | – | – | – | a molecular graph | cheap | none found |
| Assembly index a(x) | 4 | high for non-repetitive objects (disputed: ≈ compression) | ✓ | – | – | – | object, building blocks | exact is costly | proposed for AlChemy (Mathis 2024), not done |
| Assembly A (with copy numbers) | 4 | 0 when every object is unique | – | – | ✓ | – | a_i and copy numbers n_i | as above | none found |
| Hyperpath assembly (Flamm et al.) | 4 | – | – | ✓ | – | – | reaction hypergraph, building blocks | integer linear program | a formal framework for CRNs |
| Fontana & Buss organisation description | 2 and 3 | – | – | ✓ | ✓ | ✓ | grammar, laws, center, glue (by hand) | manual | AlChemy |
| K&M pointwise information of consumption | AC-specific | 0 when consumption is at chance level | – | ✓ | ✓ | – | reaction events, a null model | cheap | combinatory chemistry |

"None found" means none found in this search. The search was not exhaustive.

---

## 2. Entries

### 2.1 What "complex" means: taxonomies and definitions

#### Lloyd (2001): three questions

- **What it is.** A two-page tabulation of about forty measures under the three questions quoted in §1.1. "Degree of organization" is split into (a) "Difficulty of describing organizational structure" (effective complexity, excess entropy, effective measure complexity, sophistication, topological ε-machine size, hierarchical complexity, tree subgraph diversity, grammatical complexity, and others) and (b) "Amount of information shared between the parts of a system as the result of its organizational structure" (mutual information, algorithmic mutual information, channel capacity, correlation, stored information, organization). Logical depth, thermodynamic depth and computational complexity sit under "creation". Entropy, algorithmic complexity and Lempel–Ziv complexity sit under "description".
- **For Chemart [inference].** The three questions are three separate axes of a morphospace, not one scale. The ceiling claim is a question-3 claim (organisation), with a question-2 flavour (how hard organisations are to make). Chemart's `compressibility` answers question 1.

#### Mitchell (2009), *Complexity: A Guided Tour*, ch. 7 "Defining and Measuring Complexity"

- **What it is.** A readable survey, with sections on complexity as size, entropy, algorithmic information content, logical depth, thermodynamic depth, computational capacity, statistical complexity, fractal dimension, and degree of hierarchy. It opens with Lloyd's three questions. Size fails on the genome example: the amoeba has "about 225 times as many base pairs as humans do". On statistical complexity: "like effective complexity, statistical complexity is low for both highly ordered and random systems, and is high for systems in between". On hierarchy it presents Simon's near-decomposability and McShea's nestedness scale (§2.4). It warns that "measuring the degree of hierarchy in actual organisms can involve some subjectivity in determining what counts as a 'part' or even a 'level.'"
- **For Chemart.** `complexity_drift` (abundance-weighted mean id length) is "complexity as size". Mitchell's amoeba example is the standard warning against reading size as complexity.

#### Ladyman, Lambert & Wiesner (2013), "What is a complex system?"

- **What it is.** A philosophical review of the features attributed to complex systems and of the measures. Their proposed definition: "A complex system is an ensemble of many elements which are interacting in a disordered way, resulting in robust organisation and memory." They classify measures along two distinctions: computable versus non-computable, and statistical versus deterministic. From the abstract, they argue that "the one that best captures the qualitative notion of the order produced by complex systems is that of the Statistical Complexity." The features they examine include nonlinearity, feedback, spontaneous order, robustness and lack of central control, emergence, hierarchical organisation and numerosity. (The feature list and the definition are quoted from the authors' slide version of the paper; the abstract was checked against the published one.)
- **For the hypothesis [inference].** Two words in their definition map onto measurable things. "Memory" is what excess entropy and statistical complexity quantify. "Robust organisation" is what Fontana & Buss's self-maintenance under perturbation tests. Their definition also expects organisation to arise from *disordered* interactions, which suits artificial chemistries.

#### McShea (1996): four types of complexity

- **What it is.** A "narrow view" in which complexity is "some increasing function of the number of different types of parts or interactions". The view is purely structural: "complexity depends only on number of different parts and interactions and not on their functionality". Two dichotomies give four types:
  - "Object complexity refers to the number of different physical parts in a system, and process complexity to the number of different interactions among them."
  - "Hierarchical object complexity is the number of levels of nestedness of parts within wholes."
  - "Hierarchical process complexity is the number of levels in a causal specification hierarchy."
  - Nonhierarchical complexity "is the number of parts or interactions at a given spatial or temporal scale."

  He adds configurational complexity, the "irregularity of arrangement of parts and interactions, independent of their differentiation", and sets it aside. He rejects equating complexity with entropy: "In the narrow view, complexity is not entropy." He warns that the number of parts depends on "descriptive frame".
- **Needs.** A decision about what counts as a part type and an interaction type at a given scale.
- **Levels [inference].** N: object complexity ≈ number of species types; process complexity ≈ number of distinct reaction types. S: hierarchical object complexity ≈ number of nesting levels of organisations.
- **Limits.** Overall complexity is not defined, because the four types are "conceptually independent". Counts depend on the descriptive frame.
- **For the hypothesis.** It supplies the methodological rule quoted in §1.3: measure complexity independently of the variable you want to correlate it with. It also warns that "trends can also occur 'passively,' even if complexity is not generally advantageous". A higher maximum complexity in redundant chemistries could arise either way.

### 2.2 Measures that are largest for randomness

#### Shannon entropy and diversity

- **What it is.** `H = −Σ_i p_i log p_i` over a probability distribution. The authors' own presentation of Ladyman et al. makes the key point: Shannon entropy "cannot express the notion of randomness, order, or complexity of a single object. It can only express properties of a total set of sequences under some distribution." Crutchfield & Young (1989) state the limitation: "Periodic behavior has low information content; random, high content. What this misses, however, is the statistical simplicity of random behavior."
- **Needs.** A distribution: species abundances (T), symbol frequencies in a string (M), or a degree distribution (N).
- **In ACs.** Mathis et al. (2024) report "population entropy", computed from the species distribution. Kruszewski & Mikolov count distinct expressions. BFF's high-order entropy uses the per-byte entropy as its first term. Chemart has `shannon`, `degree_entropy` and `spectral_entropy`.
- **Limits.** Maximal for a random soup of unique molecules and zero for a single species. That makes it a measure of spread, not of organisation.
- **For the hypothesis [inference].** Use diversity as a covariate, not as an outcome. A degenerate chemistry may hold many *different* molecules doing the same job. That raises diversity without adding organisation.

#### Kolmogorov (algorithmic) complexity and compression

- **What it is.** The algorithmic information content of a string "is the length of the shortest program that will cause a given universal computer U to print out the string and then halt" (Gell-Mann & Lloyd 2003). It is uncomputable. Compressors give upper bounds, and BFF's authors describe "the well-established practice in algorithmic information theory to use Lempel-Ziv-style compressors to approximate Kolmogorov complexity."
- **Needs.** One string and a compressor.
- **Levels.** M (molecule strings), N (a serialised network; Chemart's `compressibility` compresses the canonical reaction list with zlib at level 9), T (the whole soup as one string, as in BFF).
- **Limits.**
  - Random strings score highest, which is why Gell-Mann & Lloyd call it "algorithmic randomness".
  - On short strings the compressor's fixed overhead dominates. In a check I ran (Python zlib, level 9), the 3-byte string `SKI` compresses to 11 bytes (ratio 3.7) and the 10-byte `(S(K(SI)))` to 18 bytes (ratio 1.8). Typical AC molecules are this short, so per-molecule zlib ratios measure length, not structure. **[inference, own check]**
  - For networks, the result depends on how the network is written out: species names, ordering, format. This is McShea's descriptive-frame problem.
- **For the hypothesis [inference].** The Kolmogorov complexity of a molecule says nothing about how redundant its structure-to-function map is. A redundant code can be written in incompressible strings.

#### Block Decomposition Method (Zenil et al. 2018)

- **What it is.** BDM combines the "Coding Theorem Method" for small blocks (algorithmic probability estimated from enumerating small Turing machines) with a sum over a decomposed object. The authors state it "provides efficient estimations of algorithmic complexity but that it performs like Shannon entropy when it loses accuracy", and that it applies to strings, arrays and graphs.
- **Needs.** Strings or adjacency matrices. The authors provide implementations.
- **Levels.** M (short molecule strings are exactly its intended use case), N (adjacency matrices).
- **In ACs.** None found.
- **For Chemart [inference].** It is a better per-molecule estimate of algorithmic randomness than zlib, but it still belongs to family 1.

### 2.3 Measures that vanish for both order and randomness

#### Excess entropy, also called effective measure complexity (Grassberger 1986; Crutchfield & Feldman 2003)

- **What it is.** Take the block entropy `H(L)` of length-L words of a stationary symbol sequence, the entropy rate `h_μ = lim_{L→∞} H(L)/L`, and the finite-L estimate `h_μ(L) = H(L) − H(L−1)`. Then
  - `E = Σ_{L=1}^{∞} [h_μ(L) − h_μ] = lim_{L→∞} [H(L) − h_μ L]`,
  - which equals the mutual information between the past and future halves of the sequence (Crutchfield & Feldman, Props. 7–8).

  "E measures the amount of historical information stored in the present that is communicated to the future." Crutchfield & Feldman note the same quantity appears as "stored information", "effective measure complexity" and "predictive information". Grassberger (1986) introduced it as effective measure complexity. From the abstract as indexed (I could not open the paper), he framed such measures as properties of *ensembles* of patterns rather than individual ones, closer to Shannon theory than to computational complexity. E is 0 for i.i.d. noise and log₂(p) for a period-p process.
- **Needs.** A long, stationary, discrete sequence with enough samples to estimate `H(L)` for growing L. Crutchfield & Feldman warn that "missed regularities are converted to apparent randomness," a problem that "arises particularly for small data sets".
- **Levels.** T (symbolised trajectories, e.g. which organisation or species dominates in each window). M only for long strings, which AC molecules usually are not.
- **In ACs.** None found.
- **Limits.** It needs stationarity. An evolving soup is not stationary, so the measure must be applied to windows or to the long-run regime. The result depends on how states are turned into symbols.
- **For the hypothesis [inference].** E is the "memory" in Ladyman et al.'s definition. It measures temporal structure, not composition. An organisation sitting at a noisy steady state can have low E while being compositionally rich, so E alone cannot carry the ceiling claim.

#### Statistical complexity (Crutchfield & Young 1989)

- **What it is.** Histories of a process are grouped into equivalence classes, the causal states, when they predict the same future. The resulting minimal predictive machine is the "ε-machine". Crutchfield & Young define the α-order graph complexity as the Rényi entropy of the state distribution, `C_α = (1−α)^{−1} log Σ_{v∈V} p_v^α`. Two special cases: `C_0 = log|V|` (topological) and `C_1 = −Σ_{v∈V} p_v log p_v`, "the Shannon entropy" of the causal-state distribution (usually written C_μ). They present it as "a measure of complexity distinct from and dual to the information-theoretic entropies and dimensions." Ladyman et al. endorse it.
- **Needs.** A long discrete time series and a reconstruction algorithm. It is the most data-hungry of the structure measures.
- **Levels.** T.
- **In ACs.** None found.
- **Limits.** Stationarity and a chosen symbolisation are required. AC states are multisets of strings, a huge alphabet, so they must be coarse-grained first.
- **For the hypothesis [inference].** In principle it can quantify how much computation an organisation's dynamics carry. In practice, prefer excess entropy (easier to estimate) unless Chemart's runs are very long.

#### Effective complexity (Gell-Mann & Lloyd 1996; 2003/2004)

- **What it is.** "The effective complexity (EC) of an entity" is "the length of a highly compressed description of its regularities." Formally, the regularities are an ensemble E in which the entity e is a typical member. "We introduce the AIC of the ensemble and call it Y. We then have our technical definition of effective complexity: it is the value of Y for the ensemble that is finally employed." The "total information" is `Σ = Y + I`, where I is the Shannon entropy of the ensemble, and "to within a few bits the smallest possible value of Σ is K ≡ K(e)." The recommended ensemble first minimises Σ and then, under that constraint, minimises Y (1996 abstract as indexed; the 2003 text agrees). Random strings have trivial regularities and simple strings have short ones, so both get small EC.
- **Needs.** A model class (candidate ensembles) and "some judgment" about what counts as regular. The authors say that in most practical cases the regular/random split "depends on some judgment".
- **Levels.** M, N, S (conceptually).
- **In ACs [inference].** No formal use found. BFF's authors say their high-order entropy "shares similarities with sophistication and effective complexity". Fontana & Buss describe an organisation by a grammar plus a finite set of algebraic laws; "it may take many laws to completely specify an organization", one needs 20 equations. That is an informal effective-complexity measure: a description of the organisation's regularities, independent of which molecules happen to be present.
- **Limits.** Not computable; observer-dependent.
- **For the hypothesis [inference].** This is the closest formal notion to "complexity of organisation". A practical proxy would be the size of the organisation's smallest self-maintaining generator set (Fontana & Buss's "center"), or the number of distinct reaction schemata needed to describe its closure.

#### Logical depth (Bennett 1988)

- **What it is.** "Some mathematical and natural objects (a random sequence, a sequence of zeros, a perfect crystal, a gas) are intuitively trivial, while others (e.g. the human body, the digits of π) contain internal evidence of a nontrivial causal history." Depth at significance level s is `D_s(x) = min{T(p) : |p| − |p*| < s and U(p) = x}`: "the least time required to compute it by a s-incompressible program." The measure obeys a "slow-growth law: deep objects cannot be quickly produced from shallow ones by any deterministic process". Bennett's framing of value is directly relevant to the word "redundancy": value lies "not in its information (its absolutely unpredictable parts), nor in its obvious redundancy (verbatim repetitions, unequal digit frequencies), but rather in what might be called its buried redundancy—parts predictable only with difficulty".
- **Needs.** The shortest program and its running time. It is uncomputable. Mitchell: it "does not give a practical way of measuring the complexity of any natural object of interest".
- **Levels.** M; T or S if the soup or the organisation is treated as the object.
- **In ACs.** None found.
- **Limits.** Uncomputable, and dependent on the choice of machine, though Bennett argues this dependence is modest.
- **For the hypothesis [inference].** Bennett's distinction sorts the AC measures. A soup full of verbatim copies has *obvious* redundancy, which is what BFF's high-order entropy detects. Marco's redundancy is a third thing, robustness of the structure-to-function map. In ACs the production history is recorded, so a depth-like proxy is the length of the shortest production route from the food set to a species. Flamm et al.'s hyperpath assembly formalises that (§2.6).

#### Thermodynamic depth (Lloyd & Pagels 1988), and its critique

- **What it is.** "A measure of complexity for the macroscopic states of physical systems is defined. Called depth, the measure is universal … Applied to a Hamiltonian system, the measure is equal to the difference between the system's coarse- and fine-grained entropy, a quantity that we call thermodynamic depth." The authors require "that wholly ordered and wholly random systems are not thermodynamically deep and that a complex object together with a copy is not much deeper than the object alone." Grassberger's review quotes them: the measure "must be proportional to the Shannon entropy of the set of trajectories that experiment determines can lead to the state".
- **Criticism.**
  - Crutchfield & Shalizi (1999): "Depth, therefore, is at root arbitrary," because nobody specified how to choose the macrostates. Moreover "the rate of increase in thermodynamic depth, or *dive*, is the system's reverse-time Shannon entropy rate, and so depth only measures degrees of macroscopic randomness, not structure."
  - Mitchell raises the same coarse-graining objection.
  - Grassberger (1989/2012): by this criterion "a stone (which could be a petrified plant) would be deeper than a plant which cannot be a herbified stone!"
- **Needs.** The ensemble of histories leading to a macrostate. A simulation records histories, but the choice of macrostate stays arbitrary.
- **Levels.** T.
- **For the hypothesis [inference].** Its copy property is exactly what the ceiling claim wants. After the critique, though, use it only as a secondary measure.

#### "One-hump" measures, and why vanishing at both ends is not enough (Feldman & Crutchfield 1998)

- **What it is.** The LMC measure of López-Ruiz, Mancini & Calbet (1995) multiplies entropy by a "disequilibrium", `C = H · D` with `D = Σ_i (p_i − 1/N)²`. It is zero for a perfect crystal (H = 0) and for an ideal gas (D = 0). Feldman & Crutchfield found that it "is neither an intensive nor an extensive thermodynamic variable" and "vanishes exponentially in the thermodynamic limit for all one-dimensional finite-range spin systems". Their fix makes it "a trivial function of the entropy density and hence of no use as a measure of structure or memory." Their conclusion is the test quoted in §1.2(3).
- **For Chemart [inference].** Any composite such as "entropy × (1 − entropy)" falls under this critique. A measure's hump shape alone gives it no credibility.

### 2.4 Scale, integration and hierarchy

#### TSE neural complexity (Tononi, Sporns & Edelman 1994)

- **What it is.** "CN is shown to be high when functional segregation coexists with integration and to be low when the components of a system are either completely independent (segregated) or completely dependent (integrated)." With integration `I(X) = Σ_i H(x_i) − H(X)`, the measure has three equivalent forms (as restated in Tononi et al. 1999):
  - `C_N(X) = Σ_k [⟨H(X_j^k)⟩ − (k/n) H(X)]`
  - `= ½ Σ_k ⟨MI(X_j^k ; X − X_j^k)⟩`
  - `= Σ_k [(k/n) I(X) − ⟨I(X_j^k)⟩]`

  Here ⟨·⟩ averages over all subsets of size k. In the 1999 paper, entropies are estimated "under Gaussian assumptions" from covariance matrices of linear systems driven by noise. The 1994 abstract anticipates wider use: "The approach outlined here may prove useful in analyzing complexity in other biological domains such as gene regulation and embryogenesis."
- **Needs.** A joint distribution or covariance of n components. The number of subsets grows exponentially, so in practice subsets are sampled or the parts are grouped into modules.
- **Levels.** N (linearise the kinetics around a fixed point and inject noise; this is the original procedure); T (covariance of abundance time series); S (let the parts be whole organisations).
- **In ACs.** None found.
- **For the hypothesis [inference].** This is the most direct formalisation of a network of networks. Level-1 organisations should be internally integrated and distinct from each other (segregation), and a level-2 system additionally needs integration *between* them. Two independent organisations contribute zero cross-integration. A glued level-2 system does not.

#### Degeneracy versus redundancy (Tononi, Sporns & Edelman 1999; Edelman & Gally 2001)

- **What it is.** "Degeneracy, the ability of elements that are structurally different to perform the same function, … should be distinguished from redundancy, which occurs when the same function is performed by identical elements." The measures use mutual information between subsets of the system and an output set O, under perturbation (noise injected into the subset):
  - Redundancy: `R(X;O) = Σ_j MI^P(x_j;O) − MI^P(X;O)`.
  - Degeneracy D_N(X;O) is the formal twin of C_N: replace `H(X_j^k)` with `MI^P(X_j^k;O)` in the first form above. It equals the average information shared between bipartitions of X and O.

  Findings: "degeneracy is low both for systems in which each element affects the output independently and for redundant systems in which many elements can affect the output in a similar way but do not have independent effects." Also, "networks that have been selected for degeneracy have high values of complexity." Edelman & Gally (2001) call degeneracy "a well known characteristic of the genetic code and immune systems" and argue "it is a feature of complexity at genetic, cellular, system, and population levels. Furthermore, it is both necessary for, and an inevitable outcome of, natural selection."
- **Needs.** A defined output (in Chemart, e.g. production of target species, or persistence of an organisation), a perturbation scheme (rate noise, knockouts), and the joint statistics.
- **Levels.** N, T, S.
- **In ACs.** None found. Chemart's `degeneracy` (reaction-disjoint routes from the food set, citing Edelman & Gally), `knockout_tolerance` and `synthetic_lethal_pairs` are *structural* analogues, not Tononi's information-theoretic D_N.
- **For the hypothesis [inference].** This is the key correction of vocabulary. "Biological codes are highly redundant" means *degenerate* in this literature. Identical copies lower complexity, while degeneracy is the property tied to higher complexity. When molecules are small neural networks with "tunable redundancy", check which of the two the tuning produces.

#### Complexity profile and marginal utility of information (Allen, Stacey & Bar-Yam 2017)

- **What it is.** Information is assigned to dependencies among components, and each dependency has a scale: the number of components it involves. "The complexity profile of a system A is defined as a real-valued function C_A(y) … whose value equals the total amount of information of scale y or higher in A": `C_A(y) = I({x ∈ D_A : s(x) ≥ y})`. Its properties:
  - Conservation law: `∫_0^∞ C(y) dy = S(D_A) = Σ_a σ(a) H(a)`, "independent of the way the components depend on each other".
  - `C(0) = H(A)`, the joint entropy.
  - It vanishes above the largest scale of organisation.
  - Additivity: "If a system A is the union of two independent subsystems B and C, the complexity profile of the full system is the sum of the profiles for the two subsystems".

  On redundancy and scale: "Information and scale are complementary: As information is about the degree of freedom, scale is about constraints associated with redundancy. A large-scale behavior requires redundant information among the many components engaged in that behavior." The marginal utility of information M(y) is the complementary index, the derivative of how much scale-weighted information y bits can describe. (Quotes from the arXiv version, 1409.4708.)
- **Needs.** Entropies of joint distributions over all subsets of components. It is combinatorial; the authors note that "computationally tractable approximations to the complexity profile have been developed".
- **Levels.** T and S.
- **In ACs.** None found.
- **For the hypothesis [inference].** It offers two things. The additivity property is a crisp test of level-2 integration: if the profile of two combined organisations equals the sum of their separate profiles, they merely coexist. The conservation law is the formal version of the intuition that lower-level redundancy is what lets information exist at higher scales.

#### Complexity of hierarchical structures (Huberman & Hogg 1986)

- **What it is.** Complexity of a rooted tree (a hierarchy), measured by its *diversity*. Grassberger's review summarises their starting point: "neither ordered nor completely random trees should be considered complex. This left the authors with the notion that the complexity of a tree is measured by its lack of self-similarity." Later work (Ceccatto & Huberman, as indexed) coarse-grains a tree into subtrees of fixed depth and sorts them into isomorphism classes. I could not open the 1986 paper, so its exact formula is **not reproduced here**.
- **Needs.** A rooted tree.
- **Levels.** M (λ-terms and combinator expressions are trees; Mathis et al. computed simpler tree statistics: node count, typical depth, branching factor and a "C factor", the ratio of maximum depth to the log of the node count) and S (a hierarchy of organisations drawn as a tree).
- **In ACs.** None found.
- **Criticism.** Grassberger: the quantitative measures "seem somewhat arbitrary (they are not related to the difficulty of any obvious task in a quantitative way)".
- **For the hypothesis [inference].** A level-2 system built from *different* level-1 organisations scores higher than one built from copies of a single organisation. This is the same distinction as McShea's sublevels a and b below.

#### Hierarchical structure scale (McShea 2001)

- **What it is.** "The degree of hierarchical structure of organisms—the number of levels of nesting of lower-level entities within higher-level individuals". The scale has levels, "reflecting number of layers of nestedness", and sublevels, "reflecting degree of individuation at the highest level". For a level-3 entity the sublevels are:
  - 3a: "monomorphic aggregate of level-2 entities";
  - 3b: "differentiated aggregate of level-2 entities";
  - 3c: "differentiated aggregate of level-2 entities with intermediate-level parts".

  The biological levels in Mitchell's summary run from prokaryotic cells, to eukaryotic cells as aggregates of level-1 cells, to multicellular organisms, to colonies. McShea used the scale to document "a trend in the maximum", with the tentative finding that "waiting times for transitions between sublevels may have decreased with increasing hierarchical level".
- **Needs.** Operational criteria that say what counts as an individual at each level.
- **Levels.** S.
- **In ACs.** Fontana & Buss's level 0, 1 and 2 are the AC analogue, though not an application of the scale.
- **Limits.** Mitchell notes some subjectivity in deciding what counts as a part or a level, and that "nestedness only describes the structure of an organism, not any of its functions."
- **For the hypothesis [inference].** This is the most direct way to put a number on Marco's ceiling. Record the maximum hierarchical level, with sublevels, reached by any organisation in a run, and compare chemistries by that maximum. A possible AC operationalisation:
  - level 1 = a self-maintaining organisation (Fontana & Buss; chemical organisation theory);
  - level 2 = an aggregate of level-1 organisations that is itself self-maintaining;
  - sublevel a if the constituent organisations are copies, b if they differ, c if glue or intermediate groupings exist.

#### Near-decomposability and stable intermediate forms (Simon 1962)

- **What it is.** "By a complex system I mean one made up of a large number of parts that interact in a nonsimple way." "By a hierarchic system, or hierarchy, I mean a system that is composed of interrelated subsystems, each of the latter being, in turn, hierarchic in structure until we reach some lowest level of elementary subsystem." Near-decomposability means two things: "(a) in a nearly decomposable system, the short-run behavior of each of the component subsystems is approximately independent of the short-run behavior of the other components; (b) in the long run, the behavior of any one of the components depends in only an aggregate way on the behavior of the other components." The evolutionary argument (the watchmaker parable): "complex systems will evolve from simple systems much more rapidly if there are stable intermediate forms than if there are not. The resulting complex forms in the former case will be hierarchic."
- **Needs.** Interaction strengths and a partition, e.g. within- versus between-organisation reaction fluxes.
- **Levels.** N and S.
- **In ACs.** None found.
- **For the hypothesis [inference].** Redundancy may turn level-1 organisations into Simon's "stable intermediate forms". The parable also needs those subassemblies to be combinable. Mathis et al. found AlChemy's level-1 organisations "robust against collapse" yet rarely combinable. Stability is there; composability is missing. So measure composability directly, for example with pairwise-combination experiments as Mathis did, reporting how often each outcome occurs.

### 2.5 Biological complexity: information about an environment or a function

#### Physical complexity (Adami, Ofria & Collier 2000; Adami 2002)

- **What it is.** Genomic complexity is "the amount of information a sequence stores about its environment" (2000). Physical complexity is "the amount of information that an organism stores, in its genome, about the environment in which it evolves" (2002). In practice, per-site entropies `H(i) = −Σ_j p_j(i) log p_j(i)` are taken to base D, the alphabet size, so each site's entropy lies between 0 and 1. They are estimated "by measuring substitution frequencies at each instruction across the population". Complexity is sequence length minus entropy, `C ≈ ℓ − Σ_i H(i)`. With correlations between sites (epistasis) the full form is `ℓ − H`, where H is the entropy per molecule given the environment. The logic: a locus essential to survival "will be fixed in an adapting population," whereas "inconsequential (neutral) sites will be randomized by the constant mutational load".
- **Claims.** In a fixed environment and a large population, complexity cannot decrease: "genomic complexity is forced to increase" (2000). Adami (2002) allows decreases "in co-evolving systems as well as at high mutation rates, in sexual populations, and in time-dependent landscapes".
- **Needs.** A population of aligned sequences near mutation–selection equilibrium. The Avida experiments fixed genome length; otherwise an alignment is needed.
- **Levels.** M, taken over a population of variants.
- **In ACs.** Avida (Adami et al. 2000).
- **Limits.** Needs equilibrium, alignment and a defined "environment". Correlated sites make the per-site sum an approximation.
- **For the hypothesis.** This is the measure where redundancy shows up most plainly. Neutral sites carry *zero* physical complexity: "The neutral sections that contribute only to the entropy turn out to be exceedingly important for evolution to proceed". A redundant chemistry will therefore score *lower* per molecule, consistent with the hypothesis. **[inference]** In soups without explicit lineages, a "population of a species" must first be defined. BFF's tracer tokens show a way to recover lineages.

#### Functional information (Szostak 2003; Hazen, Griffin, Carothers & Szostak 2007; Wong et al. 2023)

- **What it is.** For a function x with degree of function E_x, `I(E_x) = −log₂[F(E_x)]`, "where F(E_x) is the fraction of all possible configurations of the system that possess a degree of function ≥ E_x" (Hazen et al. 2007). Hazen et al. credit Szostak (2003) with proposing that complexity "can be quantified in the context of specific functions of the system, in contrast to prior formalisms based on genomic, sequence, or algorithmic information". Plots of I against E_x show "steps", which they read as distinct "islands" of solutions. Wong et al. (2023) propose a "law of increasing functional information: The functional information of a system will increase (i.e., the system will evolve) if many different configurations of the system undergo selection for one or more functions."
- **Needs.** A graded function, a way to sample the configuration space, and an assay for the function.
- **Levels.** M (e.g. "is a replicator", "catalyses X"), N (network functions such as self-maintenance or producing X), S (organisation-level functions).
- **In ACs.** Avida (Hazen et al. 2007), where they note that relevant commands "can be spread apart by neutral commands".
- **Limits.** Function-specific: there is no single number per system. It needs a sampling estimate of tiny fractions.
- **For the hypothesis [inference].** A redundant structure-to-function map makes F(E_x) larger and I(E_x) smaller: functions become cheaper to find. The hypothesis then predicts that redundant chemistries reach larger E_x, such as level-2 self-maintenance, because each step up costs less. A direct test is to estimate, for each chemistry, the fraction of random initial soups that reach level k. Mathis et al.'s result (level-1 organisations common, level-2 rare) reads as low information for level 1 and high information for level 2 in AlChemy.

### 2.6 Chemistry

#### Molecular complexity indices (Bertz 1981; Böttcher 2016)

- **What it is.**
  - Bertz introduced the first general graph-theoretic index. As implemented in RDKit's `BertzCT` ("Consists of a sum of two terms, one representing the complexity of the bonding, the other representing the complexity of the distribution of heteroatoms"), the bonding term is `2η log₂η − Σ_i η_i log₂η_i` over the η "connections" (pairs of adjacent bonds) grouped into symmetry classes. The heteroatom term is an entropy of the atom-type distribution. I verified the formula through the RDKit source installed in Chemart's venv, not from the original paper.
  - Böttcher's index C_m is "derived from abstracting the information content of a molecule by the degrees of freedom in the microenvironments on a per-atom basis, allowing the molecular complexity to be calculated in a simple and additive way." It "is sensitive to stereochemistry, heteroatoms, and symmetry", and "its additive character … supports direct comparisons of chemical reactions."
- **Needs.** A molecular graph.
- **Levels.** M.
- **In ACs.** None found. AC molecules are strings, trees or matrices; the analogue would be computed on the parse tree or program graph. **[inference]**
- **Limits.** Chemistry-specific and dominated by size. They say nothing about function.
- **For Chemart [inference].** Their relevance is low. Böttcher's additivity suggests one idea, though: a per-reaction change in complexity, to ask whether a network builds complexity reaction by reaction.

#### Assembly theory (Marshall et al. 2021; Sharma et al. 2023; Flamm, Merkle & Stadler 2025)

- **What it is.**
  - The molecular assembly index (MA) "is the length of the shortest of those pathways, i.e. the smallest number of joining operations requires to construct the object, where objects created in the process can subsequently be reused". The heuristic "split-branch" algorithm gives an upper bound, `MA = Σ_i (MA_i + N_i − 1)` over duplicated substructures (Marshall et al. 2021). By their probability model, "molecules with a MA of between 15 and 20 would have a chance formation of one molecule in 10^23". That is the basis of the MA ≈ 15 biosignature threshold.
  - Sharma et al. (2023) define assembly for an ensemble: `A = Σ_{i=1}^{N} e^{a_i} (n_i − 1)/N_T`. Here a_i is the assembly index of object i, n_i its copy number, N the number of unique objects and N_T the total number of objects. "Finding more than one identical copy indicates the presence of a non-random process generating the object." Objects can be strings, and a peptide chain is one of their examples.
  - Flamm, Merkle & Stadler (2025) show that "assembly pathways coincide with certain minimal hyperpaths in B-hypergraphs. This makes it possible to generalize the notion of assembly to general chemical reaction systems". The quantities are computed by integer linear programming.
- **Needs.** An object and its building blocks (a_i); copy numbers (A); a reaction hypergraph and building blocks (hyperpath version).
- **Levels.** M (a_i), T (A), N (hyperpath assembly relative to a food set; Chemart has both the network and the food set). **[inference for N]**
- **In ACs.** Mathis et al. (2024) propose "tracking the Assembly Index and the Assembly Space of λ expressions through time" in AlChemy. No application found.
- **For the hypothesis [inference].**
  1. For a single object, a_i is high for non-repetitive objects, because reuse lowers it. That places it near family 1, which is the centre of the dispute below.
  2. Assembly A counts *identical* copies. In a degenerate chemistry, functionally equivalent but non-identical variants are not pooled, so A may *undercount* selection there. Sharma et al. mention "multiple realizability", but the equation does not use it.
  3. The network (hyperpath) version measures how deep a species sits in the network's production history. That is a history measure at level 1, and a natural depth proxy for Chemart.

#### The assembly theory dispute

- **Critiques.**
  - Uthamacumaran, Abrahão, Kiani & Zenil (2024): the assembly pathway method "is an encoding scheme widely used by popular statistical compression algorithms". AT "performs similarly to other simple coding schemes and underperforms compared to system-related indexes based upon algorithmic probability". The living/non-living separation "has been reported before".
  - Abrahão, Hernández-Orozco, Kiani, Tegnér & Zenil (2024) claim to "formally prove the equivalence between Assembly Theory (AT) and Shannon Entropy via a method based upon the principles of statistical compression" of the LZ family. They add that "the assembly index is equivalent to the size of a minimal context-free grammar", and that AT does not explain selection "that could not have been arrived at using Shannon Entropy".
  - Hazen, Burns, Cleaves, Downs, Krivovichev & Wong (2024): abiotic mineral heteropolyanions reach assembly indices up to 21 ("ewingite and ilmajokite"), so "values of molecular assembly indices ≥15 do not represent unambiguous biosignatures."
  - Jaeger (2024): AT "has merit but is not nearly as novel or revolutionary as claimed. It certainly does not provide any new explanation of biological evolution or natural selection", and the paper's presentation "is starkly distorted by hype".
- **Rebuttal.** Kempes, Lachmann, Iannaccone, Fricke, Chowdhury, Walker & Cronin (2025): "Unlike computational complexity theory, which often emphasizes minimal description length via compressibility, AT explicitly focuses on the causation captured by selection". They give "mathematical examples demonstrating that the assembly index is fundamentally distinct from complexity metrics like Shannon entropy, Huffman encoding, and Lempel-Ziv-Welch compression", and "proofs showing that the assembly index belongs to a different computational complexity class compared to these measures".
- **Fair summary [assessment].** The two sides partly talk past each other. The critics target the index read as a compression of a *single object*, and the biosignature *threshold*. The proponents stress that *index and copy number together* are the evidence of selection, and that the index is a *measurable* quantity (by mass spectrometry). Hazen et al.'s minerals are the strongest empirical challenge to the threshold. For Chemart, AT is useful as a history measure on strings and, through Flamm et al., on networks. It is not a settled life detector.

### 2.7 Complexity as measured in Turing gases and other artificial chemistries

#### Fontana & Buss (1994): what they measured

- **Summary.** No scalar complexity measure. Organisations are "recognized and defined by these syntactical and functional regularities. Objects retained within an organization are characterized by a grammar and interactions among objects by algebraic relationships." An organisation is self-maintaining and has "(i) boundaries established by the invariances, (ii) strong self-repair capabilities responsible for a robustness to perturbation, and (iii) a center, defined as the smallest set kinetically persistent and self-maintaining generator set of the algebra."
- **Levels.** "Level 0 is defined by self-copying objects or simple ensembles of copying objects. Level 1 denotes a new object class, whose objects are self-maintaining organizations made of Level 0 objects, and Level 2 is defined by self-maintaining metaorganizations composed of Level 1 organizations."
- **What they counted.** Diversity, e.g. "the diversity of the system halves from about 620 different object species to slightly above 300. This reflects the fact that A is a simpler organization than B in terms of its sustained portion". The number of laws needed ("20 equations").
- **Level 2.** It shows as "glue": "The center of an L2-organization is exactly the sum of the constituent L1-centers". It is hard to obtain: "It greatly facilitates the construction of L2-organizations, if the glue, or even the constituent L1-organizations operate under different boundary conditions, in particular different collision rules".
- **Reading [inference].** Their implicit complexity is (i) the number of object species sustained and (ii) the length of the grammar-plus-laws description, an effective-complexity notion. They also admit that obtaining level 2 required special boundary conditions, which bears on the ceiling.

#### Mathis, Patel, Weimer & Forrest (2024), "Return to AlChemy"

- **Finding.** "complex, stable organizations emerge more frequently than previously expected, … these organizations are robust against collapse into trivial fixed-points, but … these stable organizations cannot be easily combined into higher order entities."
- **Measures.**
  - Number of unique expressions (level 0 shows up as a collapse towards one).
  - "Population entropy" of the species distribution (with 1000 unique expressions it equals 3, so the log is base 10).
  - Average expression length.
  - The Jaccard index of expression sets, to track similarity through time and after perturbation. Combining two level-1 organisations was classed as "Dominance", "Coexistence" (average similarity > 0.1 to both inputs) or "Mutual Destruction" (< 0.1 to both) across 455 pairs. Coexistence "rarely occurs for organizations evolved in different simulations."
  - Tree statistics of the random expression generators (nodes, typical depth, branching factor, C factor).
- **Relevance [inference].** Their combination experiment is the empirical form of the ceiling, and it is behavioural: survival of set similarity. No information-theoretic measure of integration between organisations was computed, which leaves room for TSE or the complexity profile.

#### Kruszewski & Mikolov (2021/2022), combinatory chemistry

- **Measures.** Global: "the length of the largest expression in the system", "the number of distinct expressions" (diversity), and the proportion of reducing reactions. The paper also proposes "a novel measure of emergent complexity". Because structures sustain themselves by consuming reactants, they "propose tracking reactants consumption as a proxy metric". O(x) is the share of reactants consumed as x (third argument of S-reductions) in a time window. To correct for frequency, they use "the (positive) pointwise information", `I(x) = max(log(O(x)/R(x)), 0)`. R(x) is the equilibrium relative frequency of x under a null model with only random condensation and cleavage, `x*_x ∝ e^{−b|x|}`, with b fixed by mass conservation.
- **Findings.** Diversity "explodes at the beginning … peaking very early on, then decreasing"; expression length increases.
- **Relevance [inference].** It detects autocatalytic structures at level N/T against an explicit null model, a good design pattern for any Chemart measure. It is not a hierarchy measure.

#### Agüera y Arcas et al. (2024), "Computational Life" (BFF): high-order entropy

- **Exact definition.** "we define the high-order entropy of a length n string as the difference between (1) its Shannon entropy (computed over individual tokens – i.e. bytes) and (2) its 'normalized' Kolmogorov complexity (i.e. its Kolmogorov complexity divided by n)." That is, `HOE(s) = H₁(s) − K(s)/n`. Intuitively it is "meant to capture the amount of information that can only be explained by relations between different characters."
- **Stated properties.**
  - "Given a sequence of n i.i.d. characters, its expected high-order entropy converges to 0 as n grows to infinity."
  - "Given a sequence of k i.i.d. characters with distribution D, the expected high-order entropy of the string formed by concatenating n copies of those characters converges to the Shannon entropy of D."
  - Hence "random noise will have no measurable complexity … while a soup obtained from many copies of the same string (as might arise from one 'taken over' by a self-replicator) will have substantial non-zero complexity."
- **Implementation.** K is estimated by "the compressed size of the string achieved by a state-of-the-art text compressor", specifically "brotli -q2", on the whole soup (in the case study, 2^17 tapes of 64 bytes). The units are not stated in the text I read; bits per byte is the consistent choice.
- **Behaviour.** High-order entropy "increases in the first 1000 epochs", and its jumps align with drops in "unique tokens" (tracer tokens). Self-replicators appear "in 40% of the runs within 16k epochs."
- **Reading [inference].** For a long stationary string, K/n tends to the entropy rate h_μ, so HOE approaches `H(1) − h_μ`. In Crutchfield & Feldman's notation this is the per-symbol redundancy `r(1) = h_μ(1) − h_μ`, the first term of the sum that defines excess entropy. HOE is therefore a *redundancy* measure of Bennett's "obvious" kind (verbatim repetition). It detects replication, level 0, and cannot by itself measure the complexity of organisation.

#### Avida (digital organisms)

Physical complexity (Adami et al. 2000) and functional information (Hazen et al. 2007) were both computed in Avida (§2.5). These are the only complexity measures in this report that have been applied to an artificial-life system with a function-based definition.

---

## 3. Disagreements and open questions

1. **There is no consensus measure, and even the grouping of measures is contested.** Lloyd calls them variations on a few themes. Mitchell says they have "rarely been useful" on real systems and "probably can't be captured by a single measurement scale". Ladyman et al. pick statistical complexity, while McShea deliberately uses a function-free count of part types.
2. **Randomness or structure?** Gell-Mann & Lloyd ("algorithmic randomness"), Crutchfield & Young ("the statistical simplicity of random behavior") and McShea ("complexity is not entropy") all reject randomness as complexity. Compression-based measures are still widely used as "complexity": in the AT critiques, in BFF (as one term), and in Chemart's `compressibility`.
3. **Vanishing at both ends is necessary but not sufficient** (Feldman & Crutchfield). Many proposed "one-hump" measures fail the demand to say what structure they quantify.
4. **Individual objects or ensembles?** Kolmogorov complexity, logical depth and the assembly index are properties of one object. Shannon entropy, excess entropy, statistical complexity and TSE complexity are properties of distributions (Grassberger stressed ensembles). McShea treats complexity as a property of one "specific composition and configuration". **[inference]** In Chemart, level-M measures are about individuals, level-T measures about ensembles, and level S needs both.
5. **Function-free or function-laden?** McShea defers function so that complexity can be tested against other variables. Adami, Szostak, Hazen and Wong make function (or information about an environment) the definition. Fontana & Buss define organisations by self-maintenance, itself a function. **[inference]** The redundancy hypothesis is inherently about a structure-to-function map, so function-laden measures are unavoidable for the independent variable, while McShea's argument favours a function-free measure for the dependent variable.
6. **Do copies count?** Yes for BFF's high-order entropy and for assembly A. No for thermodynamic depth, by design. Shannon diversity penalises them. For Bennett, verbatim repetition is "obvious redundancy", not value. The answer decides whether a single-replicator takeover counts as complex.
7. **Is redundancy good or bad for complexity?** It depends on the kind of redundancy. Four meanings are in play:
   - Shannon redundancy, `R = log₂|A| − h_μ` (Crutchfield & Feldman): compressibility.
   - Identical copies: Tononi's "redundancy", which lowers degeneracy and complexity.
   - Degeneracy: different structures, same function; tied to complexity (Tononi 1999; Edelman & Gally).
   - Neutrality: mutations that leave function unchanged. Counted as zero information by physical complexity and functional information, yet "exceedingly important for evolution to proceed" (Adami).

   Bar-Yam's framework adds a fifth, multivariate reading: shared information that creates large-scale structure. Marco's hypothesis concerns degeneracy and neutrality. Its intuition matches the Bar-Yam reading. It would be undermined by measures that reward identical copies. **[synthesis]**
8. **Coexistence or integration at level 2?** No AC study found here measures information shared *between* organisations. Fontana & Buss identified glue by inspection; Mathis et al. used set similarity. TSE complexity and the complexity profile (with its additivity property) are available but untested on ACs.
9. **Hierarchy levels are hard to make operational.** McShea needed "operational and consistent criteria" and still meets subjectivity (Mitchell). Fontana & Buss's grammars and centers were found by hand. **[inference]** Chemart's organisation measures (chemical organisation theory) could supply candidate level-1 units automatically; whether they match Fontana & Buss's organisations is open.
10. **The assembly theory dispute is unresolved** (§2.6). The empirical threshold is challenged (Hazen et al. 2024), and the formal status of the index is contested in both directions (Abrahão et al. 2024 versus Kempes et al. 2025).
11. **Estimation problems specific to ACs.**
    - Molecules are short strings, so zlib overhead dominates.
    - Evolving soups are not stationary, yet excess entropy and statistical complexity assume stationarity.
    - The state alphabet is enormous (multisets of strings).
    - TSE complexity and the complexity profile sum over exponentially many subsets.
    - Finite data turns missed structure into apparent randomness (Crutchfield & Feldman).

    The measures with the best theory (logical depth, effective complexity, statistical complexity) are the hardest to compute, which is Mitchell's point.
12. **Open question.** No study found here validates any complexity measure against organisational level in an artificial chemistry: does the measure rank level 0 < level 1 < level 2 across runs? Doing that validation in Chemart, using Fontana & Buss's levels as ground truth, would itself be a contribution. **[inference]**

---

## 4. References

VERIFIED means metadata confirmed (Crossref, PubMed, arXiv API or publisher page) *and* the content cited here was read in the source: full text, abstract, or an author-hosted copy. Qualifiers in brackets say what was read.

1. Abrahão, F. S., Hernández-Orozco, S., Kiani, N. A., Tegnér, J. & Zenil, H. (2024). Assembly Theory is an approximation to algorithmic complexity based on LZ compression that does not explain selection or evolution. *PLOS Complex Systems* 1(1): e0000014. https://doi.org/10.1371/journal.pcsy.0000014 — **VERIFIED** [abstract]
2. Adami, C. (2002). What is complexity? *BioEssays* 24(12): 1085–1094. https://doi.org/10.1002/bies.10192 — **VERIFIED** [abstract]
3. Adami, C., Ofria, C. & Collier, T. C. (2000). Evolution of biological complexity. *PNAS* 97(9): 4463–4468. https://doi.org/10.1073/pnas.97.9.4463 — **VERIFIED** [full text, PMC18257; equations are images, formula reconstructed from the surrounding text]
4. Agüera y Arcas, B., Alakuijala, J., Evans, J., Laurie, B., Mordvintsev, A., Niklasson, E., Randazzo, E. & Versari, L. (2024). Computational life: How well-formed, self-replicating programs emerge from simple interaction. arXiv:2406.19108 (v2, 2 Aug 2024). https://arxiv.org/abs/2406.19108 — **VERIFIED** [full text]
5. Allen, B., Stacey, B. C. & Bar-Yam, Y. (2017). Multiscale information theory and the marginal utility of information. *Entropy* 19(6): 273. https://doi.org/10.3390/e19060273 ; preprint arXiv:1409.4708 — **VERIFIED** [metadata; quotes from the arXiv full text]
6. Bennett, C. H. (1988). Logical depth and physical complexity. In R. Herken (ed.), *The Universal Turing Machine: A Half-Century Survey*, pp. 227–257. Oxford University Press. Author copy: https://web.cs.ucdavis.edu/~doty/papers/LogicalDepthAndPhysicalComplexity.pdf — **VERIFIED** [full text]
7. Bertz, S. H. (1981). The first general index of molecular complexity. *J. Am. Chem. Soc.* 103(12): 3599–3601. https://doi.org/10.1021/ja00402a071 — **VERIFIED** [metadata; formula via RDKit's `BertzCT` implementation, original not read]
8. Böttcher, T. (2016). An additive definition of molecular complexity. *J. Chem. Inf. Model.* 56(3): 462–470. https://doi.org/10.1021/acs.jcim.5b00723 — **VERIFIED** [abstract]
9. Crutchfield, J. P. & Feldman, D. P. (2003). Regularities unseen, randomness observed: Levels of entropy convergence. *Chaos* 13(1): 25–54. https://doi.org/10.1063/1.1530990 ; arXiv:cond-mat/0102181 — **VERIFIED** [full text]
10. Crutchfield, J. P. & Shalizi, C. R. (1999). Thermodynamic depth of causal states: Objective complexity via minimal representations. *Phys. Rev. E* 59(1): 275–283. https://doi.org/10.1103/PhysRevE.59.275 ; arXiv:cond-mat/9808147 — **VERIFIED** [abstract]
11. Crutchfield, J. P. & Young, K. (1989). Inferring statistical complexity. *Phys. Rev. Lett.* 63(2): 105–108. https://doi.org/10.1103/PhysRevLett.63.105 — **VERIFIED** [full text, scanned copy]
12. Edelman, G. M. & Gally, J. A. (2001). Degeneracy and complexity in biological systems. *PNAS* 98(24): 13763–13768. https://doi.org/10.1073/pnas.231499798 — **VERIFIED** [abstract]
13. Feldman, D. P. & Crutchfield, J. P. (1998). Measures of statistical complexity: Why? *Phys. Lett. A* 238(4–5): 244–252. https://doi.org/10.1016/S0375-9601(97)00855-4 ; arXiv:cond-mat/9708186 — **VERIFIED** [abstract; the LMC formula checked in López-Ruiz et al. 1995, arXiv:nlin/0205033]
14. Flamm, C., Merkle, D. & Stadler, P. F. (2025). Assembly in directed hypergraphs. *Proc. R. Soc. A* 481(2324): 20250331. https://doi.org/10.1098/rspa.2025.0331 ; arXiv:2505.22826 — **VERIFIED** [abstract]
15. Fontana, W. & Buss, L. W. (1994). "The arrival of the fittest": Toward a theory of biological organization. *Bull. Math. Biol.* 56(1): 1–64. https://doi.org/10.1007/BF02458289 ; SFI Working Paper 93-09-055 — **VERIFIED** [full text of the SFI working paper]
16. Gell-Mann, M. & Lloyd, S. (1996). Information measures, effective complexity, and total information. *Complexity* 2(1): 44–52. https://doi.org/10.1002/(SICI)1099-0526(199609/10)2:1<44::AID-CPLX10>3.0.CO;2-X — **VERIFIED** [metadata; abstract seen only as rendered by a search index]
17. Gell-Mann, M. & Lloyd, S. (2003). Effective complexity. SFI Working Paper 2003-12-068; in Gell-Mann, M. & Tsallis, C. (eds.), *Nonextensive Entropy: Interdisciplinary Applications*, Oxford University Press, pp. 387 ff. https://sfi-edu.s3.amazonaws.com/sfi-edu/production/uploads/sfi-com/dev/uploads/filer/a2/0f/a20f7840-5eb8-40a2-9d49-eb5c3456a8b9/03-12-068.pdf — **VERIFIED** [full text of the working paper; book year given as 2003 in Hazen et al. 2007's reference list, book publication usually dated 2004]
18. Grassberger, P. (1986). Toward a quantitative theory of self-generated complexity. *Int. J. Theor. Phys.* 25(9): 907–938. https://doi.org/10.1007/BF00668821 — **VERIFIED** [metadata; abstract seen only as a search-index summary]
19. Grassberger, P. (2012). Randomness, information, and complexity. arXiv:1208.3459; first published in *Proc. 5th Mexican School on Statistical Physics* (EMFE 5), Oaxtepec 1989, ed. F. Ramos-Gómez, World Scientific (1991). https://arxiv.org/abs/1208.3459 — **VERIFIED** [full text]
20. Hazen, R. M., Burns, P. C., Cleaves, H. J., Downs, R. T., Krivovichev, S. V. & Wong, M. L. (2024). Molecular assembly indices of mineral heteropolyanions: some abiotic molecules are as complex as large biomolecules. *J. R. Soc. Interface* 21(211): 20230632. https://doi.org/10.1098/rsif.2023.0632 — **VERIFIED** [abstract]
21. Hazen, R. M., Griffin, P. L., Carothers, J. M. & Szostak, J. W. (2007). Functional information and the emergence of biocomplexity. *PNAS* 104(suppl. 1): 8574–8581. https://doi.org/10.1073/pnas.0701744104 — **VERIFIED** [full text, PMC1876432]
22. Huberman, B. A. & Hogg, T. (1986). Complexity and adaptation. *Physica D* 22(1–3): 376–384. https://doi.org/10.1016/0167-2789(86)90308-1 — **VERIFIED** [metadata only; content described through Grassberger 2012 and the Ceccatto & Huberman 1988 abstract; exact formula not reproduced]
23. Jaeger, J. (2024). Assembly theory: What it does and what it does not do. *J. Mol. Evol.* 92(2): 87–92. https://doi.org/10.1007/s00239-024-10163-2 — **VERIFIED** [abstract]
24. Kempes, C. P., Lachmann, M., Iannaccone, A., Fricke, G. M., Chowdhury, M. R., Walker, S. I. & Cronin, L. (2025). Assembly theory and its relationship with computational complexity. *npj Complexity* 2: 27. https://doi.org/10.1038/s44260-025-00049-9 ; arXiv:2406.12176 (2024) — **VERIFIED** [abstract]
25. Kruszewski, G. & Mikolov, T. (2021, issue published March 2022). Emergence of self-reproducing metabolisms as recursive algorithms in an artificial chemistry. *Artificial Life* 27(3–4): 277–299. https://doi.org/10.1162/artl_a_00355 ; arXiv:2103.08245 — **VERIFIED** [full text, arXiv v3]
26. Ladyman, J., Lambert, J. & Wiesner, K. (2013). What is a complex system? *Eur. J. Philos. Sci.* 3(1): 33–67 (online 2012). https://doi.org/10.1007/s13194-012-0056-8 — **VERIFIED** [metadata and abstract; definition and feature list quoted from the authors' slide version, https://www.bristol.ac.uk/media-library/sites/eng-systems-centre/migrated/documents/ladymanslide.pdf ; the full paper was not accessible]
27. Lloyd, S. (2001). Measures of complexity: a nonexhaustive list. *IEEE Control Systems Magazine* 21(4): 7–8. https://doi.org/10.1109/MCS.2001.939938 — **VERIFIED** [full text]
28. Lloyd, S. & Pagels, H. (1988). Complexity as thermodynamic depth. *Annals of Physics* 188(1): 186–213. https://doi.org/10.1016/0003-4916(88)90094-2 — **VERIFIED** [abstract, OSTI 6857453]
29. Marshall, S. M., Mathis, C., Carrick, E., Keenan, G., Cooper, G. J. T., Graham, H., Craven, M., Gromski, P. S., Moore, D. G., Walker, S. I. & Cronin, L. (2021). Identifying molecules as biosignatures with assembly theory and mass spectrometry. *Nat. Commun.* 12: 3033. https://doi.org/10.1038/s41467-021-23258-x — **VERIFIED** [full text, PMC8144626]
30. Mathis, C., Patel, D., Weimer, W. & Forrest, S. (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos* 34(9): 093142. https://doi.org/10.1063/5.0207358 ; arXiv:2408.12137 — **VERIFIED** [full text, arXiv]
31. McShea, D. W. (1996). Perspective: Metazoan complexity and evolution: Is there a trend? *Evolution* 50(2): 477–492. https://doi.org/10.1111/j.1558-5646.1996.tb03861.x — **VERIFIED** [full text]
32. McShea, D. W. (2001). The hierarchical structure of organisms: a scale and documentation of a trend in the maximum. *Paleobiology* 27(2): 405–423. https://doi.org/10.1666/0094-8373(2001)027<0405:THSOOA>2.0.CO;2 — **VERIFIED** [abstract and introduction]
33. Mitchell, M. (2009). *Complexity: A Guided Tour*. Oxford University Press. ISBN 978-0-19-512441-5. Chapter 7, "Defining and Measuring Complexity". — **VERIFIED** [catalogue record, plus the chapter text read from a copy found online]
34. Sharma, A., Czégel, D., Lachmann, M., Kempes, C. P., Walker, S. I. & Cronin, L. (2023). Assembly theory explains and quantifies selection and evolution. *Nature* 622: 321–328. https://doi.org/10.1038/s41586-023-06600-9 — **VERIFIED** [full text, PMC10567559]
35. Simon, H. A. (1962). The architecture of complexity. *Proceedings of the American Philosophical Society* 106(6): 467–482. https://www.jstor.org/stable/985254 — **VERIFIED** [full text]
36. Szostak, J. W. (2003). Functional information: Molecular messages. *Nature* 423: 689. https://doi.org/10.1038/423689a — **VERIFIED** [metadata only; content through Hazen et al. 2007, which credits it]
37. Tononi, G., Sporns, O. & Edelman, G. M. (1994). A measure for brain complexity: relating functional segregation and integration in the nervous system. *PNAS* 91(11): 5033–5037. https://doi.org/10.1073/pnas.91.11.5033 — **VERIFIED** [abstract; formulas as restated in the 1999 paper]
38. Tononi, G., Sporns, O. & Edelman, G. M. (1999). Measures of degeneracy and redundancy in biological networks. *PNAS* 96(6): 3257–3262. https://doi.org/10.1073/pnas.96.6.3257 — **VERIFIED** [full text, PMC15929]
39. Uthamacumaran, A., Abrahão, F. S., Kiani, N. A. & Zenil, H. (2024). On the salient limitations of the methods of assembly theory and their classification of molecular biosignatures. *npj Syst. Biol. Appl.* 10: 82. https://doi.org/10.1038/s41540-024-00403-y — **VERIFIED** [abstract]
40. Wong, M. L., Cleland, C. E., Arend, D., Bartlett, S., Cleaves, H. J., Demarest, H., Prabhu, A., Lunine, J. I. & Hazen, R. M. (2023). On the roles of function and selection in evolving systems. *PNAS* 120(43): e2310223120. https://doi.org/10.1073/pnas.2310223120 — **VERIFIED** [abstract]
41. Zenil, H., Hernández-Orozco, S., Kiani, N. A., Soler-Toscano, F., Rueda-Toicen, A. & Tegnér, J. (2018). A decomposition method for global evaluation of Shannon entropy and local estimations of algorithmic complexity. *Entropy* 20(8): 605. https://doi.org/10.3390/e20080605 ; arXiv:1609.00110 — **VERIFIED** [abstract]

**Mentioned but not added as separate entries.**

- López-Ruiz, Mancini & Calbet (1995), *Phys. Lett. A* 209: 321–326, doi:10.1016/0375-9601(95)00867-5: the LMC formula, checked in arXiv:nlin/0205033.
- Ceccatto & Huberman (1988), *Physica Scripta* 37: 145, doi:10.1088/0031-8949/37/1/021: abstract read.

Both are VERIFIED.

**Count note.** The brief's named must-cite list was already about 30. The ten additions each carry one specific point:

- Gell-Mann & Lloyd 2003: the source of the verified definitions.
- McShea 2001 and Simon 1962: hierarchy.
- Crutchfield & Shalizi 1999 and Feldman & Crutchfield 1998: the key critiques.
- Tononi et al. 1999 and Edelman & Gally 2001: degeneracy versus redundancy.
- Abrahão et al. 2024: the formal claim that Kempes et al. rebut.
- Flamm et al. 2025: assembly on reaction networks.
- Zenil et al. 2018: estimating complexity for short strings.

No references are UNVERIFIED. Three were verified at metadata level only: Huberman & Hogg 1986, Szostak 2003, and Grassberger 1986 (plus its abstract as an index summary). For those, the report takes content only from other verified sources and says so.
