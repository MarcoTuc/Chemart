# Level 0, the molecule: redundancy, degeneracy, neutrality, robustness, evolvability

*Literature report B for Chemart. It covers definitions, computable measures, brittleness in computational substrates, neural networks as a substrate, and the tension between redundancy and complexity.*

**How to read this.** Every reference except one (Farmer & Belin, quoted by Ray) was checked in this session (search plus the abstract or full text). Section 6 marks each one VERIFIED or UNVERIFIED and says what was read. Verbatim quotes are in quotation marks. Anything labelled **Implication** or **(inference)** is my reasoning, not a claim from the source. Where I looked at the Neural Molecules code (`Artificial_Neural_Chemistry/neural_molecules/diagnostics/{redundancy,species}.py`), it was only to make the implications concrete. For Banzhaf & Yamamoto (2015) I read Marco's local PDF (§8.3, §11.2, ch. 17).

---

## 1. Landscape

**The words are used in overlapping ways, and the differences matter.**

- **Redundancy (many-to-one).** Many genotypes map to one phenotype. In biology it can also mean identical backup parts. It is a property of the *map*: how many genotypes a phenotype has.
- **Degeneracy.** Different parts do the same job in some contexts and different jobs in others. It is a property of *partial functional overlap that depends on context*. Tononi et al. (1999) and Whitacre (2010) both argue that degeneracy, not redundancy, is what goes with complexity and evolvability.
- **Neutrality.** A mutation is neutral if it leaves the phenotype unchanged. Neutrality is always relative to a phenotype definition, and Wagner (2005, FEBS) argues it is also relative to the context.
- **Robustness.** The fraction of mutations that are neutral. It can be measured per genotype or averaged over a phenotype. Redundancy makes robustness possible but does not guarantee it: the redundant genotypes must also be *mutational neighbours* (Hu et al. 2020; Greenbury et al. 2016).
- **Evolvability.** The ability to produce heritable phenotypic variation that is new or adaptive (Kirschner & Gerhart 1998; Payne & Wagner 2019). It is measured by how many distinct phenotypes lie one mutation away.

**The central result of the genotype–phenotype (GP) map literature.** Robustness and evolvability trade off at the level of the *genotype*. At the level of the *phenotype* they are positively correlated (Wagner 2008; Ahnert 2017). The reason is that a robust phenotype has a large, connected neutral network, and that network borders many other phenotypes.

**Five structural properties recur across GP maps** (RNA, lattice proteins, polyominoes, gene networks, linear GP) (Ahnert 2017):

1. Redundancy.
2. Strong bias: a few phenotypes take most of the genotypes.
3. Robustness far above what bias alone would predict: ρ_p ≫ f_p.
4. Many new phenotypes within a few mutations.
5. Positive robustness–evolvability correlation at the phenotype level.

**Brittleness is an old and explicit problem in artificial life.**

- Ray (1991): von Neumann machine code is "brittle, meaning that the ratio of viable programs to possible programs is virtually zero".
- Ofria et al. (2002) measured it: "over 99.7% of all nontrivial mutations are deleterious" for redcode run in Avida.
- The known fixes: remove numeric operands, address by template (Tierra/Avida), bind approximately (Stringmol's Smith–Waterman), and design degeneracy in on purpose (Clark et al. 2011).
- Conrad (1985) stated the governing trade-off: "a system cannot at the same time be effectively programmable, amenable to evolution by variation and selection, and computationally efficient" (quoted in Banzhaf & Yamamoto 2015, p. 349). Conrad (1990) is the earliest statement of Marco's hypothesis I found: "Biological structures that are characterized by a high degree of component redundancy and multiple weak interactions satisfy these conflicting pressures."

**I found no measurement of whether λ-calculus molecules are brittle.** There is no published measurement of the mutational robustness of AlChemy λ-terms, and no paper that explicitly argues AlChemy is limited *because* its molecules are brittle (§3.6).

- Fontana & Buss (1994) present λ-calculus as capturing "chemistry's diversity of equivalence classes, that many different reactants can yield the same stable product". So λ-chemistry is massively *redundant*. Whether it is *robust* is the open question.
- Real software is not as brittle as assumed: "over 30% of random mutations are neutral with respect to their test suite" (Schulte et al. 2014).

**For neural-network molecules, exact redundancy is almost trivial.**

- An irreducible single-hidden-layer tanh network is fixed by its input–output function up to hidden-unit permutations and sign flips (Sussmann 1992). That is h!·2^h = 48 weight settings for h = 3.
- So in Neural Molecules, essentially all neutrality beyond this symmetry comes from how coarsely "behavioural species" are defined, plus reducible networks with dead or duplicate units (inference from Sussmann 1992 and Farrugia-Roberts 2023).
- Weight scale changes the map itself: "The simplicity bias in P(f) becomes weaker as the width σw of the Gaussian Ppar(σw) increases" (Mingard et al. 2025, deep tanh networks).

**The redundancy–complexity tension is real and has a fairly clear resolution (§5).**

- Redundancy *in itself* pushes toward simple phenotypes and away from information:
  - neutral sites carry no information (Adami et al. 2000);
  - redundant maps favour simple outputs (Dingle et al. 2018);
  - frequent phenotypes arrive first (Schaper & Louis 2014);
  - selection for robustness can simplify the phenotype (Milano et al. 2019) or confine the population (Ancel & Fontana 2000).
- *Degeneracy*, and phenotype-level robustness over large *connected* neutral networks that border many phenotypes, is what goes with evolvability and complexity (Tononi et al. 1999; Wagner 2008; Whitacre & Bender 2010; Clark et al. 2011; Greenbury et al. 2022).

### Glossary

| Term | Definition (source) | Measured on | Not to be confused with |
|---|---|---|---|
| Redundancy (biology) | "the same function is performed by identical elements" (Tononi et al. 1999) | A set of parts | Degeneracy |
| Redundancy (GP map) | "many genotypes map to the same phenotypes" (Ahnert 2017); genotypic redundancy of a phenotype = "the total number of genotypes that map to it" (Hu et al. 2020) | The map | Robustness (redundant genotypes may not be neighbours) |
| Degeneracy | "the ability of elements that are structurally different to perform the same function or yield the same output" (Edelman & Gally 2001); functions "overlap partially… Under some conditions the functions are similar while under others they differ" (Whitacre 2010) | A set of parts × contexts | Redundancy (identical parts, identical function in every context) |
| Neutral mutation | A mutation "without phenotypic effect" (Wagner 2005 FEBS). Neutrality is "not an essential feature of a mutation": it can become non-neutral in a new environment or genetic background (same source) | A genotype + mutation operator + phenotype definition | "Nearly neutral" in fitness |
| Neutral set / network / component | All genotypes with phenotype p; the network links genotypes one point mutation apart; its connected pieces are "neutral components" (Schuster et al. 1994; Ahnert 2017) | The map | The population's current spread |
| Robustness (biological) | Systems "continue to function, survive, or reproduce when faced with mutations, environmental change, and internal noise" (Wagner 2005 FEBS) | Any level | Stability of a network's dynamics |
| Genotype robustness r_g | "The number R_G (or fraction r_G) of neutral neighbours of a genotype" (Wagner 2008) | One molecule | Phenotype robustness |
| Phenotype robustness ρ_p | r_g "averaged over all genotypes G with a given phenotype" (Wagner 2008) | The map | Genotype robustness |
| Genotype evolvability e_g | "The number E_G of different structures found in the 1-neighbourhood of a sequence" (Wagner 2008) | One molecule | Adaptive evolvability |
| Phenotype evolvability ε_p | "The number E_P of different structures found in the 1-neighbourhood of a structure P", i.e. of its whole neutral set (Wagner 2008) | The map | Genotype evolvability |
| Evolvability (biology) | "an organism's capacity to generate heritable phenotypic variation" (Kirschner & Gerhart 1998); "phenotypic variation that is both heritable and adaptive" (Payne & Wagner 2019) | Organism / lineage | Mutation rate |
| Evolvability (evol. computation) | "the ability of a population to produce variants fitter than any yet existing" (Altenberg 1994; also quoted in Banzhaf & Yamamoto 2015 §8.3.2) | Population | Robustness |
| Phenotype bias / frequency f_p | f_p = #G_p / #G, the fraction of genotypes with phenotype p; Zipf-like, "relatively few common and many rare ones" (Schuster et al. 1994) | The map | Fitness |
| Simplicity bias | P(x) ≲ 2^(−a·K̃(x) − b) (Dingle et al. 2018) | The map | Occam's razor as a prior choice |
| Brittleness | "the ratio of viable programs to possible programs is virtually zero" (Ray 1991) | Genotype space | Chaotic dynamics |
| Antiredundancy | "a hypersensitivity to mutation" that removes mutants and protects the wild type (Krakauer & Plotkin 2002) | Molecule / cell | Brittleness as a defect |
| Physical complexity | Information a sequence stores about its environment: C = ℓ − Σ_i H_i (Adami et al. 2000) | Population of sequences | Kolmogorov complexity |
| Functional equivalence class (NN) | Parameters implementing the same input–output function; for irreducible tanh networks, the orbit of the permutation and sign-flip group (Sussmann 1992; Farrugia-Roberts 2023) | Weight space | Mode connectivity (equal *loss*, not equal *function*) |

---

## 2. Concepts and measures

Most measures below need the same four ingredients:

- **(G)** a genotype encoding, meaning the molecule's internal structure string;
- **(M)** a mutation operator;
- **(Φ)** a phenotype, or species, function;
- **(S)** optionally, a sampling distribution over genotypes.

For an artificial chemistry, the natural phenotype of a molecule is **its row of the reaction table against a fixed probe set of partners**. That is its products (and optionally rates) when it reacts with each probe. This is what Neural Molecules calls the "fingerprint on a probe battery", and what Clark et al. (2011) call an "interaction function". Using the same definition of Φ for every chemistry is what makes λ-terms, BFF tapes and neural networks comparable. That last point is my inference.

**Practical consequence for Chemart (inference).** Measures that mutate a molecule (§2.4–2.12, 2.15) cannot be computed from an output reaction network alone. They need the chemistry's *reaction function* as an oracle to evaluate the mutants. Two measures work from Chemart's existing outputs:

- Clark et al.'s redundancy and degeneracy of the interaction matrix (§2.3);
- Adami's population-entropy complexity (§2.13), if the molecule strings can be aligned.

### 2.1 Redundancy vs degeneracy (verbal definitions)

**Definitions.**

- Edelman & Gally (2001): degeneracy is "the ability of elements that are structurally different to perform the same function or yield the same output… it is both necessary for, and an inevitable outcome of, natural selection."
- Tononi et al. (1999): "Because structurally different elements may produce different outputs in different contexts, degeneracy should be distinguished from redundancy, which occurs when the same function is performed by identical elements."
- Whitacre (2010): redundancy is "the coexistence of identical components with identical functionality… isomorphic and isofunctional". Degenerate components have functions that "overlap partially… Under some conditions the functions are similar while under others they differ."

**Relevance.** Marco's hypothesis is phrased in terms of *redundancy*. The literature that links this property to complexity and evolvability is almost all about *degeneracy*. The hypothesis may need to be restated as a degeneracy hypothesis: molecules whose variants keep their function with some partners and change it with others.

### 2.2 Information-theoretic redundancy and degeneracy (Tononi, Sporns & Edelman 1999)

**Definition.** A system X has n elements x_i and a set of output units O.

- `MI^P(X_j^k; O)` is the mutual information between a subset X_j^k of size k and O, measured when "that subset is injected with a fixed amount of variance (uncorrelated random noise)". ⟨·⟩ averages over all subsets of size k.
- Redundancy: `R(X;O) = Σ_{i=1..n} MI^P(x_i; O) − MI^P(X; O)`.
- Degeneracy:
  - `D_N(X;O) = Σ_{k=1..n} [ ⟨MI^P(X_j^k; O)⟩ − (k/n)·MI^P(X; O) ]`
  - equivalently `D_N(X;O) = ½ Σ_k ⟨MI^P(X_j^k; X−X_j^k; O)⟩`, where `MI^P(X_j^k; X−X_j^k; O) = MI^P(X_j^k;O) + MI^P(X−X_j^k;O) − MI^P(X;O)`
  - equivalently `D_N(X;O) = Σ_k [ (k/n)·R(X;O) − ⟨R(X_j^k;O)⟩ ]`.
- Replacing MI^P(X_j^k;O) by the entropy H(X_j^k) gives the Tononi–Sporns–Edelman neural complexity C_N(X): "the equations defining the two measures are formally identical".

**Interpretation.** Degeneracy is high when the whole system carries much information about the output, *and* small subsets carry more than a linear share. It is "low both for systems in which each element affects the output independently and for redundant systems in which many elements can affect the output in a similar way but do not have independent effects". The authors also found that "networks that have been selected for degeneracy have high values of complexity."

**Inputs.** A stochastic model of the molecule's internal elements and outputs, and a noise-injection perturbation. The paper used linear Gaussian networks (8 units, 4 outputs), where MI follows from covariance matrices.

**Computed for ACs or programs?** Not that I found. Clark et al. (2011) considered and rejected existing degeneracy measures for AC binding because of "conflation" with redundancy (§2.3).

**Feasibility.**

- NN molecules: feasible. Take the 3 hidden units as the elements and the 3 outputs as O. Inject Gaussian noise into hidden units and estimate MI by linearization or Monte Carlo. This gives an *internal* degeneracy of the molecule: how far hidden units can stand in for each other.
- λ-terms, strings, programs: no natural "elements with variance", so this does not transfer.

**Limits.** It needs a noise model; MI estimation in nonlinear systems is expensive and biased; it depends on the chosen O.

**Relevance.** It is the canonical formal distinction between redundancy and degeneracy, and the original evidence that degeneracy and complexity go together.

### 2.3 Set-theoretic redundancy and degeneracy of an interaction function (Clark et al. 2011, Stringmol)

**Definition.** Take an interaction function `f: A × B → {0,1}`, for example "a binds b", and let `B_a = { b ∈ B : f(a,b) = 1 }`.

- a_m and a_n are **redundant** ⇔ `B_{a_m} = B_{a_n}`.
- a_m and a_n are **degenerate** ⇔ `B_{a_m} ≠ B_{a_n}` and `B_{a_m} ∩ B_{a_n} ≠ ∅`.
- Group A into redundant classes Â, and B into B̂ likewise.
- Redundancy set: `R(A|f,B) = { |â| : â ∈ Â }`.
- Degeneracy set: `D(A|f,B) = { |B̂_â| : â ∈ Â }`, the number of redundant classes of B that each class interacts with.
- Both sets are rescaled by their means so that chemistries of different sizes can be compared. The mean redundancy is |A| / |Â|.

**Result.** In Stringmol the binding system (tailored Smith–Waterman, thresholded) shows a spread of both redundancy and degeneracy. "Sticky-Stringmol", in which every pair binds, has "trivial redundancy and no degeneracy". Removing degeneracy removed "emergent macro-mutations, hypercycles, sweeps and parasite evasion". The authors conclude: "degeneracy in the components of an AChem facilitates the complexity of the system as a whole".

**Inputs.** A finite set of molecules and a binary interaction function. The paper used all length-6 strings over a 5-symbol alphabet.

**Computed for ACs?** Yes. This is the only AC-native measure of degeneracy I found.

**Feasibility.** Computable for *any* Chemart chemistry directly from its reaction network outputs, with f(a,b) = 1 if a and b react. A natural extension, which is my suggestion, lets f(a,b) be the *product species*. Then redundant molecules make the same product with every partner, and degenerate ones agree with some partners and not others.

**Relevance.**

- **Implication:** Marco's redundancy ratio `n_weight_species / n_behaviour_species` has the form of Clark et al.'s *mean redundancy* |A|/|Â|, with weight cells as the elements and behaviour cells as the classes. The match is exact only if every weight cell falls inside a single behaviour cell.
- Clark et al. argue it must be paired with a degeneracy measure, and that degeneracy is the quantity that goes with complexity.
- The paper also remarks: "AlChemy (level 0) had relatively simple binding, which resulted in the collapse of the system into 'self-replicators'", and level-1 restrictions "produced more complex artifacts". This is a claim about AlChemy's *binding rule*, not about how brittle its molecules are.

### 2.4 Neutrality

**Definition.** A mutation is neutral if it does not change the phenotype (Wagner 2005 FEBS). Wagner argues neutrality is not intrinsic: "a once neutral mutation may cause phenotypic effects in a changed environment or genetic background".

**Inputs.** (G), (M), (Φ). In an AC the "environment" of a molecule is the population it reacts with. So neutrality relative to a fixed probe battery (Neural Molecules' choice) and neutrality relative to the current population are different quantities.

**Relevance.** "Neutral" is defined by Φ. Neutral Molecules defines it as "same species" under a quantizer, and its docstring already notes that species richness is resolution-dependent. The literature supports this: a continuous map has neutral networks only once phenotypes are coarse-grained (§4.2).

### 2.5 Genotype (molecule) mutational robustness

**Definitions.**

- Wagner (2008): genotype robustness `r_g = (number of neutral 1-mutants of g) / (number of 1-mutants of g)`. For point mutations over alphabet A and length L, the denominator is L(A−1).
- Ofria, Adami & Collier (2002), for Avida programs: `f = N / ((D − 1)·ℓ)`, where N is the number of neutral *or beneficial* single mutants, D the instruction-set size and ℓ the genome length. Mutations are classified as "fatal, deleterious, neutral, or beneficial". They exclude "trivial" mutations, which "affects only nonexecuted portions of code".
- Schulte et al. (2014), for real software: "the fraction of random mutations that leave a program's behavior unchanged", with behaviour measured by the program's test suite and "mutation operators… taken from genetic programming".

**Inputs.** One molecule, its full or sampled 1-mutant neighbourhood, and Φ. Cost is about L(A−1) phenotype evaluations per molecule.

**Computed for programs or ACs?**

- Avida: the evolved dominant genotypes had robustness far above the ancestor's "very low" 0.005. Robustness "seems to be strongly correlated to sequence length".
- Redcode run in Avida: "over 99.7% of all nontrivial mutations are deleterious" (Ofria et al. 2002).
- Software: ">30% of random mutations are neutral", "for mutations at both the source code and assembly instruction levels" (Schulte et al. 2014).
- Stringmol, AlChemy, BFF: not measured, as far as I found.

**Feasibility.**

| Substrate | Genotype and mutation | Neutral when… | Other outcomes to record |
|---|---|---|---|
| λ-terms | Symbol substitution on a de Bruijn string or a tree operator | Normal form, or behaviour on probes, is unchanged | Ill-formed or non-normalising terms count as "lethal" |
| BFF | 64-byte tape, byte substitution | Output tapes on probe partners are unchanged | — |
| NN weights | Continuous. Marco's `neutral_fraction` is already the right analog: the fraction of Gaussian ε-perturbations that keep the species | — | Report it as a curve over ε |

**Limits.**

- It depends on the mutation operator and on Φ's resolution. Schulte's 30% is relative to a test suite; Ofria's 99.7% is relative to competitive fitness.
- Dead code inflates it (hence Ofria's exclusion of trivial mutations).
- Step sizes for discrete and continuous substrates cannot be compared directly (§2.10 gives a fix).

**Relevance.** This is the direct operationalization of "brittle". It is the first number to measure for AlChemy λ-terms and BFF tapes.

### 2.6 Multiple mutations and epistasis (robustness curves)

**Definition.** Lenski et al. (1999) "introduced millions of single and multiple mutations into each organism and measured the effects on the organism's fitness". They compared the fitness decay with the number of mutations k against a multiplicative expectation, which tests for epistasis.

**Findings.** "The complex organisms are more robust than the simple ones with respect to the average effects of single mutations"; interactions "usually yield higher fitness than predicted".

**Inputs and feasibility.** As in §2.5, iterated to k mutations. For every substrate the curve `P_same(k)` (probability that k random mutations keep the species) is cheap to sample. For NNs the analog is `P_same(ε)`, alongside Marco's `robustness_curve` (displacement vs ε).

**Relevance.** "Change one symbol and evaluation diverges" is a claim about the *first* point of this curve. Epistasis tells you whether robustness holds up as mutations accumulate. Lenski et al.'s result that complexity goes with robustness is one data point *for* the hypothesis.

### 2.7 Genotype evolvability

**Definition.** `e_g = |{ Φ(g′) : g′ ∈ N_1(g), Φ(g′) ≠ Φ(g) }|`, the number of *distinct* new phenotypes one mutation away (Wagner 2008).

**Finding.** "genotype (sequence) robustness and evolvability share an antagonistic relationship" (Wagner 2008, RNA).

**Inputs and feasibility.** Same scan as §2.5, counting distinct species instead of neutral ones. Computable for all substrates.

**Relevance.** It separates "brittle" into two cases: mutations that *destroy* (lethal/undefined) and mutations that *vary* (many distinct new phenotypes). A λ-term whose mutants are all different but still functional would be evolvable-but-not-robust, which is not the same as brittle. That distinction is my inference.

### 2.8 Phenotype robustness and phenotype evolvability

**Definitions.**

- `ρ_p = ⟨r_g⟩` over g ∈ G_p, the neutral set of p (Wagner 2008; Ahnert 2017).
- `ε_p` = the number of distinct phenotypes other than p in the 1-mutant neighbourhood of the whole neutral set. Ahnert (2017): "the total number of different phenotypes that lie within the point-mutation neighbourhood of a phenotype".

**Findings.**

- "phenotype (structure) robustness promotes structure evolvability. A consequence is that finite populations of sequences with a robust phenotype can access large amounts of phenotypic variation while spreading through a neutral network" (Wagner 2008).
- Ahnert's review lists this positive correlation, together with the negative one at genotype level, among the shared properties of GP maps.

**Inputs.** Samples from each neutral set: inverse folding in RNA, or random sampling plus neutral random walks in general. Wagner estimated ε_p from finite samples.

**Computed for programs?** Yes, for linear GP. Hu et al. (2012) measured robustness and evolvability "at the genotypic, phenotypic, and fitness levels" over a fully enumerated compact Boolean LGP. Hu et al. (2020) found network measures of phenotype evolvability positively correlated with genotypic redundancy: "more robust phenotypes are also more evolvable".

**Feasibility.** Neutral random walks work for λ-terms, BFF and quantized NNs. Enumerating whole neutral sets is only possible for tiny genotype spaces.

**Relevance.** This is the measure that can *confirm* the hypothesis. What must be high is *phenotype* robustness together with high ε_p. High robustness of a single molecule on its own is not enough.

### 2.9 Phenotype frequency, bias and simplicity bias

**Definitions.**

- `f_p = |G_p|/|G|`, the fraction of genotype space mapping to p. It is estimated by uniform (or prior) sampling of genotypes (Schuster et al. 1994; Ahnert 2017; Hu et al. 2020, who sampled 10^9 LGP genotypes for 256 Boolean phenotypes).
- Schuster et al. (1994): "Frequencies of structures are highly non-uniform and follow a generalized form of Zipf's law".
- Simplicity-bias bound (Dingle et al. 2018): `P(x) ≲ 2^(−a·K̃(x) − b)`, with K̃ a Lempel–Ziv-based complexity estimate, `a ≈ log2(N_O) / max_x K̃(x)` and b ≈ 0 by default.
  - The bound applies to maps with "limited complexity", "redundancy" (inputs ≫ outputs), "finite size", "nonlinearity", and that are "well-behaved".
  - "many outputs fall well below the upper bound": simple outputs may still be rare.

**Inputs.** A sampling distribution over genotypes (S) and Φ.

- For NNs, f_p is the parameter-space volume, or prior probability, of the behaviour (§4.1).
- For λ-terms, f_p depends on the random-term generator. Mathis et al. (2024) characterized "the initial distribution of objects produced by two random expression generators, and their consequences on the results".

**Relevance.**

- **Implication:** a single redundancy ratio hides this distribution. Reporting the rank–frequency curve of f_p, and whether frequent species are simple (low K̃ of their fingerprint or normal form), shows whether redundancy concentrates on trivial behaviours (§5).
- Neural Molecules calibrates the ratio to 1.0 on *random* molecules. That subtracts out exactly this intrinsic bias, so the ratio measures only redundancy "the chemistry built". This is a sound choice, and conceptually close to the "excess robustness" of van Nimwegen et al. (§2.11) (inference).

### 2.10 Genetic correlations: ρ_p compared with f_p

**Definition.** In a "random GP map" null, where phenotypes are assigned to genotypes at random with the observed frequencies, "the phenotypic robustness therefore is simply ρ_p = f_p". For biological maps, "very roughly, ρ_p ∝ log f_p, so that the robustness is much larger than would be expected for the null model, in fact by several orders of magnitude for smaller f_p" (Greenbury et al. 2016). Non-neutral correlations also exist. Some reduce evolvability, since phenotype diversity next to a genotype is "lower than expected from the random model". Others enhance it.

**Inputs.** f_p (§2.9) and ρ_p (§2.8) for many phenotypes.

**Relevance. This is the key diagnostic for "redundancy vs robustness".**

- A map can be massively redundant (large f_p) and still brittle, if its redundant genotypes are scattered (ρ_p ≈ f_p).
- It also gives the cross-substrate comparison that the incommensurable mutation steps otherwise prevent. Plot ρ_p/f_p for λ-terms, BFF and NNs, each with its own native point mutation.
- To use it on NNs, discretize the weights, e.g. to a few bits each. The continuous chemistry then has a discrete genotype space with point mutations, and every GP-map measure applies unchanged. This is my suggestion, not from the literature.

### 2.11 Neutral networks: connectivity and population-level robustness

**Definitions.**

- Neutral network = the genotypes of p linked by single mutations. Its connected components matter, because a population can only diffuse within one (Schuster et al. 1994; Ahnert 2017).
- Huynen, Stadler & Fontana (1996): on neutral networks, "evolving populations split into subpopulations, which diffuse independently in sequence space". There are two mutation thresholds, "one at which genotypic information is lost and one at which phenotypic information is lost".
- van Nimwegen, Crutchfield & Huynen (1999):
  - The population's limit distribution on a neutral network is the principal eigenvector of its adjacency matrix. The population neutrality equals the spectral radius, `⟨d⟩ = ρ(G)`, when Mμ ≫ 1 (M = population size, μ = mutation rate). It equals the network mean degree d̄ when Mμ ≪ 1.
  - Excess robustness `r ≡ (⟨d⟩ − ⟨d_0⟩)/⟨d_0⟩ ≈ (ρ − d̄)/d̄`.
  - RNA test (L = 18): d̄ = 12.0 and ρ ≈ 15.7, so populations evolve about 30% more robustness than the network average, without selection for robustness.

**Inputs.** Neutral-network samples (random walks), or the population's genotypes plus neighbour scans.

**Relevance.** Two things follow.

- A soup at large Mμ will drift to the robust core of a neutral network on its own. Rising robustness over a run is therefore *expected*, and is not by itself evidence that the substrate's redundancy did the work. It needs a random-walk baseline.
- A measured increase is only meaningful as *excess* over d̄.

### 2.12 Navigability, shape-space covering and phenotype networks

**Definitions.**

- Shape-space covering (Schuster et al. 1994): "All common structures can be accessed from an arbitrary sequence by a number of mutations much smaller than the chain length."
- Navigability (Greenbury, Louis & Ahnert 2022): "the average probability that a randomly chosen phenotype pair have at least one accessible path between them, given a fitness assignment process to phenotypes". An accessible path is made of single-point mutations with monotonically non-decreasing fitness. In the random-fitness test, phenotypes get uniform random fitness and the target gets the maximum.
  - Results (bioRxiv version): "With neutral mutations allowed, navigability is almost always 1.0"; with neutral mutations disallowed it is "markedly reduced… (⟨ψ⟩ ∈ [0.38, 0.64])".
- Phenotype networks (Hu et al. 2020): nodes are phenotypes and edges are weighted by mutational connections. Degree, strength, disparity and eigenvector/PageRank centrality serve as evolvability proxies. "Neutrality is facilitated by redundancy, but not guaranteed… genotypes map to the same phenotype but are not mutationally connected."

**Inputs.** Genotype sampling plus neighbour scans, and a fitness assignment. For ACs, random fitness assignment probes structure and needs no endogenous fitness.

**Relevance.** Navigability is the GP-map quantity closest to "can the ceiling be lifted". It asks whether rare, complex phenotypes can be reached from common ones without crossing valleys. It depends on neutral mutations.

### 2.13 Physical complexity: the information cost of robustness

**Definition** (Adami, Ofria & Collier 2000).

- Per-site entropy across an equilibrated population: `H_i = −Σ_j p_j(i) log_D p_j(i)`, with logs to base D, the alphabet size.
- Complexity: `C = ℓ − Σ_i H_i`, where "The neutral sections that contribute only to the entropy turn out to be exceedingly important for evolution to proceed".
- Mutational-scan approximation: `H_i = log_D N_ν(i)`, where N_ν(i) is "the number of non-lethal substitutions" at site i. This reads as counting the wild-type symbol, so that a fully conserved site gives H_i = 0 (my reading).

**Related statement.** Ofria et al. (2002): "The robustness of a genome can be thought of as the fraction of its length that is impervious to mutations and thus carries no information… learning events… decrease robustness if the sequence length stays constant, while size increases without commensurate acquisition of information increase robustness."

**Inputs.** Either a population of alignable sequences (from a trajectory) or a mutational scan of one molecule (§2.5).

**Feasibility.**

- The population version needs only Chemart's output strings, if they are fixed-length or can be aligned.
- The scan version works for BFF and Avida-like programs directly. For λ-terms it needs a positional encoding.
- For NNs, information is not per-site. The analog is the description length of the behaviour under the weight distribution, −log₂ f_p (§4.1). That is my inference.

**Relevance.** This is the formal core of the tension (§5). At fixed molecule size, per-site robustness is information the molecule does *not* carry. Neural Molecules have a fixed 24 parameters.

### 2.14 Viability density and replicative capacity (brittleness at the scale of the whole space)

**Definitions.**

- Ray (1991): brittleness = "the ratio of viable programs to possible programs is virtually zero".
- Suzuki (2003) optimized a string-rewriting AC for "'replicative capacity', that is the occurrence ratio of self-replicating strings". He found that "a large replicative capacity assures strong connectivity between self-replicating genotypes, making the system highly evolvable".

**Inputs.** (S) and a viability predicate, e.g. "is a self-replicator" or "normalizes within the step cap".

**Feasibility.** Direct Monte Carlo for every chemistry. For AlChemy it is the fraction of random λ-terms that are non-trivial and normalize. For BFF it is the fraction of random tapes that are self-replicators or can grow into them.

**Relevance.** This is the cheapest single number for "how brittle is this chemistry's molecular code". Suzuki's result ties high density to percolating neutral networks, the §2.11 connectivity.

### 2.15 Evolvability: population and system definitions

**Kirschner & Gerhart (1998).** Evolvability is "an organism's capacity to generate heritable phenotypic variation". They name the properties that confer it: "versatile protein elements, weak linkage, compartmentation, redundancy, and exploratory behavior reduce the interdependence of components… confer evolvability on the organism by reducing constraints on change and allowing the accumulation of nonlethal variation."

**Payne & Wagner (2019).** Three themes explain evolvability: "multiple genetic and non-genetic mechanisms to generate phenotypic diversity, robustness in genetic systems, and adaptive landscape topography". They also report "mounting evidence that evolvability can evolve".

**Masel & Trotter (2010).** "Robustness to mutation allows genetic variation to accumulate in a cryptic state". The robustness–evolvability link depends on "whether recombination rates are high or low". "In both cases, the evidence supports the claim that robustness promotes evolvability."

**Draghi et al. (2010).** Robustness "can either impede or facilitate adaptation, depending on the population size, the mutation rate and the structure of the fitness landscape". Neutral diversity accelerates adaptation "as long as the number of phenotypes accessible to an individual by mutation is smaller than the total number of phenotypes in the fitness landscape".

**Relevance.** Molecule-level measures (§2.5–2.8) are proxies. Whether they turn into evolvability depends on population size, mutation rate and landscape, which belong to level 1 and above. The Draghi condition, few phenotypes reachable from one molecule and many overall, is directly checkable with e_g and the total phenotype count.

### 2.16 Summary: measure × substrate

"Map" measures need sampling of genotype space; "scan" measures need the mutational neighbourhood of given molecules.

| Measure | Data needed | λ-terms | Strings / BFF | Programs (Avida) | NN weights | From Chemart outputs alone? |
|---|---|---|---|---|---|---|
| r_g (§2.5) | Molecule + mutation op + Φ (scan) | Yes; define lethal = non-normalizing/ill-formed | Yes | Done (Ofria 2002) | Yes (ε-ball; `neutral_fraction`) | No (needs reaction oracle) |
| P_same(k), epistasis (§2.6) | As above, k-fold | Yes | Yes | Done (Lenski 1999) | Yes (P_same(ε)) | No |
| e_g (§2.7) | Scan | Yes | Yes | Yes | Yes | No |
| ρ_p, ε_p (§2.8) | Neutral-set samples | Yes (random walks) | Yes | Done in LGP (Hu 2012) | Yes after discretization | No |
| f_p, bias, simplicity bias (§2.9) | Genotype sampler + Φ | Generator-dependent | Yes | Yes | Yes (weight prior) | Partly (soup sample only) |
| ρ_p vs f_p (§2.10) | Both of the above | Yes | Yes | Yes | Yes after discretization | No |
| Excess robustness (§2.11) | Population + neighbourhood scans | Yes | Yes | Yes | Yes | No |
| Navigability (§2.12) | Samples + scans + fitness assignment | Costly | Costly | Costly | Costly | No |
| C = ℓ − ΣH_i (§2.13) | Aligned population strings, or a scan | Needs positional encoding | Yes | Done (Adami 2000) | Only via −log f_p | **Yes** (population version) |
| Viability density (§2.14) | Sampler + viability predicate | Yes | Yes | Yes | Yes | No |
| Clark R, D sets (§2.3) | Interaction matrix | Yes | Yes | Yes | Yes | **Yes** |
| TSE R, D_N (§2.2) | Noise model of internal elements | No | No | No | Yes (hidden units) | No |

---

## 3. Brittleness in computational substrates and the known fixes

### 3.1 Tierra: Ray's own statement

Ray (1991), in "An approach to the synthesis of life":

> "Von Neuman type machine languages are considered to be 'brittle', meaning that the ratio of viable programs to possible programs is virtually zero. Any mutation or recombination event in a real machine code is almost certain to produce a non-functional program. The problem of brittleness can be mitigated by designing a virtual computer whose machine code is designed with evolution in mind."

Ray adds that Farmer & Belin suggested that overcoming this brittleness and "Discovering how to make such self-replicating patterns more robust so that they evolve to increasingly more complex states is probably the central problem in the study of artificial life." (Quoted by Ray; the original is not checked, see §6.)

Ray's fixes:

- A small instruction set with no numeric operands: "The Tierran language consists of 32 instructions, which can be represented by five bits, *operands included*." He contrasts Core War redcode, whose operand-inclusive set "works out to be about 10^11 in size".
- Addressing by complementary template: "molecule A presents a template on its surface which is complementary to some surface on B"; a jump searches "for the nearest occurrence of the complementary pattern."

### 3.2 Avida: measuring and designing evolvability

**Ofria, Adami & Collier (2002).**

- Redcode "does not survive mutations (i.e., it is extremely brittle)". Ray "recognized that the brittleness of redcode is due primarily to the argumented instruction set". In Avida's redcode, "over 99.7% of all nontrivial mutations are deleterious… Those few mutations that were not deleterious were almost entirely neutral."
- They compared five instruction sets.
  - Without templates, set III was "extremely inflexible… more akin to the redcode chemistry".
  - Direct-matching versus complement-matching templates changed how often lineages "lock in" a brittle length computation.
  - An 84-instruction set lagged in fitness.
- Conclusion: differences in evolvability are "attributable mainly to their robustness to mutations and the manner in which genome-size changes occur", and "the tendency of evolvability to go hand in hand with mutational robustness or neutrality" is thought to be "universal".

**Lenski et al. (1999).** Complex digital organisms are more robust to single mutations than simple ones, and epistasis is common (§2.6).

**Wilke et al. (2001).** At high mutation rates, genotypes "located in flatter regions of the fitness surface" beat faster replicators: "survival of the flattest". This is selection for robustness at a cost to replication rate.

**Bryson & Ofria (2013).** Tested six architectural features across seven environments. "multiple argument specification and separated I/O" helped most; most other changes had minimal systematic effect (paraphrase of the arXiv abstract). Banzhaf & Yamamoto (2015, §8.3.3, Fig. 8.1) summarize this study as the broadest comparison of instruction-set evolvability.

### 3.3 Neutrality in genetic programming, and counter-evidence on "programs are brittle"

**Banzhaf (1994).**

- Binary genotypes are mapped to program phenotypes through a repair map that guarantees feasibility.
- As a result, "multiple solutions in genotype space… map into one solution in phenotype space", so neutral variants are frequent and help maintain diversity.
- The metadata is confirmed; the abstract content comes from a search snippet.

**Yu & Miller (2001).** Explicit neutrality in Cartesian GP on a Boolean benchmark; "neutrality improves evolvability". This is from a secondary summary; I could not open the abstract.

**Hu, Payne, Banzhaf & Moore (2012).** Full characterization of genotype, phenotype and fitness networks of a compact Boolean LGP, with robustness and evolvability quantified at each level.

**Hu, Tomassini & Banzhaf (2020).**

- 6.4×10^13 genotypes encode 256 phenotypes; the most common phenotype (FALSE) has more than 10^8 of 10^9 sampled genotypes.
- Network evolvability measures are positively correlated with redundancy.
- Weighted eigenvector centrality best predicts phenotype evolvability, with a remaining discrepancy "resulted by the mutational bias led by robust genotypes".

**Schulte et al. (2014): the important counterweight.**

> "Although software is often viewed as brittle, with small changes leading to catastrophic changes in behavior, our results show surprising robustness in the face of random software mutations… over 30% of random mutations are neutral with respect to their test suite. The results hold across all classes of programs, for mutations at both the source code and assembly instruction levels."

**Relevance.** "Programs are brittle" is not a general law. Measured brittleness depends on:

- the mutation operator (Tierra-style operand mutations versus GP statement operators);
- dead or unexecuted code;
- the resolution of the phenotype (a competitive-fitness difference versus passing a test suite).

Any claim that λ-terms or BFF tapes are brittle needs its own measurement, with the operator and Φ matched to the NN chemistry.

### 3.4 Approximate ("fuzzy") binding inside artificial chemistries

**Stringmol** (Hickinbotham et al. 2010, spec v0.2). A "'soft' binding process, based on Smith-Waterman alignments": "The inexact alignments of sequences allows rates of binding and execution of the microprograms to evolve."

- Complement: a template symbol's complement is "13 letters downstream"; function codes are self-complementary.
- Binding probability: `P = 0` if the alignment length λ ≤ 3, otherwise `P = min(σ, λ − m)/(λ − m)`, where σ is the alignment score and m the single-mismatch penalty.
- Why the earlier rule `(σ/λ)^λ` was replaced: "A small number of mutations away from a perfect match can have a highly detrimental effect on long binding sites… There is also little redundancy at the top end of the bind probability profile - we would like there to be many ways of achieving a bind probability of 1."

This is explicit design *for* redundancy of the binding map, the exact move Marco's hypothesis recommends at the molecule level.

**Clark et al. (2011).** Removing Stringmol's binding degeneracy removes complex system-level artifacts (§2.3).

**Hickinbotham et al. (2016)**, the "Everything's Soft" principle: "Traditional programming languages are designed to be deterministic. Work is required to build degeneracy into the programs when we wish to emulate real biochemistry. There is little facility in existing languages for the sort of evolved changes observed in the phylogeny of proteins." Their conclusion: "the inexact string matching… raises the possibility for mutation to create a new molecular species that cannot be bound to by a parasite."

**Typogenetics**, as described in Banzhaf & Yamamoto (2015, pp. 204–205). Each typoenzyme has a "'tertiary' structure, which determines to which unit it can bind", and "The initial binding position of an enzyme depends in a nontrivial way on its coding sequence". Many amino-acid sequences therefore share one binding preference: a built-in many-to-one map.

**Lock-and-key ACs** (Banzhaf & Yamamoto 2015, §11.2, pp. 234–238):

- MCS.bl binds by pattern strings.
- Stringmol binds by alignment.
- SAC uses wildcards.
- Farmer et al.'s immune network binds when "the degree of complementarity between their respective strings is above a threshold", over all alignments.
- Conrad's enzymatic computer rests on lock-key "pattern recognition".

### 3.5 What Banzhaf & Yamamoto (2015) say in §8.3 and §8.3.3

**§8.3 frames the problem.** Beyond an error threshold there is "a viability threshold (qualitative): for instance, if most mutations are lethal… then the evolutionary process will be very slow and erratic."

**§8.3.2** lists "redundancy in the representation, either in the form of gene duplication and divergence, or in the form of alternative genotypes mapping to the same phenotype" as a technique for improving evolvability.

**§8.3.3, "Evolvability of ACs".**

- "the chemical programming language used must produce viable and fertile individuals with high probability. This is an area where there is currently no firm recipe for success."
- It then summarizes Ray, Ofria et al. (the 99.7%), Bryson & Ofria, and Suzuki's string-rewriting evolvability work.
- It ends: "there remains much to be done in regard to the examination of the evolutionary potential of artificial chemistries."
- §8.3.3 does **not** discuss AlChemy's or λ-calculus's evolvability.

**Chapter 17** (p. 349) quotes Conrad's trade-off principle and adds: "by making chemical computers programmable like electronic computers, their evolvability will likely be lost, and they could end up with the same brittleness that characterizes present-day computer hardware: a change in a single bit of instruction causes a disruption in the entire computation flow."

### 3.6 Has anyone argued that λ-calculus/AlChemy molecules are brittle and that this limits AlChemy?

**Not explicitly, as far as I could find.** What exists:

1. **General arguments about programmable code.**
   - Conrad (1985; trade-off principle).
   - Conrad (1990): "Biological structures that are characterized by a high degree of component redundancy and multiple weak interactions satisfy these conflicting pressures".
   - Banzhaf & Yamamoto (2015, ch. 17).
   - Hickinbotham et al. (2016) on the lack of degeneracy in deterministic languages.

   None of them names λ-calculus or AlChemy.
2. **Statements about AlChemy that are *not* about molecular brittleness.**
   - Fontana & Buss (1994): level-0 ensembles "are typically not robust towards functional perturbations. When a small number of random objects is introduced into the system, L0-ensembles typically collapse to a single replicator." (The scanned text reads "La"; I read it as L0.) This is fragility of the *organization*.
   - Level-1 organizations have "strong self-repair capabilities responsible for a robustness to perturbation".
   - Clark et al. (2011) attribute the level-0 collapse to AlChemy's "relatively simple binding".
   - Mathis et al. (2024) find "a surprising mix of dynamical robustness and fragility": organizations are "robust against collapse into trivial fixed-points, but… cannot be easily combined into higher order entities". This is a level-2 ceiling, attributed to nothing molecular.
3. **A point against the naive version of the premise.** Fontana & Buss (1994) chose λ-calculus *because* it captures "chemistry's diversity of equivalence classes, that many different reactants can yield the same stable product". λ-chemistry is therefore highly *redundant*, with infinitely many terms per normal form or function.

   **Implication.** The defensible form of Marco's premise is not "λ has no redundancy". It is "λ's redundancy is not *mutationally connected*": ρ_p ≈ f_p, few neutral neighbours, low r_g. That is testable with §2.5, §2.10 and §2.14 and, as far as I found, has not been tested. The Neural Molecules docstring states it as a premise: "Lambda terms are brittle -- one symbol changed and the evaluation path diverges". It would become the project's first empirical result rather than an assumption.

---

## 4. Neural networks as the substrate

### 4.1 The parameter–function map is many-to-one and biased

**Valle-Pérez, Camargo & Louis (2019).** They argue that "the parameter-function map of many DNNs should be exponentially biased towards simple functions". Evidence: Boolean-function models, and larger fully connected and convolutional networks on CIFAR10 and MNIST. P(f), the probability that parameters sampled from the initialization distribution produce f, is the NN analog of phenotype frequency f_p.

**Mingard et al. (2019, arXiv).** For a single-layer perceptron without bias, P(t) = 2^(−n) for 0 ≤ t < 2^n. This gives "a strong intrinsic a-priori bias towards individual functions with low entropy".

**Mingard et al. (2025, Nat Commun).** For deep tanh networks, "The simplicity bias in P(f) becomes weaker as the width σw of the Gaussian Ppar(σw) increases". They cite results that "in the fully chaotic limit, the prior over functions is fully uniform". ReLU networks' bias "barely changes with σw".

**Implication for Neural Molecules.** The weight-scale "redundancy knob" also moves the molecules along a simplicity-bias axis. The evidence is from deep networks; for a shallow 3→3→3 molecule this has to be measured.

- Small weights: strong bias, and frequent behaviours are simple.
- Large weights: weaker bias, closer to uniform.

A sweep that shows higher neutrality at some scale could just be showing more *simple* behaviours (§5).

### 4.2 Exact symmetries and what they do to a redundancy count

**Chen, Lu & Hecht-Nielsen (1993).** Tanh MLPs have "equioutput" transformations, meaning permutations and sign flips, that "form an algebraic group isomorphic to a direct product of Weyl groups". Each optimum has "large numbers of copies" that "all lie on the same sphere". For a hidden layer of h tanh units the group has h!·2^h elements (the hyperoctahedral group): **48 for h = 3.**

**Sussmann (1992)**, as stated by Farrugia-Roberts (2023): "two irreducible parameters are functionally equivalent if and only if they are related by simple operations of exchanging and negating the weights of hidden units."

**Farrugia-Roberts (2023).** With output weights a_i, input weights b_i and biases c_i, a single-hidden-layer tanh parameter is **reducible** iff at least one of these holds:

- (i) a_i = 0 for some i;
- (ii) b_i = 0 for some i;
- (iii) (b_i, c_i) = (b_j, c_j) for some i ≠ j;
- (iv) (b_i, c_i) = (−b_j, −c_j) for some i ≠ j.

Reducible parameters are a "vanishing minority". They have "richer functional equivalence classes", which are "piecewise-linear path-connected sets". The paper also gives a canonicalization algorithm. It is stated for scalar output; my inference is that the vector-output case extends it with a_i = 0 read as a zero vector.

**Implications for the redundancy ratio (inference).**

1. **At exact functional resolution, an irreducible Neural Molecule has 48 genotypes per function and no neutral network.** Every continuous neutral direction comes from:
   - the behaviour quantizer (finite probes, projection to k = 3 dimensions, cell size);
   - reducible or near-reducible molecules with dead units, constant units, or duplicated or negated units.

   Neutrality is therefore a joint property of the substrate *and* the observer's resolution. The species docstring already concedes this: species "are quantisation cells, not clusters" and "richness is resolution-dependent".
2. **Weight species should be counted on canonicalized weights.**
   - Sort hidden units after fixing a sign convention (for example, make each unit's first nonzero outgoing weight positive).
   - Otherwise any chemistry whose copying or transformation can produce permuted or sign-flipped copies will inflate n_weight_species by up to 48 with *no* functional content.
   - A random 3-D projection of θ is not invariant under these symmetries.
3. **Compare like with like against λ-terms.** In λ-terms α-equivalence is the trivial symmetry, and it is factored out by de Bruijn indices. The NN analog, the 48-element group, should be factored out the same way.

### 4.3 Mode connectivity: equal loss is not equal function

**Garipov et al. (2018).** "optima of complex loss functions… are connected by simple curves over which training and test accuracy are nearly constant."

**Draxler et al. (2018).** Paths between minima "are essentially flat in both the training and test landscapes".

**The caveat.** Farrugia-Roberts (2023) notes that "Garipov et al. (2018) observe functional non-equivalence in low-loss paths."

**Relevance.** These results show that the set of weights with the *same performance* is connected in large networks. That is the NN analog of a connected neutral network, relative to a coarse phenotype (task loss). They do not show that the same *function* is connected. For the chemistry, what matters is connectivity of the *behaviour-species* level sets. That can be measured directly with neutral random walks in weight space (§2.11).

### 4.4 Prior work that uses NNs as molecules

**Gabor et al. (2022, Artificial Life).** Self-replicating networks as fixpoints of weight space. "backpropagation turns out to be the natural way to navigate the space of network weights and allows non-trivial self-replicators to arise naturally". They include robustness analyses, "an extensive analysis of the occurrence of fixpoint weight configurations" with their attractor basins, and AC environments of multiple networks (Crossref abstract).

**Chang & Lipson (2018, ALIFE).** "The network replicates itself by learning to output its own weights". There is "a trade-off between the network's ability to classify images and its ability to replicate".

**What neither does.** Neither measures genotype–phenotype-map robustness or evolvability in the sense of §2. Neural Molecules' neutrality diagnostics would be new there. That is a gap statement, based on the abstracts only.

---

## 5. The tension: does more redundancy mean more or less complexity?

### 5.1 Evidence that redundancy lowers complexity or information

**Information accounting.**

- Adami et al. (2000): neutral sites contribute only entropy, and C = ℓ − ΣH_i (§2.13).
- Ofria et al. (2002): robustness is "the fraction of its length that is impervious to mutations and thus carries no information".
- At fixed length, more robustness means less information per molecule.

**Bias toward simple phenotypes.**

- Dingle et al. (2018): high-probability outputs must be simple, P(x) ≲ 2^(−aK̃(x)−b).
- Valle-Pérez et al. (2019) and Mingard et al. (2025): the NN parameter–function map is simplicity-biased.
- Schaper & Louis (2014), the "arrival of the frequent": "frequent phenotypes (with larger F(p)) can fix in a population even when alternative, but less frequent, phenotypes with much higher fitness are potentially accessible. In other words, if the fittest never 'arrive' on the timescales of evolutionary change, then they can't fix."

**Selection for robustness at the expense of function.**

- Wilke et al. (2001): survival of the flattest.
- Milano, Pagliuca & Nolfi (2019): "the competition for robustness to mutations… leads to the selection of phenotypically simple but low evolvable circuits. These circuits achieve robustness by minimizing the number of functional genes rather than by relying on redundancy or degeneracy."
- Ancel & Fontana (2000): selection that reduces plasticity causes "a dramatic loss of variability (and hence a loss of evolvability) to the point of lock-in"; populations are trapped "in regions where most genetic variation is phenotypically neutral. We call this phenomenon neutral confinement."

**Antiredundancy wins in large populations.** Krakauer & Plotkin (2002): assuming a cost of redundancy, "large populations will evolve antiredundant mechanisms for removing mutants and thereby bolster the robustness of wild-type genomes; whereas small populations will evolve redundancy".

**Implication (inference).** AC soups are large, mutation-rich and dominated by replication. Brittleness there may be *selected*, not only inherited from the substrate. That is worth checking before blaming the substrate.

### 5.2 Evidence that redundancy, or degeneracy, raises complexity and evolvability

**Degeneracy goes with complexity.**

- Tononi et al. (1999): selecting for degeneracy yields high neural complexity; purely redundant systems have low degeneracy.
- Whitacre (2010) reads this as "highly redundant (non-degenerate) systems were naturally robust but never hierarchically complex". That is Whitacre's paraphrase, not Tononi's wording.

**Degeneracy, not redundancy, gives evolvability.**

- Whitacre & Bender (2010): "purely redundant systems have remarkably low evolvability while degenerate, i.e. partially redundant, systems tend to be orders of magnitude more evolvable. Surprisingly, the magnitude of observed variation in evolvability can neither be explained by differences in the size nor the topology of the neutral networks."
- Whitacre (2010), summarizing that work: "only systems with high levels of degeneracy exhibited a positive relationship between neutral network size, robustness, and evolvability", and systems of redundant proteins "were mutationally robust but greatly restricted in the number of unique phenotypes accessible from a neutral network, i.e. they were not evolvable."

**Artificial-chemistry evidence.** Clark et al. (2011): removing binding degeneracy in Stringmol removes hypercycles, macro-mutations, sweeps and parasite evasion.

**Phenotype-level robustness helps evolvability.**

- Wagner (2008): phenotype robustness promotes phenotype evolvability.
- Ahnert (2017): a positive robustness–evolvability correlation is a shared property of GP maps.
- Huynen et al. (1996): diffusion on neutral networks "enables the search of vast areas in genotype space while still preserving the dominant phenotype".
- Greenbury et al. (2022): neutral networks make landscapes navigable.
- Masel & Trotter (2010): "the evidence supports the claim that robustness promotes evolvability".

**Classic design principles.**

- Kirschner & Gerhart (1998): weak linkage, compartmentation and redundancy confer evolvability.
- Conrad (1990): "Organizations that are complex in terms of numbers of components and interactions are more likely to meet the peak-climbing condition, but less likely to meet the stability condition. Biological structures that are characterized by a high degree of component redundancy and multiple weak interactions satisfy these conflicting pressures."

**Digital organisms.**

- Lenski et al. (1999): complex organisms are more robust.
- Ofria et al. (2002): "the tendency of evolvability to go hand in hand with mutational robustness or neutrality".

### 5.3 What the literature concludes

1. **The trade-off is real at the genotype level and dissolves at the phenotype level** (Wagner 2008; Ahnert 2017). A single brittle molecule is maximally "evolvable" in the local sense, because every mutation gives something new, but it cannot keep its phenotype. A robust *phenotype* can do both, *if* its neutral network is large, connected and borders many phenotypes.
2. **Redundancy is necessary but not sufficient.** It must be mutationally connected (ρ_p ≫ f_p; Greenbury et al. 2016; Hu et al. 2020). It should be degenerate rather than pure duplication (Whitacre & Bender 2010; Clark et al. 2011). The phenotypes it favours should not all be trivially simple (Dingle et al. 2018; Milano et al. 2019).
3. **Whether robustness helps depends on population size, mutation rate and landscape.** It helps adaptation when few phenotypes are reachable from one genotype relative to the total (Draghi et al. 2010). Populations exploit neutral networks when Mμ is large (van Nimwegen et al. 1999; Huynen et al. 1996). But large populations may evolve antiredundancy (Krakauer & Plotkin 2002), and selection for robustness can lock populations in (Ancel & Fontana 2000).
4. **"More redundancy means more complexity" is not a conclusion any source supports.** The supported version is: *degenerate, connected neutrality at the phenotype level* increases evolvability, and through it the potential for complexity. Excess redundancy that buys robustness by making phenotypes simpler or sites inert lowers information and evolvability.

### 5.4 What this means for the redundancy hypothesis (inference)

To test the hypothesis the literature's way, measure four things on each chemistry, with Φ = behaviour on a shared probe battery:

- (a) r_g and viability density, the brittleness claim itself (§2.5, §2.14);
- (b) ρ_p against f_p, whether redundancy is connected (§2.10);
- (c) ε_p and e_g, whether robustness is buying access to *new* phenotypes (§2.7–2.8);
- (d) the complexity of the frequent phenotypes (K̃, or C from §2.13), whether redundancy is just buying simplicity.

The hypothesis predicts that NN molecules have higher ρ_p/f_p and higher ε_p than λ-terms or BFF tapes *without* lower phenotype complexity, and that this is what moves the organizational ceiling. Clark et al.'s degeneracy set on the reaction matrix (§2.3) is the cheapest first step, because it runs on Chemart's existing outputs.

---

## 6. References

Status key:

- **VERIFIED**: authors, year, title and venue confirmed this session, and the definition or quote taken from the abstract or full text.
- **VERIFIED (metadata)**: bibliographic details confirmed, content only from a secondary summary.
- **UNVERIFIED**: not independently confirmed.

The brief named about 35 works. The additions marked ⁺ change the conclusions: the λ-brittleness question, NN symmetries and the tension. That is why the list is longer than 30.

- Adami C, Ofria C, Collier TC (2000). Evolution of biological complexity. *PNAS* 97(9):4463–4468. doi:10.1073/pnas.97.9.4463. arXiv:physics/0005074. **VERIFIED** (full text read).
- Ahnert SE (2017). Structural properties of genotype–phenotype maps. *J R Soc Interface* 14:20170275. doi:10.1098/rsif.2017.0275. **VERIFIED**.
- ⁺Altenberg L (1994). The evolution of evolvability in genetic programming. In Kinnear KE (ed), *Advances in Genetic Programming*, MIT Press, ch. 3, pp. 47–74. Author's PDF: https://dynamics.org/Altenberg/FILES/LeeEEGP.pdf. **VERIFIED** (author's PDF read; used only for the evolvability definition).
- ⁺Ancel LW, Fontana W (2000). Plasticity, evolvability, and modularity in RNA. *J Exp Zool (Mol Dev Evol)* 288(3):242–283. doi:10.1002/1097-010X(20001015)288:3<242::AID-JEZ5>3.0.CO;2-O. **VERIFIED** (abstract).
- Banzhaf W (1994). Genotype-phenotype-mapping and neutral variation — a case study in genetic programming. In *PPSN III*, LNCS 866:322–332, Springer. doi:10.1007/3-540-58484-6_276. **VERIFIED (metadata)**; abstract content from a search snippet.
- Banzhaf W, Yamamoto L (2015). *Artificial Chemistries*. MIT Press. §8.3 (pp. 165–168), §10.5.2 (pp. 204–205), §11.2 (pp. 234–238), ch. 17 (p. 349). **VERIFIED** (read from Marco's local PDF).
- Bryson DM, Ofria C (2013). Understanding evolutionary potential in virtual CPU instruction set architectures. *PLoS ONE* 8(12):e83242. doi:10.1371/journal.pone.0083242. arXiv:1309.0719. **VERIFIED** (abstract).
- Chang O, Lipson H (2018). Neural network quine. *Proc. ALIFE 2018*, MIT Press, pp. 234–241. doi:10.1162/isal_a_00049. arXiv:1803.05859. **VERIFIED**.
- Chen AM, Lu H, Hecht-Nielsen R (1993). On the geometry of feedforward neural network error surfaces. *Neural Computation* 5(6):910–927. doi:10.1162/neco.1993.5.6.910. **VERIFIED** (abstract). The h!·2^h group order is the standard order of the hyperoctahedral group; the abstract states only "direct product of Weyl groups".
- ⁺Clark E, Nellis A, Hickinbotham S, Stepney S, Clarke T, Pay M, Young P (2011). Degeneracy enriches artificial chemistry binding systems. In *Advances in Artificial Life, ECAL 2011*, MIT Press, pp. 133–140. https://www-users.york.ac.uk/~ss44/bib/ss/nonstd/ecal11-99.htm. **VERIFIED** (full text read).
- ⁺Conrad M (1985). On design principles for a molecular computer. *Communications of the ACM* 28(5):464–480. doi:10.1145/3532.3533. **VERIFIED (metadata)**; the trade-off quote is taken from Banzhaf & Yamamoto (2015, p. 349).
- ⁺Conrad M (1990). The geometry of evolution. *BioSystems* 24:61–81. doi:10.1016/0303-2647(90)90030-5. **VERIFIED** (abstract).
- Dingle K, Camargo CQ, Louis AA (2018). Input–output maps are strongly biased towards simple outputs. *Nature Communications* 9:761. doi:10.1038/s41467-018-03101-6. **VERIFIED** (full text).
- Draghi JA, Parsons TL, Wagner GP, Plotkin JB (2010). Mutational robustness can facilitate adaptation. *Nature* 463:353–355. doi:10.1038/nature08694. **VERIFIED** (abstract).
- Draxler F, Veschgini K, Salmhofer M, Hamprecht FA (2018). Essentially no barriers in neural network energy landscape. *ICML 2018*, PMLR 80:1308–1317. arXiv:1803.00885. **VERIFIED** (abstract).
- Edelman GM, Gally JA (2001). Degeneracy and complexity in biological systems. *PNAS* 98(24):13763–13768. doi:10.1073/pnas.231499798. **VERIFIED** (abstract).
- Farmer JD, Belin A (1992). Artificial life: the coming evolution. Reprinted in *Artificial Life II*. **UNVERIFIED**: seen only in Ray (1991)'s text and reference list.
- ⁺Farrugia-Roberts M (2023). Functional equivalence and path connectivity of reducible hyperbolic tangent networks. *NeurIPS 36*, pp. 79502–79517. doi:10.52202/075280-3479. arXiv:2305.05089. **VERIFIED** (full text read).
- ⁺Fontana W, Buss LW (1994). "The arrival of the fittest": toward a theory of biological organization. *Bulletin of Mathematical Biology* 56(1):1–64. doi:10.1007/BF02458289. **VERIFIED** (SFI working-paper text 93-09-055 read).
- Gabor T, Illium S, Zorn M, Lenta C, Mattausch A, Belzner L, Linnhoff-Popien C (2022). Self-replication in neural networks. *Artificial Life* 28(2):205–223. doi:10.1162/artl_a_00359. **VERIFIED** (Crossref abstract).
- Garipov T, Izmailov P, Podoprikhin D, Vetrov D, Wilson AG (2018). Loss surfaces, mode connectivity, and fast ensembling of DNNs. *NeurIPS 2018*. arXiv:1802.10026. **VERIFIED** (abstract).
- Greenbury SF, Schaper S, Ahnert SE, Louis AA (2016). Genetic correlations greatly increase mutational robustness and can both reduce and enhance evolvability. *PLoS Computational Biology* 12(3):e1004773. doi:10.1371/journal.pcbi.1004773. **VERIFIED** (abstract + article text).
- Greenbury SF, Louis AA, Ahnert SE (2022). The structure of genotype–phenotype maps makes fitness landscapes navigable. *Nature Ecology & Evolution* 6:1742–1752. doi:10.1038/s41559-022-01867-z. Preprint doi:10.1101/2021.10.11.463990. **VERIFIED** (abstract; navigability definition and numbers from the bioRxiv full text).
- ⁺Hickinbotham S, Clark E, Stepney S, Clarke T, Nellis A, Pay M, Young P (2010). Specification of the Stringmol chemical programming language, version 0.2. Tech. Rep. YCS-2010-458, University of York. https://www.cs.york.ac.uk/library/reports/2010/YCS/458/YCS-2010-458.pdf. **VERIFIED** (full text read). The published Stringmol paper is Hickinbotham et al., "Molecular microprograms", ECAL 2009, LNCS 5777:297–304 (metadata only).
- ⁺Hickinbotham S, Clark E, Nellis A, Stepney S, Clarke T, Young P (2016). Maximizing the adjacent possible in automata chemistries. *Artificial Life* 22(1):49–75. doi:10.1162/ARTL_a_00180. **VERIFIED** (accepted manuscript read).
- Hu T, Payne JL, Banzhaf W, Moore JH (2012). Evolutionary dynamics on multiple scales: a quantitative analysis of the interplay between genotype, phenotype, and fitness in linear genetic programming. *Genetic Programming and Evolvable Machines* 13:305–337. doi:10.1007/s10710-012-9159-4. **VERIFIED (metadata)**; content from a search summary of the abstract.
- ⁺Hu T, Tomassini M, Banzhaf W (2020). A network perspective on genotype–phenotype mapping in genetic programming. *Genetic Programming and Evolvable Machines* 21:375–397. doi:10.1007/s10710-020-09379-0. **VERIFIED** (full text read).
- Huynen MA, Stadler PF, Fontana W (1996). Smoothness within ruggedness: the role of neutrality in adaptation. *PNAS* 93(1):397–401. doi:10.1073/pnas.93.1.397. **VERIFIED** (abstract).
- Kirschner M, Gerhart J (1998). Evolvability. *PNAS* 95(15):8420–8427. doi:10.1073/pnas.95.15.8420. **VERIFIED** (abstract).
- ⁺Krakauer DC, Plotkin JB (2002). Redundancy, antiredundancy, and the robustness of genomes. *PNAS* 99(3):1405–1409. doi:10.1073/pnas.032668599. **VERIFIED** (abstract).
- Lenski RE, Ofria C, Collier TC, Adami C (1999). Genome complexity, robustness and genetic interactions in digital organisms. *Nature* 400:661–664. doi:10.1038/23245. **VERIFIED** (abstract).
- Manrubia S, Cuesta JA, Aguirre J, Ahnert SE, Altenberg L, Cano AV, Catalán P, Diaz-Uriarte R, Elena SF, García-Martín JA, Hogeweg P, Khatri BS, Krug J, Louis AA, Martin NS, Payne JL, Tarnowski MJ, Weiß M (2021). From genotypes to organisms: state-of-the-art and perspectives of a cornerstone in evolutionary dynamics. *Physics of Life Reviews* 38:55–106. doi:10.1016/j.plrev.2021.03.004. arXiv:2002.00363. **VERIFIED** (abstract only; not mined for measures beyond Ahnert 2017).
- Masel J, Trotter MV (2010). Robustness and evolvability. *Trends in Genetics* 26(9):406–414. doi:10.1016/j.tig.2010.06.002. **VERIFIED** (abstract).
- ⁺Mathis C, Patel D, Weimer W, Forrest S (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos* 34(9):093142. https://pubs.aip.org/aip/cha/article/34/9/093142/3314760. arXiv:2408.12137. **VERIFIED** (abstract; venue from the publisher listing; DOI not retrieved).
- ⁺Milano N, Pagliuca P, Nolfi S (2019). Robustness, evolvability and phenotypic complexity: insights from evolving digital circuits. *Evolutionary Intelligence* 12(1):83–95. doi:10.1007/s12065-018-00197-z. arXiv:1712.04254. **VERIFIED** (abstract).
- Mingard C, Skalse J, Valle-Pérez G, Martínez-Rubio D, Mikulik V, Louis AA (2019). Neural networks are a priori biased towards Boolean functions with low entropy. arXiv:1909.11522. **VERIFIED** (abstract).
- Mingard C, Rees H, Valle-Pérez G, Louis AA (2025). Deep neural networks have an inbuilt Occam's razor. *Nature Communications* 16:220. doi:10.1038/s41467-024-54813-x. arXiv:2304.06670. **VERIFIED** (full text; σw quotes checked in the PDF).
- Ofria C, Adami C, Collier TC (2002). Design of evolvable computer languages. *IEEE Transactions on Evolutionary Computation* 6(4):420–424. doi:10.1109/TEVC.2002.802442. **VERIFIED** (full text read).
- Payne JL, Wagner A (2019). The causes of evolvability and their evolution. *Nature Reviews Genetics* 20:24–38. doi:10.1038/s41576-018-0069-z. **VERIFIED** (abstract).
- Ray TS (1992; preprint dated 1991). An approach to the synthesis of life. In Langton CG, Taylor C, Farmer JD, Rasmussen S (eds), *Artificial Life II*, Santa Fe Institute Studies in the Sciences of Complexity, Addison-Wesley, pp. 371ff. Author's text: http://tomray.me/pubs/alife2/tierra.tex. **VERIFIED** (author's text read). The SFI volume number differs between sources and is not given here.
- ⁺Schaper S, Louis AA (2014). The arrival of the frequent: how bias in genotype–phenotype maps can steer populations to local optima. *PLoS ONE* 9(2):e86635. doi:10.1371/journal.pone.0086635. **VERIFIED** (abstract).
- ⁺Schulte E, Fry ZP, Fast E, Weimer W, Forrest S (2014). Software mutational robustness. *Genetic Programming and Evolvable Machines* 15(3):281–312. doi:10.1007/s10710-013-9195-8. arXiv:1204.4224. **VERIFIED** (abstract).
- Schuster P, Fontana W, Stadler PF, Hofacker IL (1994). From sequences to shapes and back: a case study in RNA secondary structures. *Proc R Soc Lond B* 255(1344):279–284. doi:10.1098/rspb.1994.0040. **VERIFIED** (abstract).
- ⁺Suzuki H (2003). An example of design optimization for high evolvability: string rewriting grammar. *BioSystems* 69(2–3):211–221. doi:10.1016/S0303-2647(02)00138-7. **VERIFIED** (abstract).
- ⁺Sussmann HJ (1992). Uniqueness of the weights for minimal feedforward nets with a given input-output map. *Neural Networks* 5(4):589–593. doi:10.1016/S0893-6080(05)80037-1. **VERIFIED (metadata)**; the theorem as stated by Farrugia-Roberts (2023). I did not read the original.
- Tononi G, Sporns O, Edelman GM (1999). Measures of degeneracy and redundancy in biological networks. *PNAS* 96(6):3257–3262. doi:10.1073/pnas.96.6.3257. **VERIFIED** (abstract + formulas from the PMC full text, PMC15929).
- Valle-Pérez G, Camargo CQ, Louis AA (2019). Deep learning generalizes because the parameter-function map is biased towards simple functions. *ICLR 2019*. arXiv:1805.08522. **VERIFIED** (abstract).
- van Nimwegen E, Crutchfield JP, Huynen M (1999). Neutral evolution of mutational robustness. *PNAS* 96(17):9716–9720. doi:10.1073/pnas.96.17.9716. arXiv:adap-org/9903006. **VERIFIED** (full text; numbers checked).
- Wagner A (2005). *Robustness and Evolvability in Living Systems*. Princeton University Press. ISBN 0-691-12240-7. **VERIFIED (metadata)**; definitions here come from the same author's 2005 FEBS paper (next entry), not the book.
- ⁺Wagner A (2005). Robustness, evolvability, and neutrality. *FEBS Letters* 579(8):1772–1778. doi:10.1016/j.febslet.2005.01.063. **VERIFIED** (abstract).
- Wagner A (2008). Robustness and evolvability: a paradox resolved. *Proc R Soc B* 275(1630):91–100. doi:10.1098/rspb.2007.1137. **VERIFIED** (abstract + PMC text for the definitions).
- Whitacre JM (2010). Degeneracy: a link between evolvability, robustness and complexity in biological systems. *Theoretical Biology and Medical Modelling* 7:6. doi:10.1186/1742-4682-7-6. **VERIFIED** (abstract + PMC text).
- Whitacre J, Bender A (2010). Degeneracy: a design principle for achieving robustness and evolvability. *Journal of Theoretical Biology* 263(1):143–153. doi:10.1016/j.jtbi.2009.11.008. **VERIFIED** (abstract).
- Wilke CO, Wang JL, Ofria C, Lenski RE, Adami C (2001). Evolution of digital organisms at high mutation rates leads to survival of the flattest. *Nature* 412:331–333. doi:10.1038/35085569. **VERIFIED** (abstract).
- Yu T, Miller J (2001). Neutrality and the evolvability of Boolean function landscape. In *EuroGP 2001*, LNCS 2038:204–217, Springer. doi:10.1007/3-540-45355-5_16. **VERIFIED (metadata)**; the finding "neutrality improves evolvability" is from a secondary summary.
