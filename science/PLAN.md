# Science plan: neutral networks of molecules across chemistries

*A living document. Started 2026-09-23. Changes are logged at the bottom.*

## Why this exists

The literature review in `../research/measures-literature.md` reached three conclusions:

- The redundancy hypothesis is really about **neutrality and degeneracy of the structure → behaviour map**: whether swapping one symbol leaves what a molecule does unchanged.
- λ-calculus is redundant (many terms share a normal form) but probably not *local*, and nobody has measured it.
- The ceiling between level 1 (organisations) and level 2 (organisations of organisations) may be as much about missing units as about brittle molecules.

This plan measures the first two directly. It compares two artificial chemistries with real chemistry.

The unifying object is **the neutral network of a molecule, defined by what it does to, and receives from, other molecules**.

The work lives in `science/`, separate from the platform. It imports `chemart` but never modifies it.

## Decisions so far

| date | decision |
|---|---|
| 2026-09-23 | **AlChemy** runs on the Mathis group's published Rust engine, `AgentElement/functional-supercollider`. It is the reimplementation cited by Vimal, Mathis, Weimer & Forrest 2025 (arXiv:2509.03534). The older fork `colemathis/alchemy-reimplemented` is kept for reference only. |
| 2026-09-23 | **Stringmol** runs on chemart's exact Python port. Dynamics move to the upstream C++ only if Python is too slow (decided on day 1 of Phase 2). |
| 2026-09-23 | **Real chemistry:** MCM v3.3.1 (atmospheric, abiotic) and Rhea (enzymatic). |
| 2026-09-23 | **Code:** self-contained in `science/`. It reuses chemart's functions and vendors the Rust engine. |
| 2026-09-23 | **First result:** the probability that level-1 organisations collapse, as a function of mutation rate. |

## 0. Common frame

### 0.1 Terminology

- **Neutral mutation:** a one-symbol (one-atom) edit that leaves behaviour unchanged. This is Marco's "swapping one symbol doesn't change the behaviour".
- **Robustness r:** the fraction of a molecule's point mutations that are neutral.
- **Degeneracy:** different structures doing the same job, including *partial* overlap: the same outcome with some partners, a different one with others (Edelman & Gally 2001). Degeneracy describes the chemistry's map; a neutral mutant and its parent are a degenerate pair.
- **Junk against degeneracy:** a mutation can be neutral in two ways.
  - *Junk:* the site is never used. Examples: an erased subterm, a Stringmol region that is never executed or aligned, a spectator atom.
  - *Degeneracy:* the site is used, and tolerates change. This is the codon-like case, and **it is the main number**.
- **Locality:** neutral variants are each other's neighbours, so the neutral network is large and connected rather than scattered.

### 0.2 One setup for every chemistry

| | AlChemy (Rust) | Stringmol (chemart port) | MCM / Rhea |
|---|---|---|---|
| genotype | closed λ normal form (`lambda_calculus::Term` = `Var(usize)` with 1-based de Bruijn indices, or `Abs`, or `App`) | a string over 33 symbols | a molecular graph (RDKit canonical SMILES; Rhea neutralised and tautomer-canonical) |
| **uniform** operator, used for cross-chemistry comparison | one-node edit: M1, rebind a variable index | uniform substitution | one-atom edit: element swap, H↔CH₃, H↔OH, bond order; for Rhea also a stereocentre flip |
| native operator, reported separately | — | loop-neighbour substitution: the copy error, close to binding-neutral by design | — |
| oracle | the vendored Rust `collide` | `Machine.react(a, b, "possible", budget)` and `bind_probability` | the observed network. The behaviour of a molecule that is not in the network is *unknown*, not "different". |

**Behaviour Φ(m)** is m's outcome with each partner, in both roles: m acting, and m acted on. It is recorded at three resolutions:

- **class:** reacts or not, plus the failure reason (limit, filtered copy, identity, free variables, too large);
- **coarse product:** up to η for λ, the product multiset for Stringmol, the reaction class for real chemistry;
- **exact product.**

**Partners** come in three sets:

- the molecule's own organisation;
- members of *other* organisations, which is what level-2 merging needs;
- a frozen reference probe set.

Probe-saturation curves check that each set is large enough.

### 0.3 Molecule-level measures

- **r, ω, D.** r is the neutral fraction, reported against molecule size. ω is the graded overlap, the fraction of partners with the same outcome. D is the graded product distance: the mean normalised tree-edit distance between the products of m and m′ (Levenshtein distance for Stringmol).
- **Lethality**, broken down by failure reason for each partner.
- **The junk/degeneracy split.** The λ test: replace a site with a fresh free variable. The site is unused if no product contains that variable.
- **Carry-through:** is (m′)x an edit of (m)x? In AlChemy this is the only route to heredity at the organisation level.
- **P_same(k) curves:** the probability that behaviour is unchanged after k random mutations.
- **The locality ratio**, one estimator used for every chemistry:
  **L = P(same Φ | one-edit pair) / P(same Φ | size-matched random pair).**
- **The Greenbury test, ρ_p against f_p** (Greenbury et al. 2016). This compares how robust a behaviour is (ρ_p) with how common it is (f_p), and it is run only where f_p can be estimated.
  - M1 never changes a term's skeleton. So f_p comes from **exhaustively enumerating every index assignment of small skeletons** (Πdᵢ ≲ 10⁵), which gives exact f_p, ρ_p and ε_p.
  - It is not run where f_p is ill-defined: for exact-product phenotypes, for Stringmol's local sampler, or on observed networks.
- **Clark et al.'s (2011) redundancy and degeneracy sets**, computed on every observed interaction table.

### 0.4 Organisation-level measures (within one chemistry only)

- **Similarity uses our own definitions.** The engine's `jacard_index` is Σmin/(|A|+|B|), so identical soups score 0.5. We use:
  - set Jaccard, as Mathis et al. did, above an abundance cutoff;
  - abundance-weighted Σmin/N.
- **The primary collapse criterion is quench-recovery.** After T turnovers at rate μ, set μ = 0 for 20 turnovers, then classify the outcome:

  | outcome | criterion |
  |---|---|
  | **collapsed** | the reference core is no longer produced from inside the set (windowed self-maintenance), or the soup is trivial or inactive |
  | **shifted** | the soup is in a different non-trivial state that maintains itself |
  | **survived** | anything else |

  Collapse is reported as the excess over the μ = 0 rate. Similarity is a secondary measure.
- **Time is counted in productive reactions, not collisions.** One turnover is M productive reactions.
- **Evolvability is measured alongside robustness:**
  - P(shift);
  - the number of distinct self-maintaining states reached;
  - whether shifts simplify the organisation;
  - the rate of compositional drift at a sub-threshold μ, compared with μ = 0.
- **Collapse anatomy:** which of keystone loss, erosion, takeover or freezing happened.

### 0.5 Comparability rules

- **At the molecule level, compare across chemistries** using the uniform operators, L, the junk/degeneracy split, r against size, and the observed-pair estimator.
  - The observed-pair estimator uses pairs that both occur and differ by one edit.
  - Its bias is calibrated on the AlChemy and Stringmol species sets, where the full oracle is available.
- **At the organisation level, compare only within a chemistry.** μ₅₀ is never put in the same column for different chemistries, because their organisations are different kinds of thing:
  - AlChemy L1 organisations are production networks with no inheritance, since copy actions are filtered.
  - Stringmol populations are replicator ecologies (Fontana & Buss level 0) whose copy errors *are* inherited.

## 1. Study 1: AlChemy (Rust engine)

### 1.1 Engine

- **Source.** Vendor `AgentElement/functional-supercollider` (GPL-3.0) at a **pinned commit** into `alchemy/engine/`, with a `VENDOR.md`. Diff it against the older fork `colemathis/alchemy-reimplemented` and note the differences.
- **Provenance.** Mathis et al.'s 2024 paper used the original C++ code (`mathis-group/AlChemy`). The 2025 paper used the Rust engine with 8000 reduction steps, a limit of 1000 vertices, 5000–6000 expressions and 10⁶ collisions. **No published validation of the Rust engine against the original exists.**
- **Our crate**, `alchemy/mut/`, builds on the engine's `lib.rs` and provides:
  - `mutate`: M1 and N1;
  - `scan`: writes Φ as JSONL;
  - `run`: a soup with
    - product mutation, in the order reduce → mutate → re-reduce → filter;
    - replace-semantics injection (`perturb()` appends, which would grow the soup);
    - snapshots streamed to disk (`simulate_and_record` keeps them in RAM);
    - a guard on product size.
- **Orchestration and analysis** are Python scripts in `alchemy/`, calling the crate through subprocess with JSON in and JSONL out.
- **Configuration matches the paper's L1:** copy actions are barred and the reduction cutoff is 500. The Rust default `discard_identity` was not used by the paper; it is set to false and tested as a factor.

### 1.2 Operators

| operator | edit | notes |
|---|---|---|
| **M1 rebind** (primary) | change one variable's de Bruijn index to another valid binder | the result stays closed and normal |
| **N1** | insert or delete a vacuous λ, only at non-operator positions | inserting above an operator amounts to deleting the operand, so it is excluded |
| M3, M4 | leaf → subterm; application swap | deferred |

### 1.3 Experiments, in order

**A. Engine gate.** This gate blocks further work, but only on engine correctness.

- Unit tests: M1 and N1 results are closed, normal and canonical.
- Cross-check 10⁴ collisions, on random and on evolved pairs, against chemart's Python AlChemy. Normal forms are unique, so the products must be identical whenever both engines terminate.
- Benchmark collisions per second.
- A qualitative comparison with Mathis et al., reported but *not* blocking:
  - their Fig. 3B stability runs, with the arms "seed" and "perturb";
  - identity flooding at **L0**. At L1 the identity's products are filtered as copies, so flooding would only dilute the soup.

**B. Organisation library.**

- Two generators: Rust `BTreeGen`, and `fontana_generator.py` piped to stdin (Rust's `FontanaGen` is a stub).
- 75 seeds per generator, 10⁶ collisions per seed, M = 1000.
- Non-triviality and stationarity rules are declared in advance. Organisations that recur across seeds are merged.
- 20% of the library is split off as a **pilot set**.

**C. First result: collapse against mutation rate.**

- μ ∈ {0, 10⁻³, 3·10⁻³, 10⁻², 3·10⁻², 0.1, 0.3}.
- **Three arms at each μ.** Each replaces the product of a reaction with probability μ:

  | arm | the product is replaced by | what it shows |
  |---|---|---|
  | **mutant** | an M1 point mutant of the product | the test |
  | **inert** | a placeholder that never reacts | pure loss of production |
  | **random** | a term from the generator | any perturbation |

  Reading the arms:
  - mutant μ₅₀ below inert: mutants actively harm the organisation;
  - mutant μ₅₀ above inert: neutrality buffers it;
  - mutant μ₅₀ above random: mutation is local.
- 5 replicates × 100 turnovers per (organisation, μ, arm), plus 20 baseline runs per organisation. Quench-recovery judges collapse.
- Outputs:
  - P_collapse(μ) with confidence intervals;
  - μ₅₀ for each arm;
  - P(shift);
  - the anatomy of each collapse.
- The pilot set runs first and fixes the μ grid and T. The confirmatory set runs after registration.

**D. Molecular scan-lite.**

- Molecules: every library member, plus two nulls: size-matched generator terms, and first-generation products.
- All M1 and N1 edits of each molecule, against the three partner sets of 0.2.
- Outputs: r, ω, D, failure classes, the junk/degeneracy split, carry-through, L, and the exhaustive enumeration of small skeletons.
- This supplies r̄ for P1.5 and P1.6.

**E. Knock-in and knockout assays** on 20 organisations. These mirror Mathis's L1 protocol (10 molecules, 10⁵ collisions, 7 repeats).

| arm | what is done |
|---|---|
| knock-in | every copy of a member m is replaced by a mutant m′ |
| knockout | m is replaced by an inert term |
| sham | k copies of the unmutated parent are added |
| dose series | k ∈ {10, 100, 300} |

- Harm = knock-in − knockout.
- Integration is predicted without dynamics: is m′ in the closure of the organisation ∪ {m′}? The dynamics then test the prediction.
- These arms replace a single-mutant injection. At L1 a single mutant is almost always diluted out whether or not it is neutral, so that assay would have no power.

### 1.4 Predictions

These are registered after the pilot and before the confirmatory runs.

- **P1.1** M1 mutations are mostly non-neutral at exact resolution, and L is near 1 at class resolution.
- **P1.2** *(exploratory)* Members differ in r from size-matched random terms. No van Nimwegen mechanism applies, because nothing replicates.
- **P1.3** The collapse hazard rises more than 10× within half a decade of μ. μ₅₀ is stable across T ∈ {25, 50, 100, 200} turnovers and M ∈ {300, 1000, 3000}.
- **P1.4** Mutant μ₅₀ is below inert μ₅₀: brittle mutants actively harm.
- **P1.5** Across organisations, the slope of log μ₅₀ on log(1 − r̄) is −1.
  - r̄ is flux-weighted and measured against the organisation's own members, with the same operator as experiment C.
  - Rival explanation: network redundancy, i.e. the number of distinct producer pairs per member.
  - Interventional version: two operators with different r̄ applied to the same organisations give a μ₅₀ ratio of (1 − r̄₁)/(1 − r̄₂).
- **P1.6** Responses are bimodal: either no lasting change or loss of the organisation, with little neutral drift. This is the brittleness signature for evolvability.

## 2. Study 2: Stringmol

- **What it is.** Chemart's exact port of the C++ original. It is well-mixed only, so its populations are replicator ecologies (level 0), not L1 organisations. Its copy errors are inherited.
- **Oracle.**
  - `Machine.react(a, b, "possible", budget)` with copy rates set to 0, and `bind_probability(a, b)`.
  - The `$` and `?` instructions are stochastic. "Same behaviour" is therefore judged against a **parent-vs-parent null**, with an adaptive number of samples K.
  - Products are compared as multisets.
- **Scans.**
  - Molecules: `SEED_REPLICASE`, `CONFIG_REPLICASE`, `ALIFE12_SPECIES_9`, SP29–31, and members of evolved populations.
  - Operators: uniform substitution (primary) and loop-neighbour substitution (native, reported separately). Binding-neutrality and execution-neutrality are reported separately.
  - Partners: the co-evolved population, other populations, and a frozen reference set.
  - Also computed: the junk/degeneracy split (regions never aligned or executed), L, and Clark et al.'s sets.
  - Scan chosen molecules, not whole populations: one 64-mer costs about 2 CPU-hours.
- **Dynamics.**
  - The container reactor with `cell_radius` ≈ 100, sweeping `substitution_rate` with indels fixed.
  - Quench-recovery as in AlChemy, with extinction as the absorbing collapse. Time to extinction is also recorded.
  - Knock-ins as in 1.3E.
  - **Benchmark on day 1 of Phase 2.** If a sweep would take more than about a day, run dynamics on the upstream C++ (`franticspider/stringmol`, GPL-3), cross-checked against the port.
- **Predictions.**
  - **P2.1** Uniform substitutions are more often degenerate (used but tolerant) than λ M1.
  - **P2.2 (competing)** Degeneracy either (a) raises μ₅₀ relative to an inert arm, or (b) lowers survival by letting binding-competent parasites arise. Clark et al. found both effects.

## 3. Study 3: real chemistry, MCM and Rhea

The question: *what does changing one atom do to a molecule's chemistry?*

### 3.1 Framing

- **MCM is built from structure–activity rules** (Saunders et al. 2003). So locality measured in MCM is **the chemists' model of locality**, largely true by construction. It is validated against *measured* rate coefficients for one-edit pairs from McGillen et al. 2020 (ESSD 12:1203–1216): 2,765 recommended coefficients, covering 1,357 compounds with OH, 709 with Cl, 389 with NO₃ and 310 with O₃.
- **Rhea is curated physiological biochemistry.** Enzyme promiscuity is under-represented (Khersonsky & Tawfik 2010), so neutrality measured there is a lower bound.
- **Hypotheses:**
  - Reactivity is local, carried by functional groups.
  - Edits away from reactive centres keep *which* reactions happen, and are carried into the products.
  - Rates shift gradually (linear free-energy relationships; Hammett 1937).
  - Edits at reactive centres, or in small molecules, change behaviour qualitatively. These are the chemical analogue of "activity cliffs" (Maggiora 2006).
  - Enzymatic chemistry may be *less* neutral at the level of the metabolite. Enzymes recognise the whole molecule, stereochemistry included, while the degeneracy sits in enzyme sequence space.

### 3.2 Data and cleaning

**MCM v3.3.1** (`mcm_3-3-1.tar.gz` from the MCM archive, University of York):

- Files: `mcm_3-3-1_kpp_complete.eqn` for the reactions, and `mcm_3-3-1_species_complete.tsv` for SMILES and InChI.
- Size: 17,224 reactions, 5,832 species.
- Partners: OH, O₃, NO₃, NO, NO₂ and HO₂; the RO₂ pool, recognised from `RO2` in the rate expression; photolysis.
- Reaction class: the partner plus the change in functional groups, computed with RDKit SMARTS.

**Rhea release 142** (`ftp.expasy.org/databases/rhea/tsv/`):

- Files:
  - `rhea-reaction-smiles.tsv`, keeping one left-to-right reaction per master (via `rhea-directions.tsv`);
  - `rhea-chebi-smiles.tsv`;
  - `rhea2ec.tsv` and `rhea2uniprot_sprot.tsv`, joined on `MASTER_ID`.
- Size: 18,611 reactions, 15,203 compounds.
- Cleaning:
  - neutralise charges and canonicalise tautomers before enumerating edits;
  - drop currency partners (H₂O, H⁺, NAD(P), ATP, …) from behaviour profiles;
  - drop reactions that interconvert the pair itself (e.g. a methylation);
  - cluster enzymes by family.
- **Generic `*` reactions** ("any R") are kept and analysed as curated statements of neutrality.

Loaders live in `realchem/` and download into `data/`, which is git-ignored. Chemart's `toychem` helpers `canonical` and `_fragments` can be reused.

### 3.3 Analyses

1. **Gate C.** Count the one-edit pairs in each network before committing to anything else.
2. **Per pair:**
   - the Jaccard index of (partner, class) profiles;
   - shared enzymes (Rhea);
   - whether the edit is carried through into the products;
   - rate ratios (MCM, with McGillen's measured values).
3. **L = P(same | one-edit pair) / P(same | size- and annotation-matched random pair)**, which is valid on observed networks. Stratify by edit type, heavy-atom count, and the distance from the edit to the reactive site.
4. **Matched molecular pairs** (Hussain & Rea 2010) as a robustness check.

**Predictions.**

- **P3.1** Neutrality rises with molecule size and with distance from the reactive centre. This is tested on McGillen's measured rates, not only on MCM.
- **P3.2** *(descriptive)* Rhea compared with MCM at matched size.

## 4. Synthesis

- **One molecule-level table across the four chemistries:** L, r against size, the junk/degeneracy split, ω, D and lethality.
- **Organisation-level results within each chemistry:** for AlChemy, μ₅₀ against the arms, P1.5 and evolvability; for Stringmol, the same where they apply.
- **Deciding what comes next:** neural molecules through the same pipeline; spatial or compartment designs; or level-2 merger assays using the "other organisations" partner set.

## 5. Order of work

- **Phase 0 (done 2026-09-23):** this plan and `README.md`.
- **Phase 1 milestone (about 2 weeks; delivers the first result):**
  1. Vendor, pin, build and benchmark the engine; cross-check it against chemart; unit tests (Gate A).
  2. Build the organisation library, split into pilot and confirmatory sets.
  3. Run experiment C on the pilot set: a first P_collapse(μ), with the three arms.
  4. Run scan-lite on the pilot set, for r̄ and the junk/degeneracy split.
  5. Register P1.1–P1.6 in `PREREGISTRATION.md`.
  6. Run the confirmatory experiment C, then knock-ins on 20 organisations.
- **Phase 2:** Stringmol. Benchmark on day 1, then scans, then dynamics.
- **Phase 3:** real chemistry. This can run in parallel from the start, beginning with the loaders and Gate C.
- **Phase 4:** a synthesis note in `notes/`, then decide what comes next.

## 6. Compute

This machine has 12 cores and 15 GB of RAM. The estimates assume about 10⁴ collisions per second per core; the Gate A benchmark will check that.

| task | size | cost |
|---|---|---|
| organisation library | 150 runs × 10⁶ collisions | about 4 CPU-hours |
| scan-lite | about 10⁸ collisions | about 3 CPU-hours |
| experiment C | about 10⁹ collisions | about 25 CPU-hours |
| knock-ins | similar to experiment C | similar |
| Stringmol substitution scan (Python) | one 64-mer | about 2 CPU-hours |
| RDKit edit enumeration | 10⁶–10⁷ edits | minutes |

All the AlChemy work fits in under two days of wall-clock time, even if the engine is ten times slower than assumed. Neutral random walks are deferred.

## 7. Statistics and reproducibility

- **The unit of replication is the organisation.** Organisations that recur across seeds are merged. Molecule-level data get crossed random effects (organisation × species).
- **Fitting μ:**
  - model μ = 0 as a lower asymptote;
  - fit a logistic in log μ, or a frailty Cox model for time to collapse;
  - add an adaptive second stage of μ values near each transition.
- **Pilot and confirmatory sets.** The pilot (20%) sets the grid; the confirmatory set tests the registered predictions.
  - Each prediction has one primary operator, one behaviour definition and one criterion.
  - The primary tests get a Holm correction.
  - Confidence intervals come from a cluster bootstrap.
- **Sensitivity checks:** quench length, T, M, the probe sets, the reduction cutoff, `discard_identity`, and K.
- **Records.** Every run writes a manifest: engine commit, config, seed, command. Raw runs go to `runs/` (git-ignored) and summary tables to `results/`.

## 8. Layout

```
science/
  README.md  PLAN.md  PREREGISTRATION.md (per phase)
  common/      Python: phenotypes, L, stats, IO
  alchemy/     engine/ (vendored, pinned), mut/ (our crate), *.py
  stringmol/   *.py (chemart oracle); C++ fallback if needed
  realchem/    mcm/, rhea/, mcgillen/
  notes/       per-phase results notes
  data/ runs/  git-ignored
```

## 9. Risks

- **The Rust engine may disagree with chemart or with the original.** Localise the differences with the fork and PyAlChemy before trusting any result.
- **Level-1 organisations drift on their own** (Mathis et al.). This is why collapse is judged by quench-recovery against a baseline.
- **The real networks may have few one-edit pairs** (Gate C). The fallback is scans with RDKit reaction templates as the oracle. Templates have locality built in, and that has to be said.
- **Well-mixed Stringmol goes extinct.** Time to extinction is then the measure.

## 10. Optional extensions (not scheduled)

- **BFF at level 0.** Chemart has it, the hypothesis names it, and its fixed 64-byte tapes are the cleanest case for the Greenbury test.
- Neutral random walks.
- Neural molecules through the same pipeline.
- Spatial or compartment versions.

## Open questions to settle while iterating

- **The exact non-triviality and stationarity rules** for the organisation library: minimum number of species, window, and similarity threshold.
- **How the inert placeholder is implemented.** An engine-level token excluded from collisions is clean. A closed term that reacts elastically with everything may not exist.
- **The quench length and the self-maintenance window** used in quench-recovery.
- **The size of the reference probe set, and how it is drawn.**
- **Whether to add BFF at level 0** (section 10).

## Key references

The full, verified reference lists are in `../research/sources/`.

- Mathis, Patel, Weimer & Forrest (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos* 34:093142.
- Vimal, Mathis, Weimer & Forrest (2025). Prebiotic functional programs. arXiv:2509.03534.
- Fontana & Buss (1994). "The arrival of the fittest". *Bull. Math. Biol.* 56:1–64.
- Edelman & Gally (2001). Degeneracy and complexity in biological systems. *PNAS* 98:13763.
- Greenbury, Schaper, Ahnert & Louis (2016). *PLoS Comp. Biol.* 12:e1004773.
- Clark et al. (2011). Degeneracy enriches artificial chemistry binding systems. *ECAL 2011*.
- Saunders, Jenkin, Derwent & Pilling (2003). MCM v3 protocol, Part A. *ACP* 3:161–180.
- McGillen et al. (2020). Database for the kinetics of the gas-phase atmospheric reactions of organic compounds. *ESSD* 12:1203–1216.
- Khersonsky & Tawfik (2010). Enzyme promiscuity. *Annu. Rev. Biochem.* 79:471–505.
- Hussain & Rea (2010). Matched molecular pairs. *J. Chem. Inf. Model.* 50:339–348.
- Maggiora (2006). On outliers and activity cliffs. *J. Chem. Inf. Model.* 46:1535.
- Hammett (1937). The effect of structure upon the reactions of organic compounds. *J. Am. Chem. Soc.* 59:96–103.

## Changelog

- **2026-09-23:** First version.
  - Scoped with Marco: MCM and Rhea; self-contained code; first result is collapse against mutation rate.
  - Revised after an independent review of the draft. Added: quench-recovery collapse; knock-in and knockout arms; evolvability outcomes; the locality ratio L with exhaustive skeleton enumeration; the junk/degeneracy split; restated predictions; pilot/confirmatory statistics.
  - The Rust engine is the upstream `functional-supercollider`, the one cited by the group's 2025 paper.
