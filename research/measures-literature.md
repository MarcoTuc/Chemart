# Measuring complexity, organisation and hierarchy: a literature summary

*Written 2026-09-23 for Chemart. It condenses five literature reports (in
`sources/`) and the relevant chapters of Banzhaf & Yamamoto (2015). Every
reference cited here was checked in the report named next to it (A–E), which
gives the full citation and says how it was checked. Points marked
**[synthesis]** are this summary's own reasoning, not claims from the
literature.*

## Why this exists

Chemart has two purposes. One is a morphospace of artificial chemistries. The
other is to check the premise of the redundancy hypothesis behind Neural
Molecules. That premise has two parts:

- Algorithmic chemistries such as AlChemy and BFF stop at a ceiling of
  organisation because their molecules are brittle.
- Molecules whose redundancy can be tuned, such as small neural networks,
  could lift that ceiling.

The intuition is a ladder of three levels: **molecules → networks of molecules
→ networks of networks.** In the literature these are Fontana & Buss's level 0,
level 1 and level 2 (1994). Strictly, their level 0 is "self-copying objects or
simple ensembles of copying objects", so single molecules sit just below it.

The ceiling is the step from level 1 to level 2. Mathis et al. (2024) re-ran
AlChemy. They found level-1 organisations common and robust, "but that these
stable organizations cannot be easily combined into higher order entities."

To test any of this you have to measure it, and measurement is where the
vocabulary gets slippery. "Complex", "organised", "hierarchical" and
"redundant" each have several technical meanings. The measures built on them
disagree, sometimes in opposite directions. This summary sorts them out level
by level. It ends with what Chemart has, what it lacks, and what all this means
for the hypothesis.

## What is in this folder

| path | contents |
|---|---|
| `measures-literature.md` | this summary |
| `sources/A_complexity.md` | what "complex" means: the families of complexity measures (41 references) |
| `sources/B_level0_redundancy.md` | level 0: redundancy, degeneracy, neutrality, robustness, evolvability (55) |
| `sources/C_level1_organisation.md` | level 1: what an organisation is, what makes it a unit, the ceiling in the primary sources (43) |
| `sources/D_level2_hierarchy.md` | level 2: the meanings of hierarchy, higher-order entities, how to detect them (72) |
| `sources/E_open_endedness_anchors.md` | open-endedness measures, real-chemistry reference points, morphospace method (49) |
| `library/` | 51 of the PDFs collected during the research, indexed in `library/README.md`. PDFs are git-ignored. |

Research agents wrote the reports. They were required to check every
reference online (metadata plus the abstract or full text), to quote
definitions from the sources, and to label their own inferences. Read quotes
as checked facts and inferences as proposals.

## The short answer

1. **"Complex" means at least five different things** (A).
   - *Randomness* measures are highest for noise, so they rank a random soup as
     the most complex. Examples: Shannon entropy, compression, Kolmogorov
     complexity.
   - *Structure* measures vanish for both order and noise. Examples: excess
     entropy, statistical complexity, effective complexity, logical depth.
   - Others measure *scale and nesting*, *construction history*, or *function*.

   Chemart's `compressibility` and `shannon` are randomness measures, so they
   cannot test the ceiling.
2. **A level-1 "organisation" has a precise meaning** (C). It is a set of
   molecules that is closed (makes nothing outside itself) and self-maintaining
   (regenerates every member). RAF theory's version is catalytic closure from a
   food set. These are *maintenance* criteria. Being maintained is not the same
   as being a *unit* that can become part of something bigger.
3. **"Hierarchy" has four meanings, and "networks of networks" is one specific
   combination of them** (D). It means parts nested in wholes *and*
   individuated. Flow, trophic and treeness measures see who drives whom.
   Nested-module methods see graph structure. None of them sees a
   self-maintaining whole made of self-maintaining parts.
4. **In the primary sources, the ceiling is mostly about missing units, not
   missing networks** (C, D).
   - Fontana & Buss's level-1 organisations "are self-maintaining, but not
     reproducing", and "in no sense can one identify multiple instances of the
     same L1 organization in our flow reactor".
   - Their level 2 was built by the experimenter, who merged separately grown
     organisations.
   - Every robust higher level in the artificial-chemistry literature had
     space, compartments or group reproduction built in. None was well mixed.
     Banzhaf & Yamamoto say the same (§15.3).
5. **The redundancy premise needs restating** (B).
   - λ-calculus is highly *redundant*: infinitely many terms share each normal
     form, and Fontana chose it for exactly that reason.
   - What it probably lacks is *locality*. Terms that behave alike are not one
     mutation apart. Nobody has measured this.
   - The literature supports "degeneracy plus connected neutrality →
     evolvability". It does not support "more redundancy → more complexity".
6. **Redundancy has a second, deeper role: it is what makes levels** (book
   §15.3–15.4; C; A).
   - A higher-level entity is an equivalence class of lower-level states, and
     neutral, many-to-one mappings are what create those classes.
   - Fontana & Buss make this a necessary condition for organisation: the
     number of distinct products must grow more slowly than the number of
     collision histories.
   - **[synthesis]** This is the strongest theoretical link between the
     hypothesis and the three-level ladder.
7. **Measures respond to redundancy in opposite directions** (A, B, D).
   - Physical complexity and functional information *fall* as redundancy
     rises.
   - BFF's high-order entropy and assembly theory *reward* copies.
   - Redundancy pushes Rosas et al.'s emergence criteria Ψ and Δ *negative*.
   - Causal emergence measured by effective information *rises* with
     degeneracy.

   Know each measure's bias before using it in a redundancy experiment.
8. **The strict organisation count is the wrong detector for a finite soup**
   (C). In a constructive chemistry an organisation is infinite and only partly
   present at any moment. So a snapshot is "no longer closed under interaction"
   even while the organisation persists. That predicts Neural Molecules'
   0.00–0.02 organisations per run. Use windowed, flux-weighted closure,
   persistence and perturbation tests instead.
9. **No standard open-endedness measure detects a new level on its own** (E).
   They all count types at a level chosen in advance. A new level has to be
   recognised first. Then the same statistics are rerun one level up, against
   null models.
10. **The cheapest valuable additions to Chemart** are in
    [section 8](#8-what-chemart-has-and-what-to-add):
    - measure molecular locality for λ, BFF and neural molecules;
    - measure convergence;
    - compute redundancy and degeneracy of the interaction table;
    - generate level-2 candidates (organisations joined by "glue") from the
      organisations and RAFs Chemart already computes;
    - add a merge-and-knockout assay.

---

## 1. What "complex" means

### Five families of measures

Lloyd (2001) sorted about forty measures by three questions: how hard is it to
describe, how hard is it to create, and how organised is it? Report A splits
them into five families:

| family | the question | examples | scores highest for |
|---|---|---|---|
| randomness | how much information does it take to specify this exact object? | Shannon entropy, Kolmogorov complexity, compression (zlib, brotli) | pure noise |
| structure | how much pattern is there beyond noise? | excess entropy, statistical complexity, effective complexity | mixtures of order and randomness (zero at both ends) |
| scale and hierarchy | at how many scales, or nested levels, is there structure? | TSE complexity, the complexity profile, McShea's nesting levels, Simon's near-decomposability | parts that are both distinct and integrated |
| history | how much work did it take to make? | logical depth, thermodynamic depth, assembly index | objects with a long, non-trivial production history |
| function | how much information is *about* something? | physical complexity, functional information | rare configurations that work |

### Why the split matters for artificial chemistries

- **Random starts.** AlChemy and BFF start from random soups, so randomness
  measures *fall* as organisation appears. Gell-Mann & Lloyd call Kolmogorov
  complexity "algorithmic randomness" for this reason.
- **Copies.** Families disagree about copies. Shannon diversity penalises them.
  BFF's high-order entropy and the assembly measure A reward them. Thermodynamic
  depth ignores them by design: "a complex object together with a copy is not
  much deeper than the object alone" (Lloyd & Pagels 1988). So a takeover by
  one replicator looks trivial, complex or unchanged depending on the family
  you chose.
- **BFF's measure detects replication.** Its "high-order entropy" is the
  Shannon entropy per byte minus the compressed size per byte. Read closely, it
  measures verbatim repetition, what Bennett (1988) calls "obvious
  redundancy". It detects replication (level 0), not organisation (A).
- **Vanishing at both ends is not enough.** Feldman & Crutchfield (1998) showed
  that hump-shaped measures can reduce to trivial functions of entropy. A
  measure must say what structure it counts.
- **Size is not complexity.** Mitchell's amoeba has about 225 times the base
  pairs of a human. In Tierra, organism size grew while complexity did not
  (Standish 2003). Chemart's `complexity_drift` measures size.
- **AC molecules are too short to compress.** zlib's fixed overhead dominates:
  `SKI` compresses to 11 bytes. Per-molecule compression ratios therefore
  measure length (A).

### What the hypothesis needs from a complexity measure

The hypothesis makes three different claims, and each needs its own
measurement (A §1.3).

- **Level 0, "algorithmic molecules are brittle":** this is a claim about the
  map from a molecule's structure to its behaviour, not about complexity.
  Measure robustness and degeneracy directly (section 2). Expect function-based
  complexity to *drop* as redundancy rises, because neutral sites carry no
  information. That drop is predicted, not a refutation.
- **The ceiling:** the measure must be near zero for a random soup, stay low
  for a takeover by one replicator, and grow with the number of
  differentiated, interdependent parts and nested levels. Candidates:
  - McShea's hierarchy level, taking its maximum over a run;
  - TSE complexity or the complexity profile, computed on abundance dynamics;
  - the size of an organisation's smallest self-maintaining generator set
    (Fontana & Buss's "centre"), as a proxy for effective complexity;
  - the functional information of reaching level *k*, i.e. the fraction of
    random soups that get there.
- **Level 2 needs integration, not coexistence:**
  - The complexity profile of two independent subsystems is the sum of their
    profiles (Allen, Stacey & Bar-Yam 2017). Any excess over the sum is
    integration.
  - TSE complexity is zero across independent parts.
  - No study of an artificial chemistry has measured either.

One rule follows: **measure cause and effect independently** (McShea 1996). If
redundancy and organisational complexity are measured by related operations,
the test is circular.

---

## 2. Level 0: molecules

### The words

| term | meaning | measured on |
|---|---|---|
| redundancy | many structures map to one behaviour; in biology also "the same function performed by identical elements" (Tononi, Sporns & Edelman 1999) | the map |
| degeneracy | "the ability of elements that are structurally different to perform the same function or yield the same output" (Edelman & Gally 2001). Functions overlap only partly: alike in some contexts, different in others (Whitacre 2010) | parts × contexts |
| neutrality | a mutation that leaves the behaviour unchanged. It is always relative to a definition of behaviour, and to the context (Wagner 2005) | molecule + mutation operator + behaviour definition |
| robustness | the share of mutations that are neutral: per molecule (genotype robustness r_g) or averaged over everything with one behaviour (phenotype robustness ρ_p) (Wagner 2008) | molecule, or the map |
| evolvability | the number of *new* behaviours one mutation away, from one molecule (e_g) or from a whole neutral set (ε_p) (Wagner 2008) | molecule, or the map |
| locality | "neighbouring genotypes correspond to neighbouring phenotypes" (Galván-López et al. 2011, after Rothlauf 2006) | the map |

### What the literature found

- **Robustness and evolvability work against each other for one molecule and
  together for a behaviour** (Wagner 2008; Ahnert 2017). A robust behaviour
  has a large, connected neutral network, and that network borders many other
  behaviours.
- **Redundancy is not robustness.** Redundant structures also have to be
  mutational neighbours. The diagnostic compares phenotype robustness ρ_p with
  phenotype frequency f_p. A randomly scrambled map gives ρ_p = f_p. RNA
  folding and other biological maps give ρ_p ≫ f_p (Greenbury et al. 2016).
- **Brittleness is an old problem in artificial life.**
  - Ray (1991): machine code is "brittle, meaning that the ratio of viable
    programs to possible programs is virtually zero".
  - Ofria et al. (2002) measured "over 99.7% of all nontrivial mutations" as
    deleterious for redcode.
  - The known fixes are to remove numeric operands, to address by template
    (Tierra, Avida), and to bind approximately. Stringmol's designers wrote
    that "we would like there to be many ways of achieving a bind probability
    of 1."
- **The earliest statement of the hypothesis found** is Conrad (1990):
  "Biological structures that are characterized by a high degree of component
  redundancy and multiple weak interactions satisfy these conflicting
  pressures." Banzhaf & Yamamoto (§8.3.3) state the open problem: the
  chemical language "must produce viable and fertile individuals with high
  probability. This is an area where there is currently no firm recipe for
  success."
- **The only direct evidence from an artificial chemistry:** Clark et al.
  (2011) removed binding degeneracy from Stringmol and lost "emergent
  macro-mutations, hypercycles, sweeps and parasite evasion".

### λ-calculus: redundant, but probably not local

- **Evaluation is many-to-one.** `M`, `(λx.x) M` and `(λx.λy.x) M N` (for any
  `N`) all reduce to the same normal form. With Church numerals, `PLUS 2 3` and
  `PLUS 1 4` both give 5. Fontana & Buss built on this deliberately:
  "chemistry's diversity of equivalence classes, that many different reactants
  can yield the same stable product."
- **Where the redundancy lives in AlChemy.** Chemart's AlChemy, like the
  original, stores each molecule as a canonical closed normal form in de Bruijn
  notation. So:
  - the redundancy lives in the *reaction*, where many collisions give one
    product;
  - it does not live in the *molecule*. Distinct normal forms can almost always
    be told apart by some context (Böhm's theorem, 1968, a standard result not
    checked in the reports).
- **Why it can still be brittle.** Terms that behave alike are scattered, not
  neighbours. `λx.λy.x` returns its first argument. Change one variable and
  `λx.λy.y` returns the second. Delete one symbol from `λx.(x)x`
  (self-application) and you get `λx.x` (identity).
- **Untested.** Nobody has measured λ-term robustness, and no paper argues that
  AlChemy is limited *because* its molecules are brittle (B §3.6). The only
  published data point runs against the naive version: real software is not
  that brittle, with "over 30% of random mutations" neutral against a test
  suite (Schulte et al. 2014).

**[synthesis]** The testable premise is therefore: *λ-terms and BFF tapes have
low locality (low r_g, ρ_p ≈ f_p, low viability density); neural molecules have
high locality.*

### Neural-network molecules

- **Locality comes free.** Behaviour changes smoothly with the weights, so near
  neighbours behave alike.
- **Exact redundancy is tiny.** For an irreducible network with 3 hidden tanh
  units, the only weight settings with exactly the same function are the 48
  permutations and sign flips of hidden units (Sussmann 1992; Farrugia-Roberts
  2023). Neutral networks appear only once behaviour is coarse-grained into
  species. So neutrality is a joint property of the substrate and the
  observer's resolution. Two consequences:
  - canonicalise weights before counting weight species, as de Bruijn indices
    do for λ-terms;
  - report neutrality as a curve over species resolution, not as one number.
- **Weight scale also moves simplicity bias** (Mingard et al. 2025, deep tanh
  networks; it has yet to be measured for 3→3→3). A sweep that seems to vary
  redundancy may actually vary how simple the common behaviours are.
- **Capacity limits complexity too.** In Geb, the maximum complexity of
  neural-network agents grew without bound only when both network capacity and
  population size were scaled together, and then only logarithmically
  (Channon 2019; E). Fixed 24-parameter molecules have a capacity ceiling of
  their own.

### The tension: does redundancy raise or lower complexity?

**Evidence that it lowers complexity:**

- Neutral sites carry zero information (physical complexity, Adami et al. 2000).
- Many-to-one maps are biased toward simple outputs (Dingle et al. 2018), and
  so are neural-network parameter maps (Valle-Pérez et al. 2019).
- Frequent behaviours arrive first and fix, even when fitter ones exist (the
  "arrival of the frequent", Schaper & Louis 2014).
- Selection for robustness can simplify circuits (Milano et al. 2019) or lock
  populations in (Ancel & Fontana 2000).

**Evidence that degeneracy raises evolvability, and with it complexity:**

- Selecting for degeneracy yields high complexity (Tononi et al. 1999).
- "Purely redundant systems have remarkably low evolvability while degenerate,
  i.e. partially redundant, systems tend to be orders of magnitude more
  evolvable" (Whitacre & Bender 2010).
- Stringmol's degeneracy result (Clark et al. 2011, above).

**What the literature concludes** (B §5): no source supports "more redundancy →
more complexity". What is supported is that *degenerate, mutationally connected
neutrality at the level of behaviour* raises evolvability, and through it the
potential for complexity.

**Two cautions for soups:**

- Populations drift into robust regions without any selection for robustness.
  In RNA, the excess is about 30% (van Nimwegen et al. 1999). So robustness
  rising during a run needs a random-walk baseline.
- Large populations may evolve *antiredundancy* (Krakauer & Plotkin 2002).
  Brittleness can be selected, not only inherited from the substrate.

### Level-0 measures

All of these except the last two need the chemistry's reaction rule as
something you can call, to evaluate mutants. Behaviour is best defined
identically for every chemistry: **the products a molecule makes with each
member of a fixed probe set of partners** (B §2). That definition is what
makes λ-terms, BFF tapes and neural networks comparable.

| measure | what it tells you | source |
|---|---|---|
| genotype robustness r_g; P_same(k) for k mutations | the brittleness claim itself | Wagner 2008; Lenski et al. 1999 |
| viability density | the share of random molecules that do anything (normalise, replicate): the cheapest single brittleness number | Ray 1991; Suzuki 2003 |
| genotype evolvability e_g | whether mutants are *new* rather than dead | Wagner 2008 |
| ρ_p and ε_p from neutral random walks | whether redundancy is connected, and whether it borders many behaviours | Wagner 2008; Ahnert 2017 |
| ρ_p against f_p | redundant-but-scattered (ρ_p ≈ f_p) versus redundant-and-connected (ρ_p ≫ f_p). It is comparable across substrates because each uses its own point mutation | Greenbury et al. 2016 |
| rank–frequency of f_p, plus the complexity of the frequent behaviours | whether redundancy only buys simplicity | Dingle et al. 2018 |
| excess robustness over the network mean | robustness beyond what drift alone gives | van Nimwegen et al. 1999 |
| redundancy and degeneracy sets of the interaction table | *computable from Chemart's existing reaction networks* | Clark et al. 2011 |
| physical complexity ℓ − Σ H_i | information per molecule; expected to fall with redundancy | Adami et al. 2000 |

---

## 3. Level 1: networks

### What "organisation" means

| concept | criterion | captures | misses |
|---|---|---|---|
| Fontana & Buss organisation | "a kinetically self-maintaining algebraic structure", with a boundary set by invariances, self-repair, and a centre | infinite organisations only partly present; self-repair | reproduction; no algorithm |
| chemical organisation (COT) | closed + self-maintaining (a feasible flux regenerates every member) | the lattice of organisations; fixed points lie inside organisations | whether the flux is realised; small numbers; units |
| RAF | every reaction's reactants and one catalyst are made from the food set by the set itself | collective catalysis; subRAF structure; polynomial-time algorithms | stoichiometry; closure (products leak); the union of RAFs is always a RAF, so it says nothing about integration |
| autocatalytic cores | minimal stoichiometric self-amplifying motifs | growth, not just persistence | whether autocatalysis actually amplifies; units |
| autopoiesis; closure of constraints | self-produced boundary; mutually dependent constraints | individuation; levels as "separated and hierarchically nested" closures | reproduction; no general algorithm |
| chemoton; dynamic kinetic stability | reproduction of the whole; persistence through reproduction | a unit of reproduction | a design, or a stability notion, not a detector |

A closed RAF is always a COT organisation, but not the other way round
(Hordijk, Steel & Dittrich 2018).

### Maintained is not a unit

Report C collects, from the literature, what a level-1 entity needs before it
can act as a *part* at level 2:

- **Maintenance:** closure; kinetic persistence under flow; self-repair.
- **A boundary** that decides membership.
- **Individuation:** several separable instances.
- **Reproduction of the whole.**
- **Heredity with more than a few states:** "systems with limited heredity have
  only a limited evolutionary potential" (Szathmáry 2000).
- **Selectability:** evolvability appears only with many viable cores inside
  compartments that divide (Vasas et al. 2012).
- **Composability through glue.**

The published network measures, including all of Chemart's, test only
maintenance.

### Two strengths of level 2

- **Weak level 2: an ecological composite.** Two organisations persist inside
  a larger one, joined by glue that "is not self-maintaining, but which acts to
  knit the self-maintaining Level 1 sets into a higher-order self-maintaining
  entity" (Fontana & Buss). The parts keep their own self-maintenance.
- **Strong level 2: a transition in individuality.** A new unit of
  reproduction appears (Szathmáry 2015), and the parts lose autonomy (Goldsby
  et al. 2012). No AlChemy result reaches this.

Keep the two apart. They need different measures (section 4).

### Why a finite soup has no strict organisations

These points are from C §5.4.

- **Kinetic confinement.** Only a changing finite core of an infinite
  organisation is present. New products appear all the time and are "diluted
  out" (Mathis et al.).
- **Turnover.** "The expressions are continually changing" while the
  organisation persists. The invariant is a grammar or an equivalence class,
  not a fixed set of species.
- **Finite windows.** A network built from observed events holds only the
  reactions that happened to fire, so rare products break closure.
- **Small numbers.** The theorems are about ODEs, and single molecules need a
  discrete theory (Kreyssig et al. 2014).
- **Identity in continuous chemistries.** Under exact identity nearly every
  reaction makes a new species, so closure is impossible by construction. The
  count depends on the species tolerance, so report it as a curve.
- **Recursively growing structures** keep their function without any member
  recurring exactly (Kruszewski & Mikolov 2021). Any species-set criterion
  misses them.
- **Trivial organisations** inflate the count the other way: single
  self-replicators, and the empty set when decay is modelled.

### Conditions Fontana & Buss found necessary

- **Convergence.** "for organization to occur it is necessary that the total
  number of normal forms … grows with lower order than the total number of
  collision sequences." If every collision makes something new, "no network can
  ever be formed." Convergence is measurable, and it is the first place tunable
  redundancy should show an effect.
- **Copying pre-empts level 1.** If identity functions arise early, "the system
  typically does not reach Level 1." BFF is a world dominated by copying, so
  its first barrier may be copying rather than brittleness (C §4.4).
- **Food and waste conventions.** Reactants are not consumed, and inert
  products are removed. "If this assumption is relaxed, self-maintaining
  organizations fail to emerge." Kruszewski & Mikolov did get
  self-reproducing metabolisms with mass conservation, so this is a caution,
  not a verdict. It still matters for any chemistry that conserves mass,
  Neural Molecules included.

### Level-1 measures to add

| measure | definition | source |
|---|---|---|
| M1 closure in a window | of the reaction events in window W whose reactants are all in S, the fraction whose products are all in S too | B&Y §13.4, weighted by occurrence |
| M2 self-maintenance in a window | the fraction of members whose production in W is at least their loss. B&Y's printed α_S turns out to equal α_C, so it is replaced here | COT flux condition on realised fluxes |
| M3 persistence | Jaccard similarity of the set to itself after a lag, plus an abundance-weighted version and a membership half-life | Mathis et al. 2024 |
| M4 perturbation | survival when the soup is flooded with identity; incorporation of injected molecules; delete-and-regrow from a seed | Fontana & Buss; Mathis et al. |
| M5 convergence | distinct products against reaction events (over species tolerance, for continuous chemistries) | Fontana & Buss |
| M6 λ₁ per window | spectral radius of the observed catalytic graph | Jain & Krishna 2002 |
| M7 merger assay | merge two evolved organisations; classify dominance, coexistence, mutual destruction or fusion with glue; measure the glue fraction; isolation test | Fontana & Buss; Mathis et al. |
| M8 decomposition | an organisation that is the union of two smaller ones plus glue; a RAF that is the union of two subRAFs | COT lattice; Hordijk, Steel & Kauffman 2012 |

---

## 4. Level 2: networks of networks

### Four meanings of "hierarchy"

| sense | the question | measures | can it see level 2? |
|---|---|---|---|
| order, control, flow | who drives whom? | flow hierarchy (Chemart's `flow_hierarchy`), global reaching centrality, trophic coherence, treeness/feedforwardness/orderability | no |
| compositional nesting | what is part of what? | nested modules (Ravasz; hierarchical random graphs; the nested SBM), the COT lattice, P-system membranes | half: nesting without individuation, except where the parts are self-maintaining |
| individuality, levels of selection | which wholes are units that selection acts on? | multilevel Price partition, Okasha's stages, Michod's covariance | yes, if collectives reproduce |
| specification | which realm refines which? | none | no |

Beware one search trap: "networks of networks" in network science means
*interdependent* networks, which fail and percolate together (Gao et al.
2012). That meaning is unrelated.

Flow-type hierarchy measures actively penalise what level 1 is made of.
Treeness, feedforwardness and orderability collapse every strongly connected
component into a single node, which means every autocatalytic core, and they
score cycles as violations of hierarchy (D §6.1). They make good morphospace
axes, but they are not level-2 detectors.

**Working definition** (D §1). A level-2 entity is a set of at least two
level-1 organisations that is:

- **nested:** each part is self-maintaining;
- **integrated:** there is glue that no part makes alone, and the parts depend
  on each other;
- **individuated:** closed, self-repairing, and bounded as a dynamical unit;
- **more than an aggregate:** it has slower dynamics of its own.

Reproduction as a unit is an optional fifth property.

### The ceiling, as the primary sources describe it

- **Fontana & Buss (1994).** "L2-organizations are more readily attained by
  first generating L1-organizations separately, and subsequently combining
  them into the same reactor." It helps if the glue "operate[s] under
  different boundary conditions, in particular different collision rules." The
  compact version is "Centers compose, organizations not." Banzhaf & Yamamoto
  (§9.1): spontaneous level 2 "is extremely rare".
- **Mathis et al. (2024).** They merged 455 pairs of evolved level-1 soups.
  Coexistence (Jaccard similarity above 0.1 to both inputs) "rarely occurs for
  organizations evolved in different simulations."
- **What enabled higher levels where they did occur.** In each case, space,
  compartments or imposed group reproduction:
  - spiral waves of hypercycles on a lattice (Boerlijst & Hogeweg 1991);
  - Stringmol in space (Hickinbotham, Stepney & Hogeweg 2021);
  - proto-cells in a combinator chemistry on a planar graph (Speroni di
    Fenizio, Dittrich & Banzhaf 2001);
  - self-reproducing cells in the Ono–Ikegami lattice model (2000);
  - dividing compartments (Vasas et al. 2012; Takeuchi & Hogeweg 2009);
  - colonies in Avida (Goldsby et al. 2012).

  No well-mixed chemistry robustly produced level 2 from random starts (C
  §4.7).
- **Why well-mixed reactors resist level 2.** Banzhaf & Yamamoto (§15.3):
  "Generally, this is impossible in a well-mixed reactor, only one of the
  potentially alternative groups will be dominant." Their suggested route is
  selection among many reactors (Decraene, Mitchell & McMullin 2008).
  Szathmáry (2015) explains why. Two *unlike* partners, which is what Fontana
  & Buss's level 2 joins, stabilise only through "fairness in reproduction"
  of the partners as a unit. A single well-stirred reactor has no mechanism
  for that. So the ceiling is expected there however redundant the molecules
  are (D §3.2).
- **A published relative of the hypothesis.** Rasmussen et al.'s (2001)
  Ansatz: "an appropriate increase in object complexity of the primitives is
  necessary and sufficient for generation of successively higher-order
  emergent properties through aggregation." It claims, as the hypothesis does,
  that properties of the lowest objects limit which levels can form. It was
  contested. Gross & McMullin (2001) got comparable phenomena from simpler
  primitives. Dorin & McCormack (2002) built an infinite self-assembling
  hierarchy whose new properties were trivial. The lesson: loose definitions
  of "level" are trivially satisfiable.

### Detecting level 2: a pipeline

Report D (§7) assembles this pipeline from the sources. Stages A–D need only
data Chemart already produces. Stage E needs collectives that reproduce.

1. **Stage A, level-1 units.** Find COT organisations or irreducible and
   closed subRAFs, and keep those that persist and repair themselves under
   perturbation.
2. **Stage B, candidates.** Find organisations or RAFs that are joins of at
   least two incomparable units and have non-empty *glue*: species made only
   by reactions between units.
   - RAF theory has a polynomial-time test for whether a RAF is the union of
     two proper subRAFs.
   - Score each candidate on McShea's (2001) scale: nesting depth, and
     sublevels a (copies), b (differentiated parts) and c (intermediate
     parts).
   - Keep the per-run maximum, following McShea's "trend in the maximum".
3. **Stage C, individuation.**
   - Knock out each part: does the other persist, and does the glue vanish?
     That is the loss-of-autonomy test.
   - Map which perturbations the whole repairs and which destroy it, compared
     with its parts (Beer 2014).
   - Inject parasites and measure whether they invade.
   - Check that the parts stay distinct and don't dissolve into one
     organisation (McShea's "overconnectedness").
4. **Stage D, level-ness.** Let Z be the abundances of the parts and the
   glue.
   - Informational closure: I(X̃ₜ; Zₜ₊₁ | Zₜ) ≈ 0, and non-trivial
     (Bertschinger et al. 2006; Pfante et al. 2014; Rosas et al. 2024).
   - Krakauer et al.'s (2020) individuality measures, against size-matched
     random sets, because those measures never decrease with size.
   - Timescale separation (Simon, Salthe, Flack).
5. **Stage E, evolutionary individuality.** This needs reproducing
   collectives: compartments, space, or many reactors with dispersal
   ("ecological scaffolding", Black, Bourrat & Rainey 2020). Then measure the
   Price partition between and within collectives, Okasha's stages of fitness
   decoupling, and Michod's covariance.

**Acceptance rule (Wimsatt).** Call something level 2 only when the
structural, interventional and informational stages agree. Report each stage
as its own morphospace axis rather than folding them into one score.

### Pitfalls at level 2

- **Graph hierarchy is not an entity** (D §6).
- **"A new observable at the higher level" is too weak a criterion** (Dorin &
  McCormack 2002).
- **Redundancy biases the emergence measures in opposite directions.**
  - Rosas et al. (2020): "redundancy will drive Ψ and Δ more negative", and
    the criteria "double-count redundancy up to n times". They are biased
    against exactly the chemistries the hypothesis wants to build.
  - Effective-information causal emergence finds *more* emergence in
    degenerate networks (Klein & Hoel 2020).
- **Individuality scores never decrease with system size**, so control for
  size (Krakauer et al. 2020).
- **The number of possible coarse-grainings grows super-exponentially**
  (Rosas et al. 2024). Generate candidates from structure rather than
  searching all of them.
- **Noise and isolated systems are trivially "closed".** Require non-trivial
  closure.
- **"Collectively autocatalytic" does not mean "cooperative"** (Szathmáry
  2015).

---

## 5. Over time: novelty, innovation, emergence

Banzhaf et al. (2016) give the framing closest to the book:

- **variation** is novelty within the model;
- **innovation** is novelty that changes the model (a new type);
- **emergence** is novelty that changes the meta-model, and "one particularly
  important phenomenon is the emergence of a new level of organisation."

Only the *first* entity at a new level counts as emergence. A measure defined
inside a model sees only variation (Stepney 2021: standard measures are
"measuring 'Flatland'").

| measure | what it does | runs on Chemart trajectories? | registers a new level? |
|---|---|---|---|
| evolutionary activity, classes 1–3 (Bedau et al. 1998) | persistence-weighted counts of components against a neutral "shadow" run | yes, with a shadow built from `fired` | only if the components are level-N units |
| MODES change, novelty, ecology (Dolson et al. 2019) | counts and entropy of *persistent* types | yes, with persistence counted in turnovers | no. The authors "would welcome a measurement of a system's potential to produce major transitions" |
| MODES complexity | the most informative sites in any persistent genome | only with a molecular knockout | no, and it "may cause fragile genomes to appear more complex than robust ones" |
| QNN (Droop & Hickinbotham 2012) | squared positive changes in species proportions, with no shadow needed; used on Stringmol | yes, directly | no |
| cardinality leap (Sayama 2019) | cumulative unique higher-order replicating entities; molecules per replication event | partly | closest of the quantitative measures, but it still rose under a biased random control |
| passive against driven trends (McShea 1994) | a driven trend raises the minimum and biases ancestor→descendant changes | yes, from `fired` | no |
| hierarchy scale, trend in the maximum (McShea 2001) | the deepest nesting reached | via a recognizer | yes |

Two consequences for Chemart's existing measures:

- `novelty_rate` has no persistence filter and no shadow run.
- `complexity_drift` can rise passively, by spreading away from the minimum
  length, and it measures size rather than complexity.

Judge boundedness by fitting bounded and unbounded models, not by eye (MODES).
Remember that a trivial string system can pass every proposed open-endedness
criterion (Hintze 2019).

---

## 6. Anchors in real chemistry

- **Raw topology separates little.** Planetary atmospheres, the interstellar
  medium and metabolism are all small worlds (Solé & Munteanu 2004).
  Scale-free claims mostly fail strict tests (Broido & Clauset 2019; Smith,
  Kim & Walker 2021). Tracing carbon atoms instead of letting currency
  metabolites act as shortcuts makes *E. coli*'s "small world" not small
  (Arita 2004). Representation decides the answer.
- **Departure from size-matched null models separates more.** Earth's
  atmospheric network has "the most nonrandom topology" (Wong et al. 2023).
  Random samples of the same reaction universe do not reproduce biochemistry's
  scaling laws (Kim et al. 2019).
- **Real levels differ in degree, not kind.** Across individual → ecosystem →
  biosphere networks, topology shows "no sharp transition" (Smith et al.
  2021). The levels differ in how measures scale with size, and against random
  merges of genomes (Kim et al. 2019).
  - **[synthesis]** This is a real-data analogue of the three-level ladder. It
    suggests the same test for Chemart: compare a candidate level-2
    organisation against random unions of level-1 organisations of the same
    number and size.
- **Calibration points Chemart can recompute:**
  - maxRAF size against the food set in prokaryotic metabolism (Xavier et al.
    2020);
  - scope from prebiotic seed compounds (Goldford et al. 2017);
  - atmospheric networks of Earth, Mars and Venus (Solé & Munteanu 2004).

  Treat the assembly index with care: minerals reach 21, above the proposed
  biosignature threshold of 15 (Hazen et al. 2024).
- **Morphospace pitfalls** (E §4):
  - axes are hypothesis-laden and correlated through size;
  - raw values mostly measure size (van Wijk et al. 2010);
  - read only departures from a null cloud;
  - empty regions can be geometric, functional, historical or sampling
    artefacts;
  - distances mean little unless the axes are commensurable;
  - catalog entries are not independent points, since many are variants of
    one another.

---

## 7. What this means for the redundancy hypothesis

**[synthesis]** This section is this summary's reading of the five reports.

### The hypothesis, split into testable parts

- **H0, the premise (level 0).** Algorithmic molecules (λ-terms, BFF tapes)
  have low locality and neural molecules high.
  - *Measure:* r_g, ρ_p against f_p, and viability density, all with the same
    probe-set definition of behaviour.
  - *Status:* untested for any artificial chemistry.
- **H1 (level 1).** Higher locality or degeneracy gives more convergence, and
  more frequent and more robust level-1 organisations.
  - *Measure:* convergence (M5), windowed closure (M1, M2), persistence (M3),
    perturbation (M4).
  - *Mechanisms suggested by the sources:* Fontana & Buss's convergence
    condition; resistance to parasites and copiers.
- **H2 (level 2).** Higher locality or degeneracy lifts the ceiling from level
  1 to level 2.
  - *Measure:* coexistence and glue fraction in merger assays (M7), then the
    Stage C–D tests on the composites.
  - *Status:* the literature places this ceiling mainly in missing
    individuation and reproduction, which well-mixed reactors cannot supply.
  - **H2 therefore needs a factorial design:** redundancy × {well mixed,
    spatial or graph-based, dividing compartments or many-reactor
    scaffolding}.
  - *How to read the result:* if redundancy raises coexistence even in
    well-mixed merges, it matters on its own. If it helps only once units
    exist, the ceiling was structural and redundancy is a secondary factor.

### A deeper, recursive version

Several sources point the same way:

- Banzhaf & Yamamoto (§15.3–15.4), citing Wagner: "neutral mappings are a key
  factor in producing the equivalence classes of lower-level entities that
  lead to higher-level identity."
- Ellis (in B&Y §14.4.3): "multiple realizability of higher level functions,
  and consequent existence of equivalence classes of lower level variables."
- Fontana & Buss's convergence condition.
- Allen, Stacey & Bar-Yam: "A large-scale behavior requires redundant
  information among the many components engaged in that behavior."

Read together, they suggest the hypothesis is not only "robust molecules
evolve better". It is **"a level can form only on top of a level whose states
are degenerate with respect to what the higher level depends on"**. That
statement applies at every rung of molecules → networks → networks of
networks. Its measurable signature is the Stage D test (closure and timescale
separation of the higher-level variables), and M5 at the lowest rung.

### Confounds to control

- **Well-mixedness:** above all.
- **Copying pre-emption:** relevant to BFF.
- **Food and waste conventions:** relevant to mass-conserving chemistries.
- **Survival:** redundancy can change survival too. This is the lesson of
  Neural Molecules' P3, where weight scale moved survival and redundancy
  together.
- **Size:** Krakauer's measures grow with size, and so do most morphospace
  axes.
- **Species resolution:** report results as curves over it.
- **Trivial symmetries:** canonicalise before counting.
- **Population size and mutation rate:** van Nimwegen's excess robustness;
  Krakauer & Plotkin's antiredundancy; Draghi et al.'s conditions.
- **Circularity:** measure redundancy and complexity by unrelated operations
  (McShea).

---

## 8. What Chemart has and what to add

| level | what to measure | Chemart now | to add | data it needs | report |
|---|---|---|---|---|---|
| 0 | locality of structure → behaviour | `mean_structure_length`, `structure_function_mi` (crude) | r_g, e_g, ρ_p against f_p, P_same(k), viability density | the reaction rule as a callable, a probe set, a mutation operator per chemistry | B |
| 0→1 | degeneracy of the interaction table | `degeneracy` (routes from food; structural) | Clark et al.'s redundancy and degeneracy sets | reaction networks (exist) | B |
| 1 | organisation in dynamics | `organisations` (strict; up to 200 species), `max_raf_fraction`, `irreducible_rafs`, `autocatalytic_cores`, `catalytic_spectral_radius` | M1–M4 and M6 over windows | trajectories (exist); the ability to perturb a running soup | C |
| 1 | convergence | — | M5 | trajectories | C |
| 2 | candidates | — | COT join with non-empty glue; the subRAF union test; McShea score and its per-run maximum | organisations and RAFs (exist) | D |
| 2 | individuation | — | merger assay (M7); knockout and loss of autonomy; parasite challenge | the ability to merge soups and knock sets out | C, D |
| 2 | level-ness | — | informational closure; Krakauer's measures with a size control; timescale separation | replicate trajectories reduced to a few macro-variables | D |
| 2 | evolutionary individuality | — | Price partition; Okasha's stages | compartments, space, or many reactors with dispersal | D |
| time | novelty and growth | `novelty_rate`, `complexity_drift`, `turnover` | persistence filter plus neutral shadow; activity statistics; QNN; MODES change and ecology; passive-against-driven tests; boundedness by model fit | trajectories (exist) | E |
| complexity | structure, not randomness | `compressibility`, `shannon`, `degree_entropy` (all in the randomness family) | excess entropy on symbolised trajectories; TSE complexity or the complexity profile on abundance covariance | long or replicate trajectories | A |
| anchors | real chemistry | `measures.scaling` over sweeps | size-matched nulls from each chemistry's own reaction universe; random unions of level-1 units | imported real networks | E |

**The data gap.** A frame records only the well-mixed population (`state`),
the reactions fired, and chemistry-specific observables. There are no
positions and no compartment membership. The exceptions are chemistries with
their own reactors, such as autopoiesis-vmu and ono-ikegami-protocell, and P
systems, where membrane nesting is given by construction. Stage E, and any
individuation test with more than one instance per organisation, needs space,
compartments or many reactors.

**The doc gap.** The references on `docs/guide/measures.md` are still marked
"not yet checked". Most are now verified in reports A, D and E.

---

## 9. Decisions for Marco

1. **Which level 2?** Weak (coexistence plus glue), strong (a new unit of
   reproduction whose parts lose autonomy), or both as separate axes. Each
   needs different measures.
2. **Space, compartments or scaffolding.** Should Chemart's evolve face gain
   space, compartments or many-reactor scaffolding, and should Neural
   Molecules? Without one of them, level-2 individuality cannot even be
   measured, let alone reached.
3. **What the hypothesis is about.** Redundancy, degeneracy or locality? The
   sources suggest locality plus degeneracy, with redundancy as the
   many-to-one map that makes levels.
4. **The shared behaviour definition.** Products against a fixed probe set is
   the proposal. It decides whether level-0 measures can be compared across
   substrates.
5. **Validation against ground truth.** Should the complexity measures be
   validated against Fontana & Buss's levels, i.e. does a measure rank level
   0 < 1 < 2 across AlChemy runs? No study has done this (A §3, item 12). It
   would be a contribution in itself.

---

## 10. Where to start reading

✓ marks a PDF that is in `library/`.

1. Fontana & Buss 1994, "The arrival of the fittest". The three levels, the
   ceiling, and the convergence condition. ✓
2. Mathis, Patel, Weimer & Forrest 2024, "Return to AlChemy". The re-run and
   the merger experiment. ✓
3. Banzhaf & Yamamoto 2015, §8.3–8.6, §12.3–12.4, §13.3–13.4, §14.4–14.5 and
   §15. The book's own frame. (local copy)
4. Wagner 2008, "Robustness and evolvability: a paradox resolved"; Greenbury
   et al. 2016 (ρ_p against f_p).
5. Whitacre 2010 on degeneracy ✓, and Clark et al. 2011, degeneracy in
   Stringmol ✓.
6. Lloyd 2001 (two pages) and Mitchell 2009, ch. 7. The map of complexity
   measures. (Both PDFs were lost when report A's agent cleaned up; see
   `library/README.md`.)
7. McShea 2001, the hierarchy scale ✓.
8. Szathmáry 2015, "Toward major evolutionary transitions theory 2.0".
9. Rosas et al. 2024 ✓ and Krakauer et al. 2020 ✓ (preprint). Levels as
   closure; individuality.
10. Banzhaf et al. 2016, variation, innovation and emergence ✓; Dolson et al.
    2019, the MODES toolbox ✓.
11. Kim et al. 2019, levels and null models in real biochemistry.

---

## References

Full, verified citations with DOIs are in the report named after each entry.
Two entries, marked below, were checked for this summary instead.

- Adami, Ofria & Collier (2000). Evolution of biological complexity. *PNAS.* [A, B]
- Ahnert (2017). Structural properties of genotype–phenotype maps. *J R Soc Interface.* [B]
- Allen, Stacey & Bar-Yam (2017). Multiscale information theory and the marginal utility of information. *Entropy.* [A]
- Ancel & Fontana (2000). Plasticity, evolvability, and modularity in RNA. *J Exp Zool.* [B]
- Arita (2004). The metabolic world of *Escherichia coli* is not small. *PNAS.* [E]
- Banzhaf et al. (2016). Defining and simulating open-ended novelty. *Theory in Biosciences.* [E]
- Banzhaf & Yamamoto (2015). *Artificial Chemistries.* MIT Press. [A–E]
- Bedau, Snyder & Packard (1998). A classification of long-term evolutionary dynamics. *ALife VI.* [E]
- Beer (2014). The cognitive domain of a glider in the Game of Life. *Artificial Life.* [D]
- Bennett (1988). Logical depth and physical complexity. [A]
- Bertschinger, Olbrich, Ay & Jost (2006). Information and closure in systems theory. [D]
- Black, Bourrat & Rainey (2020). Ecological scaffolding and the evolution of individuality. *Nature Ecology & Evolution.* [D]
- Boerlijst & Hogeweg (1991). Spiral wave structure in pre-biotic evolution. *Physica D.* [C, D]
- Broido & Clauset (2019). Scale-free networks are rare. *Nature Communications.* [E]
- Channon (2019). Maximum individual complexity is indefinitely scalable in Geb. *Artificial Life.* [E]
- Clark et al. (2011). Degeneracy enriches artificial chemistry binding systems. *ECAL.* [B]
- Conrad (1990). The geometry of evolution. *BioSystems.* [B]
- Corominas-Murtra, Goñi, Solé & Rodríguez-Caso (2013). On the origins of hierarchy in complex networks. *PNAS.* [D, E]
- Decraene, Mitchell & McMullin (2008). Exploring evolutionary stability in a concurrent artificial chemistry. *ECCS.* [cited in B&Y §15.3]
- Dingle, Camargo & Louis (2018). Input–output maps are strongly biased towards simple outputs. *Nature Communications.* [B]
- Dittrich & Speroni di Fenizio (2007). Chemical organisation theory. *Bull Math Biol.* [C, D]
- Dolson, Vostinar, Wiser & Ofria (2019). The MODES toolbox. *Artificial Life.* [E]
- Dorin & McCormack (2002). Self-assembling dynamical hierarchies. *ALife VIII.* [D]
- Draghi, Parsons, Wagner & Plotkin (2010). Mutational robustness can facilitate adaptation. *Nature.* [B]
- Droop & Hickinbotham (2012). A quantitative measure of non-neutral evolutionary activity. *ALife 13.* [E]
- Edelman & Gally (2001). Degeneracy and complexity in biological systems. *PNAS.* [A, B]
- Farrugia-Roberts (2023). Functional equivalence and path connectivity of reducible hyperbolic tangent networks. *NeurIPS.* [B]
- Feldman & Crutchfield (1998). Measures of statistical complexity: why? *Phys Lett A.* [A]
- Fontana & Buss (1994a). "The arrival of the fittest": toward a theory of biological organization. *Bull Math Biol.* [A–D]
- Fontana & Buss (1994b). What would be conserved if "the tape were played twice"? *PNAS.* [C]
- Gao, Buldyrev, Stanley & Havlin (2012). Networks formed from interdependent networks. *Nature Physics.* [D]
- Galván-López, McDermott, O'Neill & Brabazon (2011). Defining locality as a problem difficulty measure in genetic programming. *Genetic Programming and Evolvable Machines* 12(4):365–401. doi:10.1007/s10710-011-9136-3. [checked for this summary]
- Gell-Mann & Lloyd (1996, 2003). Effective complexity. [A]
- Goldford, Hartman, Smith & Segrè (2017). Remnants of an ancient metabolism without phosphate. *Cell.* [E]
- Goldsby, Dornhaus, Kerr & Ofria (2012). Task-switching costs promote the evolution of division of labor and shifts in individuality. *PNAS.* [C, D]
- Greenbury, Schaper, Ahnert & Louis (2016). Genetic correlations greatly increase mutational robustness. *PLoS Comp Biol.* [B]
- Gross & McMullin (2001). Is it the right ansatz? *Artificial Life.* [D]
- Hazen et al. (2024). Molecular assembly indices of mineral heteropolyanions. *J R Soc Interface.* [A, E]
- Hickinbotham, Stepney & Hogeweg (2021). Nothing in evolution makes sense except in the light of parasitism. *R Soc Open Sci.* [C]
- Hintze (2019). Open-endedness for the sake of open-endedness. *Artificial Life.* [E]
- Hordijk, Steel & Dittrich (2018). Autocatalytic sets and chemical organizations. *New J Phys.* [C]
- Hordijk, Steel & Kauffman (2012). The structure of autocatalytic sets. *Acta Biotheoretica.* [C, D]
- Jain & Krishna (2002). Crashes, recoveries, and "core shifts". *Phys Rev E.* [C]
- Kim, Smith, Mathis, Raymond & Walker (2019). Universal scaling across biochemical networks on Earth. *Science Advances.* [E]
- Klein & Hoel (2020). The emergence of informative higher scales in complex networks. *Complexity.* [D]
- Krakauer, Bertschinger, Olbrich, Flack & Ay (2020). The information theory of individuality. *Theory in Biosciences.* [C, D, E]
- Krakauer & Plotkin (2002). Redundancy, antiredundancy, and the robustness of genomes. *PNAS.* [B]
- Kreyssig et al. (2014). Effects of small particle numbers on long-term behaviour in discrete biochemical systems. *Bioinformatics.* [C]
- Kruszewski & Mikolov (2021). Emergence of self-reproducing metabolisms as recursive algorithms in an artificial chemistry. *Artificial Life.* [A, C]
- Lloyd (2001). Measures of complexity: a nonexhaustive list. *IEEE Control Systems.* [A]
- Lloyd & Pagels (1988). Complexity as thermodynamic depth. *Ann Phys.* [A]
- Mathis, Patel, Weimer & Forrest (2024). Self-organization in computation and chemistry: Return to AlChemy. *Chaos.* [A–E]
- McShea (1994). Mechanisms of large-scale evolutionary trends. *Evolution.* [E]
- McShea (1996). Metazoan complexity and evolution: is there a trend? *Evolution.* [A]
- McShea (2001). The hierarchical structure of organisms. *Paleobiology.* [A, D, E]
- Milano, Pagliuca & Nolfi (2019). Robustness, evolvability and phenotypic complexity. *Evolutionary Intelligence.* [B]
- Mingard, Rees, Valle-Pérez & Louis (2025). Deep neural networks have an inbuilt Occam's razor. *Nature Communications.* [B]
- Mitchell (2009). *Complexity: A Guided Tour*, ch. 7. [A]
- Ofria, Adami & Collier (2002). Design of evolvable computer languages. *IEEE TEC.* [B]
- Pfante, Bertschinger, Olbrich, Ay & Jost (2014). Comparison between different methods of level identification. *Adv Complex Syst.* [D]
- Rasmussen, Baas, Mayer, Nilsson & Olesen (2001). Ansatz for dynamical hierarchies. *Artificial Life* 7(4):329–353. Reports D and E list only four authors, and the book dates it 2002. The five-author list and the 2001 date were checked against PubMed 11911785 for this summary. [D, E]
- Ray (1991/1992). An approach to the synthesis of life. *Artificial Life II.* [B]
- Rosas et al. (2020). Reconciling emergences. *PLoS Comp Biol.* [D]
- Rosas et al. (2024). Software in the natural world: a computational approach to hierarchical emergence. arXiv:2402.09090. [D]
- Rothlauf (2006). *Representations for Genetic and Evolutionary Algorithms*, 2nd ed. Springer. [checked for this summary]
- Sayama (2019). Cardinality leap for open-ended evolution. *Artificial Life.* [E]
- Schaper & Louis (2014). The arrival of the frequent. *PLoS ONE.* [B]
- Schulte, Fry, Fast, Weimer & Forrest (2014). Software mutational robustness. *GPEM.* [B]
- Simon (1962). The architecture of complexity. *Proc Am Phil Soc.* [A, D]
- Smith, Kim & Walker (2021). Scarcity of scale-free topology is universal across biochemical networks. *Sci Rep.* [E]
- Solé & Munteanu (2004). The large-scale organization of chemical reaction networks in astrophysics. *EPL.* [E]
- Speroni di Fenizio, Dittrich & Banzhaf (2001). Spontaneous formation of proto-cells in a universal artificial chemistry on a planar graph. *ECAL.* [C]
- Standish (2003). Open-ended artificial evolution. *IJCIA.* [E]
- Stepney (2021). Modelling and measuring open-endedness. *OEE4.* [E]
- Sussmann (1992). Uniqueness of the weights for minimal feedforward nets. *Neural Networks.* [B]
- Szathmáry (2000). The evolution of replicators. *Phil Trans B.* [C]
- Szathmáry (2015). Toward major evolutionary transitions theory 2.0. *PNAS.* [C, D]
- Takeuchi & Hogeweg (2009). Multilevel selection in models of prebiotic evolution II. *PLoS Comp Biol.* [C, D]
- Tononi, Sporns & Edelman (1994). A measure for brain complexity. *PNAS.* [A]
- Tononi, Sporns & Edelman (1999). Measures of degeneracy and redundancy in biological networks. *PNAS.* [A, B]
- Valle-Pérez, Camargo & Louis (2019). Deep learning generalizes because the parameter–function map is biased towards simple functions. *ICLR.* [B]
- van Nimwegen, Crutchfield & Huynen (1999). Neutral evolution of mutational robustness. *PNAS.* [B]
- van Wijk, Stam & Daffertshofer (2010). Comparing brain networks of different size and connectivity density. *PLoS ONE.* [E]
- Vasas, Szathmáry & Santos (2010). Lack of evolvability in self-sustaining autocatalytic networks. *PNAS.* [C]
- Vasas et al. (2012). Evolution before genes. *Biology Direct.* [C]
- Wagner (2005). Robustness, evolvability, and neutrality. *FEBS Letters.* [B]
- Wagner (2008). Robustness and evolvability: a paradox resolved. *Proc R Soc B.* [B]
- Whitacre (2010). Degeneracy: a link between evolvability, robustness and complexity. *Theor Biol Med Model.* [B]
- Whitacre & Bender (2010). Degeneracy: a design principle for achieving robustness and evolvability. *J Theor Biol.* [B]
- Wong et al. (2023). Toward network-based planetary biosignatures. *JGR Planets.* [E]
- Xavier, Hordijk, Kauffman, Steel & Martin (2020). Autocatalytic chemical networks at the origin of metabolism. *Proc R Soc B.* [E]
