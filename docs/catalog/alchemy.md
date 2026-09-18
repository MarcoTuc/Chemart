# AlChemy (lambda-calculus chemistry)

`alchemy` · *Fontana, 1991; Fontana & Buss, 1994*

*Also known as:* *lambda chemistry*, *Fontana's algorithmic chemistry*, *Turing gas*

The chemistry that gave the field its vocabulary. Molecules are lambda-expressions in normal form, and a collision applies one to the other: s1 + s2 -> s1 + s2 + normal_form(s1 applied to s2). Both reactants survive, so every molecule is simultaneously a function and a datum for other functions. Left alone, a soup usually collapses to a single self-copier; ban copying and it instead settles into self-maintaining sets of molecules that collectively produce each other. Those sets are where the notions of closure, self-maintenance and *organisation* come from.

| | |
|---|---|
| **family** | rewriting |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.1 |
| **refs** | [281], [282], [283], [284], doi:10.1007/BF02458289, doi:10.1073/pnas.91.2.757 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): closed normal forms of untyped lambda-expressions modulo renaming of bound variables; species id is the de Bruijn string (^ abstraction, (M)N application, 1-based indices, e.g. ^^2 for λx1.λx2.x1), structure the paper's standardized form λx1.λx2.x1

**R — reactions** (implicit, arity 2): s1 + s2 -> s1 + s2 + NF((s1) s2); elastic if pragmatic reduction fails or a boundary condition rejects the product

**A — reactor**: well-stirred-multiset, ode
 · *dilution:* each reactive collision adds its product and removes a random molecule; M constant (paper 5.2-5.3)

## What you get

```python
net = chemart.generate_network("alchemy", seed=1)
```

```
alchemy: 365 species, 1156 reactions, status=observed
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: analysis, final_state, notation
```

First reactions:

```
^^^^(2)4 + ^^((2)(2)(2)(1)^3)^(1)1 -> ^^^^(2)4 + ^^((2)(2)(2)(1)^3)^(1)1 + ^^^(2)^^((2)(2)(2)(1)^3)^(1)1  [mass-action k=1.0]  (x1)
^^^^((((1)1)(3)1)1)(2)^^6 + ^^(2)^1 -> ^^^^((((1)1)(3)1)1)(2)^^6 + ^^(2)^1 + ^^^((((1)1)(3)1)1)(2)^^^^(2)^1  [mass-action k=1.0]  (x1)
^^(2)1 + ^^(2)^1 -> ^^(2)1 + 2 ^^(2)^1  [mass-action k=1.0]  (x1)
^^^^((((1)1)(3)1)1)(2)^^6 + ^^^(3)(((((3)^3)(3)^4)2)2)(1)2 -> ^^^^((((1)1)(3)1)1)(2)^^6 + ^^^(3)(((((3)^3)(3)^4)2)2)(1)2 + ^^^((((1)1)(3)1)1)(2)^^^^^(3)(((((3)^3)(3)^4)2)2)(1)2  [mass-action k=1.0]  (x1)
^^(2)(((1)^(3)(1)^4)^2)^^3 + ^^(2)^1 -> ^^(2)(((1)^(3)(1)^4)^2)^^3 + ^^(2)^1 + ^^((((2)^^((2)^^^(2)^1)^1)^3)^^4)^1  [mass-action k=1.0]  (x1)
^^((2)2)^^^4 + ^(1)^^2 -> ^^((2)2)^^^4 + ^(1)^^2 + ^^^2  [mass-action k=2.0]  (x1)
^(((((1)1)(1)1)1)1)^^2 + ^^(2)^^1 -> ^(((((1)1)(1)1)1)1)^^2 + ^^(2)^^1 + ^^^^1  [mass-action k=1.0]  (x1)
^^^((1)1)^((4)((3)1)(2)2)1 + ^^^((3)^3)(((3)1)(1)(((1)2)(1)3)2)2 -> ^^^((1)1)^((4)((3)1)(2)2)1 + ^^^((3)^3)(((3)1)(1)(((1)2)(1)3)2)2 + ^^((1)1)^^((((4)2)(3)3)^3)(((((4)2)(3)3)1)(1)(((1)2)(1)((4)2)(3)3)2)2  [mass-action k=1.0]  (x1)
… and 1148 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `soup` | structural | soup: the stochastic flow reactor of Fontana & Buss 5.3, observed reactions with firing counts; closure: every reaction reachable from the seed set (the closure A* of eq. 28), complete or cut off by max_species <br>one of `soup`, `closure` |
| `M` | `int` | `100` | population | reactor size; without terms, the number of distinct random normal forms that seed the reactor (or the closure) <br>`2` … `100000` · *range:* paper: 1000 (reactor capacity, 5.3), 3000 for the L2 merger (6.4.1); book: 1000..3000 |
| `collisions` | `int` | `2000` | population | soup only: number of collisions (elastic ones included) <br>`0` … `10000000` · *range:* paper figs. 4-5: 5*10^5..6*10^5; Mathis et al. 2024: 10^5..6*10^6 |
| `terms` | `list` | `` | structural | explicit seed molecules as closed lambda terms, written λx.(M)N (λ or \) or as de Bruijn ids; reduced to normal form. Soup: M is split equally among them. Overrides the random generator <br>*range:* e.g. the L1 example-1 center [λx1.λx2.λx3.x1, λx1.λx2.λx3.λx4.x2, λx1.λx2.λx3.λx4.λx5.x3] |
| `filter` | `enum` | `none` | selection | functional boundary condition: no-copy declares elastic every collision whose product is identical to one of its two reactants (paper 6.2) <br>one of `none`, `no-copy` · *range:* none gives Level 0 (copiers, hypercycles); no-copy gives Level 1 organisations |
| `forbidden_patterns` | `list` | `` | selection | syntactic boundary conditions: regular expressions searched in the product's de Bruijn id; a match makes the collision elastic and excludes the term from the random seed <br>*range:* paper 6.2.3 bans three consecutive abstractions: ['\^\^\^'] |
| `mediator` | `str` | `` | structural | collision rule Phi of eq. 14: empty for plain application, otherwise a closed term and the product is NF(((Phi)s1)s2) <br>*range:* e.g. λf.λg.((λx1.((x1)λx2.x1)x1)f)(λx1.λx2.((x2)x1)λx3.x3)g (paper eqs. 71-73) |
| `max_steps` | `int` | `10000` | structural | pragmatic reduction time limit: beta contractions allowed before the collision is declared elastic <br>`1` … `1000000` · *range:* paper 5.3: 10,000; Mathis et al. 2024: 500 |
| `max_size` | `int` | `4000` | structural | pragmatic reduction space limit: characters of the whole term during reduction (variable names count 2, λx. 4, the parentheses of an application 2) <br>`2` … `1000000` · *range:* paper 5.3: 4000 characters |
| `max_nf_size` | `int` | `1000` | structural | largest normal form (in the same characters) accepted as a product <br>`2` … `1000000` · *range:* paper 5.3: 1000 characters |
| `p_variable` | `float` | `0.3` | population | random generator (paper 5.3 step 2): probability of a variable node (p1); p_application = 1 - p_variable - p_abstraction <br>`0.0` … `1.0` |
| `p_abstraction` | `float` | `0.4` | population | random generator: probability of an abstraction node (p2) <br>`0.0` … `1.0` |
| `max_depth` | `int` | `7` | population | random generator: nesting level at which a variable is forced <br>`0` … `50` · *range:* paper 5.3: 20; Mathis et al. 2024: 7 |
| `p_bound` | `float` | `0.8` | population | random generator: probability that a variable refers to an enclosing binder (chosen uniformly) rather than to a free name <br>`0.0` … `1.0` |
| `n_free` | `int` | `3` | population | random generator: number of free variable names; free names are bound by leading abstractions (standardization) before reduction <br>`1` … `100` |
| `max_species` | `int` | `200` | structural | closure only: species budget; organisations are usually infinite, so their closure is cut off (status truncated) <br>`1` … `100000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- L0 (no filter): diversity collapses to one self-copier (often λx1.x1) or a small ecology closed under interaction in which every object is a fixed point of some object (fig. 1, eq. 33)
- L1 (copy actions barred): self-maintaining organisations that do not close in the reactor; the simplest is the projector family A(i,j) = λx1...λxi.xj with laws A(i,j) o v = A(i-1,j-1) for j > 1 and A(i,1) o A(k,l) = A(k+i-1,l+i-1) (example 1, eqs. 34-36), whose center is its i smallest objects
- L1 example 2 (three consecutive abstractions barred): two families that act as integers under a difference operation (eqs. 50-61), center {A_-2, A_0, B_0, B_2}
- L2: two L1 organisations coupled by glue (6.4); e.g. organisation A of 6.4.2 is built by T from itself: T o T = 0010T, ..., 10T o 10T = T (fig. 6)
- syntactic convergence of survivors: the system exhausts its novelty (book 9.1, [281])
- Mathis et al. 2024: unfiltered runs with the original code sometimes also give large organisations (tens to hundreds of species) rather than a single copier

## Sources

- Fontana, W. & Buss, L. W. (1994). 'The arrival of the fittest': Toward a theory of biological organization. Bulletin of Mathematical Biology 56(1):1-64. SFI working paper 93-09-055: axioms and model assumptions (4.2, 4.4), collision rule and reaction scheme (5.1, eqs. 13-17), flow reactor equations (5.2, eq. 18), protocol (5.3), closure (5.4), Level 0 (6.1, fig. 1), Level 1 examples (6.2, eqs. 34-68), Level 2 (6.4, eqs. 71-84, fig. 6). https://sfi-edu.s3.amazonaws.com/sfi-edu/production/uploads/sfi-com/dev/uploads/filer/f2/8b/f28b8075-4e73-4696-989e-1cfbe06b2dc8/93-09-055.pdf
- Mathis, C., Patel, D., Weimer, W. & Forrest, S. (2024). Self-organization in computation & chemistry: Return to AlChemy. arXiv:2408.12137: the original random generator and its standardization of free variables (5, 5.1), L0/L1 statistics with the original code. https://arxiv.org/abs/2408.12137
- Fontana, W. & Buss, L. W. (1994). What would be conserved if 'the tape were played twice'? PNAS 91(2):757-761 (abstract only; the PMC copy is a page scan): the three generic levels. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC43028/

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Terms use de Bruijn indices, so identification modulo renaming (paper 4.4.3) is equality; the paper's standardized names (x indexed by binder occurrence) are used for the structure string, and the species id is the de Bruijn string.
- Reduction is normal-order (leftmost-outermost) beta reduction with capture-free substitution, which reaches the normal form whenever one exists (standardization theorem, paper 4.3.1). The paper instead applies Revesz's small-step axioms 4-7 and counts each one as a time unit; so max_steps counts beta contractions and the paper's 10,000 is not the same budget.
- Space is measured in characters of the paper's notation with every variable name counted as 2 characters (exact below 10 nested binders) and checked on the whole term after each contraction; a normal form above max_nf_size is elastic (paper: 4000 and 1000 characters).
- Molecules are closed terms: user terms with free variables are rejected, and random terms are standardized by prefixing an abstraction for each free variable (Mathis et al. 5.1 infer this was used in the original runs). A closed normal form always begins with λ, so the paper's standard filter (operators must be abstraction forms) needs no parameter; the 6.3 variant that admits non-abstraction products is not modelled.
- Random generator: the paper gives only p1 + p2 + p3 = 1 and a maximum nesting level of 20; Mathis et al. add that p1 and p2 vary linearly with depth and that variables are bound to a random ancestor abstraction or left free. Chemart keeps p_variable and p_abstraction constant with depth, and the defaults (0.3, 0.4, p_bound 0.8, 3 free names, depth 7 as in Mathis et al.) are Chemart choices, not published values. As in the paper, the initial objects are distinct normal forms; terms without a normal form or matching a forbidden pattern are redrawn.
- filter: none (Level 0) and no-copy (Level 1, product identical to either reactant is elastic). The v1 choice 'no-copy+no-glue (level 2)' is dropped: the paper has no such filter; its L2 organisations came from merging two L1 populations (6.4.1) or arose rarely under L1 conditions with extra syntactic filters (6.4.2). Syntactic filters are regular expressions on the de Bruijn id (the paper used regular-expression filters on its own notation), e.g. '\^\^\^' for three consecutive abstractions.
- mediator implements the generalized collision rule of eq. 14 for all collisions; the paper's 6.4.1 applies a different rule only to collisions involving glue objects, which is not modelled. The relaxed catalytic scheme of eq. 16, reduced copy efficiency (6.3) and perturbation by random injections are not modelled.
- Rates: all reactions have unit rate per ordered collision (paper 6.2.2 'Kinetics', eq. 18), so a reaction's mass-action constant is the number of ordered pairs (s1, s2), (s2, s1) giving it (1 or 2), with a constant-total outflow.
- Soup: chemart.soup draws two distinct molecules as operator and operand and removes a random molecule after inserting the product; the paper removes one before adding the product, an O(1/M) difference. Initial populations: M distinct random normal forms, or the given terms in equal copies.
- Parameters dropped: seed (an argument of generate_network) and init_generator (callable), replaced by the generator knobs. Defaults are small (M = 100, 2000 collisions); paper values are in range.
- Garbled OCR resolved by computation with this reducer: eq. 34 is A(i,j) = λx1...λxi.xj; the numerals of example 2 are A_i = a^(i+2)A and B_i = a^i B (from eq. 53 and h* = aaA), and the 'otherwise' branch of eqs. 59-61 is k = i - j - 1; law 81 of L2 example 2 is 1v1 o 1v2 = v1 o v2 (the scan drops 'o v2'; this is what the cycles of fig. 6 require); the fig. 1 terms were read as given in the tests and match the caption's B o C.

## Notes

Historically the most important AC: the notions of closure, self-maintenance and organisation come from here. The closure of an organisation is usually infinite, so method closure is normally truncated; analysis.self_maintaining reports whether every species of the returned set is produced inside it (eq. 38).

---

*Specification: `catalog/chemistries/alchemy.yaml` · generator: `chemart/chemistries/alchemy.py` · tests: `tests/chemistries/test_alchemy.py`*
