# P systems (membrane computing)

`p-systems` · *Paun, 1998*

*Also known as:* *membrane systems*, *transition P systems*

Computing with nested membranes. Objects are symbols held in the regions of a membrane tree, and rules rewrite multisets of them while targeting the products - keep them here, send them out, or push them into an inner membrane - and may dissolve the membrane they fired in. Every rule that can fire in a step does, so the parallelism is maximal. The structure computes: halting configurations of a fixed small system can yield exactly the square numbers.

| | |
|---|---|
| **family** | rewriting |
| **kind** | formalism |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.5; multicompartment SSA 4.3; biology 18.3.2 |
| **refs** | [678], [679], [680], [681], [683], [661], [718], [150], doi:10.1006/jcss.1999.1693, doi:10.1016/S0304-3975(02)00136-6 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `compartments`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): symbol objects over an alphabet O, located in the regions of a tree of nested membranes [1 [2 ]2 [3 ]3 ]1 (1 = skin); multisets written as strings (a^2 b^3 c). Species id = object@region (a@2), with env for the environment outside the skin

**R — reactions** (explicit, arity variable): multiset rewriting u -> v with target indications on the products, here (default), out (the surrounding region, the environment for the skin), in (a non-deterministically chosen inner membrane) or in_j, optionally with delta (dissolve the membrane of the region): aab -> a b_out c c_in; d -> d delta. Priority relations rho_i order the rules of a region. Chemart syntax: 'r1: 2 a + b -> a + b@out + c + c@in', 'd -> d + delta', priorities 'r1 > r3'

**A — reactor**: maximally-parallel, compartments
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("p-systems", seed=1)
```

```
p-systems: 8 species, 4 reactions, status=complete
provides: catalysts, compartments, initial-state, stoichiometry, topology
seed: 1
extras: analysis, compartments, dissolution, membranes, output_region, priorities, priority_semantics, reaction_rules, rules, system
```

First reactions:

```
d@1 + c@1 + c'@1 -> a@3
a@2 + c@2 -> c'@2
a@2 + c'@2 -> c@2
d@2 -> d@2
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `system` | `enum` | `divisibility` | structural | published system: divisibility (book eq. 9.21 and fig. 9.3: is n a multiple of k? output membrane 3 is empty iff it is), n-squared (Paun 2006 fig. 3: generates n^2, n >= 1, sent to the environment), custom (membranes, objects, rules, priorities, output) <br>one of `divisibility`, `n-squared`, `custom` |
| `n` | `int` | `7` | population | divisibility only: copies of a in membrane 2 (the number tested) <br>`1` … `1000000` · *range:* positive; the computation halts after about n/k + 3 steps, so raise max_steps for large n/k |
| `k` | `int` | `3` | population | divisibility only: copies of c in membrane 2 (the divisor) <br>`1` … `1000000` |
| `membranes` | `str` | `` | spatial | custom only: membrane structure mu <br>*range:* labelled brackets, e.g. '[1 [2 ]2 [3 ]3 ]1'; labels are letters, digits or _, 'env' and '0' are reserved |
| `objects` | `dict` | `` | population | custom only: initial multiset w_i per region <br>*range:* e.g. {'2': '7 a + 3 c + d'} or {'2': {'a': 7, 'c': 3, 'd': 1}} |
| `rules` | `dict` | `` | structural | custom only: rule set R_i per region, in the explicit reaction syntax with target suffixes <br>*range:* e.g. {'2': ['r1: a + c -> c2', 'r3: d -> d + delta'], '1': ['e -> e@out', 'x -> y@in3']}; targets @here (default), @out, @in, @in<label>; delta or δ dissolves; the optional label defaults to r<position> |
| `priorities` | `dict` | `` | structural | custom only: priority relation rho_i per region over its rule labels <br>*range:* e.g. {'2': ['r1 > r3', 'r2 > r3']}; chains 'r1 > r2 > r3' allowed; closed transitively, cycles rejected |
| `output` | `str` | `` | structural | custom only: output region i_o <br>*range:* a membrane label, or 'env' / '0' for the environment |
| `priority` | `enum` | `strong` | structural | semantics of the priority relations in the computation: strong (a lower rule is not used while a higher one is applicable; needed by the book's divisibility system) or weak (the higher rule is used maximally, the lower one on the remaining objects); Paun 2006 sec. 11 <br>one of `strong`, `weak` |
| `max_steps` | `int` | `1000` | population | step budget of the maximally parallel computation recorded in extras.analysis (does not change the network) <br>`1` … `10000000` · *range:* a halting n-squared computation with result n^2 takes 2n + 2 steps; non-halting computations stop here |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- divisibility (book 9.5, eq. 9.21): r1 and r2 alternately consume k copies of a per step, r3 dissolves membrane 2 once a is exhausted, and r4 puts one a into membrane 3 exactly when n is not a multiple of k
- N(Pi) = {n^2 | n >= 1} for the system of Paun 2006 fig. 3: every halting computation sends (n+1)^2 copies of e to the environment
- maximal parallelism (Paun 2006 sec. 4): on a^5 b^2 c^6 the rules aab -> abcc and bb -> aac give either a^3 b^2 c^10 or a^7 c^7 in one step
- strong vs weak priority (Paun 2006 sec. 11): with ff -> f > cf -> cd delta on fffc, the strong reading uses only ff -> f, the weak reading also uses cf -> cd delta on the remaining fc
- target indications: aab -> a b_out c c_in is used twice on a^5 b^2 c^6 (Paun 2006 sec. 7)
- computational universality, space-for-time trade-off via membrane division or worm objects, and polynomial solutions of NP-complete problems are reported in [679, 683] and are not reproduced

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, MIT Press, section 9.5: notation, target indications, dissolution, eq. 9.20, and the divisibility system eq. 9.21 / fig. 9.3 (primes and the in_3 subscript checked on the rendered PDF page 190, because the text extraction drops them).
- Paun, G. (2006). Introduction to Membrane Computing. In G. Ciobanu, G. Paun, M. J. Perez-Jimenez (eds.), Applications of Membrane Computing, Springer, pp. 1-42 (doi:10.1007/3-540-29937-8_1): maximal parallelism and the a^5 b^2 c^6 example (sec. 4), target indications and dissolution (sec. 7), the formal definition of a transition P system (sec. 8), computations, halting, internal and external output and the n^2 system of fig. 3 (sec. 9), strong and weak priorities with the ff -> f / cf -> cd delta example (sec. 11). Preprint used: https://cs.ioc.ee/yik/schools/win2007/paun/IntrodMC.pdf and http://bioinfo.uib.es/~recerca/BUM/paun.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Original papers not reached: Paun's TUCS Report 208 / JCSS 61 (2000) [678], the 2002 book [680] (with the yes/no version of the divisibility system), the Paun-Rozenberg guide [681] and ppage.psystems.eu (connection refused) could not be downloaded; the formal definition and semantics come from Paun's 2006 introduction, which restates them. The yes/no divisibility system of [680] is therefore not offered.
- Scope: transition P systems with symbol objects, targets here/out/in/in_j, dissolution and priorities (book 9.5). Membrane division and active membranes [679, 683], symport/antiport, promoters/inhibitors, string objects and the stochastic MGA variants [661, 718] (book 4.3, 18.3.2) are not implemented: the book gives no algorithmic detail or rate constants for them. The v1 A.reactor 'ssa' is dropped.
- Priorities: the book states rho_2 = {r1 > r3, r2 > r3} and that r3 fires 'when neither r1 nor r2 can be fired', which is Paun's strong interpretation; under the weak one r3 would fire in the first step, as soon as the c are used up, and the computation would give a wrong answer. Default priority = strong; weak is available for Paun 2006 sec. 11.
- Non-determinism: any maximal multiset of rules is allowed by the definition; Chemart samples one by repeatedly drawing a still-applicable rule with weight equal to the number of times it can still be applied and applying it a uniform number of times in [1, that], until no rule can take the remaining objects. An object sent 'in' goes to a uniformly chosen inner membrane. This is one admissible scheduler, not a published distribution.
- Step semantics (Paun 2006 secs. 7-9): all regions choose their multisets on the configuration at the start of the step; products are delivered by the membrane structure at the start of the step (objects leaving the skin go to the environment and never return); then every membrane whose region used a delta rule dissolves, and its objects and inner membranes join the nearest upper membrane not dissolved in the same step; the skin can never dissolve (a delta rule in the skin is rejected). An 'in' target with no (such) inner membrane makes the rule inapplicable.
- Halting and output: a computation halts when no rule is applicable in any region; the result is the number of objects in the output region in the halting configuration, or of objects sent to the environment when i_o = env (Paun uses i_o = 0, also accepted). A computation whose output membrane was dissolved, or that does not halt within max_steps, has result null. The terminal-object restriction T is not modelled (all objects count); catalysts are not declared separately.
- Network: every rule becomes reactions over object@region species. Because dissolution changes the tree, out/in/in_j targets resolve to every region reachable in some structure (a further ancestor when the intermediate membranes have delta rules, a deeper membrane when the membranes between can dissolve), and an 'in' with several possible membranes gives one reaction per distribution of the objects, so the network is the complete set of reactions the system can perform. Status complete.
- Facts without a formula: products of a delta rule keep the rule's own region (d -> d delta is the reaction d@2 -> d@2); the dissolution itself and the release of contents are recorded in extras.dissolution, priorities in extras.priorities (transitive closure, [higher, lower] pairs) and extras.priority_semantics. Reactions have no rate: P systems define none.
- Book notation: the text extraction of eq. 9.21 loses the primes and the subscript; the rendered page reads R_1 = {r4: dcc' -> a_in3}, R_2 = {r1: ac -> c', r2: ac' -> c, r3: d -> d delta}, w_2 = a^n c^k d, i_o = 3. The book's 'adapted from [150, 680]' system is used as printed.
- Paun 2006 fig. 3 system (n-squared): i_o is the environment, as the text says; rules r1..r3 of region 3 are a -> ab, a -> b delta, f -> ff; region 2 has b -> d, d -> de, ff -> f, cf -> cd delta; region 1 has e -> e_out, f -> f. It has no priorities: computations in which cf -> cd delta is used too early never halt because of f -> f.
- v1 parameters: Sigma (the alphabet is implied by the objects and rules), mu, w_i, R_i (callable), rho_i and i_o became the system choice plus the custom JSON/text fields membranes, objects, rules, priorities and output; n and k give the divisibility input; priority and max_steps control the recorded computation.
- extras: system, membranes (canonical structure), compartments (per membrane: parent, children, skin, dissolvable, rules, species; plus env), rules (id, region, label, text, dissolves, reaction indices), reaction_rules (rule id of each reaction), priorities, priority_semantics, dissolution, output_region, inapplicable_rules (rules no structure lets fire), and analysis (halted, steps, output, result, final structure and configuration, dissolutions, rule applications, a trace of the first 100 steps; divisibility adds n_multiple_of_k, n-squared result_is_square).

## Notes

Formally Pi = (Sigma, mu, w_1..w_m, (R_1,rho_1)..(R_m,rho_m), i_o). The module also exports parse_membranes, build_system, instantiate, select and run, and the published systems divisibility(n, k) and N_SQUARED.

---

*Specification: `catalog/chemistries/p-systems.yaml` · generator: `chemart/chemistries/p_systems.py` · tests: `tests/chemistries/test_p_systems.py`*
