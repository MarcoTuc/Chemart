# Chemical Casting Model (CCM)

`ccm` · *Kanada, 1992-1996*

*Also known as:* *chemical casting model*, *chemical computation model*, *CCM with FAM*, *annealed CCM*

Emergent problem solving from purely local rules. The whole working memory is a single candidate solution - one atom per queen, or per graph vertex - and rules rewrite small groups of atoms, with the untouched atoms acting as catalysts. A rule fires only if it does not lower a *local* order degree, so nothing ever computes a global objective. Yet the system reliably reaches global optima like an eight-queens solution or a valid map colouring, because random local improvement plus the occasional frustrated move is enough to escape traps.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 17.2.2 |
| **refs** | [438], [439], https://ieeexplore.ieee.org/document/538377/, https://www.kanadas.com/Papers/CCM-papers.html |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (explicit): atoms with a type and a state, optionally connected by links; the whole working memory is ONE point of the search space. A species is an atom in a state: q<row>=<column> for a queen (n-queens), v<i>=c<k> for a vertex with colour k (graph-coloring); structure '(queen, row=r, column=c)' or '(vertex, id=i, color=k)'. Rows, columns, vertices and colours are numbered from 0

**R — reactions** (explicit, arity variable): LHS -> RHS: the matched atoms are rewritten and the unchanged ones are catalysts. n-queens swap rule with k catalysts: q<a>=ca + q<b>=cb + k queens -> q<a>=cb + q<b>=ca + the same k queens; moving rule: q<a>=c + all other queens -> q<a>=c' + all other queens. graph-coloring: v<i>=c1 + k linked neighbours (or all of them) -> v<i>=c2 + the same neighbours. An instance reacts only if the instance order degree (sum of the pairwise LODs of the matched atoms) does not decrease

**A — reactor**: well-stirred-multiset
 · *dilution:* none - the working memory keeps one atom per queen or vertex; a reaction rewrites atom states

## What you get

```python
net = chemart.generate_network("ccm", seed=1)
```

```
ccm: 64 species, 187 reactions, status=observed
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: analysis, final_assignment, initial_assignment, instance
```

First reactions:

```
q3=3 + q4=4 + q6=6 -> q3=4 + q4=3 + q6=6  (x1)
q7=7 + q0=0 + q1=1 -> q7=0 + q0=7 + q1=1  (x1)
q6=6 + q7=0 + q1=1 -> q6=0 + q7=6 + q1=1  (x1)
q5=5 + q4=3 + q0=7 -> q5=3 + q4=5 + q0=7  (x1)
q0=7 + q3=4 + q7=6 -> q0=4 + q3=7 + q7=6  (x1)
q1=1 + q4=5 + q2=2 -> q1=5 + q4=1 + q2=2  (x1)
q0=4 + q6=0 + q2=2 -> q0=0 + q6=4 + q2=2  (x1)
q3=7 + q6=4 + q2=2 -> q3=4 + q6=7 + q2=2  (x1)
… and 179 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `problem` | `enum` | `n-queens` | structural | constraint satisfaction problem with its published rule and LOD: n-queens (HICSS-27 sec. 4, Queens_Sort) or graph-coloring (FUZZ-IEEE'95 sec. 3, Kanada 1996 sec. 3.2, USAmap_color) <br>one of `n-queens`, `graph-coloring` |
| `N` | `int` | `8` | structural | n-queens only: board size and number of queen atoms <br>`4` … `1000` · *range:* HICSS-27 fig. 7: 4..50; figs. 8-10: 8; fig. 12: 6 (the demo requires N >= 4) |
| `graph` | `enum` | `usa-mainland` | structural | graph-coloring only: usa-mainland (the 48 contiguous states, 106 borders, of FUZZ-IEEE'95 fig. 5 and USAmap_color); random (G(V, p), the model of the DSJC benchmark graphs of Kanada 1996 table 2); custom (edges) <br>one of `usa-mainland`, `random`, `custom` |
| `V` | `int` | `30` | structural | graph-coloring with graph=random: number of vertices <br>`1` … `10000` · *range:* Kanada 1996 table 2: DSJC125.x and DSJC250.x |
| `edge_probability` | `float` | `0.1` | structural | graph-coloring with graph=random: probability of each edge <br>`0.0` … `1.0` · *range:* DSJC: 0.1, 0.5, 0.9 |
| `edges` | `list` | `` | structural | graph-coloring with graph=custom: undirected edges [[u, v], ...] on vertices 0..max id |
| `colors` | `int` | `4` | structural | graph-coloring only: number of colours <br>`2` … `1000` · *range:* USA map: 4; DSJC125.1: 5-6, DSJC125.5: 17-18, DSJC250.5: 29-31 |
| `rule` | `enum` | `single-catalyst` | structural | catalysts of the rule, i.e. its locality. n-queens: column swap of two queens with 0-3 random other queens as catalysts, or variable-catalyst = move one queen to another column with all other queens as catalysts. graph-coloring: recolour a vertex with 0-3 random linked neighbours as catalysts, or all its neighbours. More catalysts: fewer reactions, but easier trapping in local maxima <br>one of `no-catalyst`, `single-catalyst`, `double-catalyst`, `triple-catalyst`, `variable-catalyst` · *range:* HICSS-27 figs. 9-10: Nc = 0..3 swap catalysts; demos: no, single, double and variable catalyst rules |
| `acceptance` | `enum` | `non-decreasing` | selection | reaction test on the instance order degree: non-decreasing IOD_before <= IOD_after (Kanada's definition and demos) or increasing IOD_before < IOD_after (the book's wording, HICSS-27 footnote 3); frustration is subtracted from IOD_before in both <br>one of `non-decreasing`, `increasing` |
| `frustration` | `bool` | `True` | thermodynamic | frustration accumulation method (FAM), the book's annealing variant; false is the original CCM |
| `f0` | `float` | `1e-05` | thermodynamic | initial frustration of every atom, restored after each reaction of the atom (must be > 0 with frustration) <br>`0.0` … `1000.0` · *range:* demos and Kanada 1996 table 2: 1e-5 (down to 1e-45 on DSJC250.5); fig. 15: 0.8 |
| `c` | `float` | `2.0` | thermodynamic | frustration growth factor: f -> c f after each failed test of an instance with unsatisfied constraints <br>`1.0` … `1000.0` · *range:* demos and Kanada 1996: 2; figs. 16-17: 1.1-4; fig. 15: 1.05 |
| `initial` | `enum` | `ordered` | population | initial working memory: ordered = all queens on the diagonal (column = row, HICSS-27 sec. 4.1, Queens_Sort) or all vertices colour 0 (FUZZ-IEEE'95 sec. 3, USAmap_color); random = a random permutation of columns (HICSS-27 sec. 4.2) or random colours <br>one of `ordered`, `random` |
| `max_tests` | `int` | `200000` | population | budget of rule tests (LHS matches); rules without catalysts never terminate and always stop here <br>`1` … `1000000000` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- convergence to global optima despite purely local evaluation: 8-queens and the USA map are solved in every run (HICSS-27 sec. 4.2, FUZZ-IEEE'95 sec. 3)
- conflicting systems: in n-queens and colouring a reaction that raises the local IOD can lower the global order degree, so GOD is not monotone (HICSS-27 fig. 8, sec. 5.2)
- rules without catalysts never stop: pure random walk (HICSS-27 sec. 5.3, FUZZ-IEEE'95)
- more catalysts: fewer reactions to a solution but more matches per reaction and easier trapping (HICSS-27 figs. 9-10, fig. 12)
- escaping local maxima: the 6-queens state of HICSS-27 fig. 12 is a GOD local maximum that the single-catalyst swap escapes
- FAM lets the variable-catalyst rule escape stationary non-solutions that trap it without frustration (FUZZ-IEEE'95, Kanada 1996 sec. 4)
- near-lock-free parallelism: access conflicts act as harmless low-intensity noise (book 17.2.2; Kanada 1996 sec. 5; not simulated)

## Sources

- Kanada, Y. & Hirokawa, M. (1994). Stochastic problem solving by local computation based on self-organization paradigm. 27th HICSS, pp. 82-91 (book [439]). Sec. 3: atoms, rules, instances, LODs (self and mutual order degrees), instance order condition IOD_before <= IOD_after; sec. 4: the N-queens swap rule with a catalyst and its LOD (figs. 4-6); sec. 5: GOD, conflicting vs cooperative systems, catalysts (figs. 8-11), local maxima (fig. 12); table 1. https://www.kanadas.com/CCM/HICSS27/probsolving-hicss.new.pdf
- Kanada, Y. (1995). Combinatorial problem solving using randomized dynamic tunneling on a production system. IEEE SMC'95, vol. 4, pp. 3784-3789 (book [438]). Sec. 2.1: atoms with type and state, links, reaction condition (sum of LODs of the atoms involved not decreased); sec. 2.2-4: CCM* (random rule composition) on 0-1 integer programming. https://www.kanadas.com/CCM/SMC95/tunneling.pdf
- Kanada, Y. (1995). Fuzzy constraint satisfaction using CCM - a local information based computation model. FUZZ-IEEE/IFES'95, pp. 2319-2326. Sec. 3: graph-colouring rule (fig. 3) and LOD, zero/two/variable catalysts (fig. 4), USA mainland map, mean degree 4.42, correct solution in every run, FAM needed by the variable-catalyst rule; table 2 (112 reactions, 4406 matchings on average). https://www.kanadas.com/CCM/FUZZ-IEEE95/fuzzycs-ccm.upd.pdf
- Kanada, Y. (1996). Constraint satisfaction by parallel optimization of local evaluation functions with annealing (manuscript, 2/4/96). Sec. 3.2: colouring rule/LOD, catalysts, termination by a number of failed trials; sec. 4: FAM (f -> c f on failure with unsatisfied constraints, reset to f0 on reaction, f0 = 1e-5, c = 2, fig. 11); sec. 6: DSJC results (table 2, figs. 16-17). https://www.kanadas.com/CCM/LargeCSP/largecsp.pdf
- Kanada, Y. (1996, JavaScript port 2025). N queens problem and sorting using CCM (Queens_Sort) and USA mainland map coloring using CCM (USAmap_color): reference programs with the no/single/double swap rules, the variable-catalyst moving rule, the recolouring rules, FAM (f0 = 1e-5, c = 2), the termination test and the USA border graph. https://www.kanadas.com/ccm/queens-sort/Queens_Sort.js , https://www.kanadas.com/ccm/coloring/USAmap_color.js , method pages https://www.kanadas.com/ccm/queens-sort/method.html , https://www.kanadas.com/ccm/coloring/method.html

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Book error: 17.2.2 credits the annealing variant with frustration to [438]. [438] (SMC'95) introduces CCM*, a randomized tunnelling by dynamic rule composition, and uses no frustration; FAM is from Kanada's SWoPP'94 paper (Japanese, not used), FUZZ-IEEE'95 sec. 3-4 and the 1996 manuscript sec. 4, which are followed here. CCM* is not implemented.
- Book wording: a rule 'may be applied if the sum of LODs on the LHS is smaller than on the RHS' (strict). Kanada's papers say the IOD must not decrease (HICSS-27 sec. 3.3, SMC'95 sec. 2.1), with the strict test as an allowed variant (HICSS-27 footnote 3). Both are offered as acceptance; the default is Kanada's non-decreasing test.
- Reaction test and FAM follow the demos exactly: an instance whose matched pairs are all satisfied never reacts (and its frustration does not grow), so solved states are stationary under the non-decreasing test; otherwise it reacts if O' >= O - F (or O' > O - F), F = summed frustration of the changed atoms (catalysts not counted, 1996 sec. 4); on failure f -> c f for each changed atom, on reaction f -> f0. The 1996 manuscript prints O - F <= O' - F' with F' after the update; its own worked example (fig. 11) and the demos compare against O - F only, which is used. The alternative of resetting f when constraints become satisfied by another reaction (1996 footnote) is not implemented.
- IOD pairs: only pairs involving a changed atom are summed (pairs of two catalysts are unchanged and cancel, 1996 sec. 3.2); the n-queens swap counts the swapped pair too, as Queens_Sort does. Rules without catalysts evaluate no pair and react unconditionally (demo swap0/change_color0; 1996 sec. 3.2), so they random-walk and only max_tests stops them.
- n-queens LOD: HICSS-27 fig. 4 defines o only by diagonals, since the swap rule keeps one queen per row and column; Queens_Sort.order also sets o = 0 for a shared column, which the moving rule needs. The demo's definition is used; under swap rules the two coincide.
- n-queens rules: swaps with 0-3 catalysts (HICSS-27 figs. 9-10 use Nc = 0..3; the demo 0..2) and the demo's variable-catalyst moving rule. Atoms are drawn uniformly without repetition, as the demos' rejection loops do. The moving rule draws a new column different from the current one (the demo can redraw the same column, which counts as a reaction with no change); this matches USAmap_color.select_color and records no empty reactions. A rule matching more than N queens is rejected (HICSS-27 sec. 5.3).
- graph-coloring rules: the demo's recolouring with one or two random neighbours as catalysts (vertices with too few neighbours are skipped as failed tests, USAmap_color.change_color1/2), extended to three, or all neighbours (variable catalyst, FUZZ-IEEE'95 fig. 4c). The new colour is drawn uniformly among the other colours (select_color); the papers draw from all colours but a same-colour draw can never change the state. The map's graph is the demo's; its mean degree 212/48 = 4.42 matches FUZZ-IEEE'95 sec. 3.
- Termination as in the demos' loopTick: with threshold T (1000 + N^2 for queens, 1000 for colouring), when T tests in a row fail the run stops if fewer than 20 T tests were made, else T doubles (the 1996 manuscript's 'maximum number of trials increased to 1/20 of the total number of LHS tests'). max_tests is a hard budget. GOD is never used by the dynamics; it is only reported (analysis).
- Published counts: FUZZ-IEEE'95 table 2 reports 112 reactions and 4406 LHS matchings on average for the USA map with the variable-catalyst rule and FAM, run with 'a slightly modified version of the rule and LOD' in SOOC-94; Chemart's demo-based rule gives about 117 reactions and 2500 tests (10 seeds). The reaction count is checked in the tests (60-200); the matching count is not, since the modification is not described.
- Initial states follow the demos: queens on the diagonal (the GOD minimum, HICSS-27 sec. 4.1) and all vertices colour 0 (FUZZ-IEEE'95); HICSS-27 measured fig. 7 from random layouts, offered as initial=random.
- Not implemented: the TSP, 0-1 knapsack/integer programming and sorting casters (HICSS-27 table 1 gives only counts; SMC'95 gives the 0-1 IPP rule but its runs need CCM*), fuzzy colouring (FUZZ-IEEE'95 sec. 4-5) and parallel scheduling (1996 sec. 5). Timings (SOOC on Macintosh, C on a Cray CS6400) are not reproducible; tests check the qualitative claims and the counts.
- Observed network: species are the initial atom states and every atom state produced; reactions are the distinct accepted rule applications (changed atoms first, then catalysts on both sides) with counts; initial_state has one of each initial atom state. extras.analysis: tests, reactions, uphill_reactions (accepted only thanks to frustration), terminated, god_initial/final/max, mod_final, solved, first_solution_reaction, god (GOD after each reaction, first 20000), mean_frustration_final; extras.instance (N, or vertices, edges, colors, mean_degree), initial_assignment and final_assignment. No rate constants: CCM has none.
- v1 parameters: the LOD callable becomes the LOD fixed by the problem; catalysts (bool) becomes the rule enum; c keeps its meaning (FAM growth factor) with f0 added.

## Notes

Kanada's model for emergent computation: a production system whose rules test only local evaluation functions and fire in random order, so that global order (solutions of constraint satisfaction problems) emerges without a global objective. CCM* (SMC'95) adds randomized rule composition; FAM adds local annealing without a global temperature.

---

*Specification: `catalog/chemistries/ccm.yaml` · generator: `chemart/chemistries/ccm.py` · tests: `tests/chemistries/test_ccm.py`*
