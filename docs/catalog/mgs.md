# MGS

`mgs` · *Giavitto & Michel, 2001*

*Also known as:* *(encore) un modele general de simulation*, *topological collections and transformations*

Rewriting where the neighbourhood relation is part of the data. A collection is not just a bag - it can be a sequence, a grid, a graph, or a group-based field - and a rule matches a *path of neighbouring* elements and substitutes for it. One formalism therefore covers sorting a sequence, a cellular automaton on a lattice, and restriction-enzyme cuts on a string, differing only in which topology the collection carries. Useful to Chemart as a way of thinking about space abstractly.

| | |
|---|---|
| **family** | rewriting |
| **kind** | formalism |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.6; biology 18.3.2 |
| **refs** | [317], [319], [320], [321], [322], doi:10.1016/S1571-0661(04)00293-2, doi:10.1007/3-540-45833-6_12 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `space`, `mass-conservation`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): a topological collection: values organised by a neighbourhood relation. Monoidal collections (set, bag: every element neighbours every other; seq: the right neighbour), group-based data fields (grid <north, east>, hexagon <east, north, northeast; east + north = northeast>) whose positions exist a priori and may be <undef>, and nested collections (a TUBE bag of DNA sequences). Species are values (set, bag), whole strands (tube) or position-qualified values value@position (seq, grid, hexagon: 3@2, true@1_4, undef@0_0)

**R — reactions** (explicit, arity variable): trans T = { pattern / guard => replacement; ... }: a rule selects a path of neighbouring elements matching the pattern (x, y: y neighbours x; x |d> y: y is the d-neighbour of x; X+: a non-empty subsequence) whose guard holds, and substitutes the replacement for it. Book: x, y / x > y => y, x sorts a sequence

**A — reactor**: graph-rewrite, maximally-parallel
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("mgs", seed=1)
```

```
mgs: 21 species, 25 reactions, status=complete
provides: initial-state, mass-conservation, space, stoichiometry, topology
seed: 1
extras: conservation, program, reaction_rules, rules, run, space
```

First reactions:

```
4@0 + 2@1 -> 2@0 + 4@1
5@2 + 1@3 -> 1@2 + 5@3
4@1 + 1@2 -> 1@1 + 4@2
2@1 + 1@2 -> 1@1 + 2@2
5@3 + 3@4 -> 3@3 + 5@4
4@0 + 1@1 -> 1@0 + 4@1
2@0 + 1@1 -> 1@0 + 2@1
4@1 + 2@2 -> 2@1 + 4@2
… and 17 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `program` | `enum` | `sort` | structural | published program: sort (book 9.6, x, y / x > y => y, x on a seq), sieve (Cohen 2003, Eratosthenes' sieve on the set 2..30), restriction-enzymes (book 9.6 and Giavitto & Michel 2002 sec. 3.1, EcoRI and Void on a TUBE of DNA strands), eden (2002 sec. 3.2, Eden growth x, <undef> / x => x, true on a 7x7 grid), bead-sort (Cohen 2003 fig. 4, beads 3, 2, 4, 2 falling north to south on a grid), turn (2002 fig. 2, a five-cell cross turned by 90 degrees), custom (collection, rules, initial) <br>one of `sort`, `sieve`, `restriction-enzymes`, `eden`, `bead-sort`, `turn`, `custom` |
| `collection` | `enum` | `default` | spatial | collection type; default keeps the program's own (sort seq, sieve set, restriction-enzymes tube, eden, bead-sort and turn grid); required for custom. The 2002 paper runs the same eden transformation on grid and hexagon <br>one of `default`, `set`, `bag`, `seq`, `grid`, `hexagon`, `tube` |
| `rules` | `list` | `` | structural | custom only: the transformation, one rule per string '[Label =] pattern => replacement' in priority order; the pattern and guard language is documented in chemart/chemistries/mgs.py (no Python evaluation) <br>*range:* e.g. ['x, y / x > y => y, x'] or ['EcoRI = X+, "G","A", Y+ => (X,"G") :: ("A",Y) :: ():TUBE'] |
| `initial` | `list` | `` | structural | initial collection; empty keeps the program's own (sort [4, 2, 5, 1, 3], sieve 2..30, the two published DNA strands, eden one true cell in the middle of a 7x7 grid, the published bead rows and cross); required for custom <br>*range:* seq/set/bag: flat list of integers, booleans or names; grid/hexagon: list of equal-length rows with null for <undef> (row 0 to the north); tube: strings or lists of atoms |
| `torus` | `bool` | `` | spatial | grid and hexagon only: close the window into a torus (the GBF with equations rows north = 0, columns east = 0) instead of a bounded window of the infinite GBF |
| `strategy` | `enum` | `maximal-parallel` | stochastic | application strategy of the run recorded in extras.run: maximal-parallel (each rule in priority order takes a random maximal set of non-intersecting matches, all substituted at once; MGS default) or asynchronous (one random match of the highest-priority applicable rule per step); the closure does not depend on it <br>one of `maximal-parallel`, `asynchronous` |
| `max_species` | `int` | `5000` | structural | species budget of the closure; growth rules on tubes or bags can have infinite closures, which are truncated <br>`1` … `1000000` · *range:* measured closures: sort 21 species / 25 reactions, sieve 29 / 52, restriction-enzymes 6 / 2, eden 97 / 164 (hexagon 97 / 234), bead-sort 22 / 4, turn 21 / 252; all under 1 s |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- book 9.6 and Cohen 2003: x, y / x > y => y, x iterated to its fixpoint sorts a sequence; in the network the sorted configuration is inert
- book 9.6 and Giavitto & Michel 2002 sec. 3.1: EcoRI cuts CCCGAATTCAA into CCCG and AATTCAA; React[fixpoint] on the tube CCCGAATTCAA, TTGAATTCGGG returns AATTCAA, TTG, CCCG, AATTCGGG
- Cohen 2003: x, y/(y mod x = 0) => [x] iterated on the set 2..n leaves the primes up to n
- Giavitto & Michel 2002 sec. 3.2: Eden growth x, <undef> / x => x, true spreads from one cell; the same transformation runs on grid and hexagon
- Giavitto & Michel 2002 fig. 2: Turn rotates a five-cell cross by 90 degrees per step
- Cohen 2003 fig. 4: bead-sort lets beads fall along north, sorting rows 3, 2, 4, 2 into 2, 2, 3, 4
- Giavitto et al. 2011: MGS also models flocking (Delaunay neighbourhoods) and epithelial growth (not implemented)

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 9.6 (topological collections; trans sort = { x, y / x > y => y, x }; the EcoRI rule and its cut of CCCGAATTCAA into CCCG and AATTCAA) and 18.3.2 (one paragraph on the biological models of [321]).
- Giavitto, J.-L. & Michel, O. (2002). Data structure as topological spaces. Unconventional Models of Computation (UMC 2002), LNCS 2509:137-150. Collection types and neighbourhoods (sec. 2.1: set, bag, seq, GBF Grid2 and Hexagon), the path pattern grammar (sec. 2.2), substitution, Newtonian and Leibnizian collections and the iteration options (sec. 2.3), the restriction enzymes program with its published input and result (sec. 3.1), the Eden model on grid and hexagon (sec. 3.2, fig. 1) and the Turn transformation (fig. 2). http://repmus.ircam.fr/_media/giavitto/export/giavitto__umc02__data_structure_as_topological_spaces.pdf
- Cohen, J. (2003). Typing rule-based transformations over topological collections. RULE 2003, ENTCS 86(2). Rule application semantics (sec. 2: rules in order, each on a maximal set of non-intersecting occurrences), the bubble sort and the Eratosthenes sieve on a set (sec. 2.1), the default rule x => [x] (sec. 4.1), the bead-sort on a grid (sec. 5.2, fig. 4) and the strategies (sec. 5.3). https://arxiv.org/abs/0912.4877
- Giavitto, J.-L., Michel, O. & Spicher, A. (2011). Interaction based simulation of dynamical system with a dynamical structure (DS)2 in MGS. SCS 2011. Guards, the linearity of patterns, the maximal-parallel and Gillespie strategies, rule priorities by declaration order, Delaunay collections. http://repmus.ircam.fr/_media/giavitto/export/giavitto_scs2011.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book's primary reference [320] (ENTCS 59(4), 2001) could not be read: ScienceDirect refuses automated downloads and the Wayback Machine and the MGS site (mgs.spatial-computing.org, certificate for lacl.fr only) were unreachable. The language was reconstructed from the 2002 UMC paper by the same authors, which restates it, Cohen 2003 and the 2011 SCS paper. The Fundamenta Informaticae 2002 and BioSystems 2003 [321] papers were not obtained either; Anabaena and cAMP models are not offered.
- Scope: set, bag, seq, the GBFs grid <north, east> and hexagon <east, north, northeast; east + north = northeast>, and tube (a bag of sequences, the only nesting needed by the published programs). Maps/records, free-group GBFs (trees), Delaunay and cellular-complex collections, rule attributes (probabilities, non-consuming rules), sub-collection patterns on sets ((x / x < 3)+ as S) and user-defined functions are not implemented.
- Pattern language: variables, integer/string/boolean literals, <undef>, ',' and |d> (d any signed sum of generators, as in Turn's |Y-X>), per-element guards x/guard (all guards are conjoined and evaluated with every variable bound, which is equivalent because guards have no side effects) and X+, X*, x+ as X on tube strands. Guards and replacements use a small expression language (or, and, not, =, ==, !=, <, <=, >, >=, +, -, *, mod, div, size) parsed by Chemart, never Python eval. Replacements accept the book form y, x, Cohen's [x] and y::x::empty_seq (Cohen prints 'empty seq'), and () for nothing.
- Positional collections: species are value@position and a replacement must have as many elements as the pattern (for GBFs this is the language's rule; Cohen 2003 sec. 3). Leibnizian length changes on seq (x => x, x) are not offered because they renumber positions; strand rewriting with length changes is available in a tube. <undef> is a species per site (undef@r_c) so that growth into an empty site is a reaction. Grids are finite windows of the infinite GBF (bounded or torus): the published figures show only part of the infinite domain.
- Tube rules match whole strands: the published program applies Restriction to a DNA string and keeps hd of the result, a TUBE holding the whole cut strand, so X+ and Y+ must cover the strand; strands with several recognition sites give one reaction per cut. A strand id is its atoms joined (CCCG) when all are single characters, otherwise (atom,atom).
- Priorities: a lower-priority rule is not applied to a sub-collection (same positions and values, or same strand) that a higher-priority rule matches; this is exact for Void after EcoRI and for default rules x => [x]. Elastic rewritings (replacement equal to the match) produce no reaction and are skipped in runs.
- Closure semantics: species-level closure like chemart.expand (every rule on every tuple of known species lying on a path of the topology); a reaction may therefore combine site values that never coexist in a single run. Sets are treated as bags whose replacement is deduplicated; the union with an equal value already present elsewhere cannot be expressed as stoichiometry.
- Runs: the MGS default T[iter=fixpoint] with the maximal-parallel strategy (Cohen 2003 sec. 2; SCS 2011 sec. 4.1) or an asynchronous alternative, stopped when no rule changes the collection or after 10000 steps (Turn never reaches a fixpoint). The run is recorded in extras.run and does not change the network; the stochastic Gillespie strategy of MGS is not offered because no rates are published.
- Default inputs the sources do not give: the book's sort has no input ([4, 2, 5, 1, 3] is used); the sieve runs on 2..30 (Cohen: 'the integers between 2 and n'); Eden starts from one true cell in the middle of a 7x7 grid. The restriction-enzymes strands, the bead rows 3, 2, 4, 2 and the Turn cross are the published ones.
- v1 parameters collection_type (enum) and transformations (callable) became collection plus program/rules. The v1 provides tag compartments was dropped: the only nesting implemented is a tube of strands, not membranes. 'catalysts' comes from rules that keep a matched element (sieve x + y -> x, Eden).
- No kinetics: MGS transformations have no rates, so reactions carry no rate law.

## Notes

Chemart can borrow its collection typing as the space abstraction. The module exports Machine (closure, step, run), Rule, Topology and evaluate for direct use; step() applies one transformation step so figure sequences can be reproduced.

---

*Specification: `catalog/chemistries/mgs.yaml` · generator: `chemart/chemistries/mgs.py` · tests: `tests/chemistries/test_mgs.py`*
