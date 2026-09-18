# Molecular Classifier System (MCS.bl)

`mcs-bl` · *Decraene, Mitchell & McMullin, 2006-2011; based on Holland's Broadcast Language (1975) and classifier systems*

*Also known as:* *MCS.b*, *Molecular Classifier System - broadcast language*

Every molecule can be the enzyme or the substrate. Strings over a small alphabet of data and control symbols act on each other by pattern matching: one string's units match a region of another and rewrite it, with both reactants surviving. Among all 4096 strings of length four exactly one is a self-replicase that copies any string - and, tellingly, it enjoys no selective advantage over non-enzymatic strings, which is the model's point about how hard it is to get selection started.

| | |
|---|---|
| **family** | bio-inspired |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.1.1 |
| **refs** | [221], [222], [224], [225], [452], [570], [395], [396], doi:10.1142/S0219525911002913 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): Broadcast devices: strings over {0, 1, *, :, #, $, %, '} (thesis glyphs diamond, reversed triangle, triangle, prime for # $ % '). Species id m_<one letter per symbol: 0 1, * S, : C, # H, $ D, % P, ' Q>, structure the raw string. A device holds zero or more broadcast units '*condition:action' (null device: none).

**R — reactions** (implicit, arity 2): enzyme + substrate -> enzyme + substrate + product (both reactants catalytic); elastic if no unit of the enzyme matches the substrate

**A — reactor**: well-stirred-multiset, compartments
 · *dilution:* single reactor: grows to n_max, then each product displaces a random molecule other than the reactants (catalytic network equation with dilution flow, thesis eq. 5.1); cell model: cells divide at a criterion and a random other cell is removed

## What you get

```python
net = chemart.generate_network("mcs-bl", seed=1)
```

```
mcs-bl: 4 species, 6 reactions, status=complete
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: seed
```

First reactions:

```
m_SD0CD1 + m_SD0CD0 -> 2 m_SD0CD1 + m_SD0CD0  [mass-action k=1.0]
m_SD0CD1 + m_SD1CD0 -> m_SD0CD1 + m_SD1CD0 + m_SD1CD1  [mass-action k=1.0]
m_SD0CD0 + m_SD1CD0 -> m_SD0CD0 + 2 m_SD1CD0  [mass-action k=1.0]
m_SD1CD0 + m_SD0CD1 -> m_SD1CD0 + m_SD0CD1 + m_SD0CD0  [mass-action k=1.0]
m_SD1CD0 + m_SD1CD1 -> 2 m_SD1CD0 + m_SD1CD1  [mass-action k=1.0]
m_SD1CD1 + m_SD0CD1 -> m_SD1CD1 + 2 m_SD0CD1  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `closure` | structural | closure: every reaction among the strings reachable from the seeds (status complete or truncated); soup: the thesis' single-reactor model (4.2.1) for `steps` collisions with mutation, returning the reactions that fired with counts (status observed) <br>one of `closure`, `soup` |
| `strings` | `list` | `['*$0:$1', '*$0:$0', '*$1:$0', '*$1:$1']` | structural | seed broadcast devices over {0,1,*,:,#,$,%,'} (thesis glyphs are accepted and translated). Default: cell type c0 of book table 11.2 (thesis table 7.1) <br>*range:* book table 11.2 cell types: c1 [*$0:$1, *$0:$0, *$#:$0, *$#:$1], c2 [*$0:$1, *$0:$0, *$#:$%0, *$#:$%1], c3 [*$0:$1, *$0:$0, *$%:$0, *$%:$1]; universal replicase [*$:$]; thesis table 5.2 replicases *$1:$1 .. *$0101:$0101 |
| `self_replication` | `bool` | `` | selection | allow reactions whose product is identical to both enzyme and substrate (s + s -> 2s + s); false makes them elastic <br>*range:* true in thesis ch. 5-6 and book table 11.1; false in the table 11.2 experiment (thesis 7.2, ACS 2011, CEC 2009) |
| `max_length` | `int` | `500` | structural | BDL_max: a reaction whose product is longer is not permitted (elastic); bounds the elongation catastrophe <br>`4` … `1000000` · *range:* thesis: BDL_max = 500 (ch. 6-7), 10^6 (5.5) |
| `max_species` | `int` | `200` | structural | closure only: species budget; elongating chemistries are cut off by it (status truncated) <br>`1` … `100000` |
| `initial_copies` | `int` | `10` | population | copies of each seed string (initial_state; initial soup) <br>`1` … `1000000` · *range:* thesis 7.2: 10 of each seed species; 5.3: 100 copies of the replicase |
| `n_random` | `int` | `` | population | extra seed strings drawn uniformly from the alphabet (each with initial_copies copies) <br>`0` … `1000000` · *range:* thesis 5.2: 100 random molecules; 5.3: 900 |
| `random_length` | `int` | `10` | population | length of the random seed strings <br>`1` … `1000000` · *range:* thesis 5.2-5.3: 10 |
| `n_max` | `int` | `1000` | population | soup only: reactor capacity; below it products are added, at it each product displaces a random molecule other than the reactants <br>`3` … `100000000` · *range:* thesis 5.2-5.3: 1000; cells in 7.2: 10^6 |
| `steps` | `int` | `10000` | population | soup only: number of collisions (elastic ones included) <br>`0` … `1000000000` · *range:* thesis 5.2: 5x10^6 collisions per run; 7.2: over 4x10^7 per cell per hour |
| `p_s` | `float` | `1e-05` | stochastic | soup only: per-symbol mutation probability of each product; a mutation is a flip to another symbol, an insertion of a random symbol after it, or a deletion, with equal probability <br>`0.0` … `1.0` · *range:* book fig. 11.3 / ACS 2011: 1e-5; thesis 5.2: 1e-3; D.1: 5e-5 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- universal replicase *$:$ copies any string; it has no selective advantage over non-enzymatic strings (thesis eqs. 5.3-5.6: zero expected growth under dilution)
- of the 8^4 = 4,096 strings of length 4 only *$:$ is a self-replicase (thesis 5.2, which misprints the count as 4^8 = 65,536)
- elongation catastrophe: the mutant *$0101:$00101 of *$0101:$0101 makes ever longer strings *$0101:$0...0101 (thesis 5.5, ALife 2008)
- succession of dominant cell types c0 -> c1 -> c2 -> c3 under selection for s1 (book table 11.2); c1-c3 are phenotypically equivalent and have 9 reactions against 6 in c0, all closed
- selection pressure can be aimed at a chosen molecular species by changing the division criterion
- hand-designed collectively autocatalytic sets *$0:$1 + *$1:$0 close into c0 (thesis appendix D.2)

## Sources

- Decraene, J. (2009). Autocatalytic Closure and the Evolution of Cellular Information Processing Networks. PhD thesis, Dublin City University: MCS reactor and mutation (4.2), broadcast language syntax and semantics with examples r1-r20 and table 4.2 (4.3), self-replicases (tables 5.1-5.2, eqs. 5.1-5.6), elongation catastrophe (5.5), cell model (6.4.2), table 7.1-7.2 experiment (7.2), appendices B and D.2. https://doras.dcu.ie/14810/1/jd-thesis-09-08-09.pdf
- Decraene, J. & McMullin, B. (2011). The evolution of complexity in self-maintaining cellular information processing networks. Advances in Complex Systems 14(1):55-75 (book ref [221]); preprint: alphabet and table 2 (operations), mutation (2.2), multilevel model (2.3, fig. 2), table 3 (= book table 11.2), table 4. http://doras.dcu.ie/16292/1/decraene-ACS-09.pdf
- Decraene, J., Mitchell, G. G. & McMullin, B. (2009). Crosstalk and the cooperation of collectively autocatalytic reaction networks. IEEE CEC 2009, 2249-2256 (book ref [225]): table I, random choice among satisfied units, reactor algorithm (fig. 5). http://doras.dcu.ie/4609/1/decraene-cec-09.pdf
- Decraene, J., Mitchell, G. G. & McMullin, B. (2008). Unexpected evolutionary dynamics in a string based artificial chemistry. ALife XI, 158-165 (book ref [224]): reactor, BDL_max rule, elongation catastrophe. http://doras.dcu.ie/4589/1/decraene-alife-final-08.pdf
- Decraene, J. (2006). The Holland Broadcast Language. Technical report ALL-06-01, DCU: specification and C++ source (BUnit.cpp: setRegexForm, setCompactForm, getSymbolValueFromSignal, reactWith). http://doras.dcu.ie/4590/1/decraene-ALL-06-01.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Symbols: the book's extracted text prints the alphabet as {3, $, %, } and table 11.2 uses 3; comparing with thesis table 7.2 (glyphs), book 3/# is the diamond (single wildcard, any suffix when last), $ the reversed triangle (prefix/suffix wildcard that transposes), % the triangle (single wildcard that transposes) and the lost fourth symbol the quote. The book's table 11.1 also lost the quotes of the activation (*$1:'*$) and inhibition (*'*0$:0$) enzymes.
- Matching semantics come from the thesis prose and examples plus the C++ source of the authors' earlier BL implementation (ALL-06-01); the MCS.bl C++ package (esignet.net/dokumente/upload/WP13...) is not archived. From the source: $ matches one or more symbols (regex .+), so it never matches an empty string; a symbol is quoted iff the preceding symbol is a quote; unquoted quotes are dropped; the $ value is positional (the other condition symbols take one character each); $ in the middle of a condition is ignored, and with $ at both ends the last is ignored; only the first % counts. The source's length arithmetic counts quote characters in the condition (a bug for quoted symbols after a leading $); here the matched string is computed from the pattern symbols.
- # as the last condition symbol matches any non-empty suffix (thesis r13-r15); after a leading $ it can only take one symbol, since $ absorbs the rest (this is the thesis' statement that # and % act the same in *$#:$0 and *$%:$0, and gives c1-c3 nine reactions each).
- Book table 11.1 / thesis table 4.2 / CEC table I / ACS table 2 activation row is inconsistent: *$1:'*$ on 0:1 gives *0: (not an active device) by the same rule the cleavage and concatenation rows use (*$1:$ on *0:1 gives *0:, so $ excludes the matched 1). The intended activation is reproduced by *$:'*$ on 0:1 -> *0:1 (tested).
- Several units: CEC 2009 (II-A) picks one satisfied unit at random. The closure lists every possible product; a reaction's rate weight is the fraction of the enzyme's satisfied units that give it.
- Elastic collisions: no unit matches; the product is empty; the product is longer than max_length (ALife 2008: such reactions 'were simply not permitted'; truncation was a later, abandoned variant); or, with self_replication false, the product equals both reactants.
- Rates: thesis eq. 5.1-5.2 (catalytic network equation, alpha = 1 per ordered enzyme-substrate pair that yields the product) with dilution flow = outflow constant-total. Each multiset reaction {a, b} -> {a, b, p} gets mass-action k = sum over the orders (a, b), (b, a) of the fraction of satisfied units giving p (1 or 2 for single-unit enzymes; one order when a = b). Soup reactions made by a mutation have no structural k and get rate None.
- Soup: two distinct molecules drawn at random, the first the enzyme (thesis 4.2.1, 5.3); below n_max the product is added, at n_max it replaces a random molecule other than the two reactants (thesis; ALife 2008 does not exclude them). Mutation (thesis 4.2.2) is applied to the product once per symbol with probability p_s, sampled as a binomial number of distinct positions; self-replication is filtered before mutation and length after. The thesis' spontaneous mutation (a fraction r_mut of the population every x timesteps, x not given) is not modelled; the table 11.2 experiment used r_mut = 0.
- Multilevel cell selection is not part of one network and is not modelled: v1 params n_cells, cell_capacity and division_target are dropped (they size a population of reactors, not the chemistry); the model is described in A.notes. Division variants: thesis 6.4.2 divides when a cell is full and overwrites a random target cell; thesis 7.2 / ACS 2011 divide at n_target copies of s_T and remove a random other cell.
- Book prose erratum (11.1.1): 'c1 consists of s1, s2, s3, s4, whereas c2 consists of s1, s2, s5, s6' is shifted by one; table 11.2 (and thesis table 7.2) give c0 = {s1..s4}, c1 = {s1, s2, s5, s6}. The book's 'one molecule (rather than the full number) was selected as cell division threshold' means one target species. The reaction counts 6 (c0) and 9 (c1) are from thesis 7.2.3 and ACS 2011 (three additional reactions).
- Thesis errata found while testing: appendix D.2 lists s4 + s2 -> s4 for s4 = *$0:$0, s2 = *$1:$0; the rule gives s2 (replication of the substrate), and the other 15 entries of that table agree; D.2 prints the species without the leading * (s1 = $0:$1) and s3 as $1:$ for *$1:$1. Section 5.5 prints s_R4 = *$0101:*$0101 with an unquoted * in the action, which would leave no unit; table 5.2 and the fig. 5.4 caption give *$0101:$0101, whose reactions reproduce 5.5 exactly (tested). Section 5.2 (and ALife 2008) count '4^8 (65,536)' molecules of length 4; with 8 symbols there are 8^4 = 4,096, of which only *$:$ self-replicates (tested).
- Open doubt: of the 15 spontaneous self-replicases of thesis table 5.1 (ALife 2008 table 2), 9 self-replicate under the published semantics (tested), but 6 have $$ in the action (e.g. *$*$:$%$%%), which by thesis example r17 doubles the substrate instead; the thesis nevertheless calls all 14 non-*$0:$0 entries equivalent to *$:$. The MCS.bl source that would settle this is not available.

## Notes

A reflexive string-rewriting chemistry: every molecule can be both enzyme and substrate. Chemart provides the molecular level (closure or single reactor); the cell level of the book's figures (31 concurrent cells racing in wall-clock time) is not a reaction network.

---

*Specification: `catalog/chemistries/mcs-bl.yaml` · generator: `chemart/chemistries/mcs_bl.py` · tests: `tests/chemistries/test_mcs_bl.py`*
