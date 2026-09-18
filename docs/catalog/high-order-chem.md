# High-order chemistry (rules as molecules)

`high-order-chem` · *Yamamoto, 2014 (PyCellChemistry); Banzhaf & Yamamoto, 2015*

*Also known as:* *HighOrderChem*, *NumberChemHO*

A chemistry whose rules are themselves molecules. The vessel holds two multisets, one of data and one of rules; each step draws a rule, counts its binding sites, draws that many data molecules, applies it, and returns everything - rule included - to the vessel. Because a rule is a molecule, its concentration matters: rules compete for substrate and can be added or removed at run time. The natural next step, rules that rewrite rules, is what the design points at without yet doing.

| | |
|---|---|
| **family** | core |
| **kind** | framework |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | appendix: Writing Your Own Artificial Chemistry in Python, 'A High-Order Chemistry' (figure 3, divrule); module list (HighOrderChem.py) |
| **refs** | https://github.com/laryamamoto/PyCellChemistry |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): two multisets: data molecules (JSON integers or lists; species n<k> or compact JSON text such as [0,3,1,2]) and rule molecules (species rule:<name>, structure the rule text, e.g. divrule(m1, m2))

**R — reactions** (implicit, arity variable - the number of binding sites of the drawn rule, plus the rule molecule itself): rule + m1 + ... + mk -> rule + products(m1, ..., mk), k = binding sites of the rule; e.g. rule:divrule + n12 + n3 -> rule:divrule + n4 + n3

**A — reactor**: well-stirred-multiset
 · *dilution:* none - products replace educts; the data population changes only if a rule returns more or fewer molecules than it binds

## What you get

```python
net = chemart.generate_network("high-order-chem", seed=1)
```

```
high-order-chem: 198 species, 154 reactions, status=observed
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: analysis, final_state, rules
```

First reactions:

```
rule:divrule + n261 + n87 -> rule:divrule + n3 + n87  (x1)
rule:divrule + n771 + n3 -> rule:divrule + n257 + n3  (x1)
rule:divrule + n981 + n3 -> rule:divrule + n327 + n3  (x1)
rule:divrule + n410 + n820 -> rule:divrule + n410 + n2  (x1)
rule:divrule + n2 + n866 -> rule:divrule + n2 + n433  (x1)
rule:divrule + n369 + n3 -> rule:divrule + n123 + n3  (x1)
rule:divrule + n512 + n64 -> rule:divrule + n8 + n64  (x1)
rule:divrule + n918 + n2 -> rule:divrule + n459 + n2  (x1)
… and 146 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `soup` | structural | soup: the book's HighOrderChem.iterate loop, observed effective reactions with counts; closure: every reaction reachable from the distinct data molecules under deterministic rules (Chemart addition) <br>one of `soup`, `closure` |
| `rules` | `dict` | `{'divrule': 4}` | structural | the rule multiset {rule: multiplicity}; this IS the chemistry. A rule is a built-in name (divrule, exchangeMachine, cutMachine, invertMachine, recombinationMachine), optionally written as its rule molecule 'divrule(m1, m2)', or an expression rule '[name:] x, y -> expr, ... [if condition]' (integers, + - * / % comparisons and or not, min max abs; see the module docstring). Rules are selected with probability proportional to multiplicity <br>*range:* NumberChemHO.py: {divrule: 4}; MolecularTSP.py: {exchangeMachine: 100, cutMachine: 100, invertMachine: 100, recombinationMachine: 1} |
| `data` | `list` | `` | structural | explicit initial data multiset as a list of JSON integers or lists; overrides the random initialisation selected by init |
| `init` | `enum` | `numbers` | structural | random data multiset when data is empty: numbers draws M integers uniformly from [minn, maxn] (NumberChemHO.py); tours draws M random tours on the ring of `cities` cities (MolecularTSP.py) <br>one of `numbers`, `tours` |
| `M` | `int` | `100` | population | number of data molecules drawn when data is empty <br>`1` … `100000` · *range:* NumberChemHO.py popsize 100; MolecularTSP.py popsize 9 |
| `minn` | `int` | `2` | structural | init numbers: lower bound of the uniform draw (inclusive) <br>≥ `1` |
| `maxn` | `int` | `1000` | structural | init numbers: upper bound of the uniform draw (inclusive) <br>≥ `1` |
| `cities` | `int` | `10` | structural | number of cities on the ring instance that the tour machines (and init tours) use; tours are permutations of 0..cities-1 <br>`3` … `1000` · *range:* MolecularTSP.py default 10 |
| `iterations` | `int` | `10000` | population | soup only: number of iterations of the algorithm (one rule drawn per iteration, elastic and idle ones included) <br>`0` … `10000000` · *range:* NumberChemHO.run: 10000; MolecularTSP.py: up to 1000 generations of ceil(M * 100 / \|rules\|) iterations |
| `max_species` | `int` | `1000` | structural | closure only: budget of data species; status truncated when exceeded <br>≥ `1` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the divrule rule set reproduces the prime-number chemistry: primes are never consumed, composites are divided, and the prime fraction grows towards 1 (book appendix: 'the results should be the same as the original NumberChem implementation')
- rule molecules are catalysts: every reaction has its rule on both sides and the rule multiset stays constant
- rules compete for substrate: each iteration one rule is chosen with probability proportional to its multiplicity
- changing the rule multiset changes the chemistry without changing the algorithm (e.g. the molecular TSP machines on tour molecules)

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries. MIT Press. Appendix, 'A High-Order Chemistry' (figure 3 and the divrule listing).
- PyCellChemistry src/HighOrderChem.py, src/NumberChemHO.py, src/MolecularTSP.py and src/artchem/Multiset.py (Yamamoto, 2014). https://github.com/laryamamoto/PyCellChemistry

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Security: the reference executes rule strings with exec. Chemart never executes a string: a rule is looked up by name in a registry of built-in functions ported from the reference examples (divrule from NumberChemHO.py and the book; exchangeMachine, cutMachine, invertMachine, recombinationMachine from MolecularTSP.py), or parsed and interpreted as an expression rule in a small documented language. Unknown names are rejected with ValueError. Subclass methods (the book's self.fold(m)) are replaced by the registry; a leading self. is accepted and ignored.
- Binding sites: for named rules nsites is computed from the rule molecule exactly as figure 3 does (len(bsite.split(','))), and a user-written rule molecule must have as many sites as the function takes (the reference would crash on a mismatch). For expression rules nsites is the number of binding-site variables.
- Multisets are lists; a random draw removes a uniformly chosen molecule, the same distribution as Multiset.expelrnd (each molecule with probability 1/total). The rule is drawn and reinjected every iteration, so the rule multiset never changes and rules are catalysts in the network (the rule species appears on both sides).
- When fewer data molecules than binding sites are present, the rule is reinjected and nothing else happens (figure 3 lines 17 and 26); such iterations count towards `iterations` and are reported as extras.analysis.idle_draws.
- divrule is written with Python 2 integer division in the book (m1 / m2 on ints); it is floor division (//), exact here because the quotient is only taken when it divides. Built-in rules collide elastically (return their educts) on molecules outside their domain (non-integers or 0 for divrule, non-tours for the machines), where the reference would raise.
- Default rules {divrule: 4}, M = 100 from [2, 1000] and 10000 iterations are NumberChemHO.py's; the book only says 'at least one molecule of divrule(m1, m2)'. With a single rule type the multiplicity does not change the dynamics.
- The observed network keeps the effective reactions (HighOrderChem.is_effective: product multiset differs from educt multiset) with firing counts; initial_state and extras.final_state include the rule molecules; extras.analysis records rule_draws per rule, and prime_fraction per generation of M iterations when the only rule is divrule on integers (NumberChemHO.nprimes).
- Tour machines use MolecularTSP.py's ring topology (TSPgraph ring=True: gridsize 2N, radius N, fully meshed roads) and its fitness with penalties; a tour molecule is the list of cities rather than the string 'fitness [tour]', since the fitness is a function of the tour. The reference's %g rounding of the stored fitness is kept (the drawn molecules' fitness is compared after '%g' formatting, the new tour's exactly). The random topology is not offered: the dedicated molecular-tsp chemistry covers the full model. List molecules must be permutations of 0..cities-1 when a machine is in the rule set.
- method closure is a Chemart addition (chemart.expand.expand, ordered reactants, each rule producing its own reaction); it is only defined for deterministic rules, so the stochastic tour machines are rejected.
- The book's caveat is kept: this is not yet truly high-order, since rules never take rules as educts and cannot rewrite rules (no typed binding sites).
- v1 parameters rules (callable) and seed (matrix) are replaced by the JSON rule multiset `rules` and the data multiset `data`/`init`; the v1 name `seed` is the generate_network seed.

## Notes

Not a chemistry but a substrate for constructive chemistries, and the closest PyCellChemistry analogue of Chemart's own react/expand/soup split. The book notes that making it truly high-order would need typed binding sites so that rules can take other rules as educts, and rule sets that keep rewritten rules functional.

---

*Specification: `catalog/chemistries/high-order-chem.yaml` · generator: `chemart/chemistries/high_order_chem.py` · tests: `tests/chemistries/test_high_order_chem.py`*
