# Gamma / gamma-calculus

`gamma` · *Banatre & Le Metayer, 1986 (INRIA RR-566), 1990, 1993; gamma-calculus: Banatre, Fradet & Radenac, 2004*

*Also known as:* *general abstract multiset manipulation*, *chemical reaction model*, *programming by multiset transformation*

The language that made multiset rewriting a programming paradigm. A program is a set of rules of the form 'pick n elements satisfying this condition and replace them with these', and the elements are just values in a bag. Computing the maximum of a multiset is one rule: take any two, discard the smaller. There are no rates and no time, only the stable state - and that state is the answer regardless of the order the comparisons happened in.

| | |
|---|---|
| **family** | rewriting |
| **kind** | formalism |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 9.2 |
| **refs** | [55], [56], [57], [58], doi:10.1016/0167-6423(90)90044-E, doi:10.1145/151233.151242, doi:10.1007/3-540-45523-X_2 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): elements of a multiset: integers, or tuples of elements such as (index, value) pairs or (position, value, sum) triples. Species id n<k> for the integer k and the canonical text (1,5) for a tuple; structure is the canonical text

**R — reactions** (implicit, arity variable): name: x1, ..., xn -> A(x1, ..., xn) <= R(x1, ..., xn): n elements satisfying the reaction condition R are replaced by the elements produced by the action A (any number, possibly none). Example (book eq. 9.7): max: x, y -> y <= x <= y. Chemart text: 'max: x, y -> y if x <= y'

**A — reactor**: well-stirred-multiset
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("gamma", seed=1)
```

```
gamma: 4 species, 10 reactions, status=complete
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: analysis, program, reaction_rules, stages
```

First reactions:

```
2 n1 -> n1
n1 + n3 -> n3
n1 + n4 -> n4
n1 + n8 -> n8
2 n3 -> n3
n3 + n4 -> n4
n3 + n8 -> n8
2 n4 -> n4
… and 2 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `program` | `enum` | `max` | structural | published program: max (book eq. 9.7; Gamma15 sec. 1), sort (exchange sort of (index, value) pairs, Gamma15 sec. 1), primes (rem o iota from {(2, n)}, Gamma15 sec. 2.1), fibonacci (add o dec1 from {n}, Gamma15 sec. 2.1), max-segment-sum (maxg o maxl on (i, x, x) triples, Gamma15 sec. 2.1), majority (maj, RULE04 fig. 1), largest-prime (primes then max on 2..10, the gamma-cn example of the York abstract), custom (rules) <br>one of `max`, `sort`, `primes`, `fibonacci`, `max-segment-sum`, `majority`, `largest-prime`, `custom` |
| `rules` | `str` | `` | structural | custom only: program text. Patterns are variables, integers or tuples; expressions use integers, tuples, + - * / (floor) mod, comparisons, and or not, true false, multiple(x, y), min, max, abs; the book's arrows and <= >= != = symbols are accepted. Never evaluated as Python <br>*range:* one reaction per line, 'name: patterns -> expressions if condition'; a line 'then' starts the next sequentially composed stage; e.g. 'gcd: x, y -> x - y, y if x > y' |
| `multiset` | `list` | `` | population | initial multiset: integers and lists (tuples of at least two elements); empty uses the named program's example input (required for custom) <br>*range:* max [3, 8, 1, 8, 4]; sort [[1, 8], [2, 3], [3, 6], [4, 1], [5, 9]]; primes [[2, 20]]; fibonacci [5]; max-segment-sum [[1, 3, 3], [2, -4, -4], [3, 5, 5], [4, -1, -1], [5, 2, 2]]; majority [1, 2, 1, 3, 1, 1, 2]; largest-prime [2, ..., 10] |
| `max_species` | `int` | `200` | structural | closure budget: species of the network <br>`1` … `100000` · *range:* the fibonacci add stage and other unbounded reducers have infinite closures and are truncated |
| `max_steps` | `int` | `100000` | population | reactions fired by the sampled execution before it is reported as not stable (for programs that do not terminate) <br>`0` … `100000000` |
| `max_states` | `int` | `5000` | population | budget of multisets visited when enumerating every stable multiset reachable under any reaction order (analysis.results) <br>`0` … `10000000` · *range:* 0 disables; primes from (2, 20) exceeds 5000 states in its iota stage |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- max (book eq. 9.7): the stable multiset is the single maximum, whatever the order of comparisons
- sort (Gamma15): exchanging ill-ordered (index, value) pairs ends with values increasing with the index; the set of indexes and the multiset of values are invariant
- primes (Gamma15 rem o iota): from (2, n) iota produces 2..n and rem leaves exactly the primes <= n
- fibonacci (Gamma15 add o dec1): the multiset expands into ones and shrinks to fib(n) (fib(0) = fib(1) = 1)
- maximum segment sum (Gamma15): the stable multiset holds the triples whose s is the maximum segment sum
- majority (RULE04): removing pairs of distinct elements leaves only the majority element, with an order-dependent multiplicity
- re-injecting elements into a stable multiset restarts the reactions and yields the result for the enlarged input (book 9.2)

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 9.2: reaction condition/action pairs, max (eq. 9.7), sort (eq. 9.8), primes (eq. 9.9), inert result, re-injection, gamma-calculus.
- Banatre, J.-P., Fradet, P. & Le Metayer, D. (2001). Gamma and the chemical reaction model: fifteen years after. Multiset Processing, LNCS 2235:17-44 (book [56]; 'Gamma15'). Sec. 1: max and exchange sort; sec. 2.1: primes = rem o iota, fibonacci = add o dec1, maximum segment sum; sec. 2.2: tropes; sec. 3.1: invariant/variant proof of sort and derivation of rem. https://pop-art.inrialpes.fr/~fradet/PDFs/Gamma15.pdf
- Banatre, J.-P., Fradet, P. & Radenac, Y. (2004). Principles of chemical programming. RULE'04, ENTCS 124(1):133-147, 2005 ('RULE04'). Fig. 1: max, primes (multiple), maj; secs. 2-4: gamma0, gamma-c, gamma-n, gamma-cn calculi and the relation of Gamma to gamma-cn. https://pop-art.inrialpes.fr/~fradet/PDFs/RULE04.pdf
- Banatre, J.-P., Fradet, P. & Radenac, Y. Higher-order chemical model of computation (workshop abstract, York). Gamma examples and the gamma-cn molecule computing the largest prime below 10. https://www.cs.york.ac.uk/nature/workshop/papers/BanatreFradetRadenac.pdf
- Banatre, J.-P., Fradet, P. & Radenac, Y. (2003). Higher-order chemistry. MolCoNet Workshop on Membrane Computing, Tarragona: gamma-abstractions gamma(P)[C].M, gamma-conversion, inertness, maximum examples. https://pop-art.inrialpes.fr/~fradet/PDFs/WMC03.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Book eq. 9.8 prints the sort condition as (i >= j) and (x >= y); read literally it exchanges pairs that are already in increasing order and so sorts decreasingly, contradicting the book's prose (the sorted set starts from index 1, reactions only on pairs out of order). Gamma15 sec. 1 and its correctness proof (sec. 3.1, variant: i > j implies x >= y) give (i > j) and (x < y); Chemart uses that. The literal book rule can be written as custom rules.
- Book eq. 9.9 (prime: x, y -> y/x <= y mod x = 0) is the book's division chemistry (prime-number-chemistry in this catalog): its prose keeps x and replaces y by y/x, while the printed action drops x. Gamma's own primes program is rem: x, y -> y <= multiple(x, y) (Gamma15 sec. 2.1, RULE04 fig. 1), which removes multiples of other elements; Chemart uses it, composed with iota as in Gamma15.
- The 1990 Science of Computer Programming paper [58] and the 1993 CACM paper were not accessible; the programs come from Gamma15 (the revised survey by the same authors), RULE04 and the York abstract. gcd and sum are not offered as named programs because no accessible source prints them; custom rules express them (e.g. gcd: x, y -> x - y, y if x > y).
- Garbled extraction: pdftotext of Gamma15 lost the <= symbol and the primes. max reads x <= y (the text says the elements are replaced by y); '-' printed as ';' in dec1 (x -> x - 1, x - 2) as in fib(x ; 1); maxl is (i,x,s),(i',x',s') -> (i,x,s),(i',x',s+x') <= i' = i + 1 and s + x' > s'; maxg is (i,x,s),(i',x',s') -> (i',x',s') <= s' > s, the only reading that keeps the larger sum. Chemart names i', x', s' as j, y, t.
- Example inputs are Chemart choices (the sources give programs, not data): max [3, 8, 1, 8, 4] (with a duplicated maximum), sort values 8 3 6 1 9, primes n = 20, fibonacci n = 5, max-segment-sum 3 -4 5 -1 2 (answer 6), majority [1, 2, 1, 3, 1, 1, 2], largest-prime 2..10 as in the York example.
- Reaction language: a small parsed expression language over integers and tuples (no floats, strings or symbols), never Python eval. '/' is floor division (Gamma15's [ ] integer part is plain grouping); multiple(x, y) is x mod y = 0 with y != 0; a repeated pattern variable requires equal elements (the gamma-calculus match unifies variables). An ill-typed expression (tuple arithmetic, comparing a tuple with <, division by zero) makes the reaction inapplicable to that tuple rather than an error.
- A reaction whose products equal its reactants (as a multiset) does not count as a reaction: it cannot change the multiset, so a multiset where only such reactions apply is stable. This matches chemart.expand's elastic collisions.
- Sequential composition is written with a 'then' line. The closure of P2 o P1 is the closure of P1 from the initial elements, then the closure of P2 seeded with every species found so far: a superset of the elements P2 can actually meet, since P2 starts from P1's stable results. Reactions keep the rule name(s) that produce them in extras.reaction_rules; a reaction produced by rules of two stages carries both names.
- Closure semantics: expand over ordered tuples of species, assuming enough copies (x, x -> x appears for a single copy of x), alternatives=True because several reactions can apply to the same tuple. Unbounded programs (fibonacci's add, custom reducers) are truncated by max_species.
- Execution: Gamma leaves the choice of tuple unspecified. The sampled run fires one reaction per step, chosen with probability proportional to the number of ordered tuples of distinct element occurrences it applies to; parallel firing of disjoint tuples is equivalent to some such interleaving. analysis.results enumerates every stable multiset reachable under all orders (depth-first over multisets, per stage, budget max_states), which is how Chemart shows that a result is order independent (deterministic true) or not (majority). A cycle without exit yields no stable multiset and is not reported separately.
- Higher-order Gamma / gamma-calculus (book [55]; RULE04, WMC03) is not generated: reactions as consumable molecules, nested solutions and inertness of sub-solutions need a term rewriting engine that the closure/soup model does not capture faithfully. Its first-order core is available: RULE04 sec. 4.2.1 shows a Gamma program is a gamma-cn molecule, and the York largest-prime example (a solution that must become inert before max is introduced) is sequential composition, offered as program largest-prime. The v1 higher_order flag is dropped.
- v1 fixes: constructive is true (actions create new elements), the reactor is a well-stirred multiset with nondeterministic choice rather than maximally parallel (parallelism is allowed, not required), the v1 callable rules and matrix initial_multiset became program/rules text and a JSON list. No rates: Gamma has no kinetics. analysis keys: final_multiset, final_values (JSON values) and stable/steps of the sampled run, stage_multisets (after each stage), results and deterministic (null when the state budget is exceeded).

## Notes

Purely qualitative: no rate constants, no time. Useful to Chemart as a front-end language for writing implicit chemistries, and the ancestor of CHAM, HOCL, P systems and Fraglets. The module exports Rule, parse_program, applicable, run, results, closure and PROGRAMS.

---

*Specification: `catalog/chemistries/gamma.yaml` · generator: `chemart/chemistries/gamma.py` · tests: `tests/chemistries/test_gamma.py`*
