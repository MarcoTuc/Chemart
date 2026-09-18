# Prime number (number-division) chemistry

`prime-number-chemistry` · *Banzhaf, Dittrich & Rauhe, 1996*

*Also known as:* *NumberChem*, *division chemistry*, *number-division chemistry*

The example that shows why a chemistry cannot always be a finite matrix. Molecules are integers; when a smaller number divides a larger one, the larger is replaced by the quotient while the divisor survives as a catalyst. Nothing ever consumes a prime, so primes accumulate as the inert residue of the reaction - the sieve of Eratosthenes as a reactor. The species set is unbounded and only ever discovered by running the rule.

| | |
|---|---|
| **family** | core |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 1 (eqs. 1.4-1.5), 2.5.2, appendix (NumberChem.py, NumberChemHO divrule) |
| **refs** | [72], [113], https://www.cs.mun.ca/~banzhaf/papers/nanotechnology7.pdf |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): natural numbers >= 2; species id n<k>, structure the decimal integer k

**R — reactions** (implicit, arity 2): s1 + s2 -> s1 + s2/s1 if s1 < s2 and s1 | s2; otherwise elastic (book eq. 2.43)

**A — reactor**: well-stirred-multiset
 · *dilution:* none - the divisor acts as a catalyst, so the molecule count stays M

## What you get

```python
net = chemart.generate_network("prime-number-chemistry", seed=1)
```

```
prime-number-chemistry: 201 species, 159 reactions, status=observed
provides: catalysts, initial-state, stoichiometry, topology
seed: 1
extras: analysis, final_state, primes
```

First reactions:

```
n261 + n29 -> n29 + n9  (x1)
n162 + n324 -> n162 + n2  (x1)
n530 + n2 -> n2 + n265  (x1)
n2 + n364 -> n2 + n182  (x1)
n250 + n750 -> n250 + n3  (x1)
n2 + n56 -> n2 + n28  (x1)
n756 + n2 -> n2 + n378  (x1)
n162 + n2 -> n2 + n81  (x1)
… and 151 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `soup` | structural | soup: the book's stochastic run, observed reactions with counts; closure: every reaction reachable from the distinct seed numbers (Chemart addition) <br>one of `soup`, `closure` |
| `M` | `int` | `100` | population | number of molecules drawn for the initial soup (ignored when numbers is given) <br>`2` … `100000` · *range:* book and paper: 100; paper fig. 6 scans 10-200 |
| `minn` | `int` | `2` | structural | lower bound of the uniform initialisation interval (inclusive) <br>≥ `2` |
| `maxn` | `int` | `1000` | structural | upper bound of the initialisation interval (inclusive); maxn >> M makes the run constructive <br>≥ `2` · *range:* book 2.5.2 and appendix: 1000; paper [72] fig. 5: 10000 |
| `numbers` | `list` | `` | structural | explicit initial multiset of integers >= 2; overrides the random draw of M numbers from [minn, maxn] <br>*range:* book fig. 2.7: each of 2..101 once |
| `iterations` | `int` | `10000` | population | soup only: number of collisions, elastic ones included (M collisions = one generation) <br>`0` … `10000000` · *range:* appendix: 10000; figs. 2.8-2.9: 20000; paper [72]: 700 generations = 700 M |
| `max_species` | `int` | `1000` | structural | closure only: species budget; the closure is finite, so it is truncated only when this is exceeded <br>≥ `1` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- primes emerge as the non-reactive set: no reaction consumes a prime, so the prime count never decreases
- the prime concentration grows sigmoidally, reaching 1.0 when M is large enough (paper [72] fig. 5: M = 100, maxn = 10000)
- initialising with each of 2..101 produces no new number types; initialising 100 numbers from [2, 1000] does (book 2.5.2, figs. 2.7-2.9)
- no prime factorisation - divisors are catalysts, so copies are not produced and only a limited number of copies of each prime remains
- phase transition in soup size: small soups (M < 40) end with non-primes left, large ones (M > 100) end all-prime, with strong run-to-run fluctuations in between (paper [72] fig. 6, 700 generations, 30 runs per M)

## Sources

- Banzhaf, W., Dittrich, P. & Rauhe, H. (1996). Emergent computation by catalytic reactions. Nanotechnology 7:307-314. Section 4.4, eq. 18, reactor algorithm II, figs. 5-6. https://www.cs.mun.ca/~banzhaf/papers/nanotechnology7.pdf
- PyCellChemistry src/NumberChem.py and src/artchem/Multiset.py (Yamamoto, 2013). https://github.com/laryamamoto/PyCellChemistry

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Collision algorithm follows the appendix code (NumberChem.py): two molecules are removed uniformly at random without replacement (Multiset.expelrnd twice, so they are distinct molecules), sorted, the larger replaced by larger/smaller if it is strictly larger and divisible, and both reinjected. chemart.soup.soup with dilution none does exactly this; the population stays M.
- The three descriptions agree on what is replaced: book eq. 2.43 and the 2.5.2 pseudocode (whose two branches are the two orderings of the draw), NumberChem.py (only m2 is overwritten) and paper [72] eq. 18 with reactor algorithm II (only s1 is replaced by s3) all replace only the dividend; the divisor is a catalyst. No source consumes both reactants. The 2.5.2 pseudocode's draw(P) does not say whether the same molecule can be drawn twice; the appendix code resolves it (no).
- Initial draw: M independent uniform integers in [minn, maxn], both ends inclusive (np.random.randint(minn, maxn+1) in NumberChem.py), so duplicates are possible. `numbers` gives an explicit multiset instead (book fig. 2.7 initialises with each of 2..101 once).
- iterations counts every collision, elastic ones included (NumberChem.py loop; book 2.6.1: a generation is M collisions). The soup output is status observed: species are every number present at the start or produced, reactions carry firing counts, initial_state is the drawn multiset, extras.final_state the final one, extras.analysis.prime_fraction the prime fraction after each generation.
- method closure is a Chemart addition, not in the book: chemart.expand.expand of the distinct seed numbers under the same rule. It is always finite (a quotient divides an existing number), so status is complete unless max_species cuts it.
- Species ids are n<k> rather than bare integers so that stoichiometry text is unambiguous (2 + 4 -> 2 + 2 reads n2 + n4 -> 2 n2); structure holds the integer. No rate constants: the book gives none.
- The v1 seed parameter is dropped (the seed is an argument of generate_network); the v1 phenomenon 'sigmoidal growth of prime concentration to 1.0' is qualified: a composite whose divisors are all absent from the soup is never divided, so the final prime fraction reaches 1.0 only for large enough soups (paper [72] fig. 6: M < 40 mostly ends in dead ends, M > 100 nearly always reaches 100%).
- Book ref [113] (Berry & Boudol, the chemical abstract machine) is cited in 2.5.2 as a previous study of the number-division chemistry; it is kept as the book gives it, but it is about CHAM, not this chemistry. The algorithm and numbers come from [72] and the appendix.

## Notes

The canonical example of an implicit/constructive chemistry, and the canonical argument for why a Chemart chemistry cannot always be a finite matrix: |S| is unbounded, and the reachable subnetwork depends on the seed multiset. Chemart exposes it both as the book's soup run and as the reachability closure of the seed numbers.

---

*Specification: `catalog/chemistries/prime-number-chemistry.yaml` · generator: `chemart/chemistries/prime_number_chemistry.py` · tests: `tests/chemistries/test_prime_number_chemistry.py`*
