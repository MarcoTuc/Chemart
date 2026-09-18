# N-economy (natural number economy)

`n-economy` · *Straatman, Banzhaf et al.; after Hanel's natural number economy*

An economy written as a chemistry of numbers. Natural numbers are commodities: primes are raw materials, composites are manufactured goods, and production reactions consume and produce them with integer stoichiometry, with 2 standing in for labour and 3 for money. Since a composite is exactly its prime factorisation, a good's *composition* is built into its name, and conservation follows from arithmetic rather than being imposed.

| | |
|---|---|
| **family** | non-chemical |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 20.3 |
| **refs** | [816], [377] |
| **provides** | `topology`, `stoichiometry` |

## Molecules, reactions, reactor

**S — molecules** (implicit): natural numbers as commodities; primes are raw materials, composites are manufactured goods; 2 = labour, 3 = money

**R — reactions** (explicit, arity variable): Production reactions with integer stoichiometry, e.g. 2x2 -> 2x5 ; 2x2 + 2x5 + 2x13 -> 3x130 ; 3x2 + 3x130 -> 6x260 ; 6x2 + 6x260 -> 6x104 ; 6x104 -> 13x2. Generally sum_i alpha_i G_i -> sum_i beta_i G_i + sum_j G_{m+j}, with beta_i <= alpha_i.

**A — reactor**: lattice-2d
 · *dilution:* consumption agents; a market prices goods

## What you get

```python
net = chemart.generate_network("n-economy", seed=1)
```

```
n-economy: 7 species, 5 reactions, status=complete
provides: stoichiometry, topology
seed: 1
```

First reactions:

```
2 2 -> 2 5
2 2 + 2 5 + 2 13 -> 3 130
3 2 + 3 130 -> 6 260
6 2 + 6 260 -> 6 104
6 104 -> 13 2
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `product_set` | `list` | `[2, 3, 5, 13, 104, 130, 260]` | structural | the goods (integers >= 2); primes are raw materials, 2 is labour and 3 is money |
| `technology` | `list` | `[{'in': {'2': 2}, 'out': {'5': 2}}, {'in': {'2'…` | structural | production processes as {in: {good: amount}, out: {good: amount}}; default is book eq. 20.6 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- noisy rise and fall of product concentrations, equilibrating after ~300 iterations (agent simulation)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Only the production system (eq. 20.6) is generated. The agent-based consumption side, the market and the grid are outside the reaction network, so agent_skills and grid are not parameters.
- The earlier catalog note claimed prime factorisation gives an exact conservation law. It does not for the book's reactions (e.g. 2x2 -> 2x5), so mass-conservation is not claimed. Species carry their factorisation as structure.

---

*Specification: `catalog/chemistries/n-economy.yaml` · generator: `chemart/chemistries/n_economy.py` · tests: `tests/chemistries/test_n_economy.py`*
