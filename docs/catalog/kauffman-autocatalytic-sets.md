# Kauffman autocatalytic sets (binary polymer model)

`kauffman-autocatalytic-sets` · *Kauffman, 1986; Farmer, Kauffman & Packard, 1986*

The original argument that metabolism could precede genes. Polymers join and split, and any polymer may catalyse any such reaction with some probability. As polymers get longer the number of possible reactions grows faster than the number of molecules, so above a critical catalysis probability a mutually-catalysing set almost surely exists - order for free, from combinatorics rather than design. This entry gives the topology; the dynamics and the formal test for such sets live in bagley-farmer and raf.

| | |
|---|---|
| **family** | origin-of-life |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 6.3.1 |
| **refs** | [446], [447], [260] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `mass-conservation`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (implicit): polymers = strings over an alphabet of size B (a, b for B = 2), up to a maximum length

**R — reactions** (implicit, arity [1, 2, 3]): condensation A + B -> AB and cleavage AB -> A + B for every split of every polymer, each optionally catalysed by a polymer E: A + B + E <-> AB + E

**A — reactor**: graph-rewrite
 · *dilution:* food set inflow keeps the system out of equilibrium

## What you get

```python
net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
```

```
kauffman-autocatalytic-sets: 62 species, 1608 reactions, status=complete
provides: catalysts, mass-conservation, stoichiometry, topology
seed: 1
extras: conservation, food_set
```

First reactions:

```
2 a -> aa
aa -> 2 a
a + b -> ab
ab -> a + b
b + a -> ba
ba -> b + a
2 b -> bb
bb -> 2 b
… and 1600 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `B` | `int` | `2` | structural | alphabet size <br>`2` … `4` |
| `max_length` | `int` | `5` | structural | longest polymer in the network (at most 512 polymers in total) <br>`2` … `8` |
| `P` | `float` | `0.05` | structural | probability that a given polymer catalyses a given condensation/cleavage pair <br>`0` … `1` · *range:* critical value P_crit ~ B^(-2L) for a firing disk of polymers up to length L (eq. 6.3) |
| `food_set` | `list` | `['a', 'b', 'aa', 'bb']` | population | the firing disk of polymers supplied from outside (book figure 6.12) |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- percolation transition at P_crit: above it, autocatalytic sets form almost surely
- larger B favours set formation, so proteins (B=20) beat RNA (B=4)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Built as the binary polymer model: all polymers up to max_length, every split of every polymer as a reversible condensation/cleavage pair, and each polymer catalysing each pair (both directions) with probability P.
- Water (H in eq. 6.2) is left out, as in the usual formulation of the polymer model. No rates: the book's analysis is static. The food set is recorded in extras.food_set, and monomer counts per letter are exact conservation laws (extras.conservation).
- The earlier parameter L (maximum food-set length) is replaced by an explicit food_set; the book's default food set is {a, b, aa, bb}.

## Notes

Topology only; see `bagley-farmer` for the dynamic version. RAF analysis (`raf`) applies directly to this network.

---

*Specification: `catalog/chemistries/kauffman-autocatalytic-sets.yaml` · generator: `chemart/chemistries/kauffman_autocatalytic_sets.py` · tests: `tests/chemistries/test_kauffman_autocatalytic_sets.py`*
