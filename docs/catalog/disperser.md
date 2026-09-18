# Chemical disperser (load balancing)

`disperser` · *Meyer & Tschudin, 2009*

A distributed algorithm that happens to be a chemistry. Each network node holds a species whose concentration is its workload, and each directed link carries a catalyst that converts a job at one end into a job at the other. Run it and the loads provably converge to the network-wide average, even while new jobs arrive. Mass conservation does the safety argument for free: the algorithm can never invent or lose work, so it cannot produce a negative load.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 17.3.1; appendix (Disperser.py) |
| **refs** | [576], [580], [941] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `compartments`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): job molecules X_i (one species per node) and catalysts C_i_j (one per directed link, located at node i)

**R — reactions** (explicit, arity 2): C_ij + X_i --k--> C_ij + X_j, for all (i,j) in E

**A — reactor**: ssa, compartments
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("disperser", seed=1)
```

```
disperser: 12 species, 8 reactions, status=complete
provides: catalysts, compartments, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: compartments
```

First reactions:

```
C1_2 + X1 -> C1_2 + X2  [mass-action k=1.0]
C1_3 + X1 -> C1_3 + X3  [mass-action k=1.0]
C1_4 + X1 -> C1_4 + X4  [mass-action k=1.0]
C2_1 + X2 -> C2_1 + X1  [mass-action k=1.0]
C2_3 + X2 -> C2_3 + X3  [mass-action k=1.0]
C3_1 + X3 -> C3_1 + X1  [mass-action k=1.0]
C3_2 + X3 -> C3_2 + X2  [mass-action k=1.0]
C4_1 + X4 -> C4_1 + X1  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `graph` | `list` | `[[1, 2], [1, 3], [1, 4], [2, 3]]` | spatial | undirected links [u, v] of the network; every link carries a catalyst in each direction |
| `k` | `float` | `1.0` | kinetic | rate coefficient, the same for all links; node degree then sets the outflow rate <br>≥ `0` |
| `catalyst_concentration` | `float` | `1.0` | population | concentration of every catalyst C_i_j, held constant and equal everywhere <br>≥ `0` |
| `initial_jobs` | `dict` | `{'4': 1000}` | population | node id (as a string) -> initial number of job molecules; other nodes start empty |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- provable convergence to the network-wide average x_bar = sum_i x_i(0)/n, stable under injection and deletion (1000 jobs on 4 nodes -> 250 each)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book's 4-node topology (figure 17.3) is only drawn. Edges 1-2, 1-3, 1-4, 2-3 are the graph consistent with the text: after injection at node 4, nodes 2 and 3 are the farthest; after injection at node 2, node 4 is.
- Catalyst concentration is not given numerically; default 1. extras.compartments maps each node to its X and outgoing C species.

## Notes

The cleanest 'chemistry as a provable distributed algorithm' example: mass conservation means the algorithm can never produce a negative load, and the ODE can be checked for convergence and stability directly.

---

*Specification: `catalog/chemistries/disperser.yaml` · generator: `chemart/chemistries/disperser.py` · tests: `tests/chemistries/test_disperser.py`*
