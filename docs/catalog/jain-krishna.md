# Jain-Krishna autocatalytic set model

`jain-krishna` · *Jain & Krishna, 1998-2002*

The life cycle of an autocatalytic set, stripped to a graph. Nodes are species and a directed edge means one catalyses another's production; concentrations follow from the graph and the least-populated node is periodically replaced at random. Nothing happens for a while - then the first small autocatalytic cycle appears and the population organises around it and grows fast. Later a new cycle can undercut the incumbent, causing a crash and a reorganisation of the core. Innovation, dominance and collapse from one rule.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 15.2.2, 15.2.3 |
| **refs** | [425], [426], [427], [428] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): m nodes of a directed graph; only the CATALYTIC action is represented, not the underlying reactions

**R — reactions** (explicit, arity 1): c_ij = 1 means X_j catalyses the production of X_i: X_j -> X_j + X_i ; x_i' = sum_j c_ij x_j - x_i sum_kj c_kj x_j

**A — reactor**: ode
 · *dilution:* normalisation sum_i x_i = 1

## What you get

```python
net = chemart.generate_network("jain-krishna", seed=1)
```

```
jain-krishna: 100 species, 22 reactions, status=complete
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: analysis
```

First reactions:

```
X59 -> X59 + X5  [mass-action k=1.0]
X30 -> X30 + X14  [mass-action k=1.0]
X38 -> X38 + X14  [mass-action k=1.0]
X31 -> X31 + X26  [mass-action k=1.0]
X46 -> X46 + X33  [mass-action k=1.0]
X2 -> X2 + X45  [mass-action k=1.0]
X52 -> X52 + X45  [mass-action k=1.0]
X52 -> X52 + X49  [mass-action k=1.0]
… and 14 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `m` | `int` | `100` | structural | number of nodes (the state-space dimension stays constant) <br>`2` … `300` |
| `p` | `float` | `0.0025` | structural | link probability, both initially and on rewiring; the phase parameter <br>`0` … `1` |
| `graph_updates` | `int` | `` | population | slow-timescale steps (least-fit node replaced) applied before the network is returned; 0 returns the initial random graph <br>`0` … `20000` · *range:* the book's run finds the first autocatalytic set at update 2854 (m = 100, p = 0.0025) |
| `self_loops` | `bool` | `` | structural | allow c_ii = 1 (direct self-replicators); the model sets c_ii = 0 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- random phase -> growth phase (once the first small ACS appears) -> organising phase
- crashes and core shifts when a stronger ACS out-competes the incumbent
- attractor = eigenvector of the Perron-Frobenius eigenvalue lambda_1 of C

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book's prose says c_ij is 'a catalytic function of species i on species j', but eq. 15.7 (x_i' = sum_j c_ij x_j ...) needs j to catalyse i. The equation is followed: c_ij = 1 becomes X_j -> X_j + X_i with rate 1, under constant-total dilution.
- The fast attractor is the Perron-Frobenius eigenvector of C, recorded as the initial state. When there is no autocatalytic set (eigenvalue 0), (C + I)^m is iterated instead, which concentrates on the ends of the longest chains. Ties for the least-populated node are broken at random.
- After a replacement the book perturbs the remaining concentrations. That has no effect when the next attractor is computed exactly, so it is not modelled. extras.analysis records the eigenvalue and the attractor's support.
- x0, the concentration given to a new node, only matters for numerical integration, so it is not a parameter.

## Notes

Comes with a published TAXONOMY OF INNOVATIONS worth encoding as an analysis output: Type A short-lived (blip, incremental); Type B long-term - independent organisation (core-shift / birth of an organisation / dormant) and modification of an existing organisation (core-enhancing / neutral).

---

*Specification: `catalog/chemistries/jain-krishna.yaml` · generator: `chemart/chemistries/jain_krishna.py` · tests: `tests/chemistries/test_jain_krishna.py`*
