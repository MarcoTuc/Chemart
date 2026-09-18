# Logistic growth ('replicate and fight')

`logistic-chemistry`

Logistic growth expressed as two reactions: X duplicates, and two X's collide to leave one. The second reaction supplies the crowding term, so the population approaches its carrying capacity without any bound being imposed. Worth keeping because the deterministic and stochastic readings genuinely differ - a stochastic run fluctuates around the capacity and can overshoot it, since nothing in the rules forbids that.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.4; appendix (Logistic.py) |
| **refs** | [621] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): one species X

**R — reactions** (explicit, arity [1, 2]): X --r--> 2X ; 2X --r/K--> X

**A — reactor**: ode, ssa
 · *dilution:* implicit, via the reverse reaction

## What you get

```python
net = chemart.generate_network("logistic-chemistry", seed=1)
```

```
logistic-chemistry: 1 species, 2 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X -> 2 X  [mass-action k=1.0]
2 X -> X  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `r` | `float` | `1.0` | kinetic | intrinsic growth rate <br>≥ `0` |
| `K` | `float` | `1.0` | population | carrying capacity; the reverse coefficient is d = r/K <br>≥ `1e-09` |
| `x0` | `float` | `0.1` | population | initial population <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- sigmoid approach to K under ODE
- stochastic runs fluctuate around and OVERSHOOT K - the reaction does not enforce a hard bound

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- K = 1 as in figure 7.2; r and x0 are not given (defaults 1 and 0.1).

---

*Specification: `catalog/chemistries/logistic-chemistry.yaml` · generator: `chemart/chemistries/logistic_chemistry.py` · tests: `tests/chemistries/test_logistic_chemistry.py`*
