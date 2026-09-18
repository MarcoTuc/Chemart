# Replication and death

`replication-death` · *textbook evolutionary dynamics*

Two reactions - birth and death - and nothing else. Its value is as a floor: exponential growth when births outpace deaths is the minimal precondition for Darwinian evolution, so this is the smallest system in which selection could mean anything at all, and the baseline any richer population model reduces to.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.1 |
| **refs** | [621] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): one species X

**R — reactions** (explicit, arity 1): X --b--> 2X ; X --d--> (nothing)

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("replication-death", seed=1)
```

```
replication-death: 1 species, 2 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X -> 2 X  [mass-action k=1.0]
X -> ∅  [mass-action k=0.5]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `b` | `float` | `1.0` | kinetic | birth (replication) rate <br>≥ `0` |
| `d` | `float` | `0.5` | kinetic | death rate <br>≥ `0` |
| `x0` | `float` | `1.0` | population | initial population <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- x(t) = x0 exp((b-d)t); exponential growth is the precondition for Darwinian evolution

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives no numbers; defaults b = 1, d = 0.5, x0 = 1 give exponential growth.

---

*Specification: `catalog/chemistries/replication-death.yaml` · generator: `chemart/chemistries/replication_death.py` · tests: `tests/chemistries/test_replication_death.py`*
