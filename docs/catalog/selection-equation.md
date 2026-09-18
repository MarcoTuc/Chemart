# Selection equation under a dilution flow

`selection-equation`

Selection in its barest form: each species replicates at its own fitness while a non-selective outflow removes everyone at the mean fitness, holding the total constant. The result is competitive exclusion - the fittest species takes over and the rest vanish. Changing the exponent on the replication term changes that conclusion, which matters for real template replication of short oligonucleotides, where growth is parabolic rather than exponential and coexistence becomes possible.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.5 |
| **refs** | [621] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `rate-law`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): n species X_1..X_n with fitnesses f_i >= 0

**R — reactions** (explicit, arity 1): X_i --f_i--> 2 X_i ; X_i --phi--> (nothing), with phi = sum_i f_i x_i (the mean fitness)

**A — reactor**: ode, ssa
 · *dilution:* non-selective outflow phi keeps sum_i x_i = 1

## What you get

```python
net = chemart.generate_network("selection-equation", seed=1)
```

```
selection-equation: 3 species, 3 reactions, status=complete
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X1 -> 2 X1  [mass-action k=1.0]
X2 -> 2 X2  [mass-action k=2.0]
X3 -> 2 X3  [mass-action k=3.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `f` | `list` | `[1.0, 2.0, 3.0]` | selection | fitness vector; its length sets the number of species |
| `c` | `float` | `1.0` | kinetic | growth exponent in the generalised form x_i' = f_i x_i^c - phi x_i; c<1 gives survival of everybody, c=1 survival of the fittest, c>1 survival of the first <br>≥ `0` |
| `x0` | `list` | `` | population | initial concentrations; empty means uniform 1/n |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- competitive exclusion
- c = 0.5 (parabolic growth) applies to template replication of short oligonucleotides and BLOCKS Darwinian evolution unless the replicators are enclosed in a dividing compartment

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The dilution flux phi is outflow = 'constant-total'. For c != 1 the growth reaction carries a 'power' rate law k x^c.
- The book gives no numbers; default fitnesses 1, 2, 3 and a uniform start.

---

*Specification: `catalog/chemistries/selection-equation.yaml` · generator: `chemart/chemistries/selection_equation.py` · tests: `tests/chemistries/test_selection_equation.py`*
