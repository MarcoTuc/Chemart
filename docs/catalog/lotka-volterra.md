# Lotka-Volterra

`lotka-volterra`

Predator and prey as three reactions: prey eats an inexhaustible resource and doubles, a predator meeting prey doubles, and predators die. Deterministically the populations circle a fixed point forever. In a small volume the same rules behave quite differently - oscillations grow erratic and one species eventually goes extinct by chance - which makes it a compact demonstration that the reactor algorithm is part of the model.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.8; appendix (Lotka.py); ecology 8.2.3 |
| **refs** | [621], [792], [37], [238], [593] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): n interacting species X_1..X_n; classically prey X1 and predator X2

**R — reactions** (explicit, arity [1, 2]): X + G --ka--> 2X + G ; Y + X --kb--> 2Y ; Y --kc--> (nothing)

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("lotka-volterra", seed=1)
```

```
lotka-volterra: 2 species, 3 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X1 -> 2 X1  [mass-action k=1.0]
X2 -> ∅  [mass-action k=1.0]
X2 + X1 -> 2 X2  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `r` | `list` | `[1.0, -1.0]` | kinetic | spontaneous growth (>0) or death (<0) rate per species; its length sets n |
| `B` | `list` | `[[0.0, -1.0], [1.0, 0.0]]` | kinetic | interaction matrix b_ij: 0 none, >0 j feeds i, <0 j eats i; both >0 mutualistic, both <0 competitive, mixed antagonistic |
| `x0` | `list` | `[5.0, 2.0]` | population | initial concentrations |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- sustained oscillations around (x, y) = (kc/kb, ka g0/kb) under ODE
- erratic oscillations and extinction under SSA in small volumes
- cyclic food chains, waves, clustering

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The general form (eq. 7.30) is the parameterisation. Each term becomes the mass-action reaction with the same ODE contribution, and a predator-prey pair b_ij = -b_ji > 0 becomes the single reaction X_i + X_j -> 2 X_i (eq. 7.26).
- Defaults are figure 7.5's predator-prey system (ka = kb = kc = 1, x0 = 5, y0 = 2, g0 = 1), with the constant grass G folded into r_1 = ka g0.
- NAV is not a parameter: it sizes a stochastic simulation, not the network.

---

*Specification: `catalog/chemistries/lotka-volterra.yaml` · generator: `chemart/chemistries/lotka_volterra.py` · tests: `tests/chemistries/test_lotka_volterra.py`*
