# Analog computation of algebraic functions with concentrations

`analog-function-crn` · *Hjelmfelt et al.; Deckard & Sauro; Dittrich et al.*

Concentrations used as analog numbers. Pick reactions whose steady state satisfies the algebraic relation you want, then read the answer off the output species once it settles: X -> X + Y together with 2Y -> nothing balances at y = sqrt(k1 x / 2 k2), so a square root falls out of mass action. The idea is elegant and its limits are honest ones - no negative numbers, you must wait for equilibrium, not every network has one, and composing two such circuits is awkward.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 17.4.3 |
| **refs** | [141], [220], [511], [386] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): one species per input (X1, X2) and the output Y; values are concentrations

**R — reactions** (explicit, arity [1, 2]): sqrt: X --k1--> X + Y ; 2Y --k2--> (nothing), giving y = sqrt(k1 x / (2 k2)); k1 = 2 k2 gives y = sqrt(x)

**A — reactor**: ode
 · *dilution:* the result is read at steady state

## What you get

```python
net = chemart.generate_network("analog-function-crn", seed=1)
```

```
analog-function-crn: 2 species, 2 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: readout
```

First reactions:

```
X1 -> X1 + Y  [mass-action k=2.0]
2 Y -> ∅  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `function` | `enum` | `sqrt` | structural | the function computed at steady state as [Y] <br>one of `sqrt`, `square`, `add`, `multiply`, `divide` |
| `x1` | `float` | `4.0` | population | first input concentration <br>≥ `0` |
| `x2` | `float` | `2.0` | population | second input concentration (add, multiply, divide) <br>≥ `0` |

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- sqrt is the book's network with k1 = 2, k2 = 1 (eqs. 17.18-17.24).
- The other functions are not in the book. They are the analogous mass-action networks, each checked against its steady state in the tests: square 2X1 -> 2X1 + Y, Y -> ; add X1 -> X1 + Y, X2 -> X2 + Y, Y -> ; multiply X1 + X2 -> X1 + X2 + Y, Y -> ; divide X1 -> X1 + Y, X2 + Y -> X2.
- extras.readout names the output species and its expected steady state.

## Notes

The book lists the shortcomings honestly: no negative numbers, must wait for steady state, not every system HAS a steady state, and composition needs disjoint species per stage, which does not scale.

---

*Specification: `catalog/chemistries/analog-function-crn.yaml` · generator: `chemart/chemistries/analog_function_crn.py` · tests: `tests/chemistries/test_analog_function_crn.py`*
