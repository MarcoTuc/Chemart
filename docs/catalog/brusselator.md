# Brusselator

`brusselator` · *Prigogine & Lefever, 1968*

The minimal oscillator, and the standard first test for anything that handles reaction networks. Two species are held constant as reservoirs, two are dynamic, and the autocatalytic step 2X + Y -> 3X provides the feedback. Push the parameter b past 1 + a^2 and the steady state loses stability: concentrations stop settling and run around a limit cycle instead. The trimolecular step is chemically unrealistic and deliberately so - it buys analytic tractability.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book` — implemented exactly as the book specifies |
| **book** | 17.4.2, 19.3; ARMS 9.4 |
| **refs** | [675] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): A, B (buffered), X, Y (dynamic), D, E (waste)

**R — reactions** (explicit, arity [1, 2, 3]): A --k1--> X ; B + X --k2--> Y + D ; 2X + Y --k3--> 3X ; X --k4--> E

**A — reactor**: ode, ssa
 · *dilution:* A and B held constant by buffering

## What you get

```python
net = chemart.generate_network("brusselator", seed=1)
```

```
brusselator: 6 species, 4 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: buffered
```

First reactions:

```
A -> X  [mass-action k=1.0]
B + X -> Y + D  [mass-action k=1.0]
2 X + Y -> 3 X  [mass-action k=1.0]
X -> E  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `a` | `float` | `1.0` | population | [A], held constant <br>≥ `0` |
| `b` | `float` | `3.0` | population | [B], held constant; sustained oscillations require b > 1 + a^2 when all k_i = 1 <br>≥ `0` |
| `k1` | `float` | `1.0` | kinetic | A -> X <br>≥ `0` |
| `k2` | `float` | `1.0` | kinetic | B + X -> Y + D <br>≥ `0` |
| `k3` | `float` | `1.0` | kinetic | 2X + Y -> 3X <br>≥ `0` |
| `k4` | `float` | `1.0` | kinetic | X -> E <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- limit-cycle oscillations for b > 1 + a^2 (all k_i = 1); damped oscillations to (X, Y) = (a, b/a) below

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Buffered species A and B are listed in extras.buffered: their concentrations stay at a and b although the reactions consume them. Initial X, Y, D, E are 0, as in figure 17.6.

## Notes

The trimolecular step 2X + Y -> 3X is physically unrealistic; the book keeps it because it simplifies the analysis. A good test that Chemart handles arity-3 reactions and the associated l_i! factor in the k->c conversion.

---

*Specification: `catalog/chemistries/brusselator.yaml` · generator: `chemart/chemistries/brusselator.py` · tests: `tests/chemistries/test_brusselator.py`*
