# Okamoto's biochemical switch

`okamoto-switch` · *Okamoto, Sakai & Hayashi, 1987-1993*

A biochemical flip-flop. Two enzyme cofactors interconvert, each conversion driven by one of two input-derived substrates, so whichever input dominates pushes the pair into the corresponding state and holds it there. The result is bistability with a sharp transition: the outputs swap within seconds of the inputs crossing, and how fast depends on how quickly inputs are converted. An early demonstration that chemistry can hold a bit.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 17.4.1 |
| **refs** | [632], [633], [634], [635], [631], [386], [387], [539] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): inputs I1, I2; substrates X1..X4; enzyme cofactors A and B (the outputs)

**R — reactions** (explicit, arity [1, 2]): I1 -> X1 ; I2 -> X3 ; A + X3 --k2--> B + X4 ; B + X1 --k1--> A + X2 ; X2 --k3--> ; X4 --k4--> ; plus I1 -> I2 to drive the inputs across each other

**A — reactor**: ode
 · *dilution:* X2 and X4 decay out

## What you get

```python
net = chemart.generate_network("okamoto-switch", seed=1)
```

```
okamoto-switch: 8 species, 7 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
I1 -> I1 + X1  [mass-action k=1.0]
I2 -> I2 + X3  [mass-action k=1.0]
A + X3 -> B + X4  [mass-action k=50000.0]
B + X1 -> A + X2  [mass-action k=50000.0]
X2 -> ∅  [mass-action k=10.0]
X4 -> ∅  [mass-action k=10.0]
I1 -> I2  [mass-action k=0.006]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `k1` | `float` | `50000.0` | kinetic | B + X1 -> A + X2 <br>≥ `0` |
| `k2` | `float` | `50000.0` | kinetic | A + X3 -> B + X4 <br>≥ `0` |
| `k3` | `float` | `10.0` | kinetic | decay of X2 <br>≥ `0` |
| `k4` | `float` | `10.0` | kinetic | decay of X4 <br>≥ `0` |
| `k_in` | `float` | `1.0` | kinetic | production of X1 from I1 and of X3 from I2 <br>≥ `0` |
| `k_conv` | `float` | `0.006` | kinetic | conversion I1 -> I2; sets when the inputs cross and hence the switching time <br>≥ `0` |
| `I1_0` | `float` | `100.0` | population | initial concentration of input I1 <br>≥ `0` |
| `I2_0` | `float` | `80.0` | population | initial concentration of input I2 <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- bistability: A and B flip within seconds of the inputs crossing
- switching time tunable via the input conversion rate
- composable into chemical neurons, logic gates, finite state machines and a Turing-universal chemical computer

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Inputs are 'controlled externally ... not affected although consumed' (book), so input reactions are written I1 -> I1 + X1 and I2 -> I2 + X3, with rate k_in (not given in the book; default 1).
- The book's I1/I2 profiles are reproduced as in its simulation: I1(0) = 100, I2(0) = 80 and an irreversible I1 -> I2. k_conv = 0.006 makes the inputs cross at about 17.5 s, as in figure 17.5 (top).
- Known discrepancy: with the reactions as printed, X1 (X3) is never lost while A (B) is saturated, so surplus substrate accumulates and the flip lags the crossing by ~18 s at the defaults, instead of 'a few seconds' in figure 17.5. The original cyclic-enzyme model [633] has reversible enzyme steps and substrate influx/efflux that the book omits; the paper was not accessible, so the book's scheme is kept.
- Other initial concentrations are the book's: X1 = X3 = 0, X2 = X4 = 8, A = 1, B = 0.

---

*Specification: `catalog/chemistries/okamoto-switch.yaml` · generator: `chemart/chemistries/okamoto_switch.py` · tests: `tests/chemistries/test_okamoto_switch.py`*
