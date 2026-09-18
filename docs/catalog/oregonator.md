# Oregonator (Belousov-Zhabotinsky)

`oregonator` · *Field & Noyes, 1974*

The reduced model of the Belousov-Zhabotinsky reaction - a real chemical oscillator, not an invented one. Five reactions couple an autocatalytic species, an inhibitor and an oxidised metal ion, with two reactants buffered; the autocatalysis fires, exhausts itself, regenerates the inhibitor and waits. Add diffusion and it becomes an excitable medium supporting target patterns and spiral waves, which is the basis of reaction-diffusion computing.

| | |
|---|---|
| **family** | wet |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 19.3.3 |
| **refs** | [276], [9], [13], [813], [857] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `space`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): A, B (buffered reactants), X (autocatalyst, HBrO2), Y (inhibitor, Br-), Z (oxidised metal ion, e.g. Ce(IV)), P, Q (products)

**R — reactions** (explicit, arity [1, 2]): A + Y --k1--> X ; X + Y --k2--> P ; B + X --k3--> 2X + Z ; 2X --k4--> Q ; Z --k5--> f Y

**A — reactor**: ode, continuous-space
 · *dilution:* A and B buffered; reaction-diffusion on a grid for spatial patterns

## What you get

```python
net = chemart.generate_network("oregonator", seed=1)
```

```
oregonator: 7 species, 5 reactions, status=complete
provides: catalysts, initial-state, rate-constants, space, stoichiometry, topology
seed: 1
extras: buffered, space
```

First reactions:

```
A + Y -> X  [mass-action k=1.0]
X + Y -> P  [mass-action k=10.0]
B + X -> 2 X + Z  [mass-action k=10.0]
2 X -> Q  [mass-action k=2.5]
Z -> Y  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `k1` | `float` | `1.0` | kinetic | A + Y -> X <br>≥ `0` |
| `k2` | `float` | `10.0` | kinetic | X + Y -> P <br>≥ `0` |
| `k3` | `float` | `10.0` | kinetic | B + X -> 2X + Z <br>≥ `0` |
| `k4` | `float` | `2.5` | kinetic | 2X -> Q <br>≥ `0` |
| `k5` | `float` | `1.0` | kinetic | Z -> f Y <br>≥ `0` |
| `f` | `float` | `1.0` | kinetic | stoichiometric factor of inhibitor regeneration; the key bifurcation parameter <br>≥ `0` |
| `A` | `float` | `1.0` | population | buffered concentration of A <br>≥ `0` |
| `B` | `float` | `1.0` | population | buffered concentration of B <br>≥ `0` |
| `D` | `dict` | `{'X': 1.0, 'Y': 1.0, 'Z': 1.0}` | spatial | diffusion coefficient per species for the reaction-diffusion setting |
| `lattice` | `enum` | `hexagonal` | spatial | grid used for reaction-diffusion (the book's patterns use a hexagonal grid) <br>one of `square`, `hexagonal` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- target patterns and spiral waves
- excitable medium supporting reaction-diffusion computing (logic gates, image processing, maze solving)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Rate constants are the book's (figure 19.20). f, the buffered concentrations A and B, and the diffusion coefficients are not given; defaults are f = 1, A = B = 1, D = 1.
- A real f cannot be an integer stoichiometric coefficient. Z -> f Y is emitted as Z -> floor(f) Y at rate k5 (1 - frac(f)) plus Z -> (floor(f)+1) Y at rate k5 frac(f), which gives exactly the same rate equations.
- Space is recorded in extras.space (lattice, diffusion); A and B are listed in extras.buffered.

---

*Specification: `catalog/chemistries/oregonator.yaml` · generator: `chemart/chemistries/oregonator.py` · tests: `tests/chemistries/test_oregonator.py`*
