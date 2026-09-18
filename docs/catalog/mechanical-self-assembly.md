# Mechanical self-assembly (Hosokawa)

`mechanical-self-assembly` · *Hosokawa, Shimoyama & Miura, 1994*

Proof that chemical kinetics does not need chemicals. Plastic triangles with magnets in two faces are shaken in a dish; they stick into dimers, trimers and larger assemblies, and the yields over time fit an ordinary mass-action reaction network. The instructive detail is where the rate constant comes from: it factorises into how often two pieces collide and the probability they are correctly oriented when they do - an encounter term times a geometric term.

| | |
|---|---|
| **family** | non-chemical |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 20.1 |
| **refs** | [404], [352], [918], [919], [69] |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `mass-conservation`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): triangular plastic bodies with permanent magnets in two of three faces; x = monomer, x2..x6 = assemblies

**R — reactions** (explicit, arity 2): x+x -> x2 ; x+x2 -> x3 ; x+x3 -> x4 ; x+x4 -> x5 ; x+x5 -> x6 ; x2+x2 -> x4 ; x2+x3 -> x5 ; x2+x4 -> x6 ; x3+x3 -> x6

**A — reactor**: ode
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("mechanical-self-assembly", seed=1)
```

```
mechanical-self-assembly: 6 species, 9 reactions, status=complete
provides: initial-state, mass-conservation, rate-constants, stoichiometry, topology
seed: 1
extras: conservation
```

First reactions:

```
2 x -> x2  [mass-action k=1.0]
x + x2 -> x3  [mass-action k=1.0]
x + x3 -> x4  [mass-action k=1.0]
x + x4 -> x5  [mass-action k=1.0]
x + x5 -> x6  [mass-action k=1.0]
2 x2 -> x4  [mass-action k=1.0]
x2 + x3 -> x5  [mass-action k=1.0]
x2 + x4 -> x6  [mass-action k=1.0]
… and 1 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n_monomers` | `int` | `100` | population | initial number of monomers x (book figure 20.3a) <br>≥ `1` |
| `agitation_rate` | `float` | `1.0` | kinetic | collision frequency; rate = agitation_rate * P_b <br>≥ `0` |
| `P_b` | `dict` | `` | kinetic | bonding probability per pair, keyed 'x+x2' etc.; pairs not listed use 1.0 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- yield curves for x4, x5, x6 matching wet-lab experiments

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives no numerical bonding probabilities ('under reasonable assumptions'); every P_b defaults to 1 and can be overridden per pair.
- The book's discrete-time update is represented by mass-action rates k = agitation_rate * P_b. extras.conservation records the monomer count (x=1 ... x6=6).

## Notes

Proof that chemical kinetics is substrate-independent, and the cleanest example of a rate coefficient that FACTORISES into an encounter term and a geometric orientation term.

---

*Specification: `catalog/chemistries/mechanical-self-assembly.yaml` · generator: `chemart/chemistries/mechanical_self_assembly.py` · tests: `tests/chemistries/test_mechanical_self_assembly.py`*
