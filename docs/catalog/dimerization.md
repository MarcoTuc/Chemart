# Reversible dimerization

`dimerization` · *textbook; used as PyCellChemistry's hello-world*

The hello-world of reaction networks: A and B bind to a dimer and it falls apart again. There is nothing to discover here, which is exactly the point - it is the fixture where you check that equilibrium comes out as k_forward/k_reverse, that a stochastic run fluctuates around the deterministic answer, and that the conversion between macroscopic and stochastic rate constants is done right.

| | |
|---|---|
| **family** | core |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book` — implemented exactly as the book specifies |
| **book** | 2.2.4, appendix (Dimer.py) |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): three named species A, B and the dimer C = A.B

**R — reactions** (explicit, arity [1, 2]): A + B --k_f--> C ; C --k_r--> A + B

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("dimerization", seed=1)
```

```
dimerization: 3 species, 2 reactions, status=complete
provides: initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
A + B -> C  [mass-action k=1.0]
C -> A + B  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `k_f` | `float` | `1.0` | kinetic | forward (association) rate coefficient <br>≥ `0` |
| `k_r` | `float` | `1.0` | kinetic | reverse (dissociation) rate coefficient <br>≥ `0` |
| `A0` | `float` | `2.0` | population | initial concentration of A <br>≥ `0` |
| `B0` | `float` | `1.4` | population | initial concentration of B <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- equilibrium K = k_f/k_r
- SSA fluctuates about the ODE equilibrium

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- NAV (N_A V) is not a parameter: it sizes a stochastic simulation, not the network. Convert rates for SSA with chemart.kinetics.k_to_c.

## Notes

The trivial baseline, and the fixture for the k (macroscopic) -> c (mesoscopic) conversion via the Wolkenhauer relation c = k / (N_A V)^(m-1) * prod_i l_i! (book appendix, eq. 4).

---

*Specification: `catalog/chemistries/dimerization.yaml` · generator: `chemart/chemistries/dimerization.py` · tests: `tests/chemistries/test_dimerization.py`*
