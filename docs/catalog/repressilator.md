# Repressilator

`repressilator` · *Elowitz & Leibler, 2000*

A synthetic oscillator that was actually built in bacteria. Three genes are arranged in a cycle so each represses the next, with repression modelled as cooperative binding and both mRNA and protein decaying. With enough cooperativity the loop cannot settle - the repression travels round the ring and the protein levels oscillate. With no cooperativity it simply relaxes, so the nonlinearity is doing the work, and the oscillation only exists inside a bounded region of parameter space.

| | |
|---|---|
| **family** | wet |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book` — implemented exactly as the book specifies |
| **book** | 19.3.2; appendix (Repressilator.py) |
| **refs** | [256], [518], [926] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): 3 genes G_i, 3 gene-protein complexes C_i, 3 mRNAs M_i, 3 repressor proteins P_i (LacI, TetR, lambda-cI)

**R — reactions** (explicit, arity n+1, 1): G_i + n P_{i-1} <-ke/kr-> C_i (cyclic: P3 represses G1) ; G_i --km--> G_i + M_i ; M_i --kp--> M_i + P_i ; M_i --mu_m--> ; P_i --mu_p-->

**A — reactor**: ode, ssa
 · *dilution:* mRNA and protein decay

## What you get

```python
net = chemart.generate_network("repressilator", seed=1)
```

```
repressilator: 12 species, 18 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
G1 + 2 P3 -> C1  [mass-action k=1.0]
C1 -> G1 + 2 P3  [mass-action k=1.0]
G2 + 2 P1 -> C2  [mass-action k=1.0]
C2 -> G2 + 2 P1  [mass-action k=1.0]
G3 + 2 P2 -> C3  [mass-action k=1.0]
C3 -> G3 + 2 P2  [mass-action k=1.0]
G1 -> G1 + M1  [mass-action k=5.0]
M1 -> M1 + P1  [mass-action k=1.0]
… and 10 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n` | `int` | `2` | kinetic | Hill coefficient (proteins per binding event); n=1 gives NO oscillations - the circuit needs a sigmoid response <br>`1` … `10` |
| `ke` | `float` | `1.0` | kinetic | repressor binding <br>≥ `0` |
| `kr` | `float` | `1.0` | kinetic | repressor unbinding <br>≥ `0` |
| `km` | `float` | `5.0` | kinetic | transcription <br>≥ `0` |
| `kp` | `float` | `1.0` | kinetic | translation <br>≥ `0` |
| `mu_m` | `float` | `0.5` | kinetic | mRNA decay <br>≥ `0` |
| `mu_p` | `float` | `0.1` | kinetic | protein decay <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- sustained oscillations inside a bounded parameter region (n=2 with the default rates)
- no oscillations for n = 1
- oscillations survive stochastic noise at low gene copy number
- implemented in real E. coli

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The SSA gene copy number m is not a parameter: it sizes a stochastic simulation (V = m/N_A), not the network.
- Binding is the elementary mass-action reaction of eqs. 19.1-19.3, not a lumped Hill rate law. Initial state is figure 19.18's: G1 = C2 = C3 = 1, everything else 0.

---

*Specification: `catalog/chemistries/repressilator.yaml` · generator: `chemart/chemistries/repressilator.py` · tests: `tests/chemistries/test_repressilator.py`*
