# Random catalytic reaction networks

`random-catalytic-networks` · *Stadler, Fontana & Miller, 1993*

The null model the more elaborate chemistries all instantiate. Species have no internal structure at all; instead, for every pair, a randomly drawn coefficient says how strongly that pair catalyses the production of a third. With a dilution flux, the system settles onto a mutually-catalysing subset and the rest dies out. Because it contains the replicator and Lotka-Volterra equations as special cases, it is the baseline against which structured chemistries have to justify their structure.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.9 |
| **refs** | [792] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow` |

## Molecules, reactions, reactor

**S — molecules** (explicit): n abstract species X_1..X_n with no internal structure

**R — reactions** (explicit, arity 2): X_i + X_j --alpha^k_ij--> X_i + X_j + X_k (both reactants catalytic)

**A — reactor**: ode
 · *dilution:* non-specific dilution flux phi

## What you get

```python
net = chemart.generate_network("random-catalytic-networks", seed=1)
```

```
random-catalytic-networks: 10 species, 42 reactions, status=complete
provides: catalysts, flow, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X1 + X2 -> X1 + X2 + X3  [mass-action k=0.7535131086748066]
X1 + X5 -> X1 + X5 + X4  [mass-action k=0.5285892632600216]
X1 + X5 -> X1 + X5 + X7  [mass-action k=0.641328169139375]
X1 + X7 -> X1 + X7 + X5  [mass-action k=0.8552269742870702]
X1 + X8 -> X1 + X8 + X2  [mass-action k=0.6457208955749478]
X1 + X9 -> X1 + X9 + X7  [mass-action k=0.6734598871529389]
X1 + X10 -> X1 + X10 + X8  [mass-action k=0.8254878133935558]
2 X2 -> 2 X2 + X7  [mass-action k=0.2624947127501015]
… and 34 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n` | `int` | `10` | structural | number of species <br>`2` … `60` |
| `density` | `float` | `0.1` | structural | probability that a given reaction X_i + X_j -> X_i + X_j + X_k exists; sparse to dense <br>`0` … `1` |
| `allow_direct_replication` | `bool` | `` | structural | allow k = i or k = j; the book's interesting case forbids direct self-replication |
| `alpha_distribution` | `enum` | `uniform` | kinetic | distribution of the rate alpha of an existing reaction: uniform on [0, 1) or exponential with mean 1 <br>one of `uniform`, `exponential` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- convergence to a fixpoint where only a mutually catalytic subset survives - later recognised as a chemical organisation
- contains the replicator and Lotka-Volterra equations as special cases

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives the reaction scheme, the ODE and the direct-replication switch, but not how alpha is drawn. Each existing reaction gets alpha from `alpha_distribution` (default uniform on [0, 1)).
- Reactions are sampled per unordered reactant pair {i, j}. In the ODE the ordered terms alpha^k_ij and alpha^k_ji describe the same reaction, so one rate per pair loses nothing. The dilution phi is outflow = 'constant-total'.

## Notes

The structural template that the matrix chemistry, the automata reaction and AlChemy all instantiate: Chemart's canonical 'null model' generator for benchmarking organisation-finding code.

---

*Specification: `catalog/chemistries/random-catalytic-networks.yaml` · generator: `chemart/chemistries/random_catalytic_networks.py` · tests: `tests/chemistries/test_random_catalytic_networks.py`*
