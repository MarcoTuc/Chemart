# Quasispecies equation

`quasispecies` · *Eigen, 1971; Eigen & Schuster*

Mutation and selection as one reaction network. Each genotype replicates at its own fitness, but replication is error-prone, so producing a neighbour is part of reproducing. What survives is therefore not the fittest sequence but a *cloud* centred on it. Push the mutation rate past a threshold and the cloud stops being held together - the population delocalises over sequence space and the information is lost. That error threshold is the model's central result.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.7; appendix (Quasispecies.py) |
| **refs** | [249], [627], [621] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow` |

## Molecules, reactions, reactor

**S — molecules** (implicit): genome sequences of length L over an alphabet of size B; |S| = B^L

**R — reactions** (explicit, arity 1): X_j --f_j q_ji--> X_j + X_i ; X_i --phi--> (nothing)

**A — reactor**: ode, ssa
 · *dilution:* phi = mean fitness, population size constant

## What you get

```python
net = chemart.generate_network("quasispecies", seed=1)
```

```
quasispecies: 64 species, 4032 reactions, status=complete
provides: catalysts, flow, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
000001 -> 000001 + 000000  [mass-action k=0.00468933447645176]
000001 -> 2 000001  [mass-action k=0.13599069981710105]
000001 -> 000001 + 000010  [mass-action k=0.00016170118884316413]
000001 -> 000001 + 000011  [mass-action k=0.00468933447645176]
000001 -> 000001 + 000100  [mass-action k=0.00016170118884316413]
000001 -> 000001 + 000101  [mass-action k=0.00468933447645176]
000001 -> 000001 + 000110  [mass-action k=5.575903063557383e-06]
000001 -> 000001 + 000111  [mass-action k=0.00016170118884316413]
… and 4024 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `L` | `int` | `6` | structural | genome length <br>`1` … `8` · *range:* the book's figure 7.4 uses L = 10, above the generator limit of 256 genomes |
| `B` | `int` | `2` | structural | alphabet size <br>`2` … `4` |
| `m` | `float` | `0.2` | stochastic | mutations per genome; per-site rate p = m/L. Deterministic error threshold at m = 1; the stochastic threshold in finite populations is LOWER <br>≥ `0` · *range:* book: 0.2, 1.0, 2.0 |
| `fitness` | `enum` | `ones-fraction` | selection | ones-fraction: share of positions holding the highest symbol (the book's choice); single-peak: 2 for the all-highest genome, 1 otherwise <br>one of `ones-fraction`, `single-peak` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- quasispecies cloud around the fittest genotype
- error catastrophe above the threshold
- no adaptation below a minimum mutation rate

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- q_ji = (p/(B-1))^d (1-p)^(L-d), with d the Hamming distance and p = m/L. Reactions with zero rate are omitted.
- The full network has up to (B^L)^2 reactions, so B^L is limited to 256 genomes. The population size N is not a parameter: it sizes a stochastic simulation. No initial state is attached: the book only says the start is random and biased toward low fitness.

---

*Specification: `catalog/chemistries/quasispecies.yaml` · generator: `chemart/chemistries/quasispecies.py` · tests: `tests/chemistries/test_quasispecies.py`*
