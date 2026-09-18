# Kauffman NK model

`nk-landscape` · *Kauffman, 1993*

A dial for how rugged a fitness landscape is. A genotype's fitness is the average of per-gene contributions, and each gene's contribution depends on itself plus K others chosen at random. K = 0 gives a single smooth peak; raising K towards N makes the landscape progressively more uncorrelated and multi-peaked. It is less a chemistry than a fitness function other models plug into, exposed through the same interface so it can be used uniformly.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | analysis |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.4.1; appendix (NKlandscape.py) |
| **refs** | [447] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow` |

## Molecules, reactions, reactor

**S — molecules** (implicit): binary genotypes of length N

**R — reactions** (explicit, arity 1): fitness w = (1/N) sum_i w_i, where w_i is drawn at random as a function of gene i and its K epistatic partners; genotypes replicate at their fitness with point mutations: X_g -> 2 X_g ; X_g -> X_g + X_h

**A — reactor**: ode
 · *dilution:* constant-total dilution (replicator-mutator dynamics)

## What you get

```python
net = chemart.generate_network("nk-landscape", seed=1)
```

```
nk-landscape: 64 species, 448 reactions, status=complete
provides: catalysts, flow, rate-constants, stoichiometry, topology
seed: 1
extras: analysis
```

First reactions:

```
000000 -> 2 000000  [mass-action k=0.5369468618454898]
000000 -> 000000 + 100000  [mass-action k=0.005423705675206967]
000000 -> 000000 + 010000  [mass-action k=0.005423705675206967]
000000 -> 000000 + 001000  [mass-action k=0.005423705675206967]
000000 -> 000000 + 000100  [mass-action k=0.005423705675206967]
000000 -> 000000 + 000010  [mass-action k=0.005423705675206967]
000000 -> 000000 + 000001  [mass-action k=0.005423705675206967]
000001 -> 2 000001  [mass-action k=0.4125110307390986]
… and 440 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N` | `int` | `6` | structural | number of genes; the network has 2^N genotypes <br>`1` … `10` |
| `K` | `int` | `2` | structural | epistatic partners per gene: K=0 unimodal, growing K makes the landscape rugged, K=N-1 fully random <br>`0` … `9` · *range:* 0..N-1 |
| `topology` | `enum` | `adjacent` | structural | partners of gene i: the next K genes (circularly) or K random other genes <br>one of `adjacent`, `random` |
| `mu` | `float` | `0.01` | stochastic | per-gene mutation probability during replication <br>`0` … `1` |

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book defines only the landscape. As an analysis entry it is exposed through the network that evolves on it: every genotype replicates at its fitness, exact copies at W (1-mu)^N and single point mutants at W mu (1-mu)^(N-1) each (multiple mutations are omitted), under constant-total dilution. The landscape, epistatic partners, local optima and global optimum are in extras.analysis.
- Index convention (book figure 18.9): the table entry for gene i is indexed by the bits of its partners followed by its own bit, e.g. j = g_k g_i for K = 1.

## Notes

Mainly a fitness function other chemistries can plug in; the network form lets it go through the same interface.

---

*Specification: `catalog/chemistries/nk-landscape.yaml` · generator: `chemart/chemistries/nk_landscape.py` · tests: `tests/chemistries/test_nk_landscape.py`*
