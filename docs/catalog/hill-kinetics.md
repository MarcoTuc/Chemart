# Hill kinetics (cooperative binding)

`hill-kinetics`

Cooperative binding, and the standard way a smooth input becomes a switch. n copies of a protein bind a gene together, so the fraction of bound genes follows P^n/(K^n + P^n) - a curve that sharpens into a step as n grows. It is an idealisation on purpose: n molecules never actually arrive simultaneously, and with a single gene copy and a handful of proteins a differential equation is the wrong tool. Useful precisely because it compresses a mechanism into one parameter.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.2.2; used by the Repressilator 19.3.2 |
| **refs** | [741], [916] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-law`, `rate-constants`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): gene G, protein P, complex C, gene product X

**R — reactions** (explicit, arity n+1): G + n P <-kf/kr-> C ; activation C -> C + X ; repression G -> G + X ; lumped rates H1 = P^n/(K^n + P^n), H2 = K^n/(K^n + P^n)

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("hill-kinetics", seed=1)
```

```
hill-kinetics: 4 species, 3 reactions, status=complete
provides: catalysts, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
G + 2 P -> C  [mass-action k=1.0]
C -> G + 2 P  [mass-action k=1.0]
C -> C + X  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `form` | `enum` | `elementary` | structural | elementary binding reactions (eqs. 18.5, 18.12, 18.15) or a single production reaction with a Hill rate law (eqs. 18.14, 18.16) <br>one of `elementary`, `lumped` |
| `regulation` | `enum` | `activation` | structural | P activates (bound gene expresses) or represses (free gene expresses) X <br>one of `activation`, `repression` |
| `n` | `int` | `2` | kinetic | Hill coefficient; better read as the DEGREE OF COOPERATIVITY than as a literal molecule count. n=1 gives no sigmoid, so no switching <br>`1` … `10` |
| `K` | `float` | `1.0` | kinetic | K^n = K_d = kr/kf, the dissociation constant <br>≥ `0` |
| `kf` | `float` | `1.0` | kinetic | binding rate; the unbinding rate is kf K^n <br>≥ `0` |
| `k_expr` | `float` | `1.0` | kinetic | expression rate of X (k1 or k2 in the book) <br>≥ `0` |
| `G0` | `float` | `1.0` | population | total gene = free + complexed <br>≥ `0` |
| `P0` | `float` | `1.0` | population | initial protein <br>≥ `0` |

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book derives both forms; `form` lets the user pick. Defaults follow figure 18.5 (G(0) = P(0) = 1) with unit rates. In the lumped form the rate dict names the regulator species ('regulator': 'P') and 'mode'.

## Notes

The book flags the idealisation explicitly: n molecules never bind simultaneously, and with one gene copy and few proteins an ODE is a poor approximation.

---

*Specification: `catalog/chemistries/hill-kinetics.yaml` · generator: `chemart/chemistries/hill_kinetics.py` · tests: `tests/chemistries/test_hill_kinetics.py`*
