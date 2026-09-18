# Computing with chemical organizations

`organization-computing` · *Matsumaru, Speroni di Fenizio, Centler & Dittrich, 2006-2011*

Computing where the answer is which species survive. Each Boolean variable gets two species, one for true and one for false, and information is carried by presence or absence rather than concentration. Gate reactions produce output species from input species, and an annihilation reaction removes any variable that holds both values at once. Run it and the set of species that persists - the organisation - is the solution. The whole computation is read off a lattice of self-maintaining sets.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 17.3.3 |
| **refs** | [550], [551], [552], [553], [554], [555], [506] |
| **provides** | `topology`, `stoichiometry`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): two species per Boolean variable: x0 for x=0 and x1 for x=1; information is encoded by PRESENCE/ABSENCE, not concentration

**R — reactions** (explicit, arity variable): XOR gate: a0+b0 -> c0 ; a0+b1 -> c1 ; a1+b0 -> c1 ; a1+b1 -> c0, plus annihilation x0 + x1 -> (nothing) for each variable. Maximal independent set: s0_j + s1_j -> ; s0_j + s0_k + ... + s0_l -> n_i s1_i (over the n_i neighbours of v_i) ; s1_j -> s0_i for all (v_j, v_i) in E.

**A — reactor**: ode, ssa
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("organization-computing", seed=1)
```

```
organization-computing: 6 species, 7 reactions, status=complete
provides: initial-state, stoichiometry, topology
seed: 1
```

First reactions:

```
a0 + b0 -> c0
a0 + b1 -> c1
a1 + b0 -> c1
a1 + b1 -> c0
a0 + a1 -> ∅
b0 + b1 -> ∅
c0 + c1 -> ∅
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `problem` | `enum` | `xor` | structural | which chemical program of the book to generate <br>one of `xor`, `maximal-independent-set` |
| `a` | `bool` | `` | population | XOR input a; its species a0 or a1 is put in the initial state |
| `b` | `bool` | `True` | population | XOR input b; its species b0 or b1 is put in the initial state |
| `graph` | `list` | `[[0, 1], [1, 2], [2, 3]]` | structural | undirected edges [u, v] of the maximal-independent-set instance |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the surviving organisation IS the answer; 15 organisations in the XOR Hasse diagram
- flip-flops, oscillators, NAND chains
- robustness: perturbing the system makes it fall back into the solution organisation, an attractor

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- No rates: the book's analysis is organisational (which species survive), not kinetic.
- The book gives no MIS instance; the default graph is a 4-vertex path. Edges are undirected, so eq. 17.7 applies in both directions. Inputs a and b only affect the XOR program.

## Notes

Its shortcoming (one species per variable state) is precisely why the output schema needs to carry the ORGANISATION LATTICE as a first-class analysis product, not just the matrix.

---

*Specification: `catalog/chemistries/organization-computing.yaml` · generator: `chemart/chemistries/organization_computing.py` · tests: `tests/chemistries/test_organization_computing.py`*
