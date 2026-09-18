# Naming game as an artificial chemistry

`naming-game-ac` · *De Beule, Hovig & Benson (book ref [212]); after Steels' naming games*

How a shared vocabulary appears, written as reactions. One species stands for a meaning, others for words and for the internal associations linking meaning to word. Reactions let a meaning produce a word, a word reinforce its association, and - crucially - a word plus a competing association convert that association to the winning one. The lateral inhibition that gives is enough for a population to converge on a single shared lexicon, reproducing what agent-based naming games show.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 16.3.3 |
| **refs** | [212], [804], [805], [806], [809] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants` |

## Molecules, reactions, reactor

**S — molecules** (explicit): one meaning molecule M, word molecules S_j and adaptor/code molecules C_j for j = 1..N_s

**R — reactions** (explicit, arity [2, 3]): M + C_j --ks--> S_j ; M + S_j --kc--> C_j ; M + S_j + C_j --ka--> 2 S_j ; M + S_j + C_k --kappa1--> M + S_j + C_j ; M + S_k + C_j --kappa2--> M + S_j + C_j (k != j)

**A — reactor**: ode
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("naming-game-ac", seed=1)
```

```
naming-game-ac: 7 species, 21 reactions, status=complete
provides: catalysts, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
M + C1 -> S1  [mass-action k=1.0]
M + S1 -> C1  [mass-action k=1.0]
M + S1 + C1 -> 2 S1  [mass-action k=1.0]
M + C2 -> S2  [mass-action k=1.0]
M + S2 -> C2  [mass-action k=1.0]
M + S2 + C2 -> 2 S2  [mass-action k=1.0]
M + C3 -> S3  [mass-action k=1.0]
M + S3 -> C3  [mass-action k=1.0]
… and 13 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N_s` | `int` | `3` | structural | number of competing names for the meaning <br>`1` … `20` |
| `ks` | `float` | `1.0` | kinetic | word production M + C_j -> S_j <br>≥ `0` |
| `kc` | `float` | `1.0` | kinetic | code production M + S_j -> C_j <br>≥ `0` |
| `ka` | `float` | `1.0` | kinetic | replication of a successful word <br>≥ `0` |
| `kappa1` | `float` | `1.0` | kinetic | replacement of a mismatching adaptor <br>≥ `0` |
| `kappa2` | `float` | `1.0` | kinetic | replacement of a mismatching word <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- lateral inhibition dynamics converging to a shared lexicon, reproducing Steels' agent-based results in ODE form

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- In the text of eqs. 16.5-16.6 the mismatch marks are lost (both read M + S_j + C_j -> M + S_j + C_j). Following the prose ('eliminate mismatching words (resp. adaptors) by replacing them with matching ones'), they are M + S_j + C_k -> M + S_j + C_j and M + S_k + C_j -> M + S_j + C_j for every k != j.
- Reactions 16.2-16.4 consume M as printed. The book gives no rate values or initial state; all rates default to 1 and no initial state is attached.

---

*Specification: `catalog/chemistries/naming-game-ac.yaml` · generator: `chemart/chemistries/naming_game_ac.py` · tests: `tests/chemistries/test_naming_game_ac.py`*
