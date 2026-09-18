# Replicator equation (evolutionary game dynamics)

`replicator-equation`

Evolutionary game theory as reactions. Strategies meet in pairs and the payoff matrix decides what happens: a positive entry means the first strategy replicates in the encounter, a negative one that it is removed, with the opponent surviving either way. Evolutionarily stable strategies are the rest points. Note the shape of the reaction - the opponent appears on both sides as a catalyst, so a net stoichiometric matrix would erase the very interaction the model is about.

| | |
|---|---|
| **family** | evolutionary-dynamics |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 7.2.8 |
| **refs** | [621], [563], [564], [792] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): n strategies X_1..X_n

**R — reactions** (explicit, arity 2): X_i + X_j --a_ij--> 2X_i + X_j (a_ij > 0, cooperative/catalytic) ; X_i + X_j --(-a_ij)--> X_j (a_ij < 0, defective) ; elastic for a_ij = 0

**A — reactor**: ode, ssa
 · *dilution:* phi = sum_ij a_ij x_i x_j

## What you get

```python
net = chemart.generate_network("replicator-equation", seed=1)
```

```
replicator-equation: 3 species, 6 reactions, status=complete
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X1 + X2 -> X2  [mass-action k=1.0]
X1 + X3 -> 2 X1 + X3  [mass-action k=1.0]
X2 + X1 -> 2 X2 + X1  [mass-action k=1.0]
X2 + X3 -> X3  [mass-action k=1.0]
X3 + X1 -> X1  [mass-action k=1.0]
X3 + X2 -> 2 X3 + X2  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `A` | `list` | `[[0.0, -1.0, 1.0], [1.0, 0.0, -1.0], [-1.0, 1.0…` | selection | n x n payoff matrix a_ij; asymmetry a_ij != a_ji encodes the game. Default: rock-paper-scissors |
| `x0` | `list` | `` | population | initial frequencies; empty means uniform 1/n |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- evolutionarily stable strategies (ESS)
- known to be equivalent to Lotka-Volterra with n+1 strategies

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- n is the size of A rather than a separate parameter. The book gives no game; the default is rock-paper-scissors. The dilution flux phi is outflow = 'constant-total', which equals eq. 7.24.
- A diagonal term a_ii is the bimolecular self-interaction 2X_i -> 3X_i (a_ii > 0) or 2X_i -> X_i (a_ii < 0).

## Notes

Note the catalytic form: X_j is on both sides. Net stoichiometry alone would erase it.

---

*Specification: `catalog/chemistries/replicator-equation.yaml` · generator: `chemart/chemistries/replicator_equation.py` · tests: `tests/chemistries/test_replicator_equation.py`*
