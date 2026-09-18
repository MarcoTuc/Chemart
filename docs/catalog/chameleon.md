# Colored chameleon chemistry

`chameleon` · *Winkler, 2007 puzzle; formulated as an AC by Banzhaf & Yamamoto*

*Also known as:* *chameleon colony puzzle*

A puzzle in chemical clothing, and the smallest useful test case in the catalog. Three colours; when two chameleons of different colours meet, both take the third colour. The total count never changes, so the system just relaxes to equal thirds. Its real value is the conservation law: besides the obvious total, the differences between colour counts are conserved *modulo 3* - an integer invariant that any real-valued nullspace computation will miss entirely.

| | |
|---|---|
| **family** | core |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book` — implemented exactly as the book specifies |
| **book** | 2.5.1 |
| **refs** | [921], [922] |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `initial-state`, `mass-conservation` |

## Molecules, reactions, reactor

**S — molecules** (explicit): three colors r, g, b

**R — reactions** (explicit, arity 2): i + j -> 2k for all i != j != k over {r,g,b}; same-color collisions are elastic

**A — reactor**: well-stirred-multiset, ode
 · *dilution:* none needed - the rule is number-conserving (2 in, 2 out)

## What you get

```python
net = chemart.generate_network("chameleon", seed=1)
```

```
chameleon: 3 species, 3 reactions, status=complete
provides: initial-state, mass-conservation, rate-constants, stoichiometry, topology
seed: 1
extras: conservation
```

First reactions:

```
r + g -> 2 b  [mass-action k=1.0]
r + b -> 2 g  [mass-action k=1.0]
g + b -> 2 r  [mass-action k=1.0]
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `M` | `int` | `5400` | population | total number of chameleons (constant) <br>≥ `3` |
| `x0` | `list` | `[0.5, 0.0, 0.5]` | population | initial colour fractions of (r, g, b), summing to 1; counts are M * x0 rounded by largest remainder |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- conserved total count
- relaxation to the 1/3,1/3,1/3 fixed point
- the mod-3 invariant that makes the original puzzle unsolvable

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Rate constants are 1 for all three reactions, as implied by the book's rate equations 2.40-2.42.
- Default x0 follows figure 2.5 (only red and blue present at the start).

## Notes

extras.conservation carries the real law (1,1,1) and the two mod-3 invariants (r - g) and (g - b), which a real-valued nullspace computation will NOT find.

---

*Specification: `catalog/chemistries/chameleon.yaml` · generator: `chemart/chemistries/chameleon.py` · tests: `tests/chemistries/test_chameleon.py`*
