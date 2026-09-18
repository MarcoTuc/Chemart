# French flag model (positional information)

`french-flag` · *Wolpert, 1969*

The textbook picture of positional information, kept as a benchmark rather than as a mechanism. Cells sit along a morphogen gradient and read their position by which concentration thresholds it exceeds, adopting one of three fates - so a smooth gradient becomes three sharp stripes. It is the standard target pattern for evo-devo experiments, and the standard criticism is that it says nothing about where the gradient came from or how the thresholds are set.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.6 |
| **refs** | [927], [928], [929], [930], [456], [34], [265], [715] |
| **provides** | `topology`, `stoichiometry`, `space` |

## Molecules, reactions, reactor

**S — molecules** (explicit): undifferentiated cells U at positions along a morphogen gradient, and their fates (blue, white, red)

**R — reactions** (explicit, arity 1): U_x -> fate_x, where the fate is chosen by the highest morphogen threshold that the concentration at x reaches

**A — reactor**: lattice-2d, continuous-space
 · *dilution:* the gradient is usually imposed, not produced

## What you get

```python
net = chemart.generate_network("french-flag", seed=1)
```

```
french-flag: 60 species, 30 reactions, status=complete
provides: space, stoichiometry, topology
seed: 1
extras: space
```

First reactions:

```
U_0 -> blue_0
U_1 -> blue_1
U_2 -> blue_2
U_3 -> blue_3
U_4 -> blue_4
U_5 -> blue_5
U_6 -> blue_6
U_7 -> blue_7
… and 22 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n_cells` | `int` | `30` | structural | cells along the gradient axis (and along the other axis in 2D) <br>`3` … `200` |
| `thresholds` | `list` | `[0.667, 0.333]` | selection | morphogen concentration thresholds; k thresholds give k+1 fates (2 thresholds: blue above the high one, white between, red below) |
| `gradient` | `enum` | `linear` | spatial | imposed gradient shape along x in [0, 1]: c = 1 - x or c = exp(-x / decay_length) <br>one of `linear`, `exponential` |
| `decay_length` | `float` | `0.3` | spatial | decay length of the exponential gradient <br>≥ `0.001` |
| `dimensions` | `enum` | `1` | spatial | 1D row of cells (Wolpert's original) or its 2D extension with the gradient along x <br>one of `1`, `2` |

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The model has no reaction network of its own. Chemart represents it as one unimolecular differentiation reaction per cell, U_x -> fate_x, with the imposed morphogen concentration per cell in extras.space.
- The book gives no numbers: thresholds default to thirds of a linear gradient, and n_segments is implied by the number of thresholds. A self-generated gradient needs a gradient-producing chemistry and is not offered.

## Notes

Kept as a BENCHMARK: the standard target pattern for evo-devo experiments. The book's criticism is that the model says nothing about how the gradient is built or maintained.

---

*Specification: `catalog/chemistries/french-flag.yaml` · generator: `chemart/chemistries/french_flag.py` · tests: `tests/chemistries/test_french_flag.py`*
