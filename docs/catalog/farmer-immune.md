# Bitstring immune system model (idiotypic network)

`farmer-immune` · *Farmer, Packard & Perelson, 1986*

An immune system as a reaction network, and one of the few chemistries that fixes its own rate constants. Antibodies are bitstrings with two sites; when one antibody's paratope matches another's epitope, the recogniser is replicated and the recognised destroyed. The crucial move is that the reaction rate is computed from the *degree* of complementarity between the strings, over all alignments - structure determines kinetics. Species below a concentration threshold are periodically culled and new ones injected, so the species set turns over.

| | |
|---|---|
| **family** | bio-inspired |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 11.2.2 |
| **refs** | [261], [262], [214], [287], [288] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): antibodies X_i as binary strings with two binding sites (epitope e, paratope p); antigens Y_j with an epitope

**R — reactions** (explicit, arity [1, 2]): Recognition of an epitope by a paratope replicates the recogniser and destroys the recognised. Matching specificity: m_ij = sum_k G( sum_n [ e_i(n+k) XOR p_j(n) ] - s + 1 ), G(x)=x for x>0 else 0, over all alignments -l_p < k < l_e.

**A — reactor**: ode, ssa
 · *dilution:* metadynamics: periodically drop species below a concentration threshold and inject new ones

## What you get

```python
net = chemart.generate_network("farmer-immune", seed=1)
```

```
farmer-immune: 12 species, 88 reactions, status=complete
provides: catalysts, rate-constants, stoichiometry, topology
seed: 1
```

First reactions:

```
X1 + X2 -> 2 X1 + X2  [mass-action k=1.0]
X1 + X4 -> X4  [mass-action k=1.0]
X1 + X5 -> X5  [mass-action k=1.0]
X1 + X6 -> 2 X1 + X6  [mass-action k=1.0]
X1 + X7 -> 2 X1 + X7  [mass-action k=1.0]
X1 + X8 -> X8  [mass-action k=1.0]
X1 + X10 -> 2 X1 + X10  [mass-action k=1.0]
X2 + X1 -> X1  [mass-action k=1.0]
… and 80 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N` | `int` | `10` | population | number of antibody species <br>`1` … `100` |
| `M` | `int` | `2` | population | number of foreign antigens <br>`0` … `20` |
| `l_e` | `int` | `8` | structural | epitope length (bits) <br>`1` … `64` |
| `l_p` | `int` | `8` | structural | paratope length (bits) <br>`1` … `64` |
| `s` | `int` | `6` | structural | matching threshold: an alignment counts when at least s bits are complementary <br>≥ `1` |
| `c` | `float` | `1.0` | kinetic | weight of the interactive part vs. spontaneous decay <br>≥ `0` |
| `k1` | `float` | `1.0` | kinetic | suppression vs. stimulation weight, antibody-antibody <br>≥ `0` |
| `k2` | `float` | `0.5` | kinetic | spontaneous antibody decay rate <br>≥ `0` |
| `k3` | `float` | `1.0` | kinetic | suppression weight, antibody-antigen <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- clonal selection
- idiotypic network memory
- constructive dynamics via metadynamics

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Strings are random bits from the seed. Each term of eqs. 11.5-11.6 is one mass-action reaction whose rate the matching strengths prescribe: stimulation X_i + X_j -> 2X_i + X_j (c m_ji), suppression X_i + X_j -> X_j (c k1 m_ij), antigen stimulation X_i + Y_j -> 2X_i + Y_j (c m), antigen removal Y_j + X_i -> X_i (k3 m), decay X_i -> (k2). A self-interaction i = j uses 2X_i on the left.
- The book says strings bind when 's > min(l_e, l_p)', which contradicts its own definition of G. Here an alignment contributes when at least s bits are complementary, which is what eq. 11.4 computes.
- The book gives no parameter values; defaults are N = 10, M = 2, 8-bit sites, s = 6, c = k1 = k3 = 1, k2 = 0.5. Metadynamics (mutation and replacement) acts between network snapshots and is not generated; mutation_ops is dropped.

## Notes

One of very few entries where the chemistry itself PRESCRIBES rate coefficients (from the structural match m_ij) rather than leaving them free.

---

*Specification: `catalog/chemistries/farmer-immune.yaml` · generator: `chemart/chemistries/farmer_immune.py` · tests: `tests/chemistries/test_farmer_immune.py`*
