# GARD (Graded Autocatalysis Replication Domain)

`gard` · *Segre, Lancet et al., 1998-2001*

*Also known as:* *Lipid World model*, *polymer GARD*, *EE-GARD*, *spatial GARD*

Heredity without a sequence. An assembly of lipid types grows by absorbing lipids from outside, and lipids already inside catalyse the entry of particular types, so the assembly's *composition* biases which lipids it takes up next. Split it in two and each half tends to reconstruct the same composition - a composome. The information is in the multiset, not in an ordered polymer. It also shows cleanly that catalysis alone cannot shift an equilibrium: without the external supply, nothing happens.

| | |
|---|---|
| **family** | origin-of-life |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 6.2.4 |
| **refs** | [764], [765], [760], [761], [762], [763], [771], [772], [773], [490], [726], [884] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `compartments`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): lipid types L_1..L_N outside (buffered) and A_1..A_N inside one assembly; the assembly's COMPOSITION (a 'compositional genome' / composome) carries the information

**R — reactions** (explicit, arity [1, 2]): L_i --k_f--> A_i ; A_i --k_b--> L_i ; L_i + A_j --beta_ij--> A_i + A_j (catalysis by lipids already in the assembly)

**A — reactor**: ssa
 · *dilution:* outside lipids buffered; the assembly splits in two when it reaches N_max molecules

## What you get

```python
net = chemart.generate_network("gard", seed=1)
```

```
gard: 200 species, 10200 reactions, status=complete
provides: catalysts, compartments, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: buffered, compartments
```

First reactions:

```
L1 -> A1  [mass-action k=0.01 crowding_capacity=100.0]
A1 -> L1  [mass-action k=0.0001]
L2 -> A2  [mass-action k=0.01 crowding_capacity=100.0]
A2 -> L2  [mass-action k=0.0001]
L3 -> A3  [mass-action k=0.01 crowding_capacity=100.0]
A3 -> L3  [mass-action k=0.0001]
L4 -> A4  [mass-action k=0.01 crowding_capacity=100.0]
A4 -> L4  [mass-action k=0.0001]
… and 10192 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N_G` | `int` | `100` | structural | number of lipid types; the network has N_G^2 catalysed joining reactions <br>`2` … `200` |
| `k_f` | `float` | `0.01` | kinetic | basal joining rate constant <br>≥ `0` |
| `k_b` | `float` | `0.0001` | kinetic | basal leaving rate constant <br>≥ `0` |
| `rho` | `float` | `0.01` | population | buffered outside concentration of every lipid type <br>≥ `0` |
| `beta_mu` | `float` | `-4.0` | kinetic | mean of ln(beta_ij); beta is the receptor-affinity-distribution rate enhancement |
| `beta_sigma` | `float` | `4.0` | kinetic | standard deviation of ln(beta_ij) <br>≥ `0` |
| `N_max` | `int` | `100` | population | assembly size at which fission occurs; also the crowding capacity of joining <br>≥ `2` |
| `catalysed_leaving` | `bool` | `` | structural | also accelerate leaving by beta_ij (the book's 'in both directions' reading), adding A_i + A_j -> L_i + A_j at k_b beta_ij |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- homeostatically stable composomes (quasi-stationary states)
- composome succession as mutations accumulate
- spontaneous symmetry breaking akin to homochirality (spatial variant)
- known evolvability limits: replication fidelity is too low

## Sources

- Predicting species emergence in simulated complex pre-biotic networks. PLoS ONE 13(2): e0192871 (2018), doi:10.1371/journal.pone.0192871 - GARD rate equation, ln(beta) ~ Normal(-4, 4), k_f = 1e-2, k_b = 1e-4, rho = 1e-2, N_max = N_G = 100
- Segre, D., Ben-Eli, D. & Lancet, D. (2000). Compositional genomes: prebiotic information transfer in mutually catalytic noncovalent assemblies. PNAS 97:4112-4117 (original model; full text not accessible, parameters taken from the 2018 paper)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives the reaction scheme but no numbers. The rate equation and all parameter values are from the 2018 PLoS ONE GARD paper; each term becomes one mass-action reaction with outside lipids buffered at rho (checked in tests).
- That equation accelerates only joining, while the book says catalysis acts 'in both directions'. The default follows the paper's equation; catalysed_leaving = true adds the book's reverse catalysis.
- The crowding factor (1 - N/N_max) is not mass action. It is recorded on every joining rate as crowding_capacity, and fission at N_max is described in extras.compartments. beta is sampled from the seed, so the network is a random instance.
- The polymer, EE and spatial variants are not generated (the former 'variant' parameter is dropped); split_size is N_max; k_i, k_minus_i and beta are replaced by the distribution parameters above. thermodynamic-consistency is not claimed, because the default catalyses only one direction.

## Notes

A rare case where information lives in a MULTISET COMPOSITION rather than a sequence. Also a clean demonstration that catalysis alone cannot shift equilibrium - the flow term is mandatory.

---

*Specification: `catalog/chemistries/gard.yaml` · generator: `chemart/chemistries/gard.py` · tests: `tests/chemistries/test_gard.py`*
