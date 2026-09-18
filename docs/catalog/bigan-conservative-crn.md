# Conservative random chemical reaction networks

`bigan-conservative-crn` · *Bigan, Steyaert & Douady, 2013*

A random reaction network built to be thermodynamically honest. Each abstract chemical gets a random free energy of formation; reactions are drawn as association, dissociation and transformation steps, and every backward rate is then computed from the forward one so detailed balance holds. The result is a network you can drive out of equilibrium with a nutrient influx and still trust. It is the easiest place to see what the thermodynamic tier costs: consistency is a constraint between forward and reverse rates, not extra data.

| | |
|---|---|
| **family** | network |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.3.2 |
| **refs** | [117], [756], [956] |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `rate-law`, `energies`, `thermodynamic-consistency`, `mass-conservation`, `flow`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): N abstract chemicals A_0..A_(N-1), each with a random Gibbs free energy of formation

**R — reactions** (implicit, arity [1, 2]): association A_i + A_j -> A_k ; dissociation A_i -> A_j + A_k ; transformation A_i -> A_j ; each kept as a forward/reverse pair

**A — reactor**: ode
 · *dilution:* external nutrient inflow drives the system out of equilibrium

## What you get

```python
net = chemart.generate_network("bigan-conservative-crn", seed=1)
```

```
bigan-conservative-crn: 10 species, 46 reactions, status=complete
provides: energies, flow, initial-state, mass-conservation, rate-constants, stoichiometry, topology
seed: 1
extras: conservation, energies, reaction_pairs
```

First reactions:

```
A1 + A3 -> A5  [mass-action k=10000.0]
A5 -> A1 + A3  [mass-action k=5852.229125020699]
A2 + A5 -> A0  [mass-action k=10000.0]
A0 -> A2 + A5  [mass-action k=62.318124632246665]
A9 -> A1  [mass-action k=100.0]
A1 -> A9  [mass-action k=8.488690652415743]
2 A1 -> A5  [mass-action k=10000.0]
A5 -> 2 A1  [mass-action k=2865.976715541064]
… and 38 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `N` | `int` | `10` | structural | number of chemical species (the paper uses 10, and 20 for checks) <br>`2` … `20` |
| `max_reactions` | `int` | `` | structural | stop after this many direct reactions; 0 builds a maximum-sized network (17-96, mean 42 direct reactions for N = 10 in the paper) <br>≥ `0` |
| `kinetics` | `enum` | `mass-action` | kinetic | saturating kinetics f = k prod [A]/(1 + [A]/K) models molecular crowding <br>one of `mass-action`, `saturating` |
| `k_avg_mono` | `float` | `100.0` | kinetic | geometric mean of monomolecular forward constants (1/s): k = k_avg 10^U(-s/2, s/2) <br>≥ `0` |
| `k_avg_bi` | `float` | `10000.0` | kinetic | geometric mean of bimolecular forward constants (1/(M s)) <br>≥ `0` |
| `s` | `float` | `` | kinetic | spread of forward rate constants, in orders of magnitude <br>≥ `0` |
| `K_avg` | `float` | `0.01` | kinetic | geometric mean saturation concentration (M) for saturating kinetics <br>≥ `0` |
| `p` | `float` | `` | kinetic | spread of saturation concentrations, in orders of magnitude <br>≥ `0` |
| `G_max` | `float` | `15.0` | thermodynamic | formation free energies are drawn as G_i/RT ~ U(0, G_max) <br>≥ `0` |
| `nutrient` | `int` | `5` | population | index of the species receiving the external flux (the paper uses A5); -1 for a closed system <br>≥ `-1` |
| `nutrient_flux` | `float` | `1.0` | population | external nutrient flux f_nu (M/s) <br>≥ `0` |
| `initial_concentration` | `float` | `0.001` | population | initial concentration of every species (the paper uses 1 mM) <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- at equilibrium one chemical dominates as a function of density (saturating kinetics)
- under nutrient flux all but one species reach constant concentration while one grows unbounded - a 'directed transformation machine' producing a putative membrane precursor

## Sources

- Bigan, E., Steyaert, J.-M. & Douady, S. (2013). Properties of random complex chemical reaction networks and their relevance to biological toy models. arXiv:1303.7439, sections 2-3

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Topology follows sec. 2.1: exchange reactions are excluded, candidates are tried in random order, and a candidate is kept while a mass vector m >= 1 with m^T S = 0 exists (checked by linear programming). Once that vector is unique, a candidate is kept iff it is orthogonal to it.
- Kinetics follow sec. 2.2 with the paper's values: k_avg_mono = 1e2 /s, k_avg_bi = 1e4 /(M s), K_avg = 1e-2 M, G/RT ~ U(0, 15), s = p = 0. The forward direction is the one with decreasing free energy, and the reverse constant is k_forward exp(-|dG|/RT) with c0 = 1 M. The paper's printed reverse formulas for mixed molecularity differ only by the c0 factor, which is 1 here.
- Forward/reverse reactions are adjacent, with their indices listed in extras.reaction_pairs. extras.energies holds G_i/RT, and extras.conservation holds the admissible mass vector (integers when rational).
- The old density parameter is not part of the network (the paper sweeps initial density to compute equilibria) and is dropped. The nutrient flux is the network inflow.

## Notes

The easiest thermodynamic entry: random generation with a conservativity check and backward rates fixed by detailed balance.

---

*Specification: `catalog/chemistries/bigan-conservative-crn.yaml` · generator: `chemart/chemistries/bigan_conservative_crn.py` · tests: `tests/chemistries/test_bigan_conservative_crn.py`*
