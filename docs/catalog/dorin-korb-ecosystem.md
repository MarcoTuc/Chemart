# Dorin & Korb virtual ecosystem chemistry

`dorin-korb-ecosystem` · *Dorin & Korb, 2007*

An ecosystem whose bookkeeping goes all the way down to atoms. Square atoms drift on a grid and bond covalently; each bond type stores or releases a definite energy, and four catalysts make or break specific bonds - one chlorophyll-like catalyst fixes sunlight into an A-B 'sugar' bond, others decompose organic and inorganic matter. Because the atom inventory is closed and energy is tracked exactly, a trophic structure of producers and decomposers appears from conservation alone rather than being modelled.

| | |
|---|---|
| **family** | spatial |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 8.2.3 |
| **refs** | [240], doi:10.1007/978-3-540-74913-4_11, https://bridges.monash.edu/articles/report/Building_Virtual_Ecosystems_from_Artificial_Chemistry/20365353 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `energies`, `mass-conservation`, `space`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (implicit): molecules = covalently bonded clusters of square virtual atoms from {A, B, C, O} plus the four catalyst atoms {K, EAB, ECC, EO}, identified by the canonical (isomorphism-invariant) bond graph

**R — reactions** (implicit, arity [1, 2]): One bond of table 1 made or broken per event, between two atoms that share a grid edge: A-B, C-C, A-O, B-O, C-O. Making a bond releases its bond energy and breaking it costs the same, so the negative bonds (A-B "sugar", C-C "biomass") are the ones that store energy. The paper's three named reactions follow: AO + BO -(K, sunlight)-> AB + 2 O (photosynthesis), O + AB -(EAB)-> A + BO + energy (respiration), C + C -(energy)-> C2 (biosynthesis).

**A — reactor**: lattice-2d
 · *dilution:* none: the atom inventory is closed, so every atom count is conserved exactly

## What you get

```python
net = chemart.generate_network("dorin-korb-ecosystem", seed=1)
```

```
dorin-korb-ecosystem: 26 species, 46 reactions, status=observed
provides: catalysts, energies, initial-state, mass-conservation, space, stoichiometry, topology
seed: 1
extras: analysis, chemistry, conservation, energies, events, final_state, space, species_encoding
```

First reactions:

```
B + O -> BO  (x93)
2 C -> C2  (x3)
A + O -> AO  (x90)
BO + C16KEAB_d33a -> B + O + C16KEAB_d33a  (x2)
C2 + C4OECC_5c53 -> 2 C + C4OECC_5c53  (x1)
BO + EO -> B + O + EO  (x12)
A + B + K -> AB + K  (x13)
AO + K -> A + O + K  (x47)
… and 38 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `width` | `int` | `18` | spatial | grid width in cells (one square atom per cell) <br>`4` … `400` · *range:* paper-scale grids are hundreds of cells wide |
| `height` | `int` | `18` | spatial | grid height in cells <br>`4` … `400` |
| `steps` | `int` | `150` | population | simulation time steps; one step is one movement phase followed by one reaction phase <br>`0` … `100000` · *range:* long runs use thousands of steps |
| `atoms` | `dict` | `{'A': 28, 'B': 28, 'O': 36, 'C': 24}` | structural | initial number of free atoms of each building block A, B, C, O |
| `catalysts` | `dict` | `{'K': 4, 'EAB': 4, 'ECC': 4, 'EO': 4}` | structural | initial number of free catalyst atoms: K chlorophyll, EAB sugar-breaker, ECC organic decomposer, EO inorganic decomposer |
| `structures` | `enum` | `all` | structural | which of the paper's organism bodies to seed: the figure 5 photosynthetic autotroph (a C box with chlorophyll and the sugar enzyme on opposite inner walls around a vacuole) and the figure 6 decomposer (a C-C body carrying ECC on a C-O spacer) <br>one of `none`, `photoautotroph`, `decomposer`, `all` |
| `light_amplitude` | `float` | `8.0` | thermodynamic | peak sunlight energy incident per atom per step; 0 switches the light off <br>≥ `0.0` |
| `light_period` | `float` | `24.0` | thermodynamic | period in steps of the sin function that governs incident sunlight <br>≥ `1.0` |
| `energy_high` | `int` | `8` | thermodynamic | the magnitude table 1 calls 'high' (the A-B sugar bond energy) <br>`2` … `1000` |
| `energy_low` | `int` | `1` | thermodynamic | the magnitude table 1 calls 'low' (C-C, A-O, B-O and C-O bond energies) <br>`1` … `999` |
| `p_low` | `float` | `0.01` | kinetic | per-step probability table 1 calls 'low' <br>`0.0` … `1.0` |
| `p_moderate` | `float` | `0.1` | kinetic | per-step probability table 1 calls 'moderate' <br>`0.0` … `1.0` |
| `p_high` | `float` | `0.9` | kinetic | per-step probability table 1 calls 'high', including every catalysed entry <br>`0.0` … `1.0` |
| `move_probability` | `float` | `0.8` | spatial | probability that a molecule attempts a one-square move in a random direction in a step <br>`0.0` … `1.0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- a trophic structure emerges from atom-level bookkeeping: chlorophyll fixes sunlight into A-B sugar bonds, the sugar enzyme releases that energy again, the released energy pays for C-C biosynthesis, and the decomposers break C-C and A-O/B-O bonds to return atoms to the inorganic pool
- sugar is made only where energy is available: with the light switched off no A-B bond is ever formed, since only sunlight can pay the high A-B bond energy
- the figure 5 autotroph runs the whole cycle on its own: with nothing else on the grid, the chlorophyll on one inner wall turns the A-O and B-O trapped in its vacuole into A-B sugar, and the enzyme on the opposite wall respires that sugar again
- atoms are never created or destroyed, and every joule is either stored in a bond, spent on one, or dissipated in the step it was released

## Sources

- Dorin, A. & Korb, K. B. (2007). Building Virtual Ecosystems from Artificial Chemistry. Advances in Artificial Life (ECAL 2007), LNCS 4648, 103-112. https://doi.org/10.1007/978-3-540-74913-4_11
- Dorin, A. & Korb, K. B. (2007). Building Virtual Ecosystems from Artificial Chemistry. Monash University Technical Report 2007/212 (the full version, used here): table 1 (the bond table), sections 2.1-2.2 (catalysts, energy neighbourhoods, sunlight, movement), sections 3.1.1-3.1.3 (photosynthesis, respiration, biosynthesis), figures 5-6 (autotroph and decomposer bodies) and the appendix (bonding rules, electron shells, known atoms and molecules). https://bridges.monash.edu/articles/report/Building_Virtual_Ecosystems_from_Artificial_Chemistry/20365353

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The paper gives no numbers: table 1 grades probabilities as low/moderate/high and bond energies as +/- low/high, and section 4 states that no run with all organism types had been performed, so nothing quantitative is published. The three levels are parameters (p_low 0.01, p_moderate 0.1, p_high 0.9) and so are the two magnitudes (energy_low 1, energy_high 8). Energies are integers so that the energy ledger balances exactly.
- Table 1 is implemented verbatim, including two entries the book's summary omits: ECC also catalyses *making* C-O bonds (the spacer of the figure 6 decomposer), and the catalysed break of A-O and B-O is credited to chlorophyll as well as to the inorganic enzyme.
- Table 1 names enzAO and enzBO separately, but section 2.1 counts four catalyst types, so one inorganic decomposer EO covers both A-O and B-O. The appendix names four enzymes EA-ED without mapping them to functions; catalysts are named for their function here (K chlorophyll, EAB, ECC, EO).
- Valences are derived, not assumed: appendix rule 2 (a bond needs a free outer-shell electron and a free slot) with shell capacities 2, 4, 8 gives A 1, B 1, O 1, C 4, K 2, enzymes 2, which is what the appendix's known molecules (A-B, A-O, B-O 'one only since O has run out of electrons', C chains, K-C-K) show.
- Only the five bonds of table 1 form or break. The appendix's K-C anchor bonds have no probabilities or energy, so bonds inside a seeded organism body are static anchors that never react ('any otherwise legal bond the Designer God doesn't want can simply be made to hurdle an absurd energy threshold').
- Energy pools are spatial, not per molecule: section 2.1 says released energy reaches 'any continuous atomic structure that contacts directly or indirectly (through intermediate neighbours) the reaction site', so the pool is the contact cluster of touching atoms. This is what lets respiration at an enzyme wall pay for biosynthesis elsewhere in the body. Atoms do not move during the reaction phase, so the clusters are fixed while reactions fire.
- Sunlight is 'incident on all atoms at a rate governed by a parameterized sin function': each cluster receives floor(light_amplitude * max(0, sin(2 pi t / light_period)) * atoms) per step, spendable only by chlorophyll-catalysed reactions. Light is spent before the bond-energy pool, and anything left in either pool at the end of a step is dissipated.
- The grid is a torus and the neighbourhood is von Neumann (atoms are squares that bond along shared edges). Movement is per molecule, rigid, one square, rejected on any collision.
- Species identity is the canonical bond graph, so lattice conformation (how a molecule is folded on the grid) is not part of it, while constitutional isomers with different bond graphs are different species. Ids are the atom formula for one- and two-atom molecules and formula plus a certificate hash beyond that.
- The network is the set of reactions that fired in one run (status observed, with counts). Catalysts appear on both sides whenever the catalysing atom belongs to a molecule that is not itself a reactant.
- Photosynthesis is three elementary bond events, not one: chlorophyll breaks A-O and breaks B-O, then makes A-B, so the freed A and B compete with the 'high' probability of re-forming A-O and B-O. The sugar yield therefore turns on a ratio table 1 leaves unquantified. In the soup, where A and B are plentiful, sugar is made steadily; an isolated figure 5 autotroph, whose vacuole holds one A and one B, makes it only in some runs (seed 5 at the default parameters is one).
- A free C can only bond where a reaction has just released energy into its contact cluster, so biomass and its decomposition concentrate around the seeded bodies and dense patches. The default therefore seeds both organism bodies (structures: all); with structures: none the soup is lively but builds almost no C-C biomass.

## Notes

Along with ToyChem and Bigan, one of the three energy-carrying chemistries in the book. The paper is a feasibility argument rather than a results paper: it specifies the chemistry and hand-designs organism bodies, and states that a run with all organism types present had not yet been done.

---

*Specification: `catalog/chemistries/dorin-korb-ecosystem.yaml` · generator: `chemart/chemistries/dorin_korb_ecosystem.py` · tests: `tests/chemistries/test_dorin_korb_ecosystem.py`*
