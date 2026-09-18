# Ono & Ikegami autopoietic protocells

`ono-ikegami-protocell` · *Ono & Ikegami, 1999-2003*

Protocells from nothing but hydrophobicity. Five particle types sit on a lattice - an autocatalyst, membrane particles, food, waste and water - with an autocatalyst that reproduces and also makes membrane material, and decay that returns everything to waste and back to food. Membrane particles repel water and so cluster; give them an orientation and the clusters become filaments that close into compartments. The enclosure retains the autocatalyst while food and waste pass through. There is no genome and no replication, only self-maintenance.

| | |
|---|---|
| **family** | origin-of-life |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 6.3.2 |
| **refs** | [637], [639], [641], [538] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `space`, `energies`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): five species A (autocatalyst, hydrophilic), M (membrane, hydrophobic, isotropic M_i or anisotropic M_a carrying one of the six lattice orientations), X (food, neutral), Y (waste, neutral), W (water, hydrophilic), one particle per cell of a hexagonal torus

**R — reactions** (explicit, arity [1, 2]): A + X -> 2 A ; A + X -> A + M ; A -> Y, M -> Y, X -> Y (decay) ; Y -> X (recycling by an external energy source)

**A — reactor**: lattice-2d
 · *dilution:* none: the lattice is closed and every cell holds exactly one particle; the external energy source enters as the recycling reaction Y -> X, whose rate X_supply is the food supply

## What you get

```python
net = chemart.generate_network("ono-ikegami-protocell", seed=1)
```

```
ono-ikegami-protocell: 5 species, 6 reactions, status=observed
provides: catalysts, energies, initial-state, space, stoichiometry, topology
seed: 1
extras: analysis, energies, space
```

First reactions:

```
A + X -> 2 A  (x158)
A + X -> A + M_a  (x160)
A -> Y  (x106)
M_a -> Y  (x61)
X -> Y  (x46)
Y -> X  (x206)
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `mode` | `enum` | `spatial` | structural | 'reactions' returns the four-reaction network on its own; 'spatial' runs the lattice and returns the reactions that fired with their counts, plus the lattice, the energies and the protocell measurements <br>one of `spatial`, `reactions` |
| `membrane` | `enum` | `anisotropic` | structural | whether membrane particles are the anisotropic M_a or the isotropic M_i; the book's central result is that only M_a sustains protocells at low food supply <br>one of `anisotropic`, `isotropic` |
| `width` | `int` | `24` | spatial | lattice width in cells (axial r); the papers use a few hundred <br>`6` … `200` |
| `height` | `int` | `24` | spatial | lattice height in cells (axial q) <br>`6` … `200` |
| `steps` | `int` | `100` | population | lattice sweeps; each sweep does six Metropolis exchange passes, one rotation pass and one chemistry pass <br>`1` … `20000` |
| `relaxation` | `int` | `12` | spatial | Metropolis exchange passes per chemistry pass; this is where the paper's separation between the mobility rates (7e-3) and the reaction rates (1e-4) enters, so particles demix long before they react <br>`1` … `200` |
| `initial` | `enum` | `random` | population | 'random' is a well-mixed start (the papers' homogeneous experiment); 'cell' prepares a membrane ring around autocatalyst and food (their cell-like experiment) <br>one of `random`, `cell` |
| `cell_radius` | `int` | `5` | spatial | radius in cells of the prepared membrane ring when initial = cell <br>`1` … `60` |
| `A_fraction` | `float` | `0.08` | population | fraction of cells seeded with autocatalyst A when initial = random <br>`0` … `1` |
| `X_fraction` | `float` | `0.3` | population | fraction of cells seeded with food X; the rest of the lattice is water W <br>`0` … `1` |
| `M_fraction` | `float` | `` | population | fraction of cells seeded with membrane M when initial = random; non-zero starts from dispersed membrane particles, so their self-assembly can be watched on its own <br>`0` … `1` |
| `P_A` | `float` | `0.03` | kinetic | probability per sweep and per neighbouring A that a food particle replicates it (A + X -> 2 A); the paper's reproduction rate of A <br>`0` … `1` |
| `P_M` | `float` | `0.03` | kinetic | probability per sweep and per neighbouring A that a food particle becomes membrane (A + X -> A + M); the paper's production rate of M <br>`0` … `1` |
| `P_decay` | `float` | `0.01` | kinetic | probability per sweep that an A, M or X particle decays to waste Y; water never decays <br>`0` … `1` |
| `X_supply` | `float` | `0.3` | kinetic | probability per sweep that waste is recycled to food (Y -> X) by the external energy source; the decisive resource-supply parameter, low values being the regime where anisotropy pays off <br>`0` … `1` |
| `anisotropy` | `float` | `2.4` | spatial | strength a of the directional part of the M_a repulsion field, F[k, o] = 1 + a cos(2 * 60deg * (k - o)): repulsion is strongest along the +-o axis and weakest on the four flanks, where it turns attractive once a > 2, so a membrane particle wants membrane neighbours along its axis and water on both flanks; 0 makes M_a behave like M_i <br>`0` … `3` |
| `repulsion` | `float` | `1.0` | thermodynamic | energy of one hydrophilic/hydrophobic neighbour pair, in units of the temperature <br>≥ `0` |
| `neutral_coupling` | `float` | `0.25` | thermodynamic | energy of a neutral/hydrophobic pair as a fraction of `repulsion`; this weak coupling is what lets X and Y cross a membrane that holds A back <br>`0` … `1` |
| `temperature` | `float` | `1.0` | thermodynamic | Metropolis temperature T in min(1, exp(-dE/T)) <br>≥ `0.01` |
| `mobility_ratio` | `float` | `1.0` | spatial | the paper's m: exchanges involving A are attempted with probability 1/m, making the autocatalyst less mobile than membrane and food; the 1D paper needs m = 5 for a cell to divide <br>≥ `1` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- membrane particles demix from water into clusters; with anisotropic M_a the clusters are filaments that close into protocells
- a closed membrane is osmotically permeable to the neutral X and Y but holds the hydrophilic A back, so a metabolism persists inside
- a protocell grows and divides while its internal metabolism continues
- isotropic M_i gives clusters that shrink and die out, anisotropic M_a gives dividing protocells; the advantage is clearest at low X supply (book figure 6.14)
- without recycling of waste into food the cell starves and collapses
- 3D: parallel membranes with tubular connections, and intertwined globular filaments (not implemented)

## Sources

- Ono, N. & Ikegami, T. (2000). Self-maintenance and self-reproduction in an abstract cell model. J. theor. Biol. 206:243-253, doi:10.1006/jtbi.2000.2121. Full text from the authors' lab archive, http://sacral.c.u-tokyo.ac.jp/pdf/ono_jtheorbiol_2000.pdf (via the Wayback Machine). Appendix A gives the reaction and mobility rules and every constant; appendix B the mean-field equations; figures 5, 7, 8 and 10 the phase diagram and the parameter values of each regime.
- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, MIT Press, section 6.3.2 and figures 6.13, 6.14 - the five-species hexagonal model, the isotropic/anisotropic repulsion fields, and the low-supply result.

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book describes the two-dimensional hexagonal model of [639]/[641]; those papers, the BioSystems study [637] and the 3D extension [538] are all paywalled and were not obtainable (no open copy per OpenAlex/Unpaywall/Semantic Scholar, and the Internet Archive was offline during this work). The accessible primary source is the authors' 1D predecessor, Ono & Ikegami (2000), whose appendix A supplies the lattice algorithm reconstructed here.
- That 1D model runs a different reaction scheme - it inserts an enzyme E (A + A -> A + E, E + R -> E + A, E + R -> E + M) and its resource/waste are R and W, whereas the book's 2D model has A act directly and names the five species A, M, X, Y, W. The implemented reaction set is the book's; the 1D constants (P_E = 0.5e-6, P_A' = 2e-6, P_R = P_W = 100e-6, P_D = 7e-3, P_R0 = 5e-3, P_R1 = 0.5 P_R0, m = 1 or 5) therefore do NOT apply to it and are recorded, with the phase-diagram points of figures 4, 7 and 10, in extras.reference_rates rather than used as rate constants.
- Consequently every reaction carries rate = None: no accessible source publishes a rate for this scheme. P_A, P_M, P_decay and X_supply are simulation probabilities per sweep whose defaults were chosen so the zero-argument run is small and fast, and are not published values.
- Ono & Ikegami coarse-grain space into sites holding 100 particles each. Here each cell of the hexagonal torus holds exactly one particle, which is what lets M_a carry the orientation the book's figure 6.13 requires, and particles move by Kawasaki exchange under the Metropolis rule instead of the paper's probability-proportional hopping. Both conserve the occupancy of every cell, which is the property the paper's swap rule exists to enforce.
- The anisotropic field is nematic, F[k, o] = 1 + a cos(2 * 60deg * (k - o)), with mean 1 over the six directions, so M_i and M_a carry the same total repulsion and differ only in how it is distributed: strongest along the +-o axis, weakest on the four flanks. The book's figure 6.13 shows a darker (more repulsive) direction around M_a and gives no formula; the 180-degree-symmetric form was chosen because a membrane has two equivalent sides. The default a = 2.4 is past the point where the flank term goes negative, which is what distinguishes a membrane from a droplet: a purely repulsive M_i minimises its water contact and collapses into a compact blob, whereas an amphiphilic M_a wants membrane neighbours along its axis and water on both flanks, giving the one-particle-thick filaments the book describes.
- Equation (A.8) of the 2000 paper appears to divide the membrane mobility by m as well, which contradicts its own section 3.3 and the caption of figure 10 ('A and E particles have lower diffusion and repulsion rates than M and R'). The prose is followed: mobility_ratio slows only A.
- Equations (B.1) and (B.2) write the spontaneous resource-to-autocatalyst term as P_A' W_0 although (A.2) makes it proportional to the local resource; that term is not part of the book's scheme, so the discrepancy is recorded but not implemented.
- What this implementation does and does not reproduce. Reproduced and tested: membrane particles demix from water; isotropic M_i collapses into compact droplets while anisotropic M_a builds thin, orientationally ordered filaments (the contrast of figure 6.14); and the metabolism collapses when the recycling of waste into food is switched off. NOT reproduced: filaments closing into protocells, the osmotic permeability of a closed membrane to X and Y while it retains A, and growth followed by division. Prepared membrane rings dissolve rather than persisting, and no enclosed compartment ever retained autocatalyst. Those three phenomena are the results of [639]/[641], whose lattice rules were not obtainable, so they are listed in `phenomena` as published claims but are not asserted by the tests.
- The v1 parameters are replaced: grid becomes width/height, rate_constants becomes the four probabilities above, repulsion_strength becomes repulsion/neutral_coupling/temperature, anisotropy becomes the enum `membrane` plus the strength `anisotropy`, and dimensions is dropped because only the 2D model is implemented ([538]'s 3D lattice was not obtainable). The v1 claims of rate-constants and flow are dropped: no rate is published, and the lattice is closed rather than a flow reactor.

## Notes

No genetic material and no evolution mechanism - a self-maintenance model rather than a replicator model, which is why the book uses it to close the protocell discussion rather than the origin-of-information one.

---

*Specification: `catalog/chemistries/ono-ikegami-protocell.yaml` · generator: `chemart/chemistries/ono_ikegami_protocell.py` · tests: `tests/chemistries/test_ono_ikegami_protocell.py`*
