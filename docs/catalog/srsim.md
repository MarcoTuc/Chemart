# SRSim (rule-based spatial simulator)

`srsim` · *Gruenert & Dittrich, 2010-2013*

*Also known as:* *spatial rule-based modelling*, *BNGL + LAMMPS*

Rule-based chemistry that refuses to forget geometry. Molecules are spheres carrying binding sites at fixed positions, diffusing in three dimensions, and rules are written over patterns as in any rule-based language - but a bimolecular rule fires only when the two sites are actually close enough and correctly oriented. Geometry therefore decides which complexes can exist at all: change a monomer's binding angle and you get rods instead of rings. The species set is never enumerated, only explored.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | framework |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 18.3.3 |
| **refs** | [351], [388], [415], [861], doi:10.1186/1471-2105-11-307 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `space`, `initial-state`, `mass-conservation`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): elementary molecules (EM): spheres with a mass, a radius and a diffusion coefficient, carrying named components (binding sites) at fixed polar coordinates (dist, theta, phi) and optional internal states. A species is a complex: a connected site graph, taken up to isomorphism, whose id and structure are its canonical text, e.g. M(a[.],b[1]),M(a[1],b[.]). Because bonds are the straight connection of two component vectors, a complex also has a shape.

**R — reactions** (implicit, arity [1, 2]): BNGL-style rules over molecule patterns, e.g. 'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ k, k_rev. A rule forms one bond, breaks one bond or changes one internal state. A bimolecular rule fires only between molecules that are geometrically compatible: their distance is within distance_tolerance of the ideal bond length d_ij + d_kl and, for every bond either molecule already carries, the angle to the partner is within angle_tolerance of the ideal angle computed from the component coordinates. A compatible pair then reacts with probability 1 - exp(-k dt).

**A — reactor**: continuous-space
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("srsim", seed=1)
```

```
srsim: 13 species, 19 reactions, status=observed
provides: initial-state, mass-conservation, rate-constants, space, stoichiometry, topology
seed: 1
extras: analysis, bngl, conservation, init, interaction_law, kinetics, model, molecule_types, reaction_rules, rules, space, species_encoding
```

First reactions:

```
2 M(a[.],b[.]) -> M(a[.],b[1]),M(a[1],b[.])  [mass-action k=25.656340004316643 rules=grow]  (x20)
M(a[.],b[.]) + M(a[.],b[1]),M(a[1],b[.]) -> M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[.])  [mass-action k=25.656340004316643 rules=grow]  (x8)
M(a[.],b[.]) + M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[.]) -> M(a[1],b[2]),M(a[.],b[1]),M(a[2],b[3]),M(a[3],b[.])  [mass-action k=25.656340004316643 rules=grow]  (x4)
M(a[1],b[2]),M(a[.],b[1]),M(a[2],b[3]),M(a[3],b[.]) -> M(a[1],b[2]),M(a[3],b[1]),M(a[2],b[4]),M(a[4],b[3])  [mass-action k=25.656340004316643 rules=grow]  (x2)
M(a[.],b[1]),M(a[1],b[.]) + M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[.]) -> M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[3]),M(a[3],b[4]),M(a[4],b[.])  [mass-action k=25.656340004316643 rules=grow]  (x3)
M(a[.],b[1]),M(a[1],b[.]) -> 2 M(a[.],b[.])  [mass-action k=0.01 rules=grow_rev]  (x3)
M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[3]),M(a[3],b[4]),M(a[4],b[.]) + M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[.]) -> M(a[1],b[2]),M(a[3],b[1]),M(a[2],b[.]),M(a[4],b[3]),M(a[5],b[4]),M(a[6],b[5]),M(a[7],b[6]),M(a[.],b[7])  [mass-action k=25.656340004316643 rules=grow]  (x1)
2 M(a[.],b[1]),M(a[1],b[.]) -> M(a[1],b[2]),M(a[.],b[1]),M(a[2],b[3]),M(a[3],b[.])  [mass-action k=25.656340004316643 rules=grow]  (x2)
… and 11 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `model` | `enum` | `polymer` | structural | polymer: one monomer with two binding sites whose ideal angle is bond_angle (the paper's fig. 3 study of geometry-determined assembly); dimerization: A + B <-> AB, the reversible reaction the paper's supplement analyses; scaffold: the paper's scaffold-protein model (figs. 5-6) with the geometry of its spass.geo and the rates of its spass_w.bngl; custom: molecule_types, rules and init <br>one of `polymer`, `dimerization`, `scaffold`, `custom` |
| `n_molecules` | `int` | `50` | population | copies of the main molecule: M (polymer), A and B each (dimerization), A (scaffold) <br>`1` … `2000` · *range:* the paper's scaffold run uses 300 A and 30 S in a 100^3 box for 100000 steps |
| `n_scaffolds` | `int` | `10` | population | scaffold model only: copies of the scaffold protein S, which has four binding sites <br>`0` … `500` · *range:* the paper uses 30 |
| `bond_angle` | `float` | `90.0` | spatial | polymer model only: ideal angle between the monomer's two binding sites, which is what decides the shape of the assemblies <br>`0.0` … `180.0` · *range:* fig. 3: 180 gives rods, 90 closed quadratic structures, intermediate angles rings and helices |
| `k_on` | `float` | `0.5` | kinetic | polymer and dimerization: microscopic binding rate k2_mic, applied to a geometrically compatible pair (probability 1 - exp(-k dt) per step) <br>≥ `0.0` |
| `k_off` | `float` | `0.01` | kinetic | polymer and dimerization: microscopic unbinding rate of a bond <br>≥ `0.0` |
| `scaffold_binding` | `float` | `1.0` | kinetic | scaffold model only: multiplier of the published scaffold association rates d2_on and d3_on, i.e. the rate ks of figure 6 <br>`0.0` … `1000.0` · *range:* figure 6 compares 1 ('w-Scf', the published d2_on = 0.05) with 0 ('wo-Scf') |
| `rate_scale` | `float` | `1.0` | kinetic | multiplies every rule rate, like SRSim's preFactBindR/preFactBreakR/preFactModifyR prefactors; scaling the whole rule set rescales time, so rate_scale x simulated time reproduces a much longer published run <br>≥ `0.0` |
| `box` | `float` | `12.0` | spatial | edge of the cubic reaction volume; the box is periodic and must be wider than twice the largest reaction distance <br>≥ `1.0` · *range:* the published examples use 80^3, 100^3 and 200x30x30 |
| `steps` | `int` | `1500` | kinetic | molecular-dynamics steps; reactions are looked for after every one of them <br>`1` … `200000` · *range:* the published runs are 100000 to 3000000 steps |
| `dt` | `float` | `0.02` | kinetic | time step of the molecular dynamics <br>≥ `1e-09` · *range:* the paper's rule of thumb: sqrt(6 D dt) should be about a tenth of a particle diameter |
| `diffusion` | `float` | `1.0` | spatial | diffusion coefficient of an elementary molecule (the scaffold S gets a fifth of it, being five times heavier) <br>≥ `0.0` · *range:* the paper estimates 8e-11 m^2/s for a haemoglobin-sized protein (Stokes-Einstein) |
| `temperature` | `float` | `1.0` | thermodynamic | k_B T of the Langevin thermostat; the friction is gamma_0 = k_B T / D <br>≥ `1e-09` |
| `distance_tolerance` | `float` | `0.5` | spatial | how far the distance of two molecules may deviate from the ideal bond length for them to react <br>≥ `0.0` · *range:* every published example uses DistanceDeviation = 0.4 or 0.5 |
| `angle_tolerance` | `float` | `30.0` | spatial | how far the angles at an already bound molecule may deviate from their ideal values for a new bond to be accepted; 180 accepts any direction <br>`0.0` … `180.0` · *range:* published AngularDeviation: 180 (scaffold), 50 (microtubule and spheres), 30 (Sierpinski tiles) |
| `refractory_time` | `float` | `1.0` | kinetic | after a bond breaks, both molecules cannot bind for this long, so that diffusion can separate them (SRSim's answer to geminate recombination); it should satisfy refractory_time x k_off << 1 <br>≥ `0.0` · *range:* the published runs use 50 to 200 time steps |
| `k_bond` | `float` | `5.0` | thermodynamic | spring constant K_d of the harmonic bond potential E_d = K_d (d - d_ij - d_kl)^2 <br>≥ `0.0` · *range:* the published input scripts use fBond = 5.0 |
| `k_angle` | `float` | `5.0` | thermodynamic | spring constant K_a of the harmonic angle potential E_a = K_a (alpha - alpha_ijk)^2 <br>≥ `0.0` · *range:* published fAngle: 5.0, 15.0, and 0.0 for the scaffold run of figure 6 (angular forces switched off) |
| `k_repulsion` | `float` | `5.0` | thermodynamic | amplitude A of the soft-sphere repulsion E_s = A [1 + cos(pi r / r_c)], r_c = r_i + r_j <br>≥ `0.0` · *range:* the published input scripts use fRepulsion = 5.0 |
| `integrator` | `enum` | `langevin` | stochastic | langevin: the paper's equation F = F_S - gamma_0 v + sqrt(2 k_B T gamma_0) xi integrated with velocities (LAMMPS fix langevin + nve); brownian: its overdamped limit, whose mean square displacement is exactly 6 D t <br>one of `langevin`, `brownian` |
| `orientation` | `enum` | `bound-only` | stochastic | how an unbound molecule, which SRSim gives no rotational orientation, is treated: bound-only accepts it (SRSim's own test, which only checks the angles of bonds it already has); sampled accepts it with the probability (1 - cos angle_tolerance)/2 that a random orientation offers the component, which makes the reactive volume of the paper's supplement exact <br>one of `bound-only`, `sampled` |
| `well_mixed` | `bool` | `` | structural | run the non-spatial control instead: molecules have no positions and a pair is geometrically compatible with probability V_react / V_reactor, which is the well-stirred limit the paper compares against BioNetGen in figure 6 |
| `molecule_types` | `list` | `` | structural | custom only: elementary species, i.e. the JSON form of SRSim's .geo file (mass, radius, diffusion and the polar coordinates of every component) <br>*range:* e.g. [{'name': 'M', 'radius': 0.5, 'mass': 1.0, 'diffusion': 1.0, 'sites': [{'name': 'a', 'dist': 1.0, 'theta': 90, 'phi': 0}, {'name': 'b', 'dist': 1.0, 'theta': 90, 'phi': 90, 'states': ['u', 'p']}]}] |
| `rules` | `list` | `` | kinetic | custom only: BNGL-style rules with microscopic rates; a bare site is free, s!1 is the bond of this rule, s!+ bound to anything, s!? untested, s~state an internal state; one rule forms one bond, breaks one bond or changes one state, over one or two elementary molecules <br>*range:* e.g. ["'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ 0.5, 0.01", "'phos' A(x~u) -> A(x~p) @ 0.1"] |
| `init` | `dict` | `` | population | custom only: how many copies of each molecule type are placed in the box at random positions <br>*range:* e.g. {'M': 50} |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- geometry decides the reachable species: with a 180-degree monomer only rods form, with 60 degrees triangles close, and with 120 degrees no small ring can close at all (paper fig. 3)
- scaffold proteins co-localise their ligands: binding A to the scaffold roughly doubles the phosphorylated A in the spatial simulation while the well-mixed simulation of the same rules shows no increase (paper figs. 5-6)
- in the fast-diffusion limit the spatial simulation reproduces well-stirred mass-action kinetics, with the equilibrium constant k2_mic V_react / k_off of the paper's supplement
- a refractory time suppresses geminate recombination; the molecules waiting it out are t_ref k_off [C] per species, as the supplement predicts
- self-assembly of filaments, rings and larger complexes from a handful of rules, with a combinatorially large (unbounded) species set that is never enumerated

## Sources

- Banzhaf, W. & Yamamoto, L. (2015). Artificial Chemistries, section 18.3.3: one paragraph on SRSim as a spatial rule-based simulator that 'takes into account the detailed spatial conformation of each protein ... and positions each molecule in space such that only neighboring molecules may react with each other'.
- Gruenert, G., Ibrahim, B., Lenser, T., Lohel, M., Hinze, T. & Dittrich, P. (2010). Rule-based spatial modeling with diffusing, geometrically constrained molecules. BMC Bioinformatics 11:307 (book [351]), doi:10.1186/1471-2105-11-307, open access: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-307 - the geometry model (fig. 2), the reaction criterion and reactive volumes (fig. 4), the Langevin/bond/angle/soft-sphere potentials, the refractory time, and the four applications (scaffold proteins figs. 5-6, microtubules fig. 7, spheric self-assembly fig. 8, DNA Sierpinski tiles figs. 9-10).
- Additional file 2 of the same paper, 'Calculation of Kinetic Parameters from Macroscopic Values': k1_mic = k1_mac; k2_mic = 3 k2_mac / (2 pi eta (3 d0^2 t_dist + t_dist^3)(1 - cos t_ang)^2); the geminate-recombination analysis phi = 1 - K2/K1 and the refractory-time corrections [A*] = t_ref k_-1 [C] and [C''] = [C0]/(1 + t_ref k_-1). https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2F1471-2105-11-307/MediaObjects/12859_2009_3764_MOESM2_ESM.PDF
- Additional file 1 (SRSimSrc.zip): the released C++ sources. RuleSys/geometry_definition.cpp computes the ideal angle between two components from their polar coordinates (unit vectors, acos of the scalar product) and reads DistanceDeviation/AngularDeviation; LammpsModules/compute_reapot_atom.cpp tests candidate pairs (distance window (d0 +/- deviDist)^2 and the routine testSiteGeo, which compares the angle to every already bound component against its ideal value); LammpsModules/fix_srsim.cpp refracts both molecules after a break.
- Additional file 3 (ExamplesSrc.zip): the published input files. 001.scaffold/spass_w.bngl and spass.geo give the scaffold model reproduced here (S with four sites at theta = 0, phi = 0/72/144/216, dist 8, radius 8, mass 5; A with sites r (theta 0, dist 1, states nP/P) and s (theta 100, dist 1), radius 1; d1_on 1e-4, d1_off 0.1, d2_on 0.05, d2_off 1e-3, d3_on 1e-5, d3_off 0.9; DistanceDeviation 0.5, AngularDeviation 180) and spass_w.in its prefactors (preFactBindR 4.7, preFactBreakR 1e-3, preFactModifyR_1 1e-3, preFactModifyR_2 4.7), fRepulsion/fBond/fAngle = 5/5/0 and the 100000-step run. Its w_scf.srsim.gdat, wo_scf.srsim.gdat, w_scf.gdat and wo_scf.gdat are the data of figure 6: in SRSim 107 phosphorylated A with the scaffold against 81 without, in BioNetGen 61 with against 79 without.
- SRSim manual (docs/Manual.pdf of additional file 1): the .geo/.tgeo file formats and the meaning of the fix srsim arguments.

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives SRSim one paragraph, so everything is reconstructed from the 2010 paper, its supplement, its sources and its example models. LAMMPS and BioNetGen are not dependencies: the molecular dynamics, the BNGL-style rule language and the rule engine are reimplemented in numpy at the scale that fits Chemart's budget.
- Rule language: a documented BNGL subset. A rule has one or two elementary molecules per side (SRSim executes mono- and bimolecular reactions) and performs exactly one action - form a bond, break a bond, or change one internal state - which is SRSim's set of modifying, binding and breaking rules. Patterns constrain one elementary molecule each, not larger subgraphs, so context conditions that span three or more molecules (used in the paper's microtubule and Sierpinski models) cannot be written; exchange rules, molecule creation and deletion are not supported by SRSim either.
- Reaction criterion taken from the released code rather than from the prose: a candidate pair must lie in the distance window (d0 - t_dist, d0 + t_dist) with d0 = d_ij + d_kl, and for each partner every *already bound* component must keep its ideal angle to the new direction within t_ang (testSiteGeo). An unbound molecule has no rotational orientation and passes the angle test; the paper's supplement instead assumes both orientations are fixed when it derives the reactive volume, so orientation='sampled' reproduces that derivation exactly by accepting an unoriented partner with probability (1 - cos t_ang)/2. At t_ang = 180 degrees, the value the scaffold model uses, the two agree.
- Kinetics: a compatible pair reacts with probability 1 - exp(-k dt) (paper, Reaction Model and Kinetics). Monomolecular rules are applied per matching molecule with the same expression instead of SRSim's fragmented Gillespie loop over the reactant-pattern vector; within one time step the two are the same distribution. Pattern matches are computed once per step and every molecule reacts at most once per step, where SRSim recomputes the pattern indices after each reaction.
- Rate constants in the network: every observed reaction carries the macroscopic constant of the rule that fired it, k_mic * V_react for bimolecular rules (a volume per time, so that k [A][B] is the well-stirred rate) and k_mic for the others, following the conversion of additional file 2. Reactions produced by several rules sum those constants and extras.reaction_rules lists the rules. No symmetry factor is applied: as in KaSim's default convention, each ordered assignment of the rule's patterns to the pair is a separate instance.
- Space: a periodic box with the minimum-image convention (the scaffold example uses p p p; the other examples use reflecting walls, which are not offered). Bonded pairs are excluded from the soft-sphere repulsion, as LAMMPS does with special_bonds, and rigid bodies (used for the monomers of the spheric self-assembly) are not implemented.
- The well-mixed control is Chemart's, not SRSim's: BioNetGen is not a dependency, so well_mixed=true keeps the same rule engine but drops the positions and makes any pair compatible with probability V_react / V_reactor, which is exactly the well-stirred limit that additional file 2 derives. This is what the figure-6 comparison uses.
- Scaffold model: the published geometry and rates, with the rate prefactors of the published LAMMPS script applied (bind x 4.7, break x 1e-3, monomolecular modify x 1e-3, bimolecular modify x 4.7). The published run is 100000 steps; rate_scale compresses it (rate_scale x simulated time is the published simulated time), because scaling every rate rescales time. Molecule numbers are scaled down at constant concentration. The absolute counts of figure 6 are therefore not reproduced; the effect is: with the scaffold the number of phosphorylated A roughly doubles in the spatial simulation and does not increase in the well-mixed one.
- Species identity: the canonical text of the complex (colour refinement, then the least breadth-first text over the refined roots, sites grouped by name and bonds numbered by appearance), so isomorphic complexes are one species and components with the same name are interchangeable, as in BNGL. The refinement is exact for the complexes these models build; texts longer than 120 characters are replaced by formula#hash to keep species ids readable.
- v1 parameters agent_rules (callable) and geometry (matrix) became the model choice plus the custom molecule_types/rules/init, in JSON and in the rule mini-language. The reactor, the tolerances and the force constants are exposed with the values the published input files use.
- The network has status 'observed': it is the set of complex-level reactions that fired, with counts. Nothing is enumerated in advance - that is the point of rule-based modelling - so the network is a sample of an unbounded species set, and the rule set itself is kept in extras.rules and extras.bngl.

## Notes

The book's answer to combinatorial explosion: never enumerate. The rule set, the molecule geometries and the reaction criterion are kept unflattened in extras (rules, bngl, molecule_types, space.reaction_criterion) so that an exporter can hand them to BioNetGen or SRSim itself, while the network reports only what actually happened in the reactor.

---

*Specification: `catalog/chemistries/srsim.yaml` · generator: `chemart/chemistries/srsim.py` · tests: `tests/chemistries/test_srsim.py`*
