## Introduction

SRSim is a simulator for proteins that bind to each other and build complexes,
in which the shape of each molecule decides what it can bind to. It was written
by Gerd Gruenert, Peter Dittrich and colleagues in Jena, and published in 2010
as "Rule-based spatial modeling with diffusing, geometrically constrained
molecules". It combines two ideas that had rarely been combined.

The first idea is *rule-based modelling*. Many cell proteins carry several
binding sites and several sites that can be modified, for example
phosphorylated. Each combination of bound and modified sites is, strictly, a
different chemical species. The 2010 paper's example is the protein p53, whose
27 phosphorylation sites alone allow 2^27 = 134,217,728 states. Nobody can list
such a network reaction by reaction. A rule-based language such as BioNetGen
instead writes *rules* over partial descriptions of molecules: "a free site `a`
of any `M` can bind a free site `b` of any other `M`", whatever the rest of the
two molecules looks like. The species are then generated only as the
simulation needs them.

The second idea is *space and geometry*. In SRSim every elementary molecule is
a small sphere that diffuses in a three-dimensional box, and its binding sites
point in fixed directions from its centre. Two molecules can bind only if they
are at the right distance and their sites face each other at the right angle.
So the same rule can build very different things. A monomer whose two sites
point in opposite directions (180°) makes straight rods. One whose sites are at
90° makes mostly closed squares. Complexes that the geometry forbids are never made,
even though the rules allow them. The authors argue that this can change the
behaviour of a whole reaction system, and that simulations without geometry
cannot show it without special treatment.

SRSim is a *framework*: a general simulator that takes a model (a rule file
plus a file giving each molecule's geometry), not a single chemistry. The
original is an extension of the molecular dynamics package LAMMPS that reads
models written in BioNetGen's language. Banzhaf and Yamamoto mention it in one
paragraph, in their chapter on modelling biological systems (book §18.3.3).
They present it as the way out of the *combinatorial explosion* of large
reaction networks: rather than enumerating every compound and reaction in
advance, compute reactions on the fly from general rules, which makes the
chemistry constructive (open-ended in its species). Its nearest neighbour in
the catalog is the [Kappa calculus](kappa-calculus.md), a formal language of
the same rule-based kind but without space. Among the spatial chemistries,
[Squirm3](squirm3.md) and the [flow artificial chemistry](flow-ac.md) also
place molecules in space, but neither gives a molecule a shape with directed
binding sites.

## How it works

### Molecules: spheres with directed sites

An *elementary molecule* (EM) is a sphere with a radius, a mass and a
diffusion coefficient (how fast it spreads by random motion). It carries named
*components*, the binding sites. Each component sits at a distance `dist` from
the centre, in a direction given by two angles, `theta` (measured from the
vertical axis) and `phi` (measured around it), as on a globe. A component may
also have an internal state, such as `u` (unmodified) or `p` (phosphorylated).

A *complex* is a set of EMs joined by bonds, each bond linking one component to
another. The complex is the species: two complexes with the same molecules and
the same bonds, in whatever order they were built, are the same species. In
Chemart a species is named by a canonical text of this site graph. For
example, a dimer of the monomer `M` with sites `a` and `b` is

```
M(a[.],b[1]),M(a[1],b[.])
```

that is, two `M`s; in the first, `a` is free (`.`) and `b` carries bond 1; in
the second, `a` carries bond 1 and `b` is free. So bond 1 joins `b` of one
molecule to `a` of the other. An internal state is written in braces, as in
`c{u}`.

Because a bond is the straight join of two component vectors, its ideal length
is the sum of the two component distances, and the ideal angle between two
bonds on one EM is the angle between its two components. A complex therefore
has a shape.

### Rules

Rules are written in a subset of BioNetGen's language. The default model has a
single rule:

```
'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ 0.5, 0.01
```

It says: a free site `a` of one `M` and a free site `b` of another may join
into bond `!1`, at rate 0.5; the bond breaks again at rate 0.01. Sites not
mentioned in a pattern are not tested, so the rule applies to lone monomers
and to the ends of long chains alike. Each rule does exactly one thing: form a
bond, break a bond, or change one internal state. These are the paper's
*binding*, *breaking* and *modifying* rules.

### The reactor: molecular dynamics plus a reaction test

Time advances in small steps of length `dt`, as in a molecular dynamics
simulation. Each step has two parts.

First the EMs move. Each one feels random kicks from the surrounding solvent
and friction (Langevin dynamics, or its simpler overdamped form, Brownian
motion), plus three forces: a spring that holds each bond at its ideal length,
a spring that holds each pair of bonds on one EM at its ideal angle, and a
soft repulsion that keeps spheres from overlapping.

Then the rule system looks for reactions. For a binding rule, two EMs that
match its two patterns must also be *geometrically compatible*:

- their distance must be within a tolerance `t_dist` of the ideal bond length;
- for every bond either EM already has, the angle between that bond and the
  direction to the new partner must be within a tolerance `t_ang` of the ideal
  angle between the two components.

An EM with no bonds yet has no fixed orientation (SRSim does not track the
rotation of single spheres), so it passes the angle test. A compatible pair
then binds with probability `1 - exp(-k dt)`, where `k` is the rule's
*microscopic* rate. Breaking and modifying rules fire with the same
probability for each molecule they match.

When a bond breaks, the two molecules are still touching and would usually
rebind at once, a problem called *geminate recombination*. SRSim's answer is a
*refractory time*: after a break, both molecules cannot bind for a while, which
gives diffusion time to separate them.

A reaction between two parts of the same complex is treated like any other
bimolecular reaction. This is how rings close.

### Worked example: a square closes

The default run puts 50 monomers `M` in a periodic box of side 12, with the two
sites at 90° and an angle tolerance of 30°. Among its reactions are these (from
`net.reactions`, with the rate annotations left out; `(xN)` is how often the
reaction fired):

```
2 M(a[.],b[.]) -> M(a[.],b[1]),M(a[1],b[.])                                 (x20)
M(a[.],b[.]) + M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[.]) -> M(a[.],b[1]),M(a[1],b[2]),M(a[2],b[3]),M(a[3],b[.])   (x4)
M(a[1],b[2]),M(a[.],b[1]),M(a[2],b[3]),M(a[3],b[.]) -> M(a[1],b[2]),M(a[3],b[1]),M(a[2],b[4]),M(a[4],b[3])   (x2)
```

Step by step:

1. Two free monomers meet at about the ideal bond length (1 + 1 = 2, give or
   take 0.5). Neither has a bond, so neither has an orientation to check, and
   they bind with probability `1 - exp(-0.5 × 0.02)`, about 1% per step while
   they stay in range. This happened 20 times.
2. A monomer joins the end of a three-molecule chain. Now the end molecule
   already has one bond, so the new partner must lie at 90° ± 30° from it.
   The chain grows with a bend at every molecule.
3. A four-molecule chain closes on itself: the free `a` of one end binds the
   free `b` of the other. The result has no free sites left (no `[.]`), so it
   is a closed ring, a square. It can only close because the chain has bent by
   about 90° at each molecule and the two ends face each other.

Rerun this with the sites at 180° and the chain stays straight: its two ends
point away from each other and can never meet.

The run ends after 1,500 steps (30 time units) with 48 bonds formed and 7
broken, three closed rings (two squares and one pentagon) and a largest complex
of ten molecules.

### Rate constants in the network

Each reaction in the network carries a well-mixed (macroscopic) rate
constant, so it can be compared with ordinary mass-action kinetics. The
conversion comes from the paper's second supplement. The *reactive volume*
`V_react` is the region around a molecule in which a partner counts as
compatible: for the distance test alone, a spherical shell of thickness
`2 t_dist` around the ideal bond length `d0`, of volume
`(8π/3)(3 d0² t_dist + t_dist³)`. The macroscopic constant is
`k_mic × V_react`. In the default run `d0 = 2` and `t_dist = 0.5`, so
`V_react ≈ 51.3` and binding reactions carry `k ≈ 25.66`, the `k` printed
next to them in the list of reactions above. When both partners' orientations
are random, the supplement narrows the volume by `(1 - cos t_ang)/2` per
partner.

## Using it

The default call above is the geometry study of the paper's figure 3: one
monomer type, one polymerisation rule, with the angle between the two sites
set by `bond_angle` (90° by default). The network has status `observed`: it
lists only the complex-level reactions that actually fired, with their counts,
because the full species set is unbounded and is never enumerated. The rule
set itself is kept unflattened in `net.extras["rules"]` and, as a BioNetGen
text, in `net.extras["bngl"]`. The summary of the run is in
`net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["events_by_rule"]     # {'grow': 48, 'grow_rev': 7}
a["largest_complex"]    # 10
a["closed_complexes"]   # 3
a["size_histogram"]     # {'1': 2, '2': 1, '3': 2, '4': 3, '5': 1, '6': 1, '7': 1, '10': 1}
```

`mean_counts` and `mean_state_counts` average the species counts and the site
states over the second half of the run. `net.extras["space"]` holds the box,
the geometry and the final position of every EM; `net.extras["kinetics"]`
holds each rule's reactive volume and macroscopic constant. A species whose
text would exceed 120 characters is named by its formula and a hash, such as
`M10#df2d2d`.

`chemart.evolve("srsim", ...)` runs the same simulation and returns it as a
trajectory: one frame per step, timed in simulated time (step × `dt`), the unit
the rate constants are given in. A frame's `state` counts the complexes and its
`fired` lists the reactions of that step. In the default run the 50 monomers
assemble into fewer, larger complexes:

```python
traj = chemart.evolve("srsim", seed=1)
[(f.t, sum(f.state.values())) for f in traj.frames[::500]]
# [(0.0, 50.0), (10.0, 17.0), (20.0, 13.0), (30.0, 12.0)]
```

**Geometry decides the shapes.** Vary `bond_angle` with a tighter tolerance
(each run takes a couple of seconds):

```python
for angle in (60.0, 90.0, 120.0, 180.0):
    net = chemart.generate_network("srsim", seed=1, bond_angle=angle, angle_tolerance=15.0)
    a = net.extras["analysis"]
    rings = sorted(s.structure.count("M(") for s in net.species if "[.]" not in s.structure)
    print(angle, a["largest_complex"], a["closed_complexes"], rings)
```

```
60.0 4 6 [3]
90.0 6 1 [5]
120.0 7 0 []
180.0 2 0 []
```

At 60° triangles close (six of them). At 120° chains grow to seven molecules
but no ring closes. At 180° the run ends with 23 dimers and no ring. Straight
rods do grow, but slowly, because an end accepts a partner only inside a
narrow cone: with the default 30° tolerance the same run reaches rods of four.
In a short run a 90° chain can also close into a pentagon, since the angle
test allows 75° to 105° and a chain can pucker out of the plane.

**Scaffold proteins.** `model="scaffold"` loads the paper's scaffold model
(figures 5 and 6) with its published geometry and rates. `scaffold_binding`
switches binding to the scaffold on (1) or off (0), and `well_mixed=True` runs
the same rules with no space. These are the settings of Chemart's test, which
compresses the published run of 100,000 steps with `rate_scale` and uses a
third of the molecules at the same concentration (about 9 seconds for four
runs):

```python
def run(**kw):
    return chemart.generate_network("srsim", seed=1, model="scaffold",
        n_molecules=100, n_scaffolds=10, box=69.3, steps=1500, dt=0.5,
        diffusion=2.0, angle_tolerance=180.0, k_angle=0.0, rate_scale=100.0, **kw)

run(scaffold_binding=1.0).extras["analysis"]["mean_state_counts"]["A.r~P"]
```

The mean number of phosphorylated `A` (`A.r~P`) comes out as 10.2 with the
scaffold and 6.2 without it in the spatial run, and 5.0 against 7.0 in the
well-mixed run.

**Checking the kinetics.** `model="dimerization"` is the reversible reaction
`A + B <-> AB` that the paper's supplement analyses. With
`n_molecules=50, box=12.0, steps=2000, dt=0.02, k_on=0.5, k_off=0.05,
refractory_time=2.0, integrator="brownian"` (about 13 seconds), seed 1 settles
at 14.9 free `A`, 14.9 free `B` and 35.1 dimers on average. `model="custom"`
takes your own `molecule_types`, `rules` and `init`; the parameter table below
shows their format.

## Results

Besides the geometry study of its figure 3, the 2010 paper presents four
applications, and says they "were engineered to test and demonstrate our approach,
rather than to deliver a highly detailed representation of a special
biological system". Units are arbitrary.

**Geometry decides which complexes exist.** Figure 3 shows one monomer with
two sites and one rule: with the sites in line it builds rods, with the sites
at 90° mostly closed squares, and with sites tilted out of the plane it can
build helices. The paper notes that telling a square from a helix properly
needs *dihedral* angles (the twist around a bond), which SRSim does not
implement; Chemart does not either. The Discussion draws a consequence for
network analysis: a set of species that is closed and self-maintaining (an
*organization*, in chemical organization theory) under the rules alone may not
be one once geometry is taken into account, and the reverse. Chemart's tests
check the geometric effect with three angles at a tolerance of 15°: triangles
close at 60°; at 120° chains grow but no ring smaller than five ever closes (a ring must
turn through at least 360°, and each corner turns by at most 75°); and at 180°
no ring closes at all.

**Scaffold proteins co-localise their ligands.** Scaffold proteins bind
several other proteins and are thought to speed up reactions among them by
holding them close together. In the paper's model, proteins `A` phosphorylate
each other when they meet and lose the phosphate over time. A large scaffold
`S` binds up to four `A`, holding unphosphorylated `A` tightly and releasing
phosphorylated `A` quickly. Bound `A` can then only move over the scaffold's
surface, so they meet, and phosphorylate each other, more often. With binding
switched on, the spatial simulation ends with 107 phosphorylated `A` against
81 with binding off (figure 6, from the data files published with the paper).
The same rules run in BioNetGen, which has no space, give 61 against 79: no
increase. The paper stresses that the faster phosphorylation on the scaffold
is not written into any rule; it emerges from the geometry. It also reports
that switching on angular forces, which push the bound `A` to one pole of the
scaffold, amplifies the effect. Chemart reproduces the effect, not the counts:
its test requires the spatial run with the scaffold to exceed the run without
by at least 30%, and the well-mixed run not to.

**Spatial kinetics agree with mass action when mixing is fast.** The rate
conversion described above assumes a well-stirred reactor, so a spatial run
with fast diffusion should reach the mass-action equilibrium. The paper's
second supplement also analyses the error
the refractory time introduces. At equilibrium the molecules waiting out their
refractory time number `[A*] = [B*] = t_ref k_-1 [C]`, where `t_ref` is the
refractory time, `k_-1` the unbinding rate and `[C]` the dimer concentration;
the supplement concludes that `t_ref k_-1` should be much smaller than 1.
Chemart's tests confirm the pieces: the reactive volume formula and its
`(1 - cos t_ang)²` factor; the equilibrium constant `[C]/([A][B]) = k_mic V_react / k_off`
in a well-mixed run (within 20% on average) and in a spatial Brownian run
(within 30%; the seed-1 run above gives 455 against a predicted 513); the number of refractory molecules
(within 25%); and the mean square displacement of free particles, `6 D t`,
behind the paper's rule of thumb for choosing the time step.

**Self-assembly and molecular machines.** Three further applications build
large structures from a handful of rules, with a species set too large to
enumerate:

- *Microtubules and motor proteins* (figure 7). 600 tubulin dimers assemble
  on a ring-shaped nucleus into a hollow tube of 13 protofilaments with a
  seam, as seen in real microtubules. Motor proteins added later walk along
  it by binding and breaking bonds and carry cargo to one end. The model has
  27 rules. The abstract describes this as transport "happening faster than
  diffusion".
- *Spheres* (figure 8). Monomers made of six EMs, held rigid, first form rings
  of five and six, then close into irregular spherical shells that mix cycles
  of five, six and seven monomers.
- *DNA Sierpinski triangles* (figures 9 and 10). A simplified version of
  Rothemund, Papadakis and Winfree's DNA tiles, which compute the XOR function
  layer by layer as they attach, reproducing the Sierpinski pattern, with an
  occasional assembly error.

Chemart does not reproduce these three. Its rule patterns constrain one EM
each, so the context conditions spanning three or more molecules used by the
microtubule and Sierpinski models cannot be written, and it has no rigid
bodies, which the spheres need (see the implementation decisions above).
Chemart's tests do check the species bookkeeping these models rely on:
complexes that differ only in the numbering of their molecules get the same
name, and every elementary molecule is conserved.

**Later work.** The book notes that SRSim was then applied to the mitotic
spindle checkpoint, "a complex network of molecular regulatory interactions
responsible for correct cell division", citing two 2013 papers: Ibrahim et
al., which describes how to build a spatial rule-based model from
experimental data and applies it to the human mitotic kinetochore, and
Tschernyschkow et al., on the inner kinetochore structure. Chemart does not
implement these models.

## Further reading

- Ibrahim, B., Henze, R., Gruenert, G., Egbert, M., Huwald, J. & Dittrich, P.
  (2013). Spatial rule-based modeling: a method and its application to the
  human mitotic kinetochore. *Cells* 2(3), 506–544. (Book [415].)
- Tschernyschkow, S., Herda, S., Gruenert, G., Döring, V., Görlich, D.,
  Hofmeister, A., Hoischen, C., Dittrich, P., Diekmann, S. & Ibrahim, B.
  (2013). Rule-based modeling and simulations of the inner kinetochore
  structure. *Progress in Biophysics and Molecular Biology* 113(1), 33–45.
  (Book [861].)
- Hlavacek, W. S., Faeder, J. R., Blinov, M. L., Posner, R. G., Hucka, M. &
  Fontana, W. (2006). Rules for modeling signal-transduction systems.
  *Science's STKE* 2006(344), re6. The review of rule-based modelling the book
  recommends. (Book [388].)
- Rothemund, P. W. K., Papadakis, N. & Winfree, E. (2004). Algorithmic
  self-assembly of DNA Sierpinski triangles. *PLoS Biology* 2(12), e424.
