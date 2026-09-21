## Introduction

A self-propelled oil droplet is a real laboratory system, not a computer
model: a drop of oil, from a fraction of a millimetre to a few centimetres
across, that swims through water on its own. Martin Hanczyc, Takashi Ikegami and colleagues described the
best-known version in 2007. They loaded a drop of nitrobenzene (an oily
solvent that does not mix with water) with oleic anhydride, a fatty-acid
precursor, and placed it in alkaline water that already contained oleate soap
(oleate is the charged form of oleic acid, the fatty acid of olive oil).
Within seconds the droplet starts to move, and it keeps moving, turning,
stopping and starting again, until its fuel runs low.

The motion comes from chemistry at the droplet's surface. Where the precursor
meets water it splits into soap molecules and releases acid. The acid changes
how strongly the oil surface pulls on itself (its *interfacial tension*), and
because the reaction never runs evenly all round the droplet, one side ends up
pulling harder than the other. That imbalance drives a flow along the surface
and a circulation inside the drop, and the drop moves. It leaves a trail of
spent soap and acid behind it, and it can follow a gradient of acidity in the
water, which its authors read as a very simple form of chemotaxis, the way a
bacterium swims towards food.

Hanczyc and Ikegami (2010) built the system to ask how much of what living
things do (sensing, moving, responding) a handful of chemicals can do without
any biological machinery. The whole system has five components including
water. In their reading the oil-water surface is at once the sensor and the
motor, and they present the droplet as a chemical basis for *minimal
cognition*. Banzhaf and Yamamoto describe it in their chapter on wet
artificial chemistries, under "Droplet-Based Protocells" (book §19.2.5): a
*protocell* is a simple chemical structure studied as a model of an early or
minimal cell. What sets this one apart is that it has "no explicit information
subsystem": no genome, no replication, nothing copied. The information, the
book says, is in which chemicals are where, and it is "processed" by the
physics that makes the droplet move.

That makes it an unusual entry. The other protocells in the catalog are
simulations built around self-maintenance or reproduction: the
[chemoton](chemoton.md), [GARD](gard.md)'s growing lipid assemblies, and the
lattice cells of [Ono and Ikegami](ono-ikegami-protocell.md) and of
[Varela, Maturana and Uribe](autopoiesis-vmu.md). The other wet entries, such
as the [Oregonator](oregonator.md) or the [repressilator](repressilator.md),
are reaction networks whose interest lies in their kinetics. The droplet's
interest lies in its movement, and its chemistry is only two reactions long.
Chemart therefore returns those reactions together with the published
measurements, and does not simulate the swimming (see below for why).

## How it works

### The chemistry: making soap and acid at the surface

The droplet holds two chemicals: nitrobenzene, the solvent, which takes no
part in the reaction, and oleic anhydride, which is two oleic acid molecules
joined through a shared oxygen. The water around it is at pH 11, strongly
alkaline, and contains 10 mM (millimoles per litre) of oleate. The default
network is the reaction that happens where the two liquids touch:

```
oleic_anhydride + H2O -> 2 oleic_acid
oleic_acid -> oleate + H
```

Step by step:

1. A molecule of `oleic_anhydride` at the surface meets water and is split
   (*hydrolysed*) into two molecules of `oleic_acid`.
2. Each `oleic_acid` gives up a proton, `H`, and becomes `oleate`, the
   soap-like *surfactant*: a molecule with an oily tail and a charged head
   that sits at an oil-water surface and lowers its tension.

The sum of the two, one anhydride giving two oleates and two protons, is the
sentence Hanczyc and Ikegami (2010) write: the precursor "is hydrolyzed at the
oil-water interface to produce more oleate and protons". Chemart keeps the two
steps separate so that the proton release, which is what moves the droplet,
is visible on its own. Each anhydride carries two oleoyl groups (the fatty
chain), and each acid or oleate carries one, so the quantity
`2·[oleic_anhydride] + [oleic_acid] + [oleate]` never changes. Chemart records
it as a conservation law.

The reactions come without rates. No source that could be read publishes a
rate constant for this hydrolysis, so every rate is `None`, and the network
cannot be integrated in time.

### From acid to motion

The movement is not a reaction, so it is not in the reaction list. Hanczyc and
Ikegami (2010) describe it as a chain of causes, which Chemart stores in order
in `net.extras["interaction_law"]["steps"]`:

1. The hydrolysis does not run evenly over the surface, so the local pH near
   the droplet is not even either.
2. Interfacial tension *rises* as the pH falls from 11, up to a maximum at
   pH 9. Their pendant-drop measurements (a method that reads tension off the
   shape of a hanging drop) put a bare nitrobenzene drop at pH 11 with no
   oleate at 27 mN/m (millinewtons per metre). At pH 11 the tension with oleate present is already low,
   so the extra oleate changes little; it is the local drop in pH that
   matters. Dyes that change colour with pH show it falling as low as pH 7
   near the droplet.
3. A difference in tension along a surface pulls liquid along it, from low
   tension to high. This is the *Marangoni effect*, the same effect behind the
   "tears of wine" on the inside of a glass. Once a random fluctuation breaks
   the symmetry, the surface flow organises a pair of circulating flows
   (*convection*) inside the droplet.
4. The circulation brings fresh precursor to one end of the droplet and
   carries products away from the other, so the imbalance is kept up instead
   of evening out. The reaction feeds the flow and the flow feeds the reaction.
5. The droplet moves along the axis of its internal flow, leaving a trail of
   expelled surfactant and acid behind it.

Two things follow. The droplet's own acid trail is a pH gradient it moves
away from, and an external pH gradient can override it and steer the
droplet. And the reaction happens only at the surface, so the droplet's size
matters: a small droplet has more surface per unit of volume and burns
through its fuel faster, while a large one deforms and moves differently.

### Why Chemart stops at the chemistry

To simulate a trajectory you would need a force law: how fast a given tension
difference pushes a droplet of a given size. The sources that could be read
give the mechanism only in words, and the numerical study by the same group
(Matsuno, Hanczyc and Ikegami 2007) has no open copy. None of them gives a
droplet speed either. Chemart therefore does not integrate any motion and
states no speed; the reactor below (`continuous-space`) describes the dish the
droplet moves in, and `extras["space"]` holds the published dish sizes and the
droplet's geometry rather than positions.

## Using it

The default call returns the Hanczyc system with the composition given by
Hanczyc and Ikegami (2010), 0.5 M oleic anhydride in nitrobenzene added to
10 mM oleate at pH 11, as a 20 µL droplet (1 µL is 1 mm³), one of the sizes
Horibe et al. (2011) used. Nitrobenzene is the sixth species
listed; it is inert and appears in no reaction. The initial state is in mM,
with the oil and water species listed together:

```python
net.initial_state
# {'oleic_anhydride': 500.0, 'oleate': 10.0, 'H': 1e-08}
net.extras["phases"]
# {'oil': ['nitrobenzene', 'oleic_anhydride'], 'aqueous': ['H2O', 'oleate', 'H'],
#  'interface': ['oleic_acid']}
```

`H` at 1e-08 mM is pH 11. The rest of what Chemart knows about the droplet is
in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["droplet"]["diameter_mm"], a["droplet"]["surface_to_volume_per_mm"]   # (3.37..., 1.78...)
a["fuel_budget"]["precursor_umol"], a["fuel_budget"]["surfactant_umol"] # (10.0, 20.0)
a["modes_at_this_volume"]      # ['directional', 'circular', 'fluctuating']
a["above_tension_maximum"]     # True
```

The droplet is treated as a sphere. The fuel budget is plain arithmetic: 20 µL
at 0.5 M holds 10 µmol of anhydride, which can make at most 20 µmol of oleate
and as many protons. `modes_at_this_volume` gives the swimming patterns Horibe
et al. (2011) saw at this size (see Results), or `None` for a size they did
not report. `above_tension_maximum` says whether the bulk pH is above the
pH 9 tension maximum, the side on which the droplet's own acid raises the
tension; with `pH=8` it is `False`. `a["tension"]`, `a["measured_behaviour"]`
and `a["sizes_and_timings"]` hold the published measurements with their
sources.

**Droplet size.** `droplet_volume_uL` changes the geometry, the fuel budget and
the reported modes. `a["size_series"]` lists the six sizes Horibe et al.
used, whatever the volume you asked for:

```python
for row in net.extras["analysis"]["size_series"]:
    print(f'{row["volume_uL"]:>4g} uL  d = {row["diameter_mm"]:.2f} mm  '
          f'S/V = {row["surface_to_volume_per_mm"]:.2f} /mm  '
          f'modes = {row["modes"]}  attraction = {row["collective_attraction"]}')
```

```
   1 uL  d = 1.24 mm  S/V = 4.84 /mm  modes = ['circular', 'fluctuating']  attraction = None
   3 uL  d = 1.79 mm  S/V = 3.35 /mm  modes = None  attraction = True
  10 uL  d = 2.67 mm  S/V = 2.24 /mm  modes = None  attraction = None
  20 uL  d = 3.37 mm  S/V = 1.78 /mm  modes = ['directional', 'circular', 'fluctuating']  attraction = True
  30 uL  d = 3.86 mm  S/V = 1.56 /mm  modes = None  attraction = None
  50 uL  d = 4.57 mm  S/V = 1.31 /mm  modes = ['vibrating', 'circular', 'fluctuating']  attraction = False
```

`S/V` is the surface-to-volume ratio, 3/r for a sphere of radius r. `None`
means the paper did not report that size. Note that Horibe et al. mixed the
anhydride with nitrobenzene 1:1 by volume, not at 0.5 M, so the fuel budgets
in this table use Chemart's default loading, not theirs.

**The other droplet chemistries.** `system` switches to three related
experiments, which differ in where the energy comes from:

```python
for sys in ["fuel-surfactant", "maze-chemotaxis", "decanol-salt"]:
    net = chemart.generate_network("self-propelled-droplets", seed=1, system=sys)
    law = net.extras["interaction_law"]
    print(f"{sys}: {len(net.reactions)} reactions; direction = {law['chemotaxis']['direction']}")
    for r in net.reactions:
        print("   ", r.to_text())
```

```
fuel-surfactant: 1 reactions; direction = None
    precursor + H2O + catalyst -> octylaniline + byproduct + catalyst
maze-chemotaxis: 2 reactions; direction = down the gradient, toward the low-pH region
    hexyldecanoic_acid + OH -> hexyldecanoate + H2O
    hexyldecanoate + H -> hexyldecanoic_acid
decanol-salt: 0 reactions; direction = None
```

- `fuel-surfactant` is Toyota et al.'s (2009) droplet of 4-octylaniline
  carrying a catalyst that splits a precursor dissolved in the water. The fuel
  comes from outside and the catalyst is not used up, so the droplet is not
  consumed. The source does not name the precursor or the second product,
  hence `precursor` and `byproduct`.
- `maze-chemotaxis` is Lagzi et al.'s (2010) droplet of dichloromethane with
  2-hexyldecanoic acid, which sheds its surfactant into alkaline water and
  swims towards acid.
- `decanol-salt` is Čejková et al.'s (2014) decanol droplet in sodium
  decanoate, steered by a gradient of salt. Nothing reacts, so the reaction
  list is empty.

For these three, the initial state reuses `precursor_M` and `surfactant_mM`
for whichever species plays the matching role, and `extras["space"]["vessels"]`
is empty. They are reconstructed from abstracts; the parameter table says what
each parameter means for each system.

## Results

**Sustained motion from an onboard fuel.** Hanczyc et al. (2007) reported
that the loaded oil droplets "showed autonomous, sustained movement through
the aqueous media", and that internal convection "created a positive feedback
loop" by bringing fresh precursor to the surface. Hanczyc and Ikegami (2010)
add that the droplet starts moving within seconds. The same paper shows why
the fuel matters: without precursor, droplets still move in an imposed pH
gradient, but stop once the tension imbalance evens out, typically in a few
seconds. Hanczyc's review (2014) sorts self-moving droplets into three kinds,
with no reactive chemistry, with onboard fuel, or with an onboard catalyst,
and says onboard-fuel droplets move "for minutes to hours", until the fuel is
exhausted or waste products slow the reaction. The book adds that externally
supplied fuel "can keep the droplets moving for as long as desired", which is
Toyota et al.'s catalyst droplet. Chemart does not simulate the motion; its
tests check the chemistry that powers it (two oleates and two protons per
anhydride, the conserved oleoyl groups, no invented rate) and that the
catalyst droplet keeps its catalyst.

**Chemotaxis, and moving away from its own waste.** The droplet follows a pH
gradient, "moving towards the highest pH" in the book's words, and an
externally imposed gradient overrides the one it makes itself. The book's
figure 19.11 shows a moving droplet leaving a low-pH trail in a 27 mm dish,
and Hanczyc (2014) notes that the droplet moves "directionally away from the
waste that it produces". Hanczyc and Ikegami (2010) read this as
sensory-motor coupling with the sensor and the motor fused into one
structure, a mechanism they say "is not known in extant living systems". The
sign of the response depends on the chemistry: Lagzi et al.'s droplets move
towards *low* pH. Chemart records the direction for each system and its tests
check both signs, and check that no direction is recorded for the salt
droplets, whose accessible abstract does not state one. The tests also check
the published tension values (27 mN/m, maximum at pH 9, local pH down to 7)
and that the mechanism is stored as an ordered chain of causes.

**Four modes of swimming, set by size.** Horibe, Hanczyc and Ikegami (2011)
tracked single droplets of 1, 3, 10, 20, 30 and 50 µL over their lifetime.
They recorded speed and turning angle every second, averaged them over 20-second
bins, and sorted the result with a self-organising map, a kind of neural
network that groups similar data without being told the groups. It found four
modes: *directional* (fast, straight), *circular*, *fluctuating* (slow,
turning a lot) and *vibrating*. Droplets of 1 µL were circular and
fluctuating, 20 µL droplets added a directional mode, and only the 50 µL
droplets vibrated, their shape and internal flow having become unstable.
Middle-sized droplets were directional early and switched between circular
and fluctuating later. Chemart's tests check this size-to-mode table. The
same paper found that droplets slow down and turn more as they age, and that
turning angle is negatively correlated with speed: droplets change direction
while stopped. With no speed published, Chemart keeps these as recorded
statements, and a test checks that no speed is asserted.

**Two droplets attract.** Placed together in one dish and filmed for an hour, two 20 µL droplets
stayed closer to each other over the first 20 minutes than two droplets
filmed in separate dishes, and 3 µL droplets did the same. The attraction
faded over time, and 50 µL droplets showed none; the authors note that
convection was weak in older droplets and unstable in larger ones. The same
paper measured how droplets bounce off a glass wall, as a coefficient of
restitution (speed after divided by speed before): without a collision it
averaged about 1, with values above 1 when the reaction flared up again.
Chemart's tests check the attraction at 3 and 20 µL but not 50 µL.

**Shape takes over in large droplets.** Droplets of about 100 µm or less keep
a spherical shape. Hanczyc and Ikegami (2010) filmed droplets of 1, 5, 10 and
30 µL and saw larger ones, up to 0.5 cm across, change shape on a timescale
of seconds, with a horseshoe shape best supporting straight motion. They
conclude that above a few hundred microns the motion shifts from
convection-driven to shape-driven. Chemart records this statement in
`sizes_and_timings` but does not test it, since it models no shape.

**Maze solving.** Lagzi et al. (2010) put an acid source at one exit of a
maze; their droplets "find the shortest path through the maze" by following
the pH gradient. Čejková et al. (2014) navigated a maze with salt-driven
decanol droplets, which could also reverse direction repeatedly and carry a
chemical cargo. Chemart returns the chemistry and the recorded direction for
both, not the maze.

**Beyond this entry.** Hanczyc et al. (2007) saw some droplets turn into
"supramolecular aggregates resembling multilamellar vesicles" as surfactant
built up. The book also lists droplet division and fusion, droplets linked by
DNA strands anchored to their surface, and droplets that extract rare-earth
metal ions. These use other droplet chemistries, and Chemart does not
implement them.

## Further reading

- Matsuno, H., Hanczyc, M. M. & Ikegami, T. (2007). Self-maintained movements
  of droplets with convection flow. In *Progress in Artificial Life*
  Lecture Notes in Computer Science 4828, 179–188. Springer. The
  group's numerical model of the propulsion (book ref [556]).
