## Introduction

This is a kinetic model of a primordial soup of polymers, built by Richard
Bagley and J. Doyne Farmer at Los Alamos and the Santa Fe Institute and
published in *Artificial Life II* (1992). It asks whether a network of
catalysed reactions among random polymers can, on its own, gather most of
the material in a flow reactor into a few specific molecules, the way a living
cell holds high concentrations of a few particular proteins while a lifeless
mixture spreads its material thinly over countless kinds. Their answer is yes,
under conditions they map out. When it happens they call the network an
**autocatalytic metabolism**, because like a metabolism it takes in food from
its surroundings and turns it into its own components.

The model starts where the *autocatalytic set* argument of Kauffman (1986) and
Farmer, Kauffman and Packard (1986) left off. There, molecules are strings of
two letters, `a` and `b`, that join (condensation) and split (cleavage), and
each string may catalyse some of these reactions. An autocatalytic set is a
group of strings in which every member is made by at least one reaction
catalysed by another member. Kauffman showed with random graphs that such sets
become almost certain when the chemistry is rich enough. But that is a
statement about the *wiring*. Bagley and Farmer added rate equations and a
flow of food, and asked what the wiring does to *concentrations*. They define
an autocatalytic metabolism as an autocatalytic set whose concentrations
differ significantly from what they would be if no reaction were catalysed.

The central observation is simple. A catalyst speeds up a reaction and its
reverse by the same factor, so it changes how fast a mixture reaches chemical
equilibrium but not where it ends up. In a closed vessel, catalysed or not,
the soup settles into the same featureless state: short strings are common,
long strings are exponentially rare, and all strings of a given length are
equally common. Catalysis can only shape the soup when something keeps it
away from equilibrium, here a steady inflow of small food molecules and a
steady outflow of everything. Too little flow and the soup relaxes to
equilibrium; too much and molecules are washed out before they react. In
between, the catalysed network can hold concentrations orders of magnitude
above the uncatalysed background.

It is a simulation model: ordinary differential equations (ODEs) for the
concentration of every polymer, plus a procedure, *metadynamics*, for keeping
the set of polymers finite. Banzhaf and Yamamoto describe it in the book's
section on autocatalytic sets in origin-of-life research (§6.3.1) and again in
the chapter on evolution (§7.3.1), as a chemistry that does not assume
evolutionary dynamics but may give rise to them. A companion paper by Bagley,
Farmer and Fontana (1992) added random fluctuations and argued that such
metabolisms can evolve.

Its nearest neighbours in the catalog share its polymers but not its
dynamics. [Kauffman autocatalytic sets](kauffman-autocatalytic-sets.md) is the
same string chemistry with the same random catalysis, but only as a graph,
with no rates. [RAF sets](raf.md) is an algorithm that finds autocatalytic sets
in such a graph, also without dynamics. The [Jain-Krishna
model](jain-krishna.md) keeps only the "who catalyses whom" graph and replaces
the weakest species over time, and [random catalytic
networks](random-catalytic-networks.md) drop molecular structure altogether.
What `bagley-farmer` adds is mass-action kinetics that respect chemical
equilibrium, a flow reactor, and catalysts that are tied up while they work.

## How it works

### Molecules are strings

The monomers are the letters `a` and `b` (Chemart allows up to four letters).
A polymer is a string of them, and strings have a direction, so `ab` and `ba`
are different molecules. The letters could stand for amino acids, making the
strings peptides, or for nucleotides; the paper says this changes the
parameters but not the form of the model. Water is a species too, written
`H`, and is held at a fixed concentration.

The paper's chemistry has no longest string. Chemart builds every string up to
`max_length` letters: with the default of 5 and two letters that is
2 + 4 + 8 + 16 + 32 = 62 polymers.

### Spontaneous reactions: join and split

Any two strings can join into their concatenation, releasing water, and any
string can be split by water at any point:

```
A + B  <->  C + H        (C is the string A followed by B)
```

Joining runs at rate `kf·A·B` and splitting at rate `kr·H·C`, where each
letter stands for a concentration (*mass action*: the rate is a constant times
the product of the reactant concentrations). The two constants are the same
for every pair of strings. Their ratio `κ = kf/(kr·H)` is the *equilibrium
constant*: at equilibrium every string satisfies `C = κ·A·B`, so a string of
length `n` has concentration `κ^(n−1)` times the product of its monomers'
concentrations. The number of strings of length `n` grows as `2^n`, so each
individual long string is very rare. The paper's example: with 20 kinds of
amino acid, one particular peptide of length 30 is roughly `20^−30` as common
as a monomer.

### Catalysis multiplies both directions

A few of the possible (reaction, catalyst) combinations are *strongly*
catalysed. Following Kauffman, the paper picks them at random: each polymer
catalyses each join/split pair with probability `p`. All other combinations
count as spontaneous. A catalyst `E` at concentration `E` multiplies both the
joining and the splitting rate by `(1 + ν·E)`, where `ν` is the *catalytic
efficiency*:

```
dC/dt = (1 + ν·E)·(kf·A·B − kr·H·C)
```

Because both directions are multiplied by the same factor, the bracket still
vanishes at `C = κ·A·B`: catalysis does not move the equilibrium.

### Saturation: catalysts are busy while they work

The formula above ignores that a catalyst and its substrates stay bound
together for a while, during which none of them can take part in anything
else. When most catalyst is bound at any moment the reaction is *saturated*.
The paper's simulations account for this with a simplified scheme (its
appendix). A catalysed join or split produces a bound complex instead of free
molecules, and every complex falls apart at the same rate `ku`. Rather than
track each complex, the paper keeps one pool per species holding everything
that species has bound up; Chemart names these pools `aa_bound`, `bbbbb_bound`
and so on.

Here are two reactions of the default network, catalysed by the string
`bbbbb`, and two of the unbinding reactions:

```
2 a + bbbbb -> aa_bound + bbbbb_bound + H     [k = ν·kf = 582153000]
aa + H + bbbbb -> 2 a_bound + bbbbb_bound     [k = ν·kr = 2242500]
aa_bound -> aa                                [k = ku = 50000]
bbbbb_bound -> bbbbb                          [k = ku = 50000]
```

In the first line two `a` molecules meet the catalyst `bbbbb` and join into
`aa`, releasing water; the new `aa` and the catalyst stay bound together. The
second line is the reverse, also through a complex. The unbinding lines return
the bound molecules to the free soup at rate `ku`. The catalytic rate
constants are large because the default `ν` is 897,000; with saturation the
catalysis is limited by how fast the complexes break up, not by `ν` alone.
With `saturation=False` Chemart uses the plain `(1 + ν·E)` form instead.

### The reactor: food in, everything out

The soup sits in a *chemostat*, a well-stirred vessel with a constant supply
and drain. Each species of the *food set* (by default `a` and `b`) flows in at
a constant flux `δ`, and every species, food or not, flows out at a rate
proportional to its concentration, with rate constant `K`. Reactions only
rearrange monomers, so the total monomer mass `m` obeys `dm/dt = (monomers fed)·δ − K·m`
and settles at a fixed value `m0` whatever the starting state. Chemart takes
`δ` and `m0` as parameters and sets `K` from them; with the defaults
`K = 17.9 × 2 / 2.0 = 17.9`. The run starts with all the mass in the food.

`δ` measures how hard the soup is driven. The paper also translates it into a
*mean reaction number* `r`, the average number of reactions a monomer takes
part in before it is washed out: infinite at equilibrium, falling as `δ`
grows.

### Why flow lets catalysis focus: the smallest example

The paper's Figure 1 shows the idea with two products. The monomers `a` and
`b` are fed in, and join into either `ab` or `ba`. Only the route to `ab` is
catalysed, by an enzyme `E` held at a fixed concentration:

```
ba + H  <->  a + b  <->(E)  ab + H
```

At steady state the ratio of the two products is (the paper's eq. 15)

```
[ab]/[ba] = (1 + β) / (1 + β/γ),    γ = 1 + ν·E,   β = K/(kr·H)
```

`γ` is the strength of catalysis and `β` measures how far the flow drives the
vessel from equilibrium. At equilibrium `β = 0` and the ratio is 1: catalysis
makes no difference. As `β` grows the ratio rises toward `γ`. Chemart builds
this network with an explicit catalytic link. With `kf = kr = H = E = 1` and
`ν = 10`, so `γ = 11`, integrating it to steady state gives exactly the
formula's values:

```
delta=0.001  K=0.001   ab/ba=1.0009 eq15=1.0009
delta=0.5    K=0.5     ab/ba=1.4348 eq15=1.4348
delta=5.0    K=5       ab/ba=4.1250 eq15=4.1250
delta=50.0   K=50      ab/ba=9.1967 eq15=9.1967
```

The example also shows that focusing needs *specific* catalysis: if `ab` and
`ba` were catalysed equally, their ratio would stay at 1 at any flow.

In a real network the enzyme is not held fixed by the experimenter. The
network has to make its own catalysts, which is where autocatalytic sets come
in: a set of strings that catalyse each other's formation keeps its catalysts
at high concentration by itself. The paper stresses that an autocatalytic set
in the graph is necessary but not sufficient; whether it becomes a metabolism
depends on the rate constants.

### Metadynamics: keeping the soup finite

With unbounded strings there are infinitely many species, and in a continuous
ODE every one of them gets a nonzero concentration as soon as `t > 0`. A real
vessel holds a finite number of molecules, so the paper sets a concentration
*threshold* equal to one molecule in the vessel. Only species above it may
react. The paper solves for the steady state of the current set of equations,
adds the species that have risen above the threshold and removes those that
have fallen below, solves again, and repeats until nothing crosses the
threshold. The result is the *metadynamical fixed point*. Chemart offers this
with `threshold > 0`; by default it keeps the whole network up to `max_length`
instead, which is the exact system the paper's thresholded and lumped version
approximates.

## Using it

The default call above builds the full network up to length 5 with the rate
constants of the paper's Table 2 (a), its best parameters for its 15-reaction
test network: `kf = 649`, `kr = 2.5`, `ν = 897,000`, `ku = 50,000`,
`δ = 17.9`, `m0 = 2`. The catalytic links, however, are not the paper's: they
are drawn at random with `p = 0.01`, which gives 121 links for seed 1 (the
emergence paper gives no value of `p`). The 125 species are the 62 polymers, water `H`
and 62 bound pools; the 696 reactions are 392 spontaneous joins and splits,
242 catalysed ones and 62 unbindings.

What you need to simulate the chemostat is on the network:

```python
net.inflow                        # {'a': 17.9, 'b': 17.9}   flux δ per food species
net.outflow                       # 17.9                     washout rate K
net.initial_state                 # {'a': 1.0, 'b': 1.0, 'H': 1.0}
net.extras["buffered"]            # ['H']                    hold water constant
net.extras["catalytic_links"][0]  # {'reaction': 'a + a <-> aa', 'catalyst': 'bbbbb', 'nu': 897000.0}
net.extras["conservation"]        # letter counts per species (monomer a, monomer b)
```

`chemart.simulate.ode` integrates the chemostat from these fields: the
inflow, the washout, and water held fixed. The rate equations are stiff (the
rate constants span more than eight orders of magnitude), so use an implicit
solver. `solver="BDF"` is the fastest at the strong flows below: about half a
second a run, against three or four for the default `"LSODA"`. Near
equilibrium (`delta` below about 1) it is the other way round: LSODA takes
under a second, and BDF does not finish within two minutes.
This helper runs the chemostat for 200 washout times and adds each species'
free and bound amounts:

```python
import chemart
from chemart import simulate

def steady_state(net):
    """Run the chemostat for 200 washout times; return free + bound totals."""
    x = simulate.ode(net, 200 / net.outflow, points=2, solver="BDF").frames[-1].state
    total = {}
    for s, v in x.items():
        if s != "H":
            total[s.removesuffix("_bound")] = total.get(s.removesuffix("_bound"), 0.0) + v
    return total
```

**Catalytic focusing.** Compare the same network with and without catalysis
(`p=0.0`) at a strong flow, `delta=1000`:

```python
cat = steady_state(chemart.generate_network("bagley-farmer", seed=1, delta=1000))
unc = steady_state(chemart.generate_network("bagley-farmer", seed=1, delta=1000, p=0.0))
for s in sorted(cat, key=lambda s: cat[s] / unc[s], reverse=True)[:4]:
    print(f"{s:6} {cat[s]:.3g} with catalysis, {unc[s]:.3g} without ({cat[s] / unc[s]:.1f}x)")
```

```
babaa  0.0506 with catalysis, 0.00231 without (21.9x)
bbaba  0.0404 with catalysis, 0.00231 without (17.4x)
baaab  0.0282 with catalysis, 0.00231 without (12.2x)
bbbbb  0.021 with catalysis, 0.00231 without (9.1x)
```

Without catalysis all strings of length 5 have the same concentration; with
it, three of them hold 30% of the total mass. The whole comparison takes about
a second.

Repeating the comparison over `δ` (seed 1), with the default solver for the
slowest flow, shows the rise and fall the paper describes, though on this random network the effect is weaker than in the
paper's hand-picked one:

| `delta` | food's share of the mass | largest boost | species boosted over 10× | their share of the mass |
|---|---|---|---|---|
| 0.01 | 0.03 | 6.1× (`a`) | 0 | 0 |
| 17.9 (default) | 0.03 | 4.8× (`baba`) | 0 | 0 |
| 100 | 0.03 | 4.3× (`bbab`) | 0 | 0 |
| 1,000 | 0.11 | 22× (`babaa`) | 3 | 0.30 |
| 10,000 | 0.79 | 509× (`babaa`) | 5 | 0.13 |
| 100,000 | 0.97 | 22,000× (`babaa`) | 16 | 0.01 |

At the highest flow most of the mass is food that leaves before reacting, so
the large boosts concern tiny amounts.

**Saturation holds mass in complexes.** The near-equilibrium row (`delta=0.01`)
is not flat, even though catalysis cannot move an equilibrium. The reason is
saturation. With `saturation=False`, catalysed and uncatalysed free
concentrations agree to within about 1% at `delta=0.01`. With the default
saturation, the free species still satisfy `C = κ·A·B`, but 1.11 of the 2.0
units of monomer mass sit in the bound pools, so every free concentration is
0.40 to 0.83 of its uncatalysed value. Read the bound pools as part of the
soup when you compare runs.

**The paper's Figure 1.** Give the link explicitly and hold the enzyme fixed,
as in the table in *How it works*:

```python
net = chemart.generate_network("bagley-farmer", max_length=2, links=["a + b <-> ab | bb"],
                               saturation=False, kf=1.0, kr=1.0, nu=10.0, delta=5.0, m0=2.0)
net.initial_state["bb"] = 1.0
net.extras["buffered"] = ["H", "bb"]     # the enzyme is held fixed, as in the paper
```

`links` takes any list of `'A + B <-> C | E'` strings and replaces the random
draw, which is how to encode a network of your own.

**Metadynamics.** A threshold returns the graph at its metadynamical fixed
point, with status `truncated`. With the paper's Figure 3 constants except
`nu=100` and `ku=1000` (the setting used in the tests), on strings up to
length 4 and the Figure 3 food set, it takes about 2 seconds:

```python
net = chemart.generate_network("bagley-farmer", seed=2, max_length=4, p=0.08,
                               kf=100.0, kr=10.0, nu=100.0, ku=1000.0, delta=100.0, m0=3.0,
                               food_set=["a", "b", "ab", "ba"], threshold=0.01)
a = net.extras["analysis"]
a["outcome"], a["metadynamics_rounds"], len(a["active_species"])   # ('fixed', 4, 20)
```

`a["active_species"]` lists the 20 species above the threshold and
`a["metadynamical_fixed_point"]` their concentrations. `outcome` is not always
`'fixed'`. With the exact Figure 3 constants (`nu=1e4`, `ku=1e4`) the same call
ends with `'cycle-merged'` after 10 rounds, and so does the default network
with `threshold=0.01` for seeds 1, 2 and 3: removing a species below the
threshold also removes the reactions that consume it, so it rises again, and
Chemart merges the repeating sets (see the implementation decisions).

## Results

**At equilibrium the soup is featureless.** The paper first derives, from the
classical theory of polycondensation, that near equilibrium concentration
falls off exponentially with length and is the same for every string of a
given length. Its Figure 5(a) confirms it in simulation at `δ = 0.01`, where
the mean reaction number is about 210,000. Chemart's tests check the exact
statement behind this: the equilibrium state `C = κ·A·B` is a steady state of
the closed vessel whatever the catalytic links (without saturation), and a
saturated closed vessel relaxes to it.

**Driving focuses the material, overdriving loses it.** At `δ = 10^5`
(mean reaction number about 33) the catalysed network stands orders of
magnitude above the background and holds most of the mass (Figure 5(b)). At
`δ = 10^7.5` (about 0.5 reactions per monomer) molecules leave before they
react and the structure fades (Figure 5(c)). The paper measures the effect by
the slope `Λ` of log concentration against length, a least-squares fit: the
flatter the profile (`Λ` near zero), the further from equilibrium. For its
15-species test network (Figure 2) at `δ = 17.9`, the network's profile is
nearly flat while the background and the equilibrium profile fall steeply
(Figure 6). As `δ` rises, `Λ` climbs from its equilibrium value, peaks close to
zero at roughly `δ = 10^2`, and drops again (Figure 7); the share of mass held
by the network peaks for `δ > 10^3` (Figure 8). The authors conclude that
interesting behaviour needs a flow of energy that is neither too small nor too
large. Chemart tests the mechanism exactly: eq. 15 on the Figure 1 network at
four flows, and that the ratio is 1 near equilibrium and large when driven.
Figures 5 to 8 themselves are not reproduced: the Figure 2 network exists only
as a drawing, and the paper lumps the background into per-length averages
where Chemart keeps every string. The table under *Using it* shows the same
rise and fall on a random network.

**A broad favourable regime.** Varying one parameter at a time around the
best values, and calling the network a metabolism when its mass exceeds that
of the food and background together, the paper finds wide ranges (Table 2).
For the Figure 2 network, which has one catalytic link per reaction, the
metabolism survives for `ν` from about `10^1.6` to beyond `10^10` and `δ` from
`10^1` to beyond `10^6`. A variant with the same reactions but 118 catalytic
links supports metabolisms over a wider regime. Figures 9 to 12 show `Λ` over
pairs of parameters. The authors note that real catalytic efficiencies and
joining rates are probably at the low end of the favourable regime, and that
polymerisation must be favoured even without catalysis, for example by energy
from pyrophosphate driven by light, a mechanism they simulated and then
replaced by an effective joining rate. Chemart's defaults are Table 2 (a); the
table itself is not reproduced, and some of its printed ranges exclude their
own optima (see the implementation decisions).

**Topology matters, and more links make it matter less.** In Figure 13 the
reactions and rates are fixed and the catalytic links are redrawn at random,
5,000 networks with up to 200 links. With few links, `Λ` stays near its
equilibrium value; at intermediate numbers it depends strongly on which links
were drawn; with many links nearly all networks depart from equilibrium alike.
With one link per reaction, changing the efficiency or the joining rate of a
single reaction can shift `Λ` a lot (Figures 14 and 15); with 118 links it
hardly matters. Chemart has the ingredients (`p`, `links`,
`nu_distribution="uniform"`) but does not reproduce these sweeps.

**Robustness to diet.** To test whether the network behaves like a
metabolism, the paper changes the food of a 22-polymer network (the endpoint of
a metadynamics run) while keeping the total inflow of mass constant (Figure 16,
Table 3). With the default food `{a, b, ab, bb}` the slope is −0.125. Two diets
leave the metabolism essentially unchanged (−0.144 and −0.160); two make it
"die", falling below the uncatalysed level (−0.384 and −0.411, against −0.359
with no catalysis). Reading the reaction graph did not predict which diets
survive. Chemart does not reproduce this experiment.

**Metadynamics reaches a unique fixed point.** The paper reports that its
deterministic metadynamics always approached a unique fixed point, whether or
not saturation or pyrophosphate were modelled, which is what made its fast
algebraic solver possible. Chemart's test checks one case that reaches a fixed
graph. On other settings, including the exact Figure 3 constants, Chemart's
graph updates cycle instead (see *Using it*); the implementation decisions
attribute this to keeping every spontaneous reaction explicit rather than
lumped.

**Evolution by jumps between fixed points.** The companion paper (Bagley,
Farmer and Fontana, 1992) adds the rare spontaneous reactions that a
deterministic model misses. A species outside the metabolism can appear by a
chance reaction and, if it catalyses its own production directly or through a
short loop, grow by orders of magnitude and shift the metabolism to a new
fixed point, where it stays until the next such event. Such events mostly
start in the *shadow*, the background strings that the metabolism's own
members can make in one step. For a single self-catalysing newcomer the paper
derives the chance of taking hold as `1 − d/c`, where `c` is its catalysed
production rate and `d` its loss rate. In preliminary "stochastic
metadynamics" runs, a metabolism grew from its food set through a series of
such jumps until it reached the preset limit of 50 species, with new species appearing and old ones going extinct
(its Figure 3). The authors call these "punctuated equilibria" and argue that
variation (fluctuations) and selection (kinetics) justify the word evolution.
Chemart does not implement this: the generator returns one network, not a
history of fixed points, and the earlier `mutation_rate` parameter was
dropped.

**What the authors and the book flag.** The authors list two main weaknesses:
the model relies on a mass flow that may not persist in nature, and the random
rule lets very short strings act as catalysts, which is unrealistic. The paper
also describes a second, ordered rule that sets catalytic strength by matching
strings, but leaves its results to a future paper; Chemart does not implement
it. The book adds that later stochastic studies (Filisetti et al., 2010, 2011)
found autocatalytic sets emerge more rarely under fluctuations than in
deterministic models and often disappear, and that side reactions remain the
main objection to autocatalytic-set theories of the origin of life.

## Further reading

- Kauffman, S. A. (1986). Autocatalytic sets of proteins. *Journal of
  Theoretical Biology* 119, 1–24.
- Farmer, J. D., Kauffman, S. A. & Packard, N. H. (1986). Autocatalytic
  replication of polymers. *Physica D*, 50–67.
- Filisetti, A., Serra, R., Villani, M., Füchslin, R. M., Packard, N. H.,
  Kauffman, S. A. & Poli, I. (2010). A stochastic model of autocatalytic
  reaction networks. *Proceedings of the European Conference on Complex
  Systems (ECCS)*.
- Filisetti, A., Graudenzi, A., Serra, R., Villani, M., De Lucrezia, D.,
  Füchslin, R. M., Kauffman, S. A., Packard, N. H. & Poli, I. (2011). A
  stochastic model of the emergence of autocatalytic cycles. *Journal of
  Systems Chemistry* 2(1).
