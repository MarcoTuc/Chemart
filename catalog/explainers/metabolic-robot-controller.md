## Introduction

This entry is a reaction network used as the control program of a small
robot. The robot's distance sensors pour "sensor substances" into a simulated
reaction vessel. The network's reactions turn them into other substances, and
the concentration of one designated "motor substance" decides what the wheels
do. Nobody writes an if-then rule such as "turn when there is a wall ahead":
whatever the robot does comes from the chemistry's own dynamics. The model
comes from the group of Wolfgang Banzhaf, one of the book's two authors, at the
University of Dortmund. It is an application of artificial chemistry, not a
model of any real cell.

The work came in two steps. Ziegler, Dittrich and Banzhaf (1998) wrote the
networks by hand. They took the bacterium *E. coli* as their model. It steers
towards food and away from toxins by a chain of signalling reactions that link
its receptors to its flagellar motors. With a hand-built enzyme-substrate
network they made a real wheeled robot seek light and avoid obstacles. Ziegler
and Banzhaf (2001) then asked whether such networks can be *evolved* instead,
since, as they wrote, all earlier approaches used "more or less manually
designed" chemistries. They encoded a network as a graph with two kinds of node,
substances and reactions, allowed only graphs that are chemically consistent,
and evolved these graphs with genetic programming (GP), an evolutionary
algorithm for structures of varying size. Fitness was measured by letting each
network drive a simulated Khepera, a small two-wheeled research robot, through
a maze.

What Chemart generates is one network of that kind: a random, chemically
consistent reaction graph wired to the robot's sensors and motor, the raw
material the evolutionary run starts from and works on. The evolution itself
and the robot are not part of it.

Banzhaf and Yamamoto describe this system in the chapter on applications, under
robots controlled by artificial chemistries (book §16.1.3). They end the section
with related work: Penner, Hoar and Jacob's chemistry of bacterial signalling,
and Lones and colleagues' *artificial biochemical networks*, which couple an
artificial genetic network to a metabolic one to control chaotic systems.
Among the catalog's other entries, [Okamoto's biochemical switch](okamoto-switch.md)
and [analog computation with concentrations](analog-function-crn.md) also
compute with concentrations of hand-designed networks, but neither is connected
to a body or evolved. [Algorithmic Chemistry GP](acgp.md) also combines
genetic programming with a chemistry, but there the evolved program is read
as a chemistry, not the other way round. For the genetic side of that later work the book
cites Banzhaf's [artificial regulatory network](arn.md), although the genetic
network of Lones et al. (2014) is a model of its own.

## How it works

### A network as a bipartite graph

In a *bipartite* graph the nodes come in two kinds and every edge joins nodes
of different kinds. Ziegler and Banzhaf use one kind for substances and one for
reactions. An edge from a substance to a reaction means "is consumed by", and
an edge from a reaction to a substance means "is produced by". Every reaction
takes one or two molecules and gives one or two, so there are four types,
numbered by the equations of the 2001 paper:

```
type 4:  s1      -> s1'
type 5:  s1 + s2 -> s1'
type 6:  s1      -> s1' + s2
type 7:  s1 + s2 -> s1' + s2'
```

A third kind of edge points from a substance to a reaction and marks the
substance as a **catalyst** (it speeds the reaction up) or an **inhibitor**
(it slows it down). Ordinary edges carry a default rate constant, written
`k_0`, meaning that the reaction happens spontaneously. Catalytic and
inhibitory edges carry a strength drawn from a range `[k_min, k_max]`.

The paper gives no formula for how a catalyst acts, only that it multiplies
the rate by a factor. Chemart follows the appendix of the 1998 paper, which
writes a catalysed reaction as two parallel mass-action reactions: the plain
one, and a copy with the catalyst on both sides. Under *mass action*, the rate
of a reaction is its constant times the product of its reactants'
concentrations. So a reaction `X -> Y` catalysed by `C` with strength κ runs
at `k_0·[X]·(1 + κ[C])`. No source gives a formula for inhibition, so Chemart
only records which substance inhibits a reaction, and how strongly.

### Material balance: no matter from nothing

Not every random graph is a possible chemistry. The paper's example is the
pair `A -> B + C` and `C -> A`. Run it in a loop and a little `A` turns into
an ever-growing amount of `B`, matter created from nothing. The requirement
that rules this out is **material balance**: it must be possible to give each
substance a positive weight (a "mass") so that every reaction has the same
total weight on both sides. For `A + B -> C + D` that means
`m_A + m_B = m_C + m_D`. For the bad pair above, `m_A = m_B + m_C` and
`m_C = m_A` force `m_B = 0`, so the pair is rejected. Ziegler and Banzhaf apply
this test after every mutation and crossover, so evolution never leaves the
space of consistent networks.

Chemart builds the graph one reaction at a time. It draws a type and random
participants, and keeps the reaction only if a positive weight vector still
exists, which it checks by linear programming. The weights it finds are part
of the output.

### Sensors in, motor out

The 2001 experiment connects the network to the Khepera's eight infrared
proximity sensors, read in pairs. Each pair feeds one substance: `b` from the
left pair, `c` from the front, `d` from the right and `e` from the back. These
four form the **input set**. The inflow of each is proportional to the larger
reading of its pair (readings run from 0, nothing near, to 1023, touching),
and is at most 10% of the reactor's volume per step. A robot boxed in on all
sides therefore receives 0.4 of a volume per step in total.

The single output, the **motor substance** `a`, drives the behaviour through a
threshold. The robot drives straight until the concentration of `a` exceeds
`a_min`. Then one motor reverses, the robot spins on the spot, and `a` is used
up at rate α for as long as it stays above the threshold. This mimics bacterial
chemotaxis, which the paper describes as runs of straight swimming broken by
"tumbling" phases whose frequency is modulated. The paper's settings are
`a_min = 0.1` and α = 0.9. All other substances are internal, and the
evolved network is free to ignore any sensor.

### What the reactor does

In the 2001 paper the reactions run in a 20 × 20 grid folded into a torus
(its edges wrap around). Each cell holds at most one molecule, molecules
diffuse by random steps into empty neighbouring cells, and a reaction fires
when its reactants sit in the same Moore neighbourhood (the eight surrounding
cells). A dilution flux removes random molecules to keep the grid from
filling. The 1998 paper instead used one well-mixed compartment described by
ordinary differential equations (ODEs). Chemart supplies the network (species,
reactions, rate constants, inflows and outflow) for either kind of simulator,
but runs neither.

### A worked example

Here is the network of the default run (its call and summary are under
*Using it in Chemart* below), reaction node by reaction node, from
`net.extras["reaction_graph"]`:

```
r1  type 5  e + s2 -> s3
r2  type 4  b -> s2
r3  type 7  b + c -> d + s2        inhibited by s3 (strength 10)
r4  type 6  e -> 2 a
r5  type 7  2 s2 -> e + s2         catalysed by d (strength 10)
r6  type 5  a + d -> s3
r7  type 7  b + e -> a + c         catalysed by s2 (strength 10)
r8  type 4  s3 -> s1               catalysed by s2 (strength 10)
```

The weights that balance it are `a = 1, b = 2, e = 2, s2 = 2, c = 3, d = 3,
s1 = 4, s3 = 4`. Check `r7`: `b + e` weighs 2 + 2 = 4, and `a + c` weighs
1 + 3 = 4. The three catalysed nodes become two network reactions each, which
is why the summary counts 11 reactions for 8 nodes. For example, `r5` is
implemented as `2 s2 -> e + s2` with `k = 1` plus
`2 s2 + d -> e + s2 + d` with `k = 10`.

Now put a wall behind the robot and nothing else. Only `e` flows in, at 0.1
per unit time. Nothing makes `b` or `s2`, so `r1` and `r7` cannot fire. Every
`e` becomes two `a` through `r4`, and `a` drains at 0.9·[a]. The motor
substance settles where production equals drain, at `2 × 0.1 / 0.9 ≈ 0.222`,
above the threshold of 0.1: this random network makes the robot spin when
something is behind it. A wall in front, on the other hand, feeds only `c`,
which is consumed only by `r3` together with `b`. So `c` piles up, `a` stays
at zero, and the robot drives into the wall. A random network is not yet a
controller. The evolutionary search's job is to find the networks that are.

## Using it

The default run is one random 8-substance, 8-reaction graph with the 2001
paper's inputs, output, α and `a_min`, and all eight sensors at 1023, the
"surrounded by obstacles" case the paper uses to state the maximum inflow.
It is a possible starting individual of the evolutionary run, not the evolved
controller: the paper only draws the evolved network's kernel and does not
list its reactions, so it cannot be rebuilt.

What to read in the result:

- `net.inflow`: the inflow of each sensor substance for the given readings
  (here 0.1 each), and `net.outflow`: α for the motor substance.
- `net.extras["reaction_graph"]`: one entry per reaction node, with its
  type (4 to 7), reactants, products, its modifier (catalyst or inhibitor
  with strength) if any, and the indices of the network reactions that
  implement it.
- `net.extras["conservation"]`: the material-balance weights, and whether
  they are unique.
- `net.extras["input_set"]`, `["output_set"]`, `["sensors"]`,
  `["actuators"]`: which substances are wired to which sensors, and the
  threshold and consumption rate of each actuator. The threshold gate is not
  a rate law, so a simulator has to apply it itself.
- Inhibition appears only as `inhibitor` and `k_inhibition` keys on a
  reaction's rate. A simulator must choose how to apply it.

**Changing what the robot sees.** Pass `sensor_readings`, by sensor name.
Unlisted sensors read 0. The following script integrates the default network
as a one-compartment ODE, with mass action, catalysts as above, inhibitors
ignored and the actuator drain applied at all concentrations. It prints the
motor substance after 200 time units for five situations:

```python
import numpy as np
from scipy.integrate import solve_ivp
import chemart

def steady_a(readings, seed=1, t_end=200.0):
    net = chemart.generate_network("metabolic-robot-controller", seed=seed, sensor_readings=readings)
    ids, R, P = net.matrices()
    S, R = (P - R).toarray(), R.toarray()
    k = np.array([r.rate["k"] for r in net.reactions])       # inhibitors ignored
    inflow = np.array([net.inflow.get(s, 0.0) for s in ids])
    outflow = np.array([net.outflow.get(s, 0.0) for s in ids])
    def f(t, x):
        v = k * np.prod(np.maximum(x, 0)[:, None] ** R, axis=0)   # mass action
        return S @ v + inflow - outflow * x
    sol = solve_ivp(f, (0, t_end), np.zeros(len(ids)), method="LSODA")
    return round(float(sol.y[ids.index("a"), -1]), 3)

print("no obstacle   ", steady_a({}))
print("wall in front ", steady_a({"front1": 1023, "front2": 1023}))
print("wall on left  ", steady_a({"left1": 1023, "left2": 1023}))
print("wall behind   ", steady_a({"back1": 1023, "back2": 1023}))
print("surrounded    ", steady_a({s: 1023 for s in ["left1", "left2", "front1", "front2",
                                                   "right1", "right2", "back1", "back2"]}))
```

```
no obstacle    0.0
wall in front  0.0
wall on left   0.121
wall behind    0.222
surrounded     0.167
```

The values for the side and back walls do not change between 200 and 400 time
units. The surrounded case is still drifting slightly (0.165 at 400), because
the dead-end substance `s1` keeps accumulating: the motor substance is the
only way out of the network. The paper's lattice reactor avoids this build-up
with its random dilution, which the generated network does not include.

**The 1998 two-motor wiring.** In the 1998 robot, two substances drove the
left and right wheels directly. Map two output substances to the
`left-motor` and `right-motor` actions:

```python
net = chemart.generate_network("metabolic-robot-controller", seed=1,
                               actuator_map={"ML": "left-motor", "MR": "right-motor"},
                               n_substances=10, n_reactions=10)
net.extras["output_set"], net.outflow   # (['ML', 'MR'], {'ML': 0.9, 'MR': 0.9})
```

This is only the wiring. The 1998 hand-designed enzyme-substrate network,
with its light sensors and constant-held species, is not generated.

**Size.** The evolved graph of the 2001 paper had nearly 50 nodes and 70
edges. `n_substances=20, n_reactions=30` gives a graph of that order (20
species, 34 network reactions) in under 0.1 s. Large graphs are slow and can
fail. Once the balancing weights are unique, a new reaction is kept only if it
happens to balance under exactly those weights, and random draws rarely do.
With `n_substances=200, n_reactions=200` the run takes about 30 s, and asking
for 400 reactions among 200 substances fails after 50,000 draws with only 238
found.

**Restricting reaction types.** `reaction_types=[4]` allows only one-to-one
conversions. Every weight is then 1, and material balance holds trivially.

## Results

**Hand-designed controllers work on a real robot (1998).** Ziegler, Dittrich
and Banzhaf compared two hand-written chemistries as information processors.
One was a polymer chemistry in which strings join end to end. The other was an
enzyme-substrate network, of the kind known as "chemical neuronal networks". For
two inputs, the polymer chemistry's steady-state output varied smoothly. The
enzyme-substrate network instead split its input space into two clearly
separated regions and acted as a logic AND gate. The authors argued that this
sigmoidal (switch-like) response corrects small fluctuations when a signal
passes from one sub-network to the next. An enzyme-substrate network of this
kind, fed by eight proximity and eight ambient-light sensors, controlled a real
robot. It sought light and avoided obstacles, let obstacle avoidance override
light seeking when they conflicted, and was robust against small random amounts
of substance injected into the reactor. The paper gives this as a qualitative
demonstration, with an example trajectory. Chemart does not reproduce it: it
generates random graphs, not the 1998 network.

**Evolved controllers avoid obstacles (2001).** Ziegler and Banzhaf evolved
networks with a (µ, λ) strategy, in which the µ = 50 best of λ = 250 offspring
become the next parents. Crossover probability was 0.7, mutation probability
0.9, and runs lasted at most 500 generations. Mutations add or delete
catalysts, inhibitors, reactions or substances, rescale a catalytic strength,
or change a reaction's type. Crossover cuts a *metabolic pathway* out of one
graph, meaning the shortest path from one substance to another together with
every participant of the reactions on it, and inserts it into another graph.
Each network was run for 1,200 reactor cycles while it drove the simulated
robot. Its fitness summed the differences between the two motors' speed and
direction over the run, so the best networks drive straight for as long as
possible and turn only briefly. The run produced networks that steer the robot
through the maze. The book's Figure 16.8 shows the robot approaching a wall,
turning around in its preferred direction (left) and driving off the other
way. Chemart implements none of this machinery. Its tests check the building blocks: every generated graph has positive
balancing weights, every node is one of the four types, the paper's example
`A -> B + C`, `C -> A` is rejected, and the sensor inflow gives 0.4 in total
when surrounded.

**The same network drives a real robot.** The network evolved in simulation
was then run on a real Khepera in a 1 m² maze. It was "completely unchanged":
only time-scale parameters were adjusted, and the sensor inflow was amplified
by an extra factor because the real robot moves faster. Photographs show it
meeting a wall and turning counterclockwise. The book's summary says instead
that evolution took place in real hardware and that "fundamentally the same
network (up to scale factors) would evolve in a simulated and the real world".
The paper describes a transfer of one evolved network, not two evolutions.
There is nothing in Chemart to reproduce here.

**The evolved kernel uses only the front and left sensors.** The paper
defines the *kernel* of a graph as the substances that inputs can influence
through paths in the graph. A node outside it never changes concentration and
can be deleted without effect. The best network's graph had nearly 50 nodes
and 70 edges, but its kernel was tiny. Only the front sensors' substance and
the left sensors' substance fed the motor substance, with one catalytic link.
The front input makes the robot turn at a wall. The left input looks
superfluous, but because the robot always turns counterclockwise, the left
sensors see the obstacle a moment into the turn. Their substance keeps `a`
above the threshold, so the turn is not cut short by a brief dip. The
authors compare this with *E. coli*, whose receptors sit mostly at its front
end. Chemart does not compute the kernel (see the implementation decisions above) and has
no evolved network to analyse.

**Genomes grow, kernels do not.** The paper's Figure 11 follows the number
of nodes per graph over the run: "the minimum length of the functional code
remains constant during the evolution, but the overall length increases". The
authors present this as the "well-known phenomena" of GP, *introns* (parts of
the genome with no effect) and *bloat* (growth of the genome without a gain in
fitness). Chemart has no evolutionary run, so this is not reproduced.

**Later work: artificial biochemical networks.** Lones, Fuente, Turner and
colleagues (2014) evolved three kinds of network as controllers: an
artificial genetic network, an artificial metabolic network, and a coupled
network in which the genetic network controls the metabolic one. The tasks
were steering trajectories of the chaotic Lorenz system, moving a
trajectory across Chirikov's standard map (a simple map that can be chaotic)
in as few steps as possible, and generating the gait of a simulated
four-legged robot. The genetic networks often, though not always, did better
than the metabolic ones. In the metabolic networks, conserving mass made
evolution less sensitive to poorly chosen parameters. Coupling helped on some
tasks and not on others. Networks built from non-linear discrete maps (such as the logistic map) gave
the best controllers on the two hardest tasks. The book describes this line of
work as combining metabolic, signalling and genetic networks. In the 2014
paper, the coupled model is a genetic network controlling a metabolic one;
signalling networks appear only in its survey of related work. These models
are not Ziegler and Banzhaf's reaction graphs, and Chemart implements none of
them.

## Further reading

- Lones, M. A., Fuente, L. A., Turner, A. P., Caves, L. S. D., Stepney, S.,
  Smith, S. L. & Tyrrell, A. M. (2014). Artificial biochemical networks:
  evolving dynamical systems to control dynamical systems. *IEEE Transactions
  on Evolutionary Computation* 18(2), 145–166. Preprint:
  <https://www.macs.hw.ac.uk/~ml355/common/papers/lones-tevc2012-ABNs.pdf>
- Lones, M. A., Tyrrell, A. M., Stepney, S. & Caves, L. S. (2010).
  Controlling complex dynamics with artificial biochemical networks. In
  Esparcia-Alcázar, A. et al. (eds.), *Genetic Programming*, LNCS 6021,
  159–170. Springer.
- Penner, J., Hoar, R. & Jacob, C. (2003). Modelling bacterial signal
  transduction pathways through evolving artificial chemistries. In
  *Proceedings of the First Australian Conference on Artificial Life*,
  Canberra.
