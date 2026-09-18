# Metabolic / artificial biochemical network robot controllers

`metabolic-robot-controller` · *Ziegler & Banzhaf, 2001; Lones, Fuente, Turner et al.*

A reaction network used as a brain. Sensors feed substances into the network, the network's own reactions process them, and a motor consumes the concentration of a designated actuator species - so the robot's behaviour *is* the network's dynamics. Hand-designed networks do light-seeking and obstacle avoidance; evolved ones transfer from simulation to real hardware with little more than a scaling factor. The appeal is that computation and its substrate are not separable here.

| | |
|---|---|
| **family** | application |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 16.1.3 |
| **refs** | [958], [959], [519], [520], [68], [660] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `mass-conservation`, `flow` |

## Molecules, reactions, reactor

**S — molecules** (explicit): abstract substances s1, s2, ...; sensor substances (b, c, d, e) receive sensor-driven inflow, the actuator substance a is read and consumed by the motor

**R — reactions** (explicit, arity [1, 2]): Bipartite directed weighted graph G = (V, E), V = R u S. Reaction nodes have the types s1 -> s1' ; s1 + s2 -> s1' ; s1 -> s1' + s2 ; s1 + s2 -> s1' + s2' (paper eqs. 4-7). Stoichiometric edges carry the default weight k_0; a special edge S -> R carries a catalytic or inhibitory weight k in [k_min, k_max]. A catalysed node becomes X -> Y [k_0] plus X + C -> Y + C [k_0 k]. The graph must fulfil material balance: positive species weights m with M(R) m = 0.

**A — reactor**: lattice-2d, ode
 · *dilution:* sensor substances flow in (eq. 30); the actuator substance is consumed at rate alpha above a_min (eq. 31); the 2001 lattice reactor also removes random molecules

## What you get

```python
net = chemart.generate_network("metabolic-robot-controller", seed=1)
```

```
metabolic-robot-controller: 8 species, 11 reactions, status=complete
provides: catalysts, flow, mass-conservation, rate-constants, stoichiometry, topology
seed: 1
extras: actuators, conservation, input_set, output_set, reaction_graph, sensors
```

First reactions:

```
e + s2 -> s3  [mass-action k=1.0 node=r1]
b -> s2  [mass-action k=1.0 node=r2]
b + c -> d + s2  [mass-action k=1.0 node=r3 inhibitor=s3 k_inhibition=10.0]
e -> 2 a  [mass-action k=1.0 node=r4]
2 s2 -> e + s2  [mass-action k=1.0 node=r5]
2 s2 + d -> e + s2 + d  [mass-action k=10.0 node=r5 catalyst=d k_catalysis=10.0]
a + d -> s3  [mass-action k=1.0 node=r6]
b + e -> a + c  [mass-action k=1.0 node=r7]
… and 3 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `n_substances` | `int` | `8` | structural | number of molecule nodes \|S\|, including the sensor and actuator substances; the others are named s1, s2, ... (the paper numbers them; the prefix keeps '2 s2' readable) <br>`2` … `500` · *range:* must include the input and output substances; the evolved graph of the paper's fig. 13 has nearly 50 nodes (molecules + reactions) and 70 edges |
| `n_reactions` | `int` | `8` | structural | number of reaction nodes \|R\|, each kept only if the graph stays materially balanced <br>`0` … `1000` |
| `reaction_types` | `list` | `[4, 5, 6, 7]` | structural | allowed reaction types by the paper's equation number: 4 s1 -> s1', 5 s1 + s2 -> s1', 6 s1 -> s1' + s2, 7 s1 + s2 -> s1' + s2' (Table 1 uses all four) |
| `p_modifier` | `float` | `0.5` | structural | probability that a reaction node carries a catalytic or inhibitory edge (at most one per node, as in mutation case 1); not published <br>`0` … `1` |
| `p_inhibitor` | `float` | `0.5` | structural | probability that a modifier edge is inhibitory rather than catalytic; not published <br>`0` … `1` |
| `k_0` | `float` | `1.0` | kinetic | default weight k_default of stoichiometric edges, the rate constant of the spontaneous reaction (1.0 in Ziegler, Dittrich & Banzhaf 1998) <br>≥ `0` |
| `k_min` | `float` | `10.0` | kinetic | lower bound of catalytic/inhibitory weights, drawn uniformly in [k_min, k_max] <br>≥ `0` |
| `k_max` | `float` | `10.0` | kinetic | upper bound of catalytic/inhibitory weights <br>≥ `0` · *range:* the 2001 paper gives no values; the default k_min = k_max = 10 is the catalytic efficiency kappa of the 1998 paper |
| `sensor_map` | `dict` | `{'b': ['left1', 'left2'], 'c': ['front1', 'fron…` | structural | input set I: each sensor substance and the Khepera proximity sensors whose maximum drives its inflow (eq. 29) |
| `sensor_readings` | `dict` | `{'left1': 1023, 'left2': 1023, 'front1': 1023, …` | population | current sensor values in [0, 1023] by sensor name; unlisted sensors read 0. The default is the paper's robot surrounded by obstacles |
| `max_inflow` | `float` | `0.1` | population | maxInflow of eq. 30: inflow of a sensor substance at maximum reading, as a fraction of the reactor volume per iteration (10% in the paper) <br>`0` … `1` |
| `actuator_map` | `dict` | `{'a': 'rotate'}` | structural | output set O: each actuator substance and its action, one of rotate (2001: switch one motor's direction above a_min), left-motor, right-motor (1998: substance drives a wheel) |
| `alpha` | `float` | `0.9` | population | consumption rate alpha of the actuator substance while it is above a_min (Table 2: 0.9) <br>`0` … `1` |
| `a_min` | `float` | `0.1` | population | actuator activation threshold a_min (Table 2: 0.1) <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- light seeking and obstacle avoidance on a Khepera robot (1998, hand-designed enzyme-substrate network)
- evolved networks transfer between simulation and a real Khepera, with only a scale factor on sensor inflow
- the evolved kernel uses only the front and left sensor substances to raise the motor substance (fig. 14)
- genome length grows during evolution (introns, bloat) while the kernel stays small
- combined genetic + metabolic + signalling networks controlling chaotic dynamical systems (Lones et al.)

## Sources

- Ziegler, J. & Banzhaf, W. (2001). Evolving control metabolisms for a robot. Artificial Life 7(2):171-190. Reaction graphs sec. 3 (eqs. 4-17), material balance sec. 3.1.1, mutation sec. 3.2.1, sensors and actuators sec. 5.1-5.2 (eqs. 29-31), Tables 1-2. http://www.cs.mun.ca/~banzhaf/papers/metabolism.pdf
- Ziegler, J., Dittrich, P. & Banzhaf, W. (1998). Towards a metabolic robot control system. In Holcombe & Paton (eds.), Information Processing in Cells and Tissues, pp. 305-318. Appendix, example A: mass-action form of catalysed reactions and k = 1.0, kappa = 10. https://users.fmi.uni-jena.de/~dittrich/p/ZDB97ipcat.ps

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The network is one individual of the GP system: a random reaction graph. The (mu, lambda) evolution, mutation, crossover, pathway extraction and the robot simulation act on many networks and are not generated; the v1 evolution parameter is dropped.
- Reaction nodes are drawn like mutation case 4 (type uniform over reaction_types, participants uniform over all molecules, with replacement). No-op reactions and duplicate nodes are rejected. The paper's 'grow' initialisation is not described further.
- Material balance (sec. 3.1.1) is tested as the existence of weights m >= 1 with M(R) m = 0, by linear programming. The paper's steps 4-6 set the independent weights to unity and check the rest, which can wrongly reject a balanced network whose positive solution needs other values. A node that would break the balance is rejected rather than penalised. Once the weight vector is unique, a node is kept iff it is orthogonal to it.
- Extraction garbling: eqs. 4-7 lost their primes (printed 's1 -> s1'); they are read as s1 -> s1', etc., since the unprimed forms would be no-ops.
- Catalysis: the 2001 paper only says a catalyst multiplies the rate constant by a factor. The 1998 appendix writes a catalysed reaction as a parallel mass-action channel w5 = k5 kappa [A][B][B] with k5 = k1, so a catalysed node becomes the plain reaction (k_0) plus a channel with the catalyst on both sides (k_0 k). The kappa glyph is lost in the extracted text and is restored from the table 'kappa 10'.
- Inhibition has no published functional form. An inhibited node is one mass-action reaction with k_0 whose rate dict records inhibitor and k_inhibition; a simulator must choose the inhibitory factor.
- Modifier edges: at most one per node, as mutation case 1 implies ('adding an inhibitor/catalyst if none is present'). The modifier is any molecule, with k ~ U[k_min, k_max]. p_modifier and p_inhibitor are not in the paper.
- Rate defaults: the 2001 paper gives no k_default, k_min or k_max. The defaults k_0 = 1.0 and k_min = k_max = 10 are the 1998 paper's k_i and kappa.
- Sensor inflow (eq. 30) is expressed per reactor iteration as a fraction of the reactor volume (cellsize = 1), so max_inflow = 0.1 gives 0.4 in total when all four sensor pairs read 1023, as the paper states. sensor_map and actuator_map replace the v1 callables.
- Actuator outflow (eq. 31) is recorded as outflow alpha. The threshold gate [a] > a_min is not a rate law and is kept in extras.actuators. Reactor size, fill level and cycle counts (Table 2) only size the simulation and are dropped. No initial state is set: the paper fills 50% of the lattice without giving the composition.
- The functional set (kernel, eq. 16) and pathways P(v, w) are analyses of a graph, not part of it, and are not computed. extras.reaction_graph lists each reaction node with its type, modifier and the indices of the network reactions that implement it.

## Notes

The paper also defines a shortest-path 'metabolic pathway' extraction, the subgraph P(v, w) with every participant of each reaction on the path, used by crossover. It would fit as a Chemart analysis on extras.reaction_graph.

---

*Specification: `catalog/chemistries/metabolic-robot-controller.yaml` · generator: `chemart/chemistries/metabolic_robot_controller.py` · tests: `tests/chemistries/test_metabolic_robot_controller.py`*
