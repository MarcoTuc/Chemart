# Bagley & Farmer autocatalytic metabolism

`bagley-farmer` · *Bagley & Farmer, 1992; Bagley, Farmer & Fontana, 1992*

What Kauffman's autocatalytic sets look like once you give them kinetics. Polymers split and join, and any polymer can catalyse a reaction, multiplying both directions of it. Run in a chemostat - food flows in, everything washes out - and mass concentrates onto a small mutually-catalysing core. The crucial detail is that catalysis alone cannot do this: at equilibrium a catalysed and an uncatalysed network reach the same place, so the focusing only appears away from equilibrium, driven by the flow.

| | |
|---|---|
| **family** | origin-of-life |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 6.3.1, 7.3.1 |
| **refs** | [49], [50] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `mass-conservation`, `flow`, `initial-state`, `thermodynamic-consistency` |

## Molecules, reactions, reactor

**S — molecules** (implicit): polymers = strings over an alphabet of 1-4 letters (a, b for the paper's binary chemistry) up to max_length; water H; with saturation, one bound pool x_bound per polymer (sum of the enzyme complexes x is part of, eq. 26)

**R — reactions** (implicit, arity [1, 2, 3]): A + B <-> C + H for every split of every polymer (kf forward, kr backward). A catalytic link E multiplies both directions by (1 + nu E) (eq. 9): A + B + E -> C + H + E at nu kf and C + H + E -> A + B + E at nu kr. With saturation (appendix, eqs. 22-29): A + B + E -> C_bound + E_bound + H at nu kf, C + H + E -> A_bound + B_bound + E_bound at nu kr, and x_bound -> x at ku.

**A — reactor**: ode
 · *dilution:* chemostat: each food species enters at flux delta, every species leaves at first-order rate K = delta * (monomers in the food set) / m0 (eqs. 10-12)

## What you get

```python
net = chemart.generate_network("bagley-farmer", seed=1)
```

```
bagley-farmer: 125 species, 696 reactions, status=complete
provides: flow, initial-state, mass-conservation, rate-constants, stoichiometry, topology
seed: 1
extras: buffered, catalytic_links, conservation, flow_law, food_set
```

First reactions:

```
2 a -> aa + H  [mass-action k=649.0]
aa + H -> 2 a  [mass-action k=2.5]
a + b -> ab + H  [mass-action k=649.0]
ab + H -> a + b  [mass-action k=2.5]
b + a -> ba + H  [mass-action k=649.0]
ba + H -> b + a  [mass-action k=2.5]
2 b -> bb + H  [mass-action k=649.0]
bb + H -> 2 b  [mass-action k=2.5]
… and 688 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `alphabet_size` | `int` | `2` | structural | number of distinct monomers m (letters a, b, c, d); the paper's simulations use 2 <br>`1` … `4` |
| `max_length` | `int` | `5` | structural | longest polymer in the network (at most 512 polymers); the paper has no hard bound and relies on the threshold instead <br>`2` … `8` |
| `food_set` | `list` | `['a', 'b']` | population | polymers supplied from outside at flux delta each <br>*range:* Figures 5-12 use {a, b}; Figure 3 uses {a, b, ab, ba} |
| `p` | `float` | `0.01` | structural | probability that a given polymer strongly catalyses a given condensation/cleavage pair (ignored when links is given) <br>`0` … `1` · *range:* Figure 13 varies the number of links from 0 to 200 on a 15-reaction network |
| `links` | `list` | `` | structural | explicit catalytic links 'A + B <-> C \| E' (e.g. 'a + b <-> ab \| bb'); when non-empty they replace the random draw |
| `nu` | `float` | `897000.0` | kinetic | catalytic efficiency: a catalyst E multiplies both rate constants by (1 + nu E) <br>≥ `0` · *range:* Table 2: 8.97e5 (network a), 5.26e5 (network b); Figure 3: 1e4 |
| `nu_distribution` | `enum` | `constant` | stochastic | constant: every link has efficiency nu; uniform: each link draws its efficiency from U(0, 2 nu) <br>one of `constant`, `uniform` |
| `kf` | `float` | `649.0` | kinetic | condensation rate constant, the same for every reaction <br>≥ `0` · *range:* Table 2: 6.49e2 (a), 3.02e4 (b); Figure 3: 1e2 |
| `kr` | `float` | `2.5` | kinetic | hydrolysis (cleavage) rate constant, the same for every reaction <br>≥ `0` · *range:* Table 2: 2.50 (a), 2.70 (b); Figure 3: 10 |
| `ku` | `float` | `50000.0` | kinetic | unbinding rate of enzyme complexes (used only with saturation) <br>≥ `0` · *range:* Table 2: 5.00e4 (a), 7.11e4 (b); Figure 3: 1e4 |
| `saturation` | `bool` | `True` | kinetic | true: catalysis through bound complexes (appendix, eqs. 22-29), as in the paper's simulations; false: the unsaturated approximation of eq. 9 |
| `H` | `float` | `1.0` | population | water concentration, held constant (Figure 3 uses H = 1) <br>≥ `0` |
| `delta` | `float` | `17.9` | population | inflow flux of each food species (mass flux driving the system from equilibrium) <br>≥ `0` · *range:* Table 2: 1.79e1 (a), 1.41e2 (b); Figure 5 uses 0.01, 1e5 and 10^7.5 |
| `m0` | `float` | `2.0` | population | total monomer concentration at the chemostat fixed point; sets the outflow rate K and the initial food concentrations <br>≥ `1e-12` · *range:* 2.0 in Figures 5-14; 3 in Figure 3 |
| `threshold` | `float` | `` | population | metadynamics concentration threshold (one molecule in the reaction vessel); 0 keeps the whole bounded network, > 0 returns the graph at the metadynamical fixed point <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- catalytic focusing: mass concentrates into a small core set, but only away from equilibrium (eq. 15; Figures 5-8)
- catalysis does not change the equilibrium: at delta = 0 every catalysed and uncatalysed reaction reaches the same ratio C / AB = kf / (kr H)
- the steady-state profile falls off exponentially with length at equilibrium and flattens for the network at intermediate delta (Figures 5-7)
- the deterministic metadynamics always reaches a unique fixed point
- metabolisms jump between fixed points under spontaneous fluctuations, a first sign of evolution ([50])

## Sources

- Bagley, R. J. & Farmer, J. D. (1992). Spontaneous emergence of a metabolism. In Langton, Taylor, Farmer & Rasmussen (eds.), Artificial Life II, SFI Studies in the Sciences of Complexity X, pp. 93-140. Addison-Wesley. Also LA-UR-91-1707. Scan: https://oms-inet.files.svdcdn.com/production/files/SpontaneousEmergence.pdf (listed at https://www.inet.ox.ac.uk/publications/spontaneous-emergence-of-a-metabolism)
- Bagley, R. J., Farmer, J. D. & Fontana, W. (1992). Evolution of a metabolism. Artificial Life II, pp. 141-158. Scan: https://oms-inet.files.svdcdn.com/production/files/EvolutionMetabolism.pdf (listed at https://www.inet.ox.ac.uk/publications/evolution-of-a-metabolism)

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Kinetics follow the paper: spontaneous A + B <-> C + H with the same kf, kr for all reactions (3.1.3), catalysis multiplying both directions by (1 + nu E) (eq. 9), and saturation through complexes that unbind at ku, aggregated into one bound pool per species (eqs. 22-29). All of it is mass action; the catalyst and nu are recorded in each catalysed rate dict.
- Eq. 24 routes the spontaneous flux through the complex CE; this looks like a typo (a spontaneous reaction uses no enzyme), so spontaneous reactions produce free species directly. Eqs. 10-11 and 25-29 agree with this reading.
- The paper keeps the reaction graph finite with a concentration threshold and lumps the spontaneous background into per-length aggregates (eqs. 19-21). Here the chemistry is bounded by max_length and every spontaneous reaction is explicit, which is the exact system those aggregates approximate. The threshold metadynamics of section 3.2 is available with threshold > 0 (default off; the paper gives no numeric threshold).
- In metadynamics a species is above threshold when its free plus bound concentration is; a reaction is in the graph when its free reactants and catalyst are above threshold (food is always above). Each round integrates the chemostat ODE of the current graph to its fixed point (t = 200/K) instead of the paper's algebraic solver. The final graph is returned with status truncated, and the fixed point, active species, rounds and outcome in extras.analysis.
- Because every spontaneous reaction is explicit, a hard threshold can make the graph updates cycle: removing a species also removes the reactions that consume it, so it rises back above threshold. The paper reports a unique fixed point for its lumped background and says nothing about cycles. Here a repeated active set merges all sets in the cycle, and the merged graph's fixed point is returned with outcome cycle-merged ('fixed' when the updates stop, 'round-limit' after 50 rounds).
- Outflow: eq. 12 writes dm/dt = N_f delta - K m. With polymers in the food set this is ambiguous, so K = delta * (total monomers in the food set) / m0, which makes m0 the true monomer mass at the fixed point and equals the paper's rule for monomer food sets. The run starts with that mass in the food set, as in the paper (m(0) = m0).
- Eq. 15: solving eqs. 13-14 gives [ab]/[ba] = (1 + beta)/(1 + beta/gamma) with beta = K/(kr H). The printed beta = delta/(m0 kr H) matches only when m0 is the mass per food species (it is off by N_f = 2 otherwise); tests use K/(kr H).
- Defaults are Table 2 (a) / Figure 6: kf = 6.49e2, kr = 2.50, nu = 8.97e5, ku = 5.00e4, delta = 1.79e1, m0 = 2.0 (the caption of Figure 6 confirms delta = 1.79e1). H = 1 comes from Figure 3, since Figures 5-14 do not state it. The paper gives no value of p, so the default p = 0.01 is a Chemart choice (about 120 links on the default 62 polymers, close to the paper's 118-link network).
- Table 2 (b) prints ranges nu in [<1, 1e4] and ku in [1e7, >1e10], which exclude their own optima (5.26e5 and 7.11e4); the two ranges look swapped. Only the optima are used.
- nu_distribution uniform draws U(0, 2 nu): the paper uses 'uniform within given maximum and minimum values' without giving them. Any polymer, monomers included, may catalyse (the random rule of 3.1.1 excludes nothing); the string-matching rule of 3.1.2 is not implemented.
- The v1 parameters firing_disk, food_influx, catalysis_probability, k_spontaneous, catalytic_enhancement and concentration_threshold became food_set, delta, p, kf, nu and threshold. mutation_rate (the stochastic metadynamics of the companion paper [50]) is dropped: it drives a simulation over fixed points rather than defining the network.

## Notes

Defines an autocatalytic metabolism dynamically: an autocatalytic set whose concentrations differ significantly from the uncatalysed case. The catalytic topology is the same as `kauffman-autocatalytic-sets`, so it is Chemart's canonical demo of "same topology, different kinetics, different answer".

---

*Specification: `catalog/chemistries/bagley-farmer.yaml` · generator: `chemart/chemistries/bagley_farmer.py` · tests: `tests/chemistries/test_bagley_farmer.py`*
