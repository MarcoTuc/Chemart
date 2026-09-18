# Arrhenius-gated collision algorithm

`energy-gated-collision` · *Banzhaf & Yamamoto, 2015 (sec. 18.3.3), extending the collision algorithm of sec. 2.3.3*

Not a chemistry but the reactor algorithm underneath many of them, catalogued so its behaviour can be checked. Draw two molecules at random, look up their reaction, and let it happen only if the collision energy clears the activation barrier. The point is what that gate buys you: the measured fraction of effective collisions comes out as the Arrhenius factor exp(-Ea/RT), so an Arrhenius plot of the simulated rate is straight and its slope recovers the activation energy you put in.

| | |
|---|---|
| **family** | systems-biology |
| **kind** | analysis |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `book+decisions` — the book left gaps; each filled choice is listed below |
| **book** | 18.3.3; algorithm from 2.3.3, energetics from 2.2.6 (eqs. 2.28-2.33) |
| **refs** | [42], [757], [780], [870] |
| **provides** | `topology`, `stoichiometry`, `rate-constants`, `rate-law`, `energies`, `thermodynamic-consistency`, `initial-state` |

## Molecules, reactions, reactor

**S — molecules** (explicit): any molecules for which a potential energy can be given and a kinetic energy computed; here a small explicit set of species, each with a potential energy in kJ/mol

**R — reactions** (implicit, arity 2): pick two random molecules; look up the rule for that unordered pair; if the collision energy reaches the activation energy Ea of the reaction, replace the reactants by the products (effective collision), otherwise leave the vessel unchanged (elastic collision). The acceptance probability is the Arrhenius factor exp(-Ea/RT) of eq. 2.32, so the reaction fires at k = A exp(-Ea/RT) without any propensity ever being summed

**A — reactor**: well-stirred-multiset
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("energy-gated-collision", seed=1)
```

```
energy-gated-collision: 4 species, 2 reactions, status=observed
provides: energies, initial-state, rate-constants, rate-law, stoichiometry, topology
seed: 1
extras: analysis, energies
```

First reactions:

```
X2 + X1 -> Y1 + Y2  [arrhenius A=1.0 Ea=2.0 T=300.0 R=0.00831451 units=kJ/mol]  (x212)
Y1 + Y2 -> X1 + X2  [arrhenius A=1.0 Ea=5.0 T=300.0 R=0.00831451 units=kJ/mol]  (x81)
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `temperature` | `float` | `300.0` | thermodynamic | absolute temperature T inside the reactor, in kelvin (eq. 2.32) <br>`1.0` … `5000.0` · *range:* the gate depends on T only through Ea/RT; RT = 2.494 kJ/mol at 300 K |
| `Ea` | `float` | `2.0` | thermodynamic | activation energy of the forward reaction of the default system, in kJ/mol; must be at least max(0, delta_G) <br>`0.0` … `500.0` · *range:* real barriers are 40-400 kJ/mol, which at 300 K would reject all but ~1e-7 of the collisions; the default keeps Ea/RT below 1 so a short run measures the Arrhenius factor |
| `delta_G` | `float` | `-3.0` | thermodynamic | free-energy change of the forward reaction, in kJ/mol: it is the potential energy of the product well, and it fixes the reverse barrier at Ea - delta_G (fig. 2.4). Negative is exergonic/spontaneous <br>`-500.0` … `500.0` |
| `A` | `float` | `1.0` | kinetic | pre-exponential factor of eq. 2.32: the rate constant a barrierless reaction would have, i.e. the collision frequency factor. It scales every rate and cancels out of the acceptance measurement <br>`0.0` … `1000000000.0` |
| `energy_model` | `enum` | `bath` | thermodynamic | bath: a thermostat draws the line-of-centres collision energy from Exponential(RT), so the acceptance is exactly exp(-Ea/RT) and the reaction heat comes from the bath. conserved: every molecule carries a kinetic energy, the gate compares Ea with the sum of the two (the book's literal wording), and E = Ek + Ep is conserved exactly (eq. 2.28); the acceptance is then (1 + Ea/RT) exp(-Ea/RT) <br>one of `bath`, `conserved` |
| `population` | `int` | `400` | population | number of molecules M in the well-stirred vessel, split over the system's initial species <br>`2` … `200000` · *range:* the chapter-2 examples use M = 100 to 100,000 |
| `steps` | `int` | `4000` | stochastic | collisions to attempt, elastic ones included; M collisions are one generation (sec. 2.6.1) <br>`0` … `10000000` · *range:* 10 generations at the default M; the acceptance estimate tightens as 1/sqrt(attempts) |
| `system` | `dict` | `` | structural | run the algorithm on a user-supplied system instead of the default one: {'energies': {species: potential energy in kJ/mol}, 'reactions': [{'id', 'reactants' (exactly two), 'products', 'Ea'}], 'initial': {species: proportion}}; empty means X1 + X2 <-> Y1 + Y2 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- the measured firing frequency of a reaction is its Arrhenius factor: effective collisions / collisions of its reactants = exp(-Ea/RT), within sampling error, over a range of Ea and T
- an Arrhenius plot of the measured rate constant is straight: ln k against 1/T has slope -Ea/R and intercept ln A, recovering the barrier that was put in (eq. 2.32)
- raising Ea or lowering T suppresses the reaction exponentially; the two enter only through Ea/RT
- at Ea = 0 the gate always opens and the algorithm degenerates to the naive collision algorithm of sec. 2.3.3, whose statistics are mass action: a second-order decay follows 1/n = 1/n0 + 2 P s / (M (M-1)) in s collisions
- forward and reverse barriers that differ by delta_G drive the vessel to the equilibrium constant K = exp(-delta_G/RT) of eq. 2.33
- the cost per collision is independent of the number of reactions (only the rules of the drawn pair are looked up), but the acceptance rate is exp(-Ea/RT): the book's caveat that the algorithm only pays off when elastic collisions are rare

## Sources

- P. Atkins and J. de Paula. Physical Chemistry. Oxford University Press, 2002 - the book's own ref. [42] for eq. 2.32, and the source of the collision-theory result used here: in a Maxwell-Boltzmann gas the fraction of collisions whose kinetic energy along the line of centres exceeds Ea is exactly exp(-Ea/RT)
- A. Slepoy, A. Thompson and S. Plimpton. A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks. J. Chem. Phys. 128:205101, 2008. https://doi.org/10.1063/1.2919546 - the book's ref. [780], the composition-rejection algorithm that 18.3.3 contrasts with this one (constant time, but still needs every reaction assigned to a propensity group)
- T. E. Turner, S. Schnell and K. Burrage. Stochastic approaches for modelling in vivo reactions. Computational Biology and Chemistry 28(3):165-178, 2004. https://doi.org/10.1016/j.compbiolchem.2004.05.001 - the book's ref. [870], the review of stochastic simulation methods cited in 18.3.3
- R. Schwartz. Biological Modeling and Simulation. MIT Press, 2008, sec. 17.4 - the book's ref. [757], the overview of methods for very large reaction networks cited in 18.3.3

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- ATTRIBUTION: the algorithm is the book's own proposal, given as three lines of prose in 18.3.3 with no pseudo-code, no parameters and no published numbers. Its refs [757], [780] and [870] are cited for *other* algorithms (a textbook overview, Slepoy's constant-time composition-rejection KMC, and a review of stochastic methods) and none of them contains this algorithm, so nothing here is reconstructed from them. The gate itself is eq. 2.32 (Arrhenius, cited by the book to Atkins [42]) and the collision loop is the chapter-2 algorithm of sec. 2.3.3.
- DISTRIBUTION OF THE COLLISION ENERGY, the one real gap: the book says to react when the kinetic energy exceeds Ea but never says how that energy is distributed, and the answer decides whether the algorithm reproduces eq. 2.32 at all. The bath model uses the standard collision-theory result (Atkins): the energy along the line of centres of a Maxwell-Boltzmann gas is Exponential with mean RT, whose tail above Ea is exactly exp(-Ea/RT). No other choice gives the Arrhenius factor.
- The book's literal wording in 2.2.6 is 'the sum of their kinetic energies' of the two molecules. Two independent Exponential(RT) energies sum to a Gamma(2, RT), whose tail is (1 + Ea/RT) exp(-Ea/RT), NOT the Arrhenius factor - it overestimates the rate by a factor 1 + Ea/RT (3.4x at the default Ea, 17x at Ea = 40 kJ/mol). Both readings are implemented and measured (energy_model), and the discrepancy is asserted in the tests rather than hidden.
- ENERGY BOOKKEEPING: the book does not say whether energy is conserved. bath is a thermostat (T fixed, kinetic energy untracked, the reaction heat booked against the bath in extras.analysis.energy_budget). conserved is a closed adiabatic vessel: each molecule carries a kinetic energy initialised i.i.d. Exponential(RT), the pair's pooled energy minus delta_G goes to the products, and E = Ek + Ep (eq. 2.28) is conserved to floating-point exactness. Being adiabatic, its temperature drifts as reaction heat accumulates, so only bath reproduces the Arrhenius factor at stationarity.
- Products share the leftover energy as uniform spacings, Dirichlet(1, ..., 1). For the two-product case this Beta(1,1) split is exactly the one that maps a Gamma(2, RT) pair energy back onto two independent Exponential(RT) energies, so a thermoneutral conserved run keeps its Maxwell-Boltzmann pool instead of drifting for a reason that is only an artefact of the splitting rule.
- In the conserved model the kinetic energies are a well-stirred pool drawn independently of which molecules the soup picked, rather than being attached to individual molecules. In a well-stirred vessel the two are statistically identical, and it keeps the observed network keyed by species instead of by (species, energy) pairs.
- The default system is the smallest one that exercises the gate at two different barriers: X1 + X2 <-> Y1 + Y2 with a forward barrier Ea and, from fig. 2.4, a reverse barrier Ea - delta_G. That relation makes it thermodynamically consistent by construction (kf/kr = exp(-delta_G/RT), eq. 2.33), and user-supplied systems that give both directions of a reaction are checked against it.
- Ea >= max(0, delta_G) is enforced: a barrier below the product well would put the transition state of fig. 2.4 under the products. It also guarantees the products can always be given a non-negative kinetic energy in the conserved model.
- Energies are in kJ/mol and R = 8.31451e-3 kJ/K/mol (the book states R = 8.31451 J/K/mol). Each rate dict carries its own T, R and units next to the required A and Ea so the exponent cannot be misread.
- Defaults are small (M = 400, 4000 collisions = 10 generations) and the barriers are a few kJ/mol: real barriers of 40-400 kJ/mol would reject all but ~1e-7 of the collisions at 300 K and measure nothing in a short run. Paper-scale values are in the params' range notes.
- One collision tests one reaction: when several rules share a reactant pair, one is drawn uniformly and gated, so a channel's measured acceptance is its own Arrhenius factor and not a competition between channels. The default system gives each pair a single channel.
- The v1 parameters Ea_fn and kinetic_energy_fn were callables, which schema v2 forbids; they are replaced by the explicit per-reaction Ea of the system, the energy_model choice and the temperature.

## Notes

Not a chemistry but the REACTOR ALGORITHM of section 18.3.3, catalogued as an analysis: `generate` runs it and reports measured against predicted kinetics in extras.analysis, so the entry is its own validation. It is a rejection method - the same trade as Slepoy's composition-rejection KMC [780], but rejecting on physics (a barrier) rather than on a propensity bound, and needing no enumeration of the reaction set. Worth shipping as a Chemart reactor backend alongside SSA / tau-leaping / next-reaction, with the elastic fraction reported so the user can see when it stops paying off.

---

*Specification: `catalog/chemistries/energy-gated-collision.yaml` · generator: `chemart/chemistries/energy_gated_collision.py` · tests: `tests/chemistries/test_energy_gated_collision.py`*
