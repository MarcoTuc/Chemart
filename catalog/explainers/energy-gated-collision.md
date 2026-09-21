## Introduction

This entry is not a chemistry. It is a *reactor algorithm*: a recipe for
simulating which reactions happen in a pot of molecules, and how often. Most
of the catalog describes what the molecules are and how they react; this entry
is about the machinery that runs such a model. It is catalogued as an
*analysis* because a run is its own test: it measures the reaction rates the
algorithm produces and compares them with the textbook formula they should
obey.

The idea is Banzhaf and Yamamoto's own, given in three sentences of book
§18.3.3 ("Algorithms for Simulating Large-Scale Reaction Networks"). The
problem they raise is cost. Gillespie's stochastic simulation algorithm (SSA),
the standard way to simulate reactions one event at a time, computes a
*propensity* (the current probability per unit time) for every possible
reaction at every step, so each step costs more as the number of reactions
grows. Faster variants such as the next reaction method, tau-leaping, or the
constant-time method of Slepoy, Thompson and Plimpton (2008) still need the
complete list of possible reactions up front. The book asks: "Is it possible to
have a stochastic reaction algorithm that takes into account the different
kinetic rates of different reactions, without having to calculate the
propensities of all reactions that might possibly occur?"

Its answer is a small change to the simplest reactor in the book, the naive
collision algorithm of §2.3.3. That algorithm picks two molecules at random
and, if a rule applies to the pair, replaces them by the products. The change
is an energy gate. Real reactions have an *activation energy* `Ea`, a barrier
the colliding molecules must climb before they can turn into products. So: pick
two molecules; if their collision energy clears the barrier, react (an
*effective* collision); otherwise put them back unchanged (an *elastic*
collision). A reaction with a high barrier then fires on only a small fraction
of the collisions of its reactants, a low-barrier one on most of them, and
different reactions run at different speeds without any rate ever being summed.
In a thermal gas that fraction is `exp(-Ea/RT)`, the factor in the Arrhenius
equation (book eq. 2.32), so the gate turns barriers into Arrhenius rates.

The book states the cost of the idea in the same paragraph: the algorithm is
effective only if elastic collisions are much rarer than effective ones, which
happens only when almost every molecule can react with almost every other, and
such networks "would most probably be of little practical use". It then moves
on to rule-based modelling, where reactions are computed on the fly from
general rules. There is no paper behind this algorithm, no published run and no
reference implementation (PyCellChemistry, the book's companion code, has no
energy-gated reactor). Chemart's version fills the gaps the three sentences
leave, and tests the result against the equations of the book's chapter 2.

Nearest neighbours: [the chameleon chemistry](chameleon.md) is the book's
worked example of the naive collision algorithm with elastic collisions (book
§2.5.1); [reversible dimerization](dimerization.md) is the catalog's other
fixture for checking that a simulation reaches the right equilibrium;
[SRSim](srsim.md) is the rule-based simulator described in the same book
section. [ToyChem](toychem.md) also attaches an Arrhenius rate
`A exp(-Ea/RT)` to each reaction, with barriers computed from its molecules'
orbital energies; here the barriers are given, and the question is how a
simulation turns them into rates.

## How it works

### Energies, barriers and temperature

Every molecule has a *potential energy* `Ep`, stored in its bonds, and a
*kinetic energy* `Ek`, from its motion; the book writes the total as
`E = Ek + Ep` (eq. 2.28). A reaction goes from a reactant state to a product
state, and the difference in energy between them is `delta_G` (ΔG, the change
in Gibbs free energy). A negative ΔG means the products lie lower: the reaction
releases energy and is called *exergonic*, or spontaneous. Between the two
states sits the barrier `Ea`.

The book's figure 2.4 draws this as a hill between two valleys. Read from the
other side, the same hill is the barrier of the reverse reaction, and it is
higher or lower by exactly ΔG: `Ea(reverse) = Ea(forward) - ΔG`. The book
derives from this that the equilibrium constant, the ratio of products to
reactants once forward and reverse rates balance, is `K = exp(-ΔG/RT)` (eq.
2.33).

`T` is the absolute temperature in kelvin and `R` the gas constant,
8.31451 J/(K·mol). Their product `RT` is the typical thermal energy per mole:
2.494 kJ/mol at 300 K. Every energy in this entry is in kJ/mol.

### The gate

One step of the algorithm is one collision:

1. Draw two different molecules at random from the well-stirred vessel.
2. Look up the rule for that unordered pair. No rule: the collision is elastic.
3. Draw a collision energy and compare it with the rule's `Ea`. Below the
   barrier: elastic. At or above: replace the two reactants by the products.

The book does not say how the collision energy is distributed, and the answer
decides whether the gate reproduces Arrhenius at all. Chemart's default,
`energy_model="bath"`, uses the standard result of collision theory (Atkins
and de Paula, the book's own reference for eq. 2.32): in a gas at temperature
`T` the energy of a collision along the line joining the two molecules'
centres is exponentially distributed with mean `RT`, so the chance it exceeds
`Ea` is exactly `exp(-Ea/RT)`. The vessel is held at a fixed temperature by an
imagined heat bath, which also absorbs the heat that reactions release.

The acceptance probability depends on `Ea` and `T` only through their ratio
`Ea/RT`. A barrier of 2 kJ/mol at 300 K is 0.80 RT and lets 45% of collisions
through; 40 kJ/mol, a realistic chemical barrier, is 16 RT and lets about one
in ten million through.

### A worked example

The default system is one reversible reaction, `X1 + X2 <-> Y1 + Y2`, with a
forward barrier of 2 kJ/mol and ΔG = −3 kJ/mol. The product well sits 3 kJ/mol
lower (Chemart puts that energy on `Y1`), so the reverse barrier is
2 + 3 = 5 kJ/mol. The vessel starts with 200 `X1` and 200 `X2` and runs 4,000
collisions. The two reactions printed below are the ones that fired, with how
often:

```
X2 + X1 -> Y1 + Y2  [arrhenius A=1.0 Ea=2.0 T=300.0 R=0.00831451 units=kJ/mol]  (x212)
Y1 + Y2 -> X1 + X2  [arrhenius A=1.0 Ea=5.0 T=300.0 R=0.00831451 units=kJ/mol]  (x81)
```

Of the 4,000 draws, 473 happened to pick an `X1` and an `X2`. The gate let 212
of them through: a measured acceptance of 0.448, against the predicted
`exp(-2/2.494) = 0.4485`. Another 629 draws picked a `Y1` and a `Y2`; 81 passed
the 5 kJ/mol barrier, 0.129 against a predicted 0.135. The other 2,898 draws
picked a pair with no rule, such as `X1 + X1` or `X1 + Y2`, and changed
nothing. So 92.7% of the collisions were elastic, which is the book's caveat
seen from the inside: most of the work of the run is spent on collisions that
do nothing. The vessel ended with 69 of each `X` and 131 of each `Y`.

Time here is counted in collisions. The book's unit is the *generation*, `M`
collisions in a vessel of `M` molecules (§2.6.1), so this run lasted 10
generations.

### The book's literal wording: `energy_model="conserved"`

Book §2.2.6 says two molecules react if "the sum of their kinetic energies" is
high enough. Chemart implements that reading too. Each molecule carries a
kinetic energy, drawn at the start from the same exponential distribution with
mean `RT`; the gate compares `Ea` with the sum of the two colliding molecules'
energies; and after an effective collision the products share the pair's
energy minus ΔG, so the total `Ek + Ep` of the closed vessel is conserved
exactly. The sum of two such energies is not exponential, and its chance of
exceeding `Ea` is `(1 + Ea/RT) exp(-Ea/RT)`: larger than the Arrhenius factor by
the factor `1 + Ea/RT`, which is 17 for a 40 kJ/mol barrier at 300 K. This
reading therefore does *not* reproduce eq. 2.32. Because no heat bath is
present, an exergonic run also heats itself up, and the acceptance drifts
further upward as it goes.

### In short

The molecules are whatever species a system lists, each with a potential
energy. The reactions are bimolecular rules with a barrier each, applied only
when the gate opens. The reactor is the closed, well-stirred vessel of §2.3.3:
two molecules in, two out, so the population never changes size.

## Using it

The default run is the worked example above. It reproduces no published
experiment, because there is none; its purpose is to put measured and predicted
rates side by side. They are in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["effective_collisions"], a["elastic_fraction"]     # (293, 0.92675)
a["final_population"]                                # {'X1': 69, 'X2': 69, 'Y1': 131, 'Y2': 131}
m = a["reactions"]["X1+X2->Y1+Y2"]
m["attempts"], m["effective"], m["acceptance"]       # (473, 212, 0.448...)
m["arrhenius_factor"]                                # 0.4485...  exp(-Ea/RT)
m["k_observed"], m["k_arrhenius"]                    # A x acceptance, A x exp(-Ea/RT)
```

For each reaction, `attempts` counts the collisions of its reactants,
`effective` those that passed the gate, and `predicted_acceptance` is what the
chosen energy model should give (the Arrhenius factor for `bath`,
`(1 + Ea/RT) exp(-Ea/RT)` for `conserved`). `energy_budget` books the energy:
`supplied_by_bath` in the `bath` model, and the kinetic, potential and total
energies with their `drift` in the `conserved` model.
`net.extras["energies"]` lists the potential energies, barriers and ΔG of
every reaction.

The prefactor `A` is carried on each rate and multiplies `k_observed`, but it
does not enter the simulation: changing it leaves the run identical. In this
algorithm every reaction's underlying collision frequency is set by how often
its reactants meet, so two reactions cannot be given different prefactors.

**An Arrhenius plot.** Measure the rate at five temperatures and fit `ln k`
against `1/T`; the slope is `−Ea/R` and the intercept `ln A`. The runs below
take about 5 seconds together:

```python
import numpy as np, chemart
R = 8.31451e-3
xs, ys = [], []
for T in (250, 300, 350, 420, 500):
    net = chemart.generate_network("energy-gated-collision", seed=1, Ea=4.0, A=2.5,
                                   temperature=T, steps=40000)
    m = net.extras["analysis"]["reactions"]["X1+X2->Y1+Y2"]
    xs.append(1 / T); ys.append(np.log(m["k_observed"]))
slope, intercept = np.polyfit(xs, ys, 1)
print(-slope * R, np.exp(intercept))     # 4.038  2.529   (put in: 4.0, 2.5)
```

Fewer collisions give noisier estimates: with 12,000 per temperature the
same seed recovers `Ea` = 4.5 and `A` = 2.9.

**Comparing the two energy models.** On a thermoneutral system (`delta_G=0.0`,
so the vessel neither heats nor cools) with `Ea=4.0` and `steps=8000`, the
`bath` model accepts 0.200 of the collisions against an Arrhenius factor of
0.201; `energy_model="conserved"` accepts 0.530, matching its predicted
0.524 and far from Arrhenius. On the default exergonic system the conserved
model keeps the total energy to within 1e-13 kJ/mol (974.28 kJ/mol before and
after), while the kinetic energy rises from 974 to 1,307 kJ/mol as reactions
release heat; the reverse reaction's acceptance, 0.51 against a predicted 0.40
at the starting temperature, shows the vessel warming up.

**Equilibrium.** Run long enough for forward and reverse to balance and read
`K = (Y1·Y2)/(X1·X2)` from `final_population`. Over seeds 0–5 with
`steps=20000`, the default (ΔG = −3) gives a mean K of 3.66 against
`exp(3/2.494) = 3.33`, and `Ea=3.0, delta_G=2.0` gives 0.40 against 0.45. The
twelve runs take about 5 seconds.

**Your own system.** Pass `system` with potential energies, bimolecular
reactions with barriers, and starting proportions. Two rules may share a
reactant pair; each collision of that pair then picks one rule at random and
tests only its barrier. Each rule is thus tried on its share of the pair's
collisions (here about half), and its measured acceptance is its own
Arrhenius factor:

```python
system = {
    "energies": {"A": 0.0, "B": 0.0, "C": -5.0, "D": 1.0, "E": 0.0},
    "reactions": [
        {"id": "fast", "reactants": ["A", "B"], "products": ["C", "E"], "Ea": 1.0},
        {"id": "slow", "reactants": ["A", "B"], "products": ["D", "E"], "Ea": 6.0},
    ],
    "initial": {"A": 1, "B": 1},
}
net = chemart.generate_network("energy-gated-collision", seed=1, system=system, steps=2000)
# fast: 170 attempts, 118 effective (0.694; Arrhenius 0.670)
# slow: 170 attempts,  10 effective (0.059; Arrhenius 0.090)
```

The generator refuses a reaction with other than two reactants, a species
without an energy, a barrier below ΔG (the top of the hill cannot lie under the
product valley), and a forward/reverse pair whose barriers do not differ by
exactly ΔG.

**Realistic barriers, and cost.** With `Ea=40.0` and 100,000 collisions (about
1 second), 49,755 collisions of `X1` with `X2` all bounce off: at an
acceptance of 1.1e-7 nothing reacts. Keep `Ea/RT` below about 5 for runs of
this size. A single run costs roughly 14 µs per collision, so a million
collisions take about 14 seconds.

## Results

There are no published results for this algorithm. The book proposes it in
three sentences, with no pseudo-code, parameters or runs, and judges it
practical only for densely connected networks that it expects to be of
little use. The works the book cites around it ([757], [780], [870]) are
about other algorithms, and the implementation notes record that none of them
contains this one. What can be checked is whether it delivers what the book says it
would: the kinetics of chapter 2. Chemart's tests check each of the claims
below; they run the default reversible system at 300 K unless stated.

**The acceptance is the Arrhenius factor.** The fraction of a reactant pair's
collisions that are effective should be `exp(-Ea/RT)`. The tests measure it at
six combinations of barrier and temperature (1 to 8 kJ/mol, 200 to 500 K),
for both the forward and the reverse reaction of each run, over 8,000
collisions, and require every measurement to lie within four binomial
standard errors of the prediction. Reproduced, in the `bath` model.

**The Arrhenius plot is straight.** Plotting `ln k` against `1/T` is how
chemists read a barrier off measured rates, and it should give back the `Ea`
and `A` that went in. The test fits five temperatures from 250 to 500 K and
requires `Ea` within 15%, `A` within 20%, and no point off the line by more
than 0.12 in `ln k`. Reproduced.

**Barrier and temperature act through `Ea/RT`, exponentially.** The test checks
that the acceptance falls steadily as `Ea` rises from 0 to 8 kJ/mol and rises
steadily with temperature from 150 to 700 K, and that the extremes differ more
than tenfold. Reproduced.

**No barrier gives back the naive algorithm and mass action.** At `Ea = 0`
every collision of a reactive pair is effective, and the only elastic
collisions are pairs with no rule: this is the §2.3.3 algorithm. Its statistics
are those of the law of mass action: for a one-way reaction
`X1 + X2 → Y1 + Y2` starting with `n0` of each reactant in a vessel of `M`
molecules, the count after `s` collisions follows
`1/n = 1/n0 + 2 P s / (M (M − 1))`, where `P` is the acceptance probability.
The tests check the zero-barrier case exactly, and invert the decay law to
recover `P` at barriers of 0, 2 and 5 kJ/mol, averaged over five runs, within
25%. The tolerance is loose because the law is a mean-field approximation
that is off by a few percent at 400 molecules. Reproduced within that
tolerance.

**Barriers that differ by ΔG give the right equilibrium.** Forward and reverse
barriers set as in figure 2.4 should drive the vessel to `K = exp(-ΔG/RT)`. The
test averages six long runs (20,000 collisions) at ΔG = −3 and +2 kJ/mol and
requires the mean within 25%. Reproduced within that tolerance; this test is
marked slow.

**Energy is conserved in the `conserved` model, and that model is not
Arrhenius.** The tests check that total energy `Ek + Ep` is conserved to
floating-point precision for ΔG of 0, −3 and +2 kJ/mol, that an exergonic run
heats up, and that the acceptance follows `(1 + Ea/RT) exp(-Ea/RT)` rather than
the Arrhenius factor. So the book's literal wording ("the sum of their kinetic
energies") overestimates the rate by the factor `1 + Ea/RT`; only the
line-of-centres energy of collision theory gives eq. 2.32.

**Cost does not depend on the number of reactions, but elastic collisions
dominate.** Each collision looks up only the rules for the drawn pair, so its
cost does not grow with the size of the reaction set, while the fraction of
useful collisions is at most `exp(-Ea/RT)` for each reactive pair. The tests
check only the second half: the default run must have more than 50% elastic
collisions (it has 92.7%). The claim about cost is true of the code but is not
measured, and Chemart does not compare this algorithm's speed with SSA or any
other method.

**What Chemart does not provide.** The algorithm runs only inside this entry,
on the systems passed to it; it is not available as a reactor for the other
chemistries in the catalog. Every reaction must be bimolecular, and all
reactions share one prefactor.

## Further reading

- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical
  reactions. *Journal of Physical Chemistry* 81, 2340–2361. The SSA that
  §18.3.3 sets out to avoid; the book treats it in chapter 4.
