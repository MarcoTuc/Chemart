## Introduction

Okamoto's biochemical switch is a small network of enzyme reactions meant to
behave like an electronic switch: two input chemicals compete, and a pair of
output chemicals reads "on" or "off" according to which input is winning. The
outputs are two forms of an enzyme *cofactor* (a helper molecule an enzyme
needs in order to work), called `A` and `B`. Almost all of the cofactor sits in
one form or the other, so the pair carries one bit: `A` high and `B` low, or
the reverse. When the balance between the inputs turns, the bit flips, and it
flips sharply.

It comes from Masahiro Okamoto, K. Hayashi and colleagues in Fukuoka, Japan
(Kyushu University, then the Kyushu Institute of Technology), who worked on it
from 1980 to 1990. Building on Robert Rosen's
"two-factor" models of neurons and biochemical automata, they first proposed an
enzymatic model coupling an excitatory and an inhibitory factor (1980). They
then studied its switching in detail as a *cyclic enzyme system*, two enzymes
that share cofactors in a cycle: as a "chemical diode" (1987), as a device
realising the on-off McCulloch-Pitts neuron equation (1988), with a way to set
the switching time by a pulse (1989), and, connected in series, as a prototype
"artificial neuronic device" (1990). Their stated aim, in the 1987 abstract,
was a "switching circuit (biochip) in a bio-computer".

The book uses the switch as its first example of *wet chemical computing
simulated in silico*: modelling on a computer a chemical computer that would
be too costly to build in the laboratory (book §17.4.1). The section places it
in a line of early ideas: Seelig and Rössler's chemical flip-flops built on
bistable reactions (early 1970s), Gánti's "fluid machines" (see the
[chemoton](chemoton.md)), and later the chemical neurons, logic gates and
finite state machines of Hjelmfelt, Weinberger and Ross, which were based on a
modified version of Okamoto's neurons and led to the result that chemical
kinetics can in principle compute anything a Turing machine can.

What Chemart holds is not the original enzyme model but the book's six-reaction
sketch of it, with no enzymes written out: a fixed, non-constructive kinetic
model to be integrated as differential equations. Its neighbours in the same
chapter are the [Brusselator](brusselator.md), an oscillator (a clock for such
circuits, §17.4.2), and [analog computation with reaction
networks](analog-function-crn.md), which reads numbers off steady-state
concentrations (§17.4.3). The switch is the catalog's only entry built to store
and flip a single bit.

## How it works

### The species

- `I1` and `I2` are the two **inputs**. The book treats them as controlled
  from outside: they feed the network but are not used up by it.
- `X1` is made from `I1`, and `X3` from `I2`. These are **substrates**, the
  molecules the cofactors work on.
- `A` and `B` are the **cofactors**, the outputs. There is one unit of
  cofactor in total (`A + B = 1` throughout).
- `X2` and `X4` are **products**, which decay away.

### The reactions

The default network, as generated above, with the name of each rate
constant (the parameter that sets it):

```
I1 -> I1 + X1      k_in   = 1       input I1 supplies substrate X1
I2 -> I2 + X3      k_in   = 1       input I2 supplies substrate X3
A + X3 -> B + X4   k2     = 50000   A converts X3, and becomes B
B + X1 -> A + X2   k1     = 50000   B converts X1, and becomes A again
X2 -> ∅            k3     = 10      product X2 decays
X4 -> ∅            k4     = 10      product X4 decays
I1 -> I2           k_conv = 0.006   the inputs drift across each other
```

Every reaction runs at *mass action*: its speed is its rate constant `k`
times the product of its reactants' concentrations, so `A + X3 -> B + X4` runs
at `k2 × [A] × [X3]` (square brackets mean concentration). The two cofactor
reactions are very fast (`k1 = k2 = 5 × 10⁴`, the values the book takes from
Okamoto et al. 1987); everything else is slow.

The first two reactions are written `I1 -> I1 + X1` so that the inputs are not
consumed, as the book asks. The last reaction is not part of the switch. It is
how the book produces a test signal: `I1` starts at 100 and `I2` at 80, and
`I1` slowly turns into `I2`, so `I1` falls, `I2` rises, and the two cross when
both reach 90. With the default rate this happens at about 17.5 s, as in the
book's figure 17.5 (top).

### Why it holds a state, and when it flips

Start where the book does, with all cofactor as `A` (`A = 1`, `B = 0`). Each
molecule of `X3` that `I2` supplies is caught at once by `A` and converted,
which turns one unit of `A` into `B`. That `B` meets `X1`, of which there is
plenty, and turns straight back into `A`. So the cofactor goes round the cycle
`A → B → A`, each turn using one `X3` and one `X1`, and `A` stays at 1 because
`B` never lasts long enough to accumulate.

The cycle can only turn as fast as the scarcer substrate arrives. While `I1`
is larger than `I2`, `X1` arrives faster than `X3`, so `X1` piles up: the
network stores the surplus `∫ k_in (I1 − I2) dt`. After the inputs cross, `X3`
arrives faster than `X1` and the stock of `X1` is drawn down. As soon as it is
empty, the `B` made from `X3` has nothing left to react with, and within a
fraction of a second all the cofactor is `B`. Now the roles are exchanged:
`B` is saturated, `X1` is used up as it arrives, and surplus `X3` piles up.

Running the default network shows exactly this (the script is under *Using
it*):

```
t=  0.0  I1=100.00 I2= 80.00  A=1.000 B=0.000  X1=   0.0 X3=   0.0 X2=8.00
t= 10.0  I1= 94.18 I2= 85.82  A=1.000 B=0.000  X1= 141.2 X3=   0.0 X2=8.58
t= 17.5  I1= 90.03 I2= 89.97  A=1.000 B=0.000  X1= 172.5 X3=   0.0 X2=8.99
t= 30.0  I1= 83.53 I2= 96.47  A=1.000 B=0.000  X1=  91.0 X3=   0.0 X2=9.64
t= 35.0  I1= 81.06 I2= 98.94  A=1.000 B=0.000  X1=  13.9 X3=   0.0 X2=9.89
t= 36.0  I1= 80.57 I2= 99.43  A=0.000 B=1.000  X1=   0.0 X3=   3.5 X2=8.23
t= 40.0  I1= 78.66 I2=101.34  A=0.000 B=1.000  X1=   0.0 X3=  86.6 X2=7.87
t=120.0  I1= 48.68 I2=131.32  A=0.000 B=1.000  X1=   0.0 X3=4490.7 X2=4.87
```

`X1` peaks at the crossing (about 172 at 17.5 s), is used up by 36 s, and
then `A` and `B` swap. The swap itself is sharp: `A` falls from 0.9 to 0.1 in
about 0.04 s. But it comes about 18 s *after* the inputs crossed, when the
surplus stored before the crossing has been paid back, not a few seconds after
as in the book's figure (see *Results*). The product `X2` tracks the speed of
the cycle, staying close to `k_in × min(I1, I2) / k3`: 8 at the start, when `I2 = 80`, and
`I1 / 10` after the flip. `X4` is made and removed at the same rates and starts
at the same value, so it is always equal to `X2`.

The *reactor* is a set of ordinary differential equations (ODEs), one per
species, that give the rate of change of each concentration. Because the two
cofactor reactions are thousands of times faster than the rest, the equations
are *stiff*, and need an implicit solver (such as SciPy's `Radau`) to be
integrated efficiently. Nothing in the network grows new species. The only
removal is the decay of `X2` and `X4`; the stored substrate has no way out, so
the substrate of the larger input keeps growing for as long as the inputs
differ.

## Using it

The default call above gives the book's figure 17.5 (top): rate constants and
initial concentrations as the book gives them for its simulation
(`X1 = X3 = 0`, `X2 = X4 = 8`, `A = 1`, `B = 0`), inputs starting at 100 and
80, and the conversion rate that makes them cross at about 17.5 s. `k_in`, the
rate at which inputs make substrate, is not in the book; its default of 1 is
Chemart's choice. The network is deterministic, so `seed` changes nothing, and
`net.extras` is empty. `net.initial_state` holds the starting concentrations.

`chemart.simulate.ode` integrates the network under mass action. This script
printed the table in *How it works*, and takes a few seconds:

```python
import numpy as np

import chemart
from chemart import simulate


def run(net, t_end, n=20001):
    # the network is stiff (k1 = k2 = 5e4), so use an implicit solver
    traj = simulate.ode(net, t_end, points=n, solver="Radau")
    ids, t, X = traj.array([s.id for s in net.species])
    return t, dict(zip(ids, X.T))


t, x = run(chemart.generate_network("okamoto-switch"), 120)
for s in (0, 10, 17.5, 30, 35, 36, 40, 120):
    i = np.argmin(abs(t - s))
    print(f"t={s:5.1f}  I1={x['I1'][i]:6.2f} I2={x['I2'][i]:6.2f}  A={x['A'][i]:.3f} "
          f"B={x['B'][i]:.3f}  X1={x['X1'][i]:6.1f} X3={x['X3'][i]:6.1f} X2={x['X2'][i]:.2f}")
```

#### Delaying the switch (figure 17.5, bottom)

The book's second experiment slows the conversion `I1 -> I2`, so the inputs
cross later and the switch flips later. The book does not give the slower
rate; its close-up shows the flip between 68 and 70 s, and a `k_conv` of
0.0015 (a quarter of the default) puts the crossing at about 70 s.
With the `run` above, this runs three conversion rates and compares the
flip with the moment the stored surplus of `X1` is used up. With `I1(0) = 100`
and `I1 + I2 = 180`, that moment is the time `T` at which
`200 (1 − e^(−k_conv T)) / k_conv = 180 T`:

```python
from scipy.optimize import brentq

for k_conv in (0.006, 0.003, 0.0015):
    t, x = run(chemart.generate_network("okamoto-switch", k_conv=k_conv), 200)
    a = x["A"]
    cross = t[np.argmax(x["I1"] <= x["I2"])]
    flip = t[np.argmax(a < 0.5)]
    width = t[np.argmax(a < 0.1)] - t[np.argmax(a < 0.9)]
    # when the stored surplus of X1, the integral of k_in (I1 - I2) dt, is used up
    empty = brentq(lambda T: 200 * (1 - np.exp(-k_conv * T)) / k_conv - 180 * T, 1, 1e4)
    print(f"k_conv={k_conv}: inputs cross {cross:.1f} s, A flips {flip:.1f} s "
          f"(0.9 -> 0.1 in {width:.2f} s), X1 peaks at {x['X1'].max():.0f}, "
          f"surplus used up at {empty:.1f} s")
```

```
k_conv=0.006: inputs cross 17.6 s, A flips 35.8 s (0.9 -> 0.1 in 0.04 s), X1 peaks at 173, surplus used up at 35.8 s
k_conv=0.003: inputs cross 35.1 s, A flips 71.5 s (0.9 -> 0.1 in 0.04 s), X1 peaks at 345, surplus used up at 71.5 s
k_conv=0.0015: inputs cross 70.2 s, A flips 143.1 s (0.9 -> 0.1 in 0.04 s), X1 peaks at 690, surplus used up at 143.0 s
```

Slowing the conversion delays the flip, as in the book. But in every case the
flip falls where the surplus runs out, at about twice the crossing time, so a
slower drift means a longer lag after the crossing (73 s at `k_conv = 0.0015`),
where the book's figure shows the flip close to the crossing. This run
takes about five seconds.

#### Held inputs, and the input rate

Setting `k_conv = 0` holds the inputs fixed. The cofactor then settles on the
side of the larger input and stays there, while the substrate of the larger
input grows without limit, at `k_in × (I1 − I2)`, 20 per second here:

```python
for I1, I2 in ((100, 80), (80, 100)):
    t, x = run(chemart.generate_network("okamoto-switch", k_conv=0, I1_0=I1, I2_0=I2), 60)
    print(f"I1={I1}, I2={I2}: A(60)={x['A'][-1]:.3f}  X1 at 30 s, 60 s: "
          f"{x['X1'][10000]:.0f}, {x['X1'][-1]:.0f}  X3 at 30 s, 60 s: "
          f"{x['X3'][10000]:.0f}, {x['X3'][-1]:.0f}")

for k_in in (0.1, 1.0, 10.0):
    t, x = run(chemart.generate_network("okamoto-switch", k_in=k_in), 60)
    print(f"k_in={k_in}: A flips at {t[np.argmax(x['A'] < 0.5)]:.1f} s")
```

```
I1=100, I2=80: A(60)=1.000  X1 at 30 s, 60 s: 600, 1200  X3 at 30 s, 60 s: 0, 0
I1=80, I2=100: A(60)=0.000  X1 at 30 s, 60 s: 0, 0  X3 at 30 s, 60 s: 599, 1199
k_in=0.1: A flips at 36.0 s
k_in=1.0: A flips at 35.8 s
k_in=10.0: A flips at 35.8 s
```

With `I2` the larger, `A` starting at 1 flips to `B` at once (within 0.03 s
in a separate run with finer output), since there is no stock of `X1` to
hold it. The second loop shows that `k_in`
barely moves the flip time (36.0 s at 0.1, 35.8 s at 1 and 10): it scales the
surplus and its repayment alike. So in this model the switching time is set by
the input signals alone. The parameter table below lists all rates and
initial inputs.

## Results

#### What the book shows

The book's figure 17.5 is an ODE simulation of the six reactions with the
parameters of Okamoto et al. (1987). In the top panel the inputs cross, and "a
few seconds after the two inputs reach equal concentration, the concentration
of A flips from one to zero, while B flips in the reverse direction"; the
close-up of `A` and `B` spans 17 to 19 s. In the bottom panel the conversion
of `I1` into `I2` is slower, the inputs cross later, and the flip moves later
with them (close-up 68 to 70 s). The point is that the switching time can be
set by acting on the inputs. The book calls this "a very rough sketch" of the
switch and refers to the original paper for the detailed analysis.

**In Chemart.** The test `test_okamoto_switch_flips_after_inputs_cross`
(in `tests/chemistries/test_w1_dynamics.py`) integrates the default network and
checks that the inputs cross at 17.5 ± 0.5 s, that `A` is above 0.99 and `B`
below 0.01 at 10 s, that `A` flips to `B` after the crossing and stays there
until 120 s, and that halving `k_conv` delays both the crossing and the flip.
These hold. What Chemart does not reproduce is the timing: the flip comes at
35.8 s, 18 s after the crossing, not a few seconds after it, and the lag grows
as the conversion slows. The reason, recorded in the implementation decisions,
is that with the reactions as printed the substrate of the dominant input is
never lost while its cofactor is saturated, so the surplus accumulates and must
be used up before the switch can flip. The original model has reversible
enzyme steps and substrate inflow and outflow that the book leaves out (the
1988 paper on the monocyclic enzyme system writes each step as a reversible
reaction with its enzyme, `X1 + B + E1 ⇌ A + E1 + X2` and
`X3 + A + E2 ⇌ B + E2 + X4`, with arrows in and out). The original paper was
not accessible, so Chemart keeps the book's scheme.

#### Bistability

The book introduces the switch through bistability: a system is bistable when
it has two stable steady states with an unstable one between them, so it
settles into one or the other. The model shows two distinct output states and
a sharp jump between them. But in the six reactions as printed there is no
steady state at all while the inputs differ, because the surplus substrate
grows without limit (the held-input run above). The state the switch is in
depends on the history of the inputs, through the stored surplus, rather than
on two stable steady states. Whether the original enzyme model has true
bistable steady states cannot be checked from the sources Chemart used.

#### What the original papers report

From the abstracts of the Okamoto papers (the full texts were not read):

- **1987, "chemical diode".** The authors related the switching mode to the
  phase difference of two sinusoidal inputs, predicted the switching time
  theoretically, and observed half-wave and full-wave rectification when the
  two inputs are in antiphase (phase difference π). The book adds that the
  paper shows the switch flipping on and off in response to a sinusoidal input.
  Chemart's generator has no oscillating inputs, so this is not reproduced.
- **1988, McCulloch-Pitts equation.** Cyclic enzyme systems have "the
  reliability of ON-OFF types of operation" (the McCulloch-Pitts neuron
  equation); the paper compares them with electronic switching circuits,
  especially on their memory ("mnemonic") mechanism. A companion 1988 paper
  found that the system works as a switch only when the initial
  concentrations of enzymes or cofactors exceed a threshold.
- **1989, turning the switch on and off.** The switching time "was inevitably
  determined in accordance with the difference in amount between two inputs",
  which the authors saw as a drawback; a pulse perturbation let them set it
  independently of the inputs. The first half matches what the Chemart runs
  show, where the flip comes when the integrated input difference returns to
  zero.
- **1990, networks.** Switches connected in series into a prototype
  "artificial neuronic device" showed, among other things, that the number of
  excited elements grows with the strength of an excitatory stimulus, that the
  signal is amplified up to a limit and attenuated as it travels, and something
  like long-term potentiation when several stimuli are repeated.

None of these experiments is reproduced by Chemart.

#### What was built on it

The book traces a line from Okamoto's switches to general chemical computing.
Hjelmfelt, Weinberger and Ross (1991) built chemical neurons from a modified
version of Okamoto's neurons, using only reversible reactions controlled by an
enzyme, connected them into neural networks and made logic gates (AND, OR,
NOT, NAND). In 1992 they built finite state machines from similar neurons,
paced by a chemical clock, with a binary decoder, an adder and a stack memory,
and argued informally that a universal Turing machine could be built this way.
Magnasco (1997) extended the argument, and the book concludes that chemical
computation can be Turing universal. None of these constructions is in
Chemart: this entry is the single switch only.

## Further reading

- Okamoto, M., Sakai, T. & Hayashi, K. (1988). Biochemical switching device:
  monocyclic enzyme system. *Biotechnology and Bioengineering* 32(4), 527–537.
  doi:10.1002/bit.260320416. The fuller cyclic enzyme model, with enzymes and
  reversible steps.
- Rosen, R. (1967). Two-factor models, neural nets and biochemical automata.
  *Journal of Theoretical Biology* 15(3), 282–297. The model Okamoto's work
  builds on (book ref. [721]).
- The PubMed records of the Okamoto papers, with abstracts: [1987](https://pubmed.ncbi.nlm.nih.gov/3689885/),
  [1988, McCulloch-Pitts](https://pubmed.ncbi.nlm.nih.gov/3382700/),
  [1988, monocyclic](https://pubmed.ncbi.nlm.nih.gov/18587751/),
  [1989](https://pubmed.ncbi.nlm.nih.gov/2720139/),
  [1990](https://pubmed.ncbi.nlm.nih.gov/2224069/).
