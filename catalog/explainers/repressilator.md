## Introduction

The repressilator is a genetic clock that was designed on paper and then
built inside living bacteria. Michael Elowitz and Stanislas Leibler
published it in *Nature* in 2000. It consists of three genes, each of which
makes a *repressor*: a protein that sits on the control region of another
gene and stops that gene from being read. The three are wired in a ring.
The first protein blocks the second gene, the second blocks the third, and
the third blocks the first. In their construction the three repressors were
LacI (from *E. coli*), TetR (from a tetracycline-resistance element) and cI
(from phage λ), none of which belongs to a natural clock.

Why should a ring of "stop" signals oscillate? Follow one protein round the
loop. When protein 1 is abundant, gene 2 is shut off, so protein 2 fades
away. With protein 2 gone, gene 3 is free, so protein 3 builds up. Protein 3
then shuts off gene 1, protein 1 fades, gene 2 is released, and the story
repeats one step further round the ring. Whether the circuit keeps going
round like this, or settles at a compromise where all three proteins sit at the same
middling level, depends on the rates. Elowitz and Leibler used a simple
mathematical model to find which experimental knobs favour oscillation, and
then built the circuit in *E. coli*, with a fluorescent protein as a readout
of its state in single cells.

Banzhaf and Yamamoto present it in their chapter on wet artificial
chemistries (book §19.3.2, "Computing with Gene Regulatory Networks") as a
well-known example of a gene regulatory network (GRN) programmed to perform
a function, and as the first synthetic oscillator built in the wet lab from
bacterial genes. What Chemart holds is not the living circuit but the
book's reaction model of it: twelve species and eighteen mass-action
reactions, meant to be integrated as differential equations or simulated
molecule by molecule. It is a small, fixed network, not a constructive
chemistry.

Its closest neighbours in the catalog are the other oscillators and the
building block it relies on. [The Brusselator](brusselator.md) is an
invented minimal chemical oscillator and [the Oregonator](oregonator.md) a
reduced model of a real chemical one, the Belousov-Zhabotinsky reaction;
the repressilator is the biological, engineered counterpart, and the book
itself names it alongside the BZ reaction as a well-known biochemical
oscillator. Its repression step is the cooperative binding described in
[Hill kinetics](hill-kinetics.md), and the circuit does not oscillate
without that cooperativity.

## How it works

### The molecules

The model has four kinds of molecule for each of the three genes, twelve
species in all:

- `G1`, `G2`, `G3`: a gene in its free state, able to be transcribed;
- `C1`, `C2`, `C3`: the same gene with repressor bound to it, silent;
- `M1`, `M2`, `M3`: messenger RNA (mRNA), the working copy read off a gene;
- `P1`, `P2`, `P3`: the repressor proteins made from that mRNA.

In the real circuit `P1`, `P2` and `P3` are LacI, TetR and cI, in that
order around the ring. The model treats all three genes alike: same rates,
different targets.

### The reactions

Each gene goes through the same cycle, and the ring is closed by which
protein binds which gene. For gene 2, whose repressor is `P1`:

```
G2 + 2 P1 -> C2      binding (rate constant ke): two P1 shut gene 2 off
C2 -> G2 + 2 P1      unbinding (kr): gene 2 is free again
G2 -> G2 + M2        transcription (km): a free gene makes mRNA
M2 -> M2 + P2        translation (kp): mRNA makes protein
M2 -> ∅              mRNA decay (mu_m)
P2 -> ∅              protein decay (mu_p)
```

Gene 3 is repressed by `P2` and gene 1 by `P3`, which closes the ring. The
gene and the mRNA appear on both sides of transcription and translation:
they are catalysts, used but not consumed. Only the free gene is
transcribed, so a gene held in `C` makes nothing. Since binding and
unbinding just move a gene between its two states, `G_i + C_i` stays at its
starting value of 1 for each gene; every other species is made and
destroyed. Decay stands in for the degradation of mRNA and protein in the
cell, and it is what keeps the populations bounded.

The `2` in the binding reaction is the **Hill coefficient** `n`: the number
of protein molecules that must bind together to silence the gene. At
equilibrium, binding and unbinding balance when `ke · G · P^n = kr · C`,
so the fraction of the gene that is free is `1 / (1 + (ke/kr) · P^n)`. That
is the Hill repression curve with a threshold `K = (kr/ke)^(1/n)`, which is
1 with the default rates. For `n = 1` it falls off gently as protein
accumulates; for `n = 2` or more it behaves more like a switch, near fully
on below the threshold and near fully off above it. The book writes binding
as this elementary reaction, and Chemart keeps it that way instead of
replacing it by a Hill rate law (see the decisions below).

### A worked example: one run of the default network

The default network starts as the book's figure 19.18 does: gene 1 free
(`G1 = 1`), genes 2 and 3 repressed (`C2 = C3 = 1`), and no mRNA or protein.
Integrating the rate equations of the eighteen reactions (the script is under
*Using it* below) gives this sequence of events.

Gene 1 starts making `M1` and then `P1`. The two complexes, meanwhile,
fall apart at rate `kr = 1`, releasing their bound protein: `C2` releases
two `P1` and `C3` two `P2`. So `P1` and `P2` appear first and peak early,
at about t = 5 to 6. From there the ring takes over, and the proteins peak
in turn (peak times up to t = 200, the start-up peaks left out):

```
P2 peaks at t = 47.9, 125.3
P1 peaks at t = 73.0, 152.0
P3 peaks at t = 98.9, 178.9
```

Each protein peaks about 26 time units after the one before, in the order
`P1`, `P3`, `P2`, `P1`, ..., each reaching a maximum near 29. The order runs
*against* the numbering of the ring, for the reason told in the
introduction: when `P1` is high it shuts off gene 2, which releases gene 3,
so the next protein to rise is `P3`. Three such steps make one period,
about 81 time units once the start-up has died away. At late times `P1`
swings between 1.2 and 29.7, and the free fraction of gene 1 between 0.001
and 0.40: most of the time each gene is almost entirely bound.

With `n = 1` the same network does not oscillate: after a long run `P1`
sits at a constant level of about 9.5. The switch-like response of cooperative binding is what keeps
the loop from settling.

The book's figures label time in seconds, but its rate constants are not
the values Elowitz and Leibler used (their model had a 10-minute protein
half-life and a 2-minute mRNA half-life, and their bacteria oscillated with
periods of hours), so read the time unit here as arbitrary.

### Two ways to run it

The formal specification below lists two reactors. The deterministic one
(`ode`) treats each species as a concentration and integrates the rate
equations, as in figure 19.18. The stochastic one (`ssa`, Gillespie's
stochastic simulation algorithm) treats them as whole molecules and fires
one reaction at a time at random, with probabilities set by the rate
constants, as in figure 19.19. For the stochastic version the book picks a
small volume `V = m / N_A`, where `N_A` is Avogadro's number, so that a
concentration of 1 corresponds to `m` molecules; with `m = 100`, the run
starts with 100 copies each of `G1`, `C2` and `C3`. Converting each rate
constant to a per-molecule constant follows Wolkenhauer et al. (2004).

## Using it

The default call above gives the network of the book's figures 19.18 and
19.19: `n = 2`, `ke = kr = kp = 1`, `km = 5`, `mu_m = 0.5`, `mu_p = 0.1`,
and the starting state in `net.initial_state`. `net.extras` is empty. The
network involves no randomness, so the `seed` changes nothing. Chemart
supplies the network, not a simulator, so each recipe below carries its own
short integrator.

#### Oscillation with n = 2, none with n = 1

```python
import numpy as np
from scipy.integrate import solve_ivp

import chemart


def simulate(net, t_end, n=20001):
    ids, R, P = net.matrices()
    R = R.toarray()
    S = (P.toarray() - R).astype(float)
    k = np.array([r.rate["k"] for r in net.reactions])
    x0 = np.array([net.initial_state[s] for s in ids])

    def f(t, x):
        return S @ (k * np.prod(x[:, None] ** R, axis=0))  # mass action

    t = np.linspace(0, t_end, n)
    sol = solve_ivp(f, (0, t_end), x0, t_eval=t, method="LSODA", rtol=1e-8, atol=1e-10)
    return t, dict(zip(ids, sol.y))


def peaks(t, y):
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    return t[i]


for n in (2, 1):
    t, x = simulate(chemart.generate_network("repressilator", n=n), 2000)
    late = t >= 1800
    P1 = x["P1"][late]
    print(f"n={n}: P1 over t 1800-2000 from {P1.min():.2f} to {P1.max():.2f}")
```

```
n=2: P1 over t 1800-2000 from 1.24 to 29.69
n=1: P1 over t 1800-2000 from 9.51 to 9.51
```

Both runs together take about four seconds.

#### Finding the edge of the oscillating region

With the same `simulate` and `peaks`, change one rate at a time and measure
the swing of `P1` (maximum minus minimum) over the last 200 time units. A
swing of zero means the run has settled:

```python
for params in ({"n": 3}, {"mu_p": 0.05}, {"mu_p": 0.5}, {"km": 0.5}, {"km": 0.2},
               {"mu_m": 2.0}, {"mu_m": 3.0}, {"kr": 50.0}, {"kr": 500.0}):
    t, x = simulate(chemart.generate_network("repressilator", **params), 2000)
    P1 = x["P1"][t >= 1800]
    swing = np.ptp(P1)
    p = peaks(t, x["P1"])
    period = f", period {np.diff(p[p > 1000]).mean():.1f}" if swing > 0.01 * P1.mean() else ""
    print(params, f"P1 swing {swing:.2f}, mean {P1.mean():.2f}{period}")
```

```
{'n': 3} P1 swing 56.33, mean 16.86, period 111.4
{'mu_p': 0.05} P1 swing 34.50, mean 14.18, period 144.7
{'mu_p': 0.5} P1 swing 11.44, mean 4.66, period 28.0
{'km': 0.5} P1 swing 1.40, mean 2.13, period 72.2
{'km': 0.2} P1 swing 0.00, mean 1.38
{'mu_m': 2.0} P1 swing 0.72, mean 2.83, period 49.5
{'mu_m': 3.0} P1 swing 0.00, mean 2.42
{'kr': 50.0} P1 swing 30.45, mean 20.87, period 63.1
{'kr': 500.0} P1 swing 0.00, mean 32.34
```

The oscillation stops when transcription is weak (`km` between 0.5 and
0.2), when mRNA decays much faster than protein (`mu_m` between 2 and 3,
against `mu_p = 0.1`), or when repression is weak (`kr = 500` raises the
threshold `K` to about 22). The period follows the protein lifetime most
closely: halving `mu_p` from 0.1 to 0.05 stretches it from 81 to 145.
This run takes about 13 seconds.

#### A stochastic run, as in figure 19.19

```python
from math import comb

import numpy as np

import chemart
from chemart.kinetics import k_to_c

net = chemart.generate_network("repressilator", seed=1)
m = 100  # gene copies; volume chosen so that N_A * V = m (book fig. 19.19)
n = {s: round(c * m) for s, c in net.initial_state.items()}  # G1 = C2 = C3 = 100
c = [k_to_c(r.rate["k"], r.reactants, volume=m, avogadro=1.0) for r in net.reactions]
print("binding c:", c[0], " transcription c:", c[6])

rng = np.random.default_rng(1)
t, events, next_sample = 0.0, 0, 0.0
samples = []
while t < 200.0:
    a = np.array([ci * np.prod([comb(n[s], k) for s, k in r.reactants.items()])
                  for ci, r in zip(c, net.reactions)])
    total = a.sum()
    t += rng.exponential(1 / total)
    while next_sample <= t and next_sample <= 200.0:
        samples.append((next_sample, n["G1"], n["P1"], n["P2"], n["P3"]))
        next_sample += 1.0
    r = net.reactions[rng.choice(len(a), p=a / total)]
    for s, k in r.reactants.items():
        n[s] -= k
    for s, k in r.products.items():
        n[s] += k
    events += 1

s = np.array(samples)
print(f"{events} events")
for j, name in ((1, "G1"), (2, "P1"), (3, "P2"), (4, "P3")):
    late = s[s[:, 0] >= 100, j]
    print(f"{name} over t 100-200: {late.min():.0f} to {late.max():.0f}")
t1 = s[:, 0]; y = s[:, 2]
hi = y > 1500
rises = t1[1:][hi[1:] & ~hi[:-1]]
print("P1 rises above 1500 at t =", rises)
```

```
binding c: 0.0002  transcription c: 5.0
275604 events
G1 over t 100-200: 0 to 40
P1 over t 100-200: 121 to 2745
P2 over t 100-200: 138 to 2768
P3 over t 100-200: 117 to 2533
P1 rises above 1500 at t = [ 60. 135.]
```

The binding constant is `2 / m²`: the `1/m²` for a reaction of three
molecules, times `2!` for the two identical protein molecules. Proteins
peak near 2,700 molecules, about `m` times the deterministic peak of 29, as
in the book's figure, and `P1` crosses 1,500 on the way up at t = 60 and
135, 75 time units apart, close to the deterministic period of 81. The
gene counts are noisier than the protein counts because there are far
fewer gene molecules. The run takes about 45 seconds; the number of
events, and so the time, grows in proportion to `m`.

The parameter table below lists the Hill coefficient and the six rate
constants.

## Results

**The circuit oscillates in living cells.** Elowitz and Leibler (2000)
built the repressilator on a low-copy plasmid in *E. coli*, with a second
plasmid carrying a green fluorescent protein (GFP) gene under a promoter
repressed by TetR, so that the cell glows when TetR is low. To bring the
repressors' lifetimes closer to that of mRNA (about 2 minutes on average in
*E. coli*), they tagged each one for rapid destruction by the cell's
proteases. They followed single cells under the microscope as they grew
into microcolonies. In one cell shown in the paper the fluorescence
oscillated with a period of around 150 minutes, "roughly threefold longer
than the typical cell-division time". At least 40% of cells were
oscillatory in each of three movies, and the peak-to-peak intervals were
160 ± 40 minutes (mean ± standard deviation, 63 intervals). Because the
period is longer than the cell cycle, the oscillator's state has to be
passed from mother to daughter cells: sibling cells stayed correlated with a
decorrelation half-time of 95 ± 10 minutes, longer than the 50 to 70 minute
division time. The oscillations were noisy, with large variations in period
and amplitude from cell to cell and over time. The cell-division cycle did
not appear to be coupled to the repressilator, but the oscillator halted
when the colony entered the stationary phase.

Chemart's model is the book's reaction scheme, not Elowitz and Leibler's
own equations, and its rates are not fitted to the experiment, so it
reproduces none of these measured numbers. The YAML lists "implemented in
real E. coli" as a phenomenon; it is a fact about the construction, not
something a simulation can check.

**Oscillation needs the right parameters.** Elowitz and Leibler's design
model has two kinds of solution: a stable steady state, or an unstable one
surrounded by sustained oscillations (a *limit cycle*). They found that
oscillations are favoured by strong promoters with efficient translation,
tight repression, cooperative repression, and comparable protein and mRNA
decay rates, and that the oscillating region grows as the Hill coefficient
increases. They also note that in such networks the period is set mainly by
protein stability. The book shows the same with its own model: figure
19.18 oscillates with the default rates, and "no oscillations occur for
n = 1". Chemart's tests check both statements. They integrate the default
network and require `P1` to swing by more than 1 between t = 300 and 400,
and they integrate the `n = 1` network to t = 2000 and require the swing
over the last 200 time units to be under 1% of the mean. The recipe above
shows the other directions (weak transcription, fast mRNA decay, weak
repression) moving the book's model out of the oscillating region, and the
period tracking the protein decay rate, in line with the paper; no test
checks those.

**Oscillations survive molecular noise.** Elowitz and Leibler also ran a
stochastic version of their model. It still oscillated, but noise shortened
the correlation time of the oscillation from infinite, in the continuous
model, to about two periods. The book's figure 19.19 makes the same point
with its own model: with 100 gene copies, a Gillespie simulation keeps
oscillating, with noise more visible in the gene counts than in the protein
counts. The stochastic run above reproduces that figure; no test checks it.
Loinger and Biham (2007), whom the book cites for a more thorough study of
noise, found that fluctuations change the range of conditions in which
oscillations appear, and their amplitude and period, and that the
deterministic and stochastic descriptions agree only when every component,
including free and bound proteins and plasmids, is present in large
numbers. In the variant most like the book's (cooperative binding, an mRNA
step, and bound repressors that are not degraded) they report that the rate
equations predict oscillations that the stochastic analysis does not show.
Their model and rates differ from the book's, and Chemart does not
reproduce their study. With the book's rates, a stochastic run of the
Chemart network with a single gene copy (`m = 1`) did still cycle, though
irregularly: 16 rises of `P1` from below 3 to above 20 in 2,000 time units,
against about 25 periods of the deterministic cycle.

**Later work.** The book describes one extension. Garcia-Ojalvo, Elowitz and
Strogatz (2004) modelled bacteria whose repressilators are coupled through
*quorum sensing*, a signalling molecule that cells release and sense. Their
modelling showed that the coupled cells synchronise, giving a population
clock more accurate than the noisy clock of a single cell. Chemart does not
include this coupled circuit.

## Further reading

- Elowitz, M. B. & Leibler, S. (2000). A synthetic oscillatory network of
  transcriptional regulators. *Nature* 403, 335–338.
  <https://doi.org/10.1038/35002125>
- Loinger, A. & Biham, O. (2007). Stochastic simulations of the
  repressilator circuit. *Physical Review E* 76(5), 051917.
  <https://arxiv.org/abs/0710.1421>
- Wolkenhauer, O., Ullah, M., Kolch, W. & Cho, K.-H. (2004). Modelling and
  simulation of intracellular dynamics: choosing an appropriate framework.
  *IEEE Transactions on NanoBioscience* 3(3), 200–207. The conversion of
  rate constants used for the stochastic run.
- Garcia-Ojalvo, J., Elowitz, M. B. & Strogatz, S. H. (2004). Modeling a
  synthetic multicellular clock: repressilators coupled by quorum sensing.
  *Proceedings of the National Academy of Sciences* 101(30), 10955–10960.
