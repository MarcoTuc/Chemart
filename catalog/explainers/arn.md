## Introduction

In a living cell, genes do not simply make proteins at a fixed rate. Some of the
proteins they make, the *transcription factors* (TFs), bind to the DNA next to
other genes and switch those genes up or down. The genes and their mutual
switching form a *gene regulatory network* (GRN). The same genome can then
behave in many different ways, and a small change in the DNA can change *when*
and *how much* a gene is expressed rather than *what* it makes.

The Artificial Regulatory Network (ARN) is Wolfgang Banzhaf's abstract model of
such a network, introduced in two papers of 2003 (book refs [67, 68]). Genes and
proteins are both bitstrings. One long random bitstring is the genome; a fixed
8-bit pattern marks where each gene starts; each gene is translated into a
32-bit protein; and every protein binds to every gene's two regulatory sites,
one that boosts the gene and one that represses it, with a strength set by how
*complementary* the two bitstrings are: the more bits differ, the tighter the
binding. Nobody draws the network. It is read off the genome, so mutating a bit
of the genome rewires the network.

According to the book (§11.3.4), the initial motivation was to give artificial
evolution a developmental step, a dynamic mapping from genotype to phenotype,
in order to make it more scalable and adaptive. Banzhaf's ECAL 2003 paper
presents the model as able to reproduce phenomena of natural regulatory
networks, and argues that ARNs are useful new genetic representations for
artificial evolution. Its showcase is
*heterochrony*, a shift in the timing of gene expression that biologists
consider important in development, which the ARN produces from one- or two-bit
changes. Later work with P. Dwight Kuo and André Leier grew genomes by
repeated whole-genome duplication and mutation and found network shapes and
small recurring wiring patterns resembling those of real regulatory networks.
The model was also used as a representation for genetic programming, the
evolution of programs.

What kind of thing it is: a simulation model. A genome defines a fixed set of
proteins and a system of ordinary differential equations (ODEs) for their
concentrations. Banzhaf and Yamamoto describe it in chapter 11, among the
network chemistries inspired by biology (§11.3.4), and return to it for GRN
modelling (§18.4.3) and the evolution of morphogenesis (§18.6.1). Its nearest
neighbours in the catalog are [random Boolean networks](rbn.md), where genes
are on/off switches wired at random rather than decoded from a sequence;
[Hogeweg's cellular Potts evo-devo model](cpm-grn-evodevo.md), where a small
Boolean GRN inside each cell drives morphogenesis; and the
[repressilator](repressilator.md), an engineered three-gene circuit modelled
with [Hill kinetics](hill-kinetics.md), whose rate laws are written by hand. What sets the ARN apart is
the sequence: every interaction strength follows from bit patterns, so the
network topology and its kinetics are both products of the genome.

## How it works

### The genome and its genes

The genome is a string of bits, 4,096 of them by default. To find the genes,
Chemart scans it for the *promoter* `01010101`. At random, any 8-bit pattern
appears at a given position with probability 2⁻⁸ ≈ 0.39%. Because the promoter
is periodic, a longer run such as `0101010101` counts as one promoter, at its
first bit.

Each gene has a fixed layout, taken from Kuo, Leier and Banzhaf (2004):

```
[inhibitor site: 32 bits][enhancer site: 32 bits][promoter: 8 bits][coding region: 5 × 32 bits]
```

The two 32-bit *regulatory sites* sit just upstream of the promoter. The
160-bit coding region follows it. There is no RNA step: the coding region is
turned straight into a protein by **majority vote**. Cut it into five 32-bit
segments and stack them; bit k of the protein is 1 if at least three of the
five segments have a 1 at position k. A promoter too close to either end of the
genome to fit a whole gene is ignored. Genes may overlap, and each is treated
independently. A random 4,096-bit genome carries about ten genes; in
general about one gene per 340 bits.

In the default run (seed 1) the genome holds seven genes, and their proteins
are named `P1` to `P7`. Each species records where its gene starts and its
three bitstrings:

```
P1 position=800 inhibitor=00010000110101001001010010111010 enhancer=10010100100001010001000010000111 protein=10100100000000011111011001001001
```

### How strongly a protein regulates a gene

Every protein `j` is compared with the enhancer site and the inhibitor site of
every gene `i`. The comparison is a bitwise XOR: a bit counts when the protein
and the site differ there. The number of such complementary bits is the
**match** `d`, between 0 and `d_max = 32`. For `P1`'s own enhancer site and
`P1`'s protein:

```
enhancer  10010100100001010001000010000111
protein   10100100000000011111011001001001
XOR       00110000100001001110011011001110   -> 14 complementary bits
```

The influence of the protein on the site is `exp(β (d − d_max))`, where β is a
positive scaling factor (1 by default). A perfect complement (`d = 32`) gives
influence 1; each missing bit divides it by e ≈ 2.7. So `P1` on its own
enhancer, at 14 bits, weighs `exp(−18)`, while the best pair in this genome,
`P2`'s protein on `P7`'s enhancer at 24 bits, weighs `exp(−8)`, about 22,000
times more. Between random bitstrings the match is spread roughly like a bell
curve around 16, so a few interactions dominate and most are weak.

The book points out that the exponential has the form of the Arrhenius law of
chemistry, the rate of a reaction as a decreasing exponential of its
activation energy: read `d_max − d` as the activation energy and β as the
thermal energy RT. A good match lowers the barrier.

### The dynamics: proteins competing on a simplex

Each protein has a concentration `c_i`. The concentrations always add up to 1,
so the state of the network is a point on a *simplex*, and a protein can only
gain share at the expense of others. For each gene `i` the model adds up the
influence of all proteins on its two sites, weighted by their concentrations:

- the enhancement `a_i = (1/N) Σ_j c_j exp(β (d^a_ij − d_max))`,
- the inhibition `h_i = (1/N) Σ_j c_j exp(β (d^h_ij − d_max))`,

where `N` is the number of genes and `d^a_ij`, `d^h_ij` are the matches of
protein `j` with gene `i`'s enhancer and inhibitor sites. Protein `i` then
grows at rate `δ (a_i − h_i) c_i` (the book divides this by the total `Σ_j c_j`,
which is 1), with δ a positive overall rate factor (1 by
default). Proteins whose genes are more enhanced than inhibited, relative to
the others, gain share. The division by N is how Banzhaf (2003) models
competition between genes for the raw material of production.

Chemart writes these equations as a reaction network. Every term
`k c_j c_i` is one mass-action reaction with `k = δ/N · exp(β (d − d_max))`:

```
P_j + P_i -> P_j + 2 P_i     enhancement: protein j catalyses the production of protein i
P_j + P_i -> P_j             inhibition: protein j catalyses the removal of protein i
```

with `2 P_i` on the left when a protein regulates its own gene. On top of these,
a non-selective outflow (Chemart's `constant-total` dilution) removes every
protein in proportion to its concentration, exactly enough to keep the total at
1; this is the flow term Φ of Banzhaf (2003). So a genome with N genes gives
2 × N × N reactions: every protein on every gene, once per site.

The first reaction of the default run is `P1` enhancing itself:

```
2 P1 -> 3 P1  [mass-action k=2.1757113921018042e-09 site=enhancer match=14]
```

With N = 7, `k = (1/7) · exp(14 − 32) = 2.18 × 10⁻⁹`. Each reaction carries its
site and match next to its rate constant.

### Where the genome comes from

Two ways of making a genome are offered.

- **Random**: every bit drawn independently, as in Banzhaf (2003).
- **Duplication and divergence**: as in Kuo and Banzhaf (2004), start from a
  random 32-bit string, copy the whole string onto its own end again and again,
  and after each doubling flip every bit with a small probability, the
  *mutation rate*. Twelve doublings give 32 × 2¹² = 131,072 bits. Genes then
  arise from copies of earlier genes, from mutations, or across the join of two
  copies, so a genome contains families of related genes.

### From weighted network to graph

The full model is a complete weighted network: every protein acts on every
gene. To study its *shape*, the papers keep only the interactions whose match
reaches a threshold, and call each surviving pair a link. The book describes
this as the static view of an ARN. At threshold 0 every link is present; at 32
only perfect complements survive. Chemart's link threshold drops the weaker
reactions from the network, which also changes its dynamics.

## Using it

The default call above builds a random 4,096-bit genome, finds its seven genes,
and writes out all 2 × 7 × 7 = 98 regulation reactions with equal starting
concentrations of 1/7. It does not reproduce a particular published genome: the
papers show examples, not seeds. `net.extras["genome"]` holds the whole genome
as a string of `0`s and `1`s, and `net.extras["analysis"]` gives the gene count,
the genome length and `d_max`: `{'genes': 7, 'genome_length': 4096, 'd_max': 32}`.
The bitstrings of each gene are in its species' `structure`, as shown above.

#### Running the dynamics

Chemart generates networks; it does not simulate them. The reactions carry
everything needed, so a few lines with SciPy integrate them. This script
rebuilds the two matrices of rate constants from `net.reactions` and applies
the constant-total outflow:

```python
import numpy as np
from scipy.integrate import solve_ivp
import chemart

net = chemart.generate_network("arn", seed=1)
ids = [s.id for s in net.species]
pos = {s: n for n, s in enumerate(ids)}
N = len(ids)

# A[i, j] and H[i, j]: rate constant of protein j on the enhancer / inhibitor of gene i
A, H = np.zeros((N, N)), np.zeros((N, N))
for r in net.reactions:
    change = {s: r.products.get(s, 0) - r.reactants.get(s, 0) for s in r.reactants}
    (i,) = [s for s, d in change.items() if d]          # the gene whose protein changes
    j = next((s for s in r.reactants if s != i), i)     # the regulating protein
    (A if r.rate["site"] == "enhancer" else H)[pos[i], pos[j]] = r.rate["k"]

def f(t, c):
    growth = (A @ c - H @ c) * c        # each reaction runs at k * c_j * c_i
    return growth - c * growth.sum()   # the constant-total outflow keeps sum(c) = 1

c0 = np.array([net.initial_state[s] for s in ids])
t = np.linspace(0, 1e8, 11)
sol = solve_ivp(f, (0, t[-1]), c0, t_eval=t, method="LSODA", rtol=1e-8, atol=1e-10)
print("time     ", "  ".join(f"{s:>5}" for s in ids))
for k in range(len(t)):
    print(f"{t[k]:9.1e}", "  ".join(f"{x:5.3f}" for x in sol.y[:, k]))
```

```
time         P1     P2     P3     P4     P5     P6     P7
  0.0e+00 0.143  0.143  0.143  0.143  0.143  0.143  0.143
  1.0e+07 0.001  0.004  0.001  0.002  0.000  0.000  0.992
  2.0e+07 0.000  0.003  0.000  0.006  0.000  0.000  0.991
  3.0e+07 0.000  0.003  0.000  0.021  0.000  0.000  0.975
  4.0e+07 0.000  0.005  0.000  0.081  0.000  0.000  0.913
  5.0e+07 0.000  0.009  0.000  0.100  0.000  0.000  0.891
  6.0e+07 0.000  0.007  0.000  0.071  0.000  0.000  0.922
  7.0e+07 0.000  0.007  0.000  0.079  0.000  0.000  0.913
  8.0e+07 0.000  0.007  0.000  0.082  0.000  0.000  0.911
  9.0e+07 0.000  0.007  0.000  0.079  0.000  0.000  0.914
  1.0e+08 0.000  0.007  0.000  0.080  0.000  0.000  0.913
```

`P7`, whose enhancer is matched in 24 bits by `P2`, takes over almost at once.
Later `P4`, which inhibits `P7` at 21 bits, rises, overshoots to 0.100 and
settles near 0.080 after a damped swing: the system ends at a point attractor
with two surviving proteins. The run takes a couple of seconds.

**Time units.** With β = δ = 1 and typical matches around 16 bits, rate
constants are around 10⁻⁸, so the interesting changes happen over 10⁶ to 10⁸
time units. Every rate constant is proportional to δ, so raising `delta`
rescales time and nothing else. Raising `beta` sharpens the contrast between
good and bad matches. The papers give no numeric values for either.

**Other genomes behave differently.** Swapping the seed changes the genome and
with it the dynamics. Over seeds 0–11 at the default length, 11 of the 12 runs
end (at t = 10⁸) with one protein holding 0.87 to 1.0 of the total, as in the
papers. Seed 9 (13 genes) does
not settle: in the same script with `seed=9`, four proteins (`P8`, `P10`,
`P12`, `P13`) take turns to dominate in a repeating cycle, with a period of
about 3.7 × 10⁷ time units that was still steady at t = 4 × 10⁸. Banzhaf (2003)
reports point attractors and damped oscillations; sustained cycles like this
one are an observation from Chemart runs, not a published result.

#### Recipes

**The genomes of Kuo and Banzhaf.** `genome_length=131072` with the default
random source gives the 131,072-bit random genomes of their Fig. 5 (342–408
genes over seeds 0–19). Such a genome has about 380 genes and therefore about
290,000 reactions, which takes about 5 seconds to generate; add a link
threshold to keep it small. For duplication and divergence, set
`genome_source="duplication-divergence"`; `genome_length` must then be 32 times
a power of two, and `mutation_rate` is the per-bit flip probability after each
doubling:

```python
net = chemart.generate_network("arn", seed=2, genome_source="duplication-divergence",
                               genome_length=32768, link_threshold=22)
len(net.species), len(net.reactions)        # (193, 1139)
```

That is the setting of Kuo, Leier and Banzhaf (2004): ten doublings at 1%
mutation. The gene count varies wildly between seeds (10 to 193 genes over
seeds 0–7), and a genome full of near-copies may have no link at all above the
threshold (seed 1: 97 genes, 0 reactions at 22 bits). Small
duplication-divergence genomes, such as the default 4,096 bits, often carry no
gene at all and give an empty network.

**Counting genes only.** To study gene numbers without building the reactions,
pass `link_threshold=33`: no match can reach it, so the network has species but
no reactions, and `net.extras["analysis"]["genes"]` is the count. Over 40
seeds of 131,072 bits, duplication and divergence at 1% gave from 15 to 2,092
genes (median 151); at 5% it gave 234 to 609, close to random genomes.
Generating the 80 genomes takes about 10 seconds.

**The connectivity transition.** Sweep `link_threshold` from 0 to 32 on one
genome and divide the number of reactions by 2N². For seed 3 with 16,384 bits
(N = 49 genes):

```
threshold  0-6    8      10     12     14     16     18     20     22     24     26+
links      1.000  0.999  0.992  0.946  0.815  0.568  0.300  0.108  0.026  0.003  0.000
```

The parameter table below lists the ranges used in the papers.

## Results

### Banzhaf (2003): rich dynamics and heterochrony from random genomes

Banzhaf's ECAL 2003 paper, the source of the dynamical model, works with random
genomes. Its Table 1 counts genes in three of them: 3 genes in 1,000 bits, 37
in 10,000 and 409 in 100,000. It also reports the best match between any
protein and any position of the genome (25, 28 and 30 bits), and shows that
matches between a protein and a random genome follow a roughly Gaussian
distribution (Fig. 1). Chemart finds fewer genes (27 on average for 10,000
bits, 288 for 100,000, over 20 seeds): Table 1's counts are close to the raw
0.39% rate, while Chemart follows the rule of the later papers, in which
overlapping promoters count once and a gene needs room for both sites and its
coding region.

Starting from equal concentrations, the paper finds that "some proteins
increase their level of concentration, then fall again, with usually one being
left over", a point attractor. Different random genomes give "remarkably
different" dynamics: a damped oscillation, a slow smooth drift, a quick
settling and a long transition in which one protein peaks and then hands over
to another (Figs. 3–4). Chemart reproduces the main pattern: its tests
integrate four random 2,048-bit genomes and check that the concentrations stay
on the simplex, that one protein ends with more than 40% of the total (and more
than twice its fair share 1/N), and that the concentrations have stopped
changing by the end. The damped oscillation appears in the default run above,
but the tests do not check for it.

The paper's central result is **heterochrony**. In the genome of its Fig. 4,
protein 7 rises above a concentration of 0.8 at time 9,000 and falls below it at
22,000. Changing the match between protein 7 and gene 4's inhibitor site by
one bit moves these times to 14,500 and 30,500; increasing the match by another
bit moves them to 30,000 and 53,000 (Figs. 5–6). Some one- or two-bit changes have
no effect at all, a neutral variation. Because influence is exponential in the
match, the shifts grow with each additional bit. The paper argues that such
networks translate small pattern changes into changes of timing, as natural
GRNs do, and could be the algorithmic "missing link" between genotypes under
constant evolutionary change and stable phenotypes. Chemart does not reproduce
this experiment: the generator draws the genome from a seed and has no way to
edit a site, and the paper's equations scale each exponential to the best
match present rather than to 32, a variant Chemart does not implement (see the
implementation decisions).

The paper ends with a first evolution experiment: a (1+1) evolution strategy
(one parent, one mutated child, keep the better) drives one protein to a
target concentration, `c_6 = 0.085` at time 100; in three runs the deviation
falls from about 10⁻² to 10⁻⁹ or below (Fig. 7). Long stagnation periods in
fitness hide continuing neutral changes in the genome. Not reproduced.

### Kuo and Banzhaf (2004): duplication and divergence

Kuo and Banzhaf's Artificial Life IX paper introduced genomes grown by whole-
genome duplication and divergence, and asked whether the regulatory networks
they encode have the shapes found in many natural networks. A network is
*scale-free* when the number of nodes with k links falls off as a power law
k^−γ, so that a few hubs have very many links; it is *small-world* when any two
nodes are a few steps apart and yet neighbours of a node tend to be linked to
each other (a high *clustering coefficient*).

- **Gene numbers.** Over 200 genomes of 131,072 bits, a 1% mutation rate gives
  gene counts from near 0 to about 3,000, which the paper fits with a power law
  of exponent γ = 0.9779 (Fig. 3). At 5% the distribution becomes narrow, 200 to
  700 genes (Fig. 4), much like fully random genomes, which carry 340 to 440
  (Fig. 5). The paper's explanation: at 5% an 8-bit promoter survives a
  duplication only 66% of the time (0.95⁸), so mutation undoes what duplication
  copies. Chemart's tests reproduce all three: random genomes of 131,072 bits
  stay within 340–440 genes over 30 seeds; at 1% the counts spread beyond 700
  and below 200 but stay under 3,000; at 5% they stay within 200–700; and the
  relative spread at 1% is more than three times that at 5%.
- **Scale-free and small-world topologies.** For each genome the paper tries
  every threshold, keeps the one whose degree distribution best fits a power
  law with γ near 2.5, and finds many duplication-divergence genomes at 1%
  mutation that qualify as scale-free; the 2006 follow-up adds that most
  networks made at 5% do not. In the majority of genomes some threshold also
  gives a small-world network. The 2004 paper argues that duplicating a whole
  genome acts like *preferential attachment*, the classic
  mechanism for scale-free networks: nodes that already have many links gain
  the most when everything is copied. Chemart generates the same networks but
  its tests do not measure degree distributions, clustering or path lengths.

### Kuo, Leier and Banzhaf (2004): evolving time series

A network's proteins have no meaning outside it, and their concentrations are
bound to sum to 1. To get a usable output signal, the PPSN 2004 paper adds one
extra pair of regulatory sites at a random place in the genome, with no gene of
its own, and uses the summed enhancement minus inhibition at those sites,
normalised to [−1, 1], as the output. A (50+100) evolution strategy (50
parents, 100 mutated children, the best 50 survive) then evolved genomes of
32,768 bits, made by ten doublings at 1% mutation, towards three target
curves: `sin(t)`, a decaying exponential `2 exp(−0.1 t) − 1` and a sigmoid. All
30 runs, 10 per target, reached a mean squared error of at most 0.0059, with the
winning networks using between 20 and 234 genes. The paper concludes that ARN
dynamics are evolvable and could serve to generate arbitrary time series.
Not reproduced: evolution acts on genomes across many networks, and Chemart
generates one network per genome. The output sites are not implemented.

### Kuo, Banzhaf and Leier (2006): thresholds, motifs and evolved dynamics

The BioSystems paper collects and extends this work.

- **The connectivity transition.** Over 200 networks, the fraction of links
  present falls from 1 to 0 as the threshold rises from 0 to 32, with "a sharp
  transition from full connectivity to no connectivity" (Fig. 4); the
  scale-free and small-world networks all come from this transition region.
  Chemart's tests reproduce it on one genome: full connectivity at threshold 0,
  none at 32, a monotone fall, more than 30% of links at 8 bits and under 1% at
  24; and they check that the graph at 22 bits is a subgraph of the graph at 21,
  as the paper's Figs. 2 and 3 are.
- **Network motifs.** *Motifs* are small wiring patterns, of three or four
  genes, that occur far more often than chance. Counting all such subgraphs in
  800 duplication-divergence networks and 800 random-genome networks, the paper
  compares them with the transcription networks of the bacterium *E. coli* and
  of yeast (*S. cerevisiae*). The duplication-divergence networks are very close
  to yeast (sum of squared errors 0.0072 for three-node and 0.0984 for four-node
  subgraph distributions), and the two most frequent natural motifs, forms of
  the *single-input module* in which one gene regulates several others, are
  well represented, while only one of them is detected in random networks. The
  paper suggests that the shape of natural networks may partly come from how
  they were created, not only from later selection. Leier, Kuo and Banzhaf
  (2007, book ref [504]) explain the preference for particular motifs by the
  duplication-divergence process. Not reproduced: Chemart does not count
  motifs.
- **Evolving target dynamics.** The paper repeats the evolution experiments
  described above, and adds that two genes suffice for an oscillation or a
  sigmoid, and one for a decaying exponential; the evolved networks use far
  more genes than that.

### Later work

The book lists further uses. Lopes and Costa's ReNCoDe (2011, 2012; book refs
[521, 522]) uses an ARN for genetic programming. For morphogenesis, Chavoya and Duthen (2008, book ref [172]) extended the
ARN with structural genes, each coding for a cell type that locks in when its
regulating protein crosses a threshold, and evolved it to produce concentric
squares and French flags (book §18.6.1); this extension is not implemented. The
book also cites the ARN [68] as the regulatory layer of controllers that
combine metabolic, signalling and genetic networks (§16.1.3). Outside the
book, Nicolau, Schoenauer and Banzhaf (2010) used a modified ARN, with added
ways of connecting inputs and outputs to it, to evolve controllers for the
pole-balancing benchmark.

### What Chemart covers

Chemart builds the network of a genome exactly as the model defines it, and its
tests check the reaction network against the book's equations 11.12–11.14 term
by term. It reproduces the published statistics of genomes (gene numbers,
connectivity versus threshold) and the settling of the dynamics into a
dominant protein. It does not reproduce the topology measures (scale-free,
small-world, motifs), the heterochrony experiment, or any of the evolution
experiments, since these act on genomes across many networks.

## Further reading

- Banzhaf, W. (2003). Artificial regulatory networks and genetic programming.
  In R. Riolo & B. Worzel (eds.), *Genetic Programming Theory and Practice*,
  chapter 4, pp. 43–61. Kluwer. The companion paper of the ECAL article (book
  ref [67]).
- Leier, A., Kuo, P. D. & Banzhaf, W. (2007). Analysis of preferential network
  motif generation in an artificial regulatory network model created by
  duplication and divergence. *Advances in Complex Systems* 10:155–172.
- Lopes, R. L. & Costa, E. (2012). The regulatory network computational device.
  *Genetic Programming and Evolvable Machines* 13(3):339–375.
- Chavoya, A. & Duthen, Y. (2008). A cell pattern generation model based on an
  extended artificial regulatory network. *BioSystems* 94(1–2):95–101.
- Nicolau, M., Schoenauer, M. & Banzhaf, W. (2010). Evolving genes to balance a
  pole. *EuroGP 2010*, LNCS 6021, pp. 196–207. Springer.
