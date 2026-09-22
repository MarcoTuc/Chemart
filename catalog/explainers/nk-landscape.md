## Introduction

The NK model is a recipe for making random *fitness landscapes* whose
ruggedness you can set with a single number. A fitness landscape assigns a
fitness (a score for how well an organism reproduces) to every possible
genotype. Picture the genotypes laid out so that genotypes one mutation apart
are neighbours, and fitness as height: evolution by small mutations then walks
uphill. If the landscape has a single smooth hill, any uphill walk reaches the
top. If it is rugged, with many peaks, a walk stops on whichever local peak it
happens to climb, and the best genotype may never be found.

Stuart Kauffman built the model, first with Levin (Kauffman and Levin, 1987)
and then in his book *The Origins of Order* (1993), to study how *epistasis*, the way the effect of one
gene depends on the state of other genes, controls this ruggedness. A genotype
is a string of N genes, each 0 or 1. Every gene makes a fitness contribution
that depends on its own state and on the states of K other genes, and the
genotype's fitness is the average of the N contributions. Because nobody knows
what these interactions really are, Kauffman simply drew each contribution at
random. With K = 0 every gene counts on its own and the landscape has one peak.
With K = N − 1 every gene depends on all the others, one mutation reshuffles
every contribution, and the landscape is completely random. The values of K in
between give a family of landscapes from smooth to rugged.

The NK model is not a chemistry in the usual sense: it has no molecules and no
reactions of its own. It is a fitness function that other evolutionary models
plug into. Banzhaf and Yamamoto present it in their chapter on modelling
biological systems (book §18.4.1), as the first of Kauffman's two binary models
of genetic regulatory networks, and describe it as "a highly idealized static
GRN model". Its dynamic sibling, in which the same N-genes-with-K-inputs
wiring is run as a network of switches, is the
[random Boolean network](rbn.md) (§18.4.2). The book also meets the NK family
in chapter 8, where the NKC model couples several NK landscapes to study
coevolution; that model is not in the catalog.

To fit Chemart's common interface, this entry wraps the landscape in a
reaction network: every genotype is a species that copies itself at a rate
equal to its fitness, with occasional point mutations. That makes it the
[quasispecies equation](quasispecies.md) of §7.2.7 with an NK fitness in place
of that entry's simple fitness functions. The landscape itself, which is what
the book defines, is delivered next to the network.

## How it works

### Building a landscape

A landscape is fixed by two random choices, made once.

**Who influences whom.** Each gene i gets K *epistatic partners*, the other
genes whose states change its contribution. Kauffman used two layouts: the K
genes that follow i along the genome, wrapping round from the end to the start
(*adjacent* neighbourhoods), or K genes picked at random (*random*
neighbourhoods).

**What each combination is worth.** Gene i's contribution depends on K + 1
bits (its partners' and its own), which can take 2^(K+1) combinations. For
each gene, and for each combination, a number is drawn uniformly between 0 and
1. The book writes this table as `M(i, j)`, where `i` is the gene and `j` is the
combination read as a binary number, partners' bits first and the gene's own
bit last.

The fitness of a genotype is then

    W = (1/N) × (w_0 + w_1 + … + w_(N−1)),   with  w_i = M(i, j_i)

(book eq. 18.22), where `j_i` is the combination that gene i actually sees.

### A worked example

Here is a small landscape, with N = 3 genes and K = 1 partner each, from
`chemart.generate_network("nk-landscape", seed=7, N=3, K=1)`. With adjacent
neighbourhoods, gene 0's partner is gene 1, gene 1's is gene 2, and gene 2's
wraps round to gene 0 (`extras["analysis"]["epistatic_partners"]` is
`[[1], [2], [0]]`). The random table is not stored on the network; redrawing
the seed's random numbers gives it (rounded):

```
          j = 00   01     10     11
gene 0     0.6251 0.8972 0.7757 0.2252
gene 1     0.3002 0.8736 0.0053 0.8212
gene 2     0.7971 0.4679 0.3030 0.2784
```

Take the genotype `011`, read left to right as genes 0, 1 and 2:

- gene 0 is `0` and its partner, gene 1, is `1`: j = `10` = 2, so w_0 = 0.7757;
- gene 1 is `1` and its partner, gene 2, is `1`: j = `11` = 3, so w_1 = 0.8212;
- gene 2 is `1` and its partner, gene 0, is `0`: j = `01` = 1, so w_2 = 0.4679.

The fitness is the mean, 0.6883, which is what the network reports:

```
000 0.5741    100 0.5001
001 0.3661    101 0.3936
010 0.8154    110 0.4673
011 0.6883    111 0.4416
```

Flip gene 1 of `011` and you get `001`. That single mutation changes two
contributions, gene 1's own and gene 0's (gene 1 is its partner), and fitness
falls from 0.688 to 0.366. This is how epistasis makes the landscape rugged:
the larger K, the more contributions one mutation redraws. A genotype is a
*local optimum* when none of its N one-mutation neighbours is fitter. Here
only `010` qualifies, and it is also the *global optimum*, the fittest genotype
of all.

### From landscape to reaction network

The book stops at the landscape. To put it through the same interface as the
other chemistries, Chemart turns it into a *replicator-mutator* network, the
standard model of a population evolving under selection and mutation. Each of
the 2^N genotypes is a species. A genotype `g` of fitness W copies itself,
and each copy has a per-gene error probability μ (`mu`):

```
g -> 2 g        rate W × (1 − μ)^N          (an exact copy)
g -> g + h      rate W × μ × (1 − μ)^(N−1)   for each h one bit away from g
```

Copies with two or more errors are left out. The reactor is a set of ordinary
differential equations (ODEs) for the concentrations, with a *constant-total
dilution*: material is washed out in proportion to what is present, so that the
total concentration stays constant. Growth faster than average then means a
rising share, which is selection. The first reaction of the default network,
`000000 -> 2 000000` with `k = 0.53695`, is exact copying of a genotype of
fitness 0.5703 (0.5703 × 0.99^6); the six mutant reactions that follow have
rate 0.5703 × 0.01 × 0.99^5 = 0.00542.

## Using it

The default network is a landscape with N = 6 genes and K = 2 adjacent
partners, so 64 genotypes, each with one copying and six mutation reactions
(448 in all), and μ = 0.01. It does not reproduce any published experiment:
none is attached to the book's section. Species names are the genotypes, gene
0 first. The landscape is in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["epistatic_partners"]   # [[1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 1]]
a["fitness"]["000000"]    # 0.5703...   fitness of every genotype, a dict
a["local_optima"]         # ['011000', '011110', '101000', '101110']
a["global_optimum"]       # '101110'    fitness 0.8266
```

Networks up to N = 10 (1,024 genotypes, 11,264 reactions) build in a fraction
of a second. `mu=0` leaves only the copying reactions.

**Ruggedness against K.** Count local optima as K grows, averaged over seeds
0–9 with N = 10:

```python
for K in (0, 1, 2, 4, 6, 9):
    counts = [len(chemart.generate_network("nk-landscape", seed=s, N=10, K=K)
                  .extras["analysis"]["local_optima"]) for s in range(10)]
    print(K, sum(counts) / 10)
```

| K | 0 | 1 | 2 | 4 | 6 | 9 |
|---|---|---|---|---|---|---|
| local optima, adjacent | 1.0 | 4.5 | 8.5 | 25.5 | 48.8 | 93.4 |
| local optima, `topology="random"` | 1.0 | 3.2 | 9.1 | 22.8 | 42.4 | 93.3 |

At K = 9 the count matches the theoretical value for a fully random landscape,
2^N/(N + 1) = 93.1 (see Results). Both rows together take about half a
minute.

**Evolving on the landscape.** The network carries its rates and the
constant-total dilution, so `chemart.simulate` integrates it as it comes:

```python
from chemart import simulate

net = chemart.generate_network("nk-landscape", seed=1)
W = net.extras["analysis"]["fitness"]
for start in ("000000", "011000"):
    traj = simulate.ode(net, 2000, x0={start: 1.0}, points=41)
    for f in traj.frames:
        if f.t in (0, 50, 200, 2000):
            top = max(f.state, key=f.state.get)
            print(f"start {start}  t={f.t:6.0f}  most common {top} ({f.state[top]:.2f})  "
                  f"mean fitness {sum(W[s] * x for s, x in f.state.items()):.3f}")
```

Starting from `000000` alone, and from the local optimum `011000` alone (about
two seconds):

```
start 000000  t=     0  most common 000000 (1.00)  mean fitness 0.570
start 000000  t=    50  most common 100100 (0.32)  mean fitness 0.747
start 000000  t=   200  most common 101110 (0.43)  mean fitness 0.786
start 000000  t=  2000  most common 101110 (0.46)  mean fitness 0.786
start 011000  t=     0  most common 011000 (1.00)  mean fitness 0.700
start 011000  t=    50  most common 011000 (0.63)  mean fitness 0.673
start 011000  t=   200  most common 101110 (0.47)  mean fitness 0.786
start 011000  t=  2000  most common 101110 (0.46)  mean fitness 0.786
```

Both populations end on the global optimum `101110`, surrounded by a cloud of
mutants that holds the mean fitness (0.786) below the peak's own 0.827. The
population does not stay on the local optimum, because in a deterministic ODE
every mutant appears at once in a tiny amount and the fittest one then grows
exponentially. The model therefore cannot show populations stuck on local
peaks, which is what the NK literature studies; for that you need a
finite-population or adaptive-walk simulation on `extras["analysis"]["fitness"]`.

The mutation rate sets how tight the cloud is. The share of `101110` at
equilibrium (the leading eigenvector of `A`) is 0.88 at μ = 0.001, 0.46 at
0.01, 0.16 at 0.05 and 0.05 at 0.2; from μ = 0.05 on, another genotype is
slightly more common than the peak. The [quasispecies](quasispecies.md) entry
is the place to study this error threshold. Above a few per cent per gene, the
missing multiple-mutation reactions also make the network a poor
approximation.

## Results

The book reports Kauffman's findings only in words: at K = 0 the landscape is
unimodal, with one peak that hill climbing easily reaches; raising K makes it
more rugged, so that single mutations can change fitness drastically; and at
K = N − 1 it is random, with no correlation between a genotype's fitness and
its neighbours'. Chemart's tests check the first claim: with N = 6 and K = 0,
every one of five seeds has exactly one local optimum. The others are not
tested, but the local-optimum counts under *Using it* show the ruggedness
growing with K.

Altenberg's (1997) review summarises the quantitative results of Kauffman
(1993), Weinberger (1991) and Fontana et al. (1993) for *adaptive walks*, the
low-mutation picture in which a population sits on one genotype and moves to a
fitter one-mutation neighbour until none is left:

- **K = 0.** One globally attractive genotype. A walk reduces the distance to
  it by one gene per step, so it takes N/2 steps on average, and neighbours'
  fitnesses are highly correlated since only one of N contributions changes.
- **K = N − 1.** A genotype is a local optimum with probability 1/(N + 1), so
  there are 2^N/(N + 1) local optima on average. Walks are short, about
  ln(N − 1) steps. As N grows, the fitness of the optimum a walk reaches falls
  toward 0.5, the mean of the whole landscape: Kauffman's *complexity
  catastrophe*. Chemart's landscapes agree with the count (93.4 local optima
  against 93.1 expected for N = 10). In a quick check over 400 random walks
  per size, walks lasted 1.6, 1.7, 2.2 and 2.4 steps for N = 4, 6, 8 and 10
  (ln(N − 1) = 1.1, 1.6, 1.9, 2.2), and the fitness reached fell slowly, from
  0.68 to 0.66. These checks are not in the tests.
- **Intermediate K.** For small K the best local optima share many genes; this
  resemblance fades as K grows, faster for random than for adjacent
  neighbourhoods. The fitness correlation between one-mutation neighbours is
  1 − (K + 1)/N, the fraction of contributions left untouched.

Two further results concern how hard the landscapes are. Weinberger showed
that the fittest genotype of an adjacent-neighbourhood landscape can be found
by dynamic programming in time proportional to 2^K × N, while with random
neighbourhoods the problem is NP-complete (as hard as the hardest search
problems) for K ≥ 3; Thompson and Wright lowered that to K = 2. So the two
layouts, which look alike to adaptive walks, differ sharply in difficulty.
Chemart finds the global optimum by checking every genotype, which limits it to
N ≤ 10.

The model has been built on widely. Altenberg lists uses in models of
epistatic gene interaction, coevolution, genome growth and Wright's shifting
balance theory. In the coevolution model (NKC, book §8.2.1), Kauffman coupled
the landscapes of several species and proposed that coevolving ecosystems sit
at an "edge of chaos" between endless change and frozen equilibria.

What Chemart does not reproduce: none of Kauffman's published figures or
tables, no adaptive walks, no NKC coupling, and no stochastic population
dynamics. Its network form keeps only single point mutations, and its full
enumeration of genotypes, like the book's own `NKlandscape.py`, which the book
warns "does not scale to large values of N and K", is limited to small N. The
tests also do not check the table's index convention against a hand-computed
value.

## Further reading

- Altenberg, L. (1997). NK fitness landscapes. In T. Bäck, D. Fogel & Z.
  Michalewicz (eds.), *Handbook of Evolutionary Computation*, section B2.7.2.
  Oxford University Press. Open copy:
  <https://dynamics.org/Altenberg/FILES/LeeNKFL.pdf>. The source of the
  quantitative results above.
- Kauffman, S. A. & Levin, S. (1987). Towards a general theory of adaptive
  walks on rugged landscapes. *Journal of Theoretical Biology* 128, 11–45.
- Weinberger, E. D. (1991). Local properties of Kauffman's N-k model, a
  tuneably rugged energy landscape. *Physical Review A* 44(10), 6399–6413.
- Hordijk, W. & Kauffman, S. A. (2005). Correlation analysis of coupled fitness
  landscapes. *Complexity* 10(6), 41–49. Cited by the book for ways of
  computing NK fitness without storing the whole table.
- Fontana, W., Stadler, P. F., Bornberg-Bauer, E. G., Griesmacher, T.,
  Hofacker, I. L., Tacker, M., Tarazona, P., Weinberger, E. D. & Schuster, P.
  (1993). RNA folding and combinatory landscapes. *Physical Review E* 47(3),
  2083–2099.
- Weinberger, E. D. (1996). NP completeness of Kauffman's N-k model, a tuneable
  rugged fitness landscape. Santa Fe Institute Working Paper 96-02-003.
- Thompson, R. K. & Wright, A. H. (1996). Additively decomposable fitness
  functions. Cited by Altenberg (1997) as "to appear".
