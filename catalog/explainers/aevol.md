## Introduction

Aevol is a simulation platform for *in silico experimental evolution*: a
population of digital organisms that replicate, mutate and compete, so that
evolutionary questions can be put to them the way a microbiologist puts them
to a bacterial culture, but with the full history of every lineage on record.
Its original authors are Carole Knibbe, Guillaume Beslon and David Parsons,
of the LIRIS computer-science laboratory at INSA-Lyon, and it is distributed as
open-source software from Inria's GitLab
(<https://gitlab.inria.fr/aevol/aevol>). Unlike [Tierra](tierra.md) or
[Avida](avida.md), whose organisms are computer programs, an aevol organism
carries something shaped like a bacterial chromosome: a circular,
double-stranded string of 0s and 1s, most of which does nothing.

The central idea is that genes must be *found* on the chromosome, as in a
real cell. Short signal sequences mark where transcription starts and stops,
and where translation starts and stops; everything between them that is read
becomes a protein, and everything else is non-coding DNA. Each protein does
not fold into a 3-D shape but into a simple mathematical object: a triangle
over an abstract axis of "biological functions" from 0 to 1. Its position says
which function it performs, its width how many neighbouring functions it
covers, and its height how well (positive heights activate a function,
negative ones inhibit it). A wide, low triangle is a versatile but weak
protein; a narrow, tall one a specialised, efficient enzyme. The organism's
phenotype is the sum of all its triangles, and its fitness is how closely that
sum matches a target curve set by the environment.

Because the genome has this realistic layout, evolution acts not only on what
the genes do but also on the chromosome itself: how many genes it carries,
how they are arranged, and how much non-coding DNA lies between them.
Mutations include not just point changes but large rearrangements
(duplications, deletions, translocations, inversions of whole segments), and
those can copy or destroy many genes at once. The founding study, Knibbe et
al. (2007), used this to ask why genomes contain so much non-coding DNA, and
found that the mutation rate, with no cost of genome size imposed, largely
sets how much of it accumulates.

Aevol is a **simulation model**. Banzhaf and Yamamoto mention it in one
paragraph of their chapter on modelling biological systems, under protein
folding (book §18.1.2), as one of the closest artificial-chemistry examples
of protein folding, next to Suzuki's [Network Artificial Chemistry](nac.md). The
book describes its DNA, RNA and proteins, the codon table and the fuzzy-set
proteins, and cites Parsons, Knibbe and Beslon (2011) on rearrangements.
Among the catalog's neighbours, the [RNA folding AC](rna-folding-ac.md) also
maps sequences to function through folding, but with real thermodynamics and
no genome; the [artificial regulatory network](arn.md) also finds genes by
scanning a bit string for promoters, but its subject is how the genes it finds
regulate each other; the [quasispecies](quasispecies.md) model treats mutation and
selection as a reaction network over whole sequences, with no genes inside
them. What sets aevol apart is that the architecture of the genome, coding
and non-coding, is itself what evolves.

## How it works

Each organism is decoded in four steps, from chromosome to fitness. Then a
population of them evolves on a grid.

### Step 1: finding the transcribed regions

The chromosome is a circle of bases, each 0 or 1, read on both strands (the
second strand is the complement, read backwards). A **promoter** is any
22-base window that differs from the fixed consensus
`0101011001110010010110` in at most 4 positions. The fewer the mismatches `d`,
the more strongly the region is expressed: its transcription level is
`1 - d/5`, so 1.0, 0.8, 0.6, 0.4 or 0.2. The RNA runs from just after the
promoter to the first **terminator**, a stem-loop: 4 bases, a 3-base loop,
then the complement of the first 4 in reverse order (`abcd***d'c'b'a'`), which
in a real RNA would fold back on itself.

### Step 2: finding genes inside the RNA

Inside a transcript, a gene starts at a **ribosome binding site**: the
Shine-Dalgarno motif `011011`, any 4 bases, then the START codon `000`. From
there the gene is read three bases at a time (a *codon*) until the STOP codon
`001`. The six other codons are the "amino acids": `M0`=`100`, `M1`=`101`,
`W0`=`010`, `W1`=`011`, `H0`=`110`, `H1`=`111`. If two transcripts contain
the same gene, it is one protein whose level is the sum of theirs.

Here is the hand-built test chromosome from aevol's own unit tests, decoded by
Chemart (72 bases: a promoter with 2 mismatches, 5 filler bases, a gene, 6
filler bases, a terminator):

```python
from chemart.chemistries import aevol as A
g = A.bits("0101010001110110010110" + "11101"
           + "011011" + "0011" + "000" + "100" + "110" + "110" + "010" + "001"
           + "110011" + "01000001101")
rnas = A.transcribe(g)
for r in rnas:
    print(r.pos, r.basal, r.start, r.length)
for p in A.translate(g, rnas, 0.033333333):
    print(p.first_aa, p.e, p.structure())
```

```
0 0.6 22 50
40 0.6 code=100110110010 codons=4 m=0 w=0 h=-1 functional=0
```

There is one transcript: its promoter is at base 0, its level is 0.6, and the
RNA starts at base 22 and is 50 bases long. There is one protein: its first
codon is at base 40, it is expressed at level 0.6, and its coding bits read
`100 110 110 010`, that is `M0 H0 H0 W0`.

The promoter has 2 mismatches, so its level is `1 - 2/5 = 0.6`. The RNA runs
from the base after the promoter to the end of the terminator, the last base
of the chromosome. The ribosome binding site starts 5 bases into it (bases
27–39), and the gene after it has four codons before the STOP.

### Step 3: folding a protein into a triangle

The codons of each kind are read, in order, as one binary number per
parameter: the `M` codons (0 or 1 each) give the triangle's position `m`, the
`W` codons its half-width `w`, the `H` codons its height `h`. The digits are a
*Gray code*, a binary encoding in which consecutive numbers differ in a
single digit. Each number is divided by its maximum and rescaled: `m` in [0, 1], `w` in
[0, `w_max`], `h` in [-1, 1].

In the example above the gene has one `M0`, two `H0` and one `W0`: every digit
is 0, so `m = 0`, `w = 0` and `h = -1`. A triangle with no width does nothing,
so this protein is *non-functional*. A protein counts only if it has at least
one codon of each kind and non-zero width and height. Change the codons to
`M1 M0 W1 W0 H1 H1` and Chemart folds them into `m = 1.0`, `w = 0.0333` (the
maximum) and `h = 0.333`: a functional activator at the right edge of the
axis.

### Step 4: phenotype, target and fitness

Every functional protein adds its triangle, scaled by its transcription level,
to the phenotype. Activators and inhibitors are summed separately, each capped
at 1, then the inhibitors are subtracted and anything below 0 is set to 0. The
environment is a **target** curve on the same axis, a sum of Gaussian bumps
capped to [0, 1]. The default target is the one of aevol's example files: a
bump of height 1.2 at 0.52 with a deep hole of -1.4 carved into it at 0.5, and
a small bump of 0.3 at 0.8. Both curves are sampled at 300 points.

The **metabolic error** is the area between the phenotype and the target. An
organism with no functional protein has a flat phenotype, and its error is
the whole area under the target, 0.1525. Fitness is
`exp(-k × metabolic error)`, with `k` the selection pressure (1000 by
default). With so large a `k` the absolute fitness values are tiny (about
10⁻⁶⁵ for the one-gene ancestor of the default run), but only ratios matter: an organism whose
error is lower by 0.003 is about 20 times more likely to be chosen as a
parent.

A random chromosome rarely carries a working gene. Over 2,000 random genomes
of 5,000 bases, Chemart finds on average 21.8 promoters, but only 0.11 genes
and 0.02 functional proteins per genome. So the run starts, as in aevol, by
drawing random genomes until one does better than a flat phenotype, and fills
every cell with copies of it.

### Mutation and the evolutionary loop

The organisms live on a grid that wraps around at the edges (a torus), one per
cell. Each generation, every cell is refilled: a parent is drawn from the
3×3 patch around the cell with probability proportional to fitness, and its
copy replaces the cell's occupant. The copy may mutate. For a genome of `L`
bases and a per-base rate `u`, the number of mutations of each kind is drawn
from a binomial distribution, about `u × L` on average. There are seven kinds:
three *local* ones (a point mutation that flips one base, a small insertion
and a small deletion of 1 to 6 bases) and four *rearrangements* that act on a
segment between two random points (duplication to a random place, deletion,
translocation, inversion). Because the number of mutations grows with genome
length, a longer genome is also a more mutable one, which is the effect the
2007 study turned on.

### What Chemart records as a reaction network

Aevol has no kinetics, so Chemart records the run as an *observed* network,
one reaction per kind of event that happened, with its count. There are two
kinds of molecule: **genotypes** (`G<length>-<hash>`, the genome as the
structure) and **proteins** (`P<codons>-<hash>`, the coding bits and the
triangle). Expression is written `G -> G + P1 + ... + Pk`, the genome acting as
a catalyst of its own proteome. Replication is `P + V -> P + O` when the
parent `P` sits in a neighbouring cell and its offspring `O` replaces the
occupant `V`, or `P -> O` when a cell's own occupant is the parent. When the
parent and the occupant have the same genotype, the first form prints as
`2 G -> G + O`.

## Using it

The default run is small: a 3×3 grid of 5,000-base genomes for 30
generations, about 3 seconds, most of it spent drawing the ancestor. It is a
demonstration of the machinery, not one of the published experiments, which
used about 1,000 organisms for 20,000 generations. The first reactions listed
above are expression events: `G5000-db4e1670`, the ancestor, is born 19 times
and each time expresses one copy of the protein `P16-1d6da27a` (16 codons).
Genotypes whose genome carries a duplicated copy of that gene, like
`G8854-60289811`, express it twice.

The outcome is in `net.extras`:

```python
a = net.extras["analysis"]
a["ancestor"]    # {'genotype': 'G5000-db4e1670', 'fitness': 2.07e-65,
                 #  'metabolic_error': 0.14894, 'functional_proteins': 1}
a["best"]        # genotype 'G9280-4ed0c2ba', metabolic_error 0.14179,
                 # genome_length 9280, rnas 36, proteins 3, functional_proteins 3
a["mutations"]   # {'duplication': 7, 'deletion': 28, 'translocation': 30, 'inversion': 35,
                 #  'switch': 31, 'small_insertion': 21, 'small_deletion': 23}
net.extras["proteins"]["P16-1d6da27a"]
                 # {'m': 0.714, 'w': 0.0265, 'h': 0.669, 'codons': 16, 'functional': True}
```

`a["per_generation"]` has one row per generation with the best and mean
metabolic error and fitness, the mean genome length, the mean number of
functional proteins and the number of distinct genotypes; it is the place to
watch evolution happen. `a["best"]["phenotype"]` and
`a["phenotypic_target"]["points"]` are the two curves, ready to plot.
`net.extras["space"]["final_grid"]` gives the genotype in each cell, and
`genotype_parent` the parent of each genotype, so lineages can be traced.

**Watching adaptation.** A bigger grid for longer shows genes being duplicated
and tuned. With an 8×8 grid for 200 generations (about 18 s per run), the
best metabolic error fell from 0.149 to 0.090 on seed 1, from 0.151 to 0.135
on seed 2 and from 0.149 to 0.092 on seed 3, while the mean number of
functional proteins rose from 1 to between 8.6 and 12.3:

```python
net = chemart.generate_network("aevol", seed=1, grid_width=8, grid_height=8,
                               generations=200, min_genome_length=100)
```

These networks have about 5,500 species and 17,000 reactions, one genotype
per distinct genome ever born.

**Set `min_genome_length`.** With the default of 1 base, a large deletion can
leave a genome of a few bases, and Chemart then stops with
`ValueError: window shape cannot be larger than input array shape` (seed 2 on
a 4×4 grid for 60 generations does this). A floor such as
`min_genome_length=100` avoids it.

**Genome size.** At the default rates genomes grow: in the runs above the mean
length went from 5,000 to between 8,300 and 9,450 bases, pressed against the
cap `max_genome_length=10000`. The cap is Chemart's (aevol's is ten million),
set to keep runs and networks small; raise it to let genome size evolve
freely. To imitate the 2007 study, set all seven rates to one value `u`, as it
did, and compare values of `u`; see *Results* for why short runs cannot show
its effect.

**Other settings.** `target` takes your own list of `[height, mean, width]`
Gaussians; `selection_pressure=0` makes selection neutral; `w_max` bounds how
broad a protein can be. Runs slow down with grid area × generations × genome
length. A 32×32 grid, aevol's default, took 12.6 s for its first 10
generations, so a minute buys only a few dozen generations.

## Results

### Non-coding DNA and the mutation rate (Knibbe et al. 2007)

The founding study asked a question from molecular evolution: why do genomes,
especially eukaryotic ones, carry so much DNA that codes for nothing? The
common view was that such DNA is a passive by-product, accumulated through
insertion-biased mutation or self-copying elements. Knibbe, Coulon, Mazet,
Fayard and Beslon proposed that non-coding DNA could instead be *indirectly*
selected, because it changes how variable an organism's offspring are: in
aevol, intergenic DNA adds targets for duplications and large deletions,
which then copy or lose genes more often. They evolved 72 asexual populations
of 1,000 organisms for 20,000 generations, with six per-base mutation rates
from 5×10⁻⁶ to 2×10⁻⁴ (the same rate for every kind of mutation) and four
strengths of selection, starting from a 5,000-base genome with a single gene.

- **The mutation rate sets genome structure.** Gene number and the amount of
  non-coding DNA both reached an equilibrium that did not depend on the
  starting genome size but did depend on the mutation rate. At a low rate the
  typical final organism had 48,708 bases, 93 genes and 93% non-coding DNA;
  at a high rate it had 560 bases, 11 genes and 19% non-coding DNA, a compact
  genome with overlapping genes "much like the viral ones". Across runs,
  non-coding DNA reached up to 97% of the genome. Log-log regressions of gene
  number and of non-coding bases on the mutation rate were all significant
  (r² from 0.87 to 0.99). No cost of genome size and no mutational bias were
  applied: without selection, genomes lost all their genes and shrank below
  100 bases.
- **One neutral offspring.** The fraction of *neutral offspring*, the
  offspring that carry no mutation or only harmless ones, evolved to about the
  same value whatever the mutation rate, close to `1/W`, where `W` is the
  expected number of offspring of the best individual. In other words, genomes
  grew or shrank until the best organism left about one unchanged offspring
  per generation. Changing either the local-mutation rate or the
  rearrangement rate alone moved genome size in the same direction (their
  Table 1). The authors read this as an indirect selection of a trade-off
  between faithful transmission and exploration of new phenotypes.
- **Random genes are rare.** To start each population, random 5,000-base
  genomes were drawn until one had a beneficial gene. It took 610 draws on
  average, one functional gene per about 3 million random bases.

Chemart reproduces the last point: over 20 seeds it takes on average 670
draws (from 24 to 2,606) to find the starting genome. It does not reproduce
the main result. Those runs needed about 1,000 organisms over thousands of
generations to reach equilibrium, which is beyond a pure-Python port; in three
runs of 36 organisms for 300 generations with every rate at 2×10⁻⁴, the
genomes had only begun to change (mean length 3,300 to 5,500 bases). The
setups also differ: the 2007 paper used aevol version 4.5, with a 28-base
promoter consensus and rank-based selection over the whole population, while
Chemart ports version 9.4.0, with a 22-base consensus and fitness-proportionate
selection in a 3×3 neighbourhood. The paper reports that it repeated its
experiments with fitness-proportionate selection and found the same kind of
relationship.

### Homologous rearrangements and evolvability (Parsons et al. 2011)

The study the book cites asks what changes when rearrangements happen, as in
real cells, mostly between similar (*homologous*) sequences instead of at
random points. Parsons, Knibbe and Beslon added alignment-driven
rearrangements to aevol. The paper could not be retrieved for this page. Its
abstract, as listed on the Inria archive, reports that homologous
rearrangements are dangerous enough at high rates to select for short
genomes, yet in successful lineages their number correlates positively with
fitness gains, so they are both dangerous and needed for evolvability. The
book summarises it as showing that such rearrangements are important to the
evolvability of genetic regulatory networks. Chemart implements only the
random-breakpoint (non-homologous) rearrangements, so it cannot reproduce this
study.

### What the tests check

Chemart's decoding is checked against the hand-computed examples in aevol's
own unit tests: promoter positions and transcription levels, RNA lengths on
both strands, gene positions, the pooling of two transcripts of one gene, the
triangle and area values of the fuzzy sets, and the three-Gaussian target
(area 0.152483). The tests also check each mutation operator, the binomial
mutation counts, and that the network's replication events rebuild the final
grid exactly. For the catalog's listed phenomena:

- *A random genome rarely carries a gene*: confirmed by the runs above (about
  22 promoters, 0.1 genes per 5,000 bases); the tests check only that the
  initial genome beats a flat phenotype.
- *Proteins range from polyvalent to specialised, width bounded by `w_max`*:
  the folding tests check the extremes (`w = 0` and `w = w_max`).
- *Fitness improves as genes are added, tuned and duplicated*: a slow test
  checks, on two seeds of a 4×4 grid for 60 generations, that the best and
  mean metabolic errors fall, and that functional proteins and genome length
  grow. The best error does not fall at every generation: selection is not
  elitist, so the best genotype can be lost. In five of six runs (4×4 for 60
  generations and 8×8 for 200, seeds 1–3) it rose at least once.
- *Homologous and non-homologous rearrangements interact*: not implemented.
- *Indirect pressure on non-coding DNA*: not reproduced (see above).

Chemart also leaves out much of what later aevol versions offer: diploid and
sexual organisms, horizontal transfer, plasmids, changing or noisy
environments, stochastic gene expression and other selection schemes; the
implementation decisions list them all.

## Further reading

- Knibbe, C., Mazet, O., Chaudier, F., Fayard, J.-M. & Beslon, G. (2007).
  Evolutionary coupling between the deleteriousness of gene mutations and the
  amount of non-coding sequences. *Journal of Theoretical Biology* 244,
  621–630. The reference the 2007 study gives for the aevol platform.
- The aevol website, with documentation and the full list of publications:
  <https://www.aevol.fr/>
