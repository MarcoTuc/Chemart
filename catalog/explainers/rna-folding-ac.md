## Introduction

This chemistry asks how a metabolism could first come under genetic control,
in an early "RNA world" where RNA molecules were both the genes and the
enzymes. Alexander Ullrich and Christoph Flamm built it (Ullrich and Flamm,
2008) as the catalytic layer of a larger simulation of metabolic evolution. Its
molecules are RNA sequences, strings over the four bases A, U, G and C. Some of
them act as ribozymes, RNA enzymes, and the model decides *which* reaction a
ribozyme catalyses by looking at its shape.

The shape comes from real physics. An RNA strand folds back on itself: stretches
of complementary bases pair up into double helices (stems), and the unpaired
stretches between them form loops. The pattern of pairs is the strand's
*secondary structure*, and standard software predicts it by finding the
structure of minimum free energy. The model then reads the reaction off that
structure: the longest loop is taken as a picture of the transition state of
the reaction the ribozyme speeds up. So there are two steps, sequence →
structure → function, the genotype → phenotype map of biology in miniature. The
first step is done by thermodynamics, the second by an explicitly artificial
rule. Flamm et al. (2010) say of it: "Without any claim of physical realism, we
interpret this cycle as an encoding of the imaginary transition state of the
catalyzed reaction."

Because the fold does the work, the chemistry inherits a well-known property of
RNA: many different sequences fold into the same structure. Many mutations are
therefore *neutral*: they change the gene but not the enzyme. Ullrich and Flamm
studied this neutrality, and then used the map inside their evolving protocells.

Banzhaf and Yamamoto describe it in the chapter on modelling biological
systems (book §18.1.1, "RNA Folding"). It is a simulation model and, in
Chemart, a generator of reaction networks. Its nearest neighbours are the
[matrix chemistry](matrix-chemistry.md), which uses folding only as a metaphor
(a binary string "folds" into a matrix that acts on other strings), and
[ToyChem](toychem.md), the graph-based chemistry whose molecules the original
ribozymes act on. Here the folding is the real, energy-based one used by RNA
biologists. Chemart keeps the folding and the structure-to-function map but,
lacking ToyChem's metabolite chemistry, lets the ribozymes cut and join RNA
strands instead.

## How it works

### From sequence to structure

Each molecule is an RNA sequence, 30 bases long by default. Chemart folds it
with the ViennaRNA package (version 2.7.2), which computes the minimum free
energy (MFE) structure by dynamic programming. The structure is written in
*dot-bracket notation*: a string as long as the sequence, where a matching pair
of brackets marks two paired bases and a dot marks an unpaired one. For
example `GGGGAAAACCCC` folds into `((((....))))`, at −5.40 kcal/mol: four G–C
pairs form a stem and the four A's close a *hairpin* loop at its end.

A structure is a set of loops. A **hairpin** loop is closed by one stem. An
**interior** loop (or a **bulge**, when one side has no unpaired bases) sits
between two stems. A **multiloop** is where three or more stems meet. The
**size** of a loop, in Chemart, is the length of the ring of bases around it:
its unpaired bases plus the two bases of each stem pair that closes it. The
hairpin above has size 4 + 2 = 6. The open ends of the strand, outside every
stem, form no ring and do not count.

### From structure to function

The second step borrows an idea from the classification of organic reactions.
Fujita (1986) superimposed the molecules before and after a reaction and kept
only the bonds that change. What remains is usually a single ring, the
*imaginary transition structure* (ITS), in which bonds that break and bonds
that form alternate. Flamm et al. (2010) note that such single-ring ITSs
describe over 90% of known reactions. They organise ITSs in a tree: level 1 is
the ring size (how many atoms take part), level 3 the arrangement of the
electron reshuffling, level 4 the atom types.

The ribozyme map reads each level off the longest loop. Its size is the ITS
size. Its stems give the arrangement (Chemart uses only how many there are).
Its unpaired bases give the "atom types". Chemart turns this into a concrete
rule:

- The loop must be **even-sized**, because a ring of 2n positions can alternate
  n broken and n formed bonds only if the count is even. Odd loops are inert.
- Its size must lie between `its_min` and `its_max` (4 and 12 by default).
- A hairpin (one stem) makes a **cleavase**, which cuts a strand in two:
  `E + X → E + Y + Z`. An interior or bulge loop (two stems) makes a
  **ligase**, which joins two strands: `E + Y + Z → E + X`. Multiloops are
  composite transition structures and are inert. If two loops tie for
  longest, the one with more stems counts.
- The ribozyme recognises its target by base pairing. It takes the first three
  bases of its loop (`min_recognition`), and its **site** is their Watson–Crick
  reverse complement, the word that would pair with them. A cleavase cuts a
  substrate just after the first occurrence of the site; a ligase joins two
  strands whose junction spells it.

The first two rules follow Fujita and Flamm et al. The choice of cleavage and
ligation as the reactions is Chemart's own: the published model rewrites ToyChem
molecule graphs, and the paper that specifies which ITS stands for which
reaction could not be obtained (see the implementation decisions). Two
consequences: the ribozyme itself is untouched, so it is a true catalyst; and
cutting or joining strands never creates or destroys a base, so the total
number of nucleotides is conserved.

### A worked example

Here is the first reaction of the default network. The catalyst is a 30-base
strand; Chemart folded it as

```
UACGCUUUCUAGCAGUUAUUCAUUCAACUC
...(((....)))((((........)))).     -2.4 kcal/mol
```

It has two hairpins. The first encloses 4 unpaired bases (size 6), the second
8 (size 10). The longest loop has size 10: even, inside 4–12, closed by one
stem. So the molecule is a cleavase. Its loop reads `AUUCAUUC`; the first three
bases are `AUU`, and the word that pairs with them is `AAU`. The reaction it
catalyses:

```
UACGCUUUCUAGCAGUUAUUCAUUCAACUC + CGUUAAUUACUCCUCCGGAAUUUGUCCUAC
  -> UACGCUUUCUAGCAGUUAUUCAUUCAACUC + CGUUAAU + UACUCCUCCGGAAUUUGUCCUAC
```

The substrate contains `AAU` at bases 5–7 (`CGUU·AAU·UAC…`), so the cleavase
cuts right after it, into a 7-base and a 23-base fragment. The fragments are
new species. Each is folded in turn and may itself be a catalyst or a
substrate. That is what makes the chemistry *constructive*: the set of species
grows as reactions run.

Most strands are not catalysts. In 1,000 random 30-base sequences, 568 were
inert (their longest loop was odd, too big, a multiloop, or too short to carry
a site), 338 were cleavases and 94 ligases.

### Two reactors

Chemart offers two ways to turn the rule into a network. In **closure** mode
it starts from a pool of random sequences and applies every ribozyme to every
possible substrate (one for a cleavage, two for a ligation), folds the
products, and repeats until nothing new appears or the species budget
`max_species` is reached. In **well-stirred** mode it simulates a pot of
molecules: at each of `steps` collisions three molecules meet at random, the
first acts as the catalyst if it is one, and at most one reaction fires. In a
cleavage the third molecule rides along unchanged. Neither mode has rate
constants: no published rates exist that Chemart could reproduce.

## Using it

The default run is closure mode on 8 random 30-base sequences. It is not a
published experiment. It completes with 21 species and 8 reactions, all of them
cleavages. `net.extras` explains the network:

```python
net.extras["functions"]["UACGCUUUCUAGCAGUUAUUCAUUCAACUC"]
# {'its_size': 10, 'stems': 1, 'loop_sequence': 'AUUCAUUC',
#  'site': 'AAU', 'reaction': 'cleavage'}
net.extras["analysis"]["catalytic_species"]   # 4
net.extras["analysis"]["cleavases"], net.extras["analysis"]["ligases"]   # (3, 1)
net.extras["energies"]["UACGCUUUCUAGCAGUUAUUCAUUCAACUC"]   # -2.4
```

`functions` lists only the catalytic species. Each species' `structure` field
holds its sequence and dot-bracket fold, separated by a space, and
`extras["conservation"]` gives the nucleotide count of every species (the
conserved quantity). The one ligase, `CGCACGCUCGUUCAGGUCCACGUUAGUCCU`, has a
hairpin and an interior loop both of size 10, so the tie rule makes it a ligase.
Its site is `AGC`, and no two strands in this default pool meet at a junction
that spells it, so no ligation appears.

**A bigger pool.** With 12 seed sequences the closure no longer finishes: it
stops at the 100-species budget, with 300 reactions, 4 of them ligations
(about a second):

```python
net = chemart.generate_network("rna-folding-ac", seed=1, pool=12, max_species=100)
# rna-folding-ac: 100 species, 300 reactions, status=truncated
```

One of those ligations joins two 30-mers into a 60-mer, because the junction
`…UUGUA|GCUAUG…` spells the ligase's site `AGC`:

```
CGCACGCUCGUUCAGGUCCACGUUAGUCCU + GCUAUGCGCUUCCAGGUUUUUAACCUUCGG + GUGGCUUGCGGAACGACAUGCUUCUUUGUA
  -> CGCACGCUCGUUCAGGUCCACGUUAGUCCU + GUGGCUUGCGGAACGACAUGCUUCUUUGUAGCUAUGCGCUUCCAGGUUUUUAACCUUCGG
```

**A pot of molecules.** `mode="well-stirred"` records only the reactions that
actually fired, each with its count (`r.count`):

```python
run = chemart.generate_network("rna-folding-ac", seed=1, mode="well-stirred",
                               pool=12, steps=400)
# rna-folding-ac: 51 species, 20 reactions, status=observed
```

In this run each of the 20 reactions fired once. `dilution="constant"` removes
molecules after every reaction to keep the population at its starting size.
Well-stirred mode needs `pool` of at least 3.

Sequence length is `seq_length`; Flamm et al. (2010) used 100-base genes, and
100-base sequences fold and map without trouble, but a closure over many of
them grows fast, so keep `max_species` modest. The ITS window and the length of
the recognition word are in the parameter table below.

## Results

**The map itself.** Ullrich and Flamm (2008) introduced the sequence →
structure → function map at the CMSB conference, inside a model of evolving
protocells whose metabolism runs on ToyChem molecules. Flamm et al. (2010)
restate it in an open-access paper, which is the specification Chemart follows.
To check that the artificial second step does not dominate, they compared the
autocorrelation of the sequence → structure map with that of the full
sequence → function map, for 100-base sequences (1,000 random reference
sequences, 1,000 mutants per mutation distance). The two curves behave much
alike, which they take to justify using an ad hoc structure-to-function rule:
when folding is dominated by neutral, essentially random structures, the
second step has little influence. Chemart reproduces the map (sizes, stem
counts, the even-size rule, the hairpin/interior split) in its tests, but not
this autocorrelation study.

**Neutral networks and nearby functions.** Ullrich and Flamm's ECAL 2009 paper
(published 2011) studied the neutrality of the map through random walks in
sequence space. The book summarises the finding: a large number of neutral
variants, and a variety of alternative functions within a short distance of
any sequence. The system is therefore robust to harmful mutations and
evolvable at the same time. Their abstract reports that, compared with
genotype–phenotype maps built on cellular automata, random Boolean networks
and other RNA-folding maps, theirs gave "the highest extent, connectivity and
evolvability of the underlying neutral network". Chemart's test takes 20
random 30-base sequences and all 90 single-base mutants of each: about 40% of
mutants keep the exact structure (725 of 1,800 in the run behind the test), and
each sequence has between 3 and 7 distinct function classes (reaction type and
ITS size, or inert) one mutation away. This is Chemart's measurement, not a
number from the papers, which studied 100-base sequences with random walks.

**Book table 18.1.** The book prints two random 60-base sequences with the
structures the Vienna RNAfold server gave. Chemart reproduces the second
exactly, at −7.10 kcal/mol. For the first, current ViennaRNA finds a different
structure at −10.90 kcal/mol; the book's structure scores −9.80 and lies within
1.5 kcal/mol of the optimum, so it is a suboptimal fold under today's energy
parameters. The tests record both facts.

**Hub metabolites.** In long evolution runs (1,000 generations) the metabolic
networks came to resemble real ones: Flamm et al. (2010) report that their
node-degree distribution follows a power law, so a few metabolites are hubs,
citing Ullrich and Flamm (2008). The book describes the same result as
small-world connectivity with highly connected metabolites. Chemart does not
reproduce it. It has no metabolites, genomes or selection, only ribozymes
acting on RNA.

**Phases of metabolic evolution.** Later evolution experiments, which the book
cites to Flamm et al. (2010) and to two papers by Ullrich et al. (2010, 2011),
split evolution into phases. In the book's summary: pathways first
grow forward, as enzymes appear that extract energy by breaking food
molecules into ever simpler parts; later, existing enzymes are recruited into
new reactions and specialise, forming complex pathways. Flamm et al. (2010)
show two 100-generation runs with 5,000-base genomes carrying 100-base genes,
one starting from ten cells fed five small molecules, the other from one cell
fed glucose. In both, enzymes and metabolites that appeared early take part in
more reactions, and older reactions carry more flux, which they note might
support the "patchwork" hypothesis of enzyme recruitment. Chemart does not
reproduce these experiments, for the same reason: they need ToyChem's
chemistry, a genome and a fitness function, which are outside this generator.

## Further reading

- Ullrich, A., Flamm, C., Rohrschneider, M. & Stadler, P. F. (2010). In silico
  evolution of early metabolism. In Fellermann et al. (eds.), *Artificial
  Life XII*, 57–64. MIT Press. (Book reference [875].)
- Ullrich, A., Rohrschneider, M., Scheuermann, G., Stadler, P. F. & Flamm, C.
  (2011). In silico evolution of early metabolism. *Artificial Life* 17(2),
  87–108. (Book reference [876].)
- ViennaRNA, the folding package used here: <https://www.tbi.univie.ac.at/RNA/>
