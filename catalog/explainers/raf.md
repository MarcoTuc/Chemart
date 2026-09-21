## Introduction

RAF theory is a way of answering one question about any network of catalysed
reactions: does it contain a part that can keep itself going? Wim Hordijk and
Mike Steel introduced it in 2004, and developed it over the following decade
with Stuart Kauffman, Jotun Hein, Joshua Smith and others. It is not a
chemistry in the usual sense. It runs no simulation and invents no molecules.
It is an **analysis**: you hand it a reaction network and a list of raw
materials, and it tells you which reactions, if any, form a self-sustaining
whole.

The idea it makes precise is the **autocatalytic set** of Kauffman (1986): a
group of molecules in which every member is made by a reaction that some
other member catalyses (speeds up). No molecule copies itself, but the group
as a whole regenerates itself, as long as it is fed simple molecules from
outside, the **food set**. RAF stands for *reflexively autocatalytic and
F-generated*, the two conditions a set of reactions must meet:

- **reflexively autocatalytic (RA)**: every reaction in the set is catalysed
  by at least one molecule that the set itself can supply;
- **F-generated**: every molecule the reactions consume can be built from the
  food set F by using reactions of the set, step by step.

A simple picture: two reactions each join two food molecules into a longer
one, and each product catalyses the *other* reaction. Neither reaction can run
on its own, since its catalyst is the other's product; together, starting from
food, they sustain each other. That pair is an RAF.

Kauffman had argued that such sets arise almost inevitably once a chemistry is
rich enough, but his argument was criticised (Lifson 1997) because it seemed
to need each molecule to catalyse a number of reactions that grows
exponentially with molecule length, which is chemically unrealistic. RAF theory
turned the question into something that can be computed. Hordijk and Steel
gave an algorithm that finds the largest RAF in any network in polynomial
time, ran it on many random networks, and found that one to two
catalysed reactions per molecule are enough, a level that grows only linearly
with system size. The book presents this in §6.3.1, "Autocatalytic Sets",
under "Formalizing Autocatalytic Sets", in the chapter on the origin of life.

In the catalog, RAF sets sit next to three entries built on the same polymer
chemistry. [Kauffman's autocatalytic sets](kauffman-autocatalytic-sets.md)
builds the random reaction graph but does not search it; this entry builds the
same kind of graph and finds the self-sustaining sets in it.
[Bagley and Farmer's autocatalytic metabolism](bagley-farmer.md) adds rates and
a flow reactor, so it asks whether such sets actually take over the mass of the
system, which RAF theory, being purely structural, cannot say. The
[Jain-Krishna model](jain-krishna.md) keeps only the "who catalyses whom" graph
and lets it evolve.

## How it works

### The input: a catalytic reaction system and a food set

RAF theory works on a **catalytic reaction system** (CRS): a set of molecule
types X, a set of reactions R, and a record C of which molecules catalyse
which reactions. A reaction turns reactants into products; a catalyst is
needed for the reaction but is not used up. Alongside the CRS comes the food
set F, the molecules assumed to be freely available in the environment.

By default Chemart uses Kauffman's **binary polymer model**. Molecules are
strings of `0` and `1` up to a maximum length `n`. Reactions are ligations,
joining two strings end to end, and their reverse, cleavages, which cut a
string in two; as in the model papers, each ligation and its cleavage count as
one reversible reaction, `A + B <-> AB`. Each molecule catalyses each reaction
independently with a small probability p. The food set is every string up to
length `t = 2`: `0`, `1`, `00`, `01`, `10`, `11`.

The papers state the amount of catalysis as the **level of catalysis**
f = p × |R|, where |R| is the number of reactions: f is the average number of
reactions one molecule catalyses. It is the natural dial, because the
published results say how large f must be for an RAF to appear.

### Closure: what the food can build

The key operation is the **closure** of the food set relative to a set of
reactions R', written `cl_R'(F)`: start from the food, apply any reaction of
R' whose reactants are all present, add its products, and repeat until nothing
new appears. Catalysts are ignored at this step; closure only asks what *can*
be made from food with these reactions.

With closure the two conditions become precise. A non-empty set of reactions
R' is F-generated if every reactant of every reaction in R' lies in
`cl_R'(F)`, and reflexively autocatalytic if every reaction in R' has at least
one catalyst in `cl_R'(F)`. Chemart follows the 2011 version of the definition
(Hordijk, Kauffman and Steel), in which a catalyst may be any molecule in the
closure, food included. The 2004 version asked for a catalyst among the
molecules the reactions themselves consume or produce. The authors changed it
because early metabolic cycles were plausibly catalysed by simple molecules
that are available from the environment but take no part in the cycle.

### The algorithm: prune until nothing changes

To find RAFs, start with all reactions, compute the closure, and delete every
reaction that has a reactant, or all of its catalysts, outside it. Deleting
reactions can shrink the closure, so repeat. When a pass deletes nothing, what
remains is either empty, meaning the network has no RAF, or the **maximal
RAF** (maxRAF): the union of all RAFs in the network, which is unique. The
worst-case running time is O(|R|² log |R|), polynomial in the number of
reactions.

A maxRAF can often be cut down further. An **irreducible RAF** (irrRAF) is an
RAF that stops being one if any reaction is removed. There can be many of them
inside one maxRAF, and they are found by a random deletion search: go through
the maxRAF's reactions in random order and remove each one whose removal still
leaves a non-empty RAF.

### A worked example

Here is the example of Hordijk and Steel (2012), handed to Chemart as a
user-supplied system. The food set is `00`, `01`, `10`, `11`, and there are
four reversible ligations, each making a string of length four:

```
r1: 00 + 01 <-> 0001    catalysed by 1011
r2: 10 + 11 <-> 1011    catalysed by 0001
r3: 00 + 11 <-> 0011    no catalyst
r4: 01 + 10 <-> 0110    no catalyst
```

The closure of the food set under all four reactions is everything: all
reactants are food, so all four products can be made. The first pruning pass
therefore keeps `r1` and `r2`, whose catalysts `1011` and `0001` are in the
closure, and deletes `r3` and `r4`, which have no catalyst at all. Recomputing
the closure with only `r1` and `r2` still gives `0001` and `1011`, so the next
pass deletes nothing, and the algorithm stops:

```
True ['r1', 'r2'] ['00', '0001', '01', '10', '1011', '11']
[{'size': 2, 'reactions': ['r1', 'r2']}]
```

These are `raf_exists`, the maxRAF's reactions, the molecules they involve,
and the irrRAFs found. The maxRAF is also irreducible: on its own, `r1` has a
closure of `00`, `01`, `10`, `11` and `0001`, which lacks its catalyst `1011`,
so a one-reaction set is not an RAF. The network Chemart returns holds the
catalysed reactions only, each written with its catalyst on both sides:

```
00 + 01 + 1011 -> 0001 + 1011
0001 + 1011 -> 00 + 01 + 1011
10 + 11 + 0001 -> 1011 + 0001
1011 + 0001 -> 10 + 11 + 0001
```

To analyse your own network, pass it as `system` in the same form. The
formal specification below restates these definitions; its "reactor" is
empty, since there are no rates and no time, and food is simply assumed to be
present.

## Using it

The default call builds one random instance of the binary polymer model with
`n = 7` and `f = 1.5` and analyses it. There are 254 molecules and 1,284
reversible reactions, of which 334 happen to be catalysed, by 373 catalysis
events in all. Each event gives two network reactions (ligation and cleavage,
catalyst on both sides), hence the 746 reactions in the summary. In the
printed reactions, `2 0 + 001 -> 00 + 001` is the ligation of two `0` into
`00`, catalysed by `001`.

The analysis is in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["raf_exists"], a["max_raf_size"]           # (True, 258)
a["max_raf_closure_size"]                    # 198  molecules the maxRAF can make
a["irreducible_rafs"][0]["size"]             # 114
a["catalysis_probability"]                   # 0.0011682242990654205  p = 1.5 / 1284
```

`max_raf_reactions` lists the maxRAF's reactions by id, such as
`0+0001<->00001`, and `max_raf_molecules` the molecules they involve. The food
set is in `net.extras["food"]`, and again in `net.extras["buffered"]` because
it is taken to be always available. `net.extras["conservation"]` gives the
number of each letter as a conserved quantity, since ligation and cleavage
never create or destroy monomers. `reactions="all"` adds the uncatalysed
ligation/cleavage pairs to the network; the analysis is unchanged.

**Finding the transition.** Sweep f and count how often an RAF exists. At
`n = 8` over 20 seeds (about 6 seconds in total):

```python
for f in (1.0, 1.2, 1.4, 1.6, 2.0):
    runs = [chemart.generate_network("raf", seed=s, n=8, f=f, irreducible_rafs=0)
            .extras["analysis"] for s in range(20)]
    hits = [r["max_raf_size"] for r in runs if r["raf_exists"]]
    print(f, len(hits) / len(runs), round(sum(hits) / len(hits)) if hits else 0)
```

```
1.0 0.15 1
1.2 0.3 311
1.4 0.9 493
1.6 1.0 601
2.0 1.0 800
```

The fraction with an RAF jumps from 0.3 to 0.9 between f = 1.2 and 1.4, close
to the published half-way point of 1.0970 + 0.0189 × 8 ≈ 1.25. The three hits
at f = 1.0 are all the same trivial one-reaction RAF, `1+0<->10`. Its reactants
and product are food, and in seed 2 its catalyst is the food molecule `11`,
which the 2011 definition allows. The papers mention this kind of triviality
(Hordijk, Kauffman and Steel 2011) and note that it is easy to exclude
afterwards: check whether the maxRAF contains a reaction that does not have
all its catalysts in the food set.

**Sampling irreducible RAFs.** `irreducible_rafs` sets how many irrRAFs to
sample. Five samples on the default network give five different irrRAFs, of
sizes 114, 101, 101, 109 and 117, inside a maxRAF of 258.

**Speed.** Finding the maxRAF is fast: `n = 10` (16,388 reactions) takes about
2 seconds and `n = 12` (81,924) under 2 seconds when no RAF is found. The
irrRAF search reruns the algorithm once per maxRAF reaction, so one irrRAF at
`n = 10` adds about 9 seconds; set `irreducible_rafs=0` for sweeps.

## Results

**A polynomial-time test.** Hordijk and Steel (2004) gave the algorithm and
proved it returns the maxRAF, or nothing when there is no RAF, in worst-case
time O(|R|² log |R|); in practice it ran in sub-quadratic time on random
instances of the polymer model. They applied it to networks of about five
million catalysed reactions. The 2011 paper simplified the algorithm to the
single pruning step described above, with the same running time. Chemart's
tests check the algorithm on small networks with known answers, including the
worked example and networks that fail only one of the two conditions, and
check that the maxRAF is a fixed point. They do not measure running time.

**Linear, not exponential, catalysis.** The central result. Applying the
algorithm to many random instances of the polymer model, Hordijk and Steel
(2004) found that the level of catalysis needed for an RAF to appear grows
only linearly with `n`, and Mossel and Steel (2005) confirmed this
analytically. Hordijk, Hein and Steel (2010) put the required level at
"between 1 and 2 reactions per molecule (on average)", which they call
"(bio)chemically quite realistic". Hordijk, Kauffman and Steel (2011) fitted
the value of f at which half the instances contain an RAF, over 100 to 1,000
instances for each n from 7 to 20: f(n) = 1.0970 + 0.0189 n (their Table 1,
case A), which gives only 1.475 at `n = 20`. The theoretical bound, for RAFs
that use every molecule, has slope 1.6339 and gives 32.678 at `n = 20`. Chemart reproduces the result: its tests find the probability of an RAF near one half on the
Table 1 line for `n` = 7, 8 and 9, find it near one half at the catalysis
probabilities that Hordijk, Smith and Steel (2015) give for `n` = 8 and 10, and
find that for every `n` from 7 to 10 (an eightfold growth in molecules) f = 1
is below the transition and f = 2 always gives an RAF.

**A sharp transition.** Steel, Hordijk and Smith (2013) computed RAF sizes at
`n = 10` over 1,000 instances per value of f. Below f = 1.20 they found no RAF
at all. Just above it RAFs appeared in 6 of 1,000 instances and became more
frequent as f rose. The first ones were already large, 1,222 reactions on
average, out of |R| = 16,388. Chemart's tests reproduce this: no RAF at f = 1.15,
RAFs of roughly the published size at f = 1.25, and larger, more frequent ones
at f = 1.45, with the test expecting about 2,000 reactions there.

**Big maxRAFs, many irreducible RAFs.** In the same figure the maxRAF grows
roughly linearly with f while the size of one irrRAF stays about constant, 624
reactions near f = 1.20. A maxRAF can contain exponentially many irrRAFs, as
Hordijk, Steel and Kauffman (2012) proved, and finding a smallest RAF is, in
general, NP-hard (Steel, Hordijk and Smith 2013). Hordijk, Smith and Steel
(2015) tested whether this happens in practice. At `n = 8` with maxRAFs of 375 reactions on average, samples of 10,000 irrRAFs
were all different, from which they conclude, with 99% confidence, that at
least hundreds of millions of irrRAFs exist; two irrRAFs share about half their
reactions on average (between 25% and 80%). Chemart's tests check that the
maxRAF grows with f while the irrRAF does not, that one irrRAF at `n = 10` is
of roughly the published size, that the search returns genuinely irreducible
sets, and that different runs find different ones. The overlap statistics and
the lower bound on the number of irrRAFs are not reproduced.

**Extensions Chemart does not implement.** Mossel and Steel (2005) added
inhibition: a molecule may block a reaction. They showed that deciding whether
there is an RAF none of whose reactions is inhibited by one of its own
molecules is NP-complete (the book says NP-hard), so the polynomial algorithm
cannot be extended to it; Hordijk and Steel (2012) showed
the problem is tractable when the number of inhibitors is small. The same 2012
paper simulated the 4-reaction example with the Gillespie algorithm. With
ligation only, the products of the RAF pair grew exponentially and those of
the two uncatalysed reactions linearly; with cleavage as well, the RAF pair
held its products at a much higher level, unless ligation was made too slow
relative to cleavage. Chemart has neither
inhibition nor dynamics. As the book reports, the authors themselves list the framework's limits: it
ignores dynamics, compartments, heredity and natural selection. An RAF is a
necessary condition for a self-sustaining chemistry, not a sufficient one.

## Further reading

- Kauffman, S. A. (1986). Autocatalytic sets of proteins. *Journal of
  Theoretical Biology* 119, 1–24. The model RAF theory formalises.
- Lifson, S. (1997). On the crucial stages in the origin of animate matter.
  *Journal of Molecular Evolution* 44, 1–8. The criticism that Kauffman's
  argument needs exponentially growing catalysis.
- Hordijk, W., Steel, M. & Kauffman, S. (2012). The structure of autocatalytic
  sets: evolvability, enablement, and emergence. *Acta Biotheoretica* 60(4),
  379–392. https://doi.org/10.1007/s10441-012-9165-1. The proof that a maxRAF
  can hold exponentially many irrRAFs.
- Mossel, E. & Steel, M. (2005). Random biochemical networks: the probability
  of self-sustaining autocatalysis. *Journal of Theoretical Biology* 233(3),
  327–336. https://doi.org/10.1016/j.jtbi.2004.10.011. The analytical linear
  bound and the NP-completeness of RAFs with inhibition (book reference [599]).
