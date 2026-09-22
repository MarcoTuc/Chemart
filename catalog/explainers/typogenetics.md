## Introduction

Typogenetics is a toy genetics made of typed letters. Douglas Hofstadter
introduced it in *Gödel, Escher, Bach* (1979), in Banzhaf and Yamamoto's words
"as a way to exemplify the formal characteristics of life". Its molecules are
strands, strings over the four letters A, C, G and T. The same strand plays two
parts. Read two letters at a time, it is a *gene*: it spells out a short program
of operations such as "copy", "cut", "move right" or "insert a T", which
Hofstadter calls a typographic enzyme (*typoenzyme*). Read as a string, it is
the *substrate* those enzymes work on. A strand's enzymes can act on the strand
that coded them, so a strand decides what happens to itself.

That self-reference is the point. Hofstadter closed his description with a
puzzle: write a strand whose enzymes, applied to it, produce two copies of it.
Such a strand is a self-replicator in the most literal sense, a description
that builds a copy of itself. The puzzle has been solved, and the
solutions constructed by Morris and by Snare share a shape: the strand is an *inverted repeat*, its second half the
complementary mirror image of its first, so that making the complementary
strand yields the original again. `CGATTCGAATCG` is one.

Hofstadter's rules leave gaps: for instance, he does not say in which order
a strand's several enzymes act. Others filled them in different ways: Morris (1989) studied the system's logic, Varetto
(1993, 1998) ran a modified version in simulation, and Snare (1999, a Monash
honours thesis) wrote a specification that reconciles the three and classified
what strands do when applied to themselves repeatedly. Chemart follows Snare,
checked against the book's worked figure and against Hofstadter's own example.

Typogenetics is a simulation model of the machine-and-tape kind, and the book
places it in chapter 10 on automata and machines (§10.5.2), between
[Laing's molecular machines](laing-molecular-machines.md) and
[polymers as Turing machines](mccaskill-polymer-tm.md). Compared with Laing's
machines it adds a genetic code, so the machine is *translated* from a tape
instead of being a different form of string; compared with
[the machine-tape chemistry](ikegami-hashimoto.md) of Ikegami and Hashimoto it
has only one kind of molecule. The book also names it, in §11.2, as one of the
early sources of lock-and-key binding in artificial chemistries, here based on
DNA base pairing, and says the [automata reaction](automata-reaction.md)
chemistry was inspired by it. Its variable-length molecules, variable number of
products and free choice of binding site make it a demanding test for any tool
that enumerates a chemistry's reactions.

## How it works

### Strands and the genetic code

The letters are called *bases*, and a base's position on a strand is a *unit*.
Bases pair up as in DNA: A with T, and C with G. A and G are *purines*, C and T
*pyrimidines*; some operations look for one class or the other.

To translate a strand, cut it into consecutive pairs of letters (*duplets*),
ignoring an odd last base, and look each duplet up in the code table. Each
duplet codes one "amino acid", which is really an operation. The duplet `AA`
codes nothing: it is punctuation that ends one enzyme and starts the next, so a
strand can code several enzymes. Hofstadter's table, as printed in book table
10.1, with each operation's *kink* (s, l or r) after the colon:

```
          A        C        G        T      <- second base
A         -      cut:s    del:s    swi:r
C       mvr:s    mvl:s    cop:r    off:l
G       ina:s    inc:r    ing:r    int:l
T       rpy:r    rpu:l    lpy:l    lpu:l
```

The fifteen operations act on the unit the enzyme is bound to:

- `cut` cuts the strand just right of the bound unit;
- `del` deletes the bound base and moves one unit right;
- `mvr`, `mvl` move one unit right or left;
- `cop` turns *copy mode* on, `off` turns it off. In copy mode, every unit the
  enzyme stands on gets its complementary base written in a second row
  opposite it, so the enzyme builds a partial double strand as it moves;
- `ina`, `inc`, `ing`, `int` insert A, C, G or T just right of the bound unit
  and move onto it;
- `swi` switches the enzyme to the other row of the double strand;
- `rpy`, `rpu` search rightwards for the nearest pyrimidine or purine, and
  `lpy`, `lpu` leftwards.

The second row reads in the opposite direction, so for an enzyme on it, right
and left are swapped. An enzyme stops when it moves off the end of a strand or
into a gap.

### Where an enzyme binds: the fold

An enzyme starts at a unit holding its *preferred base*, and that preference
comes from its shape. Draw the enzyme as a chain of links; each amino acid
turns the chain left (l), right (r) or not at all (s). The direction of the
last link relative to the first gives the base: straight on means A, turned
left C, turned right G, turned back T. Hofstadter counts the turns of every
amino acid except the first and the last. Morris and Varetto count them all,
which gives a different preference for some enzymes; Chemart offers both.

If the strand holds the preferred base at several units, Hofstadter lets the
enzyme bind any of them, and different choices give different results. Varetto
always takes the rightmost.

### A worked example: the default self-replicator

The default run starts from Snare's replicator `CGATTCGAATCG`. Its duplets
`CG AT TC GA AT CG` translate to one enzyme, `cop-swi-rpu-ina-swi-cop`. The
turns of the middle four amino acids are r, l, s, r: one net right turn, so
the enzyme prefers G. Chemart writes this as the species' structure,
`cop-swi-rpu-ina-swi-cop:G`.

The strand has G at units 1, 6 and 11 (counting from 0). Bound to the last one,
the enzyme runs as follows (a real trace from Chemart's machine; `.` is an
empty position):

```
bind   upper ............
       lower CGATTCGAATCG      enzyme on the lower row at unit 11 (G)
cop    upper ...........C      copy mode on: C written opposite the G
swi    upper ...........C      the enzyme moves to the upper row
rpu    upper GCTAAGCTTAGC      search for a purine, filling in as it goes
       lower CGATTCGAATCG      ...and runs off the end: the enzyme stops
```

The search moves leftwards along the lower strand (rightwards for the upper
row), and bases it writes in copy mode on the way do not end the search, so it
fills in the whole complementary row and falls off the end. The remaining
amino acids never run. The two rows then separate into single strands, the
upper one read in its own direction: `CGATTCGAATCG` again, because the strand
is an inverted repeat. One strand has become two.

Bound at unit 6 or unit 1 instead, the same enzyme copies only part of the
strand, and the products are `CGATTCGAATCG + CGAATCG` or `CGATTCGAATCG + CG`.
These three outcomes are the first three reactions of the default network,
printed below under *Using it in Chemart*. The fourth comes from the daughter `CGAATCG`, whose two enzymes `cop:A`
and `rpu:A` only manage to write one T opposite an A.

### Reactions and the reactor

A reaction is the whole episode: translate a strand, run its enzymes in order,
and let the double strand fall apart into *daughters*. When a strand codes
several enzymes, each binds in the stretch of strand where the previous one
stopped (Snare's rule). Double strands exist only during a reaction; the
species are always single strands. In the book's system, and in Chemart's
default `self` mode, a strand's enzymes act on the strand itself, so the
reaction is `s -> daughters`. A self-replicator is `s -> 2 s`, which is why
the default network lists `CGATTCGAATCG` as a catalyst. Chemart's `pair` mode
adds a variant in which the enzymes of one strand (the gene, which survives)
act on another.

A strand whose enzymes find nowhere to bind, or change nothing, has no
reaction. Snare calls the first kind a *dud*: `CGGC` codes `cop-inc`, which
prefers A, and there is no A in `CGGC`.

## Using it

Typogenetics has two faces. `chemart.generate_network`, printed above,
returns the *closure* of the seed strands; `chemart.evolve` runs a
*population* of strands and returns a trajectory (see *A population* below).

The default call enumerates every reaction reachable from the seed strand:
it applies each new daughter to itself in turn until nothing new appears. With `binding_tiebreak="all"`, every binding choice gives
its own reaction, so one strand can have several. The default seed closes
after four species. Each species' `structure` gives the enzymes it codes and
their preferred bases:

```python
for s in net.species: print(s.id, "|", s.structure)
# CGATTCGAATCG | cop-swi-rpu-ina-swi-cop:G
# CG | cop:A
# CGAATCG | cop:A; rpu:A
# T | None
```

`net.status` is `complete` when the closure ended by itself and `truncated`
when `max_species`, `max_length` or `max_branches` cut it off. Most strands
grow without limit through insertions, so truncation is normal: Hofstadter's
example strand `TAGATCCAGTCCATCGA` hits the 200-species budget (275 reactions,
about 5 s), and so do the four random 12-base strands that `strands=[]` draws
with `seed=1`.

**Book figure 10.6.** The book's example reaction, from Morris, needs Morris's
fold. With free binding it has four outcomes, the last being the figure's:

```python
net = chemart.generate_network("typogenetics", seed=1,
                               strands=["CGACCCAACGATTTTTCAT"], fold="morris")
# first reactions:
# CGACCCAACGATTTTTCAT -> ACCCAACGATTTTTCAT + 2 CG
# CGACCCAACGATTTTTCAT -> ATTTTTCAT + CGACCCAACG + CGT
# CGACCCAACGATTTTTCAT -> ATTTTTCAT + CGACCCAACG + CGTT
# CGACCCAACGATTTTTCAT -> ATTTTTCAT + CGACCCAACG + CGTTGGGT
```

**A population.** `chemart.evolve` puts `copies` of each seed strand in a pot
and draws a strand `steps` times; after each productive reaction a random
molecule is removed, so the population stays the same size. A draw has one
outcome, so `binding_tiebreak` picks one binding site: `rightmost` (Varetto's
rule), `leftmost` or `random`; the default `all` is read as `random` here. The
trajectory has a frame per generation, as many draws as the pot holds strands.
With the rightmost rule the replicator takes over from the dud:

```python
traj = chemart.evolve("typogenetics", seed=3,
        strands=["CGATTCGAATCG", "CGGC"], binding_tiebreak="rightmost",
        copies=50, steps=3000)
[f.state.get("CGATTCGAATCG", 0) for f in traj.frames][:7]
# [50.0, 74.0, 87.0, 95.0, 96.0, 99.0, 100.0]
net = traj.network
net.extras["final_state"]   # {'CGATTCGAATCG': 100}
# CGATTCGAATCG -> 2 CGATTCGAATCG  (x2910)
```

The frames come every 100 draws, and the replicator fills the pot after 600.
With `binding_tiebreak="random"` (or the default) and the same seed, the
replicator binds its other two G sites two times out of three, and the final
pot of 100 is `{'T': 62, 'CG': 17, 'CGATTCGAATCG': 12, 'CGAATCG': 9}`:
`CGAATCG` keeps making `T`, which codes no enzyme and accumulates. Both runs
take under two seconds.

**Other variants.** `code_table="varetto"` uses Varetto's code, which moves the
four insert operations within row G; `reaction="pair"` lets every strand's
enzymes act on every strand, itself included, which makes far more reactions:
seeding it with the replicator and the dud (with `binding_tiebreak="rightmost"`
and `max_species=30`) gives 242 reactions among 30 species.

## Results

**Hofstadter's examples.** *Gödel, Escher, Bach* works one example by hand:
the enzyme `rpu-inc-cop-mvr-mvl-swi-lpu-int`, bound to the middle G of
`TAGATCCAGTCCATCGA`, leaves `ATG` and `TAGATCCAGTCCACATCGA`. It also states that
the enzyme `rpy-ina-rpu-mvr-int-mvl-cut-swi-cop` binds C, which is exactly the
enzyme the longer product `TAGATCCAGTCCACATCGA` codes. GEB itself was not
accessible, so Chemart takes both examples as independent implementations
reproduce them, and its tests check all three facts. Two details had to be
settled to get there: the worked example's binding site is taken as given,
since neither kink-counting rule gives G for that enzyme, and an insertion
moves the enzyme onto the new base, where Snare's code would give different
products. Morris's fold, which counts every kink, gives G instead of C for
that enzyme (Snare 2.4), so Chemart's default is Hofstadter's fold.

**Morris: a logic of typogenetics.** Morris (1989) treated typogenetics as a
logic of artificial life. According to Snare, he gave a 64-base self-replicator
built as an inverted repeat, coined the word "dud", and argued that no strand
can completely erase itself. His version differs from Hofstadter's: it counts
every kink and its deletion closes the gap. The book's figure 10.6 comes from
him: two enzymes coded by `CGACCCAACGATTTTTCAT` produce `CGACCCAACG`,
`ATTTTTCAT` and `CGTTGGGT`. Chemart's tests reproduce this reaction step by step
with Morris's fold, and show it is one of the strand's possible outcomes.

**Varetto: replicators and tanglecycles.** Varetto (1993) simulated a modified
system with his own code table and rightmost binding. The book reports that he
found short self-replicators whose numbers grow exponentially, and autocatalytic
cycles of strands he called "tanglecycles". A later paper (1998), according
to Snare, studied strands under limited resources, and the book mentions
variants with a lethal factor. Snare could not verify Varetto's example
replicators (`GC`, `GGC`, `GTGC`): applying Varetto's stated rules to `GGC` gives a different daughter
from the one in his results. Chemart does not reproduce Varetto's
simulations. It offers his code table (from Snare's copy of it) and his
binding rule, but his papers, and so his lethal factor, were not accessible.

**Snare: classes of strands.** Snare (1999) followed the daughters of each
strand through repeated self-application and drew them as a graph. A *dud* has
no daughters. A *self-perpetuator* lies on a cycle, so it keeps coming back but
never in two copies; a *self-replicator* lies on more than one cycle. He
constructed replicators as inverted repeats wrapped around a copying enzyme:
`CGATTAATTAATCG` and `CGATTCGAATCG`, whose enzyme binds the rightmost G, copies
the strand and runs off its end. Any inverted repeat placed inside the second
wrapper gives another replicator, so there are infinitely many. Chemart's
tests check that both strands give `s -> 2 s` and that the dud `CGGC` has no
reaction and codes `cop-inc:A`.

Snare's systematic search found no replicators. His largest graph had 90,428
strands and 102,323 edges and contained 25,335 cycles, but only 66 were
longer than one generation, and the longest took three:
`CCGGA -> CCGGGA -> CCGGAGGA -> CCGGA`. Nor could he show that typogenetics
can compute everything a Turing machine can: he built some of the functions
of primitive recursion but not all. Chemart does not repeat the search, but its
closure builds the same kind of daughter graph from a seed. Kvasnička,
Pospíchal and Kaláb (2001) later studied replicators and hypercycles in
typogenetics; their paper was not accessible for this page.

**Several products, several outcomes.** The book stresses that molecules vary
in size, a reaction can give a variable number of products, and binding is
nondeterministic. The default network shows all three: one strand, three
alternative reactions with two products each. The tests check that the
default closure lists the three outcomes.

**Beyond the book's system.** The population run and the `pair` reaction are
Chemart additions, with no published counterpart. The tests check that in the
population the replicator ends above 90 of 100 molecules against the dud, and that in
`pair` the gene strand survives every reaction.

## Further reading

- Hofstadter, D. R. (1979). *Gödel, Escher, Bach: An Eternal Golden Braid*.
  The original description and the self-replication puzzle.
- Morris, H. C. (1989). Typogenetics: a logic for artificial life. In
  C. G. Langton (ed.), *Artificial Life*. Addison-Wesley. (The book's
  bibliography gives pages 341–368, Snare's 369–395.)
- Morris, H. C. *Typogenetics: A Logic of Artificial Propagating Entities*.
  PhD thesis, University of British Columbia (1989 in the book, 1988 in
  Snare).
- Varetto, L. (1993). Typogenetics: an artificial genetic system. *Journal of
  Theoretical Biology* 160(2), 185–205.
- Varetto, L. (1998). Studying artificial life with a molecular automaton.
  *Journal of Theoretical Biology* 193(2), 257–285.
- Kvasnička, V., Pospíchal, J. & Kaláb, T. (2001). A study of replicators and
  hypercycles by typogenetics. In *Advances in Artificial Life: ECAL 2001*,
  37–54. Springer.
- Typogenetics in Python, listed in the book's resources:
  <https://www.bamsoftware.com/hacks/geb/index.html>
