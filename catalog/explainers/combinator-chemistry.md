## Introduction

The combinator chemistry is Pietro Speroni di Fenizio's rebuilding of
[AlChemy](alchemy.md), Walter Fontana's chemistry of colliding programs, into
something closer to real chemistry. It is a simulation model. Its molecules
are small programs written as *combinators*: strings of single-letter
operators such as `S`, `K` and `W`, with parentheses for grouping. When two
molecules collide, one is applied to the other as a function to its input, and
the result is simplified until nothing more can be done to it. What comes out
is one or more new molecules. Nobody writes down a list of reactions; they
follow from what the molecules are.

Speroni di Fenizio started from AlChemy at the University of Sussex and
published the first version in 2000 ("A less abstract artificial chemistry",
Artificial Life VII). He then developed it with Wolfgang Banzhaf (ECAL 2001)
and made it the worked example of his 2007 PhD thesis in Jena, *Chemical
Organization Theory* (chapters 6 and 7). The thesis names two aims, both taken
over from AlChemy: to build a simple chemistry whose molecules can explore
qualitatively different states, and to explore whether a living system could be
built from the bottom up. The 2000 paper gets there by four changes to AlChemy,
applied one after the other: use combinators instead of λ-calculus terms (the
language AlChemy uses), make every molecule out of countable *atoms* that are
conserved, let a reaction produce several molecules instead of one, and keep
the reactor open with a small inflow and outflow.

Here is a whole reaction. The molecule `WR` meets the molecule `SKI`. Applying
the first to the second gives `WR(SKI)`. The letter `W` duplicates what follows
it, giving `R(SKI)(SKI)`, and the letter `R` then lets its second item go as a
separate molecule. The result is two copies of `SKI`: `WR` has turned one
molecule into two, using free atoms from the surroundings to build the second
copy and giving back its own `W` and `R`. This is equation 2 of the ECAL 2001
paper.

The chemistry's main use was as a test bed for *chemical organisation theory*,
the theory Speroni di Fenizio worked out with Peter Dittrich and others. An
*organisation* is a set of molecules that is **closed** (reactions among its
members produce only members) and **self-maintaining** (every member is
produced by reactions among the members). In a constructive chemistry, where
new kinds of molecules keep appearing, the organisations are the stable
states the soup can settle into, and the theory follows a run as it moves from
one organisation to another. The theory is the subject of chapters 12 and 13
of the book; chapter 12, "The Structure of Organizations", is co-authored with
Speroni di Fenizio.

In the book the chemistry itself gets one sentence, among "other rewriting
systems" at the end of §9.8. It is also cited in §7.3.1 as a variant of the
machine-and-tape reaction scheme that yields a variety of products, and in
§12.1 as an example of a chemistry whose population size varies and whose
reactions produce a multiset. Its closest neighbour in the catalog is
[AlChemy](alchemy.md): the same idea of molecules as functions, but AlChemy's
reactions keep both reactants, make exactly one product, need no raw material
and run in a pot of fixed size. The combinator chemistry adds conservation of
atoms, multi-product reactions, a variable population and a flow. The book
names two other variants of the machine-and-tape scheme next to it:
[McCaskill's polymer chemistry on tape](mccaskill-polymer-tm.md), where the
machine may modify the tape, and
[Ikegami and Hashimoto's machines and tapes](ikegami-hashimoto.md), where the
machine also changes state.

## How it works

### Combinators: operators that rearrange what follows them

A combinator is a string of *atoms* (the letters `B`, `C`, `I`, `K`, `R`, `S`,
`W`) and balanced parentheses. A parenthesised part is a *sub-combinator*,
treated as one item. The first atom of a string, or of any sub-combinator, is
an operator on the items that follow it. Each atom needs a certain number of
items, its *arity*, and then rewrites them:

```
B x y z -> x(yz)      regroups
C x y z -> x z y      swaps
I x     -> x          identity
K x y   -> x          destroys y
S x y z -> x z(yz)    copies z and regroups
W x y   -> x y y      copies y
R x y   -> x, y       releases y as a separate molecule
```

Here `x`, `y`, `z` stand for any item, an atom or a sub-combinator. Whatever
follows the last item needed stays where it is. One such rewriting is a
*reduction*. Six of the atoms are the classical combinators of mathematical
logic; `R` was designed for this chemistry, so that a molecule can release a
piece of itself. Parentheses around a single item, or at the start of a string,
are dropped, so `(S(I)((K)K))` is written `SI(KK)`.

A worked example from the thesis (§6.2.1): `BBBBBB` starts with `B` followed by
five items, so the first `B` acts on the next three, giving `B(BB)BB`. The
leading `B` again has three items after it, giving `(BB)(BB)`, which is written
`BB(BB)`. Now the first `B` has only one item after it and nothing more can be
done. A combinator that cannot be reduced is in *normal form*, and molecules
are always stored in normal form. Some combinators never get there: `WWW`
reduces to itself forever. Combinators also compute. `SKI` applied to anything
`x` gives `K x (I x)` and then `x`, so `SKI` acts as the identity; the thesis
leaves the reader to discover that `CKB` does the same.

### A collision is an application

When molecule `a` meets molecule `b`, the chemistry forms `a(b)`, `a` applied
to `b`, and reduces it to normal form. The order matters: `a` acts on `b`, not
the other way round. The normal form is the main product; anything `R` let go
along the way is a further product. So a reaction has the form

```
a + b -> c1 + c2 + ... + cn        (n ≥ 1)
```

In pure combinatory logic the order in which you reduce does not change the
normal form. Here it does, for three reasons the thesis gives: atoms are
limited, sizes are limited, and `R` acts on the soup outside the molecule. For
example `Kc(Rab)` becomes `c` either way, but `b` is released into the soup
only if `R` acts before `K` discards `(Rab)`. So the thesis fixes one order (its table 6.3):
first every `K` that can act, then every `B`, `C`, `R` and `I`, and only then
the outermost `S` or `W`. Destroyers go first and copiers last, which the
thesis chose to make reactions as likely as possible to succeed when atoms are
scarce.

### Atoms are conserved: the pool

Each atom type exists in a fixed number of copies. Atoms not bound in a
molecule sit in a *pool* of free atoms. When a reduction copies something
(`S`, `W`) the copy's atoms are taken from the pool; when an atom acts it
leaves the molecule and returns to the pool, and so do the atoms of whatever
`K` destroys. If the pool runs out at any point, the reaction does not happen.
The 2000 paper traces a collision with 20 atoms of each type in the pool:
`SKI` meets `KW` and forms `SKI(KW)`, which reduces to `K(KW)(IKW)` and then to
`KW`. `SKI` has acted as the identity, and its atoms, plus those of the copied
`KW`, went back to the pool.

In Chemart the pool is visible in the network: each atom type `X` has a
species `free:X`, and every reaction lists the free atoms it takes and returns.
The `WR` + `SKI` reaction of the introduction reads

```
WR + SKI + free:I + free:K + free:S -> 2 SKI + free:R + free:W
```

The second `SKI` is built from one free `S`, `K` and `I`, while the `W` and
`R` that acted return to the pool. Every reaction therefore conserves each
atom type, and the network carries one conservation law per atom.

### When a collision does nothing

A collision is *elastic*, meaning nothing happens, when the reduction is cut
off by one of four limits: too many reductions (`max_reductions`), a
combinator too large during the reduction (`max_react_size` atoms), a product
too large (`max_size` atoms) or too deeply nested (`max_depth` levels of
parentheses), or when the pool lacks the atoms needed. These four limits
decide which molecules can exist at all; the thesis calls them the size of
the space of combinators.

### Reactive and catalytic reactions

The thesis ran the chemistry at several *levels*, each a step away from AlChemy
(its table 6.2). L0 and L1 are Fontana and Buss's names for the two kinds of
organisation AlChemy produced, described under *Results*.

| level | reaction | atoms | `R` | flow | outcome in the thesis |
|---|---|---|---|---|---|
| 1 | catalytic | unlimited | no | one molecule removed per product | a single self-copier (L0) |
| 1 with filter | catalytic, copying banned | unlimited | no | same | infinite "ladder" organisations (L1) |
| 2 | catalytic | 2,000 per type | no | slow | the soup freezes |
| 3 | reactive | 600 per type | yes | inflow and outflow | moves between organisations |

In a **catalytic** reaction, as in AlChemy, both reactants survive,
`a + b -> a + b + c1 + ... + cn`, and a random molecule is removed for each one
added, so the population stays constant. In a **reactive** reaction the
reactants are used up, `a + b -> c1 + ... + cn`, so the number of molecules
changes with every reaction. *Fontana's filter* makes elastic any reaction
that would recreate one of its reactants; Fontana and Buss used it to stop
AlChemy from collapsing onto a single self-copying molecule.

### The reactor

The reactor is a well-stirred pot: at each step two distinct molecules are
drawn at random and the first is applied to the second. A *generation* is as
many collisions as there are molecules. A reactive soup would run down, since
most reactions turn two molecules into one, so the thesis (§6.2.3) adds a
flow. Each generation, every molecule decays back into free atoms with a small
probability (`prob_destroy`), and a new random molecule is assembled from the
pool with a probability that is 1 while the population is at most
`min_molecules` and halves for every further `half_add_prob` molecules. Random
molecules are built as in the ECAL 2001 paper: atoms, `(` and `)` drawn with
equal probability until an unmatched `)` ends the string.

## Using it

The default call above computes a *closure*: starting from the seven single
atoms, it applies every molecule to every other, adds the products, and
repeats. Most early reactions just join two atoms, such as `B + C -> BC`, and
the set of molecules grows quickly, so the closure stops at `max_species` (100)
and reports `status=truncated`. The default network is therefore a sample of
the first molecules the atoms can build, not an organisation. The `free:X`
species are the pool, and `net.extras["conservation"]` holds the seven
conservation laws, one per atom type. Species names are combinators in normal
form; each species' `structure` gives the same combinator fully parenthesised.

`method="closure"` is how you look at an organisation: give its generating
molecules as `molecules` and the closure finds everything they make.
`method="soup"` runs the reactor and returns the reactions that actually fired,
each with its count. Molecules you pass must be in normal form.

**A single reaction.** The reduction machinery can be called directly:

```python
from chemart.chemistries.combinator_chemistry import Chemistry, normal_form
Chemistry().react("WR", "SKI")   # (('SKI', 'SKI'), Counter({'S': 1, 'K': 1, 'I': 1}))
normal_form("BBBBBB")            # ['BB(BB)']
normal_form("WWW")               # None: no normal form within the limits
```

The second element of `react` is the number of free atoms of each type the
reaction needs from the pool at its peak.

**The closure of two molecules.**

```python
net = chemart.generate_network("combinator-chemistry", molecules=["WR", "SKI"])
# 9 species, 3 reactions, status=complete
# WR + SKI + free:I + free:K + free:S -> 2 SKI + free:R + free:W
# SKI + WR -> WR + free:I + free:K + free:S
# 2 SKI -> SKI + free:I + free:K + free:S
```

`SKI` applied to anything returns it unchanged, so it simply falls apart into
free atoms; `WR` applied to `SKI` copies it.

**The level-1 ladder organisations (thesis figure 6.5).** Catalytic reactions,
the six atoms without `R`, Fontana's filter, and the limits of the thesis
figure 6.4:

```python
LEVEL1 = dict(reaction="catalytic", atoms="BCIKSW", filter_reproduction=True,
              max_size=15, max_react_size=15, max_depth=7, max_reductions=100)
net = chemart.generate_network("combinator-chemistry",
                               molecules=["BKK", "BK(BKK)"], **LEVEL1)
# 15 species, 123 reactions, status=complete:
# BKK, K(BKK), ..., K(K(K(K(K(K(K(BKK)))))))  and  BK(BKK), K(BK(BKK)), ...
# BKK + BKK -> BKK + BKK + K(K(BKK))
```

The two ladders are infinite in the thesis; here they stop at `max_depth` 7.

**The type A organisation of the 2000 paper.** `K` releases instead of
destroying (`k_action="release"`), no `R`:

```python
ALPHA = "C(C(K(CKK))(WC))(C(K(CKK))(WC))"
net = chemart.generate_network("combinator-chemistry", molecules=[ALPHA],
        k_action="release", atoms="BCIKSW", max_size=100, max_react_size=100, max_depth=6)
# 18 species (12 molecules plus 6 free:X), status=complete:
# ALPHA, K, KK, K(ALPHA) ... K(K(K(ALPHA))), K(KK) ... K(K(K(K(K(KK)))))
```

**A reactive soup.** `method="soup"` with the defaults runs 20 generations of a
reactive soup of about 100 random molecules, with 200 atoms of each type. It
takes under two seconds, and with seed 1 the population falls from 103 to 5
molecules: most reactions turn two molecules into one, and the inflow adds only
one molecule per generation. The history is in `net.extras["analysis"]`
(`population`, `diversity` and `free_atoms` per generation), and the end state
in `net.extras["final_state"]` and `net.extras["final_free_atoms"]`.

Organisations take thousands of generations to appear. With 600 atoms per type
(the thesis figure 6.7) and 3,000 generations:

```python
net = chemart.generate_network("combinator-chemistry", seed=2, method="soup",
                               generations=3000, atoms_per_type=600)
a = net.extras["analysis"]
a["population"][::500]   # [100, 13, 196, 218, 226, 222, 230]
list(net.extras["final_state"].items())[:2]
# [('B(SIR)(R(B(SIR)))', 176), ('C(B(SWW)(S(SWW)))(WW)', 6)]
net.extras["final_free_atoms"]
# {'B': 108, 'C': 547, 'I': 196, 'K': 598, 'R': 0, 'S': 4, 'W': 333}
```

The population crashes, recovers once a successful molecule appears, and
settles around 220, dominated by one molecule, with the pool out of `R`. This
run takes about 35 seconds; with seed 3 the soup settles on
`C(R(WC))(R(WC))` (166 copies, the pool out of `C` and `W`), and with seed 1
the population climbs to about 270, collapses to 2 near generation 2,000 and
climbs again, in about 2.5 minutes. The thesis's own runs lasted 10,000
generations; at that scale pure Python is slow.

**A level-1 catalytic soup.**

```python
L1 = dict(method="soup", reaction="catalytic", atoms="BCIKSW",
          max_size=15, max_react_size=15, max_depth=7)
net = chemart.generate_network("combinator-chemistry", seed=4, M=100, generations=5000, **L1)
net.extras["final_state"]   # {'WK': 100}
```

Seed 4 ends with every molecule `WK`, and seed 5 with every molecule `I`: each
applied to itself gives itself, the L0 organisation. Seeds 1, 2 and 3 instead
end with 7, 42 and 2 kinds of molecule between which every collision is
elastic, because the products would exceed the size limits. Each run takes
7 to 14 seconds.

**A universal copier (thesis chapter 7).** With the limits of chapter 7,
`B(WR)(WR)` turns any molecule into four copies of it:

```python
CH7 = Chemistry(max_reductions=10_000, max_react_size=100, max_size=100, max_depth=20)
CH7.react("B(WR)(WR)", "SII")[0]   # ('SII', 'SII', 'SII', 'SII')
```

`B(WR)(WR)y` becomes `WR(WRy)`, then `R(WRy)(WRy)`, which releases one `WRy`;
each `WRy` becomes `Ryy`, which gives `y` and releases another `y`.

## Results

### AlChemy reproduced with combinators

The first level of the thesis (§6.3) repeats Fontana's experiments with
combinators: catalytic reactions, six atoms, no shortage of atoms, a pot of
fixed size. Like AlChemy, the soup collapses. In the run of figure 6.3, diversity
drops within about 50 biological generations (a clock that counts only
collisions that react) to a single molecule `A` with `A + A -> A`, the organisation Fontana called level 0 (L0);
which molecule wins changes from run to run. With Fontana's filter the soup
takes much longer to converge and ends in an infinite organisation made of
"ladders" `BKK, K(BKK), K(K(BKK)), ...`, in one run (figure 6.5) first with
three ladders at generation 137, then two at generation 405, then one at 605.
This mirrors the simplest level-1 (L1) organisation of Fontana and Buss
(1994). Section 6.3.3 analyses that organisation with organisation theory:
every ladder is an organisation by itself, so n ladders give 2ⁿ organisations
arranged as an n-dimensional cube, and since a lost ladder cannot come back the
soup slides down to a single ladder.

Chemart reproduces the organisations: its tests check that the closures of the
three seed sets of figure 6.5 are exactly the ladders, cut at depth 7, and that
every reaction in them keeps its reactants. The collapse to L0 is not checked
by the tests. In the runs shown under *Using it*, two of five seeds collapse to
a single self-copier within 5,000 generations and three freeze first.

### The 2000 paper: organisations of two kinds

The 2000 paper uses reactive reactions with `K` releasing instead of
destroying (there is no `R` yet), 2,000 atoms of each type, 300 starting
molecules and limits of 100 atoms and 20 levels of parentheses, over 180 runs.
Its findings: the number of molecules, although free to vary, settles around
an average after a transient; every run reaches a different organisation, so
the space of organisations is rich; some organisations are finite and some
infinite. The organisations fall into two classes, "nearly equally divided".
Class A organisations shrink until new random molecules arrive, then use them
to grow again, and stay near 300 molecules; class B organisations grow until
they have used up the atoms, reaching 350 to 1,000 molecules.

The class A example is generated by one molecule,
`α = C(C(K(CKK))(WC))(C(K(CKK))(WC))`, which applied to any molecule `a`
returns `α`, `K` and `a`. The released `K` then builds `K(a)`, and molecules
`Kⁿ(b)` peel off one `K` per reaction, so the organisation is
`{α, K, KK, Kⁿ(α), Kⁿ(KK)}`. In the run shown it emerged after about 210
generations and then held between 296 and 330 molecules. The class B example
is generated by four molecules, among them `γ = S(K(SSK))(K(K(SSK)))`, which
applied to `a` gives `γ, a, a`, and `δ = B(WW)(W(B(WK)W))`, which gives
`aaa(aaa)(aaa), Wa, Wa, Wa`; it ended with about 900 molecules of 4 to 6
kinds. The paper reads class A as a *metabolism* in the sense of Bagley and
Farmer: it takes in material from its environment and reassembles it into more
of itself.

Chemart's tests check the reactions of `α`, `Kⁿ(b)`, `γ` and `δ` with an
arbitrary molecule, and that the closure of `α` is exactly the class A
organisation cut at depth 6. Two published counts are not matched: the paper
gives `δ * α -> 9 α, 3 W(α)` and `δ * γ -> 17 γ, 3 W(γ)`, while Chemart gets
different numbers (the implementation decisions list them), because the counts
depend on the reduction path and the original code is not available. The
statistics of the 180 runs are not reproduced.

### ECAL 2001: metabolic and balanced organisations

Speroni di Fenizio and Banzhaf (2001) add `R` and make `K` destroy again, so
molecules can both release and destroy pieces. They ran 150 runs of 10,000
generations and 26 of 30,000, with 2,000 atoms per type and a random inflow
that never stops. Their first finding is that no organisation is totally
stable: some runs held one for many thousands of generations and then switched
to another. They name two types.

A **metabolic organisation** has few kinds of molecule, often one or two, and
cannot grow by itself, but grows on the inflow. Their example is a molecule `a`
with `a * a -> a, a`, invaded by a random `S`: `S` meets `a` three times,
building `Sa`, `Saa` and then `aa(aa)`, which falls apart into four `a`. Three
`a` went in and four came out. In the active phase such an organisation digests
nearly every molecule thrown in, the number of molecules rises, and the free
atoms fall until one atom type runs out; then a resting phase begins in which
digestion stops halfway and diversity rises. Metabolisms were very rare in an
experiment without `W`. A **balanced organisation** holds molecules that build
and molecules that destroy, splits incoming molecules into atoms and rebuilds
itself from the pool, and fluctuates around an average size. In both, the
exhaustion of one atom type often pushed the soup into a different
organisation.

Chemart's tests check the paper's worked reaction `WR * SKI -> SKI, SKI`, and
the 2000 paper's pool trace of `SKI` on `KW`. The run statistics and phases are
not reproduced or tested. The molecule that took over the seed-2 soup under
*Using it*, `B(SIR)(R(B(SIR)))`, behaves like the paper's `a`: applied to any
molecule it returns itself and that molecule, and three collisions with a
free `S` turn three copies into four, exactly the paper's pathway.

### Level 2 and level 3 of the thesis

Level 2 limited the atoms while keeping catalytic reactions, hoping that
scarcity would stop one molecule from dominating. Most runs froze instead,
with too many atoms locked up and no reaction possible; with plenty of atoms
the soup ignored the limit and reached L0. Chemart does not offer this level.

Level 3 is the reactive system with `R` and a flow. Here the soup "does not
freeze anymore": it moves from organisation to organisation. The run shown in
the thesis (figures 6.7 and 6.8, 600 atoms per type, 10,000 biological
generations) was picked from 30 for its clear transitions; population size,
diversity and free atoms change together whenever the soup changes
organisation, and the lattice of organisations shows the composition changing
between generations 5,000, 7,000 and 9,000. Another run (figure 6.2) found its
first organisation at physical generation 2,357 and moved to a second at
10,212, which uses a different set of atoms. Not every run moves: some stay in
one organisation for the whole experiment. Chemart implements this level (the
default soup) but its tests only check that atoms are conserved and the
population varies. The runs under *Using it* show the qualitative picture of a
settled population, dominated by few molecules and limited by an exhausted
atom type, and in one case a collapse and regrowth, but no systematic
comparison with the thesis has been made.

### Combinators in space

Chapter 7 of the thesis (published at ECAL 2001 with Dittrich and Banzhaf)
places the molecules on the nodes of a planar triangular graph, so that only
neighbours react. Molecules that cannot react with each other cluster into
membranes that cut the graph into separate regions, which the authors
interpret as proto-cells, some of them autopoietic (self-producing their
membrane). The chapter lists the reactions of the main molecules of one run,
for example `B(WR)(WR) + y -> 4y` and `R(x) + y -> x, y`. Chemart has no
spatial version. Its tests check the two general rules and 19 of the tabulated
reactions; the implementation decisions record that the fixed reduction order
reproduces 20 of the 25 entries of the tables, and list the five it does not.

### Limits the author reported

In §6.6 the thesis says the chemistry fell short of its aim, a self-sustaining,
progressively complexifying process, for three reasons. It is brittle: a single
change to a molecule gives a molecule with entirely different behaviour, so
there are no small steps to explore. It makes universal copiers too easily, a
molecule `a` with `a + b -> a, b, b` for any `b`, which spreads through the
soup and wipes out diversity. And the seven atoms contain two complete bases
of combinatory logic (`B, C, K, W` and `S, K, I`), so a random set of molecules
often generates every possible molecule, and the lattice of organisations then
says little. Some molecules also react with nothing and lock up atoms, which is
why the outflow is needed. The thesis concludes that a less brittle chemistry,
whose smallest self-maintaining sets need many molecules, would produce more
innovation. Chemart shows the copier (`B(WR)(WR)` above) but does not test
the loss of diversity.

### Standing

The thesis describes the λ-calculus and combinator chemistries as the first
applications of chemical organisation theory, followed later by applications
to real systems by Matsumaru, Centler, Dittrich and co-authors, and says both
AlChemy and the combinator chemistries "gave mixed results" towards their two
aims. The theory the chemistry served as a test bed for is presented in
chapters 12 and 13 of the book, and is also the basis of
[organisation-oriented computing](organization-computing.md), of which
Speroni di Fenizio is a co-author. The book says nothing further about later
use of the combinator chemistry itself.

## Further reading

- Speroni di Fenizio, P., Dittrich, P. & Banzhaf, W. (2001). Spontaneous
  formation of proto-cells in an universal artificial chemistry on a planar
  graph. In J. Kelemen & P. Sosík (Eds.), *Advances in Artificial Life*
  (ECAL 2001), LNCS 2159, 206–215. Springer. The spatial version, chapter 7 of
  the thesis.
- Fontana, W. & Buss, L. W. (1994). "The arrival of the fittest": toward a
  theory of biological organization. *Bulletin of Mathematical Biology* 56,
  1–64. The L0 and L1 organisations and the filter.
- Hindley, J. R. & Seldin, J. P. (1986). *Introduction to Combinators and
  λ-Calculus*. Cambridge University Press. The combinatory logic behind the
  atoms (book reference [383]).
