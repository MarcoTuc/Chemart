## Introduction

AlChemy is Walter Fontana's "algorithmic chemistry", introduced in 1991 and
developed with the Yale biologist Leo Buss. It is a simulation model in which the molecules are small computer programs. When two
of them collide, one is run with the other as its input, and whatever the
computation returns is a new molecule. Nobody writes a list of reactions: the
reactions follow from what the molecules are, just as in real chemistry the
products follow from the structure of the reactants.

The programs are written in the *λ-calculus* (lambda calculus), the minimal
language of functions that underlies functional programming languages such as
Lisp. Everything in it is a function, and a function can take another function
as its input. That is what Fontana and Buss wanted. They argued (Fontana & Buss
1994, "The arrival of the fittest") that two features of chemistry matter for
the origin of biological organisation. The first is *construction*: molecules
combine to make new molecules. The second is *equivalence*: many different
combinations lead to the same stable product, which is what lets reactions
close up into networks instead of branching out for ever. In the λ-calculus,
applying one expression to another builds a new expression, and simplifying it
to its *normal form* (the version that cannot be simplified any further) maps
many different combinations onto one product. Molecules are those normal
forms.

The question behind the model comes from evolutionary theory. Fontana and Buss
quote DeVries (1904): "Natural selection may explain the survival of the
fittest, but it cannot explain the arrival of the fittest." Population
genetics assumes that genes, individuals and populations already exist; it
does not say where organised entities come from. AlChemy asks what kinds of
organisation appear, without any selection, in a pot of random interacting
functions. The answer came in three levels. Left alone, the pot usually
collapses to one molecule that copies itself (Level 0). If copying is
forbidden, it settles into a set of molecules that keep producing one another,
an *organisation* (Level 1). Two such organisations can combine into a larger
one (Level 2). The vocabulary of this work is *closure* (the set produces
nothing outside itself), *self-maintenance* (every member is produced from
within the set) and *organisation* (both at once); the book's formal
definitions of these terms (§12.1) follow Fontana's papers.

Banzhaf and Yamamoto present AlChemy as the first example of a rewriting
system (book §9.1). Its closest neighbour in the catalog is
[combinator chemistry](combinator-chemistry.md), which keeps the idea but
replaces the λ-calculus by combinators, expressions without variables that are
easier to reduce. [Matrix chemistry](matrix-chemistry.md) and the
[automata reaction](automata-reaction.md) chemistry share the reaction scheme
in which both reactants survive and one acts on the other; the book notes that
experiments with the latter support AlChemy's findings. What sets AlChemy
apart is that its molecules form a universal programming language, so the
space of possible molecules and actions is unbounded.

## How it works

### Molecules are functions

A λ-expression is built from three things: a variable such as `x`; an
*abstraction* `λx.M`, the function that takes an input `x` and returns `M`;
and an *application* `(M)N`, the function `M` applied to the input `N`
(Fontana and Buss always write the operator in parentheses). Running a
function on an input is called *β-reduction*: `(λx.M)N` becomes `M` with every
`x` replaced by `N`. Repeating this until no step applies gives the normal
form.

Two simple examples. The identity `I = λx.x` returns its input unchanged. The
function `K = λx.λy.x` takes an input and returns a function that ignores its
own input and returns the first one.

The names of the variables do not matter (`λx.x` and `λy.y` are the same
function), so Chemart names molecules by position instead. In a species id,
`^` stands for a `λ`, `(M)N` is an application, and a number `k` is the
variable bound by the `k`-th enclosing `λ`, counting outwards. So `I` is `^1`
and `K` is `^^2`. This is called de Bruijn notation. Each species also carries
the paper's readable form as its `structure`, with the variables renamed
`x1, x2, ...` in the order their `λ`s appear: `^^2` is `λx1.λx2.x1`.

### A collision

When molecules `s1` and `s2` collide, `s1` is applied to `s2`, the result is
reduced to normal form, and the product is added to the pot. Neither reactant
is used up:

```
s1 + s2 -> s1 + s2 + NF((s1)s2)
```

Here is a closure Chemart computed from the two molecules `I` and `K`
(the recipe is under *Using it*):

```
2 ^1 -> 3 ^1  [mass-action k=1.0]
^1 + ^^2 -> ^1 + 2 ^^2  [mass-action k=1.0]
^^2 + ^1 -> ^^2 + ^1 + ^^1  [mass-action k=1.0]
2 ^^2 -> 2 ^^2 + ^^^2  [mass-action k=1.0]
```

Read the first reactant as the operator (`2 ^1` is `I` meeting `I`). `I`
applied to anything returns it, so `I` acting on `I` makes another `I`, and `I` acting on `K` makes another
`K`. These are *copy actions*: the product is identical to one of the
reactants. `K` acting on `I` is `(λx.λy.x)I`, which reduces to `λy.I`, that is
`λx1.λx2.x2`, a new molecule `^^1`. `K` acting on itself gives `λy.K`, the new
molecule `^^^2`. A reaction thus depends on the order of the two molecules,
and because both orders are equally likely when two molecules meet, each
ordered collision counts once. When both orders give the same product the
reaction's rate constant `k` is 2.

### When a collision does nothing

Not every λ-expression has a normal form: `(λx.(x)x)λx.(x)x` reduces to itself
for ever. Fontana and Buss therefore used *pragmatic reduction*, with a limit
on the number of reduction steps and on the size of the expression during
reduction. If a limit is exceeded, or the normal form is too large, the
collision is *elastic*: the two molecules bounce off and nothing is made.

The experimenter can also declare collisions elastic on purpose. These
*boundary conditions* are how the three levels were produced:

- the **no-copy** condition bans every collision whose product is identical to
  one of the two reactants;
- **syntactic filters** ban products that contain a given pattern; Fontana and
  Buss wrote them as regular expressions, for example to exclude three `λ`s in
  a row.

The paper also fixes one standard condition: an operator must begin with a
`λ`, since otherwise applications pile up without ever reducing. In Chemart
molecules have no free variables (every variable is bound by some `λ`), and a
normal form without free variables always starts with a `λ`, so this condition
always holds.

### The reactor

The pot is a well-stirred *flow reactor* with a fixed number `M` of molecules.
It starts with `M` distinct random normal forms. At each step two molecules
are drawn at random and collide. If the collision is reactive, the product is
added and a randomly chosen molecule is removed, which keeps the size constant
and slowly washes out anything that is not being remade. Common molecules meet
more often and are removed more often, so this is mass-action kinetics. The
book notes that Fontana and Buss ran it with `M` = 1,000 to 3,000.

The same picture has a deterministic version, the paper's equation 18, for the
relative concentrations `x_i`: `dx_i/dt = Σ_jk a^i_jk x_j x_k − x_i Φ`. Here
`a^i_jk` is 1 if `j` acting on `k` gives `i` and 0 otherwise, and `Φ` is the
total production, an outflow that keeps the concentrations summing to 1.

### Closure, self-maintenance, organisation

Instead of running the reactor one can ask what a set of molecules can make at
all. The *closure* of a set is everything reachable by repeated collisions
among its members and their products. A set is *closed* if collisions among its
members produce only members, and *self-maintaining* if every member is
produced by some collision within the set. The *center* of an organisation is
its smallest self-maintaining subset. Most organisations are infinite, so a
finite reactor holds only part of one at a time.

## Using it

The default run is a small version of the basic experiment: `M = 100` random
molecules, no filter, 2,000 collisions. It does not reproduce a published run;
the published reactors were ten times larger and ran for hundreds of thousands
of collisions. The network lists every distinct reaction that fired, with
`count` for how often. The run's history is in `net.extras`:

```python
a = net.extras["analysis"]
a["distinct_species"][:4], a["distinct_species"][-1]   # ([100, 87, 80, 63], 32)
a["elastic_pairs"]                                     # {'no_normal_form': 22}
list(net.extras["final_state"].items())[:3]
# [('^^(2)1', 13), ('^^^^^1', 7), ('^^^3', 7)]
```

`distinct_species` is the number of different molecules in the pot, sampled
every `collisions_per_sample` collisions (here every 100). Diversity has
fallen from 100 to 32. `elastic_pairs` counts the ordered pairs whose collision
was elastic, by reason: `no_normal_form` (a reduction limit was hit), `copy`
and `forbidden` (the two filters). `final_state` is the pot at the end. Its most
common molecule, `^^(2)1` or `λx1.λx2.(x1)x2`, copies any molecule that starts
with a `λ`, which here means every molecule. The third reaction of the default
network shows it at work:

```
^^(2)1 + ^^(2)^1 -> ^^(2)1 + 2 ^^(2)^1  [mass-action k=1.0]  (x1)
```

Applied to `λx1.λx2.(x1)λx3.x3`, it returns the same molecule. Given enough
time, a molecule like this takes over the pot.

**Level 0: collapse to a copier.** Run longer:

```python
net = chemart.generate_network("alchemy", seed=0, collisions=20000)
net.extras["final_state"]          # {'^1': 100}
```

All 100 molecules are the identity. Not every seed ends this way (see
*Results*).

**Level 1: ban copying.** Add `filter="no-copy"`:

```python
net = chemart.generate_network("alchemy", seed=1, collisions=20000, filter="no-copy")
names = {s.id: s.structure for s in net.species}
for sid, n in net.extras["final_state"].items():
    print(n, sid, names[sid])
```

```
27 ^^^2 λx1.λx2.λx3.x2
22 ^^^^2 λx1.λx2.λx3.λx4.x3
21 ^^2 λx1.λx2.x1
20 ^^^^^2 λx1.λx2.λx3.λx4.λx5.x4
9 ^^^^^^2 λx1.λx2.λx3.λx4.λx5.λx6.x5
1 ^^^^^^^2 λx1.λx2.λx3.λx4.λx5.λx6.λx7.x6
```

Every survivor has the same shape: a chain of `λ`s ending in the variable of
the second-to-last one. This is one family of the projector organisation
described under *Results*. Meanwhile 115 ordered pairs were refused as copy
actions (`elastic_pairs`).

**Closure instead of a reactor.** `method="closure"` returns every reaction
reachable from a seed set, given as λ-terms in `terms` (write `λ` or `\`,
applications as `(M)N`, or give de Bruijn ids). Organisations are usually
infinite, so the closure stops at `max_species` and reports `status="truncated"`.
The `I`, `K` example above is

```python
net = chemart.generate_network("alchemy", method="closure",
                               terms=["λx.x", "λx.λy.x"], max_species=6)
net.summary().splitlines()[0]      # 'alchemy: 6 species, 34 reactions, status=truncated'
```

In closure mode `net.extras["analysis"]` has `elastic` (the same counts) and
`self_maintaining`, which says whether every returned species is produced by a
collision among the returned species. Seeding with the center of an
organisation and setting `max_species` to its size tests it directly:

```python
center = ["λx1.λx2.λx3.x1", "λx1.λx2.λx3.λx4.x2", "λx1.λx2.λx3.λx4.λx5.x3"]
net = chemart.generate_network("alchemy", method="closure", terms=center,
                               filter="no-copy", max_species=3)
for r in net.reactions:
    print(r.to_text())
net.extras["analysis"]             # {'elastic': {'copy': 2}, 'self_maintaining': True}
```

```
2 ^^^3 -> 2 ^^^3 + ^^^^^3  [mass-action k=1.0]
2 ^^^^3 -> 2 ^^^^3 + ^^^3  [mass-action k=1.0]
^^^^3 + ^^^^^3 -> ^^^^3 + ^^^^^3 + ^^^3  [mass-action k=1.0]
^^^^^3 + ^^^3 -> ^^^^^3 + ^^^3 + ^^^^3  [mass-action k=1.0]
2 ^^^^^3 -> 2 ^^^^^3 + ^^^^3  [mass-action k=1.0]
```

Each of the three is made by the others. Drop the last one
(`terms=center[:2], max_species=2`) and `self_maintaining` is `False`.

**Other boundary conditions.** `forbidden_patterns` takes regular expressions
matched against the product's species id: `[r"\^\^\^"]` bans three `λ`s in a
row, the filter of the paper's second Level 1 example. `mediator` sets a
different collision rule: instead of `(s1)s2`, the product is the normal form
of `((Φ)s1)s2` for a λ-term `Φ` of your choice. The reduction limits are
`max_steps`, `max_size` and `max_nf_size`, and the random molecules are shaped
by `p_variable`, `p_abstraction`, `max_depth`, `p_bound` and `n_free` (see the
parameter table and the implementation decisions: the generator's defaults are
Chemart's choices, not published values).

**Rates.** Every reaction carries a mass-action rate constant (1 or 2, see
above) and the network has a constant-total outflow, so a closure network is
equation 18 ready to integrate with any ODE solver.

**Speed.** The default run takes about 2 seconds. Runs with `M=100` and 20,000
collisions take 1–4 seconds. A run the size of the paper's reactor,
`M=1000, collisions=100000, filter="no-copy"`, took 13 seconds and ended with
47 distinct molecules.

## Results

The results below come from Fontana and Buss's long paper (1994, Bulletin of
Mathematical Biology; section numbers refer to the Santa Fe Institute working
paper 93-09-055), summarised in their PNAS paper of the same year, and from a
2024 re-examination by Mathis, Patel, Weimer and Forrest that ran the original
code again. The PNAS abstract states the three findings as generic, expected
to reappear "if 'the tape were run twice'": hypercycles of self-reproducing
objects arise; if self-replication is inhibited, self-maintaining
organisations arise; and these can combine into higher-order organisations.

### Level 0: copiers

With no filter, the pot quickly loses diversity. Early collisions keep making
new molecules; then innovation stops, and the survivors form a small set that
is closed under interaction and in which every molecule is copied by some
molecule: for each `f` there is a `g` with `g` acting on `f`, or `f` acting on
`g`, returning `f`. "In many instances the system reduces to just one object
species that is a self-copier" (§6.1), the simplest being the identity. In
other runs a small ecology of mutual copiers survives, often a *hypercycle*,
a ring in which each member copies the next. The paper's Figure 1 shows two,
one with two members and one with three. Such ecologies are fragile: when a
few random molecules are added they typically collapse to a single copier.
The book adds that the survivors tend to look alike syntactically: the system
"quickly exhausts its potential for innovation" (§9.1).

Chemart checks this in its tests. Both ecologies of Figure 1 are closed and
satisfy the copy condition, and the caption's example actions are reproduced.
Reactors with seeds 0 and 1 (`M=100`, 20,000 collisions) end with at most two
species forming a closed copy ecology. The test file notes that not every seed
does so. Over seeds 0–19 with those settings, my runs ended with one species
(always the identity) in 7 runs and with 2 to 41 species in the others; some of
these were a projector family like the Level 1 one below.

### Level 1: organisations

Ban copy actions and the pot never closes: new molecules keep appearing. But
they are molecules that were there before and were washed out. The pot has
settled in a region of the space of molecules. Fontana and Buss describe what
an observer who knows no λ-calculus would find there: a *grammar* that
describes every molecule, a few *algebraic laws* that predict every collision,
and *self-maintenance*. They take these three properties together as their
definition of an organisation.

**Example 1, the projectors (§6.2.2).** The simplest organisation they found
consists of the molecules `A(i,j) = λx1...λxi.xj`: `i` nested `λ`s returning
the `j`-th variable. Two laws describe every collision:

- if `j > 1`, `A(i,j)` acting on anything gives `A(i−1, j−1)`;
- `A(i,1)` acting on `A(k,l)` gives `A(k+i−1, l+i−1)`.

The first law moves molecules down a diagonal of the `(j, i)` plane (the
paper's Figure 2); the second moves them back up. Each diagonal, the molecules
with a fixed `i − j`, is a *family* that is closed on its own. The center of the family with
`i − j = d` is its `d + 1` smallest members, and a population in the reactor drifts towards the
center and stays there. Different families compete neutrally, and in the long
run one remains; the Level 1 run under *Using it* ended in the family with
`i − j = 1`. For the family whose center has ten members the paper also built the
ODE (equation 18) of the center alone: with equal rates it has one stable
fixed point; when the basic cycle `A(10,1)` → `A(19,10)` → ... → `A(11,2)` →
`A(10,1)` runs five times faster than the other reactions, the dynamics become
a limit cycle.

**Example 2, integers (§6.2.3).** Banning the projectors with a filter against
three `λ`s in a row gives a different organisation: two families of molecules,
written as strings `a...aA` and `a...aB`, that behave like integers under a
difference operation (`i` acting on `j` gives `j − i` if `i ≤ j`, otherwise
`i − j − 1`, with shifts between the families). Its center has four members.

**Robustness.** Organisations repair themselves. Injecting a few random
molecules (typically three molecules in 10 copies each) usually changes
nothing; only in rare cases does the organisation absorb them, by adding new
building blocks to its grammar. Once an organisation exists, copy actions can
be allowed again without harm. The ban matters only at the start, because
copiers grow faster than organisations form. Fontana and Buss found other
ways past this race (§6.3): accepting the product of a copy action with
probability 0.75 instead of 1 "was typically sufficient to permit
organizations to develop"; starting from more complex random molecules also
helped; and one organisation arose around a replicator kept in check by
parasites it copied.

Chemart checks the projector laws for every pair with `i, k ≤ 6`, the basic
cycle, that the centers of the families with `i − j` = 1, 2 and 3 are
self-maintaining while smaller subsets are not, and that a no-copy reactor
(seed 1) settles in a single family with no copy reaction. For example 2 it
checks the laws, the integer arithmetic, and that the four-member center is
self-maintaining while three-member subsets are not. It also checks the paper's
example 4. It does not simulate injections, reduced copy efficiency, the
parasite organisation, or the limit cycle.

### Level 2: organisations of organisations

When two Level 1 organisations share a reactor, cross-collisions can make
molecules that belong to neither, which Fontana and Buss called the *glue*.
If the glue exists and stabilises both, the result is a Level 2 organisation
(§6.4). In the first example they merged two separately grown organisations,
of 54 and 41 species and 1,000 molecules each, into a reactor of capacity
3,000, with a different collision rule for the glue; both were stably
maintained. In the second, two organisations arose spontaneously in one
reactor of capacity 1,000. The larger, B, kept 650–800 species on its own.
When the simpler organisation A expanded, diversity halved from about 620 to
just over 300, and B vanished, but A and the glue regenerated it repeatedly
until B was finally excluded at about 380,000 collisions. The paper describes
A exactly: it is built from one molecule `T` and two prefixes, and `T` acting
on itself leads round a cycle of molecules back to `T` (its Figure 6).
Spontaneous Level 2 organisations are "extremely rare" (book §9.1).

Chemart checks the cycle of Figure 6 and the laws of organisation A. It does
not reproduce the merger experiment: its `mediator` applies one collision rule
to every collision, while the paper used a separate rule only for collisions
involving glue.

### The 2024 re-examination

Mathis, Patel, Weimer and Forrest (2024) recompiled the original C code, which
the Santa Fe Institute still hosts, and ran it at statistical scale. They
reproduced the key results with some surprises:

- Unfiltered (Level 0) runs, 1,000 molecules and 6 million collisions with
  100 random molecules added every million, did not always collapse. Across
  1,000 simulations the final number of distinct molecules ranged over orders
  of magnitude, with some runs keeping tens or hundreds; one kept about 380
  and was unaffected by the additions. These organisations survived replacing
  large fractions of the pot with the identity: replacing 90% destroyed only
  some of them.
- Level 1 organisations varied widely in stability. Re-running one with a
  different random seed made it drift slowly into a different stable state.
- Merging two organisations from different runs (455 pairs) rarely gave
  coexistence; usually one dominated or both were destroyed.
- The results depend on how random molecules are made. With the original
  generator, and with free variables bound by adding `λ`s in front, rich
  organisations appear; with a generator that samples expression trees more
  uniformly, no-copy runs collapse to an inert pot of identities. Leaving
  free variables unbound let one term take over regardless of filters, from
  which they infer that the original runs bound them.
- They proved that a variant based on the *typed* λ-calculus can reproduce the
  state transitions of any chemical reaction network.

Chemart does not reproduce these statistics. Its random generator follows the
original in binding free variables but keeps the probabilities constant with
depth, and its step limit counts β-reductions rather than the original's
rewriting steps, so its run statistics are not directly comparable.

### What Chemart covers, in short

The calculus, the collision rule, pragmatic reduction, both boundary
conditions, the generalised collision rule, the flow reactor and the closure
are implemented, and the published Level 0, 1 and 2 structures are checked
exactly where the paper gives them in closed form. The tests also check that
closure networks give equation 18. What is missing is the large-scale reactor
statistics (paper or 2024), the dynamics of the Level 2 experiments, and the
variants of §6.3.

## Further reading

- The original AlChemy code, recompiled, with a Docker container and analysis
  scripts: <https://github.com/colemathis/AlChemy>
- Fontana, W. (1991). Algorithmic chemistry. In C. G. Langton, C. Taylor,
  J. D. Farmer & S. Rasmussen (eds.), *Artificial Life II*, 159–210.
  Westview Press. Book reference [281], which the book cites for Fontana's
  original λ-chemistry.
- Fontana, W. & Buss, L. W. (1996). The barrier of objects: from dynamical
  systems to bounded organizations. In J. Casti & A. Karlqvist (eds.),
  *Boundaries and Barriers*, 56–116. Addison-Wesley. Book reference [284],
  which the book cites for the link to proof theory (§16.5, see
  [proof search as a chemistry](proof-ac.md)).
