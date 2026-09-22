## Introduction

Combinatory Chemistry is an artificial chemistry built on combinatory logic.
Germán Kruszewski and Tomas Mikolov introduced it at the ALIFE conference in
2020, looking for a model in which units that can evolve emerge by themselves.
In their words, such units must "(1) preserve themselves in time (2)
self-reproduce and (3) tolerate a certain amount of variation when
reproducing".

Combinatory logic is a way of writing programs with no variables at all. Every
program is built by applying three basic functions, the combinators `S`, `K`
and `I`, to each other. Together they can compute anything a computer can.
Kruszewski and Mikolov turn this into a chemistry. A molecule is a program,
and running one step of a program is a reaction. The twist is that the
chemistry conserves its atoms. Plain combinatory logic throws parts of an
expression away and duplicates others. Here nothing may appear or vanish:
whatever is thrown away is released into the soup, and whatever is copied
has to be taken from it.

The simulation starts from free atoms only, thousands of loose `S`, `K` and
`I`. Random joining and splitting make the first expressions, and some of those
can compute. After a while the soup holds structures that keep themselves
going by eating one particular kind of molecule from their surroundings. The
authors compare these cycles to metabolisms. They found three kinds: simple
*autopoietic* structures that rebuild themselves, *recursive* ones that keep
growing, and *self-reproducing* ones that double.

The model is a close relative of [AlChemy](alchemy.md), where molecules are
λ-calculus terms, and of Speroni di Fenizio's
[combinator chemistry](combinator-chemistry.md), which also uses combinators
and conserves atoms. It differs from both in three ways. It reduces one step
at a time instead of to a final *normal form*, so programs that never finish
can take part. It needs no outside rule to keep the population bounded,
because conservation does that. And it starts from single atoms rather than
from random programs. The paper appeared in 2020, after Banzhaf and Yamamoto's
book, so the book does not cover it.

## How it works

### Expressions

An expression is a tree of applications. `SII` means "apply `S` to `I`, then
apply the result to `I`"; application groups to the left, so `SII` is
`((SI)I)`. Parentheses are only written where they change the grouping:
`SII(SII)` applies `SII` to a second copy of `SII`. The *size* of an
expression is its number of atoms, so `SII(SII)` has size 6.

### What the combinators do

In combinatory logic each combinator rewrites the arguments that follow it.
Such a spot is a *redex* (a reducible expression):

- `I f` becomes `f`: the identity.
- `K f g` becomes `f`: `g` is thrown away.
- `S f g x` becomes `f x (g x)`: `x` is used twice.

An expression can contain several redexes, one inside another. When it has
none, combinatory logic calls it a *normal form*.

### Reactions that conserve atoms

Applied as they stand, these rules would destroy the combinator that fires,
lose `g`, and create a copy of `x` from nothing. Combinatory Chemistry repairs
each rule so the number of each atom stays fixed (paper eqs. 1-3; α and β
stand for the rest of the expression around the redex):

```
α(I f)β          ->  α f β          + I
α(K f g)β        ->  α f β          + g + K
α(S f g x)β + x  ->  α(f x (g x))β  + S
```

The fired combinator comes out as a free atom. `K` releases what it drops as a
molecule of its own. `S` needs a second copy of `x`, its *reactant*, from the
soup, so an `S`-redex only counts as a redex of the chemistry when a copy of
`x` is present. `K`-reactions make expressions smaller, and `S`-reactions are
the only way an expression can take another molecule into its own body.

Expressions that cannot be reduced change at random instead, by *cleavage*
(`xy -> x + y`) or *condensation* (`x + y -> xy`).

### The reactor

Algorithm 1 of the paper repeats one step: draw an expression from the soup,
with probability proportional to how many copies there are. If it has a redex,
reduce one: the paper takes one at random among the first 100 found from the
outside in. If it has none, flip a coin. Either cleave it, or join it to the
last irreducible expression drawn before it. Reductions always come first, so
the paper treats them as *auto-catalysed*: an expression that can compute does
so without needing a partner.

### A worked example: a metabolic cycle

Here is `SII(SII)` in a soup that also contains one `SII`, run with Chemart's
reactor three times:

```python
import numpy as np
from chemart.chemistries.combinatory_chemistry import Pool, Reactor

pool = Pool(["SII(SII)", "SII"])
reactor = Reactor(pool, F=1, max_reductions=100)
rng = np.random.default_rng(0)
for _ in range(3):
    print(reactor.react(next(m for m in pool.count if "(" in m), rng))
print(dict(pool.count))
```

```
[('S', ('SII(SII)', 'SII'), ('I(SII)(I(SII))', 'S'))]
[('I', ('I(SII)(I(SII))',), ('I(SII)(SII)', 'I'))]
[('I', ('I(SII)(SII)',), ('SII(SII)', 'I'))]
{'S': 1, 'I': 2, 'SII(SII)': 1}
```

The first step is an `S`-reaction: `S` applies to `I`, `I` and the argument
`SII`, which it needs twice. It takes the free `SII` from the soup and
releases `S`. The two `I`-reactions that follow each release an `I`, and the
expression is back where it started. The free `SII` has been broken down into
its atoms `S`, `I`, `I`. This is the paper's Figure 2, written
`(AA) + A => (AA) + φ(A)` with `A = SII`, where φ(A) stands for the atoms of
`A`. The structure keeps itself by eating `SII`, but it does not multiply.

### Reactant assemblage

Structures like this need a steady supply of their food, and longer foods are
rare in a soup of mostly small molecules. Simulating a soup large enough to
supply them would be expensive. The paper's *reactant assemblage*
(Algorithm 2) imitates a larger soup instead: when an `S`-reaction needs a
reactant of at most `F` atoms that is missing, it is built on the spot from
free atoms, if they are there. `F = 1` switches this off.

## Using it

The default call runs 20,000 iterations of Algorithm 1 on 1,000 atoms (334
`I`, 333 `K`, 333 `S`) without reactant assemblage. That is a tenth of the
paper's atoms and a five-hundredth of its iterations, enough to see the soup
fill up but not to see structures emerge. The network lists every distinct
reaction that fired, with its count. `net.extras["reaction_kinds"]` gives the
type of each reaction in the same order (`I`, `K`, `S`, `cleave`, `condense`,
`assemblage`), and `extras["conservation"]` holds one conservation law per
atom type.

```python
from collections import Counter
net = chemart.generate_network("combinatory-chemistry", seed=1)
a = net.extras["analysis"]
print(a["diversity"])
print(a["reductions"][-3:])
print(Counter(net.extras["reaction_kinds"]))
```

```
[3, 27, 45, 61, 63, 74, 71, 70, 70, 69, 67, 63, 70, 73, 70, 73, 77, 76, 81, 77, 72]
[0.103, 0.115, 0.114]
Counter({'condense': 463, 'cleave': 182, 'K': 146, 'I': 108, 'S': 102})
```

`extras["analysis"]` samples the soup every `record_every` iterations:
`diversity` (distinct expressions), `mean_length` (atoms per molecule),
`reductions` (the share of iterations that reduced something, since the last
sample), `free_atoms` (free `S`, `K` and `I`) and `top_reactants` (the five
molecules most eaten by `S`-reactions since the last sample, the paper's
signal of emerging structures). `extras["final_state"]` is the final soup.

### The paper's runs

The paper runs 10,000 atoms for 10 million iterations. Chemart does about
100,000 iterations a second on small expressions, and slows down as long
expressions form. This F = 1 run took 160 seconds:

```python
net = chemart.generate_network("combinatory-chemistry", seed=0,
                               n_I=3334, n_K=3333, n_S=3333,
                               iterations=10_000_000, F=1, record_every=500_000)
a = net.extras["analysis"]
print(a["diversity"])
print(a["top_reactants"][4], a["top_reactants"][20])
```

```
[3, 357, 353, 367, 355, 339, 341, 334, 332, 336, 334, 334, 322, 298, 294, 330, 300, 295, 302, 304, 288]
{'I': 1306, 'K': 1092, 'S': 893, 'SII': 447, 'KI': 173} {'I': 1152, 'K': 1102, 'S': 643, 'SII': 635, 'KI': 145}
```

Set `F` between 3 and 20 for reactant assemblage, as in the paper's
Figures 4 and 5. Runs with large `F` build very long expressions and are
slower still. The metabolic cycles can be explored without a soup, with
`redexes` and `reduce_at` from `chemart.chemistries.combinatory_chemistry`, as
the tests do.

## Results

### Diversity, length and the share of reductions

Kruszewski and Mikolov ran the soup with 10,000 atoms and `F` from 1 to 20, ten
runs each (Figure 4). Diversity explodes during the first 200,000 or so
reactions and peaks at about 300 distinct expressions. It then declines,
slowly and steadily without assemblage, and faster the larger `F` is. The mean
length of expressions grows in step, as the atoms gather into fewer and longer
expressions. Reductions make up 10% to 35% of the reactions. At high `F` they
peak and then give way to cleavage and condensation, once the free atoms
needed to assemble reactants run out. The authors read this as the system
regulating its own balance of computation and chance.

Chemart reproduces the explosion and the F = 1 behaviour. In the 10-million
iteration run above, diversity reaches 357 within the first 500,000
iterations and falls slowly to 288. Mean length rises from 1.6 to 1.9 atoms
and reductions stay near 11%. The slow test
`test_paper_scale_diversity_and_sii_autopoiesis` runs the first 2 million
iterations of that run and checks the explosion (between 200 and 500 distinct
expressions after 500,000 iterations), the share of reductions (8% to 35%) and
the `SII` signature described below. How the decline speeds up with `F` has not been
checked against the figure.

### Emergent structures

To find structures, the paper counts which molecules `S`-reactions consume,
the only reactions that incorporate one molecule into another (Figure 5). The
structures give themselves away, "Tell me what you eat and I will tell you
what you are": an expression that eats `A` usually consists of copies of `A`
side by side.

- **Simple autopoietic.** `SII(SII)` keeps its form by eating `SII` (Figure
  2, the worked example above). Chance alone would make `SII` eaten less
  often than two-atom molecules such as `KK`, yet it is eaten more, because
  `SII(SII)` structures form. With `F = 3`, which makes `SII` easy to
  assemble, they become much more numerous. A rarer one, `AAA` with
  `A = SSK`, runs the cycle `(AAA) + 2A => (AAA) + A + φ(A)`. One of the two
  `A`s is released intact, which the authors describe as an emergent
  catalyst.
- **Recursive.** With `A = S(SI)I`, `AA` eats two `A`, adds one to itself and
  breaks down the other: `(AA) + 2A => A(AA) + φ(A)` (Figure 3), then
  `A(A(AA))`, and so on. A branching variant uses `A = S(SSI)K`. The paper
  sees these from `F = 4` on. When they appear, the simple autopoietic
  structures go extinct, because assemblage makes them compete for free atoms.
  Recursive structures survive shortages better: broken in two, a piece like
  `AA` still works.
- **Self-reproducing.** With `A = SI(S(SK)I)`, `AA` eats three `A` and ends
  as two copies of itself: `(AA) + 3A => 2(AA) + φ(A)` (Figure 1). It grows
  exponentially while food lasts, but in the paper's runs the recursive
  structures outcompete it, especially at `F = 8`.

Chemart checks every one of these cycles exactly. The tests search the
reduction pathways and confirm that each structure reaches the published end
state, eating exactly the published number of `A` and releasing the published
molecules: `test_fig2_sii_sii_is_simple_autopoietic`,
`test_fig3_tail_recursive_growth`, `test_fig1_self_reproduction` and
`test_ssk_cycle_releases_a_catalyst`. In the soup, the F = 1 run above shows
`SII` eaten more than any two-atom molecule (447 against 173 for `KI` at
2 million iterations, 635 against 145 at 10 million). `SII(SII)` is present at
the end, which the slow test also checks. That run also grew long tail-recursive
chains of `S(SI)I`, even though the paper reports recursive structures only
from `F = 4`. Whether self-reproducing `SI(S(SK)I)` structures appear, and the
extinctions of Figure 5, have not been reproduced in Chemart's runs.

### Conservation

Every reaction conserves the number of `S`, `K` and `I` atoms. That is how the
chemistry bounds itself without outside intervention. The test
`test_every_reaction_conserves_atoms` checks it on every reaction of a run and
on the final soup.

### Later work

The authors followed up with a journal version, Kruszewski and Mikolov (2022),
on self-reproducing metabolisms as recursive algorithms. Their public code has
since moved to a Gillespie simulation of the soup and does not contain the
2020 Algorithm 1 implemented here. Agüera y Arcas
et al. (2024), the source of [BFF](bff.md), cite it as an example of mass
conservation giving rise to ever-growing structures and transient
self-replicators.

## Further reading

- Kruszewski, G. & Mikolov, T. (2022). Emergence of Self-Reproducing
  Metabolisms as Recursive Algorithms in an Artificial Chemistry. *Artificial
  Life*, 27(3-4), 277-299.
- Fontana, W. & Buss, L. W. (1994). What would be conserved if "the tape were
  played twice"? *PNAS*, 91(2), 757-761. The AlChemy experiments the paper
  compares itself with.
