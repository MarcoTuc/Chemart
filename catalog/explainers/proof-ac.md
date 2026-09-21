## Introduction

This entry is an automatic theorem prover built as an artificial chemistry.
The molecules are logical statements and a reaction is one step of logical
inference. Put the axioms of a theory into a reactor together with the
negation of the statement you want to prove, and let the molecules collide at
random. If a collision ever produces a contradiction, the statement is proved,
and the proof is the chain of reactions that built the contradiction.

The idea goes back to Fontana and Buss (1996). Their lambda-calculus chemistry,
[AlChemy](alchemy.md), has molecules that are functions, and in lambda calculus
a proof of a statement can be read as a term of the matching type (the
Curry-Howard correspondence, Howard 1980). Fontana and Buss suggested that
mathematical truth could therefore be mapped onto chemistry: a statement is
true if the molecule that stands for its proof can be built. Banzhaf and
Yamamoto present this in their chapter on applications (book §16.5), and then
turn to a concrete system that does not use lambda calculus at all: RESAC, a
"resolution-based artificial chemistry" by Jens Busch, described in
Busch and Banzhaf (2003) and in Busch's PhD thesis (Dortmund, 2004). In RESAC
the molecules are *clauses* of first-order logic and the only reaction is
*resolution*, the inference rule most automatic theorem provers are built on.
Chemart implements RESAC.

A small example shows the flavour. The clause "every dog howls" and the
clause "the pet is a cat or a dog" meet; because one says "dog" and the other
"not a dog, or it howls", they react and give "the pet is a cat or it howls".
Both parents stay true, and so does the product: resolution never produces
anything that does not follow from its inputs. A conventional prover chooses
carefully which pairs of clauses to combine next. RESAC does not choose: pairs
meet at random, many reactions are dead ends, and the prover relies on the
number of collisions instead of on a heuristic that picks the next step.
Busch's thesis says RESAC was not designed as a theorem prover, and it argues
that such a prover could be spread over very many parallel units, the
advantage the book also stresses.

Among the catalog's constructive chemistries, this is the one whose reactions
have a meaning fixed from outside: a product is not just a new string or a new
function, as in [AlChemy](alchemy.md) or
[combinator chemistry](combinator-chemistry.md), but a logical consequence of
its parents. Like the [Chemical Casting Model](ccm.md) and the
[molecular travelling salesman](molecular-tsp.md), it uses the chemical
metaphor to solve a problem, here the problem of finding a proof.

## How it works

### Clauses, literals and the empty clause

A *literal* is a basic statement or its negation: `dog(s)` ("s is a dog") or
`~dog(s)` ("s is not a dog"). Its arguments are *terms*: constants such as
`john`, variables such as `X`, or functions of terms such as `f(X,Y)`. Chemart
writes variables with an upper-case initial and everything else in lower case,
as Prolog does. A *clause* is an "or" of literals, and its variables are read
as "for all". So `~dog(X), howls(X)` says "for every X, X is not a dog or X
howls", which is "all dogs howl". Any statement of first-order logic can be
turned into a set of clauses; an existential "some" becomes a new constant, a
*Skolem constant*, that names the thing that exists.

The *empty clause*, written `[]`, is an "or" of nothing, which is false.
Deriving it means the clauses you started from contradict each other. This is
how a resolution prover proves a theorem: it adds the *negation* of the theorem
to the axioms and derives `[]`, a proof by contradiction called a
*refutation*.

### The reaction: binary resolution

Two clauses react if one contains a literal and the other contains the
negation of a literal that can be made identical to it by substituting terms
for variables. Finding such a substitution is called *unification*, and the
least specific one that works is the *most general unifier*. For `dog(s)` and `~dog(X)` it is `X = s`. The product, the
*resolvent*, is everything else in the two clauses with the substitution
applied:

```
cat(s), dog(s)   +   howls(X), ~dog(X)    ->    cat(s), howls(s)
```

If no literal pair unifies, the collision is *elastic*: nothing happens. Two
clauses can also have several resolvents, one for each pair of literals that
unifies. Before resolving, the variables of the two clauses are renamed apart,
since the `X` of one clause has nothing to do with the `X` of another.

Three housekeeping rules come from the thesis. Identical literals in a clause
merge, because a clause is a set. A clause containing both a literal and its
negation (a *tautology*) is always true and useless, so it is never produced.
And *factoring*: if two literals of the same sign in one clause unify, the
clause may be narrowed to the instance where they coincide. Without factoring,
resolution cannot refute some contradictory sets, such as `p(X), p(Y)` against
`~p(U), ~p(V)`. Chemart folds factoring into the reaction: the resolvents of
two clauses include those of their factors.

Species are clauses up to renaming of variables: `~has(X,Y), mice(Y)` and
`mice(B), ~has(A,B)` are one species, written with the literals sorted and the
variables renamed `X1, X2, ...` in order, as in `howls(X1),~dog(X1)`.

### A worked example: John, his pet and the mice

The default problem is a puzzle from Gini's lecture notes (1995), which Busch
uses as his first example (thesis table 6.1). All dogs howl at night. Anyone
who has any cats has no mice. Light sleepers have nothing that howls at night.
John has a cat or a dog. To prove: if John is a light sleeper, he has no mice.
Negating the goal gives "John is a light sleeper and has some mice `m`", and
the pet is named by the Skolem constant `s`. The eight start clauses are the
five axioms

```
howls(X1),~dog(X1)                            all dogs howl
~cat(X1),~has(X2,X1),~has(X2,X3),~mice(X3)    cat owners have no mice
~has(X1,X2),~howls(X2),~lightsleep(X1)        light sleepers own nothing that howls
has(john,s)                                   John has s
cat(s),dog(s)                                 s is a cat or a dog
```

and the three clauses of the negated goal: `lightsleep(john)`,
`has(john,m)` and `mice(m)`.

Chemart's default run finds this refutation (`net.extras["analysis"]["proof"]`
of the default network):

```
1. cat(s),dog(s)       + howls(X1),~dog(X1)                  -> cat(s),howls(s)
2. has(john,s)         + ~has(X1,X2),~howls(X2),~lightsleep(X1)
                                                             -> ~howls(s),~lightsleep(john)
3. cat(s),howls(s)     + ~howls(s),~lightsleep(john)         -> cat(s),~lightsleep(john)
4. cat(s),~lightsleep(john) + lightsleep(john)               -> cat(s)
5. has(john,s)         + ~cat(X1),~has(X2,X1),~has(X2,X3),~mice(X3)
                                                             -> ~cat(s),~has(john,X1),~mice(X1)
6. has(john,m)         + ~cat(s),~has(john,X1),~mice(X1)     -> ~cat(s),~mice(m)
7. mice(m)             + ~cat(s),~mice(m)                    -> ~cat(s)
8. cat(s)              + ~cat(s)                             -> []
```

Read in words: the pet is a cat or it howls (1); since John owns it, either it
does not howl or John is not a light sleeper (2); so it is a cat or John is no
light sleeper (3); John is a light sleeper, so the pet is a cat (4). The other
branch starts again from "John has `s`": if `s` is a cat, John has no mice (5);
but John has `m` (6) and `m` is a mouse (7), so `s` is not a cat. The two
conclusions contradict each other (8). Unification did real work in steps 1, 2
and 5, where the variables of a general rule were bound to `s`, `john` or `m`;
the resolvent is then a statement about those particular individuals.

### Two ways to run the chemistry

**The closure.** Chemart can compute every clause that resolution can reach,
level by level: level 0 is the start clauses, level 1 everything one
resolution away, and in general level k holds the clauses first produced by
resolving a level k−1 clause with any earlier clause. This *level saturation*
is the systematic search the thesis uses to measure the size of the search
space (thesis 6.1). Premises are kept, so each reaction reads `a + b -> a + b
+ r`: a theory only grows. For the puzzle the closure is finite (100 clauses,
402 reactions) and `[]` first appears at level 4. For theories with function
symbols the closure is usually infinite, so a species budget truncates it.

**The soup.** This is RESAC's reactor (thesis algorithms 3.1 to 3.3). The
reactor is a fixed number of slots, filled with equal numbers of copies of each
start clause; the number of copies per clause is the *multiplicity*. Each step
draws two molecules at random. If they can resolve, one resolvent is picked at
random and written into the reactor. Under *educt replacement* it overwrites
one of the two reactants (chosen by a coin flip), so the reaction is `s1 + s2
-> s_i + s3`; under *free replacement* it overwrites a random molecule anywhere
in the reactor. The reactor has no memory beyond its slots, so useful clauses
can be overwritten. The thesis also adds *inflow*, which it calls necessary in
theory because nobody knows in advance how often a clause will be used in a
proof (though large enough closed reactors find proofs too). Fresh copies of
start clauses overwrite random molecules, either after every elastic collision
(*elastic inflow*, which feeds the reactor exactly when it is unproductive) or
at a fixed rate. A run ends when the target, normally `[]`, is produced, or at
a time limit counted in collisions. The proof is then read off the history of
how the target was made; dead ends do not appear in it, as the book points out.

### Restriction strategies

Resolution provers usually forbid some pairs to keep the search small, and the
thesis offers two classic restrictions (thesis 1.3.3). Under *set of support*,
one parent must descend from the negated goal, so the search stays focused on
the goal. Under *negative resolution*, one parent must contain only negated
literals. Both keep the prover complete: if the start clauses are
contradictory, `[]` can still be derived. Busch names the combination of the two
as RESAC's preset. Clauses can also be declared *unstable* above a maximum
length, counted in symbols (predicate, function, constant and variable
occurrences), which stops the population from filling with ever longer
clauses.

## Using it

The default call above computes the closure of the puzzle, with no restriction
strategy: its 100 clauses and 402 reactions are everything resolution can
derive from the eight start clauses, and the proof is among them:

```python
a = net.extras["analysis"]
a["proved"], a["proof_length"], a["target_level"], a["longest_clause_in_proof"]
# (True, 8, 4, 3)
net.extras["goal"]         # ['lightsleep(john)', 'has(john,m)', 'mice(m)']
a["proof"][-1]             # {'reactants': ['cat(s)', '~cat(s)'], 'product': '[]'}
```

`net.extras["axioms"]` and `net.extras["goal"]` list the start clauses as
species ids. `a["proof"]` is the refutation, premises before conclusions, each
step a dict of two reactants and a product; `a["levels"]` gives every clause's
level, and `a["longest_clause_in_proof"]` counts literals. Each species'
`structure` writes its clause with ` | ` between literals. In the reactions, the two
premises appear on both sides, which is why the network reports catalysts: a
premise is used but not consumed.

**The stochastic prover (thesis fig. 6.3).** `method="soup"` runs RESAC's
reactor with the settings the thesis gives for its run-time experiment,
multiplicity 20 (160 molecules) and elastic inflow, plus Chemart's defaults of
educt replacement and a limit of 20,000 collisions.

```python
net = chemart.generate_network("proof-ac", method="soup", seed=1)
a = net.extras["analysis"]
a["proved"], a["collisions_to_proof"], a["productive_collisions"], a["inflows"]
# (True, 1077, 261, 816)
times = [chemart.generate_network("proof-ac", method="soup", seed=s)
         .extras["analysis"]["collisions_to_proof"] for s in range(20)]
sorted(times)
# [873, 1077, 1516, 1595, 1876, 2015, 2200, 2515, 3391, 3815, 3928, 4187,
#  4198, 4413, 5433, 5516, 6965, 8765, 12884, 14644]
```

Only about a quarter of the collisions were productive; each of the 816
elastic ones let a start clause flow in. The soup's network lists the
reactions that actually fired, with counts, and `net.extras["final_state"]`
holds the reactor's last contents. `a["generations"]` is the number of
collisions divided by the reactor size. The proofs differ from run to run: this
run's refutation has 7 steps and a 4-literal clause, other seeds give 8 steps.
Twenty seeds take about two seconds.

**Reactor size (thesis figs. 6.1 and 6.2).** Vary `multiplicity` with a fixed
time limit:

```python
for m in (1, 5, 20, 100):
    runs = [chemart.generate_network("proof-ac", method="soup", seed=s, multiplicity=m,
                                     max_collisions=8500).extras["analysis"]
            for s in range(12)]
    print(m, sum(r["proved"] for r in runs))
# 1 3
# 5 10
# 20 11
# 100 7
```

With one copy of each clause, 3 of 12 runs find the proof in 8,500 collisions;
with 20 copies, 11 of 12; with 100 copies the success rate falls again, because
a larger reactor needs more collisions.

**Your own theory.** Set `problem="custom"` and give the axioms in `clauses`
and the negated theorem in `goal`, with `;` between clauses and `,` between
literals. To prove that Socrates is mortal:

```python
net = chemart.generate_network("proof-ac", problem="custom",
                               clauses="~man(X), mortal(X); man(socrates)",
                               goal="~mortal(socrates)")
print(net.to_text())
```
```
man(socrates) + mortal(X1),~man(X1) -> man(socrates) + mortal(X1),~man(X1) + mortal(socrates)
mortal(X1),~man(X1) + ~mortal(socrates) -> mortal(X1),~man(X1) + ~mortal(socrates) + ~man(socrates)
mortal(socrates) + ~mortal(socrates) -> mortal(socrates) + ~mortal(socrates) + []
man(socrates) + ~man(socrates) -> man(socrates) + ~man(socrates) + []
```

The closure finds two refutations, one forward from the fact and one backward
from the goal. `target` can name any clause instead of `[]`, to ask whether it
is derivable. The `strategy` parameter applies the restrictions described
above; on the puzzle they shrink the closure from 100 clauses to 70
(`set-of-support`), 88 (`negative`) or 52 (`negative-set-of-support`), and all
three still prove it.

The remaining settings (`factoring`, `max_length`, `replacement`,
`inflow_rate` and the rest) are described in the parameter table below.

**Slow settings.** `problem="group-right-inverse"` has function symbols, so its
closure never ends: the default budget of 2,000 species fills up in about
4 s, partway through level 2, without a proof. Its soup at the default scale runs all 20,000 collisions
without a proof and takes about 40 s.

## Results

**The puzzle, proved by collisions (thesis table 6.1).** Busch's first
experiment gives RESAC the eight clauses of the puzzle. The proof it found has
seven steps, and the thesis remarks that its clauses first grow longer before
they shrink to the empty clause: it passes through a clause of five literals,
`~has(X,s), ~has(Y,Z), ~has(Y,s), ~lightsleep(X), ~mice(Z)`. The thesis also
notes that clause s5 (`cat(s), dog(s)`, two positive literals) puts the problem
outside Horn clauses (clauses with at most one positive literal), the
restricted form Prolog works with. Chemart's tests
check the published proof step by step: each of the seven clauses is a
resolvent of the two premises the thesis names, all are in the default
closure, and their lengths run 2, 3, 5, 4, 3, 1, 0. The proof Chemart itself
reports, shown above, is a different one of lowest level. The published proof
starts by resolving s5 with s1, a pair that neither set of support nor negative
resolution allows, so it cannot come from RESAC's preset strategy; this is why
Chemart's default is unrestricted.

**Resolution is sound and refutation-complete.** These two properties, which
the thesis states for the calculus (thesis 1.2.2), are what make the chemistry
a prover: every product follows logically from its parents, and a set of
clauses is contradictory exactly when `[]` can be derived. Chemart tests both
on 60 random propositional clause sets against truth tables: every product is
true in every assignment that makes its parents true, and the closure contains
`[]` exactly for the unsatisfiable sets. Further tests check that the puzzle's
five axioms alone never produce `[]` (only the negated goal makes them
contradictory), that `p(X), p(Y)` against `~p(U), ~p(V)` needs factoring, and
that the restriction strategies still refute the puzzle.

**Proof times vary strongly (thesis fig. 6.3).** Measured in collisions, the
time a run needs to find the proof differs widely between runs. Busch ran the
puzzle 30,000 times with multiplicity 20 and elastic inflow and sorted the
times into classes (1-8,500 collisions, 8,501-20,000, 20,001-35,000,
35,001-50,000 and 50,001-80,000). Numbering the classes from zero, he found a
Poisson distribution with mean 0.3601, so most runs fall in the first class.
Chemart's test runs 20 seeds with the same settings and requires at least 18
proofs within 20,000 collisions and at least 12 within 8,500. The run above
gives 20 and 17.

**Reactor size (thesis figs. 6.1 and 6.2).** In nearly 10,000 runs, with 50
runs per combination of multiplicity (up to 7,000 copies per clause) and time
limit (up to 80,000 collisions), Busch found that a larger reactor needs
significantly more collisions, while a reactor that is too small can prevent
the proof entirely: the more complex a proof, the more intermediate clauses the
reactor must hold. Chemart's test reproduces the qualitative claim: with one copy per clause, fewer of 12 seeds find the proof
within 8,500 collisions than with 20. The success-rate surface of fig. 6.1 is
not reproduced.

**Group theory: right inverses (thesis 6.6).** The second problem, from
Loveland (1978), asks for a proof that in every group each element has a right
inverse, from six clauses with function symbols. RESAC finds a proof of length
8. With 14,000 copies per clause, Busch compared 100 runs in a closed reactor
with 100 under elastic inflow (fig. 6.4). In the closed reactor the start
clauses are soon replaced by more specialised, less reactive clauses, and
productivity falls; inflow keeps productivity up at the cost of lower
diversity. The open reactors consistently found the proof faster. Chemart
provides the problem but does not reproduce this: level saturation with 2,000
species does not finish level 2, and a soup of multiplicity 20 finds no refutation in 20,000
collisions, far below the thesis's scale. A test checks only that the closure
is truncated, as it must be for an infinite one.

**Harder problems (thesis 6.7 and 6.8), not in Chemart.** The thesis also
tried the Burnside problem, which asks for a proof that in every group
where `x³ = e` for all `x` a certain commutator identity holds; a resolution
proof by Robinson and Wos (1969) takes over 138 intermediate steps. RESAC found
no proof, even in a chain of six reactors. The thesis blames the six equality
axioms, which react with almost every clause, and concludes that without
problem-specific heuristics RESAC is inferior to specialised provers here. On
Schubert's Steamroller, a puzzle about animals eating each other that cannot be
written in Horn clauses, it succeeded: from 26 clauses, with 2,600 molecules,
RESAC derived the contradiction after 69 resolutions in generation 321,
showing that foxes eat a grain-eating animal, namely birds. Neither problem is
built into Chemart, and neither has been tried through the `custom` option.

## Further reading

- Gini, M. (1995). Lecture notes, http://www-users.cs.umn.edu/~gini/5511/
  (as cited in Busch's thesis). The source of the puzzle of table 6.1.
- Loveland, D. W. (1978). *Automated Theorem Proving: A Logical Basis*.
  North-Holland. The source of the group-theory problem of table 6.2.
- Robinson, G. & Wos, L. (1969). Paramodulation and theorem-proving in
  first-order theories with equality. In B. Meltzer & D. Michie (eds.),
  *Machine Intelligence 4*, 135–150. Edinburgh University Press.
  The Burnside proof cited in the thesis.
