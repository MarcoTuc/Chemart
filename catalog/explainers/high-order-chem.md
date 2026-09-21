## Introduction

In most artificial chemistries the molecules are data and the reaction rules
are part of the program: fixed, outside the soup, and the same for the whole
run. A *high-order* chemistry puts the rules into the soup as well. A rule is
itself a molecule, with a concentration, drawn at random like any other. To run
a different chemistry you change which rule molecules are in the vessel, not
the simulator.

This entry is Lidia Yamamoto's `HighOrderChem.py`, written in 2014 for
PyCellChemistry, the Python package that accompanies Banzhaf and Yamamoto's
book. The book presents it in its appendix on writing your own artificial
chemistry, as a second way to build a *constructive* chemistry (one whose set
of molecules grows as new ones are made). The first way is to copy an existing
program and edit its reaction code. The second is a "generic high-order
chemistry able to handle externally defined reaction rules", so that "we can
reuse the underlying algorithm for various different chemistries". The
authors call it "a simplified high-order chemistry in which the reaction rules
are also molecules in a multiset".

The worked example is the [prime number chemistry](prime-number-chemistry.md):
integers collide, and when one divides the other the larger is replaced by the
quotient. On top of `HighOrderChem` that whole chemistry is one rule molecule,
`divrule(m1, m2)`, placed in the rule multiset next to a few number molecules;
the book says the results "should be the same as the original NumberChem
implementation". Swap in rule molecules that edit city tours and the same
algorithm runs a version of the [molecular travelling salesman](molecular-tsp.md).

So this is a framework rather than a model of anything: a small reactor loop
that executes whatever rules it is given. What sets it apart from its
neighbours is where the rules live. In [Gamma](gamma.md), the multiset
rewriting language, the rules are the program, fixed outside the multiset; the
book notes that in the original Gamma they "are not explicitly encoded as molecules; therefore, they cannot be
rewritten". Here the rules are in a multiset of their own, so their numbers
matter and several rules compete for the same data. But this entry does not go
the whole way. Rules never react with rules, so no rule is ever rewritten. The
book is explicit that making `HighOrderChem` "truly high-order" would need
rules that operate on rules, and that the hard part would be designing a rule
set in which rewritten rules stay functional. Chemistries in which programs do
act on programs are elsewhere in the catalog, for example
[AlChemy](alchemy.md), where every molecule is a lambda-calculus expression
and a collision applies one molecule to another.

## How it works

### Two multisets and one loop

The vessel holds two *multisets* (bags in which the same item can occur many
times): one of **data molecules** and one of **rule molecules**. In the
reference code a rule molecule is a string that looks like a function call,
such as `divrule(m1, m2)`. The names inside the parentheses are its **binding
sites**: the number of data molecules the rule needs at once. `divrule(m1, m2)`
has two; a rule written `fold(m)` would have one.

One iteration of the reactor (the book's figure 3) is:

1. draw one rule molecule at random from the rule multiset;
2. count its binding sites, `nsites`;
3. if the data multiset holds at least `nsites` molecules, take `nsites` of
   them out at random, apply the rule to them, and put the molecules it returns
   (the products) into the data multiset;
4. put the rule molecule back.

Because every rule goes back after it fires, a rule is a *catalyst*: it takes
part in the reaction and comes out unchanged. The rule multiset therefore never
changes. What a rule returns replaces what it took, so the number of data
molecules changes only when a rule returns more or fewer molecules than it
bound. A rule that returns its educts (the molecules it took in) unchanged has
made an *elastic collision*: nothing happened. A reaction whose products differ
from its educts is *effective*.

Since a rule is drawn with probability proportional to how many copies of it
are in the vessel, rule multiplicities set how often each rule acts. That is
how the machines of the molecular TSP are weighted: 100 copies each of three
single-tour machines and one copy of the two-tour recombination machine.

### A worked example: division

`divrule` takes two integers. If the larger is a multiple of the smaller, the
larger is replaced by the quotient; otherwise both come back unchanged. The
first reaction of the default run (printed below, under *Using it in Chemart*) is

```
rule:divrule + n261 + n87 -> rule:divrule + n3 + n87  (x1)
```

The rule drew the numbers 261 and 87. Since 261 = 3 × 87, 261 becomes 3, and
87, the divisor, comes back unchanged. The rule appears on both sides. `(x1)`
means this exact reaction fired once in the run. A prime can only ever be the
divisor, never the number divided, so primes are never consumed, and over time
the soup fills with them.

With the three numbers 12, 2 and 3, and taking every possible collision rather
than a random run (`method="closure"`, below), the whole reaction network is
seven reactions (`net.to_text()`):

```
rule:divrule + n2 + n12 -> rule:divrule + n2 + n6
rule:divrule + n3 + n12 -> rule:divrule + n3 + n4
rule:divrule + n6 + n2 -> rule:divrule + n3 + n2
rule:divrule + n6 + n3 -> rule:divrule + n2 + n3
rule:divrule + n6 + n12 -> rule:divrule + n6 + n2
rule:divrule + n4 + n2 -> rule:divrule + 2 n2
rule:divrule + n4 + n12 -> rule:divrule + n4 + n3
```

`2 n2` means two molecules of 2. Every composite eventually breaks down into 2s and 3s.

### What a rule can be in Chemart

The reference code turns the rule string into a line of Python and runs it with
Python's `exec`, so a rule can be any Python function. Chemart never executes a
string. A rule is one of two things:

- a **named rule** from a built-in registry, ported from the reference
  examples: `divrule` (from `NumberChemHO.py` and the book) and the four tour
  machines of `MolecularTSP.py`, `exchangeMachine` (swap two cities),
  `cutMachine` (move a segment of the tour), `invertMachine` (move and reverse
  a segment) and `recombinationMachine` (graft a segment of one tour into
  another and keep the two best of parents and child). The single-tour machines
  release the new tour only if it is shorter than the old one.
- an **expression rule** in a small language that Chemart parses and
  interprets: `[name:] x, y -> expr, expr, ... [if condition]`. The variables
  on the left are the binding sites, the expressions on the right are the
  products (none means the educts are destroyed), and the optional condition
  decides whether the rule fires. Expressions use integers, `+ - * / %`
  (`/` is integer division), comparisons, `and or not`, and `min`, `max`,
  `abs`. If the condition is false or an expression does not make sense (a
  division by zero, arithmetic on a list), the rule returns its educts: an
  elastic collision.

Data molecules are integers or lists of integers. An integer `k` is the species
`n<k>`; a list, such as a tour, is written as compact JSON, `[0,3,1,2]`. A rule
molecule is the species `rule:<name>`.

The formal specification below restates this: S is the two multisets, R is the
scheme "rule plus as many data molecules as it has binding sites goes to rule
plus products", and A is the loop above.

## Using it

The default run is `NumberChemHO.py`, the reference re-implementation of the
prime number chemistry: four copies of `divrule`, 100 numbers drawn uniformly
from 2 to 1000, and 10,000 iterations. It takes a few seconds. Of the 10,000
iterations, 180 were effective collisions, which gave the 154 distinct
reactions of the summary; the rest were elastic. The run's record is in
`net.extras`:

```python
a = net.extras["analysis"]
a["effective_collisions"], a["idle_draws"]   # (180, 0)
a["rule_draws"]                              # {'rule:divrule': 10000}
pf = a["prime_fraction"]                     # fraction of primes, every 100 iterations
pf[0], pf[10], pf[50], pf[-1]                # (0.14, 0.28, 0.92, 0.99)
net.extras["rules"]                          # {'rule:divrule': 'divrule(m1, m2)'}
```

`prime_fraction` is measured once per *generation* of M iterations (here 100)
and is only recorded when `divrule` is the sole rule. `idle_draws` counts
iterations in which the drawn rule found too few data molecules to bind.
`rule_draws` says how often each rule was drawn. `net.initial_state` and
`net.extras["final_state"]` give the counts of every molecule at the start and
end, rule molecules included; both hold 4 copies of `rule:divrule` and 100
numbers, since division returns two numbers for two. The network holds each
distinct effective reaction once, with its firing count.

**Every reaction, not one run.** `method="closure"` ignores randomness and
lists every reaction reachable from the distinct starting numbers, as in the
12, 2, 3 example above. It is a Chemart addition and only works for
deterministic rules, so the tour machines are refused.

```python
net = chemart.generate_network("high-order-chem", method="closure", data=[12, 2, 3])
print(net.summary())   # high-order-chem: 6 species, 7 reactions, status=complete
```

`data` sets the starting data multiset explicitly; `max_species` caps the
closure, which then reports `status="truncated"`.

**The same chemistry written as an expression rule.** This rule reproduces
`divrule`, and gives the same seven reactions under the rule name `rule:div`:

```python
expr = "div: x, y -> max(x, y) / min(x, y), min(x, y) if x != y and max(x, y) % min(x, y) == 0"
net = chemart.generate_network("high-order-chem", method="closure", data=[12, 2, 3],
                               rules={expr: 1})
```

**Rules competing for substrate.** Add a rule that destroys whatever it binds,
`eat: x ->`, with one copy of each rule and 3,000 iterations:

```python
kept  = chemart.generate_network("high-order-chem", rules={"divrule": 1}, seed=5, iterations=3000)
eaten = chemart.generate_network("high-order-chem", rules={"divrule": 1, "eat: x ->": 1},
                                 seed=5, iterations=3000)
```

Alone, `divrule` makes 151 effective collisions. With `eat` beside it the two
rules are drawn about equally (1,510 and 1,490 times), `eat` empties the soup,
`divrule` manages only 4 divisions, and 2,807 iterations are idle because no
data is left.

**Rule multiplicity sets the odds.** With `rules={"a: x -> x": 3, "b: x -> x": 1}`
and 4,000 iterations (seed 2), `a` was drawn 2,993 times and `b` 1,007, close to
3:1. Both rules return what they take, so the network has no reactions.

**The molecular TSP machines.** The settings of `MolecularTSP.py`: nine random
tours of ten cities on a ring, and its machine weights.

```python
rules = {"exchangeMachine": 100, "cutMachine": 100, "invertMachine": 100,
         "recombinationMachine": 1}
net = chemart.generate_network("high-order-chem", rules=rules, init="tours",
                               M=9, cities=10, iterations=3000, seed=3)
```

This takes a few seconds. The nine starting tours have lengths from 110.2 to
184.6; after 3,000 iterations all nine have length 61.8, the perimeter of the
ten-city polygon, which is the shortest tour on a ring. (The nine are different
lists, the same round trip started at different cities or run backwards.) The
reference runs up to 1,000 generations of 9 × 100 / 301, rounded up to 3,
iterations each; `iterations` here is just the total. For the full molecular
TSP, with random city layouts, use the dedicated
[molecular-tsp](molecular-tsp.md) entry.

The remaining parameters are in the table below: `M`, `minn` and `maxn` for the
random numbers, `cities` for the tours.

## Results

`HighOrderChem.py` is a teaching and tooling example, not a research model, and
the book reports no experiments with it. Its claims are about what the design
allows, and each is checked in Chemart's tests.

**Divrule reproduces the prime number chemistry.** The book's one concrete
claim is that `divrule` in `HighOrderChem` gives "the same" results as the
original NumberChem program. Chemart's tests check this two ways. With
`method="closure"`, the network from the numbers 12, 2, 3, and from a random
soup of 100 numbers, has exactly the reactions of the
[prime-number-chemistry](prime-number-chemistry.md) entry once the rule
molecule is removed. In a slower test over six random runs of the default
setting, only composites are ever consumed, the prime fraction never falls,
its final value averages above 0.95, and it is within 0.05 of the prime number
chemistry's own runs. In the default run above the fraction rises from 0.14 to
0.99.

**Rules are catalysts.** Because every rule is put back, each reaction has its
rule on both sides and the rule multiset is the same at the end as at the
start. The tests check both on the default run, and that the loop binds exactly
as many molecules as the drawn rule has binding sites (idling when there are
too few).

**Rules compete.** The book points out that with rules as molecules "we could
also look at how rules compete". Each iteration draws one rule, with
probability proportional to its multiplicity; the tests check the 3:1 ratio
above and that a substrate-destroying rule starves `divrule` of divisions.

**One algorithm, many chemistries.** Changing the rule multiset changes the
chemistry without touching the loop. The tests run the four tour machines of
`MolecularTSP.py` and check that the single-tour machines only ever release a
shorter tour, and that the best and the mean tour length both fall. The book
also suggests adding arithmetic rules to build an algebraic chemistry like
[ACGP](acgp.md); Chemart's expression rules make this possible, but there is
no published experiment to reproduce.

**What is not there.** The book's own caveat applies: this is not yet truly
high-order, since rules never take rules as educts and no rule is rewritten.
Its suggested next step, typed binding sites that would let a rule draw another
rule as an argument, is not implemented, and neither is the open research
question it raises, a rule set whose rewritten rules remain functional. The
other differences from the reference are deliberate: no string is executed
(unknown rule names are rejected rather than run, which the tests also check),
rules are limited to the built-in registry and the expression language, and
built-in rules collide elastically on molecules outside their domain where the
reference would crash. The implementation decisions below list them.

## Further reading

- Banâtre, J.-P., Fradet, P., Giavitto, J.-L. & Michel, O. (2005).
  Higher-order chemical programming style. In *Unconventional Programming
  Paradigms (UPP 2004)*, Lecture Notes in Computer Science 3566. Springer.
  The γ-calculus, which the book (§9.2, bibliography [55]) cites as a
  higher-order extension of Gamma in which programs are themselves molecules.
