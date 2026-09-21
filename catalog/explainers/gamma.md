## Introduction

Gamma is a programming language in which a program is written as a set of
chemical reactions. Jean-Pierre Banâtre and Daniel Le Métayer first proposed it in a 1986 INRIA report and published it in *Science
of Computer Programming* (1990) and *Communications of the ACM* (1993). The
book dates it from the 1990 paper. According to the book, the name stands for
"general abstract multiset manipulation".

The problem it was built for is sequentiality. Most programs fix an order of
steps that the problem itself does not need. To find the largest of a set of
numbers, a textbook program walks through an array and keeps a running
maximum; the array and the loop force an order. Banâtre and Le Métayer wanted
a notation that states only what the problem requires. Their 2001 review
describes Gamma as "a formalism for the definition of programs without
artificial sequentiality", and a later paper by Banâtre, Fradet and Radenac
(2004) explains the word: "By artificial, we mean sequentiality only imposed by
the computation model and unrelated to the logic of the program."

The answer was a chemical picture. The data are molecules floating in a
solution, which formally is a *multiset* (a set that may hold several copies of
the same element). A program is a list of reactions, each saying "if some
elements satisfy this condition, replace them with these". Reactions fire in
any order, and at the same time wherever they involve different elements. The
program ends when no reaction can fire any more, and what is left in the
solution is the answer. The maximum of a set is then one reaction: take any two
numbers and keep the larger. Whatever order the comparisons happen in, only the
maximum survives, and many comparisons can run in parallel, like the matches of
a tournament.

Gamma is a *formalism*: a programming and specification language, not a model
of real chemistry. It has no rates, no energies and no time; the chemistry is a
metaphor for "unordered and parallel". Banâtre, Fradet and Radenac (2004)
call it, "to the best of our knowledge", the first "chemical model of
computation". Banzhaf and Yamamoto present it in chapter 9, "Rewriting
Systems" (book §9.2), between the λ-calculus and the [chemical abstract
machine (CHAM)](cham.md), which Berry and Boudol built on Gamma's picture and
extended with structured molecules and membranes. Its other neighbours in the
catalog take the same multiset rewriting in different directions: [P
systems](p-systems.md) nest multisets inside membranes and fire every rule that
can fire at each step, [MGS](mgs.md) gives the collection a spatial
neighbourhood, and the [high-order chemistry](high-order-chem.md) puts the
rules themselves into the vessel as molecules. The [prime number
chemistry](prime-number-chemistry.md) is a kinetic cousin: one division rule
run as a stochastic reactor rather than as a program.

## How it works

### A multiset and a list of reactions

A Gamma program works on one multiset. Its elements can be numbers or tuples, such as `(index, value)` pairs that encode a sequence. A reaction has
three parts, written in the book as

```
name : x1, ..., xn → A(x1, ..., xn) ⇐ R(x1, ..., xn)
```

- the *pattern* `x1, ..., xn` picks n elements from the multiset;
- the *reaction condition* `R`, after the `⇐`, must hold for them;
- the *action* `A` gives the elements that replace them: any number of new
  elements, possibly none.

The book's first example (eq. 9.7) is the maximum:

```
max : x, y → y ⇐ x ≤ y
```

Take any two elements `x` and `y`; if `x ≤ y`, remove both and put back `y`.
The net effect is that the smaller one disappears.

The *Gamma operator* runs a program: while some tuple of elements satisfies
some condition, pick one and apply its action. The choice is left open. If
several disjoint tuples qualify, they may react simultaneously. When no tuple
satisfies any condition, the multiset is *stable* (the γ-calculus papers say
*inert*), and that stable multiset is the result.

Two rules keep the model simple. The *locality principle* says a condition may
only look at the elements it matched: it cannot ask "is this the largest
element of the multiset?" or "how many elements are left?". And a program can
be built from stages: in the *sequential composition* `P2 ∘ P1`, `P1` runs until
it is stable and then `P2` runs on its result.

### A worked example: the maximum

Here are two runs of Chemart's execution of `max` on the multiset
`{3, 8, 1, 8, 4}`, with different random choices of which pair reacts:

```
run 1                         run 2
(4, 8) -> 8    {1, 3, 8, 8}   (3, 8) -> 8    {1, 4, 8, 8}
(1, 8) -> 8    {3, 8, 8}      (8, 8) -> 8    {1, 4, 8}
(3, 8) -> 8    {8, 8}         (1, 4) -> 4    {4, 8}
(8, 8) -> 8    {8}            (4, 8) -> 8    {8}
```

Each step consumes two elements and returns one, so five elements take four
steps. The order differs between the runs: in the second, the two copies of 8
meet early, and 1 is eliminated by 4 instead of by 8. The end is the same,
`{8}`, because no pair of distinct elements can remain in which one is not
larger than or equal to the other. That is the claim at the heart of Gamma: the
programmer states the condition under which the job is finished (no two
elements in the wrong relation) and the local step that makes progress, and
leaves the order to the machine.

### Programs that grow and shrink the multiset

Actions may create elements as well as remove them, so the set of possible
molecules is open. The review of Banâtre, Fradet and Le Métayer (2001) computes
the primes up to n in two stages. `iota` expands the interval `(2, n)` into the
numbers 2 to n, by splitting an interval `(x, y)` in half until it becomes a
single number; `rem` then removes every number that is a multiple of another:

```
iota:  (x, y) → (x, [(x + y)/2]), ([(x + y)/2] + 1, y)  ⇐  x ≠ y
       (x, y) → x                                         ⇐  x = y
rem:   x, y → y                                           ⇐  multiple(x, y)
```

`[ ]` is the integer part and `multiple(x, y)` is true when x is a multiple of
y. One run on `(2, 10)` in Chemart splits the interval in 17 steps into
`{2, 3, ..., 10}`, and then `rem` fires five times, `(9, 3)`, `(6, 3)`,
`(4, 2)`, `(10, 2)`, `(8, 2)`, leaving `{2, 3, 5, 7}`. The review points out
the three behaviours: the first reaction makes the multiset grow, the second
keeps its size and the third makes it shrink. The Fibonacci program has the
same expand-then-reduce shape: `dec1` breaks `n` into `fib(n)` ones and `add`
sums them.

### Which results are order independent

Not every Gamma program has a single answer. The *majority* program of
Banâtre, Fradet and Radenac (2004), `maj: x, y → {} ⇐ x ≠ y`, removes pairs of
different elements. If one value occurs in more than half of the multiset, only
copies of it remain, but how many copies depends on which pairs happened to
meet. The authors present it as correct all the same, since every possible
result contains only the majority element. A Gamma program can therefore have
several possible results; `max` is order independent because it has only one.

### Higher-order Gamma and the γ-calculus

In the original Gamma the program stands outside the multiset and cannot
change. Later work removed that separation. The γ-calculus of Banâtre, Fradet
and Radenac (2004) makes the reaction rules molecules too, so a rule can be
consumed, produced or passed around like data, and a solution can sit inside
another as a single molecule that may only be taken apart once it is inert.
Chemart implements only first-order Gamma.

### What Chemart builds from a program

A Gamma program has no reaction rates, so Chemart does not simulate it as a
kinetic system. It reports three things:

- **the network**: every concrete reaction the program can perform among the
  elements it can produce, found by applying the reactions to all tuples of
  known elements until nothing new appears (the *closure*). For `max` on
  `{3, 8, 1, 8, 4}`, that is the 10 reactions `x + y → y` with `x ≤ y` among
  1, 3, 4 and 8. Chemart names the integer k `n<k>`, so `n3 + n8 → n8` is
  "3 meets 8 and 8 remains", and `2 n8 → n8` is two copies of 8 becoming one.
- **one execution**: the Gamma operator run once, one reaction per step, the
  reaction chosen at random in proportion to the number of element tuples it
  can use. Parallel firing of disjoint tuples gives the same results as some
  order of single steps.
- **every result**: a search over all reaction orders, which lists every stable
  multiset the program can reach from its input. This is how Chemart shows that
  a program is deterministic in outcome, or not.

In the formal specification below, the molecules S are the elements, the
reactions R are the program's reactions applied to concrete elements, and the
reactor A is the Gamma operator.

## Using it

The default run is the book's `max` program (eq. 9.7) on `{3, 8, 1, 8, 4}`, an
input chosen by Chemart with a repeated maximum. Species are named `n<k>` for
an integer and by their text, such as `(1,5)`, for a tuple. The results are
in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["final_values"], a["stable"], a["steps"]   # ([8], True, 4)
a["results"], a["deterministic"]              # ([{'n8': 1}], True)
net.extras["stages"]                          # [['max: x, y -> y if x <= y']]
```

`final_values` is the stable multiset of the sampled execution, `steps` the
number of reactions it fired, and `stage_multisets` the multiset after each
stage of a composed program. `results` lists every stable multiset reachable
under any order, and `deterministic` says whether there is only one.
`net.extras["reaction_rules"]` gives, for each reaction of the network, the
name of the program reaction that produced it.

### The published programs

`program` selects one of the programs printed in the sources. Their example
inputs are Chemart's choices; the sources give the programs, not data. Each
line below is a run with `seed=1` and the default input:

```python
chemart.generate_network("gamma", seed=1, program="sort")
```

| `program` | default input | stable result | reachable results |
|---|---|---|---|
| `sort` | `(i, value)` pairs with values 8, 3, 6, 1, 9 | `(1,1) (2,3) (3,6) (4,8) (5,9)` | 1 |
| `primes` | `(2, 20)` | 2, 3, 5, 7, 11, 13, 17, 19 | not computed (budget) |
| `fibonacci` | 5 | 8 (after `{1 × 8}`) | 1 |
| `max-segment-sum` | values 3, −4, 5, −1, 2 | `(5, 2, 6)`: maximum segment sum 6 | 1 |
| `majority` | 1, 2, 1, 3, 1, 1, 2 | 1, 1, 1 | 2: `{1}` and `{1, 1, 1}` |
| `largest-prime` | 2 to 10 | 7 (after `{2, 3, 5, 7}`) | 1 |

In `max-segment-sum` each element is a triple `(i, x, s)`: position, value,
and the best sum of a segment ending at position i. The segment 5, −1, 2 ends
at position 5 and sums to 6. `largest-prime` is the example of the York
abstract by Banâtre, Fradet and Radenac, computing the largest prime below 10:
the primes stage runs to stability, then `max`.

Over seeds 0 to 4, the `majority` run ended with `{1}` three times and
`{1, 1, 1}` twice.

### Writing your own program

`program="custom"` with `rules` takes Gamma text, one reaction per line, in the
form `name: patterns -> expressions if condition`; a line `then` starts the next
stage of a sequential composition. The book's symbols `→`, `⇐`, `≤`, `≥` and
`≠` are accepted. The text is parsed by a small expression language, never run
as Python. A greatest-common-divisor program:

```python
net = chemart.generate_network("gamma", program="custom",
                               rules="gcd: x, y -> x - y, y if x > y",
                               multiset=[12, 18, 30])
a = net.extras["analysis"]
a["final_values"], a["deterministic"]        # ([6, 6, 6], True)
```

Each reaction replaces the larger of two numbers by their difference, so the
multiset ends with three copies of the gcd, 6.

The book's sort rule as printed (eq. 9.8, see the decisions below) can be run
the same way, and it sorts the wrong way round:

```python
chemart.generate_network("gamma", program="custom",
    rules="sort: (i, x), (j, y) -> (i, y), (j, x) if i >= j and x >= y",
    multiset=[[1, 8], [2, 3], [3, 6], [4, 1], [5, 9]]
).extras["analysis"]["final_values"]
# [[1, 9], [2, 8], [3, 6], [4, 3], [5, 1]]
```

### Budgets and slow settings

Three parameters in the table below bound the work. `max_species` caps the
closure: `fibonacci`'s `add` stage can build ever larger sums, so its network
is cut at 200 species and 10,103 reactions (`status="truncated"`) while the
execution and its result, 8, are unaffected. `max_states` caps the search over
all orders: from `(2, 20)` the `primes` program needs more than 5,000 states,
so `results` is `None`, and the run takes about 2 seconds; with
`max_states=0` the search is skipped and `(2, 30)` runs in about 0.25 seconds.
`max_steps` stops an execution that never becomes stable: the program
`x -> x + 1` on `{0}` with `max_steps=50` reports `stable=False` after 50
steps.

To re-inject input into a stable multiset (book §9.2), call the execution
function directly: `chemart.chemistries.gamma.run(parse_program(text),
Counter(values), rng, max_steps)` returns the final multiset first (then the
step count, whether it is stable, and the multiset after each stage); add new
elements to that multiset and run again.

## Results

Gamma is a language, so what has been done with it is programs, proofs,
implementations and other formalisms, not measurements. The main account is
the review by Banâtre, Fradet and Le Métayer, "Gamma and the chemical reaction
model: fifteen years after" (2001).

**Short programs for classic problems.** The papers express maximum, sorting,
primes, Fibonacci, the maximum segment sum and majority each in one to three
reactions. According to the review, the 1993 CACM paper (which Chemart could
not obtain) gives "a longer series of examples" from string processing, graphs and
geometry. The review also names five recurring reaction shapes, the *tropes*
(transmuter, reducer, optimiser, expander, selector), from which most programs
can be built; Fibonacci, for example, is an expander, a transmuter and a
reducer composed. Chemart reproduces the six small programs, each checked by
its tests:

- `max` (book eq. 9.7) ends with the single maximum. The tests draw random
  multisets, run them with several seeds and check that the order of
  comparisons never matters, and that the search over all orders finds only
  one result.
- `sort` ends with values increasing with the index. The tests check random
  sequences, and check on every reaction of the network the invariant used in
  the review's correctness proof: the set of indexes and the multiset of values
  do not change.
- `primes` (`rem ∘ iota`) turns `(2, n)` into the numbers 2 to n and then into
  the primes up to n. The tests check this up to 30, and that `rem` consumes
  exactly the composite numbers and never a prime.
- `fibonacci` (`add ∘ dec1`) expands n into `fib(n)` ones and sums them, with
  `fib(0) = fib(1) = 1`. The tests check n = 0 to 7.
- `max-segment-sum` leaves the triples carrying the maximum segment sum. The
  tests compare with a brute-force computation on random sequences.
- `majority` leaves only the majority element, and the tests check that two
  different multiplicities, 1 and 3, can be reached from the default input.

**Correctness by local reasoning.** The review shows how Gamma programs are
proved correct. The condition under which a program stops is simply the
negation of its reaction condition: for `sort`, the absence of any pair with
`i > j` and `x < y` is exactly the statement that the sequence is sorted. What must stay true
throughout can be checked on a single reaction, and termination follows from
an ordering on elements, using a result of Dershowitz and Manna on multiset
orderings. The same reasoning runs backwards: from a logical specification of
the primes, the review derives the reaction `rem` step by step. Its authors
stress that each part of the proof works by "reducing the global reasoning"
about the whole multiset "to a local reasoning (on the elements involved in a
single reaction)". Chemart checks the sort invariant, as above, but does not do proofs.

**Implementations.** Managing Gamma's parallelism efficiently "can be a
difficult task", the review says, and it lists finding the tuples that can react, applying reactions and
detecting termination as the major problems. It reports implementations on
the Connection Machine, an Intel iPSC/2, a MasPar MP-1 (with "a very good
speed-up"), a Sequent shared-memory multiprocessor, and a reconfigurable
hardware board (PRL-DEC Perle 1) with one circuit template per trope. In
C. Creveuil's 1991 thesis, optimising a naive Gamma interpreter led back to
well-known efficient sequential algorithms, shortest paths for instance. Two
applications to image processing, one reconstructing the vascular network of
the brain from two radiographs and one growing fractal models of biological
objects, found Gamma well suited but named "the lack of efficient general
purpose implementation" as a serious drawback. Chemart's executor is a simple
interpreter with no claim to efficiency.

**Extensions and descendants.** The review describes Gamma as "a source of
inspiration" in areas its authors did not expect:

- the [chemical abstract machine](cham.md) of Berry and Boudol (1992), which
  added membranes and airlocks to describe process calculi, and which the book
  also calls an extension of Gamma;
- composition operators (sequential and parallel) and their algebraic laws;
- higher-order Gamma (Le Métayer, 1994), where programs are configurations
  that can hold other programs;
- Structured Gamma (Fradet and Le Métayer, 1998), where the multiset carries
  relations between addresses described by graph grammars, and its
  descendants: shape types for C, used to check pointer structures, and
  descriptions of software architectures, including an industrial railway
  control case study.

The γ-calculus paper (Banâtre, Fradet and Radenac, 2004) looks for "a basic
calculus containing the very essence of the chemical paradigm". Its minimal
version, γ0, has four syntax rules, reactions without conditions that consume
one molecule at a time, and is still Turing-complete; it can encode the
λ-calculus. Adding reaction conditions and multi-molecule reactions gives the
γcn-calculus, which in the authors' words "closely models most of the existing
chemical programming models"; the paper compares it with Gamma, the CHAM,
higher-order Gamma, the hmm-calculus and P systems. A companion paper, cited
by the book, uses this higher-order model to specify self-organising
("autonomic") systems: as the authors' York abstract puts it, new molecules
perturb a stable solution and its rules react until it reaches a new stable
state.

**Re-injection.** The book notes that because the reactions stay active, new
input can be injected into a stable solution and the result recomputed, which
suits systems that must keep themselves up to date. Chemart's tests run `max`
to `{9}` from `{4, 9, 2}`, add 11 and 3, and check that it restarts and ends
at `{11}`.

**How the book places it.** In its chapter on chemical computing (§17.1),
the book contrasts formal calculi such as Gamma and P systems with
assembly-level chemistries: Gamma programs can be "very intuitive for humans",
but "tend to be far from chemical reality and to resist automatic
programming". It also uses Gamma's examples to show that a random order of
reactions can still give a deterministic outcome.

**What Chemart does not reproduce.** The γ-calculus and higher-order Gamma,
where rules are molecules and solutions nest, are not implemented: they need a
term-rewriting engine that the multiset model here does not provide. The York
largest-prime example, which in the γ-calculus waits for a sub-solution to
become inert, is offered as a two-stage sequential composition instead.
Structured Gamma, the tropes as a library, parallel composition and the
temporal logic of Gamma are not implemented either. The book's printed `sort`
(eq. 9.8) sorts in decreasing order, and its printed prime rule (eq. 9.9) is
the division rule of the [prime number chemistry](prime-number-chemistry.md)
rather than Gamma's `rem`; Chemart follows the 2001 review for both, as the
decisions explain. The 1990 and 1993 papers could not be obtained, so no
program is taken from them; programs such as `gcd` are
not offered by name but are easy to write as custom rules.

## Further reading

- Banâtre, J.-P. & Le Métayer, D. (1990). The Gamma model and its discipline
  of programming. *Science of Computer Programming* 15(1), 55–77 (book [58]).
- Banâtre, J.-P. & Le Métayer, D. (1993). Programming by multiset
  transformation. *Communications of the ACM* 36(1), 98–111.
- Banâtre, J.-P., Fradet, P., Giavitto, J.-L. & Michel, O. (2005).
  Higher-order chemical programming style. *Unconventional Programming
  Paradigms (UPP 2004)*, LNCS 3566 (book [55]).
- Banâtre, J.-P., Fradet, P. & Radenac, Y. (2004). Chemical specification of
  autonomic systems. *13th International Conference on Intelligent and Adaptive
  Systems and Software Engineering*, Nice (book [57]).
- Berry, G. & Boudol, G. (1992). The chemical abstract machine. *Theoretical
  Computer Science* 96, 217–248.
- Fradet, P. & Le Métayer, D. (1998). Structured Gamma. *Science of Computer
  Programming* 31, 263–289.
