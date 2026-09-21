## Introduction

P systems, also called membrane systems, are a model of computation proposed by
Gheorghe Păun in a 1998 Turku Centre for Computer Science report, published in
2000 as "Computing with membranes". The research area that grew from it is
*membrane computing*, which in Păun's words aims "to abstract computing ideas
and models from the structure and the functioning of living cells" (Păun
2006). The picture is a cell: an outer membrane encloses it, inner membranes
enclose the nucleus and the organelles, and each compartment holds its own mix
of chemicals reacting there. A P system keeps that layout, but reduces the
chemicals to symbols and their reactions to rewriting rules.

Each compartment is a bag of symbols with its own list of rules, such as
`ab → aac` (one `a` and one `b` become two `a` and a `c`). A rule can keep its
products where they are, push them out through the membrane or into an inner
compartment, and it can dissolve its own membrane. At every tick of a global
clock, every compartment applies as many rules as it can at once. When no rule
can fire anywhere, the computation has halted, and the answer is read by
counting the symbols in a chosen compartment. The book's example decides
whether a number `n` is a multiple of `k`: the system ends with one symbol in
its output compartment (no) or with none (yes).

P systems are a *formalism*, not a single simulation model. Păun (2006)
calls membrane computing "a framework for devising compartmentalized models",
and the literature holds a very large number of variants. Chemart implements
the basic one, the *transition P system*, with symbol objects, targets for the
products, membrane dissolution and rule priorities. The field did not start as
a way to model biology: Păun lists computer-science goals, models "as elegant
(minimalistic) as possible, as powerful as possible (in comparison with Turing
machines ...), and as efficient as possible". Biological modelling came later;
Banzhaf and Yamamoto note that P systems "have also been extensively used to
model biology" (book §18.3.2).

The book presents P systems in chapter 9, "Rewriting Systems" (§9.5), among
other chemistries built on multiset rewriting. [Gamma](gamma.md) rewrites a
single bag with no compartments; what Gamma lacks, says Păun (2006), is
"distributivity". [The Chemical Abstract Machine](cham.md), which Păun calls
"the direct ancestor of membrane systems", has membranes too, but they are
sub-solutions in an algebra of concurrent processes, not cell compartments.
[ARMS](arms.md) is, in the book's words, "an example of a P system specialized
for Artificial Life experiments" (Chemart's ARMS entry leaves out the
membranes). In [brane calculi](brane-calculi.md) the membranes themselves act,
fusing and engulfing, whereas in a P system they are containers and the objects
inside do the work. What P systems add to the catalog is the nested tree of
compartments, with transport between them written into the rules, and
maximally parallel execution.

## How it works

### Membranes, regions and objects

A P system has a tree of membranes. The outermost is the *skin*; outside it is
the *environment*. A membrane with nothing inside it is *elementary*. Each
membrane encloses a *region*, the space between it and the membranes directly
inside, and both share a label. The tree is written with labelled brackets:
`[1 [2 ]2 [3 ]3 ]1` is a skin, labelled 1, holding two elementary membranes,
2 and 3.

Each region holds a *multiset* of objects: a set in which copies count. Objects
are plain symbols with no internal structure, and a multiset is written as a
string, so `a^5 b^2 c^6` means five `a`, two `b` and six `c`. Chemart writes it
in reaction syntax, `5 a + 2 b + 6 c`, and names each species by object and
region: `a@2` is an `a` in region 2, `a@env` an `a` in the environment.

### Rules and maximal parallelism

Each region has its own rules. A rule `u → v` consumes the multiset `u` and
produces `v`, like a chemical reaction: Păun's `aab → abcc` is
`2 a + b → a + b + 2 c`. A rule of the form `ca → cv`, where the object `c`
comes out unchanged, is *catalytic*, and `c` is a catalyst.

In one step a region does not apply one rule once. It chooses, at random, a
*maximal* multiset of rule applications: each rule gets a number of uses, the
objects present must suffice for all of them, and no further use can be added
from the objects left over. Păun's example (2006, §4) has the rules
`aab → abcc` and `bb → aac` on `a^5 b^2 c^6`. They compete for the two `b`.
Either the first rule is used twice, giving `a^3 b^2 c^10`, or the second once,
giving `a^7 c^7`. Using the first rule only once is not allowed, because the
leftover `aab` could still react. Over 40 seeds, Chemart gave the first outcome
21 times and the second 19 times, and nothing else.

### Targets, dissolution and priorities

Each product can carry a *target*. `here`, the default, keeps it in the region;
`out` sends it to the surrounding region, or to the environment if the rule is
in the skin, where it is lost for good; `in` sends it into a randomly chosen
membrane directly inside; `in_j` into the inner membrane labelled `j`. A rule
with an `in` target cannot be applied in an elementary membrane. Chemart writes
the target after `@`: Păun's `aab → a b_out c c_in`, placed in the skin of
`[1 [2 ]2 ]1`, is the reaction `2 a@1 + b@1 -> a@1 + b@env + c@1 + c@2`.

The special symbol `δ` (delta) in a product dissolves the membrane of the region
where the rule fired, at the end of the step. The membrane disappears with its
rules, and everything inside it, objects and inner membranes alike, joins the
region around it. The skin can never be dissolved. Dissolution is the only way
the tree changes in this basic model.

A region can also rank its rules. A *priority* `r1 > r3` means `r1` goes first,
and Păun (2006, §12) gives two readings. In the *strong* one, `r3` is not used
at all in a step where `r1` can be applied. In the *weak* one, `r1` takes all
the objects it can and `r3` may use what is left.

### Halting and the result

All regions step together. A *computation* is a sequence of steps; it *halts*
when no rule can be applied anywhere, and only a halting computation has a
result. With an *output region* inside the system, the result is the number of
objects in it at halting (none, if that membrane was dissolved on the way).
With the environment as output, it is the number of objects sent out of the
skin. Because the choices are random, one system can halt with different
results. The set of all of them, written `N(Π)`, is what the system Π
*computes*.

### Worked example: is 7 a multiple of 3?

This is the book's system (eq. 9.21, "adapted from" Calude and Păun 2000 and
Păun 2002) and Chemart's default. The tree is `[1 [2 ]2 [3 ]3 ]1`. Region 2
starts with `n` copies of `a`, `k` of `c` and one `d`, and has three rules, with
both `r1` and `r2` above `r3`:

```
r1: a c  → c'
r2: a c' → c
r3: d    → d δ
```

The skin has one rule, `r4: d c c' → a_in3`, which puts an `a` into membrane 3,
the output region. Chemart's trace for `n = 7`, `k = 3`:

```
step 1   2/r1 ×3                 region 2: a 4, c' 3, d 1
step 2   2/r2 ×3                 region 2: a 1, c 3, d 1
step 3   2/r1 ×1                 region 2: c 2, c' 1, d 1
step 4   2/r3 ×1, dissolves 2    region 1: c 2, c' 1, d 1
step 5   1/r4 ×1                 region 1: c 1         region 3: a 1
```

Each step of `r1` or `r2` uses up `k` of the `a` and turns every `c` into `c'`
or back. If `n` is a multiple of `k`, the last full step leaves all `c` or all
`c'`; otherwise, as here, a mix. Only when the `a` run out can `r3` fire (strong
priority). Membrane 2 dissolves and spills `c`, `c'` and `d` into the skin,
where `r4` needs one `c` and one `c'` together, so it fires exactly when the mix
exists. Membrane 3 ends with one `a`: 7 is not a multiple of 3. With `n = 9`,
`r4` never fires and membrane 3 stays empty.

### From a P system to a reaction network

Chemart turns every rule into reactions over the `object@region` species: a rule
`u → v` in region `i` becomes `u@i → v@target`. The default network, printed
below, has four reactions, and two things about it are easy to misread. The
skin rule reads `d@1 + c@1 + c'@1 -> a@3`, though those species exist only after
membrane 2 dissolves: Chemart resolves each target in every tree reachable by
dissolution, so the network lists every reaction the system can ever perform
(status `complete`), and an `in` with several possible membranes gives one
reaction per way of sharing out the objects. And `r3: d → dδ` becomes
`d@2 -> d@2`, which changes nothing, because dissolution has no formula as a
reaction; it is recorded in `net.extras["dissolution"]`, and the priorities in
`net.extras["priorities"]`. The reactions have no rates, since P systems define
none.

## Using it

The default run is the book's divisibility system with `n = 7`, `k = 3` and
strong priorities. The computation, one random run of the maximally parallel
semantics, is in `net.extras["analysis"]`:

```python
a = net.extras["analysis"]
a["halted"], a["steps"], a["result"]     # (True, 5, 1)
a["output"]                              # {'a': 1}   contents of membrane 3
a["n_multiple_of_k"]                     # False
a["final_structure"]                     # '[1 [3 ]3 ]1'   membrane 2 is gone
a["applications"]                        # {'2/r1': 4, '2/r2': 3, '2/r3': 1, '1/r4': 1}
```

Rules are named `region/label`, so `2/r1` is rule `r1` of region 2.
`a["trace"]` holds the first 100 steps (rules applied, membranes dissolved,
configuration), and `a["dissolved"]` when each membrane went.
`net.extras["compartments"]` gives, per membrane, its parent, children, rules,
species and whether it can dissolve.

`n` and `k` set the input. `n=9, k=3` halts after 4 steps with membrane 3 empty
and `n_multiple_of_k` true. A run takes `n/k` steps rounded up, plus one when
`n` is a multiple of `k` and plus two otherwise; `n=1000000, k=7` takes 142,860
steps (about ten seconds) and needs `max_steps` raised above its default of
1,000.

**Why the priorities must be strong.** With `priority="weak"`, `r3` may fire as
soon as `r1` has taken all the `c` it can, which happens in the first step, so
membrane 2 dissolves early and the answer is wrong:

```python
a = chemart.generate_network("p-systems", n=7, k=3, priority="weak").extras["analysis"]
a["trace"][0]["applied"]      # {'2/r1': 3, '2/r3': 1}
a["n_multiple_of_k"]          # True   -- but 7 is not a multiple of 3
```

**The square numbers.** `system="n-squared"` is the system of Păun (2006),
figure 3, which sends a square number of `e` objects to the environment (see
*Results*). Most random runs never halt, so use many seeds and a small step
budget:

```python
from collections import Counter
results = Counter()
for s in range(2000):
    a = chemart.generate_network("p-systems", system="n-squared",
                                 max_steps=40, seed=s).extras["analysis"]
    results[a["result"]] += 1
# None (did not halt): 1370;  1: 469,  4: 130,  9: 20,  16: 10,  25: 1
```

Every halting run gave a square. The 2,000 runs take about five seconds.

**Your own system.** `system="custom"` takes `membranes` (the bracket string),
`objects` (a multiset per region, as `"3 f + c"` or `{"f": 3, "c": 1}`),
`rules` (a list per region, each optionally labelled `r1:`, with targets `@out`,
`@in`, `@in3` and `delta` or `δ` for dissolution), `priorities` (per region,
such as `"r1 > r2"`) and `output` (a membrane label, or `env` or `0` for the
environment). This is Păun's priority example, `r1: ff → f` above
`r2: cf → cdδ` on `fffc`, shown as (rules applied, membranes dissolved) per
step:

```python
spec = dict(membranes="[1 [2 ]2 ]1", objects={"2": "3 f + c"},
            rules={"2": ["r1: 2 f -> f", "r2: c + f -> c + d + delta"]},
            priorities={"2": ["r1 > r2"]}, output="1")
for pr in ("strong", "weak"):
    a = chemart.generate_network("p-systems", system="custom", priority=pr,
                                 seed=1, **spec).extras["analysis"]
    print(pr, [(s["applied"], s["dissolved"]) for s in a["trace"]])
# strong [({'2/r1': 1}, []), ({'2/r1': 1}, []), ({'2/r2': 1}, ['2'])]
# weak   [({'2/r1': 1, '2/r2': 1}, ['2'])]
```

Under the strong reading `r2` waits until one `f` is left; under the weak one
it takes the `fc` that `r1` left over and dissolves the membrane at once. A
delta rule in the skin, a priority cycle or an unknown `in_j` label is rejected
with an error. How Chemart picks one maximal multiset among the many allowed is
described in the implementation decisions below: it is one admissible choice,
not a published distribution.

## Results

P systems are a formalism, so the published results are mostly theorems about
what systems can compute and how fast, plus worked examples of the semantics.
Chemart reproduces the worked examples of the book and of Păun's 2006
introduction. It does not reproduce the theorems, which concern variants it
does not implement.

**Deciding divisibility (book §9.5).** In the book's system, `r1` and `r2`
alternately consume `k` copies of `a` per step, `r3` dissolves membrane 2 once
the `a` are gone, and `r4` puts one `a` into membrane 3 exactly when `n` is not
a multiple of `k`, including when `n < k`. The book itself calls this "way too
complex a solution for a simple problem": the strength of P systems lies in
parallel computation and in modelling biology. Chemart's tests check the four
reactions and the priorities, the step-by-step traces for `n = 7, k = 3` and
`n = 2, k = 5`, and the right answer for every `n` from 1 to 12 and `k` from 1
to 5. They also check that under weak priorities the system gives the wrong
answer. The book says `r3` fires "when neither r1 nor r2 can be fired", which is
the strong reading, so strong is Chemart's default.

**Generating the squares (Păun 2006, §§8–9, figure 3).** The system has three
nested membranes. In the innermost, `a → ab` and `f → ff` run together, so each
step adds one `b` and doubles the `f`, until `a → bδ` is chosen and dissolves
the membrane. After `n` growth steps this leaves `n + 1` copies of `b` and
`2^(n+1)` of `f`. In the middle membrane the `b` become `d`, each `d` makes one
`e` per step, and `ff → f` halves the `f`. The system can halt only if the
catalytic rule `cf → cdδ` waits until a single `f` is left; otherwise leftover
`f` reach the skin, where `f → f` fires for ever. A halting run therefore sends
`(n + 1)²` copies of `e` to the environment: `N(Π) = {n² | n ≥ 1}`. It is
Păun's illustration of how non-determinism and halting together define what a
system computes. Chemart's tests check the nine reactions and, over 150 seeds,
that every halting run gives a square after `2n + 2` steps for the result `n²`,
that 1 and 4 both occur, and that some runs never halt.

**The semantics (Păun 2006, §§4, 7 and 12).** The maximal-parallelism,
target and priority examples under *How it works* and *Using it* are Păun's.
Chemart's tests reproduce all three exactly, including the final contents
after `aab → a b_out c c_in` is used twice on `a^5 b^2 c^6` (`a^3 c^8` in the
skin, `c^2` in membrane 2, `b^2` in the environment). They also check Păun's
rule that when nested membranes dissolve in the same step, their contents land
in the first surviving membrane above them.

**Computing power and hard problems (not reproduced).** Most of the literature
is about computing power. Păun (2006) reports that all the classes of P systems
he surveys, with suitable features, are *universal*, equal in power to Turing
machines, often with very few membranes: one membrane with two catalysts
suffices (Freund, Kari, Oswald and Sosík 2005). On efficiency, variants that
build an exponential workspace in linear time solve NP-complete problems such
as SAT in polynomial time, trading space for time. The book names two ways to
do it, membrane division (P systems with *active membranes*, Păun 2001; Păun,
Suzuki, Tanaka and Yokomori 2004) and replicating strings ("worm objects",
Castellanos, Păun and Rodríguez-Patón 2000), and warns that physical resources
will always limit this growth. Păun (2006) adds that systems able to divide
non-elementary membranes solve PSPACE-complete problems in polynomial time, and
that division cannot be avoided: a P system without membrane division can be
simulated by a Turing machine with only polynomial slowdown. Chemart has no
membrane division, worm objects or other such features, so none of this is
reproduced.

**Modelling biology (not reproduced).** The book notes that the maximally
parallel, clock-driven semantics is "inconsistent with biological systems",
where reactions run at different speeds. Later work gave P systems stochastic
kinetics per compartment: the multicompartmental Gillespie algorithm of
Pérez-Jiménez and Romero-Campero (2006) and its variants (book §4.3), used for
example by Romero-Campero et al. (2009) and Smaldon et al. (2008). The book also
cites Suzuki and Tanaka's model of the p53 signalling pathway, written in ARMS,
"a special kind of P system" (§18.3.2). Chemart does not implement these
stochastic variants: its reactions have no rates, and the book gives neither
algorithmic details nor rate constants for them.

## Further reading

- Păun, G. & Rozenberg, G. (2002). A guide to membrane computing. *Theoretical
  Computer Science* 287, 73–100. doi:10.1016/S0304-3975(02)00136-6.
- Păun, G. (2002). *Membrane Computing: An Introduction*. Springer. According
  to the book, it contains a version of the divisibility system with explicit
  "yes" and "no" answers.
- Freund, R., Kari, L., Oswald, M. & Sosík, P. (2005). Computationally universal
  P systems without priorities: two catalysts are sufficient. *Theoretical
  Computer Science* 330(2), 251–266.
- Smaldon, J., Blakes, J., Krasnogor, N. & Lancet, D. (2008). A multi-scaled
  approach to artificial life simulation with P systems and dissipative
  particle dynamics. *Proceedings of GECCO 2008*, 249–256.
- Castellanos, J., Păun, G. & Rodríguez-Patón, A. (2000). Computing with
  membranes: P systems with worm-objects. *Proceedings of SPIRE 2000*, 65–74.
