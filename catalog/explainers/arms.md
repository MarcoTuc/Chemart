## Introduction

ARMS, the *Abstract Rewriting System on Multisets*, is a symbolic chemistry
proposed by Yasuhiro Suzuki and Hiroshi Tanaka of Tokyo Medical and Dental
University in 1997–1998. Molecules are bare symbols such as `a`, `b` or `c`,
the reaction vessel is a *multiset* (a bag in which the same symbol can occur
many times and order does not matter), and a reaction is a rewriting rule that
takes some symbols out of the bag and puts others in: `aaa → c` removes three
`a` and adds one `c`. Symbols can also be fed in from outside, and the bag can
be given a maximal size.

Suzuki and Tanaka start from the observation that real biochemical systems are
too complex to reconstruct precisely, so they wanted an abstract model that
keeps only the essentials and can still show the temporal behaviour of
chemistry: oscillations and cycles, which they point to as common in the
emergence of life. Their Artificial Life VI paper (1998) calls ARMS "like a
chemical solution in which floating molecules can interact with each other
according to reaction rules". It extends *abstract rewriting systems* from
mathematics and computer science, where a calculation is a chain of rewriting
steps that ends in a *normal form* (a state no rule can change), and it refers
to the Chemical Abstract Machine for its intuitive meaning. What the authors
say sets ARMS apart from Fontana's λ-calculus chemistry is the focus on
"temporal aspects and the emergence of cycles".

The book's small example makes the central point. Six rules and a starting bag
`{a, b, a, a}` are applied in rule-index order, and the bag settles in three
steps into a state where nothing can happen. Apply the *same* rules to the
*same* bag in the reverse order and the bag grows without end. How the rules
are scheduled is as much a part of the model as the rules themselves. Suzuki
and Tanaka then asked which property of a rule set decides whether its runs
stop, cycle or wander, and proposed a single number for it, the order
parameter λe (lambda-e), by analogy with Langton's λ parameter for cellular
automata.

ARMS is a simulation model: a minimal formalism plus an algorithm that runs
it. Banzhaf and Yamamoto present it in chapter 9, Rewriting Systems (§9.4),
between the [Chemical Abstract Machine](cham.md) (§9.3), from which it takes
the words *heating* and *cooling*, and [P systems](p-systems.md) (§9.5), of
which they call it "an example ... specialized for Artificial Life
experiments". [Gamma](gamma.md) (§9.2) also rewrites multisets, but as a
programming language whose final state is the answer whatever order the rules
fire in; ARMS is interested in exactly the cases where order and randomness
change the outcome. Later versions of ARMS added membranes and were used to
model an ecosystem (book §8.2.3) and the p53 signalling pathway (§18.3.2);
Chemart implements the membrane-free ARMS of 1998 only.

## How it works

### Symbols, rules and the bag

An ARMS is a pair Γ = (A, R): an alphabet A of symbols and a list R of
rewriting rules `r1, r2, …`. The state is the current multiset. A rule can be
applied when every symbol on its left side is in the bag (with enough copies);
applying it removes the left side and adds the right side. A rule with an empty
left side, `∅ → a`, is an *input*: it injects a symbol from outside. The rules
carry no rate constants: the book says the kinetic constants of reactions are
modelled by different frequencies of rule application.

Rules are classed by what they do to the size of the bag, borrowing the
Chemical Abstract Machine's terms:

- a **heating** rule makes the bag bigger (`a → abba`: one symbol in, four
  out). The paper describes it as breaking "a complex molecule into smaller
  ones";
- a **cooling** rule makes it smaller (`aaa → c`), rebuilding molecules from
  smaller ones;
- a rule that keeps the size (`b → d`) is neither. Chemart calls it
  *neutral*.

A run can stop. The paper's definition of a **normal form** is a multiset to
which no rule can be applied and to which no symbol can be input without
exceeding the maximal size. It is the ARMS version of a steady state, and a run
that reaches one is said to *halt* or *terminate*.

### Two ways of running the rules

The two sources schedule rules differently, and Chemart offers both.

**Ordered (book §9.4).** The rules are processed in a fixed order, by index or
by reverse index. The book does not say exactly what one step is; Chemart
reads it as one pass through the list in which each rule is applied as many
times as its left side fits into what is left of the current bag, and every
product joins the bag only at the next step. This is a deterministic,
maximally parallel update.

**Random (Suzuki and Tanaka 1998, figs. 2 and 14).** At each iteration any
input symbols are injected (if the maximal size allows), then *one* rule is
chosen at random and applied if its left side is present and the result does
not exceed the maximal size. Only successful rewrites count as steps. The rule
can be chosen uniformly, or with a bias: a heating rule with probability `p`
and a cooling rule with probability `1 − p`. That bias is the knob of the
paper's main experiment.

### Worked example: the book's six rules

The book's rule set, which is Chemart's default network:

```
r1: 3 a -> c          (cooling)
r2: b -> d            (neutral)
r3: c -> e            (neutral)
r4: d -> 2 f          (heating)
r5: a -> 2 a + 2 b    (heating)
r6: f -> h            (neutral)
```

Start from `{a, b, a, a}` and run the ordered discipline. Output of
`chemart.chemistries.arms.run(..., trace=True)`, one line per step, with the
size of the bag last:

```
order halting
  0 {'a': 3, 'b': 1} 4
  1 {'c': 1, 'd': 1} 2
  2 {'e': 1, 'f': 2} 3
  3 {'e': 1, 'h': 2} 3
```

In step 1, `r1` comes first and uses all three `a` to make `c`, and `r2` turns
`b` into `d`. By the time `r5` is reached no `a` is left, so the heating rule
never fires. In step 2, `c` becomes `e` and `d` becomes `ff`; in step 3 both
`f` become `h`. No rule has `e` or `h` on its left side, so `{e, h, h}` is a
normal form and the run halts, exactly as in the book's eq. 9.17.

Now the reverse order, `r6` first and `r1` last:

```
reverse-order non-recurrent
  0 {'a': 3, 'b': 1} 4
  1 {'a': 6, 'b': 6, 'd': 1} 13
  2 {'a': 12, 'b': 12, 'd': 6, 'f': 2} 32
  3 {'a': 24, 'b': 24, 'd': 12, 'f': 12, 'h': 2} 74
  4 {'a': 48, 'b': 48, 'd': 24, 'f': 24, 'h': 14} 158
  5 {'a': 96, 'b': 96, 'd': 48, 'f': 48, 'h': 38} 326
```

Now `r5` reaches the three `a` before `r1` does and turns each into `abba`, so
the bag holds six `a` and six new `b`; `r2` turns the original `b` into `d`.
Step 1 is the book's eq. 9.18, `{a,b,b,a,d,a,b,b,a,a,b,b,a}`. From then on `r5`
always gets the `a` first, the number of `a` doubles every step, and `r1`, the
only cooling rule, never fires.

### Measuring order: λe and cycles

Because a bag with a maximal size has only finitely many possible states, a
run that never halts must eventually revisit states: it cycles. Suzuki and
Tanaka take the *diversity of cycles* as their measure of order. An "ordered" system yields
simple cycles, such as a limit cycle; a "disordered" one yields chaotic or
complex cycles.

Their order parameter compares the two kinds of rules actually used in a run:

```
λe = Σ heating rules used / (1 + (Σ cooling rules used − 1))
```

which is simply the number of heating applications divided by the number of
cooling applications. λe is 0 when only cooling rules are used, 1 when the two
are used equally often, and "greater than 1.0" (the paper's words) when only
heating rules are used; Chemart then reports no value. The analogy is with
Langton's λ, the fraction of a cellular automaton's rule table that does not
lead to a quiet "quiescent" state: at low λ activity dies out, at high λ it is
chaotic, and complex behaviour sits in between.

The formal specification below restates all this: S is the alphabet, R the
rule list with heating, cooling and input rules, and A the two reactor
disciplines, with the maximal size as the only bound on the population.

## Using it

The default call above builds the book's six-rule network and runs it in
index order, reproducing eq. 9.17. The network is the rule list; the run is in
`net.extras`:

```python
a = net.extras["analysis"]
net.extras["rule_kinds"]   # ['cooling', 'neutral', 'neutral', 'heating', 'heating', 'neutral']
a["class"], a["steps"], a["final_state"]   # ('halting', 3, {'e': 1, 'h': 2})
a["rule_uses"]                             # [1, 1, 1, 1, 0, 2]
a["heating_uses"], a["cooling_uses"], a["lambda_e"]   # (1, 1, 1.0)
```

`analysis["class"]` sorts the run into one of four trajectory classes:
`halting` (a normal form was reached), `periodic` (a deterministic ordered run
returned to an earlier bag; `period` holds the number of steps between the
visits and the run stops there), `recurrent` (a random run revisited bags
without halting) and `non-recurrent` (neither, within `steps`). `revisits`
counts returns to an earlier bag, `periods` lists the distinct return times,
and `kinds_of_periods` is their number, Chemart's version of the paper's
"kinds of periods". `rules`, `inputs` and `max_size` echo the system; the
input reactions `∅ → s` appear in the network labelled `input` in
`reaction_rules`.

**The reverse order.** `selection="reverse-order"` gives the expanding run.
Ten steps take the bag from 4 to 10,742 symbols; heating was used 4,600 times
and cooling never, so λe is `None`:

```python
a = chemart.generate_network("arms", selection="reverse-order", steps=10).extras["analysis"]
a["class"], a["final_size"]      # ('non-recurrent', 10742)
```

The bag roughly doubles at every step; Chemart stores only counts, so this
stays fast (the tests run 30 steps and pass 2³⁰ symbols).

**The paper's λe experiment.** `system="two-symbol"` builds the rule set of
Suzuki and Tanaka's simulation: every rule whose two sides are bags of one to
five symbols drawn from `{a, b}`, 380 rules in all (see the decisions below on
how that count was chosen), with maximal size 10 and a random starting bag of
1–10 symbols. Selection is the heating-probability bias, and a run lasts 1000
rewriting steps, as in the paper:

```python
net = chemart.generate_network("arms", system="two-symbol", seed=1)   # p = 0.5
a = net.extras["analysis"]
# class 'recurrent', steps 1000, heating_uses 508, cooling_uses 492,
# lambda_e 1.03, revisits 936, kinds_of_periods 188
```

With `p=0.0` the same seed halts after 3 steps at `{a}`, λe = 0; with `p=1.0`
it halts after 3 steps at a full bag `{a: 4, b: 6}`. One run takes a fraction
of a second; a sweep over 11 values of `p` with 20 seeds each takes about 15 s.
`rule_count` samples a random subset of the 380 rules instead of all of them.

**Ru1.** `system="ru1"` is the paper's four-rule example, `aaa → b`,
`b → a`, `b → c`, `a → bb`, with `a` input at every iteration and maximal size
4, run with uniform random selection. Every one of 20 seeds (0–19) halts, because
`c` is consumed by no rule and piles up: 16 end at `{a, a, c, c}`, 3 at
`{a, c, c, c}` and 1 at `{c, c, c, c}`. The paper's fixed rule sequence is
available through `chemart.chemistries.arms.run(..., sequence=[3, 0, 2, 1])`,
not as a parameter.

**The Brusselator.** `system="brusselator"` writes the Brusselator oscillator
(see [the Brusselator](brusselator.md)) as four rewriting rules, `A → X`,
`B + X → Y + D`, `2X + Y → 3X`, `X → E`, with `A` and `B` input continually
into a bag of at most 5000 symbols. All four rules keep the size, so the
inputs fill the bag; seed 1 with `steps=10000` halts at step 9288 with a full
bag `{B: 212, D: 2288, E: 2356, Y: 144}`, having run out of `A` and `X`.
Chemart does not produce the paper's oscillations (see Results).

**Your own rules.** `system="custom"` takes `rules`, `initial` and `inputs`.
Rules use the reaction syntax shown above; an empty side is allowed:

```python
net = chemart.generate_network("arms", system="custom",
                               rules=["a -> b", "b -> a + c", "c -> "],
                               initial={"a": 1}, selection="order")
# rule_kinds ['neutral', 'heating', 'cooling']; class 'periodic', period 2
```

The bag goes `{a} → {b} → {a, c} → {b}`: back to `{b}` after two steps.
`max_size` sets or removes (`-1`) the bound and `input_probability` makes
inputs occasional rather than certain; the table below lists every parameter.

## Results

**Rule order decides halting (book §9.4).** The book, following Suzuki and
Tanaka's 1997 paper, uses the six-rule example to show that the same rules and
the same starting bag halt in one rule order (eq. 9.17) and expand without
bound in the reverse order (eq. 9.18). Chemart reproduces both exactly; the
tests check the full halting sequence and the first expanding step, and that
the size then grows at every step.

**A step-by-step walk-through (1998, fig. 3).** The paper illustrates its
random algorithm with Ru1, maximal size 4, input `a` and the rule order `r4,
r1, r3, r2`: `r4` is refused because `a → bb` would make the bag 5 symbols,
`r1` gives `{b, a}`, then an `a` is input and `r3` gives `{c, a, a}`, then an
`a` is input and `r2` cannot apply, giving `{c, a, a, a}`. Chemart's tests
reproduce this sequence. The paper prints the initial bag as `{a, a, f, a}`,
with a symbol `f` no rule uses; the successors only follow from `{a, a, a, a}`,
which Chemart uses. The paper's figs. 4 and 5, where two other fixed orders of
the same rules give two cycles of period 3 or a halt, could not be reproduced
under any of 16 readings of the algorithm and are not claimed.

**Cycles from simple setups.** In a "simple setup" (five symbols, input `a`,
maximal size 10, six rules, varying the frequency of inputs and the randomness
of rule choice) the paper reports that cycles emerged even under simple
conditions, which it contrasts with Kauffman's networks and Fontana's
λ-calculus chemistry, both of which "need large-scale computation to generate
cyclic structures". Fusion of cycles and period doubling appeared easily once
randomness was added to the inputs (fig. 6). The paper does not give these six
rules, so Chemart does not reproduce this.

**Oscillations from the Brusselator.** To show that ARMS "works as an abstract
chemical system", the paper writes the Brusselator as rewriting rules, with
reaction rates represented by how often each rule is applied, an empty
starting bag, continual inputs of `A` and `B` and maximal size 5000. It
reports oscillations in the numbers of `X` and `Y`, of three types: quasi-stable
oscillation, unstable oscillation, and divergence followed by convergence
(figs. 9–11). The rule frequencies are not published. Chemart has the rules
(the tests check them), but with uniform selection and products `D` and `E`
that are never removed, its bag simply fills up. The oscillations are not
reproduced.

**The order parameter λe (1998, figs. 15–17).** The main experiment used the
two-symbol rule set with no inputs and maximal size 10. The heating
probability `p` was raised from 0 to 1 in steps of 0.01; at each value, 100
random starting bags were run for 1000 steps. The paper reports that below
`p = 0.1` most calculations terminated, that between 0.3 and 0.85 only a few
did, and that above 0.85 the number rose again (fig. 15). The number of
generated cycles rose with `p`, stayed around 450 between 0.3 and 0.8, and fell
again (fig. 16). The number of *kinds* of periods was highest near `p = 0.5`,
where heating and cooling rules are used equally often (fig. 17). The authors
summarise this as a table linking termination to λe, to Langton's λ and to
Wolfram's four classes of cellular-automaton behaviour: strong termination at
λe near 0 or much larger than 1 (Wolfram class I or II), weak termination and
the most complex cycling at λe near 1 (class III). They read the "edge of
chaos" as arising from a slightly biased ratio of heating to cooling rules.
The book summarises the finding as: simple dynamics for small and large λe,
cycles for intermediate values.

Chemart reproduces part of this. Its tests check that, over five seeds, every
run halts at `p = 0` (λe = 0, a bag of one symbol) and at `p = 1` (no cooling,
a full bag), that every run at `p = 0.5` is still going after 1000 steps with
λe between 0.8 and 1.25, and (a slow test over 20 seeds) that the average
kinds of periods is higher at `p = 0.5` than at 0.1 or 0.9. A sweep over 20
seeds shows the peak clearly:

```
p     halted  kinds  lambda_e
0.0    20/20    0.0  0.00
0.1     0/20  153.3  0.65
0.3     0/20  190.6  0.85
0.5     0/20  192.6  1.00
0.7     0/20  184.1  1.16
0.9     0/20  166.4  1.35
1.0    20/20    0.0  none
```

What it does not reproduce is the gradual termination curve of fig. 15. With
all 380 rules some rule can always be applied at any `p` strictly between 0 and
1, so no run halts there; runs on random subsets of 12 or 30 rules halt in
roughly the same fraction at every intermediate `p` rather than mostly near
the ends. The paper's own count of rules (30976, which its formula does not
give for two symbols) leaves its rule set uncertain; the decisions below
explain the reading Chemart uses. The absolute numbers of cycles and kinds of
periods also depend on how they are counted, which the paper does not define
precisely, so Chemart's numbers are not comparable with its figures.

**Later work.** The book describes the later ARMS as a variant of P systems,
whose membranes grow and divide to simulate artificial cells, and reports two
applications. Suzuki and colleagues used it to model a plant, herbivore and carnivore
ecosystem, finding the ecology more robust when plants produce a
herbivore-induced volatile substance that attracts carnivores (book §8.2.3,
[828], [829]). Suzuki and Tanaka (2006) expressed the core of the p53
signalling pathways in ARMS and obtained regulatory feedback loops that agree
with the biological literature (book §18.3.2, [832]). These papers are closed
access and the book gives no rules, so Chemart implements neither, nor the
membranes.

## Further reading

- Suzuki, Y. & Tanaka, H. (1997). Symbolic chemical system based on abstract
  rewriting and its behavior pattern. *Artificial Life and Robotics* 1,
  211–219 (book [830]). The source of the book's six-rule example; not
  accessible for this entry.
- Suzuki, Y., Fujiwara, Y., Takabayashi, J. & Tanaka, H. (2001). Artificial
  life applications of a class of P systems: abstract rewriting systems on
  multisets. In *Multiset Processing*, LNCS 2235, 299–346, Springer (book
  [828]). The membrane version of ARMS.
- Langton, C. G. (1991). Life at the edge of chaos. In *Artificial Life II*,
  Addison-Wesley. The λ parameter that λe is modelled on.
