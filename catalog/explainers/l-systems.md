## Introduction

An L-system, or Lindenmayer system, is a way of growing a string of symbols
by rewriting it over and over. You start from a short word, the *axiom*, and a
list of rules called *productions*, such as `a -> ab` ("replace every `a` by
`ab`") and `b -> a`. In one step *every* letter of the word is replaced at the
same time. From `a` you get `ab`, then `aba`, then `abaab`, then `abaababa`:
the words grow, and their lengths run 1, 2, 3, 5, 8, the Fibonacci numbers.

The biologist Aristid Lindenmayer introduced this kind of rewriting in 1968
as a mathematical theory of the development of simple multicellular
organisms; it was later applied to higher plants and plant organs
(Prusinkiewicz and Lindenmayer, 1990). Each letter stands for a cell, or a
larger plant part, in some state, and a production says what that part turns
into. The simultaneous replacement is the whole point. Formal grammars in the
tradition of Chomsky apply one rule at a time; L-systems apply them in
parallel because, in the words of *The Algorithmic Beauty of Plants* (ABOP),
"productions are intended to capture cell divisions in multicellular
organisms, where many divisions may occur at the same time". The difference
is not cosmetic: there are languages that context-free L-systems can generate
and context-free Chomsky grammars cannot (ABOP §1.1).

The original theory described only which parts neighbour which, not their
shape. Pictures came later, by reading the finished word as instructions for a
drawing "turtle", as in the LOGO language: `F` means draw a step forward, `+`
and `-` mean turn left or right by a fixed angle, and brackets `[` `]` save
and restore the turtle's position so that a side branch can be drawn and the
turtle can return to the stem. With this reading, one-line L-systems draw the
Koch snowflake, the dragon curve and branching plants. Banzhaf and Yamamoto
credit Lindenmayer and, later, Prusinkiewicz and Lindenmayer with pioneering
this use of turtle commands. Prusinkiewicz and Lindenmayer's book is the
source of every example in Chemart.

L-systems are a **formalism**, not a chemistry. There are no molecules
colliding, no rates, no energy and nothing is conserved: a word only grows or
changes. Banzhaf and Yamamoto give them one paragraph at the end of their
chapter on rewriting systems (book §9.8), as an example of rewriting used to
describe growth and development rather than computation. Their nearest
neighbours in the catalog share the rewriting but not the purpose.
[P systems](p-systems.md) also apply every applicable rule at once
("maximal parallelism"), but to multisets of objects in nested membranes, and
they are used to compute. [ARMS](arms.md) and [Gamma](gamma.md) rewrite an
unordered bag of symbols; an L-system word is ordered, and a letter's
neighbours can decide which rule applies to it. [MGS](mgs.md) generalises
rewriting to sequences, grids and graphs. What L-systems contribute to the
catalog is the parallel-rewriting model of development in its classic form,
with its best-known published examples.

## How it works

### Words, productions and derivation steps

An L-system has three parts: an alphabet of letters, an axiom (a nonempty
starting word) and a set of productions. A production `a -> χ` has a
*predecessor*, the single letter `a`, and a *successor*, the word `χ` (which
may be empty). One *derivation step* reads the current word left to right and
replaces every letter by the successor of its production, all at once. A
letter that has no production is replaced by itself. After `n` steps you have
the word of *derivation length* `n`. When every letter has exactly one
production and no rule looks at neighbours, the system is called a
**D0L-system**: deterministic (D) and context-free (0L, zero letters of
context). Its derivation is a single chain of words.

The first published example in ABOP models the filament of the blue-green
bacterium *Anabaena catenula*. The letters `a` and `b` are two cell states
(size and readiness to divide), and the subscripts `l` and `r` give the
cell's polarity, which decides on which side the daughter cells appear.
Chemart writes subscripts with an underscore:

```
a_r -> a_l b_r      a_l -> b_l a_r      b_r -> a_r      b_l -> a_l
```

Running it (`system="anabaena"`, four steps) gives the words ABOP prints:

```
a_r
a_l b_r
b_l a_r a_r
a_l a_l b_r a_l b_r
b_l a_r b_l a_r a_r b_l a_r a_r
```

In the second step both cells change at once: `a_l` divides into `b_l a_r`
while `b_r` simply matures into `a_r`. The filament has 1, 2, 3, 5, 8 cells.

### Branches and the turtle

Brackets turn a word into a branching structure. In the word `F[+F]F[-F]F`
the turtle draws a segment, then saves its position, turns left, draws a side
branch and returns; draws another segment, a branch to the right, and a last
segment. Brackets are never rewritten ("the brackets replace themselves"), so
a single production such as `F -> F[+F]F[-F]F` makes every segment of the
plant sprout the same pattern at the next step:

```
F
F[+F]F[-F]F
F[+F]F[-F]F[+F[+F]F[-F]F]F[+F]F[-F]F[-F[+F]F[-F]F]F[+F]F[-F]F
```

These are the first three words of ABOP figure 1.24a (`system="plant-a"`),
which is drawn after five steps with a turning angle of 25.7°. Chemart
generates the words but does not draw them.

### Context: letters that look at their neighbours

In a context-sensitive L-system a production can require particular
neighbours. `b < a -> b` means "an `a` whose left neighbour is `b` becomes
`b`"; `a > c -> ...` would look at the right neighbour; `x < a > y -> ...`
looks at both. ABOP calls these 1L-systems (one-sided context) and
2L-systems (two-sided), members of the wider class of IL-systems. ABOP
uses context to simulate interactions between plant parts, "due for example
to the flow of nutrients or hormones". A context-sensitive rule beats a
context-free rule for the same letter. The signal example of ABOP
§1.8 uses `b < a -> b` and `b -> a` from `baaaaaaaa`:

```
baaaaaaaa
abaaaaaaa
aabaaaaaa
aaabaaaaa
...
aaaaaaaab
aaaaaaaaa
```

At each step the `b` turns back into `a`, while the `a` just to its right
sees a `b` on its left and becomes `b`: the signal moves one cell to the
right per step. After nine steps it has left the string, and the word
`aaaaaaaaa` rewrites into itself for ever.

In a bracketed word, a letter's neighbour in the plant is not always its
neighbour in the text, so matching must skip over branches. A left-context
search walks down towards the root: it steps out of a branch at `[` and jumps
over complete side branches `[...]`. A right-context search walks up towards
the tips and skips side branches it does not ask for. ABOP's figure 1.29
example: the rule `BC < S > G[H]M -> X` applies to the `S` in
`ABC[DE][SG[HI[JK]L]MNO]`, skipping `[DE]` on the left and `I[JK]L` on the
right; Chemart turns the word into `ABC[DE][XG[HI[JK]L]MNO]`. A list of
*ignored* symbols (ABOP's `#ignore`, typically the turning symbols `+-`) is
skipped during matching, since turns are geometry, not plant parts.

### Stochastic L-systems

A plant drawn from a deterministic L-system is always the same plant. To get
different specimens of one species, a stochastic L-system gives a letter
several productions, each with a probability, summing to 1. Every occurrence
of the letter chooses independently. ABOP §1.7's example rewrites `F` into
`F[+F]F[-F]F`, `F[+F]F` or `F[-F]F` with probabilities 0.33, 0.33 and 0.34.
The derivation is now a tree of possible words rather than a chain.

### How Chemart turns this into a reaction network

Chemart offers two chemical readings of an L-system.

In the default reading, `reading="words"`, each **word is a species** and each
derivation step is a **reaction** turning one word into the next, with one
reactant and one product (a unimolecular reaction). The species
are named `w0`, `w1`, ... in the order they are found, starting with the
axiom `w0`, and the word itself is stored as the species' structure. A
deterministic system gives a chain `w0 -> w1 -> w2 -> ...`. A stochastic one
gives each word one reaction per distinct successor word, with the
probability of that successor (the sum over all the letter-by-letter choices
that produce it). A word that rewrites into itself, like the end of the
signal above, is not given a reaction. Probabilities are not rates, since L-systems define no
kinetics, so reactions carry no rate and the probabilities are kept on the
side. The reactor is "maximally parallel": all the parallelism lives inside
one reaction, which rewrites a whole word at once.

In the second reading, `reading="symbols"`, each **letter is a species** and
each production is a reaction. `a -> ab` becomes the reaction `a -> a + b`,
in which `a` is unchanged and so acts as a catalyst for making `b`. This
forgets the order of the letters, so it only works for context-free systems,
but it captures exactly what ABOP §1.9 uses to study growth: the matrix `Q`
whose entry `q_ij` counts the letters `a_j` in the successor of `a_i`, so
that the letter counts after `k + 1` steps are the counts after `k` steps
times `Q`.

## Using it

The default run above is ABOP's equation (1.2), `a -> ab`, `b -> a` from
axiom `a`, taken five steps: six words and five reactions. The words and the
run's bookkeeping are in the species and in `net.extras`:

```python
[s.structure for s in net.species]
# ['a', 'ab', 'aba', 'abaab', 'abaababa', 'abaababaabaab']
net.extras["analysis"]["lengths"]    # [1, 2, 3, 5, 8, 13]
net.extras["analysis"]["steps"]      # [0, 1, 2, 3, 4, 5]   derivation length of each word
net.extras["probabilities"]          # [1.0, 1.0, 1.0, 1.0, 1.0]   one per reaction
```

`extras["stationary"]` lists words that rewrite into themselves, with the
probability of doing so. `extras["turtle"]` gives the published turning angle
and derivation length of a named system, for anyone who wants to draw it. The
seed has no effect: nothing in the generator is random, because a stochastic
system is enumerated, not sampled.

**Published systems.** `system` picks one of ABOP chapter 1's examples: the
parameter table below lists them with their figure numbers. `iterations`
sets how many steps to take. For example:

```python
net = chemart.generate_network("l-systems", system="koch-island", iterations=3)
[s.structure.count("F") for s in net.species]   # [4, 32, 256, 2048]
net.extras["turtle"]                            # {'angle': 90.0, 'published_iterations': 3}
```

The quadratic Koch island of ABOP figure 1.6 starts as a square of four
segments, and each segment is replaced by eight.

**Your own system.** With `system="custom"`, give `axiom`, `productions` and
optionally `ignore` and `angle`. A production is written
`left < letter > right -> successor : probability`, where the contexts and
the probability are optional, `*` means "any context" and `ε` or nothing is
the empty successor. Symbols with subscripts are written `a_r`. Words with
subscripted symbols print with spaces between symbols. For figure 1.3, which
runs the algae rules from axiom `b`:

```python
net = chemart.generate_network("l-systems", system="custom", axiom="b",
                               productions=["a -> ab", "b -> a"], iterations=6)
[s.structure for s in net.species]
# ['b', 'a', 'ab', 'aba', 'abaab', 'abaababa', 'abaababaabaab']
```

**Stochastic systems.** The number of successors multiplies quickly. After
two steps of `stochastic-plant` there are already 301 words and 300
reactions, because `F[+F]F[-F]F` has five `F`s, each with three choices, and
all 3⁵ = 243 combinations give different words:

```python
net = chemart.generate_network("l-systems", system="stochastic-plant", iterations=2)
len(net.species), len(net.reactions), net.status   # (301, 300, 'complete')
net.extras["probabilities"][:3]   # [0.33, 0.33, 0.34]   F -> each of its three successors
```

A third step would need up to 3²⁵ combinations for a single word. Chemart
refuses to expand a word with more combinations than `max_species` and marks
the network `truncated`; the unexpanded words are listed in
`extras["analysis"]["unexpanded"]`. This is quick (under a second), but it
means only two steps of this system can be enumerated.

**Growth counts.** `reading="symbols"` gives the letter network and, in
`extras["analysis"]["counts"]`, the number of each letter after every step:

```python
net = chemart.generate_network("l-systems", reading="symbols", iterations=5)
[(r.reactants, r.products) for r in net.reactions]
# [({'a': 1}, {'a': 1, 'b': 1}), ({'b': 1}, {'a': 1})]      a -> a + b,  b -> a
net.extras["analysis"]["counts"]
# [{'a': 1}, {'a': 1, 'b': 1}, {'a': 2, 'b': 1}, {'a': 3, 'b': 2}, {'a': 5, 'b': 3}, {'a': 8, 'b': 5}]
```

Each column of the product matrix from `net.matrices()` is one row of ABOP's
`Q`. For a stochastic system the counts are expected values: for
`stochastic-plant` the expected number of `F`s grows by 3.66 per step
(0.33 × 5 + 0.67 × 3). Context-sensitive systems raise an error in this
reading.

**Budgets.** Words grow exponentially, so `max_length` (100,000 symbols by
default) stops the derivation before a word gets too long, and the network
is then marked `truncated`. With `iterations=30`, the default algae system
stops at its 24th word, 75,025 symbols long. `plant-a` with `iterations=7`
stops after derivation length 6 (39,061 symbols), which is still beyond the
5 steps of its ABOP figure. Raise `max_length` to go further. Every run
quoted on this page took well under a second.

**Direct use.** The module also exports the `LSystem` class, for rewriting
without building a network:

```python
from chemart.chemistries.l_systems import LSystem, word_text
d = LSystem("a", ["a -> ab", "b -> a"]).derivation(4)
[word_text(w) for w in d]   # ['a', 'ab', 'aba', 'abaab', 'abaababa']
```

## Results

L-systems are a formalism, so their "results" are the behaviours of
particular systems and the theory built around them. ABOP chapter 1 collects
the classic ones, and they are what Chemart reproduces. Nothing below comes
from a simulation in the chemical sense: every number is a count of symbols
in a derived word.

**Parallel rewriting and Fibonacci growth.** ABOP figure 1.3 derives `b`,
`a`, `ab`, `aba`, `abaab`, `abaababa`, `abaababaabaab` from the productions
`a -> ab`, `b -> a`, and §1.9 shows with the matrix `Q` that the number of
`a`s from axiom `a` follows the Fibonacci series 1, 1, 2, 3, 5, 8, ... .
The Anabaena filament (equation 1.1) shows the same idea as a model of real
cell division, with the five words printed in §1.2. Chemart's tests check
the figure 1.3 derivation word for word, the Fibonacci word lengths and
letter counts of the algae system up to ten steps, and the Anabaena words,
whose lengths also follow the Fibonacci numbers.

**Koch curves as L-systems.** ABOP §1.3 shows that Koch constructions, in
which every segment of a figure is replaced by a scaled copy of a
"generator", can be coded as L-systems: the initiator becomes the axiom and
the generator the successor of `F`. ABOP gives the quadratic Koch island
(figure 1.6), further islands, lakes and curves (figures 1.7 to 1.9), the
dragon curve and the Sierpiński gasket (figure 1.10). Chemart transcribes
all of them. Its tests check the segment counts: 4ⁿ segments for the Koch
curve after `n` steps (3 × 4ⁿ for the snowflake), 4 × 8ⁿ for the quadratic
Koch island, whose first derived word is also checked, 5ⁿ for the quadratic
snowflake, 2ⁿ for the dragon curve and 3ⁿ for the Sierpiński gasket. The
Koch snowflake itself is not printed as an L-system in ABOP; Chemart codes
the figure 1.1 construction as §1.3 prescribes. The figures themselves are
not reproduced, since Chemart has no turtle.

**Branching plants.** Bracketed L-systems of a single line or two draw the
plant-like structures of ABOP figure 1.24. Chemart has all six; its tests
check the first two words of figure 1.24a and that the letter reading gives
the same word lengths as the word reading for figure 1.24d.

**Growth functions.** A *growth function* gives the length of the word as a
function of the derivation length. ABOP §1.9 gives four kinds: exponential
(`F -> FF` doubles every step, 2ⁿ), Fibonacci (equation 1.2), polynomial
(`a_i -> a_i a_{i+1}`, whose letter counts form Pascal's triangle, so that
letter `a_i` grows as a polynomial of degree `i`), and a context-sensitive
2L-system (equation 1.5) whose length grows as `floor(√n) + 4`. The last
matters because the growth function of any D0L-system is, after a few steps,
a sum of polynomials times exponentials (ABOP equation 1.3, after Rozenberg
and Salomaa). Square-root growth is not of that form, and ABOP obtains it
with context. Chemart's tests check all four: 2ⁿ, the Fibonacci
counts, the Pascal triangle up to seven steps, and `floor(√n) + 4` for every
`n` up to 50. (ABOP's printed formula lost its floor brackets in text
extraction; the derivation gives the floor.) ABOP also quotes Vitányi's
argument that growth which keeps increasing but levels off, as sigmoidal
(S-shaped) growth does, cannot be obtained even with context, and introduces
the parametric L-systems of §1.10 to avoid the problem. Chemart does not implement parametric L-systems.

**Signals.** ABOP §1.8 uses context to send signals: the 1L-system that moves
a `b` along a string of `a`s, and signals in a branching structure that does
not grow (figure 1.30). An *acropetal* signal, carried by the left context,
travels from the root towards the tips and reaches every segment; a
*basipetal* signal, carried by the right context, travels from a tip towards
the root and reaches only the segments on the path. Chemart's tests check the
string signal word by word, including that it ends in a word that rewrites
into itself, and the four words of each branching signal, derived by hand
from the matching rules. They also check the figure 1.29 matching example,
and that three altered versions of its rule do not match.

**Hogeweg and Hesper's plants.** In 1974 Hogeweg and Hesper studied 3,584
patterns generated by a class of bracketed 2L-systems over the alphabet
`{0, 1}`, some of which had plant-like shapes. ABOP figure 1.31 gives five
such systems with turtle interpretation, drawn after 24 to 30 steps. Chemart
has all five; at 30 steps `hogeweg-hesper-a` reaches a word of 6,910
symbols. The tests check that the network agrees step by step with a direct
derivation and that the structure branches; they do not check the pictures.

**Different specimens of one species.** ABOP §1.7 shows that a stochastic
L-system draws plants that "look like different specimens of the same
(albeit fictitious) plant species" (figure 1.27). Chemart does not draw
specimens: it lists every possible word with its probability. Its tests
check the three first successors and their probabilities, that the
successors of every expanded word sum to 1, the 243 successors of
`F[+F]F[-F]F`, the budget cut-off, and the expected growth factor 3.66 per
step in the letter reading.

**What Chemart does not reproduce.** No turtle drawing is executed, so no
figure is reproduced, only the words behind it. Parametric L-systems (ABOP
§1.10), the three-dimensional turtle symbols, and the later chapters of ABOP
are not implemented. Lindenmayer's 1968 papers were not consulted; all
examples are checked against ABOP only. Banzhaf and Yamamoto also mention
recurrence systems (Herman, Lindenmayer and Rozenberg, 1975) and a model of
virtual cities built with L-systems (Kato et al., 1998) as context; neither
is implemented.

## Further reading

- Lindenmayer, A. (1968). Mathematical models for cellular interaction in
  development, Parts I and II. *Journal of Theoretical Biology* 18, 280–315.
  The original papers (book reference [512]; not consulted for Chemart).
- Hogeweg, P. & Hesper, B. (1974). A model study on biomorphological
  description. *Pattern Recognition* 6, 165–179. The 3,584 patterns behind
  ABOP figure 1.31.
- Szilard, A. L. & Quinton, R. E. (1979). An interpretation for DOL systems
  by computer graphics. *The Science Terrapin* 4, 8–13. As ABOP §1.3
  describes it, showed that very simple D0L-systems draw fractal curves.
- Rozenberg, G. & Salomaa, A. (1980). *The Mathematical Theory of
  L Systems*. Academic Press, New York. The formal-language theory,
  including the growth-function results cited in ABOP §1.9.
- Vitányi, P. M. B. (1986). Development, growth and time. In G. Rozenberg &
  A. Salomaa (eds.), *The Book of L*, 431–444. Springer. The argument against
  sigmoidal growth quoted in ABOP §1.9.
- The Algorithmic Botany site of Prusinkiewicz's group, with ABOP and later
  papers: <http://algorithmicbotany.org/papers/>
