## Introduction

BFF is a soup of tiny computer programs that run on each other. Blaise Agüera
y Arcas and seven colleagues at Google and the University of Chicago described
it in 2024, in a paper asking how self-replicators arise in the first place,
that is, how a world goes from "pre-life" dynamics to life. Systems such as
Tierra and Avida study evolution *after* that point: they are seeded with a
hand-written program that copies itself, the "ancestor". BFF starts from
random bytes, with no ancestor and no fitness function.

Each molecule is a program of 64 bytes, written in a variant of Brainfuck, an
esoteric programming language famous for having only eight instructions. The
authors change it so that code and data sit on the same tape and the program
can read and write itself. The acronym BFF names this family of extended
Brainfuck languages. A reaction takes two programs, glues them into one
128-byte tape and runs it. Either program can overwrite itself or the other.
Afterwards the tape is cut back into two programs.

Nothing in the system rewards copying. Yet in a sizeable fraction of runs a
program appears that copies itself into whatever it is paired with. It then
takes over the soup. The paper calls this a *state transition* and uses it to
separate pre-life from life. Its main claim is that the first replicators come
mostly from the programs' interactions and self-modification, not from lucky
random initialisation or from mutation.

The authors describe the soup as a chemistry of its own: a variant of Walter
Fontana's "Turing gas", the setting of [AlChemy](alchemy.md), with the
reaction `A + B -> A' + B'`. That is why it is in the catalog, next to
[automata reaction](automata-reaction.md), where bit-string automata act on
each other. Tierra, Avida and Core War, which start from hand-written
programs, are kept in the archive. The paper appeared in 2024, after
Banzhaf and Yamamoto's book, so the book does not cover it.

## How it works

### The language

The tape holds 128 bytes: program A followed by program B. Three positions
move over it: the instruction pointer, which starts at the first byte, and two
data heads, `head0` and `head1`, which also start at 0. Ten byte values are
instructions (paper section 2):

| byte | effect |
|---|---|
| `<` `>` | move `head0` left or right |
| `{` `}` | move `head1` left or right |
| `-` `+` | decrement or increment the byte under `head0` |
| `.` | copy the byte under `head0` to `head1` |
| `,` | copy the byte under `head1` to `head0` |
| `[` | if the byte under `head0` is zero, jump forward to the matching `]` |
| `]` | if the byte under `head0` is not zero, jump back to the matching `[` |

Every other byte does nothing and can hold data. The heads wrap around the
tape, so moving `head1` left from the first byte puts it on the last. A
program ends when the instruction pointer runs off the end, when a jump finds
no matching bracket, or after 8192 bytes have been read.

### The soup

The reactor is the paper's "primordial soup". Each *epoch*, every byte of
every program is first replaced by a random byte with a small probability
(0.024% by default), the background mutation. The programs are then shuffled
into random ordered pairs, and each pair runs once as described above. No
program is added or removed. Only running and mutation change the soup.

The paper also puts the programs on a grid (section 2.2). There each program
can only pair with one no more than two cells away along each axis.

### A worked example: a self-replicator

The paper's Figure 4 shows a replicator, cleaned up by the authors: eight
bytes of code at each end and filler in between. It is a palindrome, reading
the same backwards:

```python
from chemart.chemistries.bff import FIG4_REPLICATOR, encode, run, show

tape = bytearray(encode(FIG4_REPLICATOR) + bytes(64))   # the replicator, then 64 zeros
print(run(tape), show(tape[64:]) == show(tape[:64]))
print(repr(show(tape[:64])))
```

```
8192 True
'[[{.>]-]ĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠĠ]-]>.{[['
```

The run used all 8192 steps, and the second half now equals the first. Here is
how. The outer loop `[ ... ]` never ends. Inside it, `{` moves `head1` one byte
to the left, which from byte 0 wraps it to the end of the tape. `.` copies the
byte under `head0`, and `>` moves `head0` one byte to the right. So the read
head walks forward through A while the write head walks backward from the end
of B, and A lands in B reversed. Because A is a palindrome, the reversed copy
is exact. Chemart's interpreter matches the paper's step-by-step trace: the
first byte lands at the end of B after 4 steps, and B is complete after 256.

It does not matter what B held before. That is the paper's equation (5): a
replicator S turns any food program F into a copy of itself, `S + F -> 2 S`.
In the network this shows up as a reaction with S on both sides, the pattern
Chemart records as a catalyst.

The filler byte is a space, shown as `Ġ` (see below). With zero bytes as
filler the copy is no longer exact: the inner loop `[{.>]` stops on each zero,
the `-` after it turns the zero into 255, and both halves end with 255 where the
zeros were.

### Complexity: high-order entropy

To detect a state transition without knowing what a replicator looks like,
the paper defines *high-order entropy*. It is the Shannon entropy of the
soup's bytes, taken one at a time, minus an estimate of how far the soup can be
compressed, in bits per byte. The compression estimate stands in for the
Kolmogorov complexity. Random bytes score about 0: they can be compressed no
further than their byte frequencies allow. A soup full of copies of one
program scores high. Like the authors' code, Chemart estimates the
compression with the Brotli compressor at quality 2.

## Using it

The default call runs 128 random programs for 16 epochs, the paper's
settings except for scale: the paper uses 2^17 = 131,072 programs for up to
16,000 epochs. That is enough to see the soup react, but a replicator is not
expected to appear. Each species is a 64-byte program, printed one character
per byte: the instructions as themselves, the zero byte as `0`, and any other
byte b as the character U+0100 + b, following the authors' code. That is why
random programs look like strings of accented letters.
`extras["reaction_kinds"]` says whether each reaction is an `execution`
(`A + B -> A' + B'`) or a `mutation` (`A -> A'`).

```python
from collections import Counter
net = chemart.generate_network("bff", seed=1)
a = net.extras["analysis"]
print({k: v[-1] for k, v in a.items()})
print(Counter(net.extras["reaction_kinds"]))
```

```
{'epoch': 16, 'high_order_entropy': -0.016, 'distinct_tapes': 128, 'top_tape_count': 1, 'ops_per_run': 29.67, 'zero_bytes': 30}
Counter({'execution': 609, 'mutation': 37})
```

`extras["analysis"]` records the soup every `record_every` epochs:
`high_order_entropy`, `distinct_tapes`, `top_tape_count` (copies of the most
common program), `ops_per_run` (instructions executed per pair, no-ops not
counted) and `zero_bytes`. Replicators loop until the step limit, so
`ops_per_run` climbs towards 8192 when they spread. `extras["final_state"]`
is the final soup.

### A seeded run

The paper's "seeded" runs place one copy of the Figure 4 replicator in a
random soup (`replicators=1`). With 64 programs and 100 epochs (about five
seconds):

```python
net = chemart.generate_network("bff", seed=8, tapes=64, epochs=100,
                               replicators=1, record_every=10)
a = net.extras["analysis"]
print(a["top_tape_count"]); print(a["high_order_entropy"]); print(a["ops_per_run"])
```

```
[1, 17, 30, 34, 42, 46, 60, 56, 34, 42, 31]
[-0.0591, 1.5553, 2.5642, 2.233, 1.9313, 1.3853, 1.3162, 1.2657, 1.1285, 1.1707, 1.0997]
[0.0, 2197.94, 3907.16, 4645.09, 6057.75, 7783.66, 8192.0, 8192.0, 8192.0, 8192.0, 8192.0]
```

By epoch 60 the replicator fills 60 of the 64 slots and every pair runs to the
step limit. In this run the replicator reacted 682 times as `S + F -> 2 S`.
Many seeds do not end like this: a replicator paired second is usually
overwritten by its partner before it runs, so a single copy often dies out.

`space="grid"` with a `width` runs the 2D soup. Pure Python runs a pair of
replicators in about 2 ms, so a soup of 1,024 replicators takes about a second
per epoch. The paper's runs belong to the authors' CUDA code.

## Results

### Self-replicators emerge without being put in

The paper's central result concerns random soups of 2^17 programs run for 16,000
epochs with 0.024% background mutation. About 40% of runs show a state
transition, a jump in high-order entropy that marks a replicator taking over
(Figure 5). High-order entropy first rises during the first 1,000 epochs and
then falls again on average. The transition can come at any time, and a few
runs have it almost immediately.

Varying the mutation rate from 0 to 1% (Figure 6) shows that more mutation
speeds replicators up, and that even with no mutation at all transitions
happen about as often as with the default. So mutation alone does not explain
them. Figure 7 compares four kinds of run:

- "short": random programs, 128 epochs. A transition happens only 3 times in
  1,000, so replicators present by chance at the start are very rare.
- "seeded": one Figure 4 replicator added, 128 epochs. It takes over 22% of the
  time. Many single replicators are destroyed before they spread, for example
  because half of the time they are the second program of a pair.
- "long": random programs, 16,000 epochs, about 40%.
- "long-no-noise": no mutation and a fixed pairing pattern, even more often,
  about 50%.

From these the authors conclude that self-replicators arise mostly through
self-modification and interaction, not through random initialisation or
mutation.

Chemart reproduces the contrast between seeded and short runs at small scale.
The slow test
`test_seeded_runs_reach_a_state_transition_and_random_short_runs_do_not` runs
12 seeds of 64 programs for 100 epochs. Seeded soups reach a high-order
entropy of at least 1, the threshold the paper uses in Figure 6, in 2 to 10
of the 12 runs. Random soups never do. The emergence of replicators from
random programs over 16,000 epochs is out of reach in pure Python and is not
tested.

### How the first replicator arose

In one run the authors tag every byte with a tracer token recording where it
came from, and follow the tokens through the transition (Figures 1-3). The
number of distinct tokens drops sharply as a few take over the soup, at the
same moment as the jump in high-order entropy. The first replicator appeared
at epoch 2354 (the captions of Figures 2 and 3 give 2355 and 2354). A loop
that was incomplete without its partner copied every second byte of one
program into the other, in reverse, and so built a reversed replicator there.
Copying went on until a complete replicator sat in the second program. A
"zero-poisoning" period followed before a new family of replicators took over.

### Replicators in space

On a 240 × 135 grid, self-replicators still arise (Figure 8). Once one appears
it spreads as a wave, taking a number of epochs proportional to the grid's
side. In the well-mixed soup, a replicator typically takes over at least half
of n programs in about log n epochs. The grid
leaves room for several replicator variants to coexist and compete. Chemart's
grid follows the authors' neighbourhood exactly, which the test
`test_grid_pairs_only_neighbours` checks. The wave itself is not measured.

### What Chemart checks

- The language: the Figure 4 trace, row by row, and the exact copy
  (`test_fig4_trace_end_of_tape_b`, `test_fig4_trace_first_byte_of_tape_b`,
  `test_fig4_replicator_copies_itself_exactly`).
- Equation (5): the replicator turns 20 random food programs into copies of
  itself (`test_replicator_turns_any_food_into_itself`), and the reaction
  appears in a soup's network with the replicator on both sides
  (`test_seeded_replicator_shows_as_catalyst`).
- High-order entropy is about 0 for random bytes and above 5 bits per byte for
  a soup of copies.

### Beyond BFF

The paper also finds replicators in variants of Forth, a stack-based language,
both in soups and on one long shared tape, and in emulated Z80 and 8080
processors. It gives one counterexample. In SUBLEQ, a language with a single
instruction, hand-written replicators work when seeded, but none arose from
random soups. The authors suspect this is because the shortest possible
replicator is much longer there (60 bytes for their smallest hand-written one).
Chemart implements only BFF.

## Further reading

- Fontana, W. (1990). Algorithmic chemistry: A model for functional
  self-organization. The Turing gas that the primordial soup varies.
- The authors' code, with the BFF, Forth and SUBLEQ variants:
  https://github.com/paradigms-of-intelligence/cubff
