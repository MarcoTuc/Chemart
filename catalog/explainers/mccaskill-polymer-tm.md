## Introduction

John McCaskill's *pattern processing chemistry* is a soup of binary strings
in which every string is two things at once: a piece of tape that can be read,
and a tiny machine that can read other strings and write new ones. Two strings
react when they recognise each other. One then acts as the machine, reads the
other from end to end, and writes one or more brand-new strings, which join
the soup. Neither of the two reacting strings is used up. McCaskill set it out
in a 1988 internal report of the Max Planck Institute for Biophysical
Chemistry in Göttingen, *Polymer chemistry on tape: a computational model for
emergent genetics*, and ran it as a C program.

The report's question is how a genetic system could organise itself. Its
starting point is the quasispecies model (see
[the quasispecies chemistry](quasispecies.md)), in which replication is a
built-in event and the link between a sequence and how fast it replicates has
to be written in by hand. McCaskill wanted a level below that. In his words,
the level chosen "is definitely chemical rather than biological : for example,
replication is not an elementary event. The biology should not be defined."
Here, copying is not a rule of the model. It happens only if some string
happens to encode a machine that copies. The report's abstract puts the aim as
investigating "the self-organization of reaction pathways showing sequence
heredity, an emergent genetics", and its introduction names the questions it
hoped to reach: early biological organisation, "its limitation by errors and
parasites and its evolution".

The report's first experiments are about exactly those parasites. It gives a
19-symbol string that replicates, and a second string that the replicator's
machine will copy but that encodes no machine of its own, so it is copied for
free. It then puts both into a random soup and follows who survives.

It is a simulation model. Banzhaf and Yamamoto describe it in book §10.5.3,
"Polymers as Turing Machines", among the chemistries built from machines and
tapes. Its nearest neighbours in the catalog differ in who plays which role.
In [Laing's molecular machines](laing-molecular-machines.md), the first of
the family, strings are either passive tapes or active machines, and there is
no simulated dynamics. In [typogenetics](typogenetics.md), a strand codes for
enzymes that act back on that same strand. In
[Ikegami and Hashimoto's machine–tape chemistry](ikegami-hashimoto.md),
machines and tapes are two separate populations. In
[matrix chemistry](matrix-chemistry.md), as here, any string can act on any
other and both survive, but the action is a matrix product rather than a
program. What McCaskill adds is *specific recognition*: which strings meet is
itself written in the strings, as patterns with "don't care" positions. The
chemistry was later taken into McCaskill's special-purpose parallel hardware,
which the book and the review of Dittrich, Ziegler and Banzhaf (2001) mention
but Chemart does not model.

## How it works

### Two readings of one string

A molecule is a string of 0s and 1s, for example the report's replicator
`0001111100111010111`. The string is never read one symbol at a time.
It is read in *doublets*, pairs of neighbouring symbols, and each doublet
means something. The triplet `111` is special: it is the *initiator*, the
marker where a meaningful stretch begins. From every `111` in the string, two
things are read.

**A recognizon**, the string's recognition site. Starting at the initiator,
each doublet becomes one symbol of a pattern: `01` means 0, `10` means 1, and
`00` or `11` mean `#`, "don't care". Reading stops at a later `111`, at the
end of the string, or after `R` symbols (16 by default, the report's typical
value). The more `#`s a recognizon has, the less choosy it is. The report takes this
kind of pattern over {0, 1, #}, a "condition", from John Holland.

**A rule** of a machine. The twelve bits that start at the initiator are six
doublets, and each doublet sets one field of a rule:

| field | `01` | `10` | `00` | `11` |
|---|---|---|---|---|
| READ: tape symbol it applies to | 0 | 1 | blank | 0 or 1 |
| WRITE: symbol written | 0 | 1 | blank | opposite of the symbol read |
| STATE: machine state it applies to | 0 | 1 | 0 or 1 | 0 or 1 |
| NEXT STATE | 0 | 1 | unchanged | flipped |
| READ MOTION: read head moves | stay | right | left | right |
| WRITE MOTION: write head moves | stay | right | left | right |

The table is the report's own (section 3, item 5). A string's rules together
make a small *processor*, a machine with two internal states and two heads.
The read head walks over the other molecule, the *tape*, which is not
changed. The write head walks over a fresh, blank tape. At each step the
processor looks up the rule for the symbol under the read head and its
current state, writes, moves both heads and changes state. When no rule
applies, it halts, and each unbroken run of symbols on the written tape is
released into the soup as a new molecule. A processor that has not halted
after `max_steps` steps releases nothing (the report cuts such loops off too).

### Worked example: how the replicator copies

The replicator `0001111100111010111` has five 1s in a row at positions 3 to
7 (counting from 0), so it has three overlapping initiators, at positions 3,
4 and 5. Their twelve-bit rules are:

```
position 3: 11 11 10 01 11 01   read 0 or 1, write opposite, state 1,      next 0,    read right, write stay
position 4: 11 11 00 11 10 10   read 0 or 1, write opposite, state 0 or 1, next flip, read right, write right
position 5: 11 10 01 11 01 01   read 0 or 1, write 1,        state 0,      next flip, read stay,  write stay
```

When two rules cover the same symbol and state, the later one wins. Rule 4
replaces rule 3 in state 1, and rule 5 replaces rule 4 in state 0. What is
left is a two-step dance for every symbol of the tape:

1. In state 0, write a 1, keep both heads still, go to state 1.
2. In state 1, read the same symbol again, overwrite the 1 with the opposite
   of that symbol, move both heads right, go back to state 0.

When the read head runs off the end of the tape it sees a blank, no rule
covers a blank, and the machine halts. The written tape now holds the
*complement* of the tape, every 0 turned into 1 and every 1 into 0. So when
the replicator reads itself, it writes its complement, and when it reads its
complement, it writes itself back. From the default network:

```
2 p0001111100111010111 -> 2 p0001111100111010111 + p1110000011000101000
```

This is the replicator (processor) reading a second copy of itself (tape) and
releasing `1110000011000101000`, its complement. Species are named `p`
followed by the string. Like a DNA or RNA strand copied through its opposite
strand, the replicator reproduces in two rounds rather than one; the
implementation decisions below explain why Chemart reads the report's
"self-replicating" this way.

The parasite `0000000100100010111` has a single `111`, at its very end, so it
encodes no rule at all. The replicator copies it just the same, through its
complement:

```
p0001111100111010111 + p0000000100100010111 -> p0001111100111010111 + p0000000100100010111 + p1111111011011101000
```

### Recognition: who meets whom

Before any processing, two molecules must *collide*, and that depends on
their recognizons. In the report, each molecule drawn from the soup places
one pattern in a shared *pattern space*, with every `#` replaced by a random
0 or 1. If another molecule has already placed the same pattern there, the
two collide. The molecule that placed its pattern first becomes the
processor, the newcomer becomes the tape, and both patterns are cleared.
Chemart's closure uses the equivalent test: two molecules can collide when
they have recognizons of the same length that agree wherever neither has a
`#`.

For the published strings, the recognizons are:

```
replicator  0001111100111010111   ##10#00#   ###   #10#00#   #11   #
parasite    0000000100100010111   #
```

The last `111` of the replicator and the only `111` of the parasite both sit
in the shared tail `...0111`, which gives each a one-symbol recognizon `#`.
That is what the report requires of a parasite: "the same recognizon" as the
replicator, but no replication rule.

### What the reactor does

A reaction is always `s1 + s2 -> s1 + s2 + s3 (+ ...)`: processor and tape
both survive, and the products are new strings. In the report's soup, the
strings sit in one long array of about 2^20 symbol positions, and new
strings overwrite old ones. Chemart's soup keeps a population of whole
molecules of fixed size instead: each new string replaces a molecule chosen
at random. There are no rates. The only "kinetics" is how often molecules
meet, which comes from how common they are and how specific their
recognizons are, as the report intends. An optional error rate flips written
symbols, which acts as mutation.

## Using it

The default run is not a simulation. It is the *closure* of the two
published strings: every reaction that can happen, starting from the
replicator and the parasite, and then among everything they make, until
nothing new appears. It has seven species: the replicator, the parasite,
their two complements, and three strings that the complements write:
`1111111111111111111`, `0000000000000000000` and the one-symbol `1`. The
all-ones string is worth a look: every twelve bits of it read as "read 0 or
1, write the opposite, flip state, move both heads right", so it writes the
complement of whatever it reads, in one pass.

`net.extras["seed"]` lists the starting strings. The closure has no
counts or rates; it answers the question "what can these strings do to each
other?" Pass other strings with `strings=[...]`; any binary strings of up to
256 symbols work.

**With and without recognition.** Seeded with the replicator alone, the
closure under recognition has 4 species and 5 reactions, and the second half
of the copying cycle is missing: the replicator's complement has a single
nine-symbol recognizon, `#1####00#`, and no recognizon of the replicator has
that length, so the two never meet. Switch recognition off and every ordered
pair reacts. The cycle then closes (4 species, 12 reactions):

```python
net = chemart.generate_network("mccaskill-polymer-tm",
                               strings=["0001111100111010111"], recognition="none")
```

```
p0001111100111010111 + p1110000011000101000 -> 2 p0001111100111010111 + p1110000011000101000
```

`recognition="none"` is a Chemart addition, a well-mixed limit with no
specificity. It is not in the report.

**The soup.** `method="soup"` runs the report's collision algorithm on a
population of `population` molecules (200 by default; the report used 1,000),
inoculated with `inoculum_fraction` (10%) of each seed string and filled with
random strings of length `l` (19). `steps` is the number of molecules drawn.
The network then holds only the reactions that fired, each with its `count`,
and `net.extras` holds the rest:

```python
net = chemart.generate_network("mccaskill-polymer-tm", seed=1, method="soup")
net.extras["analysis"]             # {'steps': 5000, 'recognition_collisions': 880}
fs = net.extras["final_state"]     # molecule counts at the end, most common first
list(fs.items())[:3]               # [('p0000000000000000000', 82), ('p1111111111111111111', 26), ('p1', 20)]
fs.get("p0001111100111010111", 0)  # 0: the replicator started at 20 and died out
```

`net.initial_state` is the starting population and `net.outflow` is
`constant-total`: the population size never changes. This run has 212
species and 197 reactions and takes about 3 seconds.

**The report's run, at full size.** The report's experiment was 1.2 million
steps in a population of 1,000:

```python
net = chemart.generate_network("mccaskill-polymer-tm", seed=1, method="soup",
                               steps=1_200_000, population=1000)
net.extras["analysis"]["recognition_collisions"]   # 174065  (the report: 177104)
```

This takes about 20 seconds. The number of collisions is close to the
report's, but the outcome is not (see Results): the replicator and the
parasite both go from 100 copies to none, and the population ends up mostly
the short strings `1` (582 molecules) and `0` (253).

**Mutation.** `error_rate` (soup only) is the probability that a written
symbol is flipped. With `error_rate=0.01, steps=20000` the soup reaches 708
species.

## Results

**A self-replicating string exists.** The first question of the report's
section 4 is whether the chemistry has a self-replicator at all, and it gives
one: `0001111100111010111`, "a member of a larger family of such strings".
Chemart reproduces it, as replication through the complement, and its tests
check both halves: the replicator reading itself writes its complement, and
reading its complement writes itself. The report also estimates that about
2^17 random strings would have to be sampled to find a member of the family,
"i.e. ca 10^8", which is inconsistent (2^17 is about 1.3 × 10^5). Chemart does
not reproduce either figure: under its decoding, 10,133 of the 2^19 strings of
length 19 replicate through their complement, about one in 52 (implementation
decisions).

**Parasites exist.** The second question is the stability of the replicator
against "the simplest kind of parasites : other single strings which have
the same recognizon but do not encode the replication rule". The report's
parasite is `0000000100100010111`. Chemart reproduces it, and the tests check
that it encodes no rule, that the replicator copies it through its
complement, that it copies nothing itself, and that its recognizon collides
with the replicator's.

**The first evolution experiment.** McCaskill inoculated a random population
of strings of length 19 with 10% each of the replicator and the parasite and
followed it for 1.2 × 10^6 steps, with 177,104 recognition collisions. Both
strings first grew and gathered "quasispecies-like distributions of related
strings around them". Then a slow decline set in. The parasite's master
sequence went extinct at 750,000 steps. The master replicator nearly died
out too and was replaced by another self-replicating string; by 1.2 million
steps neither original string was left, and the replicating rule, now
carried by various sequences, was about half as common as at the start, with
related rules competing. Still, "the replicating rule turned out to be a long
term survivor in the population". McCaskill's explanation: the replicators
"mutated away from the parasitic recognizon and diversified apparently to
such an extent that the parasites could not coadapt and went extinct". This
is the result the book and the review summarise as an evolutionary arms race.
(The book tells it after the hardware versions; the report shows it came
from the 1988 software.)

Chemart does **not** reproduce this run. With its reconstructed recognizons
the replicator can meet itself but never its complement, so nothing ever
writes a new replicator, and the inoculated copies die out: in the test run
(seed 0, 20,000 steps) the replicator falls from 20 to 0. With recognition
switched off, the full copying cycle does fire in the soup (seed 1, 5,000
steps: 4 times), but in both modes the replicator still dies out on the seeds
tried, because processors in the random background keep writing short junk
such as `1`, and every product displaces a random molecule. The tests check
exactly these behaviours: the replicator reading itself fires, reading its
complement does not under recognition, the replicator declines, and the cycle
fires without recognition. The decisions trace the gap to details the report
leaves out: where rules and recognizons start relative to their initiator,
where the heads start, how the written tape is cut into strings. Chemart's
choices are the only literal ones found that make the published replicator
replicate, but the recognizon reading frame is the weakest part of the
reconstruction.

**Later work, in hardware.** McCaskill's group went on to build
special-purpose parallel computers from reconfigurable chips (FPGAs): POLYP
(Tangen, Schulte and McCaskill, 1997) and NGEN (McCaskill et al., 1997). The
NGEN paper describes the earlier work as "networks of interacting molecular
Turing machines", and credits it with showing that molecular processing can
be broken into "binary recognition events and unary processing". The review
of Dittrich et al. (2001) reports that in these later versions, on a lattice,
sets of cooperating polymers evolved that interact "in a hypercyclic
fashion", and that later still "a Chemoton-like cooperative behavior
appeared", with membrane-bounded organisations that assemble their own
membrane (compare [the chemoton](chemoton.md)). Chemart models neither the
lattice nor the hardware, so these results are out of its scope. The book
also lists Breyer, Ackermann and McCaskill (1998) and Ehricht, Ellinger and
McCaskill (1997) for this chemistry; the first simulates other chemistries in
thin films and the second is a wet-lab paper, so neither describes the
pattern processing chemistry itself.

**What else Chemart leaves out.** The report sketches an alternative
processor for readers who reject the split into processor and tape: two
strings crossed like two tapes, with no state machine, both rewritten where
they cross, and the rewritten strings "returned to the population in addition
to the initial strings". Chemart implements only the main, catalytic
processor. Note that the book's reaction scheme, `s1 + s2 -> s1 + s3`, in
which the tape is used up, matches neither processor as the report describes
it. Chemart also replaces the report's single array of symbols, in which new
strings partly overwrite old ones, with whole molecules replaced at random
(the report itself names this "Moran type of displacement" as an equivalent
way to bound the population).

## Further reading

- Tangen, U., Schulte, L. & McCaskill, J. S. (1997). A parallel hardware
  evolvable computer POLYP. *Proceedings of the 5th Annual IEEE Symposium on
  Field-Programmable Custom Computing Machines*, 238–239 (book ref [843]).
