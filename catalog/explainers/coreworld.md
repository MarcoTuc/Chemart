## Introduction

Coreworld is a computer memory turned into a test tube. Steen Rasmussen,
Rasmus Feldberg, Morten Hindsholm and Carsten Knudsen, working at the
Technical University of Denmark and at Los Alamos, built it in 1989 and
published it in *Physica D* in 1990. Their simulator is called VENUS
("Virtual Evolution in a Non-deterministic Universe Simulator"), and the
memory it simulates is the *core*: a ring of a few thousand cells, each
holding one machine instruction. Small groups of instructions that read,
copy and overwrite their neighbours play the part of molecules, and the
execution of one instruction on another plays the part of a reaction.

The starting point is Core War, the programming game described by A. K.
Dewdney in 1984, in which programs written by people fight for a shared
circular memory ([Core War](corewar.md) in this catalog). Rasmussen and
colleagues kept its instruction language, Redcode, but changed the game. The
goal was no longer to kill the opponent. It was to see whether new
functional structures would appear on their own, from a random core, with no
evolutionary path laid out in advance. They said as much in the name: Core
War's simulator is MARS, so theirs is VENUS, "a very different game". To
make the memory behave more like a chemistry, they added three things
Core War does not have: a supply of *computational resources* that
execution uses up and that flows back in at a fixed rate, a limit on how far
an instruction can reach, and constant random noise.

What came out depended on the resource supply. With little resource (a
"desert"), the core settled into simple loops and pointers stuck on single
instructions. With plenty (a "jungle"), and a long enough reach, it went
through a series of stages and in some runs ended up dominated by large,
noise-resistant structures made of copy (`MOV`) and split (`SPL`)
instructions, which the authors called cooperative "organisms". Hand-written
programs, by contrast, were too fragile to survive the noise.

Coreworld is a simulation model. Banzhaf and Yamamoto place it in their
chapter on automata and machines (book §10.6.3), between Core War (§10.6.2)
and [Tierra](tierra.md) (§10.6.4). What separates it from Core War is the
resource, the locality and the noise, which make it an open-ended experiment
instead of a tournament. What separates it from Tierra is the problem the
book stresses: in Coreworld nothing marks where one organism ends and the
next begins, because any instruction can overwrite any other and execution
pointers wander freely between them. Tierra solved this with protected
memory. The book also notes that there are no explicitly defined reactions
in Coreworld; "the execution of code belonging to another molecule/individual
can be seen as a reaction", and Chemart turns exactly those executions into
a reaction network.

## How it works

### The core and its words

The core is a ring of cells (3,584 in the paper); the address after the last
wraps round to the first. Every cell holds one *word*: an instruction with
two operands, `A` and `B`. There are ten instructions, from the paper's
Table 1:

| word | what it does |
|---|---|
| `DAT` | data only; a pointer that tries to execute it dies |
| `MOV A, B` | copy A to B |
| `ADD A, B` / `SUB A, B` | add A to B, or subtract it |
| `JMP A` | send the pointer to A |
| `JMZ A, B` / `JMN A, B` | jump to A if B is zero / is not zero |
| `DJN A, B` | decrement B, then jump to A if it is not zero |
| `CMP A, B` | skip the next word if A and B are unequal |
| `SPL B` | split: the pointer continues, and a new pointer starts at B |

Each operand is a number with an *addressing mode*. All addresses are
relative to the executing word: `$1` means "the next cell" and `$-2` "two
cells back". `#5` is the plain number 5 (immediate). `@3` is indirect: go 3
cells ahead, read the number in that cell's `B` field, and go that much
further. `<3` is the same, but that number is decreased by one first.
Where Table 1 leaves details open (evaluation order, what `ADD` does to a
whole word), Chemart follows the 1988 Core War standard; the implementation
decisions below list every choice.

### Pointers execute the words

A word does nothing until an *execution pointer* lands on it. VENUS keeps a
queue of pointers, at most `L` of them (220 in the paper). In one *core
update*, the paper's unit of time, every pointer in the queue executes the
word it is on and then moves to the next cell, unless the word sent it
elsewhere. `DAT` kills a pointer and `SPL` adds one, if the queue has room.

The simplest program in this world is the one Core War players call the
Imp, `MOV $0, $1`: copy yourself into the next cell. The pointer then moves
to that cell, finds a fresh copy, and repeats, so the Imp crawls through
memory leaving copies of itself behind (the paper's fig. 3a).

### Resources, reach and noise

These are the three additions that make VENUS a chemistry rather than a game.

- **Computational resources.** Every cell holds an amount of resource `r`,
  counted in *execs*, where one exec pays for one instruction. A pointer may
  execute only if the cells within the *resource radius* `R_res` around it
  (that is, `2 R_res + 1` cells) hold at least one exec in total. The
  execution then takes exactly one exec out of those cells, from each in
  proportion to what it holds. A pointer that cannot pay waits where it is.
  After the whole queue has run, every cell gains `Δr`, but never beyond a
  ceiling `r_max`, which is always less than one exec. So with the paper's settings a
  lone pointer can always run, but many pointers crowded into one place compete for a limited
  local supply. The paper's Table 2 defines a **jungle** as `R_res = 3`,
  `Δr = 0.5` and a **desert** as `R_res = 5`, `Δr = 0.1`, both with
  `r_max = 0.5`.
- **Operation radius.** An instruction may read, write or send a pointer no
  further than `R_opr` cells from itself. This keeps interactions local.
- **Noise.** Two kinds. When a `MOV` executes, with probability `P_mut` the
  copied word is replaced by a random one (each of the ten instructions
  equally likely, random operands). And with probability `P_point` per core
  update a new pointer appears at a random cell. Both are 0.05 in Table 2.
  The injections matter because a random core is full of `DAT`s that kill
  pointers; without them all activity could die out.

The paper updates the core **in parallel**: every pointer of an update sees
the core as it was at the start, and the writes are applied afterwards; when
two write the same cell, the pointer later in the queue wins. The book calls
this Venus I, and describes a sequential variant, Venus II, in which each
pointer's writes are visible to the next pointer at once. Chemart offers
both.

### From executions to reactions

Coreworld has no molecules in the usual sense; the "organisms" are groups of
words that happen to work together. Chemart therefore reads the chemistry off
the executions. Every time an instruction writes into the core, it records a
reaction: the reactants are the executing word, the word it copies or adds
from, and every word it overwrites; the products are the same cells after the
write. The executing word and the source come out unchanged, so they act as
catalysts. The number of cells never changes, so every reaction has as many
products as reactants.

What counts as a species is a choice. By default a species is an instruction
type (`MOV`, `SPL`, ...), which gives ten species and matches what the paper
tracks: how many cells of the core hold each instruction. With
`species="word"` a species is a whole word such as `MOV_$0_$1` (the word
`MOV $0, $1`). A write that does not change the species count, for instance
an `ADD` that only changes a number, is called *elastic*: it is counted, but
not listed as a reaction.

### A worked example: the copy loop

The paper's fig. 3b is a five-word loop that copies one instruction over and
over. This run seeds it into a random core of 1,024 cells, with the noise
off, for four updates:

```python
net = chemart.generate_network("coreworld", seed=1, seed_program="copy-loop",
                               species="word", updates=4,
                               mutation_rate=0, pointer_injection_rate=0)
```

The loop sits at cells 186 to 190, and the pointer starts on the second
word:

```
186  MOV #4, $10     the word to be copied
187  MOV $-1, @3     copy the word one back to the cell named by the DAT's number
188  ADD #1, $2      add 1 to the DAT's number
189  JMP $-2         back to 187
190  DAT #2          the counter
```

Update 1: `MOV $-1, @3` at 187 copies cell 186 to the cell found through
190 (three ahead), whose number is 2, so to cell 192. That cell held a random
word, `MOV $-249, #-51`. Chemart records

```
MOV_#4_$10 + MOV_$-249_#-51 + MOV_$-1_@3 -> 2 MOV_#4_$10 + MOV_$-1_@3
```

with the copying word and the source as catalysts. Update 2: the `ADD` turns
the counter into `DAT #3`, which is the reaction
`DAT_#0_#2 + ADD_#1_$2 -> DAT_#0_#3 + ADD_#1_$2`. Update 3: the `JMP` sends
the pointer back to 187. Update 4: the `MOV` copies to cell 193, over the
random word `JMZ $-82, $232`. Every three updates the loop adds one more copy,
which is why the paper calls such loops "very powerful in multiplying any
instruction".

The same run with the default `species="instruction"` lists only one of
these three events, `JMZ + 2 MOV -> 3 MOV`. The first copy replaced one `MOV`
with another, and the counter stayed a `DAT`, so at the level of instruction
types both were elastic.

## Using it

The default run (above) is a jungle: Table 2's resource and noise settings,
the self-replicating seed program MICE, and parallel updating (Venus I), on a
smaller core than the paper's (1,024 cells, 64 pointers, 2,000 updates,
`R_opr = 100`) so that it finishes in about two seconds. MICE is an
eight-word Core War program by Chip Wendell that the paper lists in its
appendix. The paper used it to seed the random core with copy loops; its
first copy lands 713 cells away after 18 updates.

The network is what happened in this one run: 740 distinct reactions, each
with a `count` of how often it occurred. The rest of the run is in
`net.extras`:

```python
a = net.extras["analysis"]
net.initial_state["ADD"], a["final_composition"]["ADD"]   # (106, 383)
a["final_pointers"], a["self_loop_pointers"]              # (60, 24)
a["events"]["writes"], a["events"]["elastic_writes"]      # (32286, 29672)
a["events"]["splits"], a["events"]["deaths"], a["events"]["waits"]   # (962, 980, 37258)
```

In this run the random core, about a hundred cells of each instruction, has
been partly taken over by `ADD` after 2,000 updates, the kind of
amplification of a single instruction that the paper attributes to copy
loops. Of the 60 pointers left, 24 sit on a word
that sends them back to itself. `waits` counts the times a pointer could not
execute for lack of resource.

The other fields:

- `a["series"]` samples the run about 100 times: `pointers`, `executions`
  per update, `composition` (cells per instruction) and `mean_resource`. In
  the jungle `Δr` equals `r_max`, so every cell is full again after each
  update and `mean_resource` stays at 0.5.
- `a["executed"]` counts executions per instruction; `a["final_run_lengths"]`
  counts how often each instruction occurs alone, in pairs, triples and so
  on, the statistic of the paper's fig. 4.
- `net.extras["coreworld"]` holds the final core as text (`final_core`), the
  final pointer positions and where the seed was placed (`seed_origin`).

**Desert and jungle.** Change `resource_radius` and `resource_influx` to
Table 2's values:

```python
jungle = chemart.generate_network("coreworld", seed=1)
desert = chemart.generate_network("coreworld", seed=1, resource_radius=5, resource_influx=0.1)
```

The desert executes about 14 instructions per update in the second half of
the run against 37 in the jungle, splits 156 times against 962, and its final
composition stays close to the random start (no instruction above 144 cells).

**The paper's scale.** `core_size=3584, queue_len=220` is the paper's core;
`operation_radius=3584` removes the reach limit, the setting the paper found
most interesting. Each of these runs of 5,000 updates took 5 to 7 seconds.
Run time grows with the number of updates and of active pointers; runs of
the paper's length (100,000 updates and more) were not timed for this page.

```python
net = chemart.generate_network("coreworld", seed=2, core_size=3584, queue_len=220,
                               updates=5000, operation_radius=3584)
```

Three seeds went three ways. With seed 2 the core ended with 1,889 `SPL` and
1,051 `MOV` out of 3,584 cells, having changed little after about update
3,500. With seed 3 it filled with 1,711 `DJN` and 1,294 `JMN`, and 211 of
the 220 pointers were trapped on self-looping words. With seed 1 the
composition was still near uniform (at most 419 cells of one instruction).

**Brittle programs.** With the same paper-size settings and 600 updates,
MICE left 9 intact copies of its six executable words with the noise off
(`mutation_rate=0, pointer_injection_rate=0`) and none with Table 2's noise.

**Other settings.** `mode="venus-ii"` switches to sequential updating;
`seed_program` also takes `jmp` (the paper's other start, an active
`JMP $0` that touches nothing), `imp`, `copy-loop` or `none`. The default
run with `species="word"` has 24,431 species and 24,643 reactions.

## Results

The paper is qualitative. Rasmussen and colleagues learnt to read their
cores by printing them ("each hardcopy of a core occupies several paper
meters") and naming the structures they saw. A typical simulation executed
about 22 million instructions over 100,000 core updates and took twelve
hours on an IBM personal computer. The observations below are from the 1989 Los
Alamos preprint of the *Physica D* paper, which is the version Chemart was
built from; the published version was not consulted.

**Deserts stay simple.** In a desert with a small operation radius and no
seed, the core became "relatively stable" with many pointers stuck at
*fixed points*, some simple loops and a few more complicated ones. The fixed
points are words that send the pointer back to themselves: `JMP #X`,
`JMZ #X, 0`, `JMN $0, Y` and `DJN #X, Y` with `Y` not zero. Seeding copy
loops changed the desert's instruction distribution only "a little bit". The
authors could not develop "any really viable structures" in a desert. The
book summarises this as simple cyclic structures dominating under low
resource influx. Chemart's tests check that each of these fixed points traps
its pointer, and that on Table 2's settings a jungle sustains more than 1.5
times the executions of a desert, and more splits; they do not check what a
desert core looks like after a long run.

**Jungles grow SPL-fields.** In the jungle the new feature was clusters of
dense programs driven by pointers from `SPL` instructions, which the authors
called "SPL-fields" and never found in the desert, because of their high
density of pointers. Chemart does not detect or test for these structures.

**Copy loops and phase transitions.** Loops like the one in the worked
example can flood the core with whatever instruction they copy. After about
100,000 updates in jungles with the reach limit removed, the paper reports
cores holding 2,859 `CMP` (after 110,000 updates), 1,276 `SPL` with only 42
`MOV` (182,000), 1,471 `MOV`, and 1,323 `ADD` (145,000) out of 3,584 cells,
and it calls these changes of composition phase transitions. Collapses to a
quiet core usually went with the near extinction of `MOV`, the instruction
behind both the mutations and the rearranging of programs. In one run with
`R_opr = 2000` and `Δr = 0.25`, 212 of the 220 pointers were trapped at fixed
points after 435,000 updates. Chemart's tests check that the fig. 3b loop
adds one copy every three updates (100 copies in 300 updates) and that the
network records each copy as a reaction; the floods themselves are not
tested, though short runs show them (above).

**Cooperative organisms.** With `R_opr = 100` and a MICE seed, after some
ten thousand updates the jungle often held stationary "organisms" of 20 to
100 words, mostly `SPL`, kept robust by the many pointers they produce. An
operation radius "at least of the order hundred" was needed for this. In
another run the jungle took about 100,000 updates to develop two huge
organisms of `SPL` and `MOV`, with more than 800 copies of each instruction,
over a third of the core. The authors describe the combination as extremely
stable: the `MOV`s copy `MOV`s or `SPL`s, and the `SPL`s hand out pointers
to `SPL`s or `MOV`s, so noise is damped. Such an organism does not make a
true copy of itself; it expands as it moves through the core. The authors
compared it to a cooperative, "probably autocatalytic" chemical network and
concluded that it seemed easier to evolve a "metabolic network" than a clean
genetic system. Chemart does not identify organisms, so this result is not
reproduced or tested. Its runs can end SPL- and MOV-rich, as seed 2 above
did, but nothing checks whether that core behaves like the paper's
organisms.

**Four epochs.** In the jungle with a large reach and a seed, the core went
through four successive "macroscopic states" (fig. 7): the random core with
its self-replicating seed; a core populated by copy loops; areas of
identical instructions written by those loops; and cooperative programs
taking over. The authors did not claim open-ended evolution: the process
"gets trapped after three successions", and they had "to some extent
'designed' the first succession" with the seed. Chemart does not test this
sequence.

**Engineered programs are brittle.** Programs written by people, both their
own and Core War warriors, "were too brittle to survive in the noisy VENUS
universe", while the system went on to develop its own stable structures.
Chemart reproduces this: its test places MICE in a random paper-size
jungle, and after 600 updates at least ten intact copies remain without
noise and none with Table 2's noise.

**Other parts of the paper.** The authors proposed a crude complexity
measure (the number of interacting instructions, times the number of
different instructions minus one, times the fraction of perturbations that
are not fatal), under which the Imp and MICE score close to zero and the
SPL-MOV organism scores high. They also noted that parameter sets differ:
some always lead to the same features whatever the random core, others are
extremely sensitive to it. And they discussed what their cores suggest about
computer viruses and immune systems. Chemart implements none of these.

**What Chemart checks exactly.** Beyond the phenomena, the tests reproduce
the paper's figures 3a and 3b step by step, check that MICE copies itself
713 cells away and starts a pointer on the copy in exactly the 18 updates
the paper states, and check the resource equations, the operation radius,
the parallel update rule and both noise sources. The Venus II mode was
cross-checked against the Core War simulator pMARS on 477 random programs
(see the implementation decisions).

**Later work.** Banzhaf and Yamamoto write that Redcode systems had
difficulty becoming really evolvable, and that later work took hints from
Tierra back into Redcode, citing Sperl (2011) and Vowk, Wait and Schmidt
(2004). They also credit Pargellis with demonstrating the spontaneous
emergence of a self-replicating program in Coreworld-like systems. The
sequel to Coreworld, Rasmussen, Knudsen and Feldberg's "Dynamics of
programmable matter" (Artificial Life II, 1992), is the source the book cites
for Venus II; Chemart could not obtain it.

## Further reading

- Rasmussen, S., Knudsen, C. & Feldberg, R. (1992). Dynamics of programmable
  matter. In C. G. Langton, C. Taylor, J. D. Farmer & S. Rasmussen (eds.),
  *Artificial Life II*, 211–291. Addison-Wesley. Book ref [696].
- Pargellis, A. N. (1996). The spontaneous generation of digital "life".
  *Physica D* 91(1–2), 86–96. Pargellis, A. N. (1996). The evolution of
  self-replicating computer organisms. *Physica D* 98(1), 111–127.
  Pargellis, A. N. (2001). Digital life behavior in the amoeba world.
  *Artificial Life* 7(1), 63–75. Book refs [650–652].
- Sperl, T. (2011). Taking the redpill: artificial evolution in native x86
  systems. arXiv:1105.1534. Book ref [789].
- Vowk, B., Wait, A. S. & Schmidt, C. (2004). An evolutionary approach
  generates human competitive coreware programs. In *Workshop and Tutorial
  Proceedings, Ninth International Conference on the Simulation and Synthesis
  of Living Systems*, 33–36. Book ref [894].
