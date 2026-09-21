## Introduction

Core War is a game in which two or more computer programs fight inside a
shared, simulated memory. A. K. Dewdney introduced it in his "Computer
Recreations" column in *Scientific American* in May 1984, after building a
first version with his student David Jones at the University of Western
Ontario. The programs, called *warriors*, are written by people in a small
assembly language, **Redcode**, and are run by a simulator called **MARS**
(Memory Array Redcode Simulator). Once a battle starts nobody can intervene:
in Dewdney's words, "people do not play at all". Each warrior tries to make
its opponent execute an instruction that cannot be executed, usually by
writing such instructions ("bombs") over the opponent's code. A warrior that
runs into one dies.

Two tiny warriors from Dewdney's article give the flavour. The **Imp** is a
single instruction, `MOV 0 1`: "copy this cell into the next cell". Executing
it copies it one cell ahead, where it is executed next, so the Imp crawls
through memory one cell per turn, leaving a trail of copies of itself. The
**Dwarf** is four instructions long. It sits still and drops a bomb every few
cells, looping around the whole memory. When the Imp crawls over the Dwarf,
the Dwarf's own jump lands on an Imp instruction, and the Dwarf turns into a
second Imp: neither dies, and the battle is a draw.

Core War is not a chemistry by itself: it was designed as a game, it has no
rates, and it has no mutation, so nothing evolves on its own. It is in the
catalog because it is the ancestor of a line of *artificial chemistries* built
from self-modifying machine code. Banzhaf and Yamamoto present it (book
§10.6.2) as the first *assembler automaton*, a parallel machine in which many
processes execute instructions from, and write into, one circular memory, and
they reproduce its instruction set "since it is so foundational". The entries
that follow it in the book reuse the idea: [Coreworld](coreworld.md) adds a
flowing resource, random new processes and copy errors to Core War's
instructions so that programs can arise and evolve; [Tierra](tierra.md) and
[Avida](avida.md) change the instruction set and add memory protection and
mutation to make self-replicating digital organisms. Core War is the plain
game underneath: its warriors are hand-written and its rules fixed by a
standard.

What Chemart adds is a *chemical reading* of a battle. Every instruction
sitting in memory counts as a molecule, and every executed instruction as one
reaction event that consumes and rewrites the cells it touches. The book
itself remarks, of Coreworld, that such systems have no explicitly defined
reactions and that executing another program's code "can be seen as a
reaction"; the precise event reading here is Chemart's own. Chemart's entry is
a full simulator, a port of the reference program pMARS that reproduces
pMARS's results exactly, and it turns each battle into an observed reaction
network.

## How it works

### The core, the warriors and their tasks

The memory is the **core**: a ring of `M` cells, numbered 0 to `M − 1`, where
cell `M` is cell 0 again. Every cell holds exactly one Redcode instruction.
An empty core is filled with `DAT $0, $0`. All addresses are *relative*: `1`
means "the next cell", `-2` "two cells back". A warrior cannot know where in
the core it sits.

A Redcode instruction has an *opcode* (what to do), two operands called A and
B, and, in the modern standard, a *modifier* after a dot that says which parts
of the cells the instruction reads and writes (`.I` means the whole
instruction, `.AB` means "from the A-number of the source to the B-number of
the target", and so on). Each operand carries an *addressing mode*: `#` is an
immediate number, `$` a direct address, `@` an indirect address (go to the
cell, read its B-number, and use it as a further offset from that cell).
The book's Table 10.2 lists the opcodes: `DAT` (data; executing it kills the
process), `MOV` (copy), `ADD`, `SUB`, `MUL`, `DIV`, `MOD` (arithmetic; dividing
by zero kills the process), `JMP`, `JMZ`, `JMN`, `DJN` (jumps, some
conditional), `CMP`/`SEQ`, `SNE`, `SLT` (skip the next instruction if a
comparison holds), `SPL` (split: start a second process), `NOP`, and `LDP`/`STP`
(load from and store to a small private memory, *P-space*, that survives from
one round to the next).

Each warrior owns a queue of **tasks** (processes). A task is just a program
counter, the address of the next instruction it will execute. A warrior starts
with one task. `SPL` adds a task, up to a limit; executing `DAT` or dividing by
zero removes one. A warrior with no tasks left is dead.

### The simulator

MARS loads the warriors at random positions, at least a minimum distance
apart, and runs in *cycles*. In each cycle it takes the next task from each
living warrior's queue, in turn, and executes one instruction for it; the task
then goes to the back of its queue at its new address (or two tasks do, after
`SPL`; or none, after `DAT`). So a warrior with many tasks does not run
faster: its instructions are shared out among its tasks. A **round** ends when
only one warrior is left, which wins, or after a fixed number of cycles, when
all survivors tie. A battle is several rounds, with fresh positions and the
first mover rotating.

### Molecules and reactions

Chemart reads this as a chemistry with two kinds of molecules:

- **Instructions in core.** The species name is the instruction written
  without spaces: `MOV.I$0,$1` is the Imp, `DAT.F$0,$0` an empty cell. The
  core always holds exactly `M` of them, so their total count is conserved.
- **Tasks.** One species per warrior, named from the warrior's name and a
  short hash of its code, such as `Imp_115f417a`. Its count is the number of
  tasks in that warrior's queue.

One executed instruction is one reaction. Its reactants are the task that
executes and the instructions in every cell the step reads or writes; its
products are the task or tasks that go back into the queue and the same cells
after the step. Cells the step only reads come out unchanged, so they appear
on both sides, as *catalysts*. Positions are not part of the species: the
network records *what* happened, not where. Nor is anything a rate: a
reaction happens when a task reaches the instruction, in the order the
simulator fixes.

### A worked example: Imp against Dwarf

Here are reactions from the default network (printed further down under
"Using it in Chemart"), which is a battle of Imp against Dwarf. First the Imp
replicating, the reaction that fired 2,312 times:

```
Imp_115f417a + DAT.F$0,$0 + MOV.I$0,$1 -> Imp_115f417a + 2 MOV.I$0,$1
```

The Imp's task executes the cell holding `MOV.I $0, $1`, which copies the whole
instruction (`.I`) from offset 0 (itself) to offset 1 (the next cell). The next
cell held an empty `DAT.F $0, $0`; afterwards it holds a second `MOV.I $0, $1`.
One task goes back into the queue, now pointing at the new copy. One empty cell
has become an Imp cell: that is replication, written as a reaction.

The Dwarf used in Chemart is the modern form from Karonen's *Beginners' Guide
to Redcode*:

```
ADD.AB #4, 3      ; add 4 to the B-number of the cell 3 ahead (the DAT)
MOV.I  2, @2      ; copy that DAT to the cell its B-number points at
JMP    -2         ; back to the ADD
DAT    #0, #0     ; the bomb, and the pointer
```

Its first two reactions in the default run are

```
Dwarf_20242a64 + DAT.F#0,#0 + ADD.AB#4,$3 -> Dwarf_20242a64 + DAT.F#0,#4 + ADD.AB#4,$3
Dwarf_20242a64 + DAT.F#0,#4 + DAT.F$0,$0 + MOV.I$2,@2 -> Dwarf_20242a64 + 2 DAT.F#0,#4 + MOV.I$2,@2
```

In the first, `ADD` (a catalyst, unchanged) turns the Dwarf's `DAT #0, #0`
into `DAT #0, #4`. In the second, `MOV` copies that `DAT` to the cell 4 beyond
it, an empty cell, which becomes a bomb: a copy of `DAT #0, #4`. Then `JMP`
loops back (`Dwarf_20242a64 + JMP.B$-2,$0 -> Dwarf_20242a64 + JMP.B$-2,$0`,
the task alone moving on), and the next pass bombs 4 cells further. Because
the core size 800 is divisible by 4, after 200 passes the pointer comes back to
the Dwarf's own `DAT` and the Dwarf never bombs its own code.

A task that walks into a bomb dies. In round 4 of the default run the Imp ran
into one of the Dwarf's bombs:

```
Imp_115f417a + DAT.F#0,#-304 -> DAT.F#0,#-304
```

The task is consumed and nothing replaces it; the Imp has lost. And when the
Imp gets to the Dwarf first, the Imp writes `MOV.I $0, $1` over the Dwarf's
code. The Dwarf's task then executes Imp instructions, so reactions appear in
which the Dwarf's task runs `MOV.I$0,$1`: the Dwarf has been *infected* and is
now a second Imp. Chemart counts such executions of code written by another
warrior as `foreign_executions`.

## Using it

The default call runs Imp against Dwarf for four rounds in an 800-cell core,
with random load positions. The parameters are the standard King of the Hill
(KOTH) settings of the 1994 standard, the ones used by the Internet
tournaments, with the core size, the cycle limit and the task limit divided by
10 so that a battle runs in about a second. With `seed=1`, the first three
rounds are ties (the Imp reached the Dwarf and converted it) and in the fourth
the Dwarf bombed the Imp:

```python
a = net.extras["analysis"]
[r["survivors"] for r in a["rounds"]]
# [['Imp_115f417a', 'Dwarf_20242a64'], ['Imp_115f417a', 'Dwarf_20242a64'],
#  ['Imp_115f417a', 'Dwarf_20242a64'], ['Dwarf_20242a64']]
a["rounds"][3]["death_step"]     # {'Imp_115f417a': 742}
a["wins"], a["ties"]             # ({'Imp_115f417a': 0, 'Dwarf_20242a64': 1}, {'Imp_115f417a': 3, 'Dwarf_20242a64': 3})
a["score"]                       # {'Imp_115f417a': 3, 'Dwarf_20242a64': 6}
```

The score is pMARS's default: 3 points for a win and 1 for a tie between two
warriors. Each entry of `a["rounds"]` also gives the load `positions`, the
`starter`, the number of `steps` executed and the tasks left.
`a["events"]` counts, per warrior, the instructions executed, successful
splits, task deaths, `foreign_executions` (instructions executed that another
warrior wrote), cell writes and writes over another warrior's cells. Here the
Dwarf executed 19,079 foreign instructions out of 24,371: most of its life was
spent as an Imp. `net.extras["warriors"]` holds each warrior's assembled code,
and `net.extras["final_state"]` the contents of the core at the end of all
rounds, summed; in this run 2,678 of the 3,200 cells are `MOV.I$0,$1`.

The reaction network counts each distinct event once with its firing count, so
the 771 reactions of the summary are 771 *different* events. Most of them are
Dwarf bombing steps that differ only in the pointer value of the bomb
(`DAT.F#0,#4`, `DAT.F#0,#8`, ...), which is why the species count grows with
the run: every new bomb value is a new species.

**Imp converts Dwarf.** Fix the Dwarf's load address with `positions` and run
one round:

```python
net = chemart.generate_network("corewar", seed=1, positions=[300], rounds=1)
a = net.extras["analysis"]
a["rounds"][0]["survivors"]        # ['Imp_115f417a', 'Dwarf_20242a64']   a tie
a["events"]["Dwarf_20242a64"]["foreign_executions"]     # 7700 of 8000
```

The Imp overwrites the Dwarf's three active cells (reactions such as
`Imp + MOV.I$0,$1 + JMP.B$-2,$0 -> Imp + 2 MOV.I$0,$1`), and from then on both
tasks run `Imp + 2 MOV.I$0,$1 -> Imp + 2 MOV.I$0,$1`, copying an Imp onto an
Imp, for the rest of the round.

**The imp avalanche.** The book's three-line program `SPL 2 / JMP -1 /
MOV 0, 1` is built in as `"imp-avalanche"`:

```python
net = chemart.generate_network("corewar", seed=1, warriors=["imp-avalanche"],
                               max_processes=16, max_cycles=2000, rounds=1)
```

It has a single reaction that makes tasks,
`Imp_avalanche_d7f1d8f2 + SPL.B$2,$0 -> 2 Imp_avalanche_d7f1d8f2 + SPL.B$2,$0`,
which fired 15 times: the queue grew from 1 task to the limit of 16. Each
split starts a new Imp. Because the spawning task has to share turns with
every Imp it has started, new tasks come more and more slowly: after 9 cycles
there are 4 tasks, after 25 there are 6, after 49 there are 8.

**Your own warriors.** `warriors` takes Redcode source as well as the built-in
names `imp`, `dwarf` and `imp-avalanche`. The assembler accepts the 1994
format: labels, `EQU` constants, `ORG`/`END`, expressions and predefined names
such as `CORESIZE`. It does not accept `FOR`/`ROF` blocks, so some published
warriors must be unrolled by hand. Dewdney's original 1984 Dwarf, for example,
assembles, but under the modern rules it no longer kills (see the
implementation decisions):

```python
old = "ORG 1\nDAT -1\nADD #5, -1\nMOV #0, @-2\nJMP -2\n"
net = chemart.generate_network("corewar", seed=1, warriors=[old, "imp"], rounds=4)
net.extras["analysis"]["wins"]     # {'warrior_663eb73a': 0, 'Imp_115f417a': 0}   four ties
```

**Tournament scale.** `core_size=8000, max_cycles=80000, max_processes=8000`
are the full KOTH settings. One round at that scale took about 2.5 seconds
here, both for Imp against Dwarf (a 160,000-step tie) and for Rave, a
published warrior used in the tests, against Dwarf; the Rave battle produced a
network of 22,753 species and 30,852 reactions. A KOTH-style match of 100 or
more rounds therefore takes minutes. With more than two warriors, the extra
ones are placed at random, as in pMARS.

## Results

Core War's results are those of a game rather than of a scientific study: a
set of rules, a standard, and a body of warriors and strategies. The findings
below are what the book and the sources report, with what Chemart reproduces.

**The Imp moves and leaves a trail.** Dewdney (1984) describes the Imp as
moving "through the array at a speed of one address per cycle, leaving behind
a trail of MOV 0 1 instructions". (The book prints the trail as `MOV 1 0`, a
typo.) Chemart's test runs the Imp alone for 100 cycles and checks that
exactly the first 101 cells hold `MOV.I $0, $1` and the rest are untouched.

**The Dwarf bombs the core without hitting itself.** Dewdney's Dwarf bombs
every fifth cell of an 8,000-cell core; the modern form in Karonen's guide
drops a `DAT` every fourth cell. The test runs the Dwarf alone for 200 loops in
an 800-cell core and checks that bombs lie exactly at cells 7, 11, ..., 799,
and that the Dwarf's own code is intact. Dewdney also reports that in Dwarf
against Dwarf "each program wins 30 percent of the time; in 40 percent of the
contests neither program scores a fatal hit". That figure is for his 1984
Dwarf and rules, and Chemart does not test it.

**Imp against Dwarf ends in a draw or a Dwarf win.** Dewdney predicted that if
the Imp reaches the Dwarf, "Dwarf will be subverted and become a second Imp
... the battle is a draw", and Karonen's guide adds that the Imp "won't win too
many games" because anything it overwrites becomes an Imp too. The Imp can
never win: it writes only copies of itself, never a bomb. Chemart makes this
precise. With the scaled KOTH rules in an 800-cell core, over all 601 legal
load addresses of the Dwarf, the compiled pMARS gives 104 Dwarf wins and 497
ties; Chemart gives the same outcome at every address (a test marked slow,
about 7 seconds). Other tests check four single positions against pMARS's
step counts, check that the Dwarf's task ends up executing Imp code in a
drawn round, and check that over six rounds the Imp wins none.

**The imp avalanche.** The book's `SPL 2 / JMP -1 / MOV 0, 1` "creates an
avalanche of self-replicating MOV 0 1 programs". The test checks that the task
count only grows and fills the queue to its limit, that task-making reactions
occur, and that the core ends up with more Imp cells than there are tasks.

**Death and ties.** A task dies on `DAT` or on division by zero, and a round
with several survivors at the cycle limit is a tie (Core War Guidelines, and
section 5 of the 1994 draft standard). The tests check division by zero under
several modifiers, the task limit on `SPL`, and the tie scoring.

**A conforming simulator.** Validate 1.1R, a test program distributed with
pMARS, loops for ever on a simulator that follows the standard's operand
evaluation rules and destroys itself otherwise. In Chemart it survives to the
cycle limit, and the whole core at the end matches pMARS's. The tests also
match pMARS's final core, step counts and task counts on Rave against Dwarf in
an 8,000-cell core, and on three randomly generated three-warrior battles over
four rounds, which between them use every opcode, modifier and addressing mode
and P-space. The implementation decisions report that several hundred random
battles were compared in the same way. This is the strongest claim the entry
makes: for battles under the 1994 rules, Chemart's machine is pMARS.

**What Chemart does not reproduce.** The original 1984 rules, in which any
cell whose value is zero counts as a `DAT`, are not implemented; neither are
the 1986 and 1988 standards, except as they survive in 1994, nor `FOR`/`ROF`,
`PIN` (shared P-space) or read and write limits. Evolution is outside this
entry: Core War warriors do not mutate. Evolving them was hard, because, as
Vowk, Wait and Schmidt put it, Redcode is "extremely brittle": nearly any
change to a successful program is fatal. Their 2004 system seeded random
programs with instruction statistics taken from human-written warriors,
bred them against a benchmark of strong warriors over days or weeks on a PC
cluster, and produced warriors competitive on the public tournament hills,
among them one that used two Imps and an "imp-gate", a strategy they had not
seen in human code. The book cites this and Sperl (2011) as later work that
took ideas from Tierra back to Core War. None of this is in Chemart; for
Core War instructions under mutation and a resource flow, see
[Coreworld](coreworld.md).

## Further reading

- Vowk, B., Wait, A. & Schmidt, C. (2004). An evolutionary approach generates
  human competitive Corewar programs. In *Workshop and Tutorial Proceedings,
  Ninth International Conference on the Simulation and Synthesis of Living
  Systems (ALife IX)*, 33–36 (book [894]).
  <https://corewar.co.uk/vowk/alife9ac.pdf>
- Sperl, T. (2011). Taking the redpill: artificial evolution in native x86
  systems. arXiv:1105.1534 (book [789]).
- corewar.co.uk, an archive of the standards, Dewdney's articles and warriors:
  <https://corewar.co.uk/>
