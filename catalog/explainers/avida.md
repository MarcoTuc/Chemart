## Introduction

Avida is a world of self-copying computer programs, built as a laboratory for
evolution. Each "organism" is a short program in a simple assembly-like
language, running on its own small virtual computer and living in one cell of
a square grid. The program's main job is to copy itself. Copying is slightly
error-prone, so offspring sometimes differ from their parents. Space on the
grid is limited, so every newborn eventually replaces somebody. Replication,
variation and competition are all present, so the population evolves, and
nobody tells it how.

Avida grew out of Thomas Ray's [Tierra](tierra.md), where self-copying
programs share one block of computer memory. According to Ofria, Bryson and
Wilke (2009), Chris Adami began working with Tierra in 1992. He wanted
digital organisms to evolve solutions to mathematical problems, and his idea
was to watch what numbers they read in and wrote out, and to pay them extra
CPU time when an output was, say, the sum of two inputs. Tierra proved
limiting, and in 1993 Charles Ofria and C. Titus Brown joined him to build a
new platform, Avida. Banzhaf and Yamamoto credit the system to Adami and Brown
(1994). Avida is developed at Michigan State University, and its authors and
their colleagues have used it for a long series of experiments in
evolutionary biology.

Two things set Avida apart from Tierra, and both are visible in Chemart.
First, space: organisms sit on a two-dimensional grid and a newborn goes into
a neighbouring cell, so the population size is fixed by the grid and
different lineages can hold different regions. Second, an external
environment: organisms read numbers in and write numbers out, and those that
compute certain logic functions of their inputs get more CPU time, so they
copy themselves faster. This turns Avida from a system that only evolves
faster copiers into one in which new computational abilities can evolve, and
it lets an experimenter ask how a complex ability is assembled from simpler
ones. Because every birth is recorded, the full ancestry of any organism can
be traced, with no missing links.

In the book Avida is a *simulation model* of the "assembler automata" family
(chapter 10, §10.7.1, among the chemistries based on cellular automata). As a
chemistry, a molecule is one program and a reaction is replication: the
program `s` uses up energy `E`, measured in CPU time, to make a copy `s'`,
written `s --E--> s + s'`. The book returns to Avida in its section on
ecological modelling (§8.2.3) and on evolving distributed algorithms
(§17.3.2). Its nearest neighbours in the catalog are the other programs that
copy themselves in a shared machine: [Core War](corewar.md),
[Coreworld](coreworld.md) and [Tierra](tierra.md), which have no grid and no
rewarded tasks. [Self-replicating loops](sr-loops.md) also reproduce on a
grid, but there the replicator is a pattern of cellular-automaton states
rather than a program; [aevol](aevol.md) puts digital organisms on a grid too,
but with a gene-and-protein genome instead of a program.

## How it works

### An organism is a program on a small virtual computer

An Avida genome is a circular list of instructions. Chemart uses the default
instruction set of modern Avida (the "heads" set), which has 26
instructions, so a genome can be written as a string of the letters `a` to
`z`, one letter per instruction. That string is the species name's content in
Chemart: `w` is `h-alloc`, `z` is `h-search`, `c` is `nop-C`, and so on (the
full table is in the formal specification below).

Each organism runs its genome on a virtual CPU with:

- three **registers**, `AX`, `BX` and `CX`, each holding a 32-bit number;
- two **stacks** of numbers, one active at a time;
- four **heads**, pointers into the organism's own memory: the *instruction
  head*, which marks the instruction being executed; the *read* and *write
  heads*, used to copy the genome; and the *flow head*, used for jumps and
  loops;
- an **input buffer** and an **output buffer** for exchanging numbers with the
  environment.

Instructions take no arguments. Instead, three "no-operation" instructions,
`nop-A`, `nop-B` and `nop-C`, act as modifiers. Placed right after an
instruction they change which register or head it acts on (`inc` alone
increments `BX`, `inc nop-A` increments `AX`). A run of nops also forms a
**template**, a label that other instructions can search for. Templates match
by complement: `nop-A` pairs with `nop-B`, `nop-B` with `nop-C`, and `nop-C`
with `nop-A`. Ofria et al. (2009) explain the reason for labels rather than
numeric addresses: a mutation that inserts or deletes an instruction shifts
all positions after it, but a label still finds its partner.

### How the ancestor copies itself

Every run starts from one hand-written ancestor. Chemart's default is the
ancestor that ships with Avida 2 (`default-heads.org`), 100 instructions long:

```
wzcag ccccc…ccccc zvfcaxgab
```

Five instructions at the start, 86 `nop-C` fillers that do nothing, and nine
at the end. Following the comments in Avida's own file, it works like this:

1. `w` (`h-alloc`) enlarges the organism's memory to make room for the
   offspring: from 100 to 300 lines.
2. `z c a` (`h-search nop-C nop-A`) searches for the complement template,
   `nop-A nop-B`, which is the last two instructions (`a b`), and places the
   flow head just after it: the start of the new, empty memory.
3. `g c` (`mov-head nop-C`) moves the write head there.
4. The fillers run without effect until `z` (`h-search` with no template),
   which marks the start of the copy loop.
5. The loop is `v` (`h-copy`: copy the instruction under the read head to the
   write head, and move both heads on), `f c a` (`if-label nop-C nop-A`:
   have we just copied the end label `a b`?), `x` (`h-divide`, run only if
   the answer is yes) and `g` (`mov-head`: jump back to the start of the
   loop).
6. When the whole genome has been copied, `h-divide` cuts off the copy as a
   new organism.

Chemart's `test_cpu` runs one organism alone, without mutations, until its
first division:

```python
from chemart.chemistries import avida as A
A.test_cpu(A.DEFAULT_ANCESTOR)
# {'divided': True, 'gestation': 389, 'merit': 97.0, 'executed_size': 97,
#  'copied_size': 100, 'tasks': [], 'offspring': 'wzcagcc…czvfcaxgab'}
```

The ancestor needs 389 instructions to produce an exact copy of itself: 89 of
setting up, then three for each of the 100 instructions copied. This number,
the **gestation time**, is the same as in the compiled original Avida. An
`h-divide` fails, and is simply ignored, if less than half of the parent has
been executed or less than half of the offspring's memory has been written.

### Mutations

Variation enters in two ways that the experimenter controls, as in Avida's
default configuration. Each
`h-copy` writes a random instruction instead of the right one with a small
probability (a **copy mutation**). At division, one random instruction may be
inserted into the offspring, or one deleted. A program can also copy itself
wrongly by its own logic, for instance skipping part of itself; Ofria et al.
call this an *implicit mutation*.

Here are two reactions from the default run (`seed=1`), printed further down
this page:

```
100-e994d1ba -> 2 100-e994d1ba                 (x11)
100-e994d1ba -> 100-e994d1ba + 100-eb8585d3    (x1)
```

`100-e994d1ba` is the ancestor: a species name is the genome length and the
first eight hexadecimal digits of a hash of the genome. The first line records
eleven exact copies. In the second, a copy mutation turned the `nop-C` at
position 42 (counting from 0) into `d` (`if-n-equ`), producing a new genotype of the same
length. The parent appears on both sides of each reaction because it survives
its own birth event: it acts as the catalyst of its copy.

### Merit: how CPU time is shared

The grid is simulated as a parallel computer. Time is measured in
**updates**. In one update, each living organism executes 30 instructions on
average (the time slice), but not exactly 30: each instruction is given to an
organism drawn at random, with probability proportional to its **merit**.
Twice the merit means, on average, twice the instructions, and so twice as
fast a replication. With the ancestor's 389-instruction gestation, a
generation takes about 13 updates.

Merit has two factors. The first is size: the smallest of the genome length,
the number of instructions copied into the offspring, and the number of
distinct instructions executed. The ancestor executes 97 of its 100 lines,
so its merit is 97. The second is a **bonus** for logic tasks. Each organism
has three 32-bit input numbers, and the `IO` instruction writes a register to
the output and reads the next input into it. If an output is a bitwise logic
function of the inputs, the organism has performed a **task**. There are nine
rewarded tasks, the one- and two-input logic functions, and each multiplies
the bonus by a fixed factor, at most once per gestation:

| task | computes | bonus |
|---|---|---|
| NOT | not A | ×2 |
| NAND | not (A and B) | ×2 |
| AND | A and B | ×4 |
| ORN | A or not B | ×4 |
| OR | A or B | ×8 |
| ANDN | A and not B | ×8 |
| NOR | not (A or B) | ×16 |
| XOR | A xor B | ×16 |
| EQU | A equals B, bit by bit | ×32 |

Only one instruction, `nand`, computes a logic function by itself, so every
task must be assembled from `nand`, `IO` and the instructions that move
numbers around. The rewards grow with the number of `nand` operations the
function needs (Lenski et al. 2003). Merit is set at division: the bonus a
parent earned during a gestation becomes the merit of both parent and
offspring for the next one. An organism collects its bonus for its children.

A genome that performs NOT shows the effect. Replace the start of the filler
with `yopcuyk`: `IO` reads an input `x` into `BX`, `push` and `pop nop-C` copy
it to `CX`, `nand` puts `not (x and x)`, that is `not x`, into `BX`, and the
next `IO` outputs it.

```python
NOT = "wzcagc" + "yopcuyk" + "c" * 78 + "zvfcaxgab"
A.test_cpu(NOT)
# divided True, gestation 387, merit 194.0, executed_size 97, tasks ['NOT']
```

The merit doubles, from 97 to 194.

### The grid: birth, death and competition for space

The world is a torus (a grid whose edges wrap around), 12 × 12 cells by
default. When an organism divides, its offspring goes into an empty cell among
its eight neighbours if there is one; otherwise into a random neighbour or
the parent's own cell, killing whoever lived there. An organism also dies of
old age after executing 20 times its genome length. Once the neighbourhood is
full, then, every birth is a death, and space is what organisms compete for.
Performing tasks helps only indirectly: more merit, more CPU time, faster
copying, more cells taken.

### What the network records

Chemart returns the run as the network of events between genotypes (all
organisms with the same genome form one genotype, one species). There are
four kinds of reaction:

```
P -> P + O       birth into an empty cell (O = P for an exact copy: P -> 2 P)
P + V -> P + O   birth that overwrites the organism V in a neighbouring cell
P -> O           birth into the parent's own cell: the parent dies
P -> ∅           death of old age
```

The network has no rate constants: how often each event happens follows from
the programs, the merits and the scheduler. Task completions are not
species; they are recorded in `net.extras["task_events"]`.

## Using it

The default call runs a 12 × 12 world for 150 updates from one ancestor with
Avida's default mutation rates. It reproduces no published experiment; it is a
short sample of genotype space that runs in about two seconds. In the default
run above, the ancestor divides for the first time around update 13 and the
population passes 100 organisms after about 130 updates. By update 150 there
have been 279 births, 122 of them into occupied cells, and 135 of the 144 cells
are full. 177 new genotypes have appeared and 103 are alive; the ancestor is
still the commonest, with 7 organisms. No task is performed: the ancestor
cannot do any, and 150 updates are too few to evolve one. Across seeds 1–5,
55% of births produce a new genotype, and about four in five of the
ancestor's direct mutants differ from it only in the `nop-C` filler.

What to read in `net.extras`:

- `analysis`: totals (`births`, `overwrites`, `parent_replaced`,
  `age_deaths`, `copy_mutations`, `tasks_performed`, `dominant`, the most
  common living genotype) and `per_update`, lists with one value per update:
  `organisms`, `genotypes`, `births`, `overwrites`, `age_deaths`,
  `ave_merit`, `ave_gestation`, `ave_generation`, and `tasks`, the number of
  organisms whose last gestation included each task.
- `final_state`: living organisms per genotype at the end.
- `genotype_parent`: the parent genotype of each genotype when it first
  appeared, from which lines of descent can be rebuilt.
- `task_events`: `{genotype: {task: count}}`.
- `space`: the torus and `final_grid`, the genotype in each cell.
- `instruction_set` (letter to instruction) and `environment` (the rewarded
  tasks and their multipliers).
- A species' `structure` is its genome as letters.

**Pure replication.** With the three mutation probabilities set to 0, only the
ancestor exists, and the network reduces to the four reaction kinds on one
species:

```python
net = chemart.generate_network("avida", seed=1, copy_mut_prob=0.0,
                               divide_ins_prob=0.0, divide_del_prob=0.0)
# 1 species, 4 reactions; their counts, writing P for the ancestor:
# P -> 2 P  147   P + P -> P + P  268   P -> P  44   P -> ∅  4
```

**Merit decides competition.** `generate_network` injects a single ancestor,
but the `World` class that runs it can take two. Put the ancestor and the NOT
performer in the same 12 × 12 world, without mutations:

```python
import numpy as np
from chemart.chemistries import avida as A
w = A.World(np.random.default_rng(0), width=12, height=12,
            copy_mut_prob=0.0, divide_ins_prob=0.0, divide_del_prob=0.0)
w.inject(A.DEFAULT_ANCESTOR, 0)
w.inject(NOT, 78)
w.run(400)
```

Over seeds 0–4, the ancestor holds 16 to 19 cells at update 100 against about
100 for the NOT performer, 0 to 6 at update 200, and none at update 400. Each
run takes about 4 seconds.

**Evolving tasks from scratch.** A 20 × 20 world run for 1,500 updates takes
about a minute (it executes about 18 million instructions):

```python
net = chemart.generate_network("avida", seed=1, world_x=20, world_y=20, updates=1500)
```

In three such runs (seeds 1–3), the ancestor's descendants evolved NOT and
NAND in all three, first appearing between updates 757 and 1,413, and ORN in
two. In seed 1, 112 of the 397 organisms performed NAND at the end, and the
average merit had risen from about 95 to 156. No run evolved the harder
tasks in that time.

**Other environments.** `rewards` maps task names to the exponent of the
bonus (`{"NOT": 1.0}` means ×2). `rewards={}` rewards nothing, so only
replication speed matters; `rewards={"EQU": 5.0}` is the "EQU only"
environment of Lenski et al. (2003).

**Other settings.** `birth_method="mass-action"` places offspring anywhere
in the world, like a well-stirred flask instead of a Petri dish. `age_limit=0`
turns off death of old age. `ancestor` takes any genome of 8 to 2,048
letters.

**Towards paper scale.** Lenski et al. (2003) used a 60 × 60 world,
`copy_mut_prob=0.0025` and 100,000 updates. The early growth of such a world
runs quickly: from seeds 1–3, `world_x=60, world_y=60, updates=200` gives 45,
39 and 65 organisms at update 100 and 348, 318 and 398 at update 200, in 2 to
3 seconds each. A full 60 × 60 world executes 108,000 instructions per update.
The 20 × 20 runs above managed about 250,000 instructions per second, so
100,000 updates of a full world (about 10^10 instructions) would take about
half a day per population in this pure-Python version.

## Results

### Evolution of complex features (Lenski et al. 2003)

Avida's best-known result answers an old objection to Darwin: how can a
complex feature, one that needs many parts working together, arise from
random mutation and selection? Lenski, Ofria, Pennock and Adami studied the
most complex of the nine tasks, EQU, which needs at least five `nand`
operations. Their shortest hand-written program for it was 19 instructions
long. They evolved populations of 3,600 organisms from a 50-instruction
ancestor that could copy itself but perform no task, for 100,000 updates
(15,873 ancestral generations), with 0.0025 copy mutations per instruction
and 0.05 insertions and deletions per copy.

- In the "reward-all" environment of the table above, **23 of 50 populations
  evolved EQU**. When only EQU was rewarded, **none of 50 did**
  (P < 4.3 × 10^-9, Fisher's exact test). Complex functions evolved by
  building on simpler ones, but only if the simpler ones were themselves
  rewarded.
- No particular simpler function was needed: in 36 environments that each
  left one or two simpler functions unrewarded, 124 of 360 populations (34%)
  still evolved EQU.
- The first EQU performer usually differed from its parent by a single
  mutation (19 of 23 populations; double mutations in the other four), yet
  knocking out instructions one at a time showed that 17 to 43 instructions
  (median 28) were needed to perform it. The feature was complex from the
  moment it appeared, because it was built from parts that had evolved for
  other functions.
- In the case-study population, the mutation just before EQU appeared had
  knocked out NAND and was strongly deleterious; reverting it removed EQU.
  In three of the 23 populations EQU depended on such a previously harmful
  mutation. The path to a complex feature can include backward steps.

Chemart reproduces the ingredients of this experiment but not the experiment.
Its tests check the rewards of table 1 for all nine tasks (each paid at most
once per gestation), that a NOT genome gets exactly twice the
merit, and that CPU time is proportional to merit (a 2:1 ratio of executed
instructions for a 2:1 ratio of merit). It does not run the 50-population
comparison, which needs roughly 10^10 instructions per population. The paper
used Avida 1.6 and an ancestor it does not print, so Chemart uses Avida 2's
100-instruction ancestor instead. In Avida 1.6, as the paper describes it,
the size factor of merit was the genome length; Chemart follows Avida 2's
default, the smallest of genome length, copied size and executed size.

### Other findings from the Avida group

These studies are summarised by Ofria et al. (2009) and the book. Chemart
reproduces none of them.

- **Survival of the flattest** (Wilke et al. 2001). Forty pairs of strains
  were adapted to a low and a high mutation rate, then competed. At a high
  enough mutation rate the strain adapted to it always won, although it
  replicated more slowly: its mutants lost less fitness. At high mutation
  rates selection acts on a cloud of mutants rather than on individuals,
  as quasispecies theory predicts ([quasispecies](quasispecies.md)).
- **Stable ecosystems** (Cooper and Ofria 2002). With nine depletable
  resources flowing in and 1% of unused resource flowing out, in populations
  of 2,500 organisms, 30 trials produced varied communities: nine
  coexisting specialists, a few generalists, or mixtures. The communities
  persisted after mutations were switched off. Chow et al. (2004) varied the
  resource inflow over six orders of magnitude and found multi-species
  communities only at intermediate abundance. The book (§8.2.3) adds that
  Avida populations diversify into host-parasite ecologies.
- **Evolvable languages.** Ofria, Adami and Collier (2002) studied which
  instruction sets evolve well; the book reports their finding that 99.7% of
  mutations in effective code are lethal or harmful, and that removing
  instruction arguments and protecting memory helps (§8.3.3). Bryson and Ofria
  (2013) compared a wider range of instruction sets.
- **Gradual complexity** (Ofria, Huang and Torng 2008): new complex traits
  appear mostly by combining and extending existing ones (book §8.4.4).
- The book also names studies of the growth of complexity (Adami, Ofria and
  Collier 2000), gene expression, and parallels with population dynamics in
  *E. coli*, and notes that the journal *Artificial Life* devoted an issue to
  Avida experiments in 2004.

Ofria et al. (2009) also record a warning: organisms exploit any loophole in
an experiment. With the original 16-bit inputs, guessing an answer at random
was often faster than computing it; inputs are now 32 bits.

### Evolving distributed algorithms (book §17.3.2)

With merit as a user-defined reward, Avida can also be used as an
evolutionary algorithm for engineering. The book lists distributed leader
election (Knoester, McKinley and Ofria 2007), energy saving in mobile devices
(McKinley et al. 2008), flood monitoring with sensor networks (Goings et al.)
and, in detail, synchronization and desynchronization (Knoester and McKinley
2011). In that work the instruction set gained a `flash` instruction and an
`if-recvd-flash` test, and the population was split into sub-populations
called *demes*. Whole demes were rewarded for flashing together (or as far
apart as possible), fitter demes replaced others, and each new deme was
seeded from the deme's original genome (*germline replication*). The evolved
programs adjust the timing of their flashes on receiving a neighbour's flash,
using only local information. Chemart implements neither demes nor the
extra instructions.

### What Chemart reproduces

Chemart is a port of the Avida 2.14 source for its default configuration,
checked against the compiled original.

- **The ancestor replicates exactly** without mutations: 389 instructions per
  gestation, 97 lines executed, merit 97, as in Avida's own
  `heads_default_100u` test. Tested, together with the instruction examples
  of Ofria et al. (2009, figures 2 and 3), eight genomes whose gestation,
  merit, executed size and tasks match the original's analyze mode, and the
  rules under which `h-divide` fails. The catalog's decisions report that
  626 dividing genomes out of 1,003 tested agree exactly.
- **Merit doubles for NOT** and follows table 1 for all nine tasks. Tested
  (see above).
- **The population fills the grid**, after which births kill organisms and
  space is the resource. The tests check that a default run contains births
  into empty cells and overwrites, that the population ends above 100
  organisms, and that the event counts in the network match the totals. A
  slow test checks growth on the 60 × 60 world against the original: upstream
  Avida averages 64.5 organisms at update 100 and 386 ± 29 at update 200, the
  port 64 and 400 ± 31 (from the decisions).
- **The evolution of EQU** (Lenski et al. 2003): not reproduced, for lack of
  speed; simple tasks do evolve in minute-long runs (see *Using it*).
- **Host-parasite ecologies** (§8.2.3): not reproduced. Chemart has no
  parasites, no depletable resources and one fixed reward per task.
- **Decentralised synchronization with demes** (§17.3.2): not reproduced;
  demes are not implemented.

The random numbers come from numpy, so runs match the original statistically
but not bit for bit. Point mutations, sex, other schedulers and the 77-task
environment are also left out; the implementation decisions list them.

## Further reading

- Adami, C., Ofria, C. & Collier, T. C. (2000). Evolution of biological
  complexity. *PNAS* 97, 4463–4468. (Book ref. [21].)
- Wilke, C. O., Wang, J. L., Ofria, C., Lenski, R. E. & Adami, C. (2001).
  Evolution of digital organisms at high mutation rates leads to survival of
  the flattest. *Nature* 412, 331–333.
- Cooper, T. F. & Ofria, C. (2002). Evolution of stable ecosystems in
  populations of digital organisms. *Artificial Life VIII*, 227–232. (Book
  ref. [198].)
- Chow, S. S., Wilke, C. O., Ofria, C., Lenski, R. E. & Adami, C. (2004).
  Adaptive radiation from resource competition in digital organisms.
  *Science* 305, 84–86.
- Ofria, C., Adami, C. & Collier, T. C. (2002). Design of evolvable computer
  languages. *IEEE Transactions on Evolutionary Computation* 6(4), 420–424.
  (Book ref. [628].)
- Ofria, C., Huang, W. & Torng, E. (2008). On the gradual evolution of
  complexity and the sudden emergence of complex features. *Artificial Life*
  14(3), 255–263. (Book ref. [630].)
- Bryson, D. M. & Ofria, C. (2013). Understanding evolutionary potential in
  virtual CPU instruction set architectures. *PLoS ONE* 8, e83242. (Book
  ref. [139].)
- Goings, S., Goldsby, H., Cheng, B. H. & Ofria, C. An ecology-based
  evolutionary algorithm to evolve solutions to complex problems. In
  *Artificial Life 13*, 171–177. (Book ref. [337].)
- The Avida source code: <https://github.com/devosoft/avida>
