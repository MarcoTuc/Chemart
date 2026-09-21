## Introduction

Laing's artificial molecular machines are a thought experiment from the 1970s:
a chemistry in which molecules are small computing machines. Richard Laing
described it in a series of papers between 1972 and 1977 on what he called
*artificial organisms*. Banzhaf and Yamamoto say he argued for a general theory
of living systems built on such organisms, and that the idea foreshadows the
Artificial Life theme of "life as it could be". Freitas and Merkle quote Laing
on his aim: machines "all of whose basic constituents and operations are the
sort that reasonably may be said to be possible for biological systems at the
macromolecular level". What he asked of such machines was that they compute,
and above all that they reproduce.

The molecules are strings, and they come in two forms. A **passive** string is
a tape of constituents, each in state 0 or 1, like the tape of a Turing
machine. An **active** string is a machine: a chain of instruction
constituents, folded in three dimensions so that some distant instructions
touch. The fold is what gives the machine jumps and loops, in the way Wang's
version of the Turing machine allows jumps between instructions. When a
machine meets a tape they touch at one point. The instruction at that point
acts on the tape constituent it touches (it writes a 0 or a 1, reads it, or
slides to a neighbour), and then the next instruction comes into contact. The
book's figure 10.5, reprinted from Laing (1975), draws a folded ribbon of
instructions (`w(1)`, `w(0)`, `TT`, `CT`, `L`, `R`, `NOP`) wrapped around a tape
of 0s and 1s.

Laing used this system for two claims, both made on paper, by construction. First, one active
string can play any Turing machine's read-head and a passive string its tape,
so the system can carry out any computation (Laing 1973). Second, a machine can
reproduce by **self-inspection** (Laing 1975, 1976, 1977): it reads its own
structure, writes a description of itself onto a tape, and then builds a copy
from that description. Sipper's survey of self-replication singles this out:
in Laing's model the description "is dynamically constructed concomitantly with
its interpretation", whereas in most other models the genome is given in
advance, by design or by evolution. Freitas and Merkle, who call the system
"Laing molecular tapeworms", place it among the "hybrid cellular-kinematic"
replicators: strings of automata that move and slide against each other,
rather than cells fixed in a grid as in von Neumann's cellular automaton.

It is a **formalism**, not a simulation model. The book notes that Laing did
not consider an explicit dynamics, probably because simulating it was beyond
the computers of the time. It opens the book's section on artificial
chemistries based on Turing machines (§10.5), whose other entries have all
been simulated: [Typogenetics](typogenetics.md), where strands are translated
into enzymes that act on strands; McCaskill's [polymers as Turing
machines](mccaskill-polymer-tm.md), where strings recognise each other by
pattern matching and one processes the other; and Ikegami and Hashimoto's
[machine-tape chemistry](ikegami-hashimoto.md), where a machine reading a tape
makes a new machine and a new tape. In Laing's system, as implemented here,
machines and tapes stay separate kinds: a machine only ever rewrites, extends
or cuts tapes, and no tape ever becomes a machine.

## How it works

### Tapes and machines

In Chemart a tape is written `t:` followed by its bits, and a machine `m:`
followed by its instructions separated by dots. The instructions come from the
book's figure 10.5 and from Freitas and Merkle's list of the constituents
of one of Laing's versions, built from rigid constituents:

| instruction | what it does to the tape constituent in contact |
|---|---|
| `W0`, `W1` | put it in state 0 or 1 |
| `L`, `R` | slide the contact one constituent left or right; sliding off an end recruits a new 0 constituent there |
| `CT<k>` | *conditional transfer*: if the constituent is 1, continue at the matching `TT<k>`; if it is 0, go on to the next instruction |
| `TT<k>` | *transfer target*, the place a `CT<k>` jumps to; does nothing when reached in sequence |
| `NOP` | nothing |
| `D` | *detach*: cut the tape just left of the contacted constituent, releasing the left part as a separate tape |
| `H` | halt |

The fold that brings a `CT` into contact with its `TT` is written as a shared
label: `CT1` jumps to `TT1`. Every label has exactly one target, and several
`CT`s may jump to it. Running past the last instruction also halts the machine.
Chemart describes each machine's shape in its species `structure`, listing each
instruction with its position and each `CT` with the position it jumps to.

### A reaction

A machine and a tape react; two machines or two tapes do not. The machine
attaches with its first instruction to one unit of the tape (the leftmost, by
default), runs until it halts, and lets go. The machine is unchanged, so it
acts as a catalyst. The tape is replaced by whatever the run left: one changed
tape, or several pieces if `D` cut it. In general,

```
m + t -> m + t_1 + ... + t_n
```

A run that leaves the tape as it was is not a reaction and is not listed.

Here is a reaction from the default network, whose machine adds 1 to a binary
number stored with its lowest bit on the left (so `t:11` is 3 and `t:001` is 4):

```
m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:11 -> m:CT1.W1.H.TT1.W0.R.CT1.W1 + t:001
```

Step by step, with the machine attached to the leftmost bit:

1. `CT1` reads a 1, so it jumps to `TT1`, and the machine continues after it.
2. `W0` sets the bit to 0 (1 plus 1 is 0, carry 1). `R` slides to the second bit.
3. `CT1` reads 1 again and jumps back to `TT1`; `W0` sets it to 0; `R` slides
   off the right end of the tape, which recruits a new 0 there.
4. `CT1` reads that 0 and does not jump, so `W1` writes the carry. The program
   has run out, so the machine halts. The tape now reads `001`.

On `t:0` the first `CT1` does not jump, `W1` writes a 1 and `H` stops the
machine after three instructions: `t:0` becomes `t:1`.

The two other ways a reaction can change a tape are shown by the machine
`R.D`. On `t:10` it slides to the 0 and detaches it, so one reaction gives two
products. On `t:0` it first slides off the end, recruiting a new 0, and then
cuts it off again:

```
m:R.D + t:10 -> m:R.D + t:0 + t:1
m:R.D + t:0 -> m:R.D + 2 t:0
```

### Limits

Laing's machines need not halt, and tapes can grow without bound. Chemart
therefore stops a run after `max_steps` instructions, and drops any tape longer
than `max_length`. Such a run gives no reaction, and the network's status
becomes *truncated*, meaning that reactions exist beyond what is listed.

### The reactor

Laing gave no dynamics, so Chemart offers two ways to turn the rules into a
network. The **closure** starts from the seed machines and tapes, applies every
machine to every tape, adds the products, and repeats until nothing new
appears or a limit is hit. It lists every reaction reachable, with no rates.
The **soup** is a Chemart addition: a population of molecules from which random
pairs are drawn; when a machine meets a tape they react, and whenever a
reaction increases the number of molecules, random molecules (machines or
tapes) are removed until the population is back to its starting size. Machines
are catalysts and, without `D`, a tape is replaced by a single tape, so with
the default incrementer nothing is ever removed.

## Using it

The default run is not a published experiment: it is a Chemart example. It
takes the incrementer machine above and the tape `t:0`, and computes the
closure. The result is a binary counter. Each reaction adds 1 to a tape, so the
closure holds 16 tapes, one for each number from `t:0` (0) to `t:1111` (15),
linked in counting order by 15 reactions. The next step would turn `t:1111` into
`t:00001`, five bits, longer than the default `max_length` of 4, so the status
is `truncated`. The 17th species is the machine itself. In the closure,
`net.extras["seed"]` only records the seed molecules.

Raising `max_length` extends the counter: `max_length=10` gives 1,025 species
in about 6 seconds and `max_length=11` gives 2,049 in about 13. Each extra bit
doubles the species and slightly more than doubles the time, so values much
above 12 are slow.

**Running a Turing machine.** `turing_program` compiles a two-symbol Turing
machine, given as a table from (state, symbol) to (symbol to write, move, next
state), into a Laing machine. This is Chemart's own construction, not Laing's.
Because `CT` only jumps on a 1, a jump to a state goes through two targets, one
for each symbol. The two-state *busy beaver*, the machine that writes the most
1s before halting, becomes a 27-instruction Laing machine:

```python
from chemart.chemistries.laing_molecular_machines import run, turing_program
BB2 = {("A", 0): (1, "R", "B"), ("A", 1): (1, "L", "B"),
       ("B", 0): (1, "L", "A"), ("B", 1): (1, "R", "H")}
prog = turing_program(BB2)
r = run(prog, "0", max_steps=10_000)
r.pieces, r.outcome, r.steps, r.counts["L"] + r.counts["R"]
# (['1111'], 'halt', 31, 6)
```

It writes four 1s in 31 instructions, six of which are slides, one per Turing
step. To see it as a network, pass the program as a machine:
`generate_network("laing-molecular-machines", machines=[".".join(prog)])` gives
the single reaction `m:CTA1.W1.CTA0...R.H + t:0 -> m:... + t:1111`.

**Cutting tapes.** `machines=["R.D"], tapes=["10"], max_length=2` gives the
three-reaction network of the detach example above, with status `complete`.

**Binding anywhere.** `binding="all"` lets the machine attach to every unit of
the tape, and lists each outcome as a separate reaction. With `machines=["W1"],
tapes=["000"], max_length=3` the tape `t:000` has three reactions, to `t:100`,
`t:010` and `t:001`, and the closure is all eight 3-bit tapes with 12
reactions. `binding="random"` draws one attachment point per pair instead.

**A soup.** `method="soup"` draws `steps` random pairs from `copies` copies of
each seed molecule:

```python
net = chemart.generate_network("laing-molecular-machines", seed=2, method="soup",
                               copies=20, steps=400, max_length=8)
net.summary()                # 18 species, 16 reactions, status=observed
net.extras["final_state"]
# {'m:CT1.W1.H.TT1.W0.R.CT1.W1': 20, 't:011': 4, 't:0001': 3, 't:0101': 3,
#  't:111': 2, 't:1001': 2, 't:1101': 2, 't:0011': 2, 't:00001': 1, 't:1111': 1}
```

Each reaction carries a `count` of how often it fired. Here each of the 20
tapes has been counted up from 0 to somewhere between 6 and 16, depending on
how often it happened to meet a machine.

## Results

Laing's results are mathematical constructions, and his papers were not
accessible for this entry (see the implementation decisions). What follows is
what the book and the two secondary sources, Sipper (1998) and Freitas and
Merkle (2004), report.

**Universal computation.** Laing (1973) proved that his artificial organisms
can perform any computation. As quoted by Sipper, Laing (1977) puts it this
way: "one string (an active string) can be designed to play the part of any
Turing machine finite-state read-head, and another string (a passive string)
can be designed to play the part of a Turing machine tape". Chemart reproduces the claim by example rather than by Laing's
proof: its tests compile three busy beavers with `turing_program` and check
their known results. The two-state machine writes 4 ones in 6 moves (Rado
1962); of the three-state machines, one writes 6 ones in 14 moves and another
takes 21 moves and writes 5 ones (Lin and Rado 1965). Other tests check the
incrementer on the numbers 1 to 39 and that the default closure is exactly the
16-tape counter.

**The instruction semantics.** Freitas and Merkle quote the rigid-constituent
version: a slide off the end of the tape means a new 0 "must be
synthesized or otherwise be recruited from the environment", and the detachment
molecule causes "the severing of the structural attachment that molecule has
with its predecessors". Chemart's tests check both: sliding off either end adds
a 0, and a detach cuts `1111` into `11` and `01`, so that one reaction yields
several product molecules. They also check how `CT` reads the tape and that a
non-halting machine is stopped by `max_steps`.

**Self-reproduction by self-inspection.** This is the subject of Laing's
papers from 1975 to 1977, and Chemart does not reproduce it. Freitas and Merkle quote Laing's procedure
for the rigid-constituent version: the machine first "produces a string of molecules which is a
description of itself"; reading that description, it determines each required
constituent and synthesises it at the end of the string; it ends with its
description followed by a copy of itself, which the detachment molecule then
cuts free; finally an "active" status is passed to the offspring, which folds
into shape by self-assembly. Laing (1977), in the passage Sipper quotes, valued
the capacity "to explore its own structure and produce a complete description
of it". Doing this needs instructions that turn a passive constituent into an
instruction and make it active: the "synthesize" constituent of the
rigid-constituent version, and the activation, activate-and-detach and conversion primitives of his 1977 version,
whose tape also has a third, null state. No accessible source says which state
or type these produce, so Chemart leaves them out. As a result machines are
never built, only tapes, and nothing in Chemart reproduces. Laing's
compartmentalised artificial organisms are not modelled either.

**Standing.** Sipper's fifty-year overview of self-replication treats Laing's
machines as one of its string-based models and notes that Laing later took
part in NASA's 1980 study, which examined a "seed" factory on the moon that
would replicate itself from local material. Freitas and
Merkle report that Laing thought an electronic engineer could build components
with the required properties. The book presents the work as a sketch: the
steps toward life's self-replication "were sketched", and there is no dynamics
to compare with. It does not describe later work that built directly on
Laing's machines.

## Further reading

- Rado, T. (1962). On non-computable functions. *Bell System Technical
  Journal* 41(3), 877–884. The busy beaver problem and the two-state result.
- Lin, S. & Rado, T. (1965). Computer studies of Turing machine problems.
  *Journal of the ACM* 12(2), 196–212. The three-state busy beavers.
- Laing, R. (1975). Some alternative reproductive strategies in artificial
  molecular machines. *Journal of Theoretical Biology* 54, 63–84.
  <https://doi.org/10.1016/S0022-5193(75)80055-5>
