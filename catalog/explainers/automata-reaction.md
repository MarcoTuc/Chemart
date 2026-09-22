## Introduction

The automata reaction is an artificial chemistry built by Peter Dittrich and
Wolfgang Banzhaf at the University of Dortmund and published in *Artificial
Life* in 1998 as "Self-evolution in a constructive binary string system". Its
molecules are 32-bit binary words. When two of them collide, the first is read
as a tiny program, the program runs on the second as its input, and whatever it
writes is a new, third word. Both colliding words survive. Every word is
therefore two things at once: a machine that can process other words, and data
that other words can process.

The question the authors asked was whether evolution can arise in such a soup
on its own. Evolutionary algorithms and most artificial-life systems supply
evolution from outside: a mutation operator changes individuals, a fitness
function scores them, a selection step picks the winners. Dittrich and Banzhaf
wanted a system with none of these, in which, in their words, "every variation
is carried out by the individuals" and "selection pressure is generated
implicitly through interaction among the individuals". They called what they
observed *self-evolution*. Along the way the soup shows a range of simpler
behaviours, depending on its size: takeover by a single word, a small closed
group of words that keep regenerating each other, and a long exploratory phase
ended by the arrival of self-copying words.

The design choice that sets it apart is speed. Banzhaf and Yamamoto place it at
the start of their section on "Artificial Chemistries Based on von Neumann
Machines" (book §10.6.1): the reaction was designed "respecting the medium",
here the ordinary computer, so 32-bit words fit a processor register and
populations of about a million molecules can be run for thousands of
generations. The paper says the machine was inspired by Hofstadter's
[Typogenetics](typogenetics.md), where strings also act on strings. It extends
Banzhaf's earlier binary-string chemistry, in which strings are folded into
matrices instead of programs ([matrix chemistry](matrix-chemistry.md)). Its
findings are compared with Fontana's lambda-calculus chemistry
[AlChemy](alchemy.md), whose closed, self-maintaining organizations it also
produces. Unlike Ikegami and Hashimoto's [machine-tape chemistry](ikegami-hashimoto.md)
there is no separate population of machines: any word can play either role.
And unlike the other von Neumann chemistries that follow it in the book,
[Coreworld](coreworld.md) and [Tierra](tierra.md), programs do not run in a
shared memory and copy themselves; each collision is one short, bounded
computation that yields one new word.

It is a simulation model: a set of molecules, a reaction rule and a reactor, to
be run and watched.

## How it works

### Molecules are words, collisions are short programs

A molecule is one of the 2³² possible 32-bit words. Chemart names it by its
hexadecimal value with a `w` in front: `w1e1ca260` is the word `0x1e1ca260`.

A collision takes an ordered pair of words, the *operator* s1 and the
*operand* s2, and computes a product s3:

```
s1 + s2  ->  s1 + s2 + s3
```

The computation happens on a small machine with two 32-bit registers. The
*operator register* holds s1 and can only be read. The *IO register* starts as
a copy of s2, can be written, and holds s3 when the machine stops. Each register
has a pointer to one bit, and both pointers start at bit 0, the least
significant bit (the right-hand end of the word as written).

The program comes from s1 itself. The 32 bits are cut into eight 4-bit
groups (*nibbles*), read from the least significant one, and each nibble is
looked up in a table of 16 instructions:

- **logic instructions** (`ID`, `NOT`, `AND`, `OR`, `EXOR`, `EQ`) combine the
  bit under the operator pointer with the bit under the IO pointer and write
  the result into the IO register; `ID` simply copies the operator bit. Then
  the pointers move one step.
- **`MOV`** moves the pointers. If a 4-bit *pattern* has been set with
  `SETP <pattern>` (which uses the next nibble as its argument), `MOV` keeps
  moving, up to 32 steps, until the pattern is found in the register.
- **copy mode** (`CPON`, `CPOFF`): while it is on, every pointer step first
  repeats the last logic operation, so one `MOV` can rewrite a whole run of
  bits.
- `TDIR` reverses the direction of movement, `TMM` switches between moving
  both pointers, only the operator pointer or only the IO pointer, `UNSETP`
  clears the pattern, `NOP` does nothing, and `STOP` halts the machine.

There are no jumps and no loops. The program ends after its eighth nibble or at
`STOP`, so every collision finishes in a bounded number of steps and always
yields a product. There are two *code tables*: they agree on 15 codes and differ
only at `1010`, which is `NOT` in table 1 and `EQ` (bit equality) in table 2.

### Three collisions

The first reaction of the default run (shown below) is the simplest kind of
change:

```
wb7362806 + w13410911 -> wb7362806 + w13410911 + w13410910
```

The operator `b7362806` begins, from the right, with the nibbles `6`, `0`, `8`:
`CPON`, `ID`, `STOP`. `ID` copies bit 0 of the operator (a 0) over bit 0 of the
operand (a 1), and `STOP` ends the program. The product is the operand with its
last bit changed: `13410911` becomes `13410910`.

A word whose first nibble is `8` starts with `STOP`, so it returns every operand
unchanged. `9a716b98` is such a word, and it is one of the paper's examples of
a *passive replicator*, a word s1 for which `s1 + s2 => s2`:

```
9a716b98 + 9878fa02 => 9878fa02
```

The opposite is an *active replicator*, which writes a copy of itself whatever
the operand: `s1 + s2 => s1`. The paper's example `1e1ca260` folds into
`ID, CPON, SETP 1010, EXOR, MOV, ID, MOV` and, as Chemart confirms,

```
1e1ca260 + 078f21ff => 1e1ca260
```

Because the pointers always start at the least significant bit and move toward
the higher bits unless the program says otherwise, most changes happen at the
right-hand end of the operand. Products therefore usually look like their
parents, which the paper names as "a prerequisite for evolution".

### The reactor

The population, the *soup*, is a multiset of M words, initially drawn at
random. The reactor algorithm repeats two steps: pick two words at random
without removing them, compute their product, and put the product in place of a
randomly chosen word. The soup size stays M, and the words pushed out form the
*dilution flux*. M collisions make one *generation*, so runs of different sizes
can be compared.

A collision can also be *elastic*: nothing is inserted. The paper uses this
through a *filter*, a condition the product must pass. Its filter f1
(paper eq. 3) rejects any product equal to one of the two reactants, which
switches off exact replication, both passive and active. With the filter,
words can no longer multiply by copying; they persist only if other words keep
producing them.

For comparison the paper also runs a trivial *AND reaction*, in which the
product is the bitwise AND of the two words.

### What is measured

The paper follows the soup with population-level ("macroscopic") measures, and
Chemart reports three of them for every generation:

- **diversity**: the number of distinct words in the soup divided by M; 1 means
  every word is different, 1/M means a single word fills the soup;
- **productivity**: the fraction of collisions that inserted a product (always
  1 without a filter);
- **innovativity**: the fraction of collisions that produced a word never seen
  before in the run.

Two further notions describe groups of words. A set of words is *closed* if
every collision among them produces a word of the set, and *self-maintaining*
if every word of the set is produced by some collision within it. A closed,
self-maintaining set is an *organization*, the notion taken from Fontana's
work. The paper shows them as *reaction tables*: row i, column j gives the
product of word i acting on word j.

In the formal specification below, `A_s1(s2)` is the product of the machine
folded from s1 run on s2, and "elastic if the filter rejects s3" is filter f1.
The *catalytic network equation* named in the decisions is the paper's
differential-equation description of the same reactor for large M: each
reaction `i + j -> i + j + k` has rate constant 1 if it exists and passes the
filter, 0 otherwise.

## Using it

The automata reaction has two faces. `chemart.generate_network("automata-reaction")`,
printed above, returns a *closure*: every reaction reachable from 10 random
words, cut off at 50 species (`status=truncated`); it is mainly a tool to
check published organizations (below). `chemart.evolve("automata-reaction")`
runs the paper's reactor and returns a trajectory: a frame per generation of M
collisions, with the contents of the soup, and at the end the network of every
distinct reaction that fired.

The default run is the reactor on 1,000 random words for 10 generations with
code table 1 and no filter, 10,000 collisions in about a second. It is not a
published experiment, but it already shows the trend of the paper's larger
runs:

```python
traj = chemart.evolve("automata-reaction", seed=1)
[round(len(f.state) / 1000, 3) for f in traj.frames]    # diversity
# [1.0, 0.827, 0.734, 0.725, 0.706, 0.691, 0.678, 0.665, 0.648, 0.623, 0.604]
traj.series("productivity")   # [None, 1.0, 1.0, ...]   every collision yields a product
traj.series("innovativity")
# [None, 0.599, 0.484, 0.472, 0.424, 0.373, 0.341, 0.303, 0.299, 0.278, 0.245]
net = traj.network
print(net.summary())
# automata-reaction: 4818 species, 9872 reactions, status=observed
list(net.extras["final_state"].items())[:3]
# [('w94c116e2', 100), ('w94c116e0', 13), ('w94c116e3', 11)]
```

Frame `t` counts generations. Each frame's `state` is the soup at that moment,
so diversity is the number of words in it divided by M. Productivity and
innovativity are the frame's `observables`; they describe the generation that
ended at the frame, so the first frame, the initial soup, has none.
`final_state` is the soup at the end, most frequent word first. The top word
here, `94c116e2`, is an active replicator: it returned itself on 100 of 100
random operands we tried. Each distinct reaction that happened is recorded
once in `net.reactions`, with how often it fired and a rate constant `k` that
counts how many of the two orders of the pair give that product:

```
wb7362806 + w13410911 -> wb7362806 + w13410911 + w13410910  [mass-action k=1.0]  (x1)
```

The species' `structure` field is the word as a 32-character bit string, most
significant bit first.

**Extinction by the AND reaction** (paper Fig. 2). The all-zero word is produced
by many pairs and copies itself with any partner, so it takes over:

```python
traj = chemart.evolve("automata-reaction", seed=0, mechanism="and", M=1000, generations=10)
[round(len(f.state) / 1000, 3) for f in traj.frames]
# [1.0, 0.835, 0.467, 0.212, 0.092, 0.038, 0.02, 0.011, 0.005, 0.003, 0.001]
list(traj.network.extras["final_state"].items())   # [('w00000000', 1000)]
```

**A small soup settles into an organization** (paper Fig. 3, M = 100):

```python
traj = chemart.evolve("automata-reaction", seed=0, M=100, generations=150)
traj.network.extras["final_state"]
# {'w66847fcb': 34, 'w66847fc9': 33, 'w66847fca': 19, 'w66847fc8': 14}
```

Diversity falls from 1.0 to 0.14 by generation 50 and 0.04 by 100, and no new
word appears in the last ten generations. Four words that differ only in their
last two bits remain, like the paper's four (`7240a7ef`, `7240a7ea`,
`7240a7eb`, `7240a7ee`). To check that such a set is closed, use
`generate_network`, which computes every reaction reachable from the given
words instead of running a soup:

```python
net = chemart.generate_network("automata-reaction",
                               words=["66847fcb", "66847fc9", "66847fca", "66847fc8"])
print(net.summary())
# automata-reaction: 4 species, 12 reactions, status=complete
```

`status=complete` means no collision leads outside the four words. On random
words (`n_seeds` of them) the closure grows until `max_species` stops it
(`status=truncated`).

**Exploration ended by active replicators** (paper Fig. 4, M = 10⁴). This run
takes about 45 seconds:

```python
traj = chemart.evolve("automata-reaction", seed=1, M=10000, generations=60)
[round(len(f.state) / 10000, 3) for f in traj.frames][::5]     # generations 0, 5, ..., 60
# [1.0, 0.697, 0.687, 0.669, 0.597, 0.283, 0.032, 0.007, 0.005, 0.004, 0.004, 0.003, 0.003]
traj.series("innovativity")[5::5]                              # generations 5, 10, ..., 60
# [0.3876, 0.3046, 0.2583, 0.199, 0.0492, 0.0008, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```

For about 20 generations the soup stays diverse and keeps producing new words;
then diversity collapses between generations 20 and 30. The two most frequent
words at the end, `616614e2` (4,684 copies) and `c16614e2`, are active
replicators.

**The evolution setting** (paper §5.4) is `code_table=2` and
`forbid_exact_replication=True`. At M = 10⁴ for 60 generations (about 55 s) the
soup is still in its exploratory phase: diversity 0.81, productivity about 0.64,
innovativity about 0.18 at the end, and 152,730 species in the
observed network. The paper's runs used M = 10⁵ or 10⁶ for 1,000 to 7,000
generations, that is 10⁸ collisions or more. Chemart's pure-Python machine
manages roughly 10,000 collisions a second, so these would take hours; the
YAML notes that paper-scale runs need a vectorised or compiled machine.

**Single reactions.** The machine itself can be called directly:

```python
from chemart.chemistries.automata_reaction import automata, disassemble
disassemble(0x101d662b)            # ['EXOR', 'SETP 0110', 'CPON', 'NOP', 'MOV', 'ID', 'MOV']
f"{automata(0x101d662b, 0xcb5d5c2d):08x}"   # 'cb5d5606'
```

This is the paper's worked example (its Table 1); `disassemble(w, 2)` uses code
table 2. Pass `words=[...]` to start a soup from chosen words instead of `M`
random ones.

## Results

**Static properties of the reaction.** The paper's Table 2 lists eight
reactions: typical ones, one using `TDIR`, two passive and two active
replications. It observes that products resemble their parents and that
changes appear mostly at the right edge of the operand, "sometimes on the
left", because the pointers start at bit 0 moving left. It reports that about
30% of random strings are passive replicators and about 0.004% active
self-replicators, without saying how many operands each string was tested on.
Chemart reproduces all of Table 2 and the program listings of Figs. 3, 4 and 10
exactly. With its machine, 37% of 20,000 random strings copied one random
operand unchanged and 20% copied each of three; the decisions record about
0.005% active replicators for eight operands. The test checks only that the
passive fraction lies between 15% and 45% and the active fraction below 0.2%.
The resemblance of products to parents is not tested statistically.

**Four kinds of soup behaviour.** Section 5 of the paper shows "typical" runs,
since averaging many runs was not feasible then:

- *Extinction* (Fig. 2, AND reaction, M = 10⁴): diversity drops to 1/M as
  `00000000`, which "is able to replicate with every other string", fills the
  soup. Chemart's test reproduces this at M = 1,000: at least 990 of 1,000
  words are `00000000` after 10 generations.
- *Emergence of an organization* (Fig. 3, M = 100, code table 1): diversity and
  innovativity drop until a closed, self-maintaining set of four strings
  dominates; drift then slowly thins it. The tests check that the published
  four-string set is closed and gives exactly the published reaction table,
  and that soups with M = 100 run for 150 generations end with
  diversity at most 0.1, no new words in the last ten generations, and the two
  most frequent survivors differing in at most 8 bits.
- *Exploration and innovation* (Fig. 4, M = 10⁴): after a quick initial drop,
  an explorative phase (generations 2 to 40) of high diversity, innovativity
  and productivity ends "around gen. 40" when active replicators appear and
  take over. The tests check that the paper's four replicators (`1e1ca260` and
  relatives) form the published reaction table and copy themselves on random
  operands; the soup run itself is not tested, but the run under *Using it*
  shows the same shape, with the collapse near generation 25.
- *Evolution* (Figs. 5 to 9, M = 10⁵ and 10⁶, code table 2, filter f1): the
  headline result. Code table 2 was chosen because without `NOT` "it is harder
  to construct a self-maintaining organization" under the filter. An
  exploratory phase (generations 0 to 110) gives way to a period (110 to 200)
  in which "many different species are competing for space" and are replaced
  by "better" ones (Fig. 7 follows representative strings at M = 10⁶). The
  resulting organization contains "approximately 50000 different string types
  where only about 3000 are present at a specific time", and it keeps
  developing over 7,000 generations by absorbing newly invented strings. The
  authors liken the long quasi-stable periods interrupted by rapid change to
  punctuated equilibrium. Reaction tables of the 90 most frequent strings
  (Figs. 8, 9) show the strings becoming closely coupled and very similar, with
  fewer elastic collisions. Chemart does not reproduce these runs: they are far
  beyond the speed of its machine (see *Using it*). What is tested is the
  static part: the published eight-string block of the generation-3,400 table
  (Fig. 9) is recomputed exactly, and the filter is checked to make every
  replication elastic.

**Cross-over without a cross-over operator** (Fig. 10). In the evolution runs
the authors found strings that insert part of their own sequence into any
operand, for example `1e64a24e + 00000000 => 0000024e`, so that eight such
strings form a closed organization of recombinants. Recombination, normally an
operator imposed by an evolutionary algorithm, appears as a reaction. Chemart's
tests reproduce this reaction and two others from the figure exactly.

**The authors' reading.** The paper concludes that "complex forms of evolution
can take place in systems without any explicit variation operator, like
mutation or recombination, and without any explicit (artificial) selection".
Compared with Bagley, Farmer and Fontana's model of evolving metabolisms, it
needs only rate constants 0 and 1, only catalysed reactions and a much simpler
simulation algorithm, and it can hold far more distinct objects (10⁶ against
75); its drawbacks are a finite soup, the difficulty of very different rates,
and fixed-length strings. Banzhaf and Yamamoto summarise the observations as
dominance of self-replicators, formation of catalytic networks, increasing
complexity when exact replication is forbidden, and syntactic similarity within
emergent organizations, which "support the findings of Fontana et al.,
especially the phenomenon of syntactic and semantic closure". Chemart tests the
closure of the published organizations; it does not test semantic closure in
Fontana's broader sense.

**Lazy and eager replicators.** The book also lists "the spontaneous transition
from lazy to eager replicators", where a lazy replicator copies whatever it
meets (`s1 + s2 -> s2`, the passive replicator above) and an eager one copies
itself (`s1 + s2 -> s1`). The book gives no separate source for it; the 1998
paper describes the soup evolving "towards active replication" (Fig. 4).
Chemart's M = 10⁴ run shows active replicators taking over, but no test checks
it.

**Later work.** Dittrich, Ziegler and Banzhaf (1998) used the same chemistry,
with code table 2 and M = 10⁵, to develop a *mesoscopic* analysis: classifying
and clustering the population between single-molecule tracing and whole-soup
averages. In a series of 100 runs they found behaviour "ranging from early
stabilization with short transients to complex, oscillating dynamics". In run
A4-23 they traced a sudden dip in productivity around generation 700 to a new
cluster of strings that had acquired the ability to replicate the left side of
a string, not only the right. That paper states its filter as
`s1 ≠ s2 and s1 ≠ s3`, which differs from eq. 3 of the 1998 journal paper
that Chemart implements. Banzhaf, Dittrich and Eller (1999) placed binary
strings on a spatial topology; Chemart does not include that variant (see the
decisions).

**Not reproduced.** The distance distribution complexity (DDC), the paper's
fourth macroscopic measure (the entropy of the distribution of Hamming
distances between strings), is not computed. The C source's code tables 3 and
4 and its binding-preference option, unused in the paper, are not exposed.

## Further reading

- Dittrich, P., Ziegler, J. & Banzhaf, W. (1998). Mesoscopic analysis of
  self-evolution in an artificial chemistry. In C. Adami, R. K. Belew,
  H. Kitano & C. E. Taylor (eds.), *Artificial Life VI*, pp. 95–103. MIT Press.
  <https://users.fmi.uni-jena.de/~dittrich/p/DZB98alife6.ps>
- Banzhaf, W., Dittrich, P. & Eller, B. (1999). Topological interactions in a
  binary string system. *Physica D* 125, 85–104. The spatial variant.
- Bagley, R. J., Farmer, J. D. & Fontana, W. (1992). Evolution of a metabolism.
  In C. G. Langton, C. Taylor, J. D. Farmer & S. Rasmussen (eds.),
  *Proceedings of the Workshop on Artificial Life (ALIFE '90)*, pp. 141–158.
  Addison-Wesley. The model the paper compares itself with.
- Peter Dittrich's publication list, with preprints and the autoreac source:
  <https://users.fmi.uni-jena.de/~dittrich/>
