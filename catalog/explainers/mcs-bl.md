## Introduction

The Molecular Classifier System (MCS.bl, "bl" for Broadcast Language) is a
string chemistry in which every molecule is a tiny if-then rule. When two
strings meet, one acts as the enzyme and the other as the substrate. If the
substrate matches the pattern in the enzyme's *condition*, the enzyme writes a
new string according to its *action*, and both reactants survive. A string
that reads "anything ending in 0 → the same thing ending in 1" is an enzyme
that flips the last symbol of its substrate; applied to a string that is itself
such a rule, it turns one enzyme into another. Because every string is at once
a possible rule and a possible substrate, the chemistry is *reflexive*: the
molecules rewrite each other.

It was built by James Decraene, George Mitchell and Barry McMullin at Dublin
City University between 2006 and 2011, within the European ESIGNET project
("Evolving Cell Signalling Networks in silico"). Their subject was the
networks of signalling molecules by which a cell processes information, and a
specific question about them: in a primitive cell with no genes, a network can
only be inherited if it keeps regenerating its own members, that is, if it is
*collectively autocatalytic* (a concept from Kauffman) or *self-maintaining*. Can such
networks evolve to do a useful job while never losing that property? The
strings come from John Holland's *broadcast language*, a rule notation from
his 1975 book *Adaptation in Natural and Artificial Systems*, related to his
later classifier systems. Holland's classifier systems keep rules and messages
apart; the MCS removes that separation, since in a cell every molecule can be
both (Decraene's thesis, 2009).

MCS.bl is a simulation model with two levels. At the molecular level, strings
collide at random in a reactor of fixed capacity and are copied with occasional
errors. At the cellular level, many such reactors ("cells") grow, divide and
replace one another, so that cells whose networks work better take over. The
work produced two kinds of result. First, a set of negative, "unexpected"
findings in a single reactor: a universal copier string gains nothing from
copying everything, and a well-designed replicator is destroyed by its own
mutants, which grow ever longer (the *elongation catastrophe*). Second, with
cells and selection aimed at one chosen molecule, a hand-designed
four-molecule network evolved into a slightly faster one, and two networks
merged into larger ones.

Banzhaf and Yamamoto present it first among the string-based chemistries of
their chapter on bio-inspired chemistries (book §11.1.1), before
[Stringmol](stringmol.md) and [SAC](sac.md). Like those two, it is a
reflexive string rewriting system. What sets it apart is that its molecules are
pattern-matching rules rather than little programs executed by a machine
(Stringmol) or operators on a stack of strings (SAC), and that its main
questions are the ones raised by [AlChemy](alchemy.md) and
[Tierra](tierra.md), whose experiments it set out to mirror: when do
self-replicators and self-maintaining organisations arise and last? Its
self-maintaining sets are those of [Kauffman's autocatalytic
sets](kauffman-autocatalytic-sets.md), but made of rewriting rules instead of
polymers that join and cut one another.

## How it works

### Molecules: strings of eight symbols

A molecule, called a *broadcast device*, is a string over an alphabet of eight
symbols:

- `0` and `1`, the information;
- `*`, which starts a rule, and `:`, which separates its condition from its
  action;
- three wildcards, written here `#`, `$` and `%` (the thesis and papers draw
  them as a diamond, a reversed triangle and a triangle);
- the quote `'`, which makes the next symbol literal.

A rule, a *broadcast unit*, runs from a `*` to the next `*`, and has the form
`*condition:action`. A string can hold several rules, or none; a string with no
complete rule is inert (a *null device*) and can only be a substrate. Symbols
outside any rule are ignored, like junk DNA.

### A reaction: match the condition, write the action

The enzyme's condition must match the *whole* substrate. `0`, `1` and quoted
symbols must match themselves. `$` matches one or more symbols at the start or
the end of the substrate; `#` and `%` match exactly one symbol (a final `#`
matches any non-empty tail). The action is then written out symbol by symbol:
`0`, `1` and quoted symbols are copied, and a `$` or `%` in the action is
replaced by whatever it matched in the condition. That is how a string carries
part of its substrate into its product (the papers call it *transposition*).
If the condition fails, the collision is *elastic*: nothing is made.

The operation is always

```
enzyme + substrate -> enzyme + substrate + product
```

so both reactants act as catalysts, and order matters: `a` acting on `b` is a
different reaction from `b` acting on `a`. The book's Table 11.1 shows what
this can express: copying, cutting a string (cleavage), extending it
(concatenation), turning an inert string into a rule (activation), or a rule
into an inert string (inhibition).

The shortest copier is `*$:$`. Its condition `$` matches any whole string, and
its action `$` writes it back, so it copies whatever it meets, itself included.
Out of the 8⁴ = 4,096 strings of length four, it is the only one that copies
itself. A longer string such as `*$0101:$0101` is a *specific* replicase: it
copies only strings that end in `0101`, which includes itself.

### A worked example: the seed network c0

The default network is the hand-designed seed of the published evolution
experiment, a set of four rules that Decraene and McMullin call `c0`:

```
s1 = *$0:$1    s2 = *$0:$0    s3 = *$1:$0    s4 = *$1:$1
```

Each reads "a string ending in x → the same string ending in y". Take `s1`
acting on `s3` (`*$1:$0`). The condition `$0` needs a string ending in `0`:
`*$1:$0` does, and `$` matches its first five symbols, `*$1:$`. The action
`$1` writes those five symbols and then `1`, giving `*$1:$1`, which is `s4`.
Now `s3` acting on `s1` (`*$0:$1`, which ends in `1`): `$` takes `*$0:$`, the
action `$0` gives `*$0:$0`, which is `s2`. And `s1` acting on `s2` gives `s1`
again: `s2` is copied with its last symbol changed to `1`.

In Chemart's species names each symbol becomes a letter or digit (`*` → `S`,
`$` → `D`, `:` → `C`, `#` → `H`, `%` → `P`, `'` → `Q`), so `s1 = *$0:$1` is
`m_SD0CD1`. The six reactions of `c0`, printed in the default output below,
are, written as enzyme on substrate:

```
s1 on s2 -> s1        s3 on s1 -> s2
s1 on s3 -> s4        s3 on s4 -> s3
s2 on s3 -> s3        s4 on s1 -> s1
```

In reaction-network form each is `enzyme + substrate -> enzyme + substrate +
product`; the first, for instance, is printed as
`m_SD0CD1 + m_SD0CD0 -> 2 m_SD0CD1 + m_SD0CD0`.

Every one of the four strings is made by some reaction among the four, and
nothing else is made: the set is *closed* and collectively autocatalytic.
Two of the strings could copy themselves (`s2` acting on another `s2` gives
`s2`, and likewise `s4`), but in this experiment self-replication was switched
off, so those collisions are elastic and the set keeps itself going only
collectively.

### The reactor

In the single-reactor model (thesis §4.2.1), two molecules are picked at
random, the first as the enzyme. If they react, the product is added; once the
reactor holds its maximum number of molecules, each new product replaces a
random molecule other than the two reactants. This random removal acts as a
*dilution flow*, and in the continuous limit the reactor follows the
*catalytic network equation* of Stadler, Fontana and Miller (1993; see
[random catalytic networks](random-catalytic-networks.md)): each species grows
at the rate its producing reactions fire, minus its share of the removals.
Products longer than a maximum length are not allowed.

Each product can be miscopied: every symbol, with a small probability `p_s`,
is replaced by another symbol, has a random symbol inserted after it, or is
deleted, with equal odds. Longer strings therefore suffer more errors.

### Cells and multilevel selection

The published evolution experiments add a second level. Each cell is a
separate reactor that only grows. When it meets a *division criterion* (it is
full, or it holds enough copies of a chosen *target* species), half its
molecules, drawn at random, go to a daughter cell, and a random other cell is
removed so the number of cells stays fixed. Cells whose networks reach the
criterion sooner divide more often and spread. Changing the criterion changes
what is selected, without saying how to achieve it. Random sharing at division
can also lose a rare species altogether, a *cellular mutation*. In the
experiments each cell ran on its own processor and cells raced in real time,
so a network that took less computer time per reaction also had an advantage.

Chemart implements the molecular level (the rule language, the closure of a
set of strings and the single reactor with mutation), not the cells. The
formal specification below summarises the rules and the reactor.

## Using it

The default call returns the closure of `c0`: the four strings and their six
reactions, with `status=complete` meaning no new string can be made.
`net.initial_state` gives 10 copies of each seed string, the amount the thesis
uses for `c0` (§7.2.2). Each reaction's `k` is a collision weight: for each of
the two orders in which the pair can meet, the fraction of the enzyme's
matching rules that give this product, summed. It is 1 for every reaction
here, and 2 when each string, acting on the other, makes the same product.
Integrated with the constant-total outflow the network declares, these
reactions are the catalytic network equation described above.

**The evolved cell types.** Book Table 11.2 lists the dominant cell types of a
run; pass their strings (the thesis glyphs are also accepted):

```python
net = chemart.generate_network("mcs-bl", strings=["*$0:$1", "*$0:$0", "*$#:$0", "*$#:$1"])
print(net.summary())
# mcs-bl: 4 species, 9 reactions, status=complete
```

This is `c1`: still closed, with nine reactions against six for `c0`. The
wildcard `#` lets `s5 = *$#:$0` and `s6 = *$#:$1` act on strings ending in
either symbol. `c2` and `c3` give the same nine reactions once `#` and `%`
are treated as the same symbol.

**The universal copier.** Self-replication must be allowed:

```python
net = chemart.generate_network("mcs-bl", strings=["*$:$", "0101"], self_replication=True)
# 2 species, 2 reactions:
#   m_SDCD + m_SDCD  -> 3 m_SDCD              (*$:$ copies itself)
#   m_SDCD + m_0101  -> m_SDCD + 2 m_0101     (*$:$ copies 0101)
```

**The elongation catastrophe.** The replicase `*$0101:$0101` and its one-symbol
mutant `*$0101:$00101` generate longer and longer strings. The maximum length
stops the closure: with `max_length=16` it ends at five strings of length 12
to 16 and nine reactions; with the default `max_length=500` the species budget
cuts it off (`max_species=200`, 10,100 reactions, `status=truncated`).

**The single reactor.** `method="soup"` runs the reactor for `steps`
collisions and returns the reactions that fired, with how often
(`r.count`); `net.extras["final_state"]` holds the final population and
`net.extras["analysis"]` counts collisions, productive ones and mutant
products. The default soup, 10,000 collisions of `c0` with `p_s=1e-5` and
capacity 1,000, takes about 2.5 s; 4,160 collisions were productive, and the
population ended at 300, 294, 210 and 196 copies of `s1`, `s3`, `s2`, `s4`.

To repeat the thesis' universal-copier experiment (§5.3: 100 copies of `*$:$`
among 900 random strings of length 10, capacity 1,000, no mutation), repeat
the string in the list and set one copy each:

```python
net = chemart.generate_network("mcs-bl", method="soup", seed=0,
        strings=["*$:$"] * 100, n_random=900, initial_copies=1,
        self_replication=True, p_s=0.0, steps=200000)
net.extras["final_state"].get("m_SDCD", 0)     # 11
```

After 200,000 collisions (7 to 11 s) the copier had fallen from 100 to 11, 15
and 0 molecules for seeds 0, 1 and 2.

For the catastrophe itself (§5.5), seed 100 copies of `*$0101:$0101` among 900
random strings with `p_s=1e-3` and 100,000 collisions (7 to 11 s): for seeds 0
and 1 the replicase was extinct by the end, and the mean string length had
risen from about 10 to 391 and 380 symbols, against a limit of 500. Runs grow
slower as strings lengthen; the thesis ran 5 million collisions per run, far
beyond a quick call here.

`n_random` and `random_length` add random seed strings; `max_length`,
`n_max`, `p_s` and `steps` are in the parameter table below.

## Results

### Self-replicators are rare and gain nothing (thesis ch. 5; ALife 2008)

**Spontaneous replicators.** Decraene seeded 30 runs with 100 random strings
of length 10 and let them react, with mutation, for 5 million collisions each
(thesis §5.2). Fifteen different self-replicating strings appeared, but none
ever reached more than a single molecule, and since in MCS.bl a string needs
to meet a second copy of itself to replicate, no replication was seen. The
thesis explains part of the rarity: of all strings of length four, the shortest
that can hold a rule, only `*$:$` copies itself. (The thesis and the ALife 2008
paper give that number of strings as "4^8 (65,536)"; with eight symbols it is
8⁴ = 4,096.) Chemart's tests check by enumeration that `*$:$` is the only one.

**The universal copier has no advantage.** Placed by hand, 100 copies among 900
random strings, `*$:$` did not take over: averaged over 30 runs it steadily
declined (§5.3). The thesis explains it with the catalytic network equation.
In the best case for the copier, where every other string is inert, the copier
makes copies of itself and of the others in exactly the proportion they
already have, so both have zero expected growth (eqs. 5.3–5.6) and only random
drift remains; any side reaction among the other strings tips the balance
against it. A copier that copies everything cannot select itself. Chemart's
tests check that the network of `*$:$` and one inert string has zero growth
rate at every composition, and the runs above show the decline.

**Specificity helps.** Replicases that copy only strings ending in `1`, `01`,
`101` or `0101` took over the reactor more often the longer their tag (thesis
§5.4, Figs. 5.2–5.3), because fewer random strings share the tag and can
parasitise them. Chemart can run these experiments (the strings are in the
parameter table) but does not test this result.

**The elongation catastrophe.** Seeded with the specific replicase
`*$0101:$0101` at 10% of the population, with mutation, the replicase first
filled the reactor and then collapsed; string length rose sharply and the
reaction rate fell until reactions stopped (§5.5, Fig. 5.4; ALife 2008). The
cause is the mutant `*$0101:$00101`. The replicase copies it, but the mutant
acting on the replicase also makes the mutant, and two mutants make a still
longer string `*$0101:$000101`, and so on. Each step is fitter at the
molecular level, but long strings are hit by more mutations, which break
their rules or make them so specific that nothing matches. Five changes to
the model (limits on length, on the wildcard, on the supply of symbols) did not
cure it (ALife 2008). Chemart's tests check this chain of reactions exactly,
and that the maximum length bounds it; the soup runs above show the
replicase's extinction and the growth of string length.

**Cells cure it.** With 32 cells of capacity 1,000, each seeded with 250
copies of `*$0101:$0101` and 250 random strings, cells infected by elongating
strings stopped dividing and were displaced, and no evolved cell suffered the
catastrophe (ALife 2008, 5 runs of at least 50 million collisions per cell).
Chemart does not model cells.

### Evolving a self-maintaining network (ACS 2011; thesis §7.2; book Fig. 11.3)

This is the experiment the book describes. Thirty-one cells, each on its own
processor, started from the hand-designed network `c0` and ran for one hour of
real time. Each cell could hold 10⁶ molecules, divided as soon as it held 200
copies of the target `s1 = *$0:$1`, and miscopied each product symbol with
probability 10⁻⁵; self-replication was disallowed. So a cell's fitness was how
fast its network made `s1` while keeping every species it needed.

In the reported run, 1,235 distinct cell types appeared, many of them not
self-maintaining, but only self-maintaining networks lasted: strings outside
the closed set were diluted away by divisions. Three takeovers happened, from
`c0` to `c1`, `c2` and `c3` (book Table 11.2), each in 4 to 9 seconds. The
authors measured the time a newborn cell needs to divide, over 4,500
incubations per type: 3.94×10⁻² s for `c0`, 3.65×10⁻² s for `c1`, 3.52×10⁻² s
for `c2` and 3.65×10⁻² s for `c3`. They concluded that the first takeover was
selection (the new strings `s5` and `s6` make `s1` faster), while `c1`, `c2`
and `c3` do the same thing (the wildcards `#` and `%` act alike in these
conditions, and the `%` in `c2`'s actions is ignored), so the later takeovers
were neutral drift. The durations alone could not tell selection from drift: a
model of drift in a population of 31 (the Moran process) predicts a takeover
in about 945 reproductions, with a standard deviation of about 506, and the
observed takeovers took 1,400 to 3,150, which the authors judge compatible
with drift. The network grew from 6
to 9 reactions with the same number of species and string length, which the
authors call "only a very limited growth of complexity". In ten more runs,
four ended in `c1` or `c3`, four in those plus extra strings that did not
help, and two in `c0` variants with no difference in behaviour.

Chemart reproduces the molecular side: its tests check that `c0`, `c1`, `c2`
and `c3` are closed with 6, 9, 9 and 9 reactions, the six reactions of `c0`,
and that `c1`, `c2` and `c3` have the same reactions once `#` and `%` are
identified. The takeovers, gestation times and cell statistics need the cell
model and are not reproduced. The general point, that selection can be aimed
at a chosen species by choosing the division criterion, is a property of the
cell model and so is not reproduced either.

A related table in the thesis (appendix D.2) shows that just `*$0:$1` and
`*$1:$0` make the other two strings of `c0` and nothing else; Chemart's tests
check this closure and the table's reactions, one of which the thesis prints
wrongly (see the implementation decisions).

### Crosstalk: merging networks (CEC 2009; ACS 2011; thesis §7.3)

Two self-maintaining networks put in one cell compete, and one displaces the
other, unless their molecules can react with each other (*crosstalk*); then
they can cooperate (CEC 2009). In ACS 2011 the authors seeded cells with `c1`
and a second evolved network `c4`, and let cells divide only once they held
200 copies each of both targets. First a combined network `c5` of 12 species
and 55 reactions dominated; then a mutant `c6` with three species fewer (9
species, 32 reactions) displaced it, and later `c7` (12 species, 66
reactions). Mean gestation times fell from 1.03 s to 0.80 s and 0.74 s. The
authors read this as networks of higher complexity that do both tasks at once,
arising by the merging of networks. Chemart can compute the reactions of a
mixed set of strings, but does not reproduce these experiments.

### What Chemart does not cover

Besides the cell level, discussed above, the thesis' spontaneous-mutation mechanism is left out (the table
11.2 experiment did not use it). The MCS.bl source code is not archived, so
the rule semantics were rebuilt from the thesis' twenty worked examples and the
source of Decraene's earlier broadcast-language implementation. One doubt
remains: six of the fifteen spontaneous replicators listed by the thesis do
not copy themselves under these rules (see the implementation decisions).

## Further reading

- Holland, J. H. (1975). *Adaptation in Natural and Artificial Systems*.
  University of Michigan Press. The broadcast language (in the ACS paper's
  citation of the 1992 edition, pp. 143–152).
- Stadler, P. F., Fontana, W. & Miller, J. H. (1993). Random catalytic
  reaction networks. *Physica D* 63(3–4), 378–392. The catalytic network
  equation used for the reactor.
- Decraene, J., Mitchell, G. G. & McMullin, B. (2006). Evolving artificial
  cell signaling networks using molecular classifier systems. BIONETICS 2006
  (book ref. [222]). The earliest of the MCS papers the book cites.
