## Introduction

Stringmol is an artificial chemistry in which every molecule is a short
program, and a reaction is what happens when two of these programs meet, stick
together and one of them runs. It was built by Simon Hickinbotham, Edward
Clark, Susan Stepney and colleagues at the University of York between 2009 and
2011. Their model was bacterial evolution. In living cells, as the ALife XII
paper puts it, "the phenotype includes the genotype-reading structures", and
Stringmol keeps that property: the molecules that copy strings are strings
themselves. For the first experiments the group took its cue from the RNA
world, in which one kind of molecule is both the template that gets copied and
the machine that copies it.

A Stringmol molecule is a string of letters and a few punctuation-like symbols,
for example the 64-symbol "seed replicase" that most experiments start from:

```
OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB
```

The capital letters are *template codes*: they do nothing themselves, but
decide where molecules stick together and where a program jumps. The seven
symbols `$ > ^ ? = % }` are *function codes*, the instructions. Two molecules
bind where a stretch of one roughly complements a stretch of the other; then
one of them, the enzyme, executes its string, copying and cutting symbols,
until it reaches the end instruction and the two let go. Bound to a copy of
itself, the seed replicase builds a third copy: two molecules in, three out.

Copying occasionally writes a wrong symbol, and because binding and jumps
rely on inexact matching, a mutated molecule often still works, only
differently. Running the chemistry for millions of time steps in one
well-mixed container, the York group watched replicases replace one another,
parasites wipe the population out, and pairs of molecules that could no longer
copy themselves survive by copying each other. The paper calls Stringmol an
"open-ended chemical system": one in which new kinds of molecule, with new
behaviours, can keep arising.

Banzhaf and Yamamoto describe Stringmol among the bio-inspired string
chemistries (book §11.1.2), after the [Molecular Classifier System](mcs-bl.md),
where strings also act as enzyme or substrate. Like
[Typogenetics](typogenetics.md) and [Laing's molecular
machines](laing-molecular-machines.md), it is an *automata chemistry*: its
molecules are programs. What distinguishes it is that binding is decided by a
sequence alignment of the kind used in bioinformatics, and that the enzyme
works on its partner through movable pointers. The digital-organism systems
Tierra and Avida, also self-copying programs, are archived in this catalog.

## How it works

### Symbols and complements

There are 33 symbols: the 26 template codes `A`–`Z` and the seven function
codes. Every template code has a *complement* 13 letters further on in the
alphabet, wrapping round: `A` pairs with `N`, `B` with `O`, `O` with `B`. Each
function code is its own complement. This is why the seed replicase needs two
different binding regions: a stretch of letters cannot bind to an identical
stretch, only to its complement.

Matching is graded rather than all-or-nothing. The symbols are arranged in a
fixed circle, and a symbol that sits next to the exact complement on that circle
counts as a near miss, one further away as a worse miss. To decide whether two
strings bind, Chemart, like the original, takes the complement of the
initiating string and aligns it against the partner with the
*Smith-Waterman algorithm*, the standard method in bioinformatics for finding
the best-matching pair of substrings of two sequences, allowing mismatches and
gaps. The alignment gives a score `s` over an aligned length `l`, and the
probability of binding is

```
P(bind) = min(s, l − 1.124) / (l − 1.124)      (0 if l ≤ 2)
```

so a long, nearly perfect match binds almost surely.

### Enzyme, substrate and pointers

Once bound, the molecule whose aligned stretch starts further from the
beginning of its string becomes *active* (the enzyme); the other is *passive*
(the substrate). The active molecule carries four pointers, each of which can
point into its own string or into its partner's:

- **I**, the instruction pointer: the symbol being executed;
- **F**, the flow pointer: a bookmark that I and the other pointers can jump to;
- **R**, the read pointer, and **W**, the write pointer, used by copying.

All four start at the beginning of the aligned stretch on each molecule. The
instructions are:

| code | name | what it does |
|---|---|---|
| `$` | search | move F to just after the best complementary match of the letters that follow |
| `>` | move | move R or W (when followed by `B` or `C`) to F; otherwise send I to just after F |
| `^` | toggle | switch F (or I, R, W with `A`, `B`, `C`) to the other string |
| `?` | if | skip the next instruction if R has run off the end of its string (or, with letters after it, with a probability set by how well they match at R) |
| `=` | copy | write the symbol under R at W, then advance both |
| `%` | cleave | cut the string at F; the tail becomes a new molecule |
| `}` | end | finish: the two molecules separate |

Template codes met by I are skipped. Searches and `if` tests use the same
inexact matching as binding, and succeed with a probability, so a slightly
changed label still finds its target most of the time.

### A worked example: the seed replicase copies itself

Two seed replicases, `R + R`. The complement of the initiator's first 18
symbols, `OOGEOLHHHRLUEUOBBB`, matches the partner's symbols 15 to 33,
`BBBRBXUUUDYGRHBLROO`, with score 9.73. The partner's stretch starts later, so
the partner becomes active, with binding probability
`min(9.73, 18 − 1.124)/(18 − 1.124) = 0.576`. Its pointers start at its own
position 15 and at position 0 of the substrate. Tracing the reaction step by
step in Chemart:

1. I walks over the template codes from position 15 to 35, doing nothing (21
   instructions).
2. At 36, `$BLUBO` searches for the complement of `BLUBO`, which is `OYHOB`,
   and finds it at the very end of the enzyme's own string: F = 64, the end.
3. `^B` toggles R onto the substrate, where it sits at position 0; `>C` moves W
   to F, the end of the enzyme.
4. `$=?>` is the copy loop. A `$` with no letters after it puts F on itself
   (position 46), so `>` will send I back to the `=` just after it. The loop
   is: `=` copies one symbol from the substrate to the end of the enzyme, `?`
   checks whether R has run off the substrate, and `>` jumps back.
   The trace (condensed) shows the enzyme growing by one symbol per turn:

   ```
   step 33  '='  R on substrate at 1   W at 65   enzyme length 65
   step 36  '='  R on substrate at 2   W at 66   enzyme length 66
   ...
   step 222 '='  R on substrate at 64  W at 128  enzyme length 128
   ```

5. With R past the end, `?` skips the jump. `$BLUBO` searches for `OYHOB` again
   and finds the first match, which now marks the border between the original
   enzyme and the new copy: F = 64.
6. `%` cuts there. The 64-symbol copy becomes a free molecule; `}` ends the
   reaction.

The event takes 232 instructions and yields `R + R → 3 R`. The copy starts
where the substrate's aligned stretch starts, here its first symbol; binding at
another site gives a truncated copy, as in the Results.

### Copying errors

Each `=` can mutate. With a small probability (the *substitution rate*) it
writes a neighbour of the read symbol on the circle instead of the symbol
itself; with a much smaller one (the *indel rate*) it either writes an extra
random symbol or skips one. Mutation happens only here: a new species is always
a faulty copy.

### The container

The reactor is one well-mixed container with a clock. At every time step each
molecule is visited once, in random order. A visited molecule may *decay*
(vanish, taking its partner with it if bound) with a fixed probability. If
energy is left, an unbound molecule may meet another unbound one, with a chance
that grows with their number and with the ratio of molecule to container size,
and then binds with `P(bind)`; a bound pair executes one instruction. Each bind
and each instruction costs one unit of energy, and a fixed amount is added
after every step. Energy caps the chemistry per step, decay removes molecules
continually, and the population settles where the two balance.

## Using it

The default run puts 100 seed replicases in the container for 3,000 steps at
the settings of the specification and the ALife XII paper (25 energy units
per step, decay probability 1/65² per step, substitution rate 10⁻⁵). The only
reaction is self-copying, `R + R → 3 R`, completed 101 times; no mutant
appeared. Species names are the sequences with the function codes written as
lower-case letters (`$ s`, `> m`, `^ t`, `? i`, `= c`, `% x`, `} e`);
`species.structure` holds the real sequence. The rest of the run is in
`net.extras`:

```python
ex = net.extras
ex["final_state"]           # {seed: 98}   unbound molecules at the end
ex["decayed"]               # {seed: 65}   unbound molecules that decayed
ex["aborted"][0]["count"]   # 6    reactions cut short because a molecule decayed
len(ex["in_progress"])      # 13   complexes still bound at the end
ex["analysis"]["population"][-1], ex["analysis"]["energy"][-1]   # (123, 49174)
```

The population grew from 100 to 123 while unspent energy piled up: with the
default container radius molecules rarely meet, so binding, not energy, limits
the chemistry. `ex["analysis"]` samples population, distinct strings and energy
about 500 times; the string count includes half-built copies, one per bound
complex, which is why it ends at 13. `ex["epochs"]` records every change of the
most abundant species (the papers' sweeps), and `ex["active_counts"]` how often
each reactant was the enzyme. Reactions and extras balance exactly against the
initial and final populations.

**The ALife XII cascade.** `method="closure"` does not simulate a container:
it takes a set of molecules, reacts every ordered pair once in isolation (the
first molecule initiates the bind), copies exactly, and repeats on the products.
Given species 9 of the ALife XII paper and its single-point mutant 29, it finds
the paper's cascade within seconds:

```python
from chemart.chemistries import stringmol as sm
sp9 = sm.ALIFE12_SPECIES_9
sp29 = sp9.replace("$BLUBO^", "$BLUBP^")            # the single point mutation
net = chemart.generate_network("stringmol", method="closure",
                               molecules={sp9: 1, sp29: 1}, max_species=8)
```

```
complete 8 species, 23 reactions
9 + 9 -> 9 + 9 + 9
9 + 29 -> 9 + 29 + 29
29 + 9 -> 30 + 9
29 + 29 -> [129-mer] + 29
30 + 9 -> 30 + 9 + 31
```

(first five reactions, renamed with the paper's species numbers). Order
matters: the alignment starts from the initiator's complement, so `9 + 29` and
`29 + 9` bind at different sites and do different things.

**Mutation and collapse.** Raise the substitution rate into the range of the
authors' configuration files and shrink the container so that molecules meet
often, and a run shows sweeps and a collapse (this one takes 15–30 seconds):

```python
net = chemart.generate_network("stringmol", seed=1, cell_radius=100, steps=40000,
                               substitution_rate=1e-3, indel_rate=1e-4)
ex = net.extras
print(net.summary().splitlines()[0])
print("extinct:", ex["extinct"], "after", ex["time_steps"], "steps")
print("population:", ex["analysis"]["population"][::50])
length = {s.id: len(s.structure) for s in net.species}
for t, sid in ex["epochs"][:3]:          # time, length and name of each new leader
    print(t, length[sid], sid)
```

```
stringmol: 106 species, 371 reactions, status=observed
extinct: True after 34593 steps
population: [100, 278, 419, 452, 201, 79, 32, 12, 4]
0 64 OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROOREsBLUBOtBmCscimssBLUBOxeOYHOB
8160 22 tBmCscimssBLUBOxeOYHOB
8560 42 OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROOREsOYHOB
```

After about 8,500 steps the seed replicase loses first place to a 42-symbol
molecule: its binding regions followed directly by `$OYHOB`, with the copy
program deleted. In isolated reactions this molecule is copied when a replicase
acts on it but copies nothing itself, and it is shorter, so it is copied faster.
It spreads, the replicases decline, and the container empties. The ALife XII
paper describes this end in words; this run is a small, fast analogue at 100
times its mutation rate, not a reproduction.

`method="soup"` draws pairs from a constant-size population and runs each
reaction to completion at once, with no energy or decay. Container runs of the
papers' length (about 350 molecules for 10⁵ to 10⁷ steps) are far too slow in
pure Python.

## Results

**A replicase that copies itself.** The technical specification (version 0.2,
2010) defines the chemistry and gives a hand-designed seed replicase; the ALife
XII paper's seed is 65 symbols long and takes 240 time steps to build a copy.
Chemart's tests check that the specification's seed, the seed of the authors'
configuration files and ALife XII species 9 each give `R + R → 3 R`, and that
the closure of the seed is exactly that one reaction.

**Invasion when rare.** In the ECAL 2009 paper (published 2011) a replicase
`r` with imperfectly matching binding regions (binding probability 0.293
between two copies) met a one-mutation variant `m` whose sites match perfectly
(probability 1, also with `r`). A single `m` added to an equilibrium population
of `r` took over in 88 of 100 runs. That paper used an earlier binding formula,
`(s/l)^l`, and length-dependent decay; Chemart implements version 0.2 and does
not repeat the experiment.

**Diversity from a monoculture.** The ALife XII paper (Hickinbotham et al.,
2010) is the main study. Starting from seed replicases with a substitution
probability of 10⁻⁵ per copied symbol, 25 energy units per step and a
population of about 350 molecules, it ran 1,000 trials, each until no molecule
was left. A new species appeared every 18,700 time steps on average, and most
were lost again quickly. The authors classified what they saw by inspecting the
population plots:

- *Extinction.* No trial kept going for ever. The modal extinction time was
  750,000 time steps, with about 40 new species produced on average by then.
  The usual cause was a parasite, a molecule that cannot copy itself but that
  the replicase copies faster than itself.
- *Characteristic sweeps.* A mutant drives the dominant species out in less
  than 50,000 time steps; these are the main cause of change.
- *Drift*, a neutral mutant rising by chance, in 92 trials; *slow sweeps* in
  52; *rapid sweep sequences*, in which a mutant triggers a cascade of new
  dominant species, in 31.
- *Sub-populations* of more than 50 molecules beside the dominant species, in
  nearly every run. In 26 trials one outlasted a sweep, which the authors read
  as sub-populations living off the dominant species as non-lethal parasites.
- *Hypercycles* in 30 trials: two or more species that cannot copy themselves
  but copy each other. Emergent hypercycles occurred in 8 trials, spontaneous
  ones in 15 and multispecies ones in 14.

The paper traces one hypercycle, in trial 277, to a single mutation. Species
29 differs from 9 only in `$BLUBO` becoming `$BLUBP`. The damaged search lands
one symbol short, so when 29 copies 9 it writes the copy over its own last
symbol, and the second search no longer finds the cut site: the result is
species 30, 29 with a whole copy of 9 attached. When 30 meets 9 the binding
site shifts and the copy skips the first 14 symbols of 9, giving the short
species 31. Species 31 lived as a sub-population for about 5,750,000 steps
before becoming one partner of a hypercycle that lasted about 3 million steps.
Chemart's tests reproduce each of these reactions exactly, in isolated
reactions and in the closure, including the exact sequences of 30 and 31; the
specification's own `$BLUBP` mutant of the seed, which also makes a nearly
double-length molecule, is tested too.

Chemart does not reproduce the statistics of the 1,000 trials, which need
hundreds of thousands to millions of steps per trial. Its tests check the
ingredients instead: copy mutations move only to neighbouring symbols on the
circle, energy influx and decay hold a population steady (with no energy
nothing reacts and the population decays), and the enzyme role goes to the
later binding site. The machine is a line-by-line port of the authors' C++
program, checked step by step against the compiled original on seven container
runs. One difference from the book: its Fig. 11.5, taken from the
specification, arranges the symbols on the circle in a different order from
the authors' program (the complements are the same). Chemart follows the
program, which produced the published runs; the implementation decisions list
this and the other discrepancies.

## Further reading

- Eigen, M. & Schuster, P. (1977). The hypercycle: a principle of natural
  self-organization. Part A: Emergence of the hypercycle.
  *Die Naturwissenschaften* 64, 541–565. The source of the hypercycle concept
  used in the ALife XII paper.
