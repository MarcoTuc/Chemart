## Introduction

This chemistry composes short pieces of music. Its molecules are fragments of
music: single notes, runs of notes, one-bar melodies with a chord written under
them, and sequences of bars. Its reactions are the steps a composer might
follow: extend a melody by one note, give a bar a chord, fix a note that clashes
with that chord, join two bars whose chords follow each other naturally, and
cut a finished four-bar phrase out of a longer sequence. Put a few thousand
such molecules in a well-stirred pot, let them collide at random, and phrases
in a fixed, simple style come out: a smooth melody of eighth notes in C major,
one chord per bar, and a chord progression that follows textbook cadences.

It was built by Tomoya Miura and Kazuto Tominaga at the Tokyo University of
Technology and published in 2006 (GWAL-7, the German Workshop on Artificial
Life). Their question was not really about music. They observed that running
an artificial chemistry by random molecular collisions is like running a
*nondeterministic algorithm*, a procedure with choice points where any of
several options may be taken: picking which molecules collide is picking an
option. Music composition, with a choice at every note and chord, was their
test case for using an artificial chemistry as a platform to write down and
execute such algorithms. The randomness a composer's program would otherwise
have to call explicitly comes for free from the chemistry.

The point the authors stress is that there is no fitness function. Genetic
algorithms for music generate many candidates and score them against the rules
of music theory, discarding the bad ones. Here the rules of the style *are* the
reactions, so every phrase that comes out obeys them by construction, and
fragments that would break the style are repaired or reused rather than thrown
away. The price, which the authors name as the main limitation, is that the
designer must already know how to compose in the style.

It is a simulation model, an application of a general formalism: the
chemistry underneath is Tominaga's pattern-matching-and-recombination chemistry,
catalogued separately as [Tominaga's stackable strings](tominaga-stacked-strings.md).
That entry is the general-purpose notation; this one is one particular set of
65 rules and an initial pot written in it. Banzhaf and Yamamoto give it one page
in their chapter on applications (book §16.4, "Music Composition Using
Algorithmic Chemistries"), next to [the naming game](naming-game-ac.md) and
[proof search](proof-ac.md). Among the catalog's applications it is closest in
spirit to [the Chemical Casting Model](ccm.md), where local rules also replace a
global algorithm; but CCM still scores candidates locally to solve an
optimisation problem, while here nothing is ever scored. The book also mentions
a later polyphonic system by Tominaga and Setomoto (2008); Chemart implements
only the 2006 homophonic one.

## How it works

### Molecules are stacks of lines

A molecule (the paper says *object*) is one line of *elements*, or several
lines stacked and aligned column by column. It is written line by line as
`displacement#elements/`, where the displacement says how many columns a line is
shifted relative to the first. Some molecules from a real run:

```
0#E4/                                   a seed: a single note
0#E4F4E4/                               a short melody
0#C4D4C4D4E4F4E4D4/0#SFFFFFFF/          one bar with a chord under it
0#Degree/0#D4E4/                        a "step from D4 up to E4" token
```

The elements are:

- **notes** `C4 D4 E4 F4 G4 A4 B4 C5`: the C major scale over one octave, each
  an eighth note (there are no other note lengths);
- **chords** `C Dm Em F G Am Bm`: the triad on each degree of the scale (`Bm`
  stands for the diminished B chord);
- **chord functions** `T S D D2`: tonic, subdominant, dominant and second
  dominant, the harmonic role of a chord (see below);
- **control elements** `Degree Chord Avoid Start Stop Dummy`, which mark the
  "reagent" molecules and finished phrases.

A bar is a two-line molecule. The first line holds its eight notes. The second
line starts with the chord's function, under the first note, and then repeats
the chord name seven times, so that every later note has its chord directly
beneath it. In `0#C4D4C4D4E4F4E4D4/0#SFFFFFFF/` the chord is F acting as
subdominant (`S`). The column alignment is what lets a rule ask "is there an F4
standing over a C chord?".

### Reactions are pattern rewrites

A rule has patterns on both sides. A pattern looks like a molecule with
*wildcards*: `<1>` matches one element, `<*1>` at the start of a line matches
any run of elements (possibly empty) extending leftwards, `<1*>` at the end of a
line does the same rightwards, and `<1..7>` is shorthand for `<1><2>...<7>`.
The molecules matched on the left are consumed, and the right side is built
from the matched pieces. Each wildcard appears exactly once on each side, so
nothing is created or destroyed: every element that goes in comes out, just as
atoms do in a chemical reaction. Rules take one or two molecules.

### The five steps of the composition, as reactions

The paper first writes the composition as a five-step procedure and then turns
each step into a group of rules, which it numbers by the examples it prints
(the formal specification below uses those numbers). Every example below is a reaction that fired
in the default Chemart run (`seed=1`).

**1. Grow a melody by steps (the paper's rule 1).** The pot holds 1,400 `Degree` tokens, one kind
for each move to the neighbouring note up or down. A melody ending in G4 meets
a "G4 to F4" token and takes the new note; the spent token is left behind:

```
0#G4/ + 0#Degree/0#G4F4/ -> 0#G4F4/ + 0#Degree/0#G4/
```

Because the only tokens are steps between neighbouring notes, a melody can
never jump (musicians call this *conjunct motion*). There are 8 such rules, one
per note.

**2. Cut off a bar and give it a chord (rules 2–5).** Once a melody is at least nine notes
long, a `Chord` token can cut its first eight notes off as a bar and write its
chord underneath. The chord must contain the bar's first note: here the melody
starts on C, and F (F, A, C) contains C:

```
0#C4D4C4D4E4F4E4D4E4F4E4/ + 0#Chord/0#SFFFFFFF/
    -> 0#C4D4C4D4E4F4E4D4/0#SFFFFFFF/ + 0#E4F4E4/ + 0#Chord/
```

The leftover `0#E4F4E4/` keeps growing. There is one chord token per chord and
function, eleven in all, because some chords have two functions (Am is both
tonic and subdominant). This gives 38 rules.

**3. Repair clashing notes (rules 6 and 7).** Each chord except F has an *avoid pitch*, a note
that sounds dissonant against it (F4 against C, C4 and C5 against G, and so
on). Since the chord was chosen only from the first note, the rest of the bar
may contain avoid notes. A two-stage repair swaps one for a placeholder
`Dummy`, then swaps the `Dummy` for whatever note an `Avoid` token happens to
carry:

```
0#B4C5B4C5B4A4B4C5/0#DGGGGGGG/ + 0#Avoid/0#Dummy/
    -> 0#B4C5B4C5B4A4B4Dummy/0#DGGGGGGG/ + 0#Avoid/0#C5/
0#B4C5B4C5B4A4B4Dummy/0#DGGGGGGG/ + 0#Avoid/0#D4/
    -> 0#B4C5B4C5B4A4B4D4/0#DGGGGGGG/ + 0#Avoid/0#Dummy/
```

The last C5 of this G-chord bar became a D4. The replacement is random, so it
can be another avoid note, which a later collision will replace again. This is
also the only way a melody gets a jump. 9 rules.

**4. Join bars along cadences (rule 8).** Chords are grouped by function, and some
successions of functions are what the paper calls typical *cadences*,
progressions that sound natural: tonic–dominant–tonic,
tonic–subdominant–tonic, tonic–second dominant–dominant–tonic. The paper turns this into six allowed transitions:
T→D2, T→D, T→S, D→T, D2→D and D2→T. A rule for each lets a sequence whose last
bar has the first function absorb a sequence whose first bar has the second.
Here a tonic Em bar takes a second-dominant F bar:

```
0#E4F4G4A4G4F4E4D4/0#TEmEmEmEmEmEmEm/ + 0#C4D4C4D4C4D4E4D4/0#D2FFFFFFF/
    -> 0#E4F4G4A4G4F4E4D4C4D4C4D4C4D4E4D4/0#TEmEmEmEmEmEmEmD2FFFFFFF/
```

**5. Cut out phrases (rules 9 and 10).** A phrase must start on a tonic bar. When a sequence of
at least five bars starts with one, a `StartStop` token cuts off its first four
bars and caps them with `Start` and `Stop`; the rest goes back into the pot.
If a sequence starts on a non-tonic bar, one of three rules detaches that bar
so that it can be used elsewhere. In the default run, a six-bar sequence
starting on a G (dominant) bar lost that bar, and the remaining sequence, which
now began on a tonic Em bar, yielded the run's phrase:

```
0#StartE4F4E4A4E4F4E4D4D4E4D4E4F4E4F4G4D4F4D4E4F4G4A4B4G4A4G4F4E4F4E4D4Stop/
1#TEmEmEmEmEmEmEmD2DmDmDmDmDmDmDmDBmBmBmBmBmBmBmTEmEmEmEmEmEmEm/
```

(one molecule, shown on two lines). Read bar by bar, that is Em, Dm, Bm, Em, with
the functions tonic, second dominant, dominant, tonic: the textbook cadence
T–D2–D–T. The chord line is displaced by 1 because `Start` now occupies the
first column.

The rule counts add up to the paper's 65: 8 + 38 + 9 + 6 + 1 + 3.

### The steps interleave

Nothing orders these steps. Melodies grow while other bars are being given
chords, repaired, joined and cut, and several sequences grow in the same pot at
once. What makes the whole thing work is supply: 1,826 starting molecules of 43
kinds, in counts the authors call "somewhat arbitrary", chosen so that a few
phrases appear in reasonable time. A finished phrase normally takes no further part: the `Start` and `Stop` caps
shift its columns so that none of the growing, joining or cutting rules fits
it. The exception is repair: a phrase that still contains an avoid note, or a
pending `Dummy`, can still be acted on by the repair rules.

### The reactor

The paper's simulator picked molecules at random and let them collide, and the
authors deliberately report only the products, not rates or time courses: they
say it lacks the theoretical foundation to simulate quantitative behaviour.
Chemart does the same with its own sampling scheme. At each step it picks a
rule with probability proportional to the number of ways its reactants could
be chosen from the pot, then picks the reactants in proportion to their copy
numbers. No rule has a rate constant. The run stops when a set number of
phrases has been finished, or after a fixed number of steps.

## Using it

The default call runs the published system (all 65 rules, the published
starting pot) until one phrase is finished. With `seed=1` that takes 893
reactions, 705 of them growing melodies, and about 2 seconds. The 772 species
and 719 reactions in the summary are the distinct molecules that appeared and
the distinct reactions that fired; each reaction carries its firing count. The
first reactions listed are all step 1: seeds and short melodies growing.

The music is in `net.extras["analysis"]["phrases"]`, one entry per finished
phrase:

```python
p = net.extras["analysis"]["phrases"][0]
p["chords"], p["cadence"]   # (['Em', 'Dm', 'Bm', 'Em'], ['T', 'D2', 'D', 'T'])
p["music"]
# 'phrase: Em(T) E4 F4 E4 A4 E4 F4 E4 D4 | Dm(D2) D4 E4 D4 E4 F4 E4 F4 G4 |
#  Bm(D) D4 F4 D4 E4 F4 G4 A4 B4 | Em(T) G4 A4 G4 F4 E4 F4 E4 D4'
```

`p["notes"]` is the 32-note melody and `p["object"]` the molecule itself. The
jumps E4→A4 and D4→F4 are repaired avoid notes. This phrase still has F4s
under Em, its avoid pitch, so the repair rule could still act on it. Also in
`analysis`: the design totals (`elements` 25, `initial_object_kinds` 43,
`initial_objects` 1826, `rules` 65) and `failed_collisions`, draws that could
not react. `net.extras["rules"]` lists the 65 rules in the paper's notation,
`net.extras["reaction_rules"]` names the rule behind each reaction,
`net.extras["final_state"]` is the pot at the end, and
`net.extras["conservation"]` holds one conservation law per element kind. Each
species' `structure` field gives its musical reading (`notes: ...`,
`bars: ...` or `phrase: ...`).

**Several phrases in one run.** The paper's runs collected three to five
phrases. Set `phrases=0` to run a fixed number of steps instead of stopping at
the first phrase:

```python
net = chemart.generate_network("music-ac", seed=5, steps=3000, phrases=0)
[p["chords"] for p in net.extras["analysis"]["phrases"]][:5]
# 10 phrases in about 3.5 s; the first five:
# [['Em', 'F', 'Em', 'Em'], ['C', 'Em', 'Am', 'Dm'], ['Am', 'Dm', 'Em', 'Am'],
#  ['Am', 'Bm', 'C', 'Em'], ['Em', 'Dm', 'G', 'Em']]
```

The fifth, Em, Dm, G, Em, is the chord sequence of the paper's Figure 3(2).
Longer runs add little. With `seed=1, steps=20000, phrases=0` (about 10 s) all
1,400 step tokens are used up, 167 bars are made, but only 12 phrases come out.
Most of the remaining collisions are a futile cycle: a non-tonic bar is joined
to a sequence and then detached again (8,897 and 8,795 firings of the joining
and detaching rules).

**Without the repair step.** `avoid_notes=False` drops the 9 repair rules,
leaving 56. Melodies are then purely stepwise:

```python
net = chemart.generate_network("music-ac", seed=5, avoid_notes=False)
# phrase: Em(T) E4 F4 G4 A4 B4 A4 G4 F4 | F(D2) A4 B4 C5 B4 C5 B4 C5 B4 |
#         Bm(D) D4 E4 F4 E4 F4 G4 A4 B4 | Am(T) A4 G4 A4 G4 F4 E4 D4 E4
```

**Allowing subdominant to tonic.** The paper lists T–S–T as a typical cadence
but gives no rule for S→T (see Results). `cadences="with-s-t"` adds it (66
rules). With `seed=5, steps=3000, phrases=0` it gives 15 phrases in about
1.5 s, among them Am, Am, Em, Bm with functions T, S, T, D, which the published
rules cannot produce.

**Other sizes.** `notes_per_bar` and `bars_per_phrase` resize the bar and the
phrase and rewrite the rules to match; `notes_per_bar=4, bars_per_phrase=2`
gives two-bar phrases of four notes per bar. `copies` multiplies the whole
starting pot. Runs of a few thousand steps take seconds; the time grows with
the number of distinct molecules, so budgets of tens of thousands of steps take
ten seconds or more.

## Results

**Phrases in a defined style, without a fitness function.** The 2006 paper's
result is qualitative. The system "successfully generated musical phrases"
without human intervention: a typical run gave three to five four-bar phrases
in about an hour on a 1.8 GHz PowerPC G5. The paper prints four of them, Figure
1 (chords C, G, C, G) and Figure 3 (Am, Am, C, G; Em, Dm, G, Em; C, Bm, Am, F).
Its discussion claims every phrase observes the style by construction: a smooth
melody, little dissonance from avoid notes, and a chord progression following
the cadence rules, with each phrase starting on a tonic bar. There is no
listening test or other evaluation of musical quality. Chemart's tests check
the style on every phrase its runs produce: four bars of eight notes, a chord
in each bar that contains the bar's first note and has the stated function, a
tonic first bar, and a function sequence made only of allowed transitions.
They also check the design totals (43 kinds of starting molecule, 1,826
molecules, 65 rules), the starting counts, the paper's worked examples of
matching and recombination, and that seven of the printed rules are generated
character for character. Chemart reproduces the phrases, not the timings: its
runs take seconds, and the paper's hour on its own simulator is not a
meaningful comparison.

**Figure 3(1) contradicts the published rules.** The paper's six cadence rules
are what its total of 65 requires, but they omit S→T. Checking every way of
assigning functions to the chords of the printed phrases, Figure 1, Figure 3(2)
and Figure 3(3) are derivable, but Figure 3(1), Am, Am, C, G, is not: it needs
T, S, T, D. The paper does list T–S–T among its typical cadences, but its
rules do not allow it, and the text does not say how Figure 3(1) was produced. Chemart
keeps the published six by default and offers the seventh as
`cadences="with-s-t"`; the tests check both derivability results.

**Finished phrases are final products.** The paper says a phrase capped with
`Start` and `Stop` "is no longer recombined by any rule", and that it collected
these as the output, since the system has no stopping condition. The tests
confirm that no rule matches the Figure 1 phrase. The claim is not quite
general, though: a phrase that still carries an avoid note can still be
repaired, as the default run's phrase shows. The paper promises only "little"
dissonance, so this is not an error in it, but the phrase-cutting rule does not
wait for repairs to finish. Chemart's runs also show phrases cut while a
`Dummy` was pending (one appears in the `with-s-t` run above), so a phrase's
note list can contain `Dummy`.

**Jumps come from repairs.** Each phrase in the paper's figures has a jump in
an otherwise stepwise melody, which the authors attribute to avoid notes being
replaced by random notes. The tests confirm both halves: with the repair rules
off, every interval inside a bar is a single step; with them on, a longer run
(`seed=5`, 3,000 steps) produces jumps and fires the repair rules.

**Bars are reused, not discarded.** The paper's Figure 2 walks through a
sequence that ends up starting on a dominant bar, which cannot begin a phrase,
so the bar is cut off and reused. The authors present this as an advantage over
genetic algorithms, where unfit individuals are simply dropped. The tests check
the detaching rule on such a sequence, and that it fires in a longer run.

**Elements are conserved.** Every rule conserves elements, which the paper
builds in from the start. The tests verify one conservation law per element
kind on the network of a run.

**Randomness and parallelism for free.** The authors' broader claim is that
the chemistry's own nondeterminism supplies the choices a composing program
would have to make explicitly, and that phrases grow in parallel without any
parallel programming. Chemart shows several phrases being built in one pot
(10 in the 3,000-step run above), which the slow test checks; the rest of the
claim is an argument about programming, not a measurement.

**Later work.** Tominaga and Setomoto (2008) extended the approach to
polyphonic phrases: two melodies, upper and lower, both generated by the
chemistry. Book Figure 16.11 shows three phrases from each system. The paper
is not openly available. The one open description, in Fernández and Vico's
2013 survey of algorithmic composition, says it encoded polyphonic compositions
in the strings and rules of counterpoint in the reactions, and that "the
aesthetical value of the resulting compositions varied widely". Tominaga's lab
page (now only on the Wayback Machine) also offered samples with various note
durations, non-repeating phrases and generated beats, with the beat method "to
be published soon". Chemart implements none of this: without the 2008 paper
there is no rule set to reconstruct, so the `mode` parameter has only the
homophonic choice.
