## Introduction

The *naming game* is a model of how a group comes to share a word for
something without anyone deciding it. Luc Steels and colleagues at the VUB
Artificial Intelligence Lab in Brussels studied it in the late 1990s with
populations of software agents and then robots. Two agents meet; the speaker
names an object; the hearer tries to understand. In Steels' description
(Steels 2000), when a game succeeds "the scores of the associations that were
used go up and those of competing associations go down", and when it fails the
scores of the associations used go down. This reward-the-winner,
punish-the-rivals rule is what he calls *lateral inhibition*. Repeated over many games, one word
comes to dominate and its rivals die out: the population has a shared lexicon
(a common vocabulary).

De Beule, Hovig and Benson (2011) rewrote the simplest version of this game as
a chemical reaction network. Their paper, in the journal *Biosemiotics*, is
about codes in biology in general: it proposes a formal model of how a
convention can emerge among the users of a code, with examples from language
and from the immune system. In their chemistry, a *meaning* molecule `M`, *word*
molecules (`S1`, `S2`, … for each candidate name) and *adaptor* or *code*
molecules (`C1`, `C2`, …, one per name, standing for the association that links
the meaning to that word) float in a well-stirred reactor. Words make adaptors,
adaptors make words, a word that meets its own adaptor is copied, and a word
that meets a rival's adaptor makes one of the two switch sides. The authors
turned the reactions into ordinary differential equations (ODEs) and, as the
book reports, showed that they give the same lateral-inhibition behaviour as
the agent-based naming games, converging to a shared lexicon.

It is a small, fixed **simulation model**: a handful of species and mass-action
reactions, with no new kinds of molecule ever created. The book presents it in its
chapter on applications of artificial chemistries (section 16.3.3, "Language
Cognition"), as one of the few cases where an artificial chemistry has been
used for language, and points out that such models can bridge ODE-based and
agent-based models of language. It is the only entry in the catalog from the
book's section on language. Its closest relatives are the generic selection models: in
[the selection equation](selection-equation.md) and
[the replicator equation](replicator-equation.md), too, competing species
exclude each other, but there the competition is written in directly as
fitness or payoff, whereas here it has to come out of the word–adaptor
reactions. Its siblings in chapter 16, [music composition](music-ac.md) and
[proof search](proof-ac.md), are applications of a very different kind.

Chemart implements the five reaction schemes of the book. The book gives the
schemes only, with no rate values, no initial state and no numerical results, and its printing
of the last two schemes is garbled (see *How it works*). As the runs below
show, with all rates equal (Chemart's defaults) the network does not converge
at all; convergence to one name appears only when the two "switching" rates
differ in a particular direction.

## How it works

### Species

There is one meaning, `M`, and `N_s` candidate names for it (three by
default). For each name `j` there are two species:

- `Sj`, a word: the name as it circulates, as a sign;
- `Cj`, an adaptor (or code): the association that maps the meaning `M` to the
  word `Sj`, the chemical counterpart of an entry in an agent's lexicon.

The quantity that tells how popular name `j` is, is the total amount of its
molecules, `Sj + Cj`. A shared lexicon means one name holds everything and the
others are at zero.

### Reactions

The book prints five reaction schemes (its eqs. 16.2–16.6), each with its own
rate constant. For name `j` (the default network, `N_s = 3`, shown for name 1):

```
M + C1 -> S1        rate ks      a code produces its word for the meaning
M + S1 -> C1        rate kc      a word produces the code for the meaning
M + S1 + C1 -> 2 S1 rate ka      a word that meets its own code is copied
```

and, for every rival name `k ≠ j`, two reactions that resolve a mismatch
between a word of one name and a code of another:

```
M + S1 + C2 -> M + S1 + C1   rate kappa1   the code switches to the word's name
M + S1 + C2 -> M + S2 + C2   rate kappa2   the word switches to the code's name
```

These two lines are the first rival pair of the default network; `S1 + C3`
gives two more, and the same holds for each name. `M` is a reactant in every
reaction. In the first three it is used up, as the book prints them; in the
mismatch reactions it comes out again, acting as a catalyst.

The mismatch reactions need a word of explanation. In the book's text both
eqs. 16.5 and 16.6 read `M + Sj + Cj -> M + Sj + Cj`, which does nothing; the
indices that marked the mismatch were lost. The book's prose says these
reactions "eliminate mismatching words (resp. adaptors) by replacing them with
matching ones". Chemart reads that as: when word `Sj` meets a rival code `Ck`,
either the code is replaced (`kappa1`) or the word is (`kappa2`). So the same
encounter has two possible outcomes. (The book's "resp." ordering would pair
16.5, `κ1`, with replacing the word; Chemart's `kappa1` replaces the code. Only
the labels differ.)

The reactor is deterministic mass action: each reaction runs at its rate
constant times the product of its reactant concentrations, and concentrations
change continuously by the resulting ODEs.

### What the reactions conserve

Two things follow straight from the reaction list, and they decide what the
network can do.

**Each name's total is untouched by the first three reactions.** `M + C1 -> S1`,
`M + S1 -> C1` and `M + S1 + C1 -> 2 S1` only turn name 1's codes into name 1's
words or back. Only the mismatch reactions move material between names.

**The two mismatch outcomes pull in opposite directions.** In the encounter
`S1 + C2`, the `kappa1` outcome gives name 1 one molecule from name 2, and the
`kappa2` outcome gives name 2 one molecule from name 1. Adding up all
encounters, the total of name `j` changes at the rate

```
d(Sj + Cj)/dt = (kappa1 − kappa2) · M · (Sj · C_tot − Cj · S_tot)
```

where `S_tot` and `C_tot` are the total amounts of words and of codes. So:

- if `kappa1 = kappa2`, every name keeps exactly the amount it started with,
  and no lexicon forms;
- if `kappa1 > kappa2`, a name grows when its share of the words exceeds its
  share of the codes; in the run below this snowballs until one name holds
  everything;
- if `kappa2 > kappa1`, the sign flips; in the run below the names end up
  with equal amounts.

There is a second limit. `M` is used up by the first three reactions and never
made. When it runs out, every reaction stops, because each needs `M`. To see
the long-run competition, `M` has to be held constant, as if the meaning were
kept on offer. The book does not say which the original paper did.

## Using it

The default call above builds the network for three names: 7 species and
21 reactions, `3·N_s + 2·N_s·(N_s − 1)` in general, all with rate constant 1.
No initial concentrations are attached (`net.initial_state` is `None`) and
`net.extras` is empty; the book gives neither, so you supply them. Species
names follow the book: `M` the meaning, `Sj` the words, `Cj` the codes.

`chemart.simulate.ode` integrates the mass-action ODEs. The script below
starts from a mixture in which names 1 and 2 each hold 0.4 and name 3 holds
0.2, split unevenly between words and codes. To hold `M` fixed it lists `M`
among the network's buffered species (`net.extras["buffered"]`), which the
simulator keeps constant. It prints each name's total over time, in four
settings: defaults with `M` used up, defaults with `M` held fixed, and
`kappa1` or `kappa2` doubled with `M` held fixed. It takes about a second.

```python
import chemart
from chemart import simulate

x0 = {"M": 1.0, "S1": 0.3, "C1": 0.1, "S2": 0.25, "C2": 0.15, "S3": 0.2}
times = [0, 1, 5, 20, 100, 400]

def report(label, net, hold_M):
    if hold_M:
        net.extras["buffered"] = ["M"]                  # the meaning is kept on offer
    traj = simulate.ode(net, 400, x0=x0, points=401)    # a frame every time unit
    print(label)
    for t in times:
        x = traj.frames[t].state
        names = "  ".join(f"name{j}={x.get(f'S{j}', 0) + x.get(f'C{j}', 0):.3f}" for j in (1, 2, 3))
        print(f"  t={t:<4} M={x.get('M', 0):.3f}  {names}")

report("defaults, M consumed", chemart.generate_network("naming-game-ac"), False)
report("defaults, M held", chemart.generate_network("naming-game-ac"), True)
report("kappa1=2, M held", chemart.generate_network("naming-game-ac", kappa1=2.0), True)
report("kappa2=2, M held", chemart.generate_network("naming-game-ac", kappa2=2.0), True)
```


```
defaults, M consumed
  t=0    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=1    M=0.339  name1=0.400  name2=0.400  name3=0.200
  t=5    M=0.004  name1=0.400  name2=0.400  name3=0.200
  t=20   M=0.000  name1=0.400  name2=0.400  name3=0.200
  t=100  M=-0.000  name1=0.400  name2=0.400  name3=0.200
  t=400  M=-0.000  name1=0.400  name2=0.400  name3=0.200
defaults, M held
  t=0    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=1    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=5    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=20   M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=100  M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=400  M=1.000  name1=0.400  name2=0.400  name3=0.200
kappa1=2, M held
  t=0    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=1    M=1.000  name1=0.401  name2=0.388  name3=0.212
  t=5    M=1.000  name1=0.406  name2=0.391  name3=0.203
  t=20   M=1.000  name1=0.428  name2=0.404  name3=0.169
  t=100  M=1.000  name1=0.647  name2=0.329  name3=0.023
  t=400  M=1.000  name1=1.000  name2=0.000  name3=0.000
kappa2=2, M held
  t=0    M=1.000  name1=0.400  name2=0.400  name3=0.200
  t=1    M=1.000  name1=0.399  name2=0.413  name3=0.187
  t=5    M=1.000  name1=0.396  name2=0.408  name3=0.196
  t=20   M=1.000  name1=0.382  name2=0.390  name3=0.229
  t=100  M=1.000  name1=0.342  name2=0.343  name3=0.315
  t=400  M=1.000  name1=0.333  name2=0.333  name3=0.333
```

How to read it:

- **Defaults, `M` used up.** The meaning is gone by `t = 5` and everything
  stops. The words and codes of each name have reshuffled (not shown), but the
  name totals never moved.
- **Defaults, `M` held.** The reactions keep running for ever, but with
  `kappa1 = kappa2` each name keeps its 0.4, 0.4 and 0.2: no shared lexicon.
- **`kappa1 = 2`, `M` held.** Winner-take-all. At the start, name 2 has more
  of its molecules as codes (0.15 of 0.4) than as words, so it loses; name 3,
  all words, gains at first. By `t = 400` name 1 holds everything. This is the
  convergence to a single name that the naming game is about, and it needs
  code replacement to beat word replacement.
- **`kappa2 = 2`, `M` held.** The opposite: the three names are driven to a
  third each.

The other rates, `ks`, `kc` and `ka`, never move material between names. They
act through the split of each name between words and codes, which the formula
above shows sets how fast a name gains or loses, and so can change which name
wins; whether the outcome is one winner or a tie is set by the sign of
`kappa1 − kappa2`. To try more names, pass `N_s` (up to 20), extend `x0`, and
change the `(1, 2, 3)` in `report`. Every run here is fast.

## Results

**The published result.** The book reports one finding of De Beule, Hovig and
Benson (2011): turned into ODEs, the reactions show "the same lateral
inhibition behavior that had been reported in previous work", Steels (2000),
"converging to a shared lexicon". In Steels' agent-based experiments this is a
winner-take-all outcome: after a struggle between alternatives, one word
dominates. The book gives no figure, rates, initial state or numbers for the
chemical version, and the paper itself is behind a paywall and was not
available for this page, so the conditions under which the authors saw
convergence could not be checked.

**What Chemart reproduces.** Chemart's reading of the book's reactions shows
convergence to one name only in part of parameter space: with the meaning held
constant and `kappa1 > kappa2`, as in the run above. With the default equal
rates, each name's total is exactly conserved and no lexicon forms; with `M`
consumed as printed, the dynamics stop before any competition. So the network
*can* produce the published behaviour, but the defaults do not, and whether
the settings that do match the paper's is unknown. The likely sources of the
gap are the two gaps in the book: the lost indices in eqs. 16.5–16.6, which
Chemart had to reconstruct, and the missing rates and initial state.

Chemart has no chemistry-specific tests for this entry. Only the library-wide
contract tests run on it (the default network builds, the same seed gives the
same network, bad parameters are rejected). The runs on this page are not
checked by any test.

**Context in the book.** The book presents the model as an illustration of the
potential of artificial chemistries for language dynamics, noting that much
remains to be explored. It also mentions that De Beule, Hovig and Benson see
the combination of artificial chemistries with Steels' Fluid Construction
Grammar, a computational formalism for grammar, as a promising way to study
how grammatical codes form. The book describes no later work built on the
chemical naming game.

## Further reading

- Steels, L. (2000). Language as a complex adaptive system. In *Parallel
  Problem Solving from Nature PPSN VI*, Lecture Notes in Computer Science 1917,
  17–26. Springer. (Book ref [809]; open copy at
  <https://langev.com/pdf/steels00languageAs.pdf>.) Section 4
  describes the agent-based naming game and its lateral-inhibition rule.
- De Beule, J., Hovig, E. & Benson, M. (2011). Introducing dynamics into the
  field of biosemiotics.
  *Biosemiotics* 4(1), 5–24.
  <https://doi.org/10.1007/s12304-010-9101-1>
