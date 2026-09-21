## Introduction

Brane calculi are a family of formal languages, introduced by Luca Cardelli
(Microsoft Research) in a 2004 paper, for describing what biological membranes
*do*: engulf other membranes, fuse with them, pinch off vesicles, and pump
molecules across themselves. "Brane" is, in Cardelli's words, "a common
abbreviation for 'membrane' in physics". A brane calculus is a *process
calculus*, the kind of notation computer scientists use to describe concurrent
programs (the π-calculus, Mobile Ambients), adapted so that its terms are
nested membranes and its steps are membrane rearrangements.

The central idea is that computation happens *on* the membrane, not inside it.
Cardelli starts from the biology: cell membranes are "not just containers:
they are coordinators and active sites of major activity", with proteins
embedded in them that act on both sides at once. So in a brane calculus each
membrane carries a small program of actions, and a reaction is what happens
when an action on one membrane meets a matching action on another. A virus
membrane that says "be engulfed" meets a cell membrane that says "engulf", and
the virus ends up inside a new vesicle inside the cell.

The worked example that runs through the paper is the life cycle of the
Semliki Forest virus: the cell swallows the virus, a cellular compartment (the
endosome) fuses with the swallowed vesicle, the virus exploits this to release
its core into the cell, the cell's machinery copies the viral RNA and builds
new virus parts, and new viruses bud out through the cell membrane. Each of
these stages becomes a short sequence of reactions. The paper also shows
molecular pumps and ion channels, and small examples in which a molecule
released by one membrane triggers another membrane to engulf it.

It is a **formalism**: a notation with precise rules, not a simulation model.
It has no rates or concentrations; it says which rearrangements are possible,
not how fast they happen. Banzhaf and Yamamoto mention it in one paragraph of
§9.7, "Other Formal Calculi Inspired by a Chemical Metaphor", next to the
[Kappa calculus](kappa-calculus.md). The comparison they draw is with
[P systems](p-systems.md): both work on a tree of nested membranes, but P
systems compute with objects *inside* the compartments, while brane calculi
"emphasize operations on membranes and performed by membranes". In a brane
calculus the membranes themselves are what change: they engulf, fuse and
split, while the molecules between them only cross membranes when a membrane
action moves them.

## How it works

### Systems, membranes and actions

A configuration is written as text. `s[P]` is a membrane whose surface is `s`
and whose contents are `P`. Contents are a *system*: membranes and free
molecules, separated by commas, in any order. A surface is a *brane*: a set of
actions separated by `|`, also in any order. Two operators complete the
language:

- `a.s` means "perform action `a`, then the surface becomes `s`". An action is
  used up when it fires.
- `!` means "an unlimited supply": `!a` is an action that can fire any number
  of times, and `!m` is an unlimited supply of molecule `m`. A reaction uses
  one copy taken from the supply and leaves the supply itself in place.

`0` stands for nothing: the empty system or the empty surface. `[]` alone is
an empty membrane with contents omitted.

Because both commas and bars are unordered, `a|b[P, Q]` and `b|a[Q, P]` are the
same configuration. Cardelli calls this *structural congruence*: it captures
that membranes are a two-dimensional fluid (proteins drift freely on them)
floating in a three-dimensional fluid (the contents drift freely too).

Chemart writes Cardelli's notation in ASCII. Actions come in pairs, an action
and its *co-action* (the paper marks it with ⊥, Chemart with the prefix `co`).
An optional subscript, `phago_n` with `cophago_n`, says which action pairs with
which; actions without a subscript match each other.

### The reactions

There are six membrane reactions. The first three are Cardelli's basic set:

- **Phago** (phagocytosis, "eating"): a membrane with `phago` next to a
  membrane with `cophago(r)` is engulfed by it and wrapped in a new membrane
  `r`, supplied by the co-action.
- **Exo** (exocytosis): a membrane with `exo` inside a membrane with `coexo`
  fuses with it; its contents are released *outside* the outer membrane, and
  the two surfaces merge.
- **Pino** (pinocytosis, "drinking"): a membrane with `pino(r)` creates an
  empty vesicle `r[]` inside itself.

The other three are their mirror images, based on fusion and fission:

- **Mate**: two sibling membranes with `mate` and `comate` fuse into one,
  with both contents together.
- **Bud**: a membrane with `bud` inside a membrane with `cobud(r)` leaves
  through it, wrapped in a new membrane `r`.
- **Drip**: a membrane with `drip(r)` releases an empty vesicle `r[]` to its
  outside.

A seventh reaction moves molecules. A **bind&release** action
`p1(p2)=>q1(q2)` takes the molecules `p1` from outside the membrane and `p2`
from inside it, all at once, and puts `q1` outside and `q2` inside. Empty parts
can be left out: `ATP=>ADP,Pi(H+,H+)` takes one ATP from outside, puts ADP and
Pi (phosphate) back outside and two protons inside. This is how pumps,
channels and catalysts are written. Conservation of mass is not built in: the
modeller must write it into the action.

The formal rules under *Formal specification* below write each reaction
with variables. In `phago_n.s|s0[P]`, `phago_n` is the action that fires, `s`
is what the surface becomes after it, `s0` is the rest of the surface (which
is carried along unchanged) and `P` is the contents. Phago, for instance, reads
`phago_n.s|s0[P], cophago_n(r).t|t0[Q] -> t|t0[r[s|s0[P]], Q]`: the eaten
membrane `s|s0[P]` ends up inside the eater, wrapped in `r`.

Any other name on a membrane is an inert placeholder; any other name in a
system is a molecule. Cardelli uses such names (`X`, `Z`, `Nucleus`) for the
parts of a cell that a model leaves unspecified.

### Bitonality

Six membrane reactions share a property Cardelli calls *bitonality*. Colour
the space inside a membrane light and the space between it and the next one
dark, alternating as you go deeper. Phago, Exo, Pino, Mate, Bud and Drip never
move anything onto a background of the other colour: whatever is inside a cell
stays on an inside-coloured background, and outside material can come in only
wrapped in another membrane. "Bitonality is common in cellular-scale living
systems," Cardelli writes, and although not universal it "inspires a
collection of basic reactions that are biologically implementable". Bind&release
is the exception: it carries molecules across a membrane.

### A worked example: the virus gets in

The default run is the infection stage (Cardelli's Figure 7). Its starting
configuration, as printed by Chemart's `layout`, is a virus next to a cell:

```
membrane !coexo|!cophago(mate)
  membrane !coexo|!comate
  Z
membrane phago.exo
  membrane X|!bud
    vRNA
```

The second top-level membrane is the virus: an envelope with the surface
`phago.exo`, around the nucleocapsid (the capsid `X|!bud` holding the viral RNA
`vRNA`). The first is the cell: its outer membrane can engulf (`!cophago(mate)`)
and fuse (`!coexo`) without limit, and it contains an endosome
(`!coexo|!comate`) and the unspecified rest of the cell, `Z`. Chemart finds
exactly three reactions:

```
!coexo|!cophago(mate)[!coexo|!comate[],Z],phago.exo[X|!bud[vRNA]]
  -> !coexo|!cophago(mate)[!coexo|!comate[],Z,mate[exo[X|!bud[vRNA]]]]      phago
  -> !coexo|!cophago(mate)[!coexo|!comate[exo[X|!bud[vRNA]]],Z]             mate
  -> !coexo|!cophago(mate)[!coexo|!comate[],X|!bud[vRNA],Z]                  exo
```

1. **Phago.** The virus's `phago` meets the cell's `cophago(mate)`. The virus
   is now inside the cell, wrapped in a new vesicle whose surface is `mate`;
   its own envelope has used up `phago` and now reads `exo`.
2. **Mate.** The vesicle's `mate` meets the endosome's `comate`, and the two
   fuse. The virus is now inside the endosome.
3. **Exo.** The virus envelope's `exo` meets the endosome's `coexo`. The
   envelope fuses with the endosome membrane and releases its contents on the
   far side, which is the inside of the cell. The nucleocapsid `X|!bud[vRNA]`
   is now free in the cell, next to an empty endosome.

This is the trick Cardelli describes: the virus uses the cell's own transport
route and a fusion step to put its core where the cell's contents are. Nothing
more can happen, so the last configuration is terminal.

### What Chemart turns into a network

The calculus is not a chemistry of many interacting molecules, so Chemart maps
it onto a reaction network in one particular way. A **species is a whole
configuration**: the entire nest of membranes, written in a canonical text form
in which the unordered parts are sorted and `!` is simplified, so that two
congruent configurations get the same name. Each **reaction is one step**,
from one configuration to the next, labelled with the rule that fired. Named
definitions such as `membrane` are expanded in the species names, which is why
they are long.

The reactor, the part of an artificial chemistry that decides which reactions
happen, does not simulate anything. Starting from the initial configuration,
Chemart lists every configuration reachable by any sequence of steps, and
every step between them. This is the *reduction graph*. A reaction can fire
anywhere, at any depth of nesting. When an unlimited supply is involved (`!`),
the graph can be infinite, and `max_species` caps the number of configurations
explored.

## Using it

The default call reproduces Cardelli's Figure 7, the infection stage walked
through above. Its 4 species are the 4 configurations and its 3 reactions the
3 steps. Besides the network, `net.extras` holds:

```python
net.extras["reaction_rules"]   # ['phago', 'mate', 'exo']  rule of each reaction, in order
net.extras["main"]             # canonical text of the starting configuration
net.extras["program"]          # the definitions it was built from
a = net.extras["analysis"]
a["infection (fig. 7)"]        # {'configuration': '!coexo|!cophago(mate)[...]', 'reached': True}
a["terminal"]                  # configurations with no way out (only for a complete closure)
```

`analysis` names each configuration the paper derives and says whether the
closure reaches it. `species[i].structure` gives the indented membrane tree
of a configuration, as shown above.

**The published systems.** `system` picks one of Cardelli's examples; the
table lists them all. Each closes in a fraction of a second:

```python
for system in ["viral-reproduction", "eat-me", "seek-and-store", "plant-vacuole"]:
    net = chemart.generate_network("brane-calculi", system=system)
    print(system, len(net.species), net.extras["reaction_rules"])
```

```
viral-reproduction 4 ['exo', 'bud', 'phago']
eat-me 4 ['bind&release', 'bind&release', 'phago']
seek-and-store 5 ['bind&release', 'pino', 'bind&release', 'mate']
plant-vacuole 5 ['bind&release', 'bind&release', 'bind&release', 'bind&release', 'bind&release']
```

**Viral replication** has no end. The RNA copies itself without limit, so the
reduction graph is infinite and is cut off at `max_species`. The analysis
reports the largest numbers of free nucleocapsids and envelope vesicles in any
configuration explored. They grow slowly with the budget:

| `max_species` | reactions | max nucaps | max envelope vesicles | time |
|---|---|---|---|---|
| 1,000 (default) | 2,795 | 2 | 4 | 0.5 s |
| 5,000 | 16,567 | 3 | 6 | 3 s |
| 20,000 | 74,678 | 4 | 8 | 12 s |

Larger budgets are slow and add little: the graph is explored breadth first,
so most of the budget goes on the many orders in which the same few copies can
be made.

**Your own system.** With `system="custom"`, `program` takes definitions
`name := term`, one per line, with `#` comments; `main` is the starting
configuration. The full syntax is in the module docstring of
`chemart/chemistries/brane_calculi.py`. For example, Cardelli's encoding of
Mate by one Phago and two Exo:

```python
prog = """
# Mate encoded by Phago and two Exo (Cardelli 2004, sec. 3.2)
left  := phago_n.exo_m.s[P]
right := cophago_n(coexo_m.exo_k).coexo_k.t[Q]
main  := left, right
"""
net = chemart.generate_network("brane-calculi", system="custom", program=prog)
for r, rule in zip(net.reactions, net.extras["reaction_rules"]):
    print(rule, ":", r.reactants, "->", r.products)
```

```
phago : {'cophago_n(coexo_m.exo_k).coexo_k.t[Q],phago_n.exo_m.s[P]': 1} -> {'coexo_k.t[Q,coexo_m.exo_k[exo_m.s[P]]]': 1}
exo : {'coexo_k.t[Q,coexo_m.exo_k[exo_m.s[P]]]': 1} -> {'coexo_k.t[P,Q,exo_k|s[]]': 1}
exo : {'coexo_k.t[P,Q,exo_k|s[]]': 1} -> {'s|t[P,Q]': 1}
```

The end point `s|t[P,Q]` is what one Mate step gives: one membrane with both
surfaces and both contents. The module also exports its parser and one-step
function (`parse_system`, `parse_program`, `steps`, `reachable`, `layout`,
`membranes`, `molecule_parities`) for working with terms directly.

## Results

Cardelli's paper is a proposal of a notation, and its results are
derivations: sequences of reactions showing that the notation can express a
biological process. It contains no simulations and no numbers. Chemart's tests
replay each derivation step by step and check that each configuration has
exactly the one successor the paper gives.

**Viral infection (Figure 7).** Virus and cell reduce by Phago, Mate and Exo to
a cell with the nucleocapsid free inside it, `membrane[nucap, cytosol]` in the
paper's shorthand. Chemart reproduces this exactly: the three steps are the
only ones possible, and the paper's end configuration is the only terminal one.

**Viral reproduction (Figure 8).** A vesicle lined with the viral envelope
protein fuses with the cell membrane by Exo, which puts the envelope protein
on the cell surface; the nucleocapsid then buds out through that patch,
wrapped in a new envelope, as a new virus. Chemart reproduces both steps. It
also finds a third step the paper does not discuss: the same cell membrane can
engulf the new virus again. The virus is then stuck in a vesicle, because the
example's unspecified cell contents `Z'` hold no endosome to fuse with.

**Nucleocapsid replication (§4.6).** The middle of the life cycle. A trigger
molecule in the cell makes the capsid push its RNA out. The RNA is copied; a
fictitious translation membrane turns it into capsid proteins that assemble,
by a Drip, into an empty capsid that then takes up an RNA, which is a new
nucleocapsid; and the endoplasmic reticulum, the cell's membrane factory,
turns it into envelope vesicles. Cardelli's claim is that
`nucap ∘ cytosol` reduces to configurations with any numbers of nucleocapsids,
envelope vesicles and spare RNAs. Chemart's tests check each of the three
paths separately, and check that the truncated graph contains configurations
with at least two nucleocapsids and with nucleocapsids, vesicles and capsid
residue together. The paper writes RNA copying as a plain chemical reaction;
Chemart uses the paper's own §4.2 recipe for a reaction, an empty membrane
whose bind&release action performs it.

**Mate, Bud and Drip reduce to Phago, Exo and Pino (§3.2).** Cardelli shows
that each fusion or fission reaction can be done by three of the basic ones,
calling this "just a test of expressive power"; in practice all six should be
primitives. Chemart offers all six directly, and its tests run the three
encodings and check that each lands where the direct reaction does. Two of the
paper's printed derivations contain typographical errors, corrected as the
decisions describe.

**Bitonality.** Cardelli notes that the six membrane reactions preserve the
nesting parity of every subsystem. Chemart's tests check this on every
membrane reaction of the viral-infection and viral-reproduction networks and
of a custom Pino network: the nesting parity of every molecule is unchanged,
and the number of membranes changes as the rule requires (up by one for Phago,
Pino, Bud and Drip, down by one for Exo and Mate). The other published systems
are not checked this way.

**Pumps, channels and molecular triggers (§4.4, §4.5).** The plant vacuole is
a storage compartment in plant cells. Its membrane carries three replicated
bind&release actions: a proton pump that uses ATP to push two protons in, a
channel that lets chloride in, and an antiporter that swaps a sodium ion
outside for a proton inside. The paper gives the membrane but no
surroundings; Chemart places one ATP, one Cl⁻ and one Na⁺ outside. The tests
check that pumping must come first, that channel and antiporter can then fire
in either order, and that the single end state holds H⁺, Cl⁻ and Na⁺ inside
with ADP, Pi and one H⁺ outside. In *Eat Me*, a membrane releases a molecule
that makes a neighbour engulf it; in *Seek and Store*, a cell recognises a
nutrient, takes it in and stores it in an internal vesicle. Both derivations
are reproduced step by step.

**Not offered.** Cardelli's paper also sketches protein complexes (§4.7) and
extensions (§5): communication between membranes, choice, restriction of
names, and "atonal" In/Out transport that breaks bitonality. Chemart
implements none of these. Nor does it add rates; the paper only suggests that
rates in the style of the stochastic π-calculus "should yield quantitative
modeling".

**Later work.** Busi and Gorrieri (2006) compared the two basic calculi as
computing devices. They showed that a fragment of Phago/Exo/Pino with only
Phago and Exo is Turing powerful, by encoding a Random Access Machine, while
for Mate/Bud/Drip it is decidable whether every computation terminates, so no
faithful such encoding exists there. As they explain, in Mate/Bud/Drip the
maximum nesting depth of membranes cannot grow during a computation, while in
Phago/Exo/Pino it can. These names, "Phago/Exo/Pino" and "Mate/Bud/Drip", are how
the two calculi are known in that literature. The book also cites Fellermann
et al. [269] for applying brane calculi to the transport of molecular cargo
by synthetic containers; that reference is a two-page abstract with no open
copy, and its follow-up (Fellermann and Cardelli, 2014) uses a different
formalism, which Chemart does not implement.

## Further reading

- Busi, N. & Gorrieri, R. (2006). On the computational power of Brane Calculi.
  *Transactions on Computational Systems Biology VI*, LNCS 4220, pp. 16–43,
  Springer. Authors' copy: <http://www.cs.unibo.it/~gorrieri/Papers/tcsm06.pdf>
- Fellermann, H. & Cardelli, L. (2014). Programming chemistry in DNA-addressable
  bioreactors. *Journal of the Royal Society Interface* 11(99): 20130987.
