## Introduction

The Chemical Abstract Machine (CHAM, or "cham" in the original paper) is a
way of describing how concurrent programs run, using a chemical solution as
the model of a computer's state. Gérard Berry and Gérard Boudol introduced it
at the POPL conference in 1990 (journal version 1992). Computer science has
*abstract machines* for sequential computation, such as the Turing machine or
the SECD machine for the λ-calculus. For concurrent programming, Berry and
Boudol wrote, "the situation is much less clear": the models it had were
either too weak (Petri nets, communicating automata) or were specification
languages rather than machines. The cham was their proposal for that missing machine.

The picture comes from the Gamma language of Banâtre and Le Métayer
([Gamma](gamma.md), the section just before this one in the book). The state
of a system is a solution in which molecules float about; a "magical
mechanism" stirs it, and molecules that meet may react according to rules.
Any number of reactions can happen at once, as long as they involve
different molecules, so parallelism is built in rather than added. What the
cham adds to Gamma, in the paper's own summary, is "the structure of molecules
as terms and the notion of a subsolution". Molecules are no longer plain
values but expressions of a formal language, and a whole solution can be
wrapped in a *membrane* and float, as a single molecule, inside a bigger one.

The language Berry and Boudol put into the solution is a *process calculus*:
a small algebra of programs that communicate over named channels. The best
known is Robin Milner's CCS (Calculus of Communicating Systems). In CCS the
program `a.P` waits to perform action `a` and then behaves as `P`; the
program `~a.Q` offers the complementary action (in print, a with a bar over
it); and when the two run side by side they synchronise and continue as `P`
and `Q`. In a cham, `a.P` and `~a.Q` are two *ions* floating in a solution, carrying
complementary valences much as ions carry opposite charges. When they meet,
the valences vanish and the bodies `P` and `Q` are released. Communication has become a chemical reaction, and the
non-determinism of concurrent programs (which of several partners a process
talks to) has become which molecule happens to collide with which.

The cham is a *formalism*, not a simulation model: there are no rates, no
concentrations and no time. Its contribution is conceptual. The paper shows
that the familiar process calculi can be written as chams, with the same
observable behaviour, and that the rules become simple local rewrites
instead of the inference rules normally used to define such languages.
Banzhaf and Yamamoto present it in their chapter on rewriting systems (book
§9.3), as a descendant of Gamma. It is the most "programming-language" member
of that chapter. [Gamma](gamma.md) rewrites bags of plain values;
[ARMS](arms.md) (§9.4) borrows the cham's heating and cooling vocabulary for
symbol strings; [P systems](p-systems.md) and [brane calculi](brane-calculi.md)
also nest membranes, but P systems rewrite symbols with maximal parallelism
and brane calculi make the membranes themselves act. [SOAS](soas.md) uses the
cham to specify self-assembling production lines.

## How it works

### Molecules, solutions and membranes

A cham has three ingredients: molecules, solutions and rules.

- **Molecules** are terms of an algebra. For the process-calculus chams they
  are CCS programs: `0` (the process that does nothing), prefixes such as
  `a.m` and `~a.m`, parallel composition `p|q`, restriction `m\a` (the channel
  `a` is private to `m`), and so on. In the examples below, `P` and `Q` stand
  for arbitrary processes; Chemart treats them as constants that never act.
- The full TCCS language used in the paper's §4 adds more operators, which
  appear in the formal specification below: relabelling `m[b/a]` (the
  channel `a` is renamed `b`), the internal sum `p (+) q` (the program picks
  one branch by itself), the external sum `p [] q` (the branch is picked by
  whichever side first communicates with the environment), and recursion
  `fix_X(X=...)`. The external sum is handled with pair molecules `<m, S>`,
  tagged `l:` or `r:` to remember which side offered an action. `tau.m`, a
  prefix with CCS's internal action τ, can be written, but only the paper's
  CCS machine uses it, and Chemart does not offer that machine.
- A **solution** is a multiset of molecules, written `{| m1, m2, ... |}`.
  Being a multiset, it has no order: `{| P, Q |}` and `{| Q, P |}` are the
  same solution. This is how the cham gets for free the rule that `p|q`
  behaves like `q|p`, which ordinary semantics has to state and prove.
- A solution is itself a molecule. Wrapping a solution in the **membrane**
  `{| ... |}` makes it one molecule of the enclosing solution, so solutions
  nest to any depth. Reactions inside a membrane happen locally.

A molecule that cannot be broken down further and carries a communication
capability, like `a.m`, is an **ion**, and the capability `a` is its
**valence**.

### Three kinds of rules

Each particular cham is given by a short list of *specific rules* of the
form `m1, ..., mk -> m1', ..., ml'`. Berry and Boudol sort them into three
classes:

- **Heating** rules break one molecule into simpler ones. The basic one is
  `p|q ⇀ p, q`: a parallel program floating in the solution falls apart into
  its two components, now separate molecules.
- **Cooling** rules are the inverse, rebuilding a compound molecule. Most
  heating rules are paired with a cooling rule, written together with a
  double arrow `⇌`. Heating and cooling are *structural*: they rearrange how
  a program is written, not what it can do.
- **Reaction** rules are irreversible and do the actual computing. The CCS
  reaction is `a.m, ~a.n -> m, n`: two complementary ions meet, the valences
  vanish, and the bodies are released.

To run a program `p`, put it alone in a solution, `{| p |}`, and heat it: the
parallel components fall apart until ions are exposed, and ions with opposite
valences react. The paper's first example is `a.b.0 | ~a.0 | ~b.0`. Heating
twice gives the three ions `a.b.0, ~a.0, ~b.0`; `a` reacts with `~a`, leaving
`b.0, 0, ~b.0`; then `b` reacts with `~b`, leaving `0, 0, 0`. A *cleanup*
rule, "0 evaporates when heated", removes the inert zeros and leaves the
empty solution.

### Four general laws

The specific rules say *what* can change. Four laws, common to every cham,
say *where* they may be applied (paper §3.1):

1. **Reaction law**: an instance of a rule's right-hand side can replace the
   matching instance of its left-hand side.
2. **Chemical law**: a reaction can happen in any solution regardless of the
   other molecules present. If `S -> S'`, then `S` together with `S''` becomes
   `S'` together with `S''`.
3. **Membrane law**: a subsolution inside a membrane evolves freely,
   whatever molecule it sits in.
4. **Airlock law**: `{| m |} ⊎ S ⇌ {| m ◁ S |}`, where `⊎` is multiset union.
   Any one molecule `m` can be pulled out of a solution and attached, as an
   **airlock** `m ◁ S`, to a membrane holding the rest. This is reversible.

Rules apply only to molecules floating in a solution, never inside a molecule
that is not a solution. Several rules may fire in parallel if no molecule is
used twice, but the paper notes that this does not add expressive power: a
non-conflicting parallel step is equivalent to any sequence of its single
steps. When more than one rule applies, the choice is non-deterministic.

### Why airlocks are needed

Restriction is where membranes earn their keep. In `(~a.P | Q)\b` the channel
`b` is private, but the inner process should still be able to offer its
action `~a` to the outside, since `a` is not `b`. If the body were a single
ion there would be a simple rule for this, the *restriction ion*
`(α.m)\a ⇌ α.(m\a)` provided α is neither `a` nor `~a`: the valence moves
out through the restriction. Here the body is compound, so first it is put
in its own local solution by the *restriction membrane* rule
`m\a ⇌ {| m |}\a`, where it can heat and react on its own. To get one of its
ions out, the airlock law isolates the ion next to the rest of the
subsolution, and a *heavy ion* rule `(α.m) ◁ S ⇌ α.(m ◁ S)` moves the valence
to the outside of the airlock, carrying the rest along. Berry and Boudol
built it this way so that every step can be undone: nothing is lost about
which subsolution the ion came from.

### A worked example: the default network

The default program in Chemart is the paper's communication example (§2.3),
`{| a.0 | (~a.P | Q)\b |}`. Chemart writes the co-name `ā` as `~a` and the
airlock `◁` as `<|`. Restriction binds tightest, so `~a.(P<|{|Q|})\b` means
`~a.((P ◁ {|Q|})\b)`. These are reactions from the default network, in the
order of the paper's derivation:

```
a.0|(~a.P|Q)\b -> a.0 + (~a.P|Q)\b                      heating  parallel
(~a.P|Q)\b -> {|~a.P|Q|}\b                              heating  restriction membrane
{|~a.P|Q|}\b -> {|Q,~a.P|}\b                            heating  parallel (inside the membrane)
{|Q,~a.P|}\b -> {|~a.P<|{|Q|}|}\b                       heating  airlock
{|~a.P<|{|Q|}|}\b -> {|~a.(P<|{|Q|})|}\b                heating  heavy ion
{|~a.(P<|{|Q|})|}\b -> (~a.(P<|{|Q|}))\b                cooling  restriction membrane
(~a.(P<|{|Q|}))\b -> ~a.(P<|{|Q|})\b                    heating  restriction ion
a.0 + ~a.(P<|{|Q|})\b -> 0 + (P<|{|Q|})\b               reaction reaction
0 -> ∅                                                  heating  inaction cleanup
(P<|{|Q|})\b -> {|P<|{|Q|}|}\b                          heating  restriction membrane
{|P<|{|Q|}|}\b -> {|P,Q|}\b                             cooling  airlock
```

Step by step: the top-level parallel falls apart into `a.0` and the
restricted process. The restricted process opens a membrane, and inside it
`~a.P | Q` falls apart into two molecules. The airlock pulls `~a.P` out next
to a membrane holding `Q`, and the heavy ion rule moves the valence `~a` to
the front. Now the membrane holds a single ion, so it is cooled away, and the
restriction ion rule moves `~a` through the restriction (allowed because `a`
is not `b`). The result, `~a.((P ◁ {|Q|})\b)`, is an ion floating in the top
solution next to `a.0`, and the two react. The released `0` evaporates, the
remainder opens its membrane again, and cooling the airlock puts `P` and `Q`
back together inside the restriction: `{| {|P, Q|}\b |}`. The ion `a.0` has
talked to a process hidden behind a private channel, and the private channel
is still private.

### What Chemart computes

Chemart does not run the cham as a simulation. It computes the *closure* of
the initial solution: every molecule that can appear floating in the top
solution, and every single-rule step between them. Each step is one
reaction of the network, labelled with its class (heating, cooling or
reaction), the rule that produced it, and how many membranes deep it
happened. The full structural relation is infinite (cooling can rebuild
`p|q`, `q|p`, `(p|q)|r` and so on without end), so the exploration follows a
fixed policy: heating, reactions, airlocks inside membranes, and the cooling
steps needed to undo them. The paper's three derivations all lie inside this
closure. A reaction such as `a.0 + X -> 0 + Y` reads: in any solution that
contains the molecules `a.0` and `X`, they can be replaced by `0` and `Y`
(the chemical law); the network lists molecules, not whole solutions. Since the cham has
no kinetics, reactions carry no rates.

## Using it

The default call above builds the communication example with the complete
TCCS machine: 18 molecules and 30 single steps, computed in milliseconds.
Every reversible step appears twice, once as heating and once as cooling,
which is why the first eight reactions come in pairs. The initial
state is `{'a.0|(~a.P|Q)\b': 1.0}`. Species ids are the molecules as text,
with no spaces and with the members of each solution sorted.

What each reaction is, and which membranes each molecule contains, are in
`net.extras`:

```python
net = chemart.generate_network("cham", seed=1)
for r, lab in zip(net.reactions, net.extras["reaction_rules"]):
    print(r.to_text(), lab)
# a.0 + ~a.(P<|{|Q|})\b -> 0 + (P<|{|Q|})\b
#     {'class': 'reaction', 'rule': 'reaction', 'reversible': False, 'depth': 0}
# {|P<|{|Q|}|}\b -> {|P,Q|}\b
#     {'class': 'cooling', 'rule': 'airlock', 'reversible': True, 'depth': 1}

net.extras["compartments"][r"{|P,Q|}\b"]
# [{'depth': 1, 'context': 'restriction \\b', 'molecules': ['P', 'Q']}]
```

`extras["machine"]` names the rule set, `extras["rules"]` lists its 18 rules
as text, and `extras["airlock_law"]` says whether the airlock law is on.

**The other published examples.** Set `program`:

| `program` | paper | molecules / steps |
|---|---|---|
| `heavy-ion-communication` (default) | §2.3, communication through a membrane | 18 / 30 |
| `ccs-minus-execution` | §2.2, the execution example above | 7 / 7 |
| `nondeterministic-choice` | §2.2, one ion and two possible partners | 6 / 3 |
| `external-sum` | §4.2, the external sum `[]` | 19 / 29 |

The execution example is small enough to print whole:

```python
net = chemart.generate_network("cham", program="ccs-minus-execution", machine="ccs-minus")
# a.b.0|~a.0|~b.0 -> a.b.0|~a.0 + ~b.0     heating  parallel
# a.b.0|~a.0 + ~b.0 -> a.b.0|~a.0|~b.0     cooling  parallel
# a.b.0|~a.0 -> a.b.0 + ~a.0               heating  parallel
# a.b.0 + ~a.0 -> a.b.0|~a.0               cooling  parallel
# a.b.0 + ~a.0 -> b.0 + 0                  reaction reaction
# 0 -> ∅                                   heating  inaction cleanup
# b.0 + ~b.0 -> 2 0                        reaction reaction
```

The non-deterministic example shows that the network holds *both* outcomes
as separate reactions, `a.0 + ~a.b.0 -> 0 + b.0` and
`a.0 + ~a.c.0 -> 0 + c.0`. Which one happens is not decided: there is no
reactor that picks.

**Machines.** `machine="ccs-minus"` is the small machine of the paper's §2
(parallel, reaction, restriction membrane and ion, heavy ion).
`machine="tccs"`, the default, is the complete machine of §4.3, which adds
relabelling, the two sums and recursion. On the default program both give
the same network; on `external-sum`, `ccs-minus` has no rule for `[]` and
stops after 3 molecules and 2 steps. `cleanup=False` drops the rules that
evaporate `0` and empty membranes.

**Switching off the airlock.** With `airlock=False` the default program stops
at 5 molecules and 6 steps, all structural. The restricted process opens its
membrane and falls apart inside it, but no ion ever gets out, so no reaction
happens. This is the paper's reason for introducing airlocks.

**Your own programs and rules.** `program="custom"` takes a `solution` in the
same notation, and `machine="custom"` takes `rules` as text, one per line,
with `<=>`, `=>`, `~>` and `->` for heating-with-cooling, heating, cooling and
reaction, and typed metavariables (`?p` an agent, `?m` any molecule, `?a` a
name). The parameter table below gives the full syntax; nothing is evaluated
as Python.

```python
net = chemart.generate_network("cham", program="custom", solution="a.0, ~a.0, ~a.0",
                               machine="custom", rules="reaction: ?a.?m, ~?a.?n -> ?m, ?n")
[r.to_text() for r in net.reactions]      # ['a.0 + ~a.0 -> 2 0']
```

Recursion can make the closure infinite. The recursive process
`fix_X(X=(a.0|X)\b)` never closes: with the default `max_species=500` it
stops at 500 molecules and 2,111 steps, `status="truncated"`, after about
1.7 seconds.

**Single steps.** `chemart.chemistries.cham.transitions(solution, rules)`
returns every one-step successor of a solution with its label, including the
cooling and airlock moves the closure leaves out. On
`{| a.0, ~a.b.0, ~a.c.0 |}` with the `ccs-minus` rules it gives 20
successors: the two reactions, 6 coolings that rebuild a parallel pair, and
12 airlocks.

## Results

The cham's results are in the 1990 paper, and they are statements about
programming languages rather than measurements. Chemart reproduces the worked
derivations; the theorems it does not attempt.

**CCS as chemistry (paper §2.2).** The execution of `a.b.0 | ~a.0 | ~b.0`
heats to three ions, reacts twice, and the cleanup rule leaves the empty
solution. The authors contrast this with the usual definition of CCS by
*structural operational semantics* (SOS), in which each step is a proof by
inference rules and simplifications such as `p|0 = p` must be proved
separately. In the cham, structure is handled once by heating, reactions are
chained while the solution stays hot, and a `0` simply evaporates. Chemart's
tests check every step of this derivation with the `ccs-minus` rules and
check that each is a reaction of the generated network.

**Non-determinism (§2.2).** In `{| a.0, ~a.b.0, ~a.c.0 |}` the ion `a.0` can
react with either partner, giving `{| b.0, ~a.c.0 |}` or `{| ~a.b.0, c.0 |}`
after cleanup. The tests check that exactly these two reactions are possible.

**Communication through a membrane (§2.3).** This is the example of the
worked derivation above, where membranes, airlocks and heavy ions let a
restricted process communicate with the outside while keeping every
structural step reversible. The tests check all eleven steps, the default
network's size (18 molecules, 30 steps) and that with airlocks switched off
no reaction can happen. The paper prints this derivation with two slips:
the second-last solution shows the membrane as `\c` instead of `\b`, and the
evaporation of the released `0` is not listed. Chemart follows the corrected
version (see the decisions).

**The full TCCS calculus and the external sum (§4).** TCCS is De Nicola and
Hennessy's variant of CCS, which separates internal steps from steps visible
to an observer. Its complete cham adds rules for relabelling, internal sum,
external sum and recursion; with the heavy ion and the three cleanup rules,
Chemart's `tccs` machine has 18 rules. The external sum `p [] q` (a choice
settled by the environment) needs a pair molecule `<m, S>` with a membrane
per summand and tags `l:` and `r:` recording which side offered the action;
the paper concludes that designing a cham "has much to do with standard
programming" and that the external sum "is not very natural in concurrent
abstract machines". In its example `{| a.~b.0 | ((~a.0 | b.0) [] Q) |}` the
left summand exposes `~a`, reacts, the left projection discards `Q`, and
cooling the airlock gives `{| ~b.0, 0, b.0 |}`. The tests check all nine
steps and each TCCS rule separately. The printed rule summary omits the heavy
ion rule, which this example needs, and the direction of several arrows
cannot be read in the scan; Chemart's `tccs` machine includes the heavy ion
rule and settles the arrows from the prose (see the decisions).

**Reversibility (§2.2, §3.2).** Most heating rules are paired with an
inverse cooling rule; the one-way rules are the reactions (including the
internal sum), the cleanup rules, the two projections and the unfolding of
recursion. The tests check that every reversible step along the three derivations can be undone, and that every
reversible reaction in the default network appears with its inverse.

**Equivalence with TCCS (§4.4).** The paper's main theorem says that, as
far as visible transitions are concerned, the solution `{| p |}` "can do
whatever the term p can do, and it cannot do more", and that the cham
"differs from the original TCCS calculus only in the number of internal
steps involved in computations". Chemart does not check this; it is a
proof, not a computation.

**What Chemart does not offer.** The CCS cham of §4.5, where an internal
communication produces a `τ`-ion that is only consumed when an outside
observer accepts it, is left out: that step is an observation, not a rule,
so the closure would stop at the first `τ`-ion. The paper found this
simulation of CCS "rather unsatisfactory", and showed weak bisimulation (the
same observable behaviour, ignoring internal steps) between `p` and `{| p |}`.
Also left out is the paper's §5, a new higher-order concurrent calculus
extracted from the cham, which the paper calls the γ-calculus (not the
γ-calculus of the Gamma family). It contains the lazy λ-calculus and can
define a "parallel or" that no λ-term can, but it has a different molecule
algebra. The chams of Inverardi and Wolf (1995), the book's other reference,
are not included either: that paper is closed access and was not used.

**Later work.** The book mentions several later systems that use or draw on the cham.
Suzuki and Tanaka, in the book's words, "demonstrated that it can be
utilized to model chemical systems", and defined [ARMS](arms.md) (§9.4), an
abstract rewriting system on multisets. MGS (§9.6, [mgs](mgs.md)) integrates
elements of the cham, cellular automata, P systems and L-systems in one rule
system. Self-organising assembly systems ([SOAS](soas.md), §20) use the cham
to specify rules by which industrial robots assemble themselves into a
production line. The book's chapter on networking also lists the cham, with
Gamma, MGS and P systems, as a relative of the Fraglets language
([fraglets](fraglets.md)). Inverardi and Wolf (1995) used cham specifications
to describe and analyse software architectures.

## Further reading

- Banâtre, J.-P. & Le Métayer, D. (1990). The Gamma model and its discipline
  of programming. *Science of Computer Programming* 15(1), 55–77 (book
  reference [58]). The language the cham grew out of.
- Milner, R. (1989). *Communication and Concurrency*. Prentice Hall. The
  standard text on CCS (book reference [589]).
- De Nicola, R. & Hennessy, M. (1987). CCS without τ's. *TAPSOFT '87*,
  Lecture Notes in Computer Science 249, 138–152. The TCCS calculus.
- Berry, G. & Boudol, G. (1989). The chemical abstract machine. INRIA
  research report RR-1133: <https://inria.hal.science/inria-00075426>
