## Introduction

Michael Conrad, one of the pioneers of molecular computing, argued from the
1980s on that molecules compute in a way digital machines do not. An enzyme
recognises its target by *shape*: the target fits a pocket of the enzyme the
way a key fits a lock, and the better the fit, the tighter they bind. Conrad
saw this lock-and-key binding as a form of pattern recognition, a task that
conventional computers find hard, done in a single physical step. He described
such a computation as a system "falling" downhill in free energy towards its
answer (book §11.2.1, citing Conrad 1992). This entry is Chemart's version of
that idea, and of the laboratory experiment that later put it to work.

The idea comes in two forms. The first is a design Conrad sketched, the
*self-assembly processor* of the book's Figure 11.10: each input signal
releases a molecule of a particular shape, the molecules assemble into a
complex, and an adaptor molecule that recognises the shape of that complex
starts a signalling cascade that amplifies it into an output. That processor
was never built. The second is a working experiment, *Enzymatic Computing* by
Zauner and Conrad (2001) at Wayne State University. They used one real
enzyme, malate dehydrogenase (MDH), a standard enzyme of cell metabolism. The
inputs were two bits, each coded as the presence or absence of salt (magnesium
or calcium ions) in the test tube. The output was how much product the enzyme
made in a fixed time. With the right amounts of salt, the single enzyme
computed the exclusive-or (XOR): output 1 when exactly one input is 1.

The XOR matters because it is the classic function that a single threshold
element cannot compute. A threshold element adds up its inputs and fires above
some level; the artificial neurons of a one-layer perceptron and, as the paper
puts it, transistors work this way. Minsky and Papert showed that no such
element can do the XOR. The enzyme could, because its activity first rises and
then falls as the salt increases: one dose of salt speeds it up, a double dose
slows it down again. The authors' conclusion is that an enzyme is more than a
summing element: it is a context-sensitive pattern recogniser that responds
both to the molecules it binds and to its chemical surroundings.

The book places Conrad's computer in its chapter on bio-inspired chemistries,
under "Lock-and-Key Artificial Chemistries" (§11.2.1), as a *wet* artificial
chemistry, meaning one made of real molecules. It was not proposed as an
artificial chemistry; the book reads it as one in hindsight, next to
[Typogenetics](typogenetics.md) and the [bitstring immune system
model](farmer-immune.md), as an origin of the lock-and-key matching that
string chemistries like [MCS.bl](mcs-bl.md), [Stringmol](stringmol.md) and
[SAC](sac.md) use. Unlike those, nothing here is a string that evolves: the
chemistry is one enzyme, a handful of small molecules, and a fixed recipe. It
also differs from the catalog's other wet entries, [Adleman's DNA
computation](dna-hpp.md) and the [DNA automaton](dna-automaton.md), which
compute with DNA sequences rather than with an enzyme's response to its
milieu. The book returns to Conrad in §17.1.2 for a principle he stated
alongside this work, the trade-off between programmability, evolvability and
efficiency, described under Results.

## How it works

### The experiment Chemart models

MDH carries out one reaction:

```
L-malate + NAD+ -> oxaloacetate + NADH + H+
```

L-malate is the *substrate*, the molecule the enzyme converts. NAD+ is a
*coenzyme*, a helper molecule that takes up the hydrogen removed from malate
and leaves as NADH. NADH absorbs ultraviolet light at 339 nm and NAD+ does
not, so a spectrophotometer (an instrument that measures how much light a
sample absorbs) can watch NADH build up in the cuvette, the small transparent
test tube. The *response* is the NADH absorbance read a fixed time after the
enzyme is added.

Magnesium (Mg2+) and calcium (Ca2+) ions take no part in this reaction. They
are the *milieu*: they change how active the enzyme is. Zauner and Conrad
attribute such effects to the enzyme's conformational dynamics, the way its
shape (*conformation*) shifts in response to what surrounds it, but they did
not model this effect
("our purpose here was not to investigate mechanism"); they measured it.

### Coding two bits in salt

Each input bit arrives as a portion of solution. A 1-signal contains a fixed
amount of salt, the *encoding*; a 0-signal contains none. With two input
lines, the cuvette receives zero, one or two portions of salt, so the enzyme
sees one of three *milieu states*:

- **a**: inputs `00`, no salt;
- **b**: inputs `01` or `10`, one portion;
- **c**: inputs `11`, two portions.

The enzyme cannot tell `01` from `10`, so only operations that treat the two
inputs alike (commutative ones) can be computed. Writing `r(a)`, `r(b)` and
`r(c)` for the enzyme's response in each state, an operation is a choice of
which states must give a high output. The XOR needs state b high and states a
and c low. Its *signal strength* is the margin between the lowest response
that must be high and the highest response that must be low:

```
Δs(XOR) = r(b) − max(r(a), r(c))
```

If `Δs` is positive, a single threshold placed in the gap separates the
patterns, and the operation works with that encoding. A response that only
rises (or only falls) with salt puts state b between a and c, so `Δs(XOR)` can
never be positive. A response that peaks at one portion and drops at two can
make it positive. That is the whole trick.

### The shapes: Chemart's lock-and-key layer

The book describes Conrad's computer through lock-and-key recognition, but the
2001 paper gives no shapes. Chemart therefore adds its own encoding (see the
implementation decisions below): every small molecule carries a 6-bit *key*,
and the enzyme carries three 6-bit *locks*, one for each binding site: the
coenzyme site, the substrate site and the *modulator* site where the ions sit.
A key fits a lock where its bits are the inverse of the lock's bits. A
*ligand*, any small molecule that binds the enzyme, binds a site when at least `match_bits` of the six positions are complementary
(five by default).

| molecule | key | site | lock | complementary positions |
|---|---|---|---|---|
| NAD+ | `110100` | coenzyme | `001011` | 6 of 6 |
| L-malate | `101010` | substrate | `010101` | 6 of 6 |
| Mg2+ | `111000` | modulator | `000111` | 6 of 6 |
| Ca2+ | `011000` | modulator | `000111` | 5 of 6 |
| NADH | `110111` | coenzyme | `001011` | 4 of 6 |
| oxaloacetate | `011010` | substrate | `010101` | 4 of 6 |

These scores come from `net.extras["interaction_law"]["match_scores"]`. The
choice is made so that the right molecules bind, the products do not rebind,
and calcium is a slightly worse fit than magnesium, mirroring the paper's
observation that magnesium gives the stronger signal. Nothing in the paper
says the ions' effect is a matter of shape fit; it is Chemart's way of
expressing the book's lock-and-key criterion.

### A worked example: the default network

The default network is the tabletop device with magnesium as the only signal
carrier and both inputs set to 1 (milieu state c). Its twelve reactions, from
`net.to_text()`, are one binding step for the ion and one catalytic cycle for
each form of the enzyme:

```
MDH + Mg -> MDH_Mg
MDH_Mg -> MDH + Mg

MDH_Mg + NAD -> MDH_Mg_NAD
MDH_Mg_NAD -> MDH_Mg + NAD
MDH_Mg_NAD + MAL -> MDH_Mg_NAD_MAL
MDH_Mg_NAD_MAL -> MDH_Mg_NAD + MAL
MDH_Mg_NAD_MAL -> MDH_Mg + OAA + NADH + H
```

and the same cycle for the free enzyme `MDH`. Species names list what is bound
to the enzyme: `MDH_Mg_NAD_MAL` is the magnesium-loaded enzyme holding NAD+
and malate. Step by step:

1. A magnesium ion fits the modulator lock (6 of 6) and binds reversibly. The
   enzyme is now in a different *conformer*, `MDH_Mg`.
2. NAD+ fits the coenzyme lock and binds, then malate fits the substrate lock
   and binds. Both steps are reversible: the book stresses that lock-and-key
   bindings are typically transient.
3. The full *ternary* complex (enzyme, NAD+ and malate together) turns
   over, that is, completes the reaction, releasing oxaloacetate (`OAA`), NADH and a
   proton (`H`). This is the one irreversible step, and the only place NADH,
   the readout, is made. The proton is taken up by the buffer.
4. NADH and oxaloacetate fit their old locks in only 4 of 6 positions, below
   the threshold of 5, so they do not bind again.

The two conformers run the same cycle. What should distinguish them is *how
fast* they run it, and that is exactly what the experiment measured and never
turned into rate constants. Every reaction in the network therefore has no
rate. The network records which complexes can form and which quantities are
conserved (the enzyme, the NAD+/NADH pool, the malate/oxaloacetate pool, and
each ion, free plus bound), not how much NADH appears.

### Where the computation lives

Because the response levels are not modelled, Chemart computes the
classification separately, from Table 1 of the paper. For the chosen
operation it lists the truth table, the milieu states that must be active and
the signal-strength formula. Functions in the module apply those formulas to
any response levels you supply.

The reactor is a single batch assay: everything is mixed in one cuvette and
read once. The formal specification below summarises the molecules (shapes and
enzyme states), the reactions (binding, release and turnover) and this reactor.

## Using it

The default run is the magnesium-only device of the paper's Figures 8 and 9,
set to the XOR, with input `11`. The concentrations are those of the published
recipe after mixing: two 0.8 mL signal portions and 0.5 mL of enzyme solution
make 2.1 mL, so one 1-signal of 190 mM MgCl2 gives 72.38 mM Mg2+ in the
cuvette, and `11` gives twice that.

```python
net = chemart.generate_network("conrad-enzymatic", seed=1)
net.initial_state
# {'NAD': 1.29, 'MAL': 5.41, 'Mg': 144.76}
a = net.extras["analysis"]
a["truth_table"]            # {'00': 0, '01': 1, '10': 1, '11': 0}
a["milieu_state"]           # 'c'
a["milieu_mM"]              # {'Mg': {'a': 0.0, 'b': 72.38, 'c': 144.76}}
a["signal_strength"]        # 'r(b) - max(r(a), r(c))'
a["linearly_separable"]     # False
```

The enzyme itself is absent from `initial_state`: the paper doses it as a
volume of commercial suspension, so its molar concentration is unknown. The
full protocol is in `net.extras["experiment"]`, and the published device
result (135 presentations, all correct) in
`net.extras["analysis"]["published_classification"]`. The seed has no effect:
the network is deterministic.

**Changing the input.** `inputs` only changes how much salt is in the cuvette;
the reactions stay the same. `inputs="00"` sets `Mg` to 0, `"01"` and `"10"`
to 72.38.

**Magnesium and calcium together.** The paper's example of a mixed encoding
is 20 mM Mg2+ plus 40 mM Ca2+ per 1-signal, which realises the XOR. With both
ions there are four conformers (free, Mg-bound, Ca-bound, both), and the
network grows to 19 species and 28 reactions:

```python
net = chemart.generate_network("conrad-enzymatic", seed=1,
                               encoding_mM={"Mg": 20.0, "Ca": 40.0}, inputs="01")
net.summary()          # conrad-enzymatic: 19 species, 28 reactions, status=complete
net.initial_state      # {'NAD': 1.29, 'MAL': 5.41, 'Mg': 20.0, 'Ca': 40.0}
```

**Stricter recognition.** `match_bits=6` admits only perfect complements, so
calcium no longer binds and any encoding that uses it is refused:

```
ValueError: encoding_mM: Ca does not recognise the modulator site at match_bits=6
(its key is complementary in 5 of 6 positions); lower match_bits or drop Ca
```

**The other operations.** `operation` selects one of the six commutative
two-bit operations. It changes only the analysis, not the network:

```
AND   {'00': 0, '01': 0, '10': 0, '11': 1}  active ['c']       r(c) - max(r(a), r(b))
OR    {'00': 0, '01': 1, '10': 1, '11': 1}  active ['b', 'c']  min(r(b), r(c)) - r(a)
XOR   {'00': 0, '01': 1, '10': 1, '11': 0}  active ['b']       r(b) - max(r(a), r(c))
NAND  {'00': 1, '01': 1, '10': 1, '11': 0}  active ['a', 'b']  min(r(a), r(b)) - r(c)
NOR   {'00': 1, '01': 0, '10': 0, '11': 0}  active ['a']       r(a) - max(r(b), r(c))
NXOR  {'00': 1, '01': 0, '10': 0, '11': 1}  active ['a', 'c']  min(r(a), r(c)) - r(b)
```

**Testing a response of your own.** The classification functions take
response levels for the three milieu states. The levels below are made up, a
response that peaks at one portion of salt; they are not measurements:

```python
from chemart.chemistries.conrad_enzymatic import signal_strength, implementable, classify
r = {"a": 0.2, "b": 0.9, "c": 0.4}
signal_strength("XOR", r), implementable("XOR", r)   # (0.5, True)
signal_strength("OR", r),  implementable("OR", r)    # (0.2, True)
signal_strength("AND", r), implementable("AND", r)   # (-0.5, False)
classify("XOR", r, threshold=0.6)   # {'00': 0, '01': 1, '10': 1, '11': 0}
```

With this response, both the OR and the XOR work, and the threshold decides
which one you get: the same situation the paper describes for its encodings.
Every call is instantaneous.

## Results

### The enzyme's response surface

Zauner and Conrad first mapped how MDH responds to the two ions. They ran 36
assays, six levels of MgCl2 against six of CaCl2, up to about 133 mM of each
ion in the 3 mL cuvette, and read the NADH absorbance 300 s after starting the
reaction. Interpolating between the points gave a *response surface* (Figure
1). Its key feature is qualitative: along both the magnesium and the calcium
axis, activity rises to a maximum and then falls. The authors call this convex,
or strictly nonmonotonic, and report that the shape holds over much of the
reaction course after 120 s. They note that ionic strength was already known
to suppress mitochondrial MDH and that magnesium had been reported to activate
the reaction in blood serum, and suggest the peak comes from these two effects opposing each other.

Chemart does not reproduce the surface. It was published as a picture, not as
numbers, so no absorbance appears anywhere in the entry.

### Which logic operations one enzyme can do

The paper's Table 1 turns the surface into logic. Of the 16 operations on two
bits, 6 ignore one or both inputs and 4 are not commutative, which leaves six:
AND, OR, XOR, NAND, NOR and NXOR. For each, the table lists the milieu states
that must be active and the signal-strength formula. Scanning the surface with
these formulas, the authors found encodings with positive signal strength for
AND, OR and XOR (Figure 5) and for NAND and NOR (Figure 6). NXOR was not
implementable with the convention that a high response is an active output.
Two encodings are named: 20 mM Mg2+ per 1-signal gives an OR, and adding 40 mM
Ca2+ to the 1-signal turns it into an XOR without moving the threshold. Figure
5C also shows that both ions work for the XOR alone or together, and that the
signal is stronger with magnesium alone than with calcium alone.

Chemart's tests check that its truth tables and signal-strength formulas are
exactly those of Table 1, on random response levels. To check the "five of
six" result, the tests use an invented response shaped like the published one
(a rising stimulation times a falling suppression) at many encodings: exactly
AND, OR, XOR, NAND and NOR appear, NXOR never does. The reason NXOR fails is
Chemart's reading, not a statement of the paper: NXOR needs the middle state b
to be *lower* than both a and c, and a response that rises then falls never
dips in the middle. The stronger magnesium signal is represented only by
calcium's weaker shape fit, which is an encoding choice, and the tests check
the scores.

### Turning an inseparable problem into a separable one

Figure 7 of the paper states the central argument. With a response that only
rises, a high threshold gives the AND and a low one the OR, but the XOR needs
two thresholds, fire above the lower one and not above the higher. An element
whose response is nonmonotonic removes the need for the second threshold: it
groups the inputs so that one threshold separates them. MDH has the required
property, so it transforms the linearly inseparable XOR into a separable
problem. Chemart's tests reproduce this argument on the formulas: no rising or
falling response makes the XOR (or NXOR) implementable, every such response
makes AND/OR (or NOR/NAND) implementable, and a response peaking at one
portion of salt makes the XOR work.

### The tabletop device

To show the XOR working "in a device context", the authors built a small
apparatus (Figure 8): a flow cuvette in a spectrophotometer, hand-operated
syringe pumps for the 0-signal, 1-signal and enzyme solutions, and a flushing
line to clean the cuvette between runs. Injecting the enzyme starts a timer,
and the computer reads the absorbance at a fixed time as the classification.
Using MgCl2 as the signalling substance, they presented 135 patterns in mixed
order: 45 of `00`, 46 of `01`/`10` and 44 of `11`. With a reading 10 s after
the start, one threshold separated all 135 correctly, "though in a few
instances the signal strength was low" (Figure 9). Earlier readings still
classified, but less reliably, and the authors note that more enzyme would
raise the signal strength.

Chemart stores these counts and its default network is this device's recipe.
The corresponding test checks the stored numbers and that a hand-set response
(high in state b, low in a and c) classifies the 135 patterns correctly; it
does not simulate the device, which would require the measured responses.

### What the authors did and did not claim

The reaction conditions were far from those in a cell: the magnesium levels
are probably higher than any plausible change inside a cell, and the pH was
about 10. The authors therefore draw no conclusion about enzymes in living
cells, and they do not expect enzymes to compete with electronics as logic
gates. Their point is that proteins can supply *nonlinear* input-output
transforms that are precisely reproducible, since every copy of a protein is
identical, unlike the statistical aggregates of solid-state devices. Networks
of such components could then compute a given function with far fewer parts
than a network of logic gates. They call the two-bit experiment "a stepping
stone" and say that response surfaces for more complex signal patterns are
needed to judge the approach.

### Conrad's ideas in the book

The book credits Conrad's enzymatic computing with contributing to the
foundations of wet molecular computing (§11.2). It reports two further ideas of
his, from the 1985 and 1992 papers, which Chemart does not model:

- **Quantum parallelism.** Conrad held that free-energy minimisation, as in
  protein folding, lock-key docking or conformational change, uses quantum
  superposition to search for the minimum-energy state in parallel.
- **Programmability versus evolvability.** Because the relation between a
  protein's sequence and its shape is hard to predict, such computers would be
  programmed by learning or evolution rather than by writing instructions.
  Conrad stated this as a general principle (quoted in §17.1.2): "a system
  cannot at the same time be effectively programmable, amenable to evolution by
  variation and selection, and computationally efficient." The book argues it
  still holds, pointing to programs made by genetic programming, which look
  nothing like those written by people, and warning that making chemical
  computers programmable like electronics may cost them their evolvability.

The book also cites a later chapter by Conrad and Zauner (2003) that calls the
approach *conformation-based computing*. Chemart has not read Conrad's own
papers, which are paywalled, so nothing on this page is attributed to them
directly.

## Further reading

- Zauner, K.-P. & Conrad, M. (1996). Simulating the interplay of structure,
  kinetics, and dynamics in complex biochemical networks. *Proc. German
  Conference on Bioinformatics (GCB'96)*, 336–338. Book reference [952].
- Zauner, K.-P. & Conrad, M. (1998). Conformation-driven computing: simulating
  the context-conformation-action loop. *Supramolecular Science* 5(5–6),
  791–794. Book reference [953]; the book (§2.3.3) cites these two as examples
  of hybrid simulations, in which single macromolecules are simulated
  explicitly and small molecules as concentrations.
- Arkin, A. & Ross, J. (1994). Computational functions in biochemical reaction
  networks. *Biophysical Journal* 67, 560–578. The simulation study that,
  according to Zauner and Conrad, showed an enzymatic XOR to be possible in
  theory.
- Minsky, M. L. & Papert, S. (1969). *Perceptrons: An Introduction to
  Computational Geometry*. MIT Press. The result that a single-layer
  perceptron cannot compute the XOR.
