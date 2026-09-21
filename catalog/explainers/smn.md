## Introduction

The String Metabolic Network (SMN) is an artificial chemistry built to ask
why real metabolisms are organised the way they are. It was proposed by
Naoaki Ono, Yoshi Fujiwara and Kikuo Yuta at the European Conference on
Artificial Life in 2005, in a paper titled "Artificial metabolic system: an
evolutionary model for community organization in metabolic networks".

The question comes from network biology. A metabolism can be drawn as a
graph: chemicals and reactions are nodes, and links join each reaction to
the chemicals it uses and makes. Studies of real organisms around 2000 found
that these graphs tend to be scale-free or small-world (Jeong et al. 2000;
Wagner and Fell 2001) and often *modular* (Ravasz et al. 2002; Parter et al.
2007). Modular means the graph falls into clusters, or *communities*: groups
of chemicals and reactions tightly linked among themselves, with few links
between groups. Where does that structure come from? One way to find out is
to evolve metabolisms from scratch in a world simple enough to control, and
see whether the same structure appears.

SMN supplies that world. Its chemicals are strings of letters such as `adbg`
or `eefg`, and it has only two kinds of reaction, both reversible: two
strings can be joined end to end (and a string can be cut in two), or two
strings can swap their tails. Every reaction is catalysed by its own enzyme,
and an organism is defined by the enzymes it owns: its genome is simply a
list of the reactions it can carry out. Evolution changes the genome by
copying an enzyme and altering the copy, and selection favours organisms
that build long strings. Ono and colleagues evolved organisms this way and
compared the community structure of their metabolic networks with that of
real organisms.

SMN is a simulation model. Banzhaf and Yamamoto describe it in their chapter
on modelling biological systems, under biochemical pathways (book §18.3.1),
as one of two artificial chemistries used to study how metabolic networks
evolve. The other is the bond-number chemistry of Hintze and Adami, whose
reactions also swap the tails of two molecules; its catalog entry,
`bnc-cell`, is archived. Other string chemistries in the catalog ask
different questions: in
[Kauffman's autocatalytic sets](kauffman-autocatalytic-sets.md) and the
[Bagley-Farmer metabolism](bagley-farmer.md), polymers join and split and
any polymer may catalyse a reaction, and the question is whether a
self-sustaining set arises at all. In SMN the catalysts are not molecules
but the genome, and the question is about the *shape* of the network that
evolution builds.

## How it works

### Molecules, reactions, enzymes

A molecule is a string over the eight letters `a` to `h`. The letters stand
for nothing in particular; they are the atoms of the chemistry, and a string
is a compound built from them. There are two reaction classes (book
equations 18.17 and 18.18):

```
ligation / cleavage:   A + B   <->  AB
recombination:         AB + CD <->  AD + CB
```

Here `A`, `B`, `C` and `D` are pieces of strings. Read forwards, a
*ligation* joins two strings end to end; read backwards, the same reaction
is a *cleavage* that cuts one string in two. A *recombination* cuts each of
two strings at one point, the *recombination point*, and swaps the tails.
The book's examples are `adbg + ef <-> adbgef` and
`facbaha + eefg <-> facbfg + eeaha`: in the second, `facbaha` is cut after
`facb` and `eefg` after `ee`, and the tails `aha` and `fg` change places.
The strings that go into a reaction are its *educts*.

Each reaction is catalysed by a dedicated enzyme, and the book treats the
enzyme and its reaction as the same thing. An organism's genome is its list
of enzymes. Its *compound set* is every string that takes part in those
reactions: the initial metabolites it starts life with, plus everything its
enzymes can make from them. The metabolism of one organism is the network of
its enzymes' reactions over its compound set, and that network is what one
Chemart run produces.

### A worked example

The default run (`seed=1`, shown in *Using it in Chemart* below) starts from
six random initial metabolites, `eg`, `abg`, `bcg`, `cg`, `d` and `ea`, and a
genome of twelve random enzymes. Its fifth enzyme, `E4`, is a
recombination:

```
b|cg + e|g   gives   bcg + eg <-> bg + ecg
```

The bars mark the recombination points. `bcg` is cut into `b` and `cg`, `eg`
into `e` and `g`, and the tails swap: `b` + `g` makes `bg`, and `e` + `cg`
makes `ecg`. Both products are new compounds. The next enzyme, `E5`, is the
ligation `abg + ecg <-> abgecg`, which uses `ecg`, a product of `E4`. This
is how the compound set grows: each enzyme acts on compounds that are already
there, and its products become raw material for the enzymes after it. The
twelve enzymes of the default run turn six initial strings into 21
compounds.

Count the letters on both sides of `E4`: one `b`, one `c`, one `e` and two
`g`s on each side. Every ligation, cleavage and recombination only
rearranges letters, never creates or destroys them, so the number of each
letter is conserved. The total number of letters, which the book calls the
organism's *mass*, is conserved by the reactions too.

### Flow, mass and fitness

Reactions do change the *number* of molecules. A ligation turns two
molecules into one; a cleavage turns one into two. The organism also has a
regulated inflow and outflow of substances, which keeps the total amount of
substance in the cell constant. Put the two together: if the organism turns
short strings into long ones, it holds the same number of molecules but more
letters. Its mass grows. The book gives this as the reason mass is a good
fitness measure: a larger mass means the organism is able to make longer
compounds, and the kinetic coefficients are set so that longer compounds are
harder to make.

### Mutation and evolution

The book's mutation operator works on the genome. Pick one enzyme at random,
copy it, replace one of the copy's educts with another compound already in
the compound set, choose a new recombination point where there is one, and
add the copy to the genome. The book compares this to a gene duplication
followed by a mutation. The genome only grows, and each mutant can create
new compounds for later enzymes to use.

Evolution in the paper repeats three steps: mutate the genomes of the
organisms, run each organism's metabolism for a fixed number of time steps,
and use each organism's mass as its fitness to select individuals for
the next generation.

Chemart builds one organism's network, with or without mutations, but not the
evolutionary loop. The book gives no kinetic coefficients, only that longer
products are harder to make, so the network carries no rates, and without
rates a metabolism cannot be run to measure its mass. The formal
specification below summarises the model; the *Implementation decisions*
list what Chemart had to choose where the book is silent.

## Using it

The default run builds a random organism: six distinct initial metabolites
of one to three letters, each at amount 1.0, and a genome of twelve enzymes
drawn one at a time, each acting on compounds already present, half
ligations and half recombinations on average. The book does not say how
initial genomes are made, so this is a Chemart choice. No mutations
are applied by default. Every enzyme contributes two directed reactions,
forward and backward, which is why twelve enzymes give 24 reactions. All
rates are `None`.

The genome is in `net.extras["genome"]`, one record per enzyme:

```python
net.extras["genome"][4]
# {'enzyme': 'E4', 'kind': 'recombination', 'educts': 'b|cg + e|g',
#  'reaction': 'bcg + eg <-> bg + ecg', 'reactions': [8, 9], 'parent': None}
```

`educts` is the enzyme in the notation Chemart uses for input, `reactions`
holds the indices of its two directed reactions in `net.reactions`, and
`parent` names the enzyme it was copied from, for mutants. The other extras
are `initial_metabolites`, `conservation` (one conservation law per letter),
and three text fields recording the flow law, the missing kinetics and the
book's fitness definition. `net.outflow` is `"constant-total"`, the
regulated flow.

**The book's examples.** Pass your own genome as a list of enzymes, `"A + B"`
for a ligation and `"A|B + C|D"` for a recombination. With no
`initial_metabolites`, the educts of the genome become the initial
metabolites:

```python
net = chemart.generate_network("smn", genome=["adbg + ef", "facb|aha + ee|fg"])
print(net.to_text())
```

```
adbg + ef -> adbgef
adbgef -> adbg + ef
facbaha + eefg -> facbfg + eeaha
facbfg + eeaha -> facbaha + eefg
```

**Growing a genome by mutation.** `mutations` applies the book's
duplicate-and-mutate operator that many times, without selection. Twenty
mutations take the default organism from 21 to 43 compounds:

```python
net = chemart.generate_network("smn", seed=1, mutations=20)
net.summary().splitlines()[0]   # 'smn: 43 species, 62 reactions, status=complete'
for e in net.extras["genome"][12:14]:
    print(e["enzyme"], e["parent"], e["educts"], "=>", e["reaction"])
```

```
E12 E8 e|aeg + ecgcgb|cg => eaeg + ecgcgbcg <-> ecg + ecgcgbaeg
E13 E0 abgecg + bcg => abgecg + bcg <-> abgecgbcg
```

`E12` is a copy of `E8` (`e|aeg + c|gcg`) whose second educt was replaced by
`ecgcgbcg` and given a new recombination point. The 32 enzymes give 62
rather than 64 reactions because two enzymes catalyse the same reaction and
share it. This is the book's control experiment, random mutation without
selection, as a network you can analyse.

**Checking conservation and computing mass.** Each entry of
`extras["conservation"]` is a vector that the stoichiometry leaves
unchanged, and `letter_mass` computes the book's fitness measure for any
state. Continuing with the mutated network:

```python
import numpy as np
from chemart.chemistries.smn import letter_mass
ids, R, P = net.matrices()
S = (P - R).toarray()
law = net.extras["conservation"][0]                     # letter a
v = np.array([law["vector"].get(s, 0) for s in ids])
np.abs(v @ S).max()                                      # 0
letter_mass(net.initial_state)                           # 13.0: six initial strings, 13 letters
```

To make the metabolism run you must supply rates yourself; the book gives
none. `alphabet_size` narrows the alphabet, `max_length` caps compound
length (the book has no cap), and the other population settings are in the
parameter table below. Every run here takes well under a second.

## Results

**Evolved metabolisms look like natural ones.** The main finding, as the
book reports it, concerns community structure. A community is a module:
the network is split into groups by an algorithm that minimises the number
of links between groups, and counting the nodes in each group gives the
distribution of community sizes. Figure 18.6 of the book, reprinted from the
paper, compares two such distributions. The left plot is the average over
100 organisms with large metabolic networks, taken from a database of
genomes and metabolic pathways; the right plot is the average over 30 runs
of the SMN evolutionary algorithm. Both follow a power law, and the two
curves are very similar. Chemart does not reproduce this result: it needs
the paper's kinetic coefficients and the evolutionary loop with selection,
neither of which the book describes in enough detail. The paper itself is
closed access, and no open copy was found.

**Selection is necessary.** The authors also report, according to the book,
that networks grown by random mutation without selection pressure do not
become modular. Chemart provides the control network itself, through the
`mutations` parameter, and its tests check that each mutant is a copy of an
earlier enzyme with exactly one educt replaced by a compound already in the
compound set. It does not measure modularity, so the finding is not
reproduced either.

**Letters are conserved.** Ligation, cleavage and recombination only
rearrange letters. Chemart's tests check that every per-letter vector in
`extras["conservation"]`, and the total letter count, is left unchanged by
every reaction of a network with 40 random enzymes and 60 mutations.

**Constant amount, growing mass.** Because the regulated flow holds the
amount of substance constant, the organism's letter mass can only grow by
making longer compounds, which is why the book uses mass as fitness. Chemart
represents the flow as a constant-total outflow. A test gives a small
three-enzyme organism temporary mass-action rates (the book has none),
integrates it, and checks that the total amount stays at its initial value
while the letter mass rises as ligation builds longer strings.

The tests also check the book's two example reactions, that every enzyme is
a reversible ligation or recombination of the right form, that
`alphabet_size` and `max_length` bound the compounds, and that malformed
genomes are rejected.

**Related work.** The book contrasts SMN with the model of Hintze and Adami
(2008). SMN organisms live in a fixed environment, set by their initial
metabolites and the regulated flow; Hintze and Adami evolved metabolisms in
environments that change during an organism's life, and found that networks
evolved in dynamic environments evolved more slowly and were less modular but more robust
to environmental noise. That model is catalogued as `bnc-cell`, which is
archived.

## Further reading

- Jeong, H., Tombor, B., Albert, R., Oltvai, Z. N. & Barabási, A.-L. (2000).
  The large-scale organization of metabolic networks. *Nature* 407, 651–654.
- Wagner, A. & Fell, D. A. (2001). The small world inside large metabolic
  networks. *Proceedings of the Royal Society B* 268(1478), 1803–1810.
- Ravasz, E., Somera, A. L., Mongru, D. A., Oltvai, Z. N. & Barabási, A.-L.
  (2002). Hierarchical organization of modularity in metabolic networks.
  *Science* 297(5586), 1551–1555.
- Parter, M., Kashtan, N. & Alon, U. (2007). Environmental variability and
  modularity of bacterial metabolic networks. *BMC Evolutionary Biology*
  7, 169.
- Hintze, A. & Adami, C. (2008). Evolution of complex modular biological
  networks. *PLoS Computational Biology* 4(2), e23.
