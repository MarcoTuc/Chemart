## Introduction

In 1994 Leonard Adleman, a computer scientist at the University of Southern
California, solved a small graph problem in a test tube. He encoded a map of
seven "cities" and fourteen one-way "roads" in short strands of DNA, let the
strands stick together and join into every possible route at once, and then
used standard laboratory methods to throw away every route that was not the
answer. What was left was a single kind of DNA molecule whose sequence spelled
out the solution. The book describes his paper in *Science* as the first
successful wet-lab computation with DNA strands, and as the experiment with
which wet DNA computing "took off" (§19.3.1).

The problem was the **directed Hamiltonian path problem** (HPP): given a
graph whose edges have a direction, a start vertex and an end vertex, is there
a path from start to end that enters every other vertex exactly once? The
problem is NP-complete, and, as the paper notes, every known algorithm for it
needs exponential time in the worst case. Adleman's paper opens by recalling Feynman's 1959 talk on "sub-microscopic"
computers and says its aim is to explore "the possibility of computing
directly with molecules". The idea is brute force made cheap by chemistry: a
tube holds around 10^13 copies of each strand, so all candidate routes are
built in parallel, and the hard part becomes picking the right molecule out of
the mixture.

This is a **wet-lab experiment**, not a simulation model. Chemart's entry is
an in-silico account of it: the strands as species with their real (or
randomly drawn) sequences, the joining step as reactions, and the four
selection steps as an analysis of which molecules each step keeps. It
reproduces the logic of the experiment and the published sequences and band
lengths; it does not simulate the chemistry's kinetics or its errors.

The book presents the experiment in chapter 19 (wet artificial chemistries),
section 19.3.1 on DNA computing, next to the DNA automaton of Shapiro and
Benenson ([dna-automaton](dna-automaton.md)). The two are a contrast: in
Adleman's computation a person carries out each step at the bench, while the
automaton runs by itself once its parts are mixed. Adleman's problem is a
special case of the travelling salesman problem, which the book also treats
with a purely in-silico chemistry, the [molecular TSP](molecular-tsp.md); there
the molecules are candidate tours improved by collisions, not DNA strands
filtered by a chemist.

## How it works

### The code: vertices, edges and splints

DNA is a chain of four bases, A, C, G and T, with a direction: one end is
called 5' and the other 3', and sequences are written from 5' to 3'. Two
strands stick together (**hybridise**) when they are **complementary**: A
faces T, C faces G, and the strands run in opposite directions. An enzyme,
**DNA ligase**, can then join two strands that lie end to end on a third one
into a single strand.

Adleman used these three facts to write the graph in DNA:

- Each vertex `i` gets a random 20-base code word `O_i`. Chemart calls it `Oi`,
  for example `O1`.
- Each edge `i -> j` gets an **edge oligonucleotide** (a short synthetic
  strand, "oligo" for short): the last 10 bases of `O_i` followed by the first
  10 bases of `O_j`. Chemart calls it `Ei-j`. Because the halves come from
  different ends of the code words, `E2-3` is not the same strand as `E3-2`,
  so the direction of the road is kept.
- Each vertex also gets the complement of its code word, `Ō_i`, written
  `Obari` in Chemart. It acts as a **splint**: its two halves pair with the
  end of an edge arriving at vertex `i` and the start of an edge leaving it,
  holding them end to end so ligase can join them.
- The two ends are special: edges leaving the start vertex carry the whole
  code word of the start, and edges entering the end vertex carry the whole
  code word of the end. A finished path of `n` vertices is then exactly
  `20 × n` bases long, which is 140 for Adleman's seven-vertex graph.

The paper chose random 20-mers because random code words are unlikely to share
long stretches (which could bind by mistake) or fold back on themselves, and
because a 10-base overlap between splint and edge is stable at room
temperature.

### The five steps

The paper states the method as a five-step algorithm and carries each step out
with a standard technique:

1. **Generate random paths.** 50 picomoles (pmol) of each edge oligo and of
   each splint (except those of the start and end vertices) go into one tube
   with ligase for 4 hours at room temperature. Strands join into paths of
   every length through the graph.
2. **Keep the paths from start to end.** The polymerase chain reaction (PCR)
   copies DNA exponentially, but only the stretch between two short
   "primer" strands that bind at its ends. Using `O_0` and `Ō_6` as primers
   multiplies only the molecules that begin at vertex 0 and end at vertex 6.
3. **Keep the paths with exactly seven vertices.** Gel electrophoresis sorts
   DNA by length; Adleman cut out the 140 base-pair (bp) band.
4. **Keep the paths that enter every vertex.** The strands are separated and
   passed over magnetic beads carrying `Ō_1`; only strands containing `O_1`
   stick. Repeating this with `Ō_2` to `Ō_5` keeps the paths that enter all
   five inner vertices. Together with step 3, a path of seven slots that enters
   all seven vertices enters each exactly once.
5. **Read the answer.** If anything is left, the answer is yes. Adleman also
   read out *which* path it was, by **graduated PCR**: tube `i` uses primers
   `O_0` and `Ō_i`, so it copies the stretch from the start to wherever vertex
   `i` sits. A vertex in position `k` of the path (counting the start as 0)
   gives a band of `20 × (k + 1)` bp. The Hamiltonian path `0->1->2->3->4->5->6`
   gives 40, 60, 80, 100, 120 and 140 bp in lanes 1 to 6.

Only step 1 makes new molecules. Steps 2 to 5 select among them.

### What Chemart turns into a network

Chemart's reactions are the joining events of step 1. Each one takes a path
that ends at vertex `v`, the splint of `v`, and an edge leaving `v`, and makes
a path one vertex longer. The first reaction of the default network (the list
under *Using it in Chemart*, below) is

```
E0-1 + Obar1 + E1-2 -> P0-1-2
```

With seed 1 the strands are

```
E0-1    CGTTAATTACTCCTCCGGAA TTTGTCCTAC          (all of O0, then the first half of O1)
E1-2                        ACTACCTAGC TATCGGATCG (second half of O1, first half of O2)
O1                TTTGTCCTACACTACCTAGC
Obar1             GCTAGGTAGTGTAGGACAAA            (the complement of O1, written 5'->3')
P0-1-2  CGTTAATTACTCCTCCGGAATTTGTCCTACACTACCTAGCTATCGGATCG   (50 bases)
```

Where `E0-1` ends and `E1-2` begins, the two halves of `O1` meet, and
`Obar1` pairs with exactly those 20 bases. Ligase seals the gap and the
product `P0-1-2` is the path 0 -> 1 -> 2. Species named `P` followed by
vertices are such paths; a two-vertex path is just its edge oligo. A path
molecule still carries its splints, so the product weighs 50 + 20 = 70
nucleotides, the same as the 30 + 20 + 20 of the three reactants: every
reaction conserves nucleotides, and `net.extras["conservation"]` holds the
weights.

Because every reaction needs three molecules (path, splint, edge), the
molecules-per-reaction count in the specification below is 3. The reaction
set is built by closure: start from the edges and keep extending every path by
every edge leaving its last vertex. Adleman's graph has cycles (for example
1 -> 2 -> 1), so this never ends on its own. Chemart stops at a maximum number
of vertex slots per path, which defaults to 7, the length step 3 keeps, and
marks the network `truncated`. No rate constants were published, so reactions
carry no rates, and the initial state is the ligation mix in pmol.

Steps 2 to 5 are not reactions. Chemart applies them to the set of paths it
built and records what each one keeps in `net.extras["analysis"]`.

## Using it

The default call is Adleman's experiment: his graph (figure 1 of the paper),
start vertex 0, end vertex 6, 20-base code words, and the three code words
`O2`, `O3` and `O4` that figure 2 prints. The other four code words are drawn
at random from the seed, as Adleman drew his. The network has 797 species
(14 vertex code words and splints, 14 edge oligos and 769 longer paths) and
769 reactions.

The five steps are in `net.extras["analysis"]["steps"]`:

```python
import chemart
net = chemart.generate_network("dna-hpp", seed=1)
a = net.extras["analysis"]
for s in a["steps"]:
    print(s["step"], s["operation"], s.get("kept", s.get("molecules")))
print(a["gel_survivors"])
print(a["hamiltonian_paths"], a["solution_species"])
```

```
1 ligation (T4 DNA ligase) 783
2 PCR 5
3 agarose gel 2
4 affinity purification on magnetic beads 1
5 PCR and graduated PCR readout 1
[[0, 1, 2, 3, 4, 5, 6], [0, 3, 2, 3, 4, 5, 6]]
[[0, 1, 2, 3, 4, 5, 6]] ['P0-1-2-3-4-5-6']
```

Of the 783 path molecules up to seven slots (131 of which start at vertex 0),
PCR keeps the five that run from 0 to 6: `E0-6` (40 bp), `P0-3-4-5-6`
(100 bp), `P0-1-3-4-5-6` (120 bp), and two of 140 bp. The gel keeps those
two: the Hamiltonian path and `0->3->2->3->4->5->6`, which enters vertex 3
twice and misses vertex 1. The first bead round, on `Obar1`, removes the
second, and one molecule is left. `kept_after` in step 4 lists the count after
each bead round.

The graduated-PCR readout of every gel survivor is in
`a["graduated_pcr"]`; an empty list is an empty lane:

```
P0-1-2-3-4-5-6 {'1': [40], '2': [60], '3': [80], '4': [100], '5': [120], '6': [140]}
P0-3-2-3-4-5-6 {'1': [], '2': [60], '3': [40, 80], '4': [100], '5': [120], '6': [140]}
```

The sequences are in `species.structure`, and `net.extras["encoding"]` keeps
the figure 2 strands for comparison. `a["experiment"]` records the protocol
as text, including the number of copies that 50 pmol stands for (about
3.0 × 10^13).

**Instances with no answer.** The paper names two variants of its graph that
have no Hamiltonian path. Both run in Chemart:

```python
from chemart.chemistries.dna_hpp import ADLEMAN_GRAPH
no23 = chemart.generate_network("dna-hpp", seed=1,
                                graph=[e for e in ADLEMAN_GRAPH if e != [2, 3]])
moved = chemart.generate_network("dna-hpp", seed=1, v_in=3, v_out=5)
```

PCR still keeps 3 and 9 paths respectively, but none survives the gel, so
`hamiltonian_paths` is empty.

**Your own graph.** Pass `graph` as a list of `[u, v]` edges, with `vertices`,
`v_in` and `v_out`. Set `sequences={}` to draw every code word at random (the
default's three published code words only make sense on Adleman's graph), and
set `max_path_vertices` to the number of vertices. A three-vertex example,
which has no cycle and so gives a complete network:

```python
net = chemart.generate_network("dna-hpp", seed=1, graph=[[0, 1], [1, 2], [0, 2]],
                               vertices=3, v_in=0, v_out=2,
                               max_path_vertices=3, sequences={})
print(net.summary().splitlines()[0])
print([r.to_text() for r in net.reactions])
```

```
dna-hpp: 10 species, 1 reactions, status=complete
['E0-1 + Obar1 + E1-2 -> P0-1-2']
```

**Longer walks.** Raising `max_path_vertices` builds the longer walks that the
real tube also makes and the gel throws away. The answer does not change, but
the network roughly doubles with each extra slot on Adleman's graph: 797,
1,590, 3,177 and 6,352 species for 7 to 10 slots. All of these run in well
under a second. `vertices` is capped at 14 and `max_path_vertices` at 12.

The seed changes only the random code words, never the answer: seed 2 gives a
different `O0` and the same Hamiltonian path.

## Results

**One molecule, the right path.** Adleman ran the five steps on the graph of
figure 1. The ligation product in step 1 appeared on the gel as a smear with
striations, as expected from paths of many lengths, and the PCR of step 2 gave
dominant bands for paths from vertex 0 to vertex 6. Graduated PCR on the final
product of step 4 (figure 3C) showed the bands of the path
`0->1->2->3->4->5->6`, which the figure 1 legend states is the graph's unique
Hamiltonian path. Chemart reproduces this: its tests check that the five
steps leave exactly the molecule `P0-1-2-3-4-5-6`, that brute-force search
agrees there is only this one path, and that each step keeps a subset of the
one before.

**The 140 bp band.** With 20-base code words and whole code words at both
ends, a seven-vertex path is 140 bp, and that is the band Adleman cut out.
The tests check that every gel survivor is a 140 bp path from 0 to 6, and
that the Hamiltonian path is the ligation of its six edge oligos.

**Graduated PCR as a printer.** The paper works out the lanes for three paths
of 140 or 120 bp: 40, 60, 80, 100, 120, 140 bp for the Hamiltonian path; 40,
nothing, 60, 80, 100, 120 for `0->1->3->4->5->6`, which skips vertex 2; and
nothing, 60, "80bp/40bp", 100, 120, 140 for `0->3->2->3->4->5->6`, which
passes vertex 3 twice. The tests reproduce all three patterns exactly.

**Imperfect selection in the tube.** Graduated PCR of the gel band (figure 3B)
showed the superimposed bands of all three paths above. The paper notes that
the 120 bp path was not expected and suggests the excised band was
contaminated with 120 bp molecules, adding that such contamination "does not
persist through Step 4". Chemart's gel is exact, so the 120 bp path never
passes it; the tests check that it is 120 bp long and absent from the gel
survivors, and that the double-passage path survives the gel and is removed by
the first bead round. The leak itself is not simulated.

**Graphs with no path.** The paper states that removing the edge 2 -> 3, or
moving the start and end to vertices 3 and 5, leaves no Hamiltonian path
(since, for the second, no edge enters vertex 0). The tests run both variants
and check that the selection ends empty while PCR still keeps paths. A further
test runs 25 random graphs of 4 to 6 vertices and checks that the five steps
find exactly the Hamiltonian paths that brute force finds.

**The published code.** Figure 2 prints `O2`, `O3`, `O4`, the edge oligos
`O2->3` and `O3->4`, and the splint `Ō3`. The tests check that the generator
rebuilds the two edge oligos and the splint from the three code words by the
halving rule, and that every splint pairs base by base with its code word.

**Cost, speed and energy.** The paper counted the resources. The number of
different oligos grows linearly with the edges and the number of lab
procedures linearly with the vertices, but the quantity of each oligo needed
to form a Hamiltonian path with high probability "should grow exponentially
with the number of vertices". The seven-vertex run took about 7 days, with
step 4 alone taking a full day at the bench. Adleman estimated that step 1
performed about 10^14 ligations, and that 10^20 or more seemed plausible with
larger quantities, faster than the supercomputers of the time; that one
ligation costs the energy of one ATP molecule, so about 2 × 10^19 operations
per joule against a thermodynamic limit of 34 × 10^19 and at most 10^9 for
supercomputers; and that DNA stores about 1 bit per cubic nanometre. He also
noted that for this graph "sub-attomol quantities" would probably have
sufficed instead of the 50 pmol used. Chemart records the copy number that
50 pmol stands for (tested at about 3 × 10^13, as the paper says), but none of
these estimates is computed.

**What followed.** The book summarises what came next. Lipton (1995) extended
the procedure to the satisfiability problem (SAT), and other DNA algorithms
followed for graph colouring, maximal clique and the knight problem. The
exponential growth in DNA proved decisive: Hartmanis (1995) showed that a
200-city instance solved by Adleman's method would need more than the weight
of the Earth in DNA, and the field concluded that NP-hard problems were not
DNA computing's "killer application". It moved instead toward autonomous
molecular computation, self-assembly and diagnosis, the direction of the
[DNA automaton](dna-automaton.md). The book also lists the experiment's
limits: its procedures were slow, labour-intensive and error-prone. Chemart
states the scaling claim as text only (in `a["experiment"]["scaling"]`); it
is not derived from a run, and the tests do not check it.

**A note on the book's encoding.** Section 19.3.1 describes the edges as
complemented halves of the vertex code words. That is not the paper's scheme,
where the edges are uncomplemented halves and the complements are the
splints. Chemart follows the paper; the implementation decisions give the
details.

## Further reading

- Lipton, R. J. (1995). DNA solution of hard computational problems.
  *Science* 268(5210), 542–545. The extension to SAT.
- Amos, M. (2005). *Theoretical and Experimental DNA Computation*. Springer.
  A book-length account, recommended by Banzhaf and Yamamoto as a tutorial.
- Amos, M. (2009). DNA Computing. In R. A. Meyers (ed.), *Encyclopedia of
  Complexity and Systems Science*, vol. 4, 2089–2104. Springer. A shorter
  review, also recommended by the book.
