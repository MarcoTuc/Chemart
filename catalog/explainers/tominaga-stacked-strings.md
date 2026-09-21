## Introduction

Kazuto Tominaga and his colleagues at the Tokyo University of Technology built
this chemistry as a *modelling language* for real biochemistry. Most artificial
chemistries are virtual worlds, run to see what lifelike behaviour emerges.
Tominaga's aim was different: to write down known molecular processes, such as
DNA computers or the pathway from DNA to protein, precisely enough that a
simulator can execute them, and to check on paper or on a computer that a
proposed mechanism does what it should.

A molecule is a string of symbols, or a **stack** of strings, each line shifted
left or right against the first. That one step beyond plain strings is what
makes biology writable. A DNA double strand is two stacked lines of bases; a
sticky end is a lower line that sticks out past the upper one; an enzyme sitting
on DNA is a third line placed above the stretch it binds. Reactions are
*recombination rules*: patterns with wildcards that pick out molecules of a
given shape, remove them, and put together new molecules from the matched
pieces. The simplest example, from the 2007 paper, joins copies of the string
`AB` into `ABAB`, `ABABAB` and so on, with a helper molecule `CD` that lies
across the joint like a splint and is released afterwards:

```
ABAB          AB           ABABAB
   CD    +          ->        CD
```

The chemistry was first described in a 2004 technical report, and two journal
papers show what it can express. Tominaga et al. (2007) prove that it can
emulate any Turing machine using only rules with one or two reactants, which
matters to them because real elementary reactions are almost all unimolecular
or bimolecular, and then model two landmark DNA computers: Adleman's
Hamiltonian-path experiment and Benenson's DNA automaton. Tominaga et al. (2009)
model biochemical pathways: DNA replication, transcription, translation and the
oxidation of fatty acids, and add a procedure that reasons backwards over a
model.

It is a formalism with a simulator, not a system studied for emergent
dynamics: there are no rates, energies or space, only which reactions can happen.
Banzhaf and Yamamoto mention it in one paragraph of their chapter on modelling
biological systems (book §18.3.2), in the section after the
[String Metabolic Network](smn.md).
Its closest relatives in the catalog are the laboratory systems it models,
[the DNA automaton](dna-automaton.md) and [Adleman's DNA computation](dna-hpp.md),
and the [music-composition chemistry](music-ac.md), which reuses the same stacked
notation for notes and chords. Unlike string chemistries such as
[SAC](sac.md) or [AlChemy](alchemy.md), where the symbols themselves carry
functions, here a symbol means nothing on its own: all behaviour is in the rules.

## How it works

### Molecules: stacked lines with displacements

A molecule (the papers say *object* in 2007 and *v-molecule* in 2009) is written
line by line as `displacement#elements/`. The displacement says how far a line
is shifted relative to the first line, in element widths. So

```
0#ABAB/3#CD/      is      ABAB
                             CD
```

because `CD` starts three positions to the right of the start of `ABAB`. An
element is a capital letter optionally followed by lower-case letters, so
`Rp` (RNA polymerase), `Cap`, `Coa` or `Fad` are single elements, just like the
bases `A`, `C`, `G`, `T`. How much an element stands for is up to the modeller:
a base, a ten-base block of DNA, a whole enzyme. A line may even be empty, as a
placeholder that fixes where a partner strand will go.

### Patterns and rules

A pattern looks like a molecule but may contain wildcards:

- a digit such as `1` matches exactly one element;
- `*1` at the start of a line, or `1*` at the end, matches any run of zero or
  more elements.

The pattern `0#*1AB/1#CD/` matches any two-line molecule whose upper line ends
in `AB` and whose lower line is exactly `CD`, with `C` under the `B`. It matches
`0#ABAB/3#CD/` (with `*1` = `AB`) but not `0#ABAB/1#CD/`, where `CD` sits under
the wrong letters. When counting displacements in a pattern, a sequence
wildcard counts as length zero, which is why the pattern says `1#` where the
molecule says `3#`.

A **recombination rule** `lhs -> rhs` has patterns on both sides. When
molecules in the pool match the left-hand patterns, they are removed, and the
right-hand patterns, filled in with whatever the wildcards matched, are added.
Every wildcard on the left must appear exactly once on the right. Here is rule
(2) of the `AB` example applied in Chemart:

```
0#*1AB/1#CD/ + 0#AB2*/ -> 0#*1ABAB2*/1#CD/
0#ABAB/3#CD/ + 0#AB/   ->  0#ABABAB/3#CD/
```

A system also has **sources**, which supply a given molecule without limit, and
**drains**, which remove any molecule matching a pattern. The published dynamics
is deliberately unspecified: at each step, apply some rule, operate some source
or operate some drain, in any order. There is no spatial structure and no rate
constant. The 2007 paper argues this suits DNA computing, which happens in a
well-stirred test tube, and that rates belong to whichever simulator one builds
on top.

### Worked example: one step of the DNA automaton

The default system is Benenson's DNA automaton as written in the 2007 paper. The
automaton reads a word over the letters a and b and has two states, S0 and S1;
it starts in S0, a `b` flips the state, and it accepts (ends in S0) when the word
has an even number of b's. In the laboratory the word is a DNA double strand,
and the enzyme Fok I does the reading: it binds at the sequence `GGATG` and cuts
the upper strand 9 bases and the lower strand 13 bases further on, leaving a
four-base single-stranded sticky end. Fok I itself is written as `FFFFFFFFF`.

This is the input for the word `abb` (`X` stands for arbitrary bases that Fok I
ignores; `CTGGCT` is the symbol a, `CGCAGC` is b, `TGTCGC` the terminator):

```
GGATGXXXXXXXCTGGCTCGCAGCCGCAGCTGTCGCX
CCTACXXXXXXXGACCGAGCGTCGGCGTCGACAGCGX
```

Rule (3) binds Fok I five positions to the right of the start of `GGATG`, and
rule (4) cuts. These are the two products of the run, drawn with Chemart's `layout`
and placed side by side:

```
     FFFFFFFFF                               GGCTCGCAGCCGCAGCTGTCGCX
GGATGXXXXXXXCT                  +                GCGTCGGCGTCGACAGCGX
CCTACXXXXXXXGACCGA
```

Rule (5) releases Fok I from the spent head on the left. The piece on the right
now starts with an overhanging `GGCT`. Its sequence encodes "state S0, symbol a",
and only one of the four transition molecules has the complementary overhang
`CCGA` on its lower strand, the one for S0 reading a and staying in S0:

```
GGATGTAC
CCTACATGCCGA
```

Rule (6) sticks the two together and joins the strands:

```
GGATGTACGGCTCGCAGCCGCAGCTGTCGCX
CCTACATGCCGAGCGTCGGCGTCGACAGCGX
```

The transition molecule brought a new `GGATG`, so Fok I can cut again. Because
the spacer between `GGATG` and the input differs between transition molecules,
the next cut falls at a position that encodes the new state together with the
next symbol. For `abb` the sticky ends run `GGCT` (S0, a), `CAGC` (S0, b), `CGCA`
(S1, b) and finally `TCGC`, the terminator read in state S0. That end fits the
S0 detector `0#X/0#XAGCG/`, and rule (10) forms the *reporter*
`0#XTCGCX/0#XAGCGX/`, the molecule whose presence signals acceptance. The trailing
`X` tag on each input survives into its reporter, so a pool with several inputs
shows which of them were accepted.

### What Chemart computes

Chemart turns a system into a reaction network in one of two ways. The default,
`method="closure"`, applies every rule to every combination of known species
until nothing new appears, and returns all reachable reactions. Sources become
reactions `∅ -> s` and drains `s -> ∅`. With `method="soup"` Chemart instead
samples one run of the nondeterministic process from the published copy numbers
and returns the reactions that fired, with how often. That sampler is a Chemart
addition, since the papers do not describe their simulators.

## Using it

The default call above is the 2007 automaton with its published pool: 100 Fok I,
20 of each transition molecule, 20 detectors, and 10 each of the inputs `abb` and
`aba`. Its closure is complete at 47 species and 59 reactions. The first eight
reactions shown are Fok I binding (rule 3) to every molecule that contains
`GGATG`, then cutting the two inputs (rule 4). The answer is in
`net.extras["analysis"]`:

```python
net.extras["analysis"]
# {'reporters': {'abb': ['S0'], 'aba': []}, 'accepted': ['abb']}
```

`net.extras["reaction_rules"]` gives each reaction's rule number, as in the
papers, and `net.extras["rules"]` lists the rules themselves. `aba` leaves no
reporter because it ends in S1 and the published pool contains only the S0
detector. The closure also holds side reactions the rules allow, such as a
spent Fok I head ligating to a sticky end; they do not change which words are accepted.

**Other words.** Any words over a and b can be given; `s1_detector=True` adds the
S1 detector implied by rule (11), so rejected words leave a reporter too:

```python
net = chemart.generate_network("tominaga-stacked-strings",
        words=["", "a", "b", "ab", "bb", "bab", "abba", "babb"], s1_detector=True)
net.extras["analysis"]["reporters"]
# {'': ['S0'], 'a': ['S0'], 'b': ['S1'], 'ab': ['S1'], 'bb': ['S0'],
#  'bab': ['S0'], 'abba': ['S0'], 'babb': ['S1']}      182 species, 285 reactions
```

**The biochemical pathways.** `system="transcription"` gives the 2009
transcription model (30 species, 28 reactions), and its analysis reports the
expected messenger RNA and whether it is made:
`{'mrna': '0#CapCGCAAUGCUGAGCUAGUUUU/', 'mrna_found': True}`. The DNA is set
with `dna`. `system="fatty-acid-oxidation"` starts from a fatty acyl CoA with
`carbons` carbons (default 10, as in the paper; it must be even) and returns
`{'acyl_coa_carbons': [2, 4, 6, 8, 10], 'acetyl_coa': True}`: every shorter chain
and acetyl CoA are reached.

**Your own rules.** `system="custom"` takes `rules`, `pool`, `sources` and
`drains` as text. This is the `AB` example as a single sampled run:

```python
net = chemart.generate_network("tominaga-stacked-strings", system="custom",
    rules=["0#*1AB/ + 0#CD/ -> 0#*1AB/1#CD/",
           "0#*1AB/1#CD/ + 0#AB2*/ -> 0#*1ABAB2*/1#CD/",
           "0#*1ABAB2*/1#CD/ -> 0#*1ABAB2*/ + 0#CD/"],
    pool={"0#CD/": 1}, sources=["0#AB/"], drains=["0#ABABAB1*/"],
    method="soup", steps=200, seed=1)
net.extras["final_state"]      # {'0#AB/': 4, '0#AB/1#CD/': 1}
```

Over the 200 steps it built chains up to `ABABABABAB`, and the drain removed the
chains of three or more `AB`s. The same system is built in as
`system="ab-concatenation"`. Its closure, and that of
`system="adleman-hamiltonian-path"`, is infinite, so they stop at `max_species`
with status `truncated`. The Adleman closure is slow: 2,000 species (6 s here)
do not yet contain the answer molecule, 5,000 species (about a minute) do.

## Results

**Computational universality (2007).** Tominaga et al. show how to build, from
any deterministic or nondeterministic Turing machine, a system that emulates it
using rules with three reactants, and then give a mechanical conversion to rules
with one or two reactants only. The chemistry with unary and binary rules is
therefore as expressive as a Turing machine. Chemart does not implement the
construction, which is generic rather than a system of its own.

**Adleman's experiment (2007).** Adleman (1994) solved a Hamiltonian-path
problem, finding a route through a seven-node graph that visits every node once,
by letting DNA strands for nodes and edges hybridise (pair up with
complementary strands). The model gives each
ten-base block an element: node i is a lower strand `L_i1 L_i2`, and the edge
from i to j is an upper strand made of the complements of i's second half and
j's first half. Two rule families express hybridisation and ligation (joining
strand ends). The paper ran it with 100,000 objects of each kind, reduced from
the experiment's 3 × 10^13 strands per kind, and obtained molecules for the only
answer, the path 0 → 1 → 2 → 3 → 4 → 5 → 6. Chemart's tests check the paper's
intermediate molecules `0#U_22U_31/1#L_31L_32/` and
`0#U_22U_31U_32U_41/1#L_31L_32/`, and build the answer molecule by applying the
rules along the path. The closure is infinite because walks can revisit nodes.

**Benenson's DNA automaton (2007).** Written with the base sequences of the real
experiment, the model reproduces the automaton step by step, and in the authors'
simulator "detectors that indicate that the automaton accepted the input abb
were generated". Chemart reproduces it: the tests check the published input
molecules character for character, the Fok I cut at 9 and 13 bases, the unique
choice of transition molecule by each sticky end, and the reporters. The default
network accepts `abb` and not `aba`; with the S1 detector, eight further words
end in S0 exactly when they have an even number of b's; and a sampled run of the
published pool also accepts `abb` while keeping all 100 Fok I molecules.

**Transcription (2009).** Seven rules turn a DNA double strand into messenger
RNA. RNA polymerase binds just after the promoter `TATATT`, the sequence that marks
where transcription starts (standing in for a TATA box), a `Cap` element starts the RNA, rules (14)-(17) add the nucleotide
complementary to each DNA base (`U` opposite `A`), and on reaching `UUUU` the
polymerase lets go. Chemart's tests follow these rules base by base and check
that the paper's example releases `0#CapCGCAAUGCUGAGCUAGUUUU/`. The companion
models of DNA replication (22 rules) and translation (53 rules) are not in
Chemart: the paper prints only 8 and 6 of their rules, and the rest would have
to be invented.

**Fatty acid oxidation (2009).** Four rules model the four steps of the cycle
that shortens a fatty acyl CoA by two carbons, releasing one acetyl CoA and
producing FADH2, NADH and H+. Sequence wildcards stand for the rest of the
carbon chain, so the same four rules handle a chain of any length, and a C10
acyl CoA ends as five acetyl CoA. The paper keeps every atom except a dummy
vacancy element `X`, needed because lines must be continuous. Chemart's tests
check the byproducts of each step, that the C10 network reaches every shorter
chain and acetyl CoA, and that every reaction conserves atoms apart from `X`.

**Reasoning backwards (2009).** The 2009 paper adds a search procedure that asks
whether a given set of molecules can be produced, by applying rules in reverse
from the goal. It is sound, since a successful search yields a real reaction
path, but not complete, and the authors point out that no always-terminating
procedure could be both, because the chemistry can express Turing machines.
Asked which acyl CoA yields three acetyl CoA, it answers with a six-carbon
chain (their Figure 3). Chemart does not implement the procedure. Its tests
replay Figure 3 forwards and find that the molecule as printed lacks the
extra hydrogen at the start of the carbon line that the paper's own acyl CoA
form has: rules (25)-(28) reproduce every
intermediate of the figure, but the last step then gives `0#HO/0#CC/0#HSCoa/`
instead of a third acetyl CoA. With the corrected six-carbon molecule the path
gives three acetyl CoA.

**What the authors concluded.** Because elements have no built-in functions,
Fok I can be a plain run of `F`, and systems combine by taking the union of
their rules, which the 2009 paper
contrasts with Fontana and Buss's lambda-calculus chemistry (the basis of
[AlChemy](alchemy.md)), where changing one molecule's form changes its
function. The authors also list the limits: the models are
qualitative and "not necessarily biologically accurate", and a stack of strings
cannot represent rings, trees or hairpins, although they note that pathways
such as glycolysis, the citric acid cycle and the urea cycle can still be written
when the molecules have no variable parts.

## Further reading

- Tominaga, K. (2004). A formal model based on affinity among elements for
  describing behavior of complex systems. Technical report UIUCDCS-R-2004-2413,
  Department of Computer Science, University of Illinois at Urbana-Champaign.
  The first description of the chemistry.
- Adleman, L. M. (1994). Molecular computation of solutions to combinatorial
  problems. *Science* 266(5187), 1021–1024.
- Benenson, Y., Paz-Elizur, T., Adar, R., Keinan, E., Livneh, Z. & Shapiro, E.
  (2001). Programmable and autonomous computing machine made of biomolecules.
  *Nature* 414, 430–434.
