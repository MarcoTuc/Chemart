## Introduction

SAC, the *string-based artificial chemistry* of Hideaki Suzuki and Naoaki Ono
(2002–2003), is a chemistry in which every molecule is a short piece of text
that can also be run as a program. When two strings meet, the first is read as
a set of search-and-replace rules and applied to the second. The first string
comes out of the collision unchanged, so it acts as a catalyst, an enzyme; the
second is rewritten and possibly cut into pieces. Some strings are written so
that their rules do nothing: they are inert data, the *genes*. Others are
active: the *enzymes* that copy and translate those genes.

Suzuki and Ono used this language to build a working replication system out of
twelve strings: six genes and six enzymes. Three enzymes form a *copier*, which
makes an exact copy of any gene, and three form a *constructor*, which reads a
gene and writes out the enzyme it encodes. Since the six genes encode exactly
these six enzymes, the twelve strings together make more of themselves.
Because the constructor can translate any gene whatsoever, the authors call it
universal.

The question they then asked was biological. A primitive cell, or
*protocell*, keeps its genetic information in many separate replicating
molecules, and those replicate at different speeds. When the cell divides and
its molecules are shared out at random, a slow replicator can be missing from a
daughter, which is then defective. Suzuki's slides call this *segregational
instability*. With SAC strings enclosed in cells that grow and divide, Suzuki
(2003) compared three ways of organising the genome: six independent genes, one
chromosome carrying all the genes, and independent genes whose replication is
regulated and whose copies are pulled apart by a spindle before division.

Banzhaf and Yamamoto describe SAC in their chapter on bio-inspired chemistries
(book §11.1.3), next to two other string chemistries in which strings act on
strings: the [Molecular Classifier System](mcs-bl.md) and
[Stringmol](stringmol.md). What sets SAC apart is its explicit gene/protein
split, inherited from molecular biology: a gene is its protein with every
character switched off. SAC is a simulation model, and Chemart implements the
string chemistry inside a single cell; the population of dividing cells is
described but not simulated (see *Results*).

## How it works

### Strings as programs

A SAC string uses 19 symbols (book table 11.4). The digits `0 1 2 3` are plain
data. The others give the string its meaning as a program:

- `'x` deletes the next symbol and `"x` creates it; any other symbol is simply
  matched and kept.
- `!`, `?` and `%` match any single character, and `*` any substring (the
  shortest one that fits).
- `/A/B/` swaps two substrings.
- `\x` *suppresses* x: it matches the character x literally, whatever its usual
  function.
- `&` and `$` separate rules that run one after another. After a rule, `&`
  goes on to the next rule only if the rule matched, `$` only if it did not.
- `.` is a cut mark: after rewriting, the operand is cut at every `.` that is
  not suppressed.
- `M` and `E` are membrane material, and `L` and `R` are spindle tags; both
  appear only in the third model.

The book's example is the string `\0&\!$0'0"1*0'1"2`, which reads as three
rules, written `[\0]&[\!]$[00*01 → 01*02]`. The first checks for a `0`, the
second for a `!`; only if the string contains a `0` and no `!` does the third
rule run. It replaces the second character after a leading `00` by `1`, and the
`1` of the next `01` by `2`. In Chemart:

```
\0&\!$0'0"1*0'1"2  +  002301  ->  \0&\!$0'0"1*0'1"2  +  012302
```

The first string, the *operator*, is put back as it was; the second, the
*operand*, is rewritten. If the rules change nothing, the collision is elastic
and nothing happens. The string printed in the book lacks the `&` and the `\`
before `!`, and does not work as described; Chemart follows the spelling in Suzuki's slides (see the
decisions below).

### Genes, enzymes and markers

A gene is the enzyme it encodes with a `\` in front of every character, framed
by the start marker `0000` and the terminator `0011`. Suppression makes the
gene inert: read as an operator it is a single rule that only matches, so it
never rewrites anything. Translation is the removal of the backslashes.

The enzymes work by walking a marker through a gene. Here is a real copy, made
by applying the three copier strings of the default cell to the shortest gene
(50 characters, which encodes one of the constructor's strings). Decoded, the
three copier strings are:

```
[0000*0011 → 00110202*00110202]        start: plant two markers 0202
[0202\%*0202 → \%0202*\%0202]          move: step one character, copy it to the end
[0011*02020011*0202 → 0000*0011.0000*0011]   finish: restore the frames, cut
```

The run:

```
0   0000\/\\\0 ... \2\2\00011
 -> 00110202\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\000110202
1   -> 0011\/0202\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011\/0202
2   -> 0011\/\\0202\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011\/\\0202
...
21  -> 0011\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\002020011\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00202
22  -> 0000\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011
     + 0000\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011
```

In step 0 the first copier string turns the start marker into `0011` and plants
the marker `0202` twice: in front of the gene body and after it. In each of
steps 1 to 21 the second copier string takes the suppressed character just
after the first marker (a `\` followed by any character, matched by `\\%`),
moves it in front of the marker and writes a copy of it before the second
marker. The copy grows at the tail, one character per collision. When the first
marker reaches the terminator, the third string's pattern fits at last: it
restores the frames `0000 … 0011` and inserts a `.`, and the cut releases two
identical genes. Every step here is a separate collision with a copier string,
and the copier strings are never changed.

The constructor works the same way with the marker `0220`, but its moving
string drops the backslash as it copies. Its last step on the same gene gives
back the gene and the enzyme it encodes:

```
22  -> 0000\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011  +  /\0\2\2\0/\\%/*"%0220
```

A string matches a pattern leftmost and with the shortest possible `*`. That is
what stops a marker at the next character rather than the last one; the
sources do not spell this out, and Chemart chose it because greedy matching
breaks the published cycles.

### Three models

- **Model (i)**, independent genes: the twelve strings above. Every gene
  replicates on its own, and the genes have different lengths (50 to 110
  characters in Chemart's strings), so they copy at different speeds.
- **Model (ii)**, one chromosome: all eight genes joined into one string by
  the separator `01010110`, closed by the terminator `3300`, with a copier and
  a constructor of four strings each. The whole genome is copied in one pass.
- **Model (iii)**, spindle and membrane: seven independent genes, a copier of
  four strings and a constructor of three. A finished copy comes out as two
  genes tagged `L…` and `R…`. The copier refuses to copy a tagged gene, and at
  division the tags send one copy to each daughter. The last copying step also
  releases the membrane seed `EM` and the breeder `"M\E\M`, which decodes to
  `[EM → MEM]`, so every collision between them adds one `M`.

### The cell level

In the papers each set of strings lives in a cell. Strings meet in random
pairs, the first as operator. A cell divides when its number of strings has
doubled (in model iii, only once a membrane string with ten `M`, that is
`MMMMMMMMMMEM`, exists), and its strings are shared equally between the
daughters. The population holds at most 100 cells of at most 50 strings, and
each new cell replaces an old one. None of this is part of the reaction network
that Chemart builds: it acts on whole cells, not on reactions between two
strings.

## Using it

The default call builds the *closure* of model (i): starting from the twelve
strings, Chemart applies every string to every other string, adds whatever
comes out, and repeats until nothing new appears. The result is complete, with
410 species and 410 reactions, and takes under a second. The species are all
the intermediate states of every copy and translation. Species names spell the
string with letters (`\` is `b`, `'` is `d`, `"` is `c`, `*` is `s`, `/` is
`p`, `%` is `z`, `.` is `t`; the full table is in the specification above), so
`s_00d0d0c1c1c0c2c0c2s0011c0c2c0c2` is the copier's first string. Each species'
`structure` holds the raw string. In `net.extras`, `roles` names the gene,
copier and constructor strings of the seed, and
`analysis["reproduced_seed"]` lists the seed strings that some reaction
produces; in the default run it holds all twelve.

**A cell as a soup.** `method="soup"` simulates random collisions in one
growing cell instead:

```python
net = chemart.generate_network("sac", seed=3, method="soup", steps=20000)
a = net.extras["analysis"]
a["population_size"][0], a["population_size"][-1]   # (12, 34)
sum(r.count for r in net.reactions)                  # 850 effective collisions
```

This took about two seconds. `population_size` records the number of strings
after every `chunk_steps` collisions (12 here, the initial size); it never
falls, since no reaction consumes a string. `final_state` holds the final
contents. At the end of this run none of the six genes was in its resting form
`0000…0011`: each was part-way through a copy or a translation. With
`dilution="constant"`, random strings are removed to keep the cell at its
initial size, a Chemart stand-in for the papers' size limit.

**The other models.** `organisation="single-chromosome"` or
`"spindle-membrane"` seeds the ancestral cell of model (ii) or (iii). Their
closures never end, so they always stop at `max_species`: with
`max_species=300`, model (iii) gives 563 reactions in 0.6 s and model (ii) 397
reactions in 1.8 s. For model (iii), `analysis["max_membrane_M"]` reports the
longest membrane string.

**Your own strings.** `strings` replaces the ancestral cell. The book's example
closes after one reaction:

```python
net = chemart.generate_network("sac", strings=[r"""\0&\!$0'0"1*0'1"2""", "002301"])
net.summary()        # sac: 3 species, 1 reactions, status=complete
```

To try the language directly, `chemart.chemistries.sac` exports `decode`,
`rewrite`, `react` and `ancestral_cell`.

## Results

**Catalysis, and a universal constructor.** Because the operator always comes
back unchanged, "since the operations are catalytic, a system of this type will
produce more and more of itself" (book §11.1.3). In Chemart every reaction of
the default closure has its operator as a catalyst, and the soup of model (i)
grows from 12 strings to about three times that in 20,000 collisions (the tests
check at least 24). Suzuki's slides state that the constructor "is 'universal'
in that it can create any operator string". Chemart's tests check that the
model (i) copier turns each of the six genes into two exact copies, that the
constructor turns each gene into itself plus its enzyme, and that it also
translates a gene for the book's example string, which is not part of the
system. The closure of model (i) produces all twelve seed strings, including
reactions that release two copies of a gene.

**Segregational instability.** Suzuki (2003) ran each model in a population of
dividing cells. The book's figure 11.9 and the slides show the outcome: in
model (i) the initial growth of replicators comes to a halt and the cell
population stagnates. The slides conclude that cellular selection alone
"is not powerful enough to conserve genetic information": six chromosomes "with
diverse lengths (50 to 110) were not conserved on account of segregational
instability". Models (ii) and (iii) "succeeded in conserving genetic information
stably". The slides' answer is that the genome "can be conserved if all genes
are linked on the same chromosomes (like prokaryotes) or if the replication of
genes are strictly regulated synchronously with the cell division cycle (like
eukaryotes)". Chemart does not reproduce this result,
because it does not simulate cells, division or selection; the published
settings are recorded in `net.extras["cell_level"]`.

**The single chromosome.** The tests check that the model (ii) copier copies the
whole chromosome exactly, and that the constructor releases all eight enzymes
from it. There is a catch in the published strings: at the end of a copy, the
copier's single-character mover (needed to step over the unsuppressed gene
separators) can also carry the marker past the terminator before the cutting
rule acts, and the copy then never closes. The copy is exact only if the cutting
rule acts first. Chemart gives no rule priority, so its closures and soups of
model (ii) include these runaway copies; whether the paper prevents them is
unknown.

**Regulated replication and the membrane.** For model (iii) the tests check the
mechanism of the slides: copying a gene yields the tagged copies `L…` and `R…`
plus the membrane seed and breeder, the copier leaves tagged genes alone, and
ten collisions of the breeder grow `EM` into `MMMMMMMMMMEM`, the ten-`M`
membrane that permits division. In a soup, however, the model (iii) cell does
not get that far. Two of the copier's by-products, `''M"\'"\'MMMM` and
`''\""\'"\'M\\\E\\\M`, contain doubled prefixes whose meaning no accessible
source explains; in Chemart's reading they react with nothing, and one copier
string releases the first of them from every gene it meets. With seed 1, the
cell holds 79 strings after 20,000 collisions, 60 of them this by-product, and
no membrane has grown; after 100,000 collisions it holds 124 strings, 105 of
them the by-product, and no tagged copy has yet appeared. This is a limit of
the reconstruction, not a published result.

**What the sources leave open.** Neither journal paper (Ono and Suzuki 2002;
Suzuki and Ono 2002) nor the text of Suzuki (2003) could be obtained. Every
string in Chemart comes from Suzuki's ECAL 2003 slides, transcribed from their
images, and the decoding rules were reconstructed so that all the published
examples and all three ancestral cells work. The decisions above list each
choice.

## Further reading

- Ono, N. & Suzuki, H. (2002). String-based artificial chemistry that allows
  maintenance of different types of self-replicators. *Journal of Three
  Dimensional Images* 16, 148–153. (Book ref. [642]; not consulted.)
- Suzuki, H. & Ono, N. (2002). Universal replication in a string-based
  artificial chemistry system. *Journal of Three Dimensional Images* 16(4),
  154–159. (Book ref. [826]; not consulted.)
