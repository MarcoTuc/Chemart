# Automata reaction (32-bit binary string chemistry)

`automata-reaction` · *Dittrich & Banzhaf, 1998*

*Also known as:* *Dittrich-Banzhaf automata chemistry*, *constructive binary string system*

The cheapest constructive chemistry that still does something interesting. Molecules are 32-bit words, and a collision treats the first word as a finite automaton that consumes the second and emits a third: s1 + s2 -> s1 + s2 + A_s1(s2). Both reactants survive, so every string is both machine and tape. About a third of random strings turn out to be passive replicators, products tend to resemble their parents, and none of it is designed in - it falls out of 32 bits and one rule.

| | |
|---|---|
| **family** | automata |
| **kind** | generator |
| **constructive** | yes — the species set grows at run time |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 10.6.1 |
| **refs** | [233], doi:10.1162/106454698568521 |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): 32-bit words; species id w<8 hex digits>, structure the 32-character bitstring (most significant bit first)

**R — reactions** (implicit, arity 2): s1 + s2 -> s1 + s2 + s3 with s3 = A_s1(s2); elastic if the filter rejects s3 (paper eq. 3)

**A — reactor**: well-stirred-multiset
 · *dilution:* each inserted product replaces a random molecule; M constant (paper's reactor algorithm)

## What you get

```python
net = chemart.generate_network("automata-reaction", seed=1)
```

```
automata-reaction: 4818 species, 9872 reactions, status=observed
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: analysis, final_state
```

First reactions:

```
wb7362806 + w13410911 -> wb7362806 + w13410911 + w13410910  [mass-action k=1.0]  (x1)
w54694c3a + w8b89e43c -> w54694c3a + w8b89e43c + w8b89e43f  [mass-action k=1.0]  (x1)
w262802e1 + w8b169674 -> w262802e1 + 2 w8b169674  [mass-action k=1.0]  (x1)
w3700032a + w89d94737 -> w3700032a + w89d94737 + w89d9473b  [mass-action k=1.0]  (x1)
wedb4ccb2 + w3407a349 -> wedb4ccb2 + w3407a349 + w3407a34b  [mass-action k=1.0]  (x1)
w27d6c977 + w0e4da006 -> w27d6c977 + w0e4da006 + w0e4da005  [mass-action k=1.0]  (x1)
we9eff343 + w4295bf0f -> we9eff343 + w4295bf0f + w7295bf0f  [mass-action k=1.0]  (x1)
wc9f97bcf + wf8c7532f -> wc9f97bcf + wf8c7532f + wf8c75369  [mass-action k=1.0]  (x1)
… and 9864 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `method` | `enum` | `soup` | structural | soup: the paper's reactor algorithm, observed reactions with firing counts; closure: every reaction reachable from the distinct seed words (Chemart addition, e.g. to check that a published organization is closed) <br>one of `soup`, `closure` |
| `mechanism` | `enum` | `automata` | structural | automata: the automata reaction (paper 3.2); and: the paper's reference reaction s3 = s1 AND s2 (paper 3.1, fig. 2) <br>one of `automata`, `and` |
| `code_table` | `enum` | `1` | structural | instruction table (paper fig. 1, right): table 1 maps 1010 to NOT, table 2 maps it to EQ; the other 15 codes are identical <br>one of `1`, `2` · *range:* paper figs. 3-4 use table 1, figs. 5-10 (evolution) use table 2 |
| `forbid_exact_replication` | `bool` | `` | selection | filter f1 of paper eq. 3: a collision whose product equals one of its reactants is elastic, so exact replication is disabled <br>*range:* paper section 5.4 (evolution): true |
| `M` | `int` | `1000` | population | soup size: number of random 32-bit words in the initial soup (ignored when words is given) <br>`2` … `1000000` · *range:* paper: 100 (fig. 3), 10^4 (figs. 2, 4), 10^5 (figs. 5-6), 10^6 (fig. 7) |
| `generations` | `int` | `10` | population | soup only: run length in generations of M collisions each, elastic collisions included <br>`0` … `100000` · *range:* paper: 10 (fig. 2), 140 (fig. 3), 280 (fig. 4), 1000 and 7000 (figs. 5-6) |
| `words` | `list` | `` | structural | explicit initial multiset (soup) or seed set (closure) of 32-bit words as hex strings or integers; overrides the random draw of M words <br>*range:* e.g. the fig. 3 organization [7240a7ef, 7240a7ea, 7240a7eb, 7240a7ee] |
| `max_species` | `int` | `200` | structural | closure only: species budget; the closure of random words is usually cut off by it (status truncated) <br>≥ `1` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- many random strings are passive replicators (s1 + s2 => s2; paper 3.2: about 30%, measured here 20-38% depending on how many operands are checked) and very few are active self-replicators (s1 + s2 => s1; paper: about 0.004%, measured here about 0.005%)
- products resemble their parents; variation happens mostly at the right edge of the operand (paper table 2)
- AND reaction: extinction, the soup is taken over by 00000000 (paper fig. 2, M = 10^4)
- M = 100: diversity and innovativity drop until a closed self-maintaining organization of a few strings dominates, then drift (fig. 3: 7240a7ef, 7240a7ea, 7240a7eb, 7240a7ee)
- M = 10^4: an explorative phase of high diversity and innovativity ends (about generation 40) when active replicators appear and take over (fig. 4)
- M = 10^5-10^6 with code table 2 and exact replication forbidden: self-evolution without mutation operator or fitness function, punctuated quasi-stable organizations of closely coupled similar strings (figs. 5-9)
- spontaneous emergence of cross-over: 1e64a24e + 00000000 => 0000024e (fig. 10)
- syntactic and semantic closure of emergent organizations (confirming Fontana); lazy -> eager replicator transitions (book 10.6.1)

## Sources

- Dittrich, P. & Banzhaf, W. (1998). Self-evolution in a constructive binary string system. Artificial Life 4(2):203-220. Preprint (draft April 20, 1998): reactor algorithm (sec. 2), eqs. 1-3, automata reaction (sec. 3.2, fig. 1, tables 1-2), figs. 2-10, appendix (registers and instructions). https://users.fmi.uni-jena.de/~dittrich/p/DB97alife.ps
- Dittrich, P. (1997). autoreac-1.0, ANSI C source code of the automata reaction (paper ref. [8], named there as the precise formal specification): autoReac.c runN, instructionTable01/02. https://users.fmi.uni-jena.de/~dittrich/software/autoreac-1.0.tar

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book gives only the idea (10.6.1). The machine is a line-by-line port of autoReac.c, which the paper's appendix names as the formal specification; compiled with 32-bit words it and the port agree on 2015 random pairs for each code table, and both reproduce paper table 2, the fig. 3, 4 and 10 reaction tables and the fig. 9 block (tests).
- Where the paper's appendix and the C source disagree the source wins: TMM cycles both -> operator pointer only -> IO pointer only -> both (autoReac.c movMode++), whereas the appendix prints both -> IO pointer -> operator pointer.
- C details the paper does not spell out, kept as in autoReac.c: pointers are realised by rotating the registers, and the IO register is rotated back at the end (correctResult); pointers start at bit 0 (least significant) and 'left' moves to higher bits; instructions are read from the least significant nibble; a SETP in the last nibble gets pattern 0000; a logic instruction sets the last-ALU register, runs, and then makes one pointer step, which with copy mode on applies that ALU operation once more before moving (so EXOR/NOT/EQ cancel on that bit); CPON keeps the last ALU operation (initially ID); a MOV with a pattern checks the operator register, or the IO register in IO-pointer-only mode, after each of at most 32 steps.
- autoReac.c also contains code tables 3 and 4, a binding-preference option and a no-correction switch; the paper uses none of them (only tables 1 and 2), so they are not exposed. The v1 params word_len and opcode_width are dropped: the machine is defined only for 32-bit words and 4-bit instructions.
- The v1 topology param (2d-grid) and the 'space' capability are dropped: the spatial variant is Banzhaf, Dittrich & Eller (1999, Physica D 125:85-104), not [233], and lattice-2d is removed from A.reactor.
- Reactor: chemart.soup.soup with dilution constant runs generations x M collisions. The paper replaces a random one of the M molecules by s3; soup removes a random one of the M + 1 molecules after insertion, so s3 itself is removed with probability 1/(M + 1). soup also draws two distinct molecules, while the paper does not say whether s1 and s2 may be the same molecule. Both differences are O(1/M).
- Rates: paper eq. 1 (catalytic network equation) with eq. 2 (k = 1 if s1 + s2 => s3 exists and passes the filter, 0 otherwise). Each network reaction is a multiset {a, b} -> {a, b, s3}; its mass-action k counts the ordered pairs (a, b), (b, a) that produce s3 (1 or 2; 1 when a = b). outflow constant-total is the dilution term of eq. 1.
- Observed network: species are all initial words and all products; initial_state is the random soup; extras.final_state the final one (most frequent first); extras.analysis records the paper's macroscopic measures per generation: diversity Div = distinct types / M (from generation 0), productivity (inserted products / M) and total innovativity (never-seen products / M). The paper's DDC is not computed.
- mechanism and is the paper's reference reaction (3.1), included for fig. 2; code_table does not apply to it. method closure (chemart.expand.expand over ordered pairs) is a Chemart addition.
- Paper 3.2 says about 30% of random strings are passive and about 0.004% active self-replicators, without saying how many operands were tried. With the (C-verified) machine, 38% of random strings copy one random operand, 20-25% copy each of three, and about 0.005% (54-58 per 10^6) return themselves for eight operands; the test only checks the passive fraction is of that order.
- Paper erratum: the fig. 4 program listing prints the binary code of CPON as 0100 (that is TDIR); the string 1e1ca260 has nibble 0110 there, and the listed mnemonics are what the machine runs.

## Notes

The best cost/benefit constructive chemistry in the book: 32-bit words, no loops, a bounded number of steps per reaction, large populations, and it reproduces the qualitative results of AlChemy at a fraction of the compute. The observed network of a short soup run is a sample of a 2^32-species chemistry; paper-scale runs (M = 10^6) need a vectorised or compiled machine.

---

*Specification: `catalog/chemistries/automata-reaction.yaml` · generator: `chemart/chemistries/automata_reaction.py` · tests: `tests/chemistries/test_automata_reaction.py`*
