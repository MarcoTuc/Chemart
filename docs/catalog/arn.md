# Artificial Regulatory Network (ARN)

`arn` · *Banzhaf, 2003*

A gene regulatory network grown from one bitstring. Scanning the genome for a promoter pattern marks genes; each gene's coding region is collapsed by majority vote into a 32-bit protein, and the same gene carries enhancer and inhibitor sites. How strongly protein j regulates gene i is set purely by how complementary the protein is to that site - bit-matching stands in for molecular affinity. Concentrations then follow an exponential response on the simplex. Because the whole network is encoded in a sequence, mutation and duplication act on the regulatory topology itself.

| | |
|---|---|
| **family** | network |
| **kind** | generator |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.3.4; GRN modelling 18.4.3; morphogenesis 18.6.1 |
| **refs** | [67], [68], [73], [478], [479], [480], [504], [521], [522], [172] |
| **provides** | `topology`, `stoichiometry`, `catalysts`, `rate-constants`, `flow`, `initial-state`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (implicit): one bitstring genome; each gene i (found by scanning for the promoter 01010101) yields a protein P_i: a 32-bit string obtained by bitwise majority voting over its 5 x 32-bit coding region

**R — reactions** (implicit, arity 2): Layout of a gene: [inhibitor site 32][enhancer site 32][promoter 8][5 x 32 coding bits]. d^a_ij (d^h_ij) = popcount(protein_j XOR enhancer (inhibitor) site of gene i). Enhancement P_j + P_i -> P_j + 2 P_i, k = delta/N exp(beta (d^a_ij - d_max)); inhibition P_j + P_i -> P_j, k = delta/N exp(beta (d^h_ij - d_max)); 2 P_i on the left when j = i.

**A — reactor**: ode
 · *dilution:* concentrations kept on the simplex sum_i c_i = 1 (outflow constant-total)

## What you get

```python
net = chemart.generate_network("arn", seed=1)
```

```
arn: 7 species, 98 reactions, status=complete
provides: catalysts, flow, initial-state, rate-constants, stoichiometry, topology
seed: 1
extras: analysis, genome
```

First reactions:

```
2 P1 -> 3 P1  [mass-action k=2.1757113921018042e-09 site=enhancer match=14]
2 P1 -> P1  [mass-action k=1.1878981701479541e-07 site=inhibitor match=18]
P2 + P1 -> P2 + 2 P1  [mass-action k=2.1757113921018042e-09 site=enhancer match=14]
P2 + P1 -> P2  [mass-action k=1.1878981701479541e-07 site=inhibitor match=18]
P3 + P1 -> P3 + 2 P1  [mass-action k=2.9445051749122256e-10 site=enhancer match=12]
P3 + P1 -> P3  [mass-action k=1.607645353132273e-08 site=inhibitor match=16]
P4 + P1 -> P4 + 2 P1  [mass-action k=4.3700331500260826e-08 site=enhancer match=17]
P4 + P1 -> P4  [mass-action k=1.0832229182731294e-10 site=inhibitor match=11]
… and 90 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `genome_source` | `enum` | `random` | structural | random: every bit drawn at random (Banzhaf 2003). duplication-divergence: a random 32-bit string doubled repeatedly, each doubling followed by point mutations (Kuo & Banzhaf 2004); this is what yields scale-free gene-number distributions and small-world topologies <br>one of `random`, `duplication-divergence` |
| `genome_length` | `int` | `4096` | structural | genome length L_G in bits; for duplication-divergence it must be word_bits * 2^k, k being the number of duplications <br>`1` … `1048576` · *range:* Banzhaf 2003 uses 1000, 10000, 100000; Kuo & Banzhaf use 131072 (12 duplications) and 32768 (10 duplications). A random genome carries about genome_length/340 genes |
| `mutation_rate` | `float` | `0.01` | structural | probability that each bit of the genome flips after each whole-genome duplication (duplication-divergence only) <br>`0` … `1` · *range:* Kuo & Banzhaf use 0.001, 0.01 and 0.05 |
| `promoter` | `str` | `01010101` | structural | bit pattern that marks the start of a gene; overlapping matches and periodic extensions (0101010101) count as one promoter at the first bit |
| `word_bits` | `int` | `32` | structural | length of a regulatory site, of each of the 5 coding segments and of a protein; also the maximum match d_max <br>`1` … `64` |
| `beta` | `float` | `1.0` | kinetic | positive scaling of the match in the exponential; the RT analogue <br>≥ `0` |
| `delta` | `float` | `1.0` | kinetic | positive overall rate scaling (time unit) <br>≥ `0` |
| `link_threshold` | `int` | `` | structural | keep only regulatory interactions with at least this many complementary bits; 0 keeps every interaction, i.e. the full ODE model <br>`0` … `64` · *range:* Kuo & Banzhaf draw graphs at 21 and 22 |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- random genomes of length 131072 carry 340-440 genes (Kuo & Banzhaf 2004 Fig. 5)
- duplication-divergence at 1% mutation gives a broad power-law-like gene-number distribution (0-3000 genes at L_G = 131072, Fig. 3); at 5% it resembles random genomes (200-700, Fig. 4)
- sharp transition from full to no connectivity as the link threshold rises from 0 to 32 (Kuo et al. 2006 Fig. 4)
- small-world and scale-free topologies under duplication-divergence
- from equal concentrations, a few proteins dominate and the dynamics settles into point attractors or damped oscillations (Banzhaf 2003 Figs. 3-4)
- heterochrony: one- or two-bit changes in a regulatory site shift expression onset and offset (Banzhaf 2003 Figs. 5-6)
- network motifs matching natural GRNs
- evolvable to target concentration dynamics (sinusoid, exponential, sigmoid)

## Sources

- Banzhaf, W. (2003). On the dynamics of an artificial regulatory network. Advances in Artificial Life, ECAL 2003, LNAI 2801, pp. 217-227 (book [68]). Promoter, gene length, XOR matching, eqs. 1-3, Table 1, Figs. 3-6. http://www.cs.mun.ca/~banzhaf/papers/ecal2003_final.pdf
- Kuo, P. D. & Banzhaf, W. (2004). Small world and scale-free network topologies in an artificial regulatory network model. Artificial Life IX, pp. 404-409 (book [478]). Promoter rules, duplication-divergence, gene-number histograms Figs. 3-5. http://www.cs.mun.ca/~banzhaf/papers/Genome_final.pdf
- Kuo, P. D., Leier, A. & Banzhaf, W. (2004). Evolving dynamics in an artificial regulatory network model. PPSN VIII, LNCS 3242, pp. 571-580 (book [480]). Fig. 1 gene layout, eqs. 1-2 with u_max = 32, initial concentrations 1/N. http://www.cs.mun.ca/~banzhaf/papers/ARNOptimizeNew.pdf
- Kuo, P. D., Banzhaf, W. & Leier, A. (2006). Network topology and the evolution of dynamics in an artificial genetic regulatory network model created by whole genome duplication and divergence. BioSystems 85:177-200 (book [479]). Section 2, Fig. 4 (edges vs threshold), L_G = 2^12 x 32. https://www.cs.mun.ca/~banzhaf/papers/biosystems85_2006_177.pdf

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- Gene layout follows Kuo, Leier & Banzhaf (2004) Fig. 1: inhibitor site, then enhancer site, immediately upstream of the promoter; the 160 coding bits start right after the promoter, at any bit offset. Banzhaf (2003) instead starts the gene 'at the next integer'. Genes whose sites or coding region fall outside the genome are dropped; genes may overlap and are treated independently (Kuo et al. 2006).
- d is the number of complementary bits (XOR popcount), as in the book's definition and all papers; the book's later phrase 'decreasing exponential on the number of matching (equal) bits' and its claim that beta (d - d_max) is 'zero otherwise' are slips: the exponent is 0 only at d = d_max, and influence grows with complementarity.
- Overlapping promoter matches collapse into the first one: a match is a new promoter only if no match starts in the preceding len(promoter) - 1 bits, which realises 'overlapping promoters or periodic extensions are not allowed'.
- Eq. 11.12 divides by sum_j c_j but does not by itself keep the sum at 1; the papers impose it with a flow term Phi (Banzhaf 2003 eq. 3). On that simplex each ODE term is one mass-action reaction (enhancement P_j + P_i -> P_j + 2 P_i, inhibition P_j + P_i -> P_j, k = delta/N exp(beta (d - d_max))) and Phi is exactly outflow = constant-total; tested at equation level.
- d_max = word_bits = 32, as in the book and Kuo et al. Banzhaf (2003) eqs. 1-2 instead rescale u_max to the best match present; that variant is not implemented.
- No source gives numeric beta or delta; defaults are beta = delta = 1. Initial concentrations are equal, 1/N (Banzhaf 2003 'equal concentration'; Kuo et al. 2004 '1/#genes').
- Duplication-divergence: a random word_bits string is doubled k times and after each doubling every bit of the whole genome flips with probability mutation_rate (consistent with Kuo & Banzhaf's 66% survival of an 8-bit promoter per duplication at 5%: 0.95^8 = 0.66). genome_length replaces the duplication count (L_G = 32 x 2^k), and segments_per_gene is fixed at the papers' 5.
- Default genome: random, 4096 bits (about 11 genes, a few hundred reactions). Small duplication-divergence genomes often carry no gene at all, which gives an empty network. link_threshold (the book's static-topology threshold, an edge when d >= threshold, so threshold 0 is fully connected as in Kuo et al. 2006 Fig. 4) drops weaker interactions and thereby changes the dynamics.
- The evolution experiments (fitness on target dynamics, output sites of Kuo et al. 2004) and the morphogenesis extension of [172] act on genomes between networks and are not generated. The whole genome is kept in extras.genome.

## Notes

The book presents the ARN dynamics as a normalised exponential-response ODE rather than a reaction network; on the simplex it is exactly mass action with constant-total dilution, and the bitstring match prescribes every rate constant.

---

*Specification: `catalog/chemistries/arn.yaml` · generator: `chemart/chemistries/arn.py` · tests: `tests/chemistries/test_arn.py`*
