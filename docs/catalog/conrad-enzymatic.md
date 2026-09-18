# Conrad's lock-and-key enzymatic / self-assembly processor

`conrad-enzymatic` · *Conrad, 1985-1992; wet realisation by Zauner & Conrad, 2001*

*Also known as:* *enzymatic computing*, *conformation-based computing*, *lock-key paradigm*

Computation by molecular shape recognition. Ligands carry keys, an enzyme carries locks, and binding happens when shapes are complementary. Signal ions binding a modulator site switch the enzyme between conformations, each with different catalytic behaviour - so a single enzyme, in a milieu of two ion species, classifies input patterns by how much product it makes. Conrad's argument was that this is what molecules are good at and digital logic is not: recognition of a whole pattern in one step, rather than a program over bits.

| | |
|---|---|
| **family** | bio-inspired |
| **kind** | wet |
| **constructive** | no — fixed species set |
| **fidelity** | `reconstructed` — built from the original papers listed below |
| **book** | 11.2.1; also 17.1.2 (emergence vs. programmability tradeoff) |
| **refs** | [193] Conrad, On design principles for a molecular computer, CACM 28(5):464-480, 1985, doi:10.1145/3532.3533, [194] Conrad, Molecular computing: the lock-key paradigm, Computer 25(11):11-20, 1992, doi:10.1109/2.166400, [196] Conrad & Zauner, Conformation-based computing: a rationale and a recipe, in Sienko et al. (eds), Molecular Computing, MIT Press 2003, ch. 1, doi:10.7551/mitpress/4739.003.0003, [954] Zauner & Conrad, Enzymatic computing, Biotechnol. Prog. 17(3):553-559, 2001, doi:10.1021/bp010004n, [722] Rosen, Recent developments in the theory of control and regulation of cellular processes, Int. Rev. Cytology 23:25-88, 1968 |
| **provides** | `topology`, `stoichiometry`, `initial-state`, `mass-conservation`, `sequence-structure-function` |

## Molecules, reactions, reactor

**S — molecules** (explicit): Molecules are 6-bit shapes: ligands (Mg2+, Ca2+ signal carriers, NAD+, L-malate, and the released NADH and oxaloacetate) carry a key, malate dehydrogenase carries three locks (coenzyme, substrate and modulator sites). One species per enzyme conformer (set of ions bound at the modulator site) and per binding state of its two catalytic sites.

**R — reactions** (implicit, arity [1, 2]): lock-and-key: a ligand binds a site when at least match_bits of the 6 shape positions are complementary (bit-inverted). Signal ions bind the modulator site and select a conformer; every conformer binds NAD+ then L-malate reversibly and turns the ternary complex over to oxaloacetate + NADH + H+.

**A — reactor**: well-stirred-multiset
 · *dilution:* none

## What you get

```python
net = chemart.generate_network("conrad-enzymatic", seed=1)
```

```
conrad-enzymatic: 12 species, 12 reactions, status=complete
provides: initial-state, mass-conservation, stoichiometry, topology
seed: 1
extras: analysis, buffered, conservation, experiment, interaction_law
```

First reactions:

```
MDH + NAD -> MDH_NAD
MDH_NAD -> MDH + NAD
MDH_NAD + MAL -> MDH_NAD_MAL
MDH_NAD_MAL -> MDH_NAD + MAL
MDH_NAD_MAL -> MDH + OAA + NADH + H
MDH_Mg + NAD -> MDH_Mg_NAD
MDH_Mg_NAD -> MDH_Mg + NAD
MDH_Mg_NAD + MAL -> MDH_Mg_NAD_MAL
… and 4 more
```

## Parameters

| name | type | default | role | what it does |
|---|---|---|---|---|
| `operation` | `enum` | `XOR` | structural | which of the six commutative two-bit operations of Table 1 the milieu encoding is meant to realise; fixes the truth table, the active milieu states and the signal-strength formula in extras['analysis'] <br>one of `AND`, `OR`, `XOR`, `NAND`, `NOR`, `NXOR` |
| `inputs` | `enum` | `11` | population | the two-bit input pattern presented; each 1-signal adds one encoding unit of every carrier ion to the cuvette, giving milieu state a (00), b (01/10) or c (11) <br>one of `00`, `01`, `10`, `11` |
| `encoding_mM` | `dict` | `{'Mg': 72.38}` | population | concentration (mM) that one 1-signal contributes for each signal carrier; the keys choose the carriers, and hence which enzyme conformers exist. Valid keys: Mg, Ca <br>*range:* published encodings: {Mg: 20} realises OR, {Mg: 20, Ca: 40} realises XOR; the response surface spans 0-133 mM of each ion, i.e. encodings up to ~67 mM |
| `match_bits` | `int` | `5` | structural | minimum number of complementary positions (out of 6) for lock-and-key recognition; 6 admits only perfect complements and so drops Ca2+ as a signal carrier <br>`5` … `6` |
| `malate_mM` | `float` | `5.41` | population | L-malate in the cuvette (mM); 7.1 mM in 2 x 0.8 mL of signal solution diluted into 2.1 mL <br>≥ `0` |
| `nad_mM` | `float` | `1.29` | population | NAD+ in the cuvette (mM); 5.4 mM in 0.5 mL of enzyme solution diluted into 2.1 mL <br>≥ `0` |

## Published phenomena

What the literature reports this model produces. Whether the generator reproduces each one is recorded in the decisions below.

- lock-and-key recognition: only shape-complementary key/lock pairs bind, and Ca2+ is the weaker key of the two signal carriers
- context-sensitive pattern classification: a single enzyme classifies two-bit input patterns by its catalytic activity
- a nonmonotonic (convex) response to a single signal ion turns the linearly inseparable XOR into a linearly separable problem, solvable with one threshold
- AND, OR, XOR, NAND and NOR are all realisable at suitable encodings; NXOR is realisable at none, because a unimodal response has no interior minimum
- the tabletop device classified all 135 presented two-bit patterns correctly (45 x 00, 46 x 01/10, 44 x 11)

## Sources

- K.-P. Zauner and M. Conrad. Enzymatic Computing. Biotechnology Progress 17(3):553-559, 2001. doi:10.1021/bp010004n. Full text used: https://eprints.soton.ac.uk/259138/1/ZaunerKP01EnzComp.pdf
- W. Banzhaf and L. Yamamoto. Artificial Chemistries. MIT Press, 2015, sections 11.2.1 and 17.1.2 (Figure 11.10, adapted from Conrad 1992).
- M. Conrad. On design principles for a molecular computer. CACM 28(5):464-480, 1985. doi:10.1145/3532.3533 (paywalled; not obtained).
- M. Conrad. Molecular computing: the lock-key paradigm. Computer 25(11):11-20, 1992. doi:10.1109/2.166400 (paywalled; not obtained).
- M. Conrad and K.-P. Zauner. Conformation-based computing: a rationale and a recipe. In Sienko et al. (eds), Molecular Computing, MIT Press, 2003, ch. 1, pp. 1-31. doi:10.7551/mitpress/4739.003.0003 (paywalled; not obtained).

## Decisions

Every gap, ambiguity or erratum in the sources, and how Chemart resolved it. Read this before quoting a number from this entry.

- The book's Figure 11.10 processor is a schematic, so the generated network is the in-silico network of the wet realisation the book points to ([954]): the signal ions are the input signals, the ion-loaded conformer is the recognising complex, its active site is the adaptor that binds NAD+ and L-malate, and NADH is the output. The amplifying signalling cascade of the figure has no counterpart in the 2001 experiment and is not generated.
- The 6-bit keys and locks are a Chemart encoding of the lock-and-key criterion of book 11.2.1, not published data. NAD+, L-malate and Mg2+ are perfect complements of the coenzyme, substrate and modulator locks; Ca2+ is one bit off, which reproduces the paper's observation that the signal strength is stronger for MgCl2 alone than for CaCl2 alone; the released NADH and oxaloacetate are 2 bits off their locks and do not rebind at match_bits >= 5.
- The paper is deliberately phenomenological ('our purpose here was not to investigate mechanism'), so the mechanism is the minimal ordered scheme consistent with the assay: reversible coenzyme then substrate binding on every conformer, one irreversible turnover, and modulator-ion binding on the free enzyme only. No kinetics are published, so every reaction rate is None.
- Response levels are measured, not modelled: Figure 1 is an interpolated surface picture and Figure 9 a bar chart, so no absorbance value is reproduced. What the module computes is Table 1's classification law: the milieu states that must be active, and the signal strength Delta s = min(r over active states) - max(r over the rest), which reproduces all six of the table's formulas.
- Defaults are the tabletop device of Figures 8-9 with MgCl2 alone: 190 mM MgCl2 in a 0.8 mL 1-signal, 7.1 mM L-malate in each 0.8 mL signal and 5.4 mM NAD+ in 0.5 mL of enzyme solution give, in the 2.1 mL cuvette, 72.38 mM Mg2+ per 1-signal, 5.41 mM L-malate and 1.29 mM NAD+.
- The enzyme is dosed as 6 uL/mL of a commercial suspension, so its molar concentration is not published and it is absent from initial_state; the protocol string is kept in extras['experiment']. The proton released by the turnover is listed in extras['buffered'] (glycine/NaOH, pH 10).
- The Results sentence 'The activity first increases to a maximum and then decreases as the concentration of ion decreases' contradicts itself and Figure 1; it is read as '... as the concentration of ion increases', consistent with the surface being described as convex (strictly nonmonotonic) along both ion axes.
- Conrad's own papers [193], [194] and [196] are paywalled and were not obtained, so the conceptual reading (computation as downhill free-energy minimisation, programmability traded against evolvability) is taken from book 11.2.1 and 17.1.2 and no number is attributed to them. The v1 entry's origin 'Conrad, 1974-1992' is corrected to the works the bibliography actually cites, 1985 and 1992.
- Readout times are not parameters (they do not change the network); 300 s for the response surface and 10 s for the device are recorded in extras['experiment'].

## Notes

Wet: the network returned is the in-silico network of the published experiment, not a laboratory protocol. Its lasting contribution to Chemart is the emergence-vs-programmability tradeoff of 17.1.2: "a system cannot at the same time be effectively programmable, amenable to evolution by variation and selection, and computationally efficient."

---

*Specification: `catalog/chemistries/conrad-enzymatic.yaml` · generator: `chemart/chemistries/conrad_enzymatic.py` · tests: `tests/chemistries/test_conrad_enzymatic.py`*
