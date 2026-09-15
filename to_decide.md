# To decide

Open questions parked for later. Each one blocks a small part of the
catalog, never the library core.

## 1. Generative vs algorithmic: hybrid rulings

**Criterion (agreed 2026-09-06).** Is the chemistry a *fixed reaction rule*
applied to inert elements of the set (**generative**, e.g. prime-number
chemistry, Gamma, P systems), or is the reaction rule *encoded in the
molecule's own symbolic structure*, so the result is computed on
interaction (**algorithmic**, e.g. AlChemy / lambda chemistry, matrix
chemistry, Typogenetics)? This is the book's "reflexive" notion
(§11.1.1: a string "can act as both an operator and an operand"). It is
orthogonal to `constructive`.

Provisional split: algorithmic 24 · hybrid 8 · generative 66.

Algorithmic: `alchemy` `combinator-chemistry` `reflexive-ac`
`matrix-chemistry` `automata-reaction` `typogenetics` `stringmol` `mcs-bl`
`sac` `laing-molecular-machines` `mccaskill-polymer-tm` `ikegami-hashimoto`
`fraglets` `nac` `proof-ac` `urdar` `corewar` `coreworld` `tierra` `avida`
`high-order-chem` `dna-automaton` `bondable-ca`† `rbn`†

† `subsymbolic` sub-tag (§10.7.3): the operator emerges from the atom's
internal dynamics rather than being read from a sequence.

### Hybrids: the pattern

In most hybrids the *reaction* is a fixed external rule, but *which*
reaction fires is encoded in a molecule acting as an enzyme.

| entry | ambiguity | leaning |
|---|---|---|
| `squirm3` | Base rules R1–R8 are an external table, but the catalyst atom's state encodes the whole reaction via `i = 2(2(\|T\|(…)+b₁)+b₂+S`, read off a base-4 genetic string | ? |
| `smn` | Ligation/cleavage/recombination are fixed; the enzyme selecting them is encoded in the genome | ? |
| `bnc-cell` | Same as `smn` | ? |
| `rna-folding-ac` | Sequence → fold → catalytic function is a genuinely encoded operator, but the product is computed by external ToyChem rewriting | algorithmic |
| `arn` | Complementarity binding and the response law are fixed; the interaction *network* is encoded in the genome | ? |
| `aevol` | Gene → protein → mathematical function expressing metabolic activity | out of scope (framework) |
| `molecular-tsp` | Machines are molecules operating on strings, but their behaviour is hardwired by type (E/C/I/R), not read from a sequence | generative |
| `acgp` | Named "Algorithmic Chemistry GP", but registers are inert and the instruction multiset is external | generative |

## 2. Sub-tags inside each class (proposed, not agreed)

- `rule_supplied` inside generative-with-structure: Gamma, CHAM, ARMS,
  P systems, MGS, κ-calculus, Brane calculi, L-systems take a rule set as
  their primary argument. Also possibly `chemical-casting-model`, `soas`.
- `sampled` vs `fixed` inside generative: `random-catalytic-networks`,
  `bigan-conservative-crn`, `jain-krishna` draw a topology from a
  distribution; `brusselator`, `repressilator` are written down.

## 3. Low-confidence reconstructions to revisit

These generators work, are tested against every published fact that could be
found, and record their choices in `decisions`. Their key semantics rest on
inference, though, because the primary source was unavailable. Upgrade them if
the sources turn up.

| id | what is uncertain | what would settle it |
|---|---|---|
| `mccaskill-polymer-tm` | Rule framing, head positions and initial state were found by searching for the one decoding that makes the report's published replicator work. The report's evolution run (replicator survives, parasite extinct) is **not** reproduced. | McCaskill's original C program, or the later NGEN / Breyer-Ackermann-McCaskill papers in full |
| `sac` | String-rewriting semantics (`*` shortest match, prefixes, `&` / `$` sequencing) are inferred from Suzuki's slides and the book's example; none of the three papers was obtainable. | Suzuki & Ono papers [642], [823], [826] |
| `smn` | Built from the book alone, so there are no rate constants. | Ono, Fujiwara & Yuta, ECAL 2005 (LNAI 3630) |
| `typogenetics` | The Varetto code table rests on one source (Snare's thesis). | Varetto 1993, Morris 1989 |
| `laing-molecular-machines` | Laing's main result, self-reproduction by self-inspection, is **not** modelled: no readable source defines the synthesize/activate/convert instructions. Built from the book plus secondary descriptions. | Laing 1977 on U. Michigan Deep Blue, or his 1977 Binghamton dissertation (both open but blocked automated download; fetch by hand) |
| `urdar` | Invasion results reproduce; published diversity/efficiency *numbers* do not (trends do). The `on-gain` transform rule is inferred from those results. | The authors' Java program (math.chalmers.se/~torbjrn/Urdar) |
| `combinator-chemistry` | 7 of the published reactions (thesis ch. 7 tables, 2000 paper) release a different number of copies; the count depends on the unpublished reduction order. | Speroni di Fenizio's simulator source |

## 4. Scope

Agreed: frameworks and analyses are marked out of scope (not deleted);
wet chemistries become a small section of *given* topologies.
Undecided: whether `tierra` / `avida` / `corewar` (framework-like but
algorithmic) ever get wrappers.
