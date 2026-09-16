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
| `nac` | Only NAC's **passive layer** (graph rewiring, polarity demixing) is implemented, from Suzuki's archived JSAI 2004 slides. The active layer — node programs, polymerase/helicase/splitase/replicase, centrosome-driven division, chain folding — is absent, not invented, because every paper specifying it is closed. Same author and same problem as `sac`. | Suzuki [824], [825], [827]; the ALife IX 2004 NAC paper |
| `typogenetics` | The Varetto code table rests on one source (Snare's thesis). | Varetto 1993, Morris 1989 |
| `laing-molecular-machines` | Laing's main result, self-reproduction by self-inspection, is **not** modelled: no readable source defines the synthesize/activate/convert instructions. Built from the book plus secondary descriptions. | Laing 1977 on U. Michigan Deep Blue, or his 1977 Binghamton dissertation (both open but blocked automated download; fetch by hand) |
| `urdar` | Invasion results reproduce; published diversity/efficiency *numbers* do not (trends do). The `on-gain` transform rule is inferred from those results. | The authors' Java program (math.chalmers.se/~torbjrn/Urdar) |
| `evolve-series` | **The weakest entry in the mart.** None of the six EVOLVE publications was obtainable (all closed, no repository copy, archive.org down all session), so the function table, matter/energy economy, matching semantics and population rules are all reconstruction from one book paragraph plus abstracts. Whether Conrad's actual machinery resembles it at all is unknown. | Conrad & Pattee 1970; Conrad & Strizich 1985; Rizki & Conrad 1985; Brewster & Conrad 1999 — any full text at all |
| `ono-ikegami-protocell` | The book's headline phenomena are **not** reproduced: membrane filaments never close into protocells, no closed membrane retains autocatalyst (`enclosed_A = 0` in every run), and no growth/division occurs. Only the accessible 1D predecessor's published constants and the book's reaction set are tested. The anisotropy field `F[k,o]` — likely the thing that makes closure work — is the implementer's third attempt and is a reconstruction. | Ono & Ikegami [639], [641] (paywalled); [538] for the 3D extension |
| `reflexive-ac` | The product construction (reflexive composition of two finite-state machines) is read off Salzberg & Sayama's 2025 restatement, not off the 2007 paper, which is paywalled; that paper's own elastic rules, reactor and results are unverified. Two figures of the 2025 paper (12, 13) are not reproduced. | Salzberg 2007, BioSystems 87:1-12, and his two 2006 papers |
| `combinator-chemistry` | 7 of the published reactions (thesis ch. 7 tables, 2000 paper) release a different number of copies; the count depends on the unpublished reduction order. | Speroni di Fenizio's simulator source |

## 4. Record design: per-reaction metadata

Several generators need facts attached to single reactions that are not a
rate law: L-system successor probabilities, which rule produced a reaction
(kappa-calculus, CHAM, MGS, Gamma), threshold gates (chemoton), inhibitors
(metabolic-robot-controller). Today they live in aligned lists under
`extras` (e.g. `extras.probabilities[i]` for reaction i) or as extra scalar
keys on the rate dict. Proposal from the l-systems implementer: an optional
`Reaction.extra: dict | None`. It changes the record for every chemistry, so
decide once all waves are in: adopt and migrate the aligned lists, or keep
the current convention and document it.

## 5. Scope

Note: the repository has no `LICENSE` file and `pyproject.toml` declares no
license, so nothing is in conflict today — but see §6 before publishing.

Agreed: frameworks and analyses are marked out of scope (not deleted);
wet chemistries become a small section of *given* topologies.
Settled by implementation: `tierra`, `avida`, `corewar` and `coreworld` are
not wrappers but minimal faithful re-implementations, and their `kind` moved
from `framework` to `generator`. `aevol` and `high-order-chem` are still
`framework`.

## 6. Licensing

`tests/chemistries/test_corewar.py` embeds two Redcode warriors verbatim from
the pMARS distribution — Validate 1.1R and Rave, both by Stefan Strack, GPL-2
— with attribution in the file. They are what make the pMARS cross-checks
readable (Validate is *the* MARS conformance program).

The repository currently declares no license at all, so there is no conflict
yet. Before Chemart is published, decide one of:

1. license Chemart GPL-2-or-later (simplest, but it is the strongest copyleft
   of anything vendored so far);
2. keep the warriors in a separate `tests/fixtures/pmars/` directory with its
   own GPL-2 notice, and license Chemart itself permissively;
3. drop the two warriors and keep only the recorded pMARS core hashes, losing
   readability in two tests.

Also worth a pass at the same time: every chemistry ported from upstream
source (`stringmol`, `tierra`, `avida`, `corewar`, `coreworld`,
`high-order-chem`) should say in `sources` which upstream license its port
derives from.

## 7. Conventions the genome-carrying chemistries disagree on

Raised by the `aevol` implementer; they affect `tierra`, `avida`, `aevol`,
`stringmol` and `squirm3` alike, so decide once rather than per chemistry.

**7a. How much structure to store.** Every genotype currently carries its
full genome as `Species.structure`, so `aevol`'s default network is ~1 MB of
JSON and a long run is several MB. Options: keep it (self-contained records,
but large), store genomes only for the ancestor and the final population, or
store a hash plus a lookup table in `extras`. Note the record is meant to be
JSON that an LLM can read, which argues against multi-MB defaults.

**7b. Which non-replication events are reactions.** `aevol` puts gene
expression (`G -> G + P1 + … + Pk`) in the network alongside the replication
events, so its protein species actually react. `avida` and `tierra` keep
analogous non-replication events (task completions, instruction execution) in
`extras` instead. Both readings are defensible; they should not coexist
unexamined, since `provides` and any cross-chemistry comparison depend on it.

## 8. A catalysed Michaelis-Menten rate law?

Raised by the `isologous-diversification` implementer. Kaneko & Yomo's
eq. (1) has the enzyme-proportional saturating term

    e1 · x(j) · x(m) / (1 + x(m)/x_M)

which is not in `chemart.kinetics.RATE_LAWS`. It was emitted as mass-action
`k=e1` (its correct dilute limit) with `saturation`, `saturated_species` and
`x_M` as extra scalar keys — the documented fallback, but a degradation: the
published law *is* known here, it simply has no vocabulary entry.

Proposal: add `catalysed-michaelis-menten` with params `(k, K)` plus an
enzyme species key, propensity `k·[E]·[S]/(1 + [S]/K)`. This shape recurs
across the enzyme-kinetics literature, so it likely serves more than one
entry.

Two things to settle before adopting, both cross-cutting:
1. **Survey first.** Several chemistries may have degraded similar laws to
   extra keys. Adding the law without migrating them leaves the catalog
   inconsistent in a way `provides: rate-law` would misreport.
2. **`tests/chemistries/odes.py` needs the same change.** It unpacks
   michaelis-menten as `(i, _), = reac`, so it cannot integrate *any* MM
   reaction that has a second (catalyst) reactant. Latent today because
   every chemistry degrades instead of emitting one; it becomes a real bug
   the moment this law lands.

Deliberately not adopted mid-wave: the rate-law vocabulary is part of the
record format, and changing it while 90+ entries are already written is a
migration, not an addition.
