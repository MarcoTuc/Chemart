# Fidelity and trust

Chemart reproduces published models. Some are specified completely in the book;
some required reconstruction from papers that were hard to obtain; a few rest on
abstracts. Every entry records which — and that record is the most important
metadata in the catalog when you are about to state a result.

## The three fidelity levels

```python
chemart.describe_chemistry("brusselator")["fidelity"]
```

| level | count | meaning |
|---|---|---|
| `book` | 4 | implemented exactly as Banzhaf & Yamamoto specify |
| `book+decisions` | 33 | the book left gaps; each filled choice is listed in `decisions` |
| `reconstructed` | 61 | built from the original papers, cited in `sources` |

`reconstructed` is **not** a warning. It is usually the stronger label, because
a primary paper is more precise than a survey chapter. Several of these were
cross-checked against compiled upstream simulators:

- **`avida`** against Avida 2.14.0 in analyze mode — 626 dividing genomes, zero
  mismatches on gestation, merit, executed length and tasks;
- **`corewar`** against pMARS 0.9.2 over roughly 1250 battles, agreeing on the
  whole core;
- **`coreworld`** against pMARS 0.9.4, instruction by instruction;
- **`tierra`** reproduces the ancestor's published 839/813 gestation counts and
  the 6.02 genebank's 827/809.

What matters is not the label but the `decisions` list underneath it.

## `decisions`: read this before quoting a number

```python
d = chemart.describe_chemistry("hbcb-psd")
d["fidelity"]     # 'book+decisions'
d["sources"]      # what was actually read
d["decisions"]    # every gap and how it was resolved
```

`decisions` records, one line each: which sources could not be obtained, which
rules were inferred, which published numbers are reproduced and which are not,
and any errata found in the book or the papers.

This is where the honesty lives. An entry whose default parameters are the
implementer's choice says so. If you report such a number as the paper's, you
are reporting something the catalog explicitly told you it wasn't.

## Errata found along the way

Verification against the original sources caught errors in every layer, and each
is recorded in the relevant entry's `decisions`:

- **In the catalog** — wrong authors for `bondable-ca` (Hatcher, Banzhaf & Yu,
  not "Gross") and `music-ac` (Miura, not "Miyamoto"); the wrong book section
  for `self-propelled-droplets`; `constructive` wrong for `hbcb-psd`.
- **In the book** — a misfiled reference; a figure caption listing C4H10 while
  plotting C6H10; `dna-hpp`'s edge encoding printed with the strands
  complemented and the halves swapped; an algorithm in §2.2.6 contradicting its
  own equation 2.32 by a factor of 3.4.
- **In a published paper** — Hutton (2002) states `5n² + 21n + 4` rules where
  expanding R1–R8 gives `+24`, which is the 188 the same paragraph reports.

## Entries that are thin

`to_decide.md` in the repository keeps the current list — roughly a dozen whose
key semantics rest on inference because the primary source was unavailable:

- **`evolve-series`** — the weakest entry. None of its six publications could be
  obtained, so the function table, matter/energy economy and population rules
  are reconstruction from one book paragraph plus abstracts.
- **`hbcb-psd`** — primary paper unreachable behind a publisher challenge; the
  hierarchy size, bond energies and every numeric result are parameters.
- **`ono-ikegami-protocell`** — the book's headline phenomena (membrane closure,
  osmotic permeability, division) are **not** reproduced.
- **`nac`** — only the passive layer is implemented; the active layer is absent
  rather than invented.
- **`sac`** — string-rewriting semantics inferred from the author's slides and
  the book's worked example; none of the three papers was obtainable.
- **`reflexive-ac`**, **`mccaskill-polymer-tm`**, **`laing-molecular-machines`**,
  **`typogenetics`**, **`urdar`**, **`combinator-chemistry`**,
  **`rna-folding-ac`** — each with a specific documented gap.

## Negative results are reported, not tuned away

Where a published result could not be reproduced, the entry says so:

- **`evolve-series`** — enzyme efficiency does *not* improve from near-optimal
  founders. Selection on the match is real (a scrambled genome leaves zero
  descendants), but drift wins. Parameters were not adjusted until the expected
  result appeared.
- **`ca-embedded-particles`** — one published lookup table does not reproduce
  its published behaviour under any convention tried, including all 480
  single-digit corrections.
- **`self-propelled-droplets`** — no droplet speed is asserted, because no
  accessible source gives one for this chemistry. The widely quoted 6 mm/s
  belongs to a different system.

## Entries implementing only part of their scope

Distinct from low fidelity: sound for what they cover, but one documented half is
missing because its source is unobtainable. `music-ac` implements the homophonic
system but not the polyphonic; `nac` the passive layer but not the active;
`rna-folding-ac` the genotype–phenotype map but not the evolutionary results;
`isologous-diversification` stage 1 but not stages 2–4;
`laing-molecular-machines` the machine chemistry but not self-inspection.

None of these is a defect to fix by guessing. They mark where obtaining one PDF
would materially extend the library.

## How to phrase a result honestly

- **From a `book` entry** — "the Brusselator as given in Banzhaf & Yamamoto".
- **From `reconstructed` with matching tests** — name what was reproduced:
  "reproduces the published gestation counts from Ray (1991), appendix C".
- **From a thin entry** — "Chemart's reconstruction; the primary paper was not
  obtainable, and the entry records these values as implementation choices."

If someone asks "what does model X predict?", check whether the specific number
they want is one the entry reproduces or one it chose. `decisions` and
`phenomena` together answer that in under a minute, and the answer is often "the
published claim is recorded but not asserted by the tests" — which is itself the
useful reply.
