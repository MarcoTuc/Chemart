# Fidelity and trust

Chemart reproduces published models. Some are specified completely in the
book; some required reconstruction from papers that were hard to get; a few
rest on abstracts. Every entry records which, and that record is the most
important metadata in the catalog when you are about to state a result.

## Contents
- [The three fidelity levels](#the-three-fidelity-levels)
- [decisions: read this before quoting a number](#decisions-read-this-before-quoting-a-number)
- [Entries that are known to be thin](#entries-that-are-known-to-be-thin)
- [Entries that implement only part of their scope](#entries-that-implement-only-part-of-their-scope)
- [How to phrase a result honestly](#how-to-phrase-a-result-honestly)

## The three fidelity levels

```python
chemart.describe_chemistry("brusselator")["fidelity"]
```

| level | count | meaning |
|---|---|---|
| `book` | 4 | implemented exactly as Banzhaf & Yamamoto specify |
| `book+decisions` | 33 | the book left gaps; each filled choice is listed in `decisions` |
| `reconstructed` | 61 | built from the original papers, cited in `sources` |

`reconstructed` is not a warning — most of these are the *better* entries,
because the primary paper is more precise than a survey chapter. Several were
cross-checked against compiled upstream simulators: `avida` against Avida
2.14.0 on 626 dividing genomes with zero mismatches, `corewar` against pMARS
over ~1250 battles, `coreworld` against pMARS instruction by instruction.

What matters is not the label but the `decisions` list underneath it.

## decisions: read this before quoting a number

```python
d = chemart.describe_chemistry("hbcb-psd")
d["fidelity"]     # 'book+decisions'
d["sources"]      # what was actually read
d["decisions"]    # every gap and how it was resolved
```

`decisions` records, one line each: which sources could not be obtained, which
rules were inferred, which published numbers are reproduced and which are not,
and any errata found in the book or the papers. Entries routinely say things
like "the primary paper is paywalled; the hierarchy size and bond energies are
parameters, not published values".

This is where the honesty lives. An entry whose default parameters are the
implementer's choice will say so; if you report that number as the paper's, you
are reporting something the catalog explicitly told you it wasn't.

## Entries that are known to be thin

The repository's `to_decide.md` keeps the current list. At the time of writing
it names roughly a dozen whose key semantics rest on inference because the
primary source was unavailable, including:

- **`evolve-series`** — the weakest entry: none of its six publications could
  be obtained, so the function table, matter/energy economy and population
  rules are reconstruction from one book paragraph plus abstracts.
- **`hbcb-psd`** — primary paper unreachable behind a publisher challenge; the
  hierarchy size, bond energies and every numeric result are parameters.
- **`ono-ikegami-protocell`** — the book's headline phenomena (membrane
  closure, osmotic permeability, division) are **not** reproduced.
- **`nac`** — only NAC's passive layer (graph rewiring, polarity demixing) is
  implemented, from archived slides; the active layer is absent rather than
  invented.
- **`sac`** — the string-rewriting semantics are inferred from the author's
  slides and the book's worked example; none of the three papers was
  obtainable.
- **`reflexive-ac`**, **`mccaskill-polymer-tm`**, **`laing-molecular-machines`**,
  **`typogenetics`**, **`urdar`**, **`combinator-chemistry`**,
  **`rna-folding-ac`** — each with a specific documented gap.

Read `to_decide.md` for the current table with what would settle each one.

Where a published result could not be reproduced, entries say so instead of
tuning parameters until it appeared — `evolve-series` reports that enzyme
efficiency does *not* improve from near-optimal founders;
`ca-embedded-particles` reports that one published lookup table does not
reproduce its published behaviour under any convention tried.

## Entries that implement only part of their scope

Distinct from low fidelity: these are sound for what they cover, but one
documented half is missing because its source is unobtainable. `music-ac`
implements the homophonic system but not the polyphonic one; `nac` the passive
layer but not the active; `rna-folding-ac` the genotype-phenotype map but not
the evolutionary results; `isologous-diversification` stage 1 but not stages
2-4; `laing-molecular-machines` the machine chemistry but not self-inspection.
`to_decide.md` has the table.

None of these is a defect to fix by guessing. They mark where obtaining one PDF
would materially extend the library.

## How to phrase a result honestly

A small amount of care costs nothing and keeps the library trustworthy:

- **From a `book` entry**: "the Brusselator as given in Banzhaf & Yamamoto §…"
- **From `reconstructed` with matching tests**: "reproduces the published
  gestation counts from Ray (1991) appendix C" — name what was reproduced.
- **From a thin entry**: "Chemart's reconstruction; the primary paper was not
  obtainable, and the entry records these values as implementation choices."

If a user asks "what does model X predict?", check whether the specific number
they want is one the entry reproduces or one it chose. `decisions` and
`phenomena` together answer that in under a minute, and the answer is
frequently "the published claim is recorded but not asserted by the tests",
which is itself the useful reply.
