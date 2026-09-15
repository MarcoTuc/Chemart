# Chemart catalog schema

One YAML file per family; each file is a mapping `chemistries: [ ... ]`.
Every entry uses the fields below. Fields marked (req) must be present.

The catalog is the **only specification of each generator's parameters**:
`chemart.describe_chemistry` builds its JSON Schema from it, and
`chemart.generate_network` validates arguments against it.

**Schema versions.** An entry is on **v2** exactly when its generator module
`chemart/chemistries/<id with underscores>.py` exists; the validator then
enforces the v2 rules marked below. Entries without a module are still v1
(prose-level parameters, `generator_ready`). Entries migrate to v2 one by one
as they are implemented.

## Identity
| field | meaning |
|---|---|
| `id` (req) | stable kebab-case key, unique across the whole catalog; module name is the id with `-` → `_` |
| `name` (req) | canonical name |
| `aliases` | other names used in the literature |
| `origin` | `author, year` of the original proposal |
| `book` (req) | section(s) of Banzhaf & Yamamoto (2015) describing it |
| `refs` | bracketed citation numbers as they appear in the book's bibliography, plus DOIs/URLs where known |

## Classification
| field | meaning |
|---|---|
| `family` (req) | `core`, `rewriting`, `automata`, `bio-inspired`, `origin-of-life`, `evolutionary-dynamics`, `network`, `spatial`, `application`, `systems-biology`, `wet`, `non-chemical` |
| `kind` (req) | what the thing *is*: `generator`, `formalism`, `framework`, `analysis`, `wet` |
| `constructive` (req) | `true` if the species set is open/unbounded and grows at run time; `false` if S is fixed and enumerable up front |

## The (S, R, A) triple — Banzhaf & Yamamoto §2.3
| field | meaning |
|---|---|
| `S` (req) | `repr`: molecule representation; `definition`: `explicit`\|`implicit`; `notes` |
| `R` (req) | `definition`: `explicit`\|`implicit`; `arity`: reaction order(s); `scheme`: canonical reaction template; `notes` |
| `A` (req) | `reactor`: `well-stirred-multiset`, `ode`, `ssa`, `lattice-2d`, `continuous-space`, `compartments`, `graph-rewrite`, `maximally-parallel`, `sequential-vm`, ...; `dilution`: how population size is bounded; `notes` |

## Parameters
`params:` list of

```yaml
- name: N
  type: int          # v2: int | float | bool | str | enum | list | dict
  default: 4         # v2: required, and must satisfy the spec below
  min: 1             # optional, int/float only
  max: 100           # optional, int/float only
  choices: [1, 2]    # required for enum
  range: "perfect square; the book also uses 9, 16, 25"   # optional prose
  role: structural   # structural | kinetic | thermodynamic | population | spatial | stochastic | selection
  meaning: "string length; also fixes |S| = 2^N - 1"      # v2: required
```

v2 rules:
- **JSON values only.** Functions become named choices
  (`alpha_distribution: enum [uniform, lognormal]`); rule sets become text
  or JSON (`rules: str` in a documented mini-language, or `list`/`dict`).
- **Small defaults.** `generate_network(id)` with no arguments must finish
  in seconds. Paper-scale values go in `range`.
- **No `seed` parameter.** The seed is an argument of `generate_network`.
- v1 types `matrix`, `callable`, `seed` are not allowed.

`role` is what a knob *does*, and is what lets Chemart offer coherent
scaling across chemistries (e.g. "scale every structural knob down until
|S| < 500").

## Output capability tier — what a generator can actually hand you
`provides:` a list drawn from

- `topology` — species list + who-reacts-with-whom (always present)
- `stoichiometry` — separate reactant/product multiplicities
- `catalysts` — reactions where a species appears on both sides and must be kept
- `rate-constants` — the chemistry itself prescribes k's (not just "pick some")
- `rate-law` — a non-mass-action propensity is part of the definition (MM, Hill, saturating, Arrhenius, swarm force law)
- `energies` — per-species free energies / bond energies
- `thermodynamic-consistency` — reverse rates constrained by ΔG, detailed balance holds
- `mass-conservation` — an atom/mass vector m with Sᵀm = 0
- `flow` — inflow/outflow (food set, dilution flux, CSTR) is part of the definition
- `space` — positions, lattice, or diffusion coefficients
- `compartments` — nested membranes/cells
- `initial-state` — the chemistry prescribes a seed multiset
- `sequence-structure-function` — molecules carry an internal structure that the reaction rule reads

For implemented entries, the tags a generated `Network` computes from its
content must be a subset of `provides` (checked by `tests/test_contract.py`).

## Provenance and status
| field | meaning |
|---|---|
| `generator_ready` | v1 only: `yes` (fully specified in the book), `partial` (needs a design decision), `no` (concept only / wet only) |
| `fidelity` | v2, required: `book` (implemented exactly as specified in the book), `book+decisions` (gaps filled; each listed in `decisions`), `reconstructed` (built from the original papers listed in `sources`) |
| `sources` | citations/DOIs actually used for the implementation (required for `reconstructed`) |
| `decisions` | every gap, ambiguity or erratum in the book or papers, and how Chemart resolved it (required for `book+decisions`) |
| `reference_impl` | existing code we can port or check against (e.g. PyCellChemistry module) |
| `phenomena` | what the chemistry is known to produce — used for regression tests |
| `notes` | anything else |
