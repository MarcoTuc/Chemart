# Chemart library: implementation plan (all 98 chemistries)

> **Status: complete.** All 98 chemistries are implemented (4 `book`,
> 33 `book+decisions`, 61 `reconstructed`), the catalog validates with 0
> problems, and the full suite passes with no xfails. The plan below is kept
> as the record of how it was built; the text is as approved, so its
> forward-looking phrasing ("there is no generator code yet", the wave
> schedule) describes the starting point, not the current state. Open
> questions that outlived the build are in `to_decide.md`.

## Context

Phase 1 produced a catalog of 98 artificial chemistries (`catalog/chemistries/*.yaml`, loaded by
`chemart/catalog.py`). There is no generator code yet. Goal: a **mart**. Every chemistry is
reachable through one standardized interface that is easy for people and LLMs: discoverable,
JSON in/out, zero-argument defaults, actionable errors. Categorization is postponed (`to_decide.md`).
All 98 get implemented.

Decisions (2026-09-15):
- Concept-only / partial entries are **reconstructed from the original papers**.
- Frameworks (Tierra, Avida, Core War, aevol, CPM evo-devo, SRSim): **minimal Python re-implementations**.
- LLM access: Python API + JSON CLI + JSON-Schema tool definitions (no MCP server).
- Heavy scientific dependencies are **required**. Wheels are verified for Python 3.12 on Linux x86_64: ViennaRNA 2.7.2, RDKit 2026.3.6.
- **uv only, environment in `./.venv`.** Use `uv add` / `uv add --dev`, `uv sync`, `uv run`; commit `uv.lock`. No pip, no conda interpreter, no `UV_PROJECT_ENVIRONMENT`.

## Revision: structural problems in the previous draft and their fixes

| # | problem | fix |
|---|---|---|
| 1 | **Parameters had two sources of truth**: a Python `Params` dataclass per chemistry *plus* the catalog YAML, glued together by a sync test. | The **catalog YAML is the only parameter spec**. The validator and `describe()` read it; chemistry code receives an already-validated `p`. |
| 2 | **The catalog can't drive a JSON interface today**: of 428 params, 43 are `callable`, 95 are free-form `matrix`, 288 have no default, and ranges are prose. | **W0 normalizes the parameter vocabulary**: JSON types only (`int, float, bool, str, enum, list, dict`), a mandatory `default`, and machine-readable `min`/`max`/`choices`. Callables become named choices (`alpha_distribution: enum[uniform, lognormal]`) or rule sets written as text/JSON. Prose `range` is kept as a comment. |
| 3 | **Premature abstraction**: 4 generation modes, 4+ base classes and 11 engine modules designed before any chemistry existed (a "shared VM" for Tierra *and* automata reaction, whose instruction sets have nothing in common). | **One contract: a module-level function `generate(p, rng) -> Network`.** No base classes, no mode dispatch. Only two helpers are written up front, because 20+ chemistries need them: `expand()` (closure, book §12.6) and `soup()` (well-stirred event recording). Anything else is extracted into `chemart/helpers/` only when a second chemistry needs the same code. |
| 4 | **Universal budgets** (`max_species`, `steps`, …) **duplicated chemistry-specific size knobs** already in the catalog (`M`, `iterations`, `max_length`, `reduction_budget`). | No universal budgets. Each chemistry's own size parameters are its budget, with defaults small enough for a run of a few seconds. Paper-scale values are recorded in the parameter's `range` note. |
| 5 | **`seed` in two places** (a universal argument and 6 catalog params). | `seed` is a universal argument of `generate`, removed from catalog params. |
| 6 | **Provenance scattered over 4 places**: `docs/sources/<id>.md`, catalog `sources`, docstrings, `BOOK-ERRATA.md`. **`fidelity` overlapped `generator_ready`.** | Provenance lives only in the **catalog entry**: `fidelity` (replaces `generator_ready`), `sources` (citations/DOIs), `decisions` (each gap or book/paper erratum and its resolution). `BOOK-ERRATA.md` content moves into `matrix-chemistry.decisions` and the file is deleted. Docstrings just point at the catalog id. |
| 7 | **Oversized record**: a typed RateLaw class hierarchy, `kind: crn/non-crn`, 8 optional typed blocks, a hand-written JSON Schema file mirroring the dataclasses, a stored `provides` that could drift from content. | A slim record (below). Rates are plain dicts from a fixed vocabulary. Chemistry-specific payload (space, compartments, energies, analysis results) goes under one `extras` dict. `provides` is **computed** from content. No separate schema file; a round-trip test is the format contract. `kind` is dropped: an empty reaction list is legal and `summary()` says so. |
| 8 | **Surplus API surface**: `get()`, `react/expand/run` objects; the name `chemart.list` shadows a builtin; API names ≠ tool names. | Exactly three functions, named like the tools: `list_chemistries()`, `describe_chemistry(id)`, `generate_network(id, seed=None, **params)`. CLI and tool definitions are thin generated wrappers over them. |
| 9 | **The catalog rename** (`Chemistry` → `CatalogEntry`) was needed only because of the base class. | Not needed any more; `chemart/catalog.py` keeps its names. |
| 10 | **Plan documents overlapped**: `docs/PLAN.md`, README "Next", OUTPUT-SCHEMA §3, this plan. | This plan replaces `docs/PLAN.md`. README "Next" points to it. OUTPUT-SCHEMA keeps the rationale; its §3 is replaced by a short description of the implemented record. |
| 11 | **Research could stall** on paywalled or lost papers, with no rule for what happens then. | Fallback: implement from the book, set `fidelity: book+decisions`, and record the missing source in `decisions`. Never block a wave on one paper. |
| 12 | **Dependencies were front-loaded** (networkx listed without a user). | W0 adds only `numpy`, `scipy`, `pyyaml`, dev `pytest`. `rdkit` / `ViennaRNA` are added with `uv add` in the wave that first uses them (still required deps, just added when needed). |

## Interface

```python
import chemart
chemart.list_chemistries()                 # [{id, name, summary, fidelity}] x 98
chemart.describe_chemistry("matrix-chemistry")
    # {id, name, summary, book, sources, decisions, params: JSON Schema from the catalog}
net = chemart.generate_network("matrix-chemistry", seed=0, N=4)   # no params needed
net.summary()          # species/reaction counts, status, provides
net.to_text()          # "s1 + s6 -> s1 + s6 + s4" per line, rate appended if present
net.to_dict()          # JSON-ready; Network.from_dict(d) inverts it
species_ids, R, P = net.matrices()          # scipy.sparse; S = P - R
```
CLI: `uv run chemart list | describe <id> | generate <id> -p N=4 --seed 0 --format json|text|summary`.
`chemart.tool_definitions()` returns the three functions as JSON-Schema tool specs.

Guarantees: zero-argument runs finish in seconds; invalid params raise an error naming the parameter, its allowed values and its meaning (from the catalog); the same seed gives the same network; all values are JSON (species `structure` is a string).

## Network record (`chemart/network.py`, one dataclass file)

```python
Species(id: str, structure: str | None = None)
Reaction(reactants: dict[str,int], products: dict[str,int],
         rate: dict | None = None,   # {"law": "mass-action"|"michaelis-menten"|"hill"|"saturating"|"arrhenius", ...params}
         count: int | None = None)   # firings, when the network was observed from a simulation
Network(chemistry: str, params: dict, seed: int | None,
        status: "complete" | "truncated" | "observed",
        species: list[Species], reactions: list[Reaction],
        initial_state: dict | None = None, inflow: dict | None = None, outflow: dict | None = None,
        extras: dict = {})
# computed: provides, matrices(), summary(), to_text(), to_dict()/from_dict()
```
`extras` has reserved keys: `space`, `compartments`, `energies`, `conservation`, `analysis`, `interaction_law`.
The computed `provides` reads the reserved keys (plus rates, catalysts via reactants ∩ products, inflow/outflow, initial_state). Any other key is free-form.
Wolkenhauer k→c conversion (book eq. A.4) is a function in `chemart/kinetics.py`, the only other core module.

## Code layout

```
chemart/
  __init__.py        re-exports the three functions + tool_definitions
  api.py             list/describe/generate: param validation from catalog, registry by naming convention
  catalog.py         (exists) + param-vocabulary validation, fidelity/sources/decisions fields, status command
  network.py         record above
  kinetics.py        rate-law vocabulary check, k->c
  expand.py          closure: expand(react, seed_species, arity, max_species) -> (species, reactions, status)
  soup.py            soup(react, initial, steps, rng, dilution) -> observed reactions with counts
  cli.py             argparse over api
  chemistries/<id_with_underscores>.py   each defines generate(p, rng) -> Network
  helpers/           created only on second use
tests/
  test_catalog.py    (exists)
  test_contract.py   parametrized over all catalog ids
  chemistries/test_<id>.py   published-number regression tests
```
Registry = convention: id `matrix-chemistry` → module `chemart.chemistries.matrix_chemistry`, imported lazily. No decorators.

## Work plan

### W0: foundation (no chemistries yet)
1. `uv add numpy scipy`, `uv add --dev pytest`, drop the `dev` optional-dependency table, add `[project.scripts] chemart`. Delete `main.py`.
2. Catalog schema v2 in `catalog/SCHEMA.md` + `catalog.py` validator: the JSON parameter vocabulary (fix 2), `seed` removed (fix 5), `fidelity/sources/decisions` replacing `generator_ready` (fix 6). The YAML is migrated **per entry in its wave**, not all 428 params up front. No version flag: the validator enforces v2 exactly for entries that have a `chemart/chemistries/<id>.py` module, so migration and implementation cannot drift apart.
3. `network.py`, `kinetics.py`, `expand.py`, `soup.py`, `api.py`, `cli.py`, `tool_definitions()`.
4. `test_contract.py` + the coverage report (`python -m chemart.catalog status`: implemented / v2-migrated / fidelity per wave).
5. Replace `docs/PLAN.md` with this plan; move BOOK-ERRATA content into the catalog; update README.

### Per-chemistry procedure (identical in every wave)
1. Read the book section (`pdftotext -layout` into the job tmp dir) and, for `partial`/`no` entries or when the book is thin, the original papers from `refs` (paper-lookup skill). Fallback per fix 11.
2. Migrate the catalog entry to v2: JSON params with small defaults, `fidelity`, `sources`, `decisions`. Fix wrong facts.
3. Write `chemart/chemistries/<id>.py`: `generate(p, rng) -> Network`.
4. Add a regression test reproducing published numbers (book or paper) where any exist; otherwise the catalog's `phenomena` as property tests.

### Waves (ordered by increasing research and mechanism complexity)

| wave | chemistries |
|---|---|
| **W1 written-down networks (25)** | dimerization, chameleon, brusselator, oregonator, repressilator, michaelis-menten, hill-kinetics, logistic-chemistry, lotka-volterra, replication-death, selection-equation, replicator-equation, quasispecies, okamoto-switch, naming-game-ac, organization-computing, analog-function-crn, disperser, chemoton, gard, mechanical-self-assembly, nuclear-reaction-networks, n-economy, flow-ac, french-flag |
| **W2 random topologies (12)** | random-catalytic-networks, bigan-conservative-crn, jain-krishna, kauffman-autocatalytic-sets, bagley-farmer, rbn, nk-landscape, metabolic-robot-controller, arn, acgp, farmer-immune, ecolab |
| **W3 string & number constructive (17)**: first use of `expand`/`soup` | matrix-chemistry, prime-number-chemistry, automata-reaction, smn, bnc-cell, mcs-bl, sac, mccaskill-polymer-tm, ikegami-hashimoto, typogenetics, stringmol, fraglets, molecular-tsp, urdar, proof-ac, laing-molecular-machines, tominaga-stacked-strings |
| **W4 rewriting systems (14)**: rule-set params as text/JSON, defaulting to the published worked example | alchemy, combinator-chemistry, reflexive-ac, gamma, arms, cham, soas, ccm, p-systems, brane-calculi, mgs, l-systems, kappa-calculus, high-order-chem |
| **W5 spatial, CA, agents (14)** | squirm3, swarm-chemistry, autopoiesis-vmu, ono-ikegami-protocell, dorin-korb-ecosystem, sr-loops, bondable-ca, ca-embedded-particles, isologous-diversification, cpm-grn-evodevo, srsim, nac, evolve-series, social-communication-ac |
| **W6 digital organisms (5)**: network = observed replication/parasitism events | coreworld, corewar, tierra, avida, aevol |
| **W7 physico-chemical, analyses, concept-only, wet (11)**: `uv add rdkit ViennaRNA` here | toychem, synthon, rna-folding-ac, hbcb-psd, energy-gated-collision, music-ac, conrad-enzymatic, raf, dna-automaton, dna-hpp, self-propelled-droplets |

Non-generator entries follow the same contract. Analyses (`raf`, `nk-landscape`, `energy-gated-collision`, `ca-embedded-particles`) generate their model instance and put the analysis result in `extras`. Wet entries return the in-silico network of the published experiment. Entries with no reactions (e.g. swarm chemistry) return species plus `extras` describing the interaction law.

Tests already verified for W3 (this session): matrix chemistry table 3.3; 4 self-replications / 76 replications at N=4; at N=9, foldings 1–4 give 14/12028, 122/21310, 18/11822, 94/16830; table 3.5 row 165; closure of s1..s15 at N=9 = s1..s27 minus s20..s23.

At the end of each wave: full suite green, then commit (agreed: one commit per wave, plus infrastructure commits).

### Parallel implementation (from W2 on)

Remaining chemistries are implemented by parallel agents, one per chemistry,
8 per batch, in the shared working tree:
- The catalog is split into one YAML file per chemistry, so each agent owns
  exactly `chemart/chemistries/<mod>.py`, `catalog/chemistries/<id>.yaml` and
  `tests/chemistries/test_<mod>.py`. Shared files are integrator-only; agents
  request changes in their report.
- Every agent follows [`IMPLEMENTING.md`](IMPLEMENTING.md) and must pass its
  definition of done: `python -m chemart.catalog validate --only <id>`,
  `CHEMART_ONLY=<id> pytest tests/test_contract.py tests/chemistries/test_<mod>.py`,
  and `chemart generate <id>`.
- After each batch the integrator reviews the reports (decisions, sources),
  applies requested shared changes, runs the full suite, regenerates
  `docs/CATALOG.md`, and commits.

## Verification

- `uv sync` builds `./.venv`; `uv run which python` resolves inside it; after W7, `uv run python -c "import RNA; from rdkit.Chem import rdEHTTools"` works.
- `uv run pytest`. `test_contract.py`, for every implemented id: no-arg `generate_network` succeeds in < 5 s; `from_dict(to_dict())` is lossless; the same seed gives an equal network; unknown or out-of-range params raise `ValueError` naming the parameter; every reaction's species ids exist; `matrices()` shapes match; `rate` laws are in the vocabulary; the catalog `provides` claims are ⊇ the computed `provides`.
- A final gate test: implemented ids == catalog ids (98) and every entry is on schema v2.
- `uv run python -m chemart.catalog validate` → 0 problems.
- CLI: `uv run chemart generate matrix-chemistry -p N=4 --format json | uv run python -m json.tool`; `uv run chemart describe gamma`.
- `tool_definitions()`: each schema validates, and calling through it gives the same result as the direct function.
