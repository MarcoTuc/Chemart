# Chemart — a mart of artificial chemistries

Chemart makes every artificial chemistry in the literature available behind one
interface, so you can generate its reaction network, look at it, simulate it,
and compare it with any other — without reimplementing 98 papers.

All 98 chemistries from the field's reference survey — Wolfgang Banzhaf & Lidia
Yamamoto, *Artificial Chemistries* (MIT Press, 2015) — are catalogued **and
implemented**, each with tests reproducing published results wherever any exist.

```python
import chemart

chemart.list_chemistries()                      # all 98
net = chemart.generate_network("brusselator", seed=1)
print(net.to_text())
# A -> X  [mass-action k=1.0]
# B + X -> Y + D  [mass-action k=1.0]
# 2 X + Y -> 3 X  [mass-action k=1.0]
# X -> E  [mass-action k=1.0]
```

📖 **Full documentation:** a browsable site with a page for every chemistry.
Run `uv run mkdocs serve` and open <http://127.0.0.1:8000> — see
[Documentation](#documentation) for publishing it.

---

## What is an artificial chemistry?

An artificial chemistry (AC) is a model of a chemistry that need not be real.
Banzhaf & Yamamoto define one as a triple **(S, R, A)**:

| | |
|---|---|
| **S** | the set of possible **molecules** — numbers, bitstrings, λ-terms, graphs, machine programs, RNA sequences |
| **R** | the **reaction rules** that say what happens when molecules meet |
| **A** | the **algorithm** — the reactor: a well-stirred multiset, an ODE system, a lattice, a virtual machine |

That definition is broad on purpose, and it is why the field is hard to survey:
the Brusselator (three coupled differential equations) and Tierra (self-replicating
assembly programs competing for CPU time) are both artificial chemistries, and
they share almost nothing operationally.

Chemart's bet is that they can still share an **output**. Whatever a chemistry is
internally, it can hand back a chemical reaction network: which species exist,
which reactions consume and produce them, and — when the model defines them —
the rates. That common output is what makes 98 incomparable models comparable.

## How Chemart works

### The catalog is the specification

Each chemistry is one YAML file in [`catalog/chemistries/`](catalog/chemistries/).
It records the (S, R, A) triple, the parameters with their types, defaults,
meanings and published ranges, what the generator can supply, the book section
and papers behind it, and — importantly — every decision taken where the sources
were incomplete.

The catalog is the **only** parameter specification. There is no second copy in
Python that could drift from it: `describe_chemistry` builds a JSON Schema from
the YAML, and `generate_network` validates your arguments against the same YAML.

### One function per chemistry

Each chemistry is a module defining exactly one function:

```python
def generate(p, rng) -> Network
```

`p` holds the validated parameters; `rng` is a seeded NumPy generator. No base
classes, no registration, no plugin system — the id `matrix-chemistry` resolves
to `chemart.chemistries.matrix_chemistry` by naming convention, imported lazily.
Adding a chemistry means adding three files and touching nothing else.

### Everything returns the same record

```python
Network(
    species,        # [Species(id, structure)] - structure carries the bitstring,
                    #   λ-term, genome or fold when the molecule has one
    reactions,      # [Reaction(reactants, products, rate, count)]
    status,         # "complete" | "truncated" | "observed"
    initial_state, inflow, outflow,
    extras,         # space, compartments, energies, conservation, analysis, ...
    chemistry, params, seed,     # provenance: enough to reproduce it exactly
)
```

It is plain JSON throughout, so `Network.from_dict(net.to_dict())` round-trips
exactly — a property tested for all 98. From there:

```python
net.summary()                    # counts, status, what it provides
net.to_text()                    # one reaction per line
ids, R, P = net.matrices()       # sparse; net stoichiometry S = P - R
net.provides                     # capability tags computed from content
```

### `status` is the field people skip

It changes what the network *means*, and comparing networks across different
statuses is how you get nonsense:

- **`complete`** — the whole defined network, or a closure that finished.
- **`truncated`** — a closure stopped by a size budget. There is more chemistry
  you are not seeing; raise the budget.
- **`observed`** — the reactions that *fired in one simulated run*, each with a
  firing `count`. A reaction's absence tells you about that run, not about the
  chemistry.

Half of these chemistries are **constructive**: the species set is open and
grows as reactions produce new molecules. For those, the interesting object is
a closure, and `truncated` is a normal answer rather than a failure.

### Provenance is part of the data

Reproducing published models honestly means admitting where the sources ran out.
Every entry carries a `fidelity`:

| | count | meaning |
|---|---|---|
| `book` | 4 | implemented exactly as the book specifies |
| `book+decisions` | 33 | the book left gaps; each filled choice is in `decisions` |
| `reconstructed` | 61 | built from the original papers, cited in `sources` |

`reconstructed` is usually the *stronger* label — a primary paper is more precise
than a survey chapter, and several of these were cross-checked against compiled
upstream simulators (Avida against Avida 2.14.0 on 626 genomes with zero
mismatches; Core War against pMARS over ~1250 battles).

What matters is the `decisions` list underneath. Where a paper was unobtainable,
the entry says so and records its parameters as choices rather than inventing
published values. [`to_decide.md`](to_decide.md) collects the thin entries in one
place, along with cases where a published result could **not** be reproduced —
those are reported, not tuned away.

## What's in the box

| | |
|---|---|
| chemistries catalogued and implemented | 98 |
| constructive (open, growing species set) | 51 |
| carry their own rate constants or rate law | 44 |
| carry energetics / thermodynamic consistency | 11 |
| declare a conservation law | 31 |
| define space | 21 |
| define compartments | 10 |

By kind: 80 generators, 7 formalisms, 4 analyses, 4 wet, 3 frameworks.
Browse them all in [`docs/CATALOG.md`](docs/CATALOG.md) or on the docs site.

## Quickstart

The environment is managed with [uv](https://docs.astral.sh/uv/); run everything
through `uv run` so you get the project's own `./.venv`.

```bash
uv sync
uv run chemart list                     # every chemistry, as JSON
uv run chemart describe brusselator     # metadata + parameter JSON Schema
uv run chemart generate brusselator --format text
uv run chemart generate matrix-chemistry -p N=4 --seed 0 --format summary
```

```python
import chemart

# Discover
chemart.list_chemistries()
chemart.describe_chemistry("gard")

# Generate; every parameter is optional and defaults run in seconds
net = chemart.generate_network("gard", seed=0)

# Same seed, same network - always
assert chemart.generate_network("gard", seed=0).to_dict() == net.to_dict()
```

### For LLM agents

The whole interface is three functions, exposed as JSON-Schema tool specs:

```python
chemart.tool_definitions()               # list/describe/generate, ready to register
chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
```

There is also a bundled Claude Code skill at `.claude/skills/chemart/` with
task-oriented references and runnable examples.

## Documentation

The docs site is [MkDocs](https://www.mkdocs.org/) + Material, with a page for
every chemistry generated from the catalog so the site cannot drift from the
library.

```bash
uv run python tools/gen_catalog_pages.py     # regenerate the 98 chemistry pages
uv run mkdocs serve                          # http://127.0.0.1:8000
uv run mkdocs build                          # static site into site/
```

### Publishing to GitHub Pages

This repository has no remote yet. Once it has one:

```bash
git remote add origin git@github.com:<you>/<repo>.git
git push -u origin HEAD
```

Then in the repository's **Settings → Pages**, set *Source* to **GitHub
Actions**. The bundled [workflow](.github/workflows/docs.yml) regenerates the
catalog pages, builds with `--strict` and publishes on every push to `main` or
`chemart-library`; the site lands at `https://<you>.github.io/<repo>/`.

`uv run mkdocs gh-deploy` also works and pushes a built site to a `gh-pages`
branch, but the workflow is preferable — it regenerates the per-chemistry pages
from the catalog first, so the site cannot fall behind the library.

## Repository layout

```
chemart/
  api.py           list_chemistries / describe_chemistry / generate_network
  catalog.py       loader, validator, index generator
  network.py       the Network / Species / Reaction record
  kinetics.py      the rate-law vocabulary, and k -> c conversion
  expand.py        closure of a constructive rule (the generating operator)
  soup.py          well-stirred run that records which reactions fired
  cli.py           the `chemart` command
  chemistries/     98 modules, one `generate(p, rng)` each
  helpers/         explicit reaction syntax, parameter validation
catalog/
  SCHEMA.md        field definitions
  chemistries/     98 YAML entries - the specification
docs/              the documentation site
tests/
  test_contract.py the contract every chemistry must satisfy
  chemistries/     per-chemistry tests reproducing published results
```

## Testing

```bash
uv run pytest                              # fast tests
uv run pytest -q -m "slow or not slow"     # everything, including slow
uv run python -m chemart.catalog validate  # catalog invariants
```

The contract in `tests/test_contract.py` checks, for **every** chemistry, that a
zero-argument call finishes in under five seconds, that the record round-trips
through JSON, that the same seed reproduces the same network, that the
capabilities it computes are covered by what the catalog claims, and that bad
parameters raise an error naming the parameter. Per-chemistry tests then
reproduce published numbers — tables, counts, closures, steady states.

## Open questions

Implementation is finished; what remains is judgement, collected in
[`to_decide.md`](to_decide.md):

- the generative vs algorithmic classification of the chemistries;
- **licensing** — the Core War tests vendor two GPL-2 warriors from pMARS, and
  Chemart itself declares no license yet;
- two record conventions the genome-carrying chemistries disagree on;
- a proposed catalysed Michaelis-Menten rate law, and the `arrhenius` gap in the
  ODE helper;
- the entries whose primary sources could not be reached.

## Source

All section and page references are to Banzhaf & Yamamoto (2015). Bracketed
numbers in a `refs:` field are that book's bibliography entries.
