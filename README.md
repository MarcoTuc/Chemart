```text
     °       H     H       O
    ∘  °      ╲   ╱       ╱ ╲
     °         C═C       H   H      ▗▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▖
    ° ∘       ╱   ╲                 ▐▓▓▓▓░░░░▓▓▓▓░░░░▓▓▓▓░░░░▓▓▓▓░░░░▓▓▌
    ┌─┐      H     H                ╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱
    │ │      ░░░░░░╗░░╗  ░░╗░░░░░░░╗███╗   ███╗ █████╗ ██████╗ ████████╗
    │ │     ░░╔════╝░░║  ░░║░░╔════╝████╗ ████║██╔══██╗██╔══██╗╚══██╔══╝
   ╱ ° ╲    ░░║     ░░░░░░░║░░░░░╗  ██╔████╔██║███████║██████╔╝   ██║
  ╱ ∘ ° ╲   ░░║     ░░╔══░░║░░╔══╝  ██║╚██╔╝██║██╔══██║██╔══██╗   ██║
 ╱≈≈≈≈≈≈≈╲  ╚░░░░░░╗░░║  ░░║░░░░░░░╗██║ ╚═╝ ██║██║  ██║██║  ██║   ██║
 ╰───────╯   ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝
            └───────────── chem ──────────────┘
                                    └────────────── mart ──────────────┘

              ╔═════════════════════════════════════════╗
              ║ [A] [B] [X] [Y] [λ] [RNA] [01101] [●─●] ║
              ╚═════════════════════════════════════════╝

                one stop shop for artificial chemistry
                   chemistries in stock · open 24/7
```

# Chemart — the one stop shop for artificial chemistry

Chemart makes the artificial chemistries of the literature available behind one
interface, so you can generate a chemistry's reaction network, look at it,
simulate it, and compare it with any other — without reimplementing the papers.

The chemistries come from the field's reference survey — Wolfgang Banzhaf &
Lidia Yamamoto, *Artificial Chemistries* (MIT Press, 2015) — and every one is
catalogued **and implemented**, with tests reproducing published results
wherever any exist. Entries that are not chemistry, such as the artificial-life
systems, are kept in an [archive](#the-archive).

```python
import chemart

chemart.list_chemistries()                      # the whole catalog
net = chemart.generate_network("brusselator", seed=1)
print(net.to_text())
# A -> X  [mass-action k=1.0]
# B + X -> Y + D  [mass-action k=1.0]
# 2 X + Y -> 3 X  [mass-action k=1.0]
# X -> E  [mass-action k=1.0]
```

🏪 **Chemart Hub:** share your own chemistries and networks, and load anyone's
by name: `chemart.generate_network("ada/my-chem", trust_remote_code=True)`.
See [Chemart Hub](#chemart-hub).

📖 **Full documentation:** <https://marcotuc.github.io/Chemart/>, a browsable
site with a page for every chemistry. To browse it locally, run
`uv run mkdocs serve` and open <http://127.0.0.1:8000>; see
[Documentation](#documentation).

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
the Brusselator (three coupled differential equations) and AlChemy (λ-terms that
react by applying one to the other) are both artificial chemistries, and they
share almost nothing operationally.

Chemart's bet is that they can still share an **output**. Whatever a chemistry is
internally, it can hand back a chemical reaction network: which species exist,
which reactions consume and produce them, and — when the model defines them —
the rates. That common output is what makes incomparable models comparable.

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
exactly — a property tested for every chemistry. From there:

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

Many of these chemistries are **constructive**: the species set is open and
grows as reactions produce new molecules. For those, the interesting object is
a closure, and `truncated` is a normal answer rather than a failure.

### Provenance is part of the data

Reproducing published models honestly means admitting where the sources ran out.
Every entry carries a `fidelity`:

| | meaning |
|---|---|
| `book` | implemented exactly as the book specifies |
| `book+decisions` | the book left gaps; each filled choice is in `decisions` |
| `reconstructed` | built from the original papers, cited in `sources` |

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

Browse the chemistries in [`docs/CATALOG.md`](docs/CATALOG.md) or on the docs
site. Both are generated from the catalog, and count what it holds: how many
chemistries, how many are constructive, how many carry their own kinetics or
energetics, and each one's fidelity.

### The archive

Some entries are kept out of the catalog: the artificial-life systems (Tierra,
Avida, Core War, Coreworld, Swarm Chemistry), since Chemart is about chemistry,
and others set aside on review. Each is marked `archived:` in its YAML and keeps
its generator, tests and page. `list_chemistries()`, `chemart list` and the LLM
tools leave them out; `list_chemistries(include_archived=True)` and
`chemart list --all` include them, and `generate_network(id)` still runs any of
them. Deleting the `archived:` line brings an entry back.

## Quickstart

The environment is managed with [uv](https://docs.astral.sh/uv/); run everything
through `uv run` so you get the project's own `./.venv`.

```bash
uv sync
uv run chemart list                     # every chemistry, as JSON (--all adds the archive)
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

### A notebook to start from

[`examples/simulating-chemistries.ipynb`](examples/simulating-chemistries.ipynb)
simulates a given network (the Brusselator, integrated as rate equations) and
a generator (the prime-number chemistry, whose run generates its network).
Open it with `uv run --group notebooks jupyter lab examples/`.

### For LLM agents

The whole interface is three functions, exposed as JSON-Schema tool specs:

```python
chemart.tool_definitions()               # list/describe/generate, ready to register
chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
```

There is also a bundled Claude Code skill at `.claude/skills/chemart/` with
task-oriented references and runnable examples.

## Chemart Hub

The catalog is the book. The **Chemart Hub** is the rest of the shop: a website
where people publish chemistries and reaction networks, and a client built into
`chemart` that loads them by name, the way `transformers` loads models.

```python
import chemart

chemart.generate_network("chemart/brusselator")        # the official shelf: the built-ins
chemart.load_network("ada/brusselator-b35")            # a shared network: plain data
chemart.generate_network("ada/hypercycle-lite",        # a shared chemistry: runs its code,
                         revision="f2634dc6b1d4",      # so read it, pin it,
                         trust_remote_code=True)       # and say so

net.push_to_hub("you/my-network")                      # share a network
```

```bash
uv run chemart login                 # paste a token from the hub's Settings → Tokens
uv run chemart new my-chem           # a working generator skeleton
uv run chemart check my-chem         # the contract, run locally
uv run chemart push my-chem          # -> you/my-chem
uv run chemart search autocatalysis --provides rate-constants
```

There are two kinds of repo. A **chemistry** holds a catalog entry plus the
`generate(p, rng)` that builds it, so anyone can run it with their own
parameters. A **network** holds one reaction network as JSON. Every push is a
content-addressed commit, so `revision=` pins exactly what runs, and every
generated network records `namespace/name@commit` as its provenance. Pushes are
held to the same contract as the built-in catalog: the entry validates, defaults
run in seconds, and the same seed gives the same network. The hub never runs
uploaded code. You run it, on your own machine, and only with
`trust_remote_code=True`.

The server is in this repository (`hub/`, the `chemart-hub` package):

```bash
uv sync --all-packages
uv run --package chemart-hub chemart-hub init
uv run --package chemart-hub chemart-hub create-user you --admin
uv run --package chemart-hub chemart-hub serve           # http://127.0.0.1:8000
```

Full guide: [`docs/hub.md`](docs/hub.md).

## Documentation

The docs site is [MkDocs](https://www.mkdocs.org/) + Material, with a page for
every chemistry generated from the catalog so the site cannot drift from the
library.

```bash
uv run python tools/gen_catalog_pages.py     # regenerate the chemistry pages
uv run mkdocs serve                          # http://127.0.0.1:8000
uv run mkdocs build                          # static site into site/
```

### Publishing to GitHub Pages

The site is published at <https://marcotuc.github.io/Chemart/>. The
[workflow](.github/workflows/docs.yml) regenerates the catalog pages, builds
with `--strict` and deploys on every push to `main` or `chemart-library`, so
merging is publishing. `uv sync` installs `mkdocs` by default (the `docs`
dependency group is in `default-groups`).

## Repository layout

```
chemart/
  api.py           list_chemistries / describe_chemistry / generate_network
  catalog.py       loader, validator, index generator
  network.py       the Network / Species / Reaction record
  kinetics.py      the rate-law vocabulary, and k -> c conversion
  expand.py        closure of a constructive rule (the generating operator)
  soup.py          well-stirred run that records which reactions fired
  contract.py      the checks every generator must pass (tests, `chemart check`, push)
  cli.py           the `chemart` command
  chemistries/     one module per chemistry, each defining `generate(p, rng)`
  helpers/         explicit reaction syntax, parameter validation
  hub/             Chemart Hub client: ids, cache, trust gate, push/load (stdlib only)
hub/               the Chemart Hub server (package `chemart-hub`)
  src/chemart_hub/ FastAPI app, SQLite + blob store, web UI, `chemart-hub` CLI
  tests/           API, web, and end-to-end client <-> server tests
catalog/
  SCHEMA.md        field definitions
  chemistries/     one YAML entry per chemistry - the specification
  explainers/      the prose of each chemistry's documentation page
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
uv sync --all-packages && uv run --package chemart-hub pytest hub/tests   # the hub
```

The contract (`chemart/contract.py`, run by `tests/test_contract.py`) checks, for **every** chemistry, that a
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
