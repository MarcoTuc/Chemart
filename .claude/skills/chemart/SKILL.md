---
name: chemart
description: Generate, inspect, simulate and compare chemical reaction networks from Chemart — a library holding all 98 artificial chemistries in Banzhaf & Yamamoto's *Artificial Chemistries* (MIT Press, 2015), from the Brusselator and Oregonator to AlChemy, Tierra, Avida, matrix chemistry, RAF sets, P systems and DNA computing. Use this skill whenever the user wants an artificial chemistry or an abstract CRN: naming any catalogued chemistry, asking for a reaction network to analyse or simulate, comparing chemistries by what they supply (kinetics, energetics, conservation laws, space, compartments), exploring constructive or open-ended chemistries, or adding a new chemistry to the catalog. Also use it when someone asks for "a toy CRN", "an artificial chemistry", a model of autocatalysis / self-replication / origin-of-life / protocells, or wants to drive Chemart from an LLM tool loop. Prefer it over hand-rolling a reaction network from memory — the catalogued ones are sourced and tested, and a hand-rolled Brusselator will not match the book.
allowed-tools: Read Write Edit Bash
compatibility: Requires the Chemart repository and its uv environment. Run everything through `uv run` from the repo root (Python 3.12+, numpy/scipy/pyyaml, plus rdkit and ViennaRNA for a few entries). No network access needed for the built-in catalog; the optional Chemart Hub (shared chemistries and networks, ids like `ada/my-chem`) talks to `$CHEMART_HUB_URL`.
metadata:
  version: "1.0"
---

# Chemart — a mart of artificial chemistries

Chemart turns 98 published artificial chemistries into one uniform,
JSON-friendly interface. Each one generates a **chemical reaction network**:
species, reactions with reactant/product stoichiometry, optional rate laws,
and the provenance to reproduce it.

The point is comparability. The Brusselator and Tierra have nothing in common
as models, but both hand back the same record, so you can hold them side by
side, count conservation laws, simulate the ones with kinetics, and export any
of them the same way.

## Setup

Everything runs through `uv` from the repository root — the project pins its
own `./.venv`, and `pip` or a system Python will get the wrong interpreter.

This skill ships **inside** the Chemart repository, at
`<repo>/.claude/skills/chemart/`, so if you are working in Chemart you are
already in the right place — the repo root is three levels up from this file.

```bash
cd <repo>                    # the directory containing chemart/ and catalog/
uv sync                      # once
uv run chemart list          # sanity check: 98 entries
```

From another directory, point `uv` at the project instead of cd-ing:
`uv run --project <repo> chemart list`

## The whole interface is three functions

```python
import chemart

chemart.list_chemistries()                      # id, name, summary, fidelity, implemented
chemart.describe_chemistry("brusselator")       # metadata + JSON Schema of the parameters
net = chemart.generate_network("brusselator", seed=1)   # -> Network
```

`generate_network` takes `seed` plus any catalogued parameter as a keyword.
Every parameter is optional; defaults are chosen so a zero-argument call
finishes in seconds. The same seed always gives the same network.

The CLI mirrors it exactly, which is usually the fastest way to look:

```bash
uv run chemart list
uv run chemart describe matrix-chemistry
uv run chemart generate matrix-chemistry -p N=4 --seed 0 --format summary
uv run chemart generate brusselator --format text     # one reaction per line
uv run chemart generate brusselator --format json     # the full record
```

And for tool loops: `chemart.tool_definitions()` returns the three functions
as JSON-Schema specs, `chemart.call_tool(name, arguments)` executes them and
returns JSON-ready data. See `examples/llm_tools.py`.

## The network record, briefly

```python
net.summary()      # counts, status, capability tags, extras keys
net.to_text()      # "2 X + Y -> 3 X  [mass-action k=2.0]"
net.to_dict()      # plain JSON; Network.from_dict(d) round-trips exactly
ids, R, P = net.matrices()      # scipy.sparse; net stoichiometry S = P - R
net.provides       # what this network actually contains, computed from content
```

Three things here surprise people, so check them before you draw conclusions:

**`status` changes what the network means.** `complete` is the whole defined
network or a finished closure. `truncated` means a closure hit a size budget —
there is more chemistry you are not seeing, and raising the budget parameter
will find it. `observed` means these are the reactions that *fired in a
simulation*, each carrying a firing `count`; absence of a reaction is evidence
about that run, not about the chemistry.

**Zero reactions is legal.** `swarm-chemistry` defines motion rather than
transformation, and is the only entry whose default network has no reactions:
`net.reactions == []`, with the law in `net.extras["interaction_law"]`.
`summary()` says so rather than looking broken. Code that walks reactions
should handle the empty case.

**`net.provides` is computed; the catalog's `provides` is a claim.** The
contract guarantees computed ⊆ claimed, so an entry may advertise a capability
its *default* parameters don't exercise. `rbn` is the live example: it claims
`mass-conservation`, but its default network declares no conservation laws.
Trust `net.provides` for what you actually hold.

Full field-by-field spec: `references/network-record.md`.

## Common tasks

| You want to… | Go to |
|---|---|
| find the right chemistry out of 98 | `references/finding-chemistries.md` |
| understand every field of the record | `references/network-record.md` |
| simulate the dynamics / plot time courses | `references/simulating-dynamics.md` + `scripts/simulate.py` |
| compare chemistries, find conservation laws, study closures | `references/comparing-and-analysing.md` + `scripts/survey.py` |
| judge how much to trust an entry | `references/fidelity-and-trust.md` |
| add a new chemistry, or change one | `references/extending-chemart.md` |

Runnable examples live in `examples/` — start with `examples/quickstart.py`.

## Bundled scripts

Two things every user would otherwise rewrite:

```bash
# Integrate any network that carries rate constants, and print or plot it
uv run python .claude/skills/chemart/scripts/simulate.py brusselator --t-end 40 --seed 1

# Survey the catalog by capability, family, fidelity
uv run python .claude/skills/chemart/scripts/survey.py --provides rate-constants energies
uv run python .claude/skills/chemart/scripts/survey.py --family origin-of-life --verbose
```

`simulate.py` exists because the library deliberately ships **no simulator** —
Chemart's job is to hand you a correct network, not to be a solver. The
repository's ODE helper lives in `tests/` and is not importable as
`chemart.odes`, so this script is the supported way to integrate one. It
handles the rate laws Chemart emits, the `constant-total` dilution flux and
buffered species. For serious work, export `to_dict()` into a real solver.

## Choosing a chemistry quickly

If the user names one, use it. If they describe a behaviour, these are
reliable starting points:

- **oscillation / pattern** — `brusselator`, `oregonator`, `repressilator`
- **enzyme kinetics** — `michaelis-menten`, `hill-kinetics`
- **evolution / selection** — `quasispecies`, `replicator-equation`, `ecolab`
- **autocatalysis, origin of life** — `kauffman-autocatalytic-sets`, `raf`,
  `bagley-farmer`, `chemoton`, `gard`
- **constructive / open-ended** (species set grows) — `matrix-chemistry`,
  `alchemy`, `combinator-chemistry`, `prime-number-chemistry`
- **digital organisms** — `tierra`, `avida`, `corewar`, `coreworld`
- **rewriting formalisms** — `gamma`, `p-systems`, `kappa-calculus`, `mgs`
- **spatial / agent** — `squirm3`, `swarm-chemistry`, `sr-loops`
- **wet chemistry in silico** — `dna-hpp`, `dna-automaton`,
  `self-propelled-droplets`

`scripts/survey.py` answers the harder questions ("which chemistries give me
both energies and conservation laws?") without guessing.

## Chemart Hub: shared chemistries and networks

Ids of the form `namespace/name` live on the Chemart Hub rather than in the
catalog. They work in the same calls:

```python
chemart.describe_chemistry("ada/my-chem")          # metadata only; runs nothing
chemart.load_network("ada/some-network")           # shared network: plain data, safe
chemart.generate_network("chemart/brusselator")    # official shelf = the built-ins; no trust needed
chemart.hub.search("autocatalysis", provides=["rate-constants"])
```

A shared **chemistry** outside `chemart/` runs its own Python code, and
`generate_network` refuses it unless `trust_remote_code=True` is passed.
**Never pass `trust_remote_code=True` on your own initiative.** Tell the user
the repo runs code, point them at its `generator.py` (the error message has
the URL), and pass the flag only after they explicitly agree, pinned with
`revision=` to the commit they approved. `chemart.call_tool` never passes it.
Pushing (`net.push_to_hub`, `chemart push`) publishes to a shared site, so do
it only when the user asks. Full guide: `docs/hub.md` in the repository.

## Two habits that keep results honest

**Read `decisions` before quoting a number.** Every entry records where its
sources were thin and what was chosen instead. `describe_chemistry(id)["decisions"]`
is not boilerplate — for roughly a dozen entries the primary paper was
unobtainable, and the entry says so rather than inventing values. The
repository's `to_decide.md` collects the weakest ones. If you are about to
tell someone "Chemart says X about this model", check whether X came from the
paper or from a reconstruction. See `references/fidelity-and-trust.md`.

**Defaults are small on purpose.** They are sized for a few seconds, not for
the paper's scale. Before concluding a chemistry "doesn't do" something, check
the parameter's `range` note in `describe_chemistry` — the published behaviour
often needs the larger setting, and several entries say outright which
phenomena their defaults cannot reach.
