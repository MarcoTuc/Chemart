---
name: chemart
description: Generate, inspect, simulate, evolve, measure and compare chemical reaction networks from Chemart — a library of the artificial chemistries in Banzhaf & Yamamoto's *Artificial Chemistries* (MIT Press, 2015), from the Brusselator and Oregonator to AlChemy, matrix chemistry, RAF sets, P systems and DNA computing, with artificial-life systems such as Tierra and Avida kept in an archive. Use this skill whenever the user wants an artificial chemistry or an abstract CRN: naming any catalogued chemistry, asking for a reaction network to analyse or simulate, comparing chemistries by what they supply (kinetics, energetics, conservation laws, space, compartments) or by structural and dynamical measures, running Turing gases (AlChemy, BFF, combinatory chemistry) in evolutionary time, exploring constructive or open-ended chemistries, or adding a new chemistry to the catalog. Also use it when someone asks for "a toy CRN", "an artificial chemistry", a model of autocatalysis / self-replication / origin-of-life / protocells, or wants to drive Chemart from an LLM tool loop. Prefer it over hand-rolling a reaction network from memory — the catalogued ones are sourced and tested, and a hand-rolled Brusselator will not match the book.
allowed-tools: Read Write Edit Bash
compatibility: Requires the Chemart repository and its uv environment. Run everything through `uv run` from the repo root (Python 3.12+, numpy/scipy/pyyaml, plus rdkit and ViennaRNA for a few entries). No network access needed for the built-in catalog; the optional Chemart Hub (shared chemistries and networks, ids like `ada/my-chem`) talks to `$CHEMART_HUB_URL`.
metadata:
  version: "1.1"
---

# Chemart — a mart of artificial chemistries

Chemart turns published artificial chemistries into one uniform,
JSON-friendly interface. Each one generates a **chemical reaction network**:
species, reactions with reactant/product stoichiometry, optional rate laws,
and the provenance to reproduce it.

The point is comparability. The Brusselator and AlChemy have nothing in common
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
uv run chemart list          # sanity check: a JSON list of the catalog
```

From another directory, point `uv` at the project instead of cd-ing:
`uv run --project <repo> chemart list`

## Three types, two faces

Every chemistry is one of three **types** (`describe_chemistry(id)["type"]`):

- **given** — a reaction network written down (Brusselator, Michaelis–Menten).
  Choose rates and an initial state, then simulate it.
- **generator** — an algorithm computes the network from its arguments
  (Kauffman sets, RAF, random catalytic networks). Once built it is treated as
  given; its measures change with the arguments.
- **gas** — a Turing gas: structured molecules react by a procedure (AlChemy,
  BFF, combinatory chemistry). The network is not known in advance; the soup
  *evolves*, and its measures change in evolutionary time.

A chemistry's module has one or two **faces**, listed in
`describe_chemistry(id)["faces"]`: `generate` returns a network, `evolve` runs
its process frame by frame in its own `clock` (collisions, epochs,
generations…). Many gases have both: `generate` gives a network directly
(often a closure), `evolve` runs the soup. A gas with only `evolve` still
works with `generate_network`, which runs the process to the end and returns
what fired.

## The interface

```python
import chemart
from chemart import measures, simulate

chemart.list_chemistries()                      # id, name, type, summary, fidelity, implemented
chemart.describe_chemistry("brusselator")       # metadata, type, faces, clock, JSON Schema of params
net = chemart.generate_network("brusselator", seed=1)   # -> Network

traj = simulate.ode(net, t_end=40)              # rate equations -> Trajectory
path = simulate.ssa(net, t_end=40, volume=100, seed=1)  # one Gillespie path

run = chemart.evolve("alchemy", seed=1)         # a gas in its own time -> Trajectory
for frame in chemart.evolve_frames("bff", seed=1): ...  # the same, live

chemart.measure(net)                            # {measure: value}, cheap tier by default
measures.over(run, ["shannon", "n_species"], window=5)  # measures in evolutionary time
```

`generate_network` and `evolve` take `seed` plus any catalogued parameter as a
keyword. Every parameter is optional; defaults are chosen so a zero-argument
call finishes in seconds. The same seed always gives the same network or run.
`describe_chemistry` gives each face's parameters (`params`, `evolve_params`);
one passed to the wrong face raises an error that says where it belongs.

The CLI mirrors it exactly, which is usually the fastest way to look:

```bash
uv run chemart list
uv run chemart describe matrix-chemistry
uv run chemart generate matrix-chemistry -p N=4 --seed 0 --format summary
uv run chemart generate brusselator --format json     # the full record
uv run chemart simulate brusselator --t-end 40        # also --method ssa --volume 100
uv run chemart evolve alchemy --seed 1 --track shannon --track n_species
uv run chemart measure raf --seed 1 --cost moderate --why
```

For tool loops, `chemart.tool_definitions()` returns six JSON-Schema tools —
`list_chemistries`, `describe_chemistry`, `generate_network`,
`simulate_network`, `evolve_chemistry`, `measure_network` — and
`chemart.call_tool(name, arguments)` executes them, returning JSON-ready data.
The simulate and evolve tools cap their output (12 species, 40 time points).
See `examples/llm_tools.py`.

For a person rather than a model, the **simulation pit** is a local web app:
pick a chemistry, set parameters, rates and initial state, run ODE, SSA or
evolve, and watch concentrations and measures live.

```bash
uv sync --all-packages       # once: the pit lives in the hub package
uv run chemart-hub pit       # opens http://127.0.0.1:8765, local only
```

## The network record, briefly

```python
net.summary()      # counts, status, capability tags, extras keys
net.to_text()      # "2 X + Y -> 3 X  [mass-action k=2.0]"
net.to_dict()      # plain JSON; Network.from_dict(d) round-trips exactly
ids, R, P = net.matrices()      # scipy.sparse; net stoichiometry S = P - R
net.provides       # what this network actually contains, computed from content
```

Four things here surprise people, so check them before you draw conclusions:

**`status` changes what the network means.** `complete` is the whole defined
network or a finished closure. `truncated` means a closure hit a size budget —
there is more chemistry you are not seeing, and raising the budget parameter
will find it. `observed` means these are the reactions that *fired in a
simulation*, each carrying a firing `count`; absence of a reaction is evidence
about that run, not about the chemistry.

**A gas has two kinds of network.** Its generated network is either a closure
(`complete`/`truncated`) or the record of one run (`observed`), and the two
answer different questions. Compare gases by population measures over a run
first; network measures of an observed network mix the chemistry with how it
was sampled.

**Zero reactions is legal.** `swarm-chemistry` (archived) defines motion rather
than transformation, and is the only entry whose default network has no reactions:
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
| find the right chemistry, or read about the archive | `references/finding-chemistries.md` |
| understand every field of the record | `references/network-record.md` |
| simulate a network, or evolve a gas | `references/simulating-dynamics.md` (`chemart.simulate`, `chemart.evolve`) |
| measure and compare chemistries, find conservation laws, study closures | `references/comparing-and-analysing.md` + `scripts/survey.py` |
| judge how much to trust an entry | `references/fidelity-and-trust.md` |
| add a new chemistry, or change one | `references/extending-chemart.md` |

Runnable examples live in `examples/` — start with `examples/quickstart.py`.

## Surveying and measuring

```bash
# Survey the catalog by type, capability, family, fidelity (bundled script)
uv run python .claude/skills/chemart/scripts/survey.py --type gas --verbose
uv run python .claude/skills/chemart/scripts/survey.py --provides rate-constants energies

# Measures of a chemistry's network, a saved network, or a saved run
uv run chemart measure kauffman-autocatalytic-sets --seed 1 --cost moderate
uv run chemart evolve bff --seed 1 --format json > run.json && uv run chemart measure run.json
```

Measures are registered by section (size, stoichiometry, graph, organisation,
growth, kinetics, dynamics, robustness, information) and by cost: `cheap` runs
by default, `moderate` and `exponential` on request, with node limits that
`force=True` lifts. A measure that does not apply is left out, never filled
in; `measures.applicable(obj)` says why. `docs/guide/measures.md` in the
repository has the full list with meanings and references.

## Choosing a chemistry quickly

If the user names one, use it. If they describe a behaviour, these are
reliable starting points:

- **oscillation / pattern** — `brusselator`, `oregonator`, `repressilator`
- **enzyme kinetics** — `michaelis-menten`, `hill-kinetics`
- **evolution / selection** — `jain-krishna`, `random-catalytic-networks`
- **autocatalysis, origin of life** — `kauffman-autocatalytic-sets`, `raf`,
  `bagley-farmer`, `chemoton`, `gard`
- **constructive / open-ended** (species set grows) — `matrix-chemistry`,
  `alchemy`, `combinator-chemistry`, `prime-number-chemistry`
- **Turing gases** (run with `chemart.evolve`) — `alchemy`, `bff`,
  `combinatory-chemistry`, `rbn-world`, `typogenetics`, `stringmol`
- **rewriting formalisms** — `gamma`, `p-systems`, `kappa-calculus`, `mgs`
- **spatial / agent** — `squirm3`, `sr-loops`, `flow-ac`
- **wet chemistry in silico** — `dna-hpp`, `dna-automaton`,
  `self-propelled-droplets`

Digital organisms (`tierra`, `avida`, `corewar`, `coreworld`) and `swarm-chemistry`
are artificial life, kept in the archive: run them by id only when asked for.

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
`revision=` to the commit they approved. `chemart.call_tool` never passes it, and the tools never read local files.
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
the paper's scale (a BFF default soup is 128 programs, not 131,072). Before concluding a chemistry "doesn't do" something, check
the parameter's `range` note in `describe_chemistry` — the published behaviour
often needs the larger setting, and several entries say outright which
phenomena their defaults cannot reach.
