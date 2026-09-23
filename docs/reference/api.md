# API and CLI reference

## Python API

### `chemart.list_chemistries()`

Every chemistry in the catalog.

```python
[{"id": "acgp", "name": "...", "summary": "...",
  "implemented": True, "fidelity": "reconstructed", "type": "generator"}, ...]
```

### `chemart.describe_chemistry(chemistry)`

Everything known about one chemistry, JSON-ready.

| key | contents |
|---|---|
| `id`, `name`, `aliases`, `origin` | identity |
| `intuition` | plain-language explanation of how the chemistry works and what the idea behind it is — the place to start when you meet an entry for the first time |
| `family`, `kind`, `constructive` | classification |
| `type` | `given`, `generator` or `gas`: where the network comes from ([How it works](../concepts.md#three-types-of-chemistry)) |
| `faces`, `clock` | what can be called: `["generate"]`, `["evolve"]` or both; the unit of evolve time |
| `implemented`, `fidelity` | status |
| `book`, `refs`, `sources` | provenance |
| `decisions` | every gap in the sources and how it was resolved |
| `molecules`, `reactions`, `reactor` | the (S, R, A) triple |
| `provides` | capability tags the chemistry claims |
| `phenomena` | what the model is known to produce |
| `notes` | anything else |
| `params` | JSON Schema of the parameters of `generate_network` |
| `evolve_params` | JSON Schema of the parameters of `evolve` (with an evolve face only) |

`phenomena` is the most underused field: it is the list of behaviours the
literature reports, which is usually what someone actually wants when they ask
for "a chemistry that does X".

### `chemart.generate_network(chemistry, seed=None, *, revision=None, trust_remote_code=False, **params)`

Runs the generator and returns a [`Network`](record.md). Raises `ValueError`
for an unknown chemistry (suggesting close ids), an unknown parameter (listing
the valid ones) or an out-of-range value (naming the parameter and its
meaning).

`chemistry` may also be a [Chemart Hub](../hub.md) id, `namespace/name`,
optionally pinned as `namespace/name@<commit>`. `revision` selects the commit
(`"main"`, a full id or a unique prefix of at least 7 hex digits). A hub
chemistry that ships its own code raises `ValueError` unless
`trust_remote_code=True`; the official `chemart/` repos run the built-in
generators and need no trust.

### `chemart.evolve(chemistry, seed=None, *, every=1, **params)`

Runs the chemistry's process and returns a
[`Trajectory`](record.md#the-trajectory-record) with `method="evolve"` and the
chemistry's `clock`. Only chemistries with an evolve face can be evolved; the
error says what to call instead. `every` keeps one frame in `every`, adding up
the reactions fired in between. See [Evolving a chemistry](../guide/evolving.md).

`chemart.evolve_frames(chemistry, seed=None, *, every=1, **params)` is the same
run as a generator of frames, for live use; its return value is the observed
network.

### `chemart.simulate`

`ode(net, t_end, *, rates=None, x0=None, points=200, solver="LSODA", seed=None, temperature=None, gas_constant=None, fill_only=False)`
and `ssa(net, t_end, *, rates=None, x0=None, volume=1, avogadro=1, points=200, seed=None, max_events=10**6, ...)`
return a `Trajectory`; `assign(net, rates=, x0=, rng=, fill_only=)` returns the
rated network. A network that cannot be simulated raises `NotSimulable`,
naming what is missing. See [Simulating dynamics](../guide/simulating.md).

### `chemart.measure(obj, names=None, *, cost="cheap", food=None, seed=0, force=False)`

Measures of a network, a trajectory or a population state, as `{name: value}`.
Measures that do not apply are left out; `chemart.measures.applicable` says
why. The module `chemart.measures` also has `describe`, `track`, `over`,
`zscores`, `sweep` and `scaling`; see [Measures for comparison](../guide/measures.md).

### `chemart.load_network(repo_id, revision=None)`

Loads a network shared on the hub (a network repo). It is plain data, so no
trust is needed. `Network.push_to_hub(repo_id, title=, description=, tags=,
license=, message=)` publishes one.

### `chemart.hub`

`search(query, repo_type=, family=, provides=[...], tag=, author=, sort=)`,
`repo_info(id)`, `snapshot_download(id, revision=)`, `push_generator(folder,
id=None)`, `push_network(net, id)`, `check(folder)`, `new(folder)`,
`login(token=None)`, `logout()`, `whoami()`. Errors are `HubError` subclasses:
`HubConnectionError`, `HubAuthError`, `RepoNotFoundError`,
`HubConflictError`, `HubValidationError` (whose `problems` lists every
reason a push was refused). See [Chemart Hub](../hub.md).

### `chemart.tool_definitions()` and `chemart.call_tool(name, arguments)`

The public functions as JSON-Schema tool specs, and an executor for them.
`call_tool` returns JSON-ready data — `generate_network` comes back as the plain
record dict rather than a `Network` object.

```python
tools = chemart.tool_definitions()
[t["name"] for t in tools]
# ['list_chemistries', 'describe_chemistry', 'generate_network',
#  'simulate_network', 'evolve_chemistry', 'measure_network']

chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
chemart.call_tool("simulate_network", {"chemistry": "brusselator", "method": "ssa", "volume": 100})
chemart.call_tool("evolve_chemistry", {"chemistry": "alchemy", "seed": 1, "track": ["shannon"]})
chemart.call_tool("measure_network", {"chemistry": "raf", "cost": "moderate"})
```

The simulate, evolve and measure tools keep their answers small enough for a model's context:
`simulate_network` and `evolve_chemistry` return at most 12 species and 40
time points per series (`chemart.api.TOOL_SPECIES`, `TOOL_POINTS`), and
`evolve_chemistry` returns the size of the observed network rather than the
network. `measure_network` returns the values and, for the measures left out,
why they do not apply.

The `chemistry` argument is an **enum** of every implemented id (of every
chemistry with an evolve face, for `evolve_chemistry`), so a model cannot ask
for one that does not exist. A tool call never runs code from the hub:
`call_tool` refuses `trust_remote_code` among the parameters. It never reads
files either: `simulate_network` refuses rates or states given as a path.

!!! tip "Keep records small when feeding them back to a model"
    Chemistries that carry genomes or sequences in `Species.structure` can
    produce large records — `aevol`'s default is about 1.2 MB of JSON. Drop
    structures, or summarise, before returning one to a model context.

### Also exported

`chemart.Network`, `chemart.Species`, `chemart.Reaction`, `chemart.Trajectory`,
`chemart.Frame` — see [The network record](record.md).

## Command line

```bash
uv run chemart list
uv run chemart describe <chemistry> [--revision REV]
uv run chemart generate <chemistry> [-p NAME=VALUE ...] [--seed N]
                                    [--format summary|text|json]
                                    [--revision REV] [--trust-remote-code]
uv run chemart simulate <chemistry> [-p NAME=VALUE ...] [--seed N] [--method ode|ssa]
                                    [--t-end T] [--points N] [--volume V]
                                    [--rates SPEC] [--x0 SPEC] [--fill-only]
                                    [--species S ...] [--format table|csv|json]
uv run chemart evolve <chemistry> [-p NAME=VALUE ...] [--seed N] [--every N]
                                  [--track MEASURE ...] [--window N]
                                  [--species [S ...]] [--format table|csv|json]
uv run chemart measure <chemistry | file.json> [-p NAME=VALUE ...] [--seed N]
                                  [--names M ...] [--cost cheap|moderate|exponential]
                                  [--force] [--why] [--format table|json]
```

`-p` may be repeated; values are parsed as JSON when possible and taken as
strings otherwise. `--rates` and `--x0` take a number, a JSON distribution or
table, or a `.json`/`.csv` file. `evolve` prints the population, the
chemistry's observables and each `--track` measure per frame (`--window 0`:
the cumulative network); `--format json` writes the whole trajectory, which
`measure` reads back. `measure --why` also lists the measures that do not apply,
with the reason. Errors go to stderr and exit with status 2 (status 3 for
errors talking to the hub).

The [simulation pit](../hub.md#the-simulation-pit) is a local app for the same
work:

```bash
uv run chemart-hub pit [--port 8765] [--no-browser]
```

Chemart Hub commands:

```bash
uv run chemart login [--token T]      # save a token for $CHEMART_HUB_URL
uv run chemart whoami | logout
uv run chemart search [QUERY] [--type generator|network] [--provides TAG ...] [--tag T] [--author A]
uv run chemart download <ns/name> [--revision REV]     # prints the local folder
uv run chemart new <folder> [--id ID] [--name NAME]    # a generator skeleton
uv run chemart check <folder>                          # the contract, locally
uv run chemart push <folder> [ns/name] [-m MESSAGE]
uv run chemart push-network <file.json> <ns/name> [--title --description --license --tag]
```

## Catalog tooling

```bash
uv run python -m chemart.catalog validate           # invariants; 0 problems expected
uv run python -m chemart.catalog validate --only <id>
uv run python -m chemart.catalog status             # implementation + fidelity counts
uv run python -m chemart.catalog index              # regenerate docs/CATALOG.md
uv run python -m chemart.catalog show <id>          # human-readable entry
uv run python -m chemart.catalog query --provides rate-constants
uv run python -m chemart.catalog query --family origin-of-life
```

### From Python

```python
from chemart.catalog import load

for c in load():                  # parses once and caches; each call gets its own list
    if c.constructive and "rate-constants" in c.provides:
        print(c.id, c.fidelity, [p.name for p in c.params_by_role("structural")])
```

`Chemistry` carries `id`, `name`, `family`, `kind`, `type`, `clock`, `constructive`, `provides`,
`fidelity`, `sources`, `decisions`, `phenomena`, `book`, `refs`, `params`, plus
the properties `implemented`, `module`, `tiers` and the method
`params_by_role(role)`.

Parameter objects carry `name`, `type`, `default`, `min`, `max`, `choices`,
`meaning`, `role`, `range` and `face` (`generate`, `evolve`, or None when both
faces use it).

## Helpers

| import | for |
|---|---|
| `chemart.helpers.explicit.network` | building a network from written-down reactions |
| `chemart.helpers.params` | validating list/dict parameters (`vector`, `square_matrix`, `edges`, `apportion`) |
| `chemart.expand.expand` | closure of a constructive rule |
| `chemart.soup.stir` | well-stirred soup as a generator of frames, for an evolve face |
| `chemart.soup.Tally` | counts of the reactions that fired, overall and per frame |
| `chemart.soup.soup` | `stir` run to the end, returning what fired |
| `chemart.kinetics.RATE_LAWS` | the rate-law vocabulary |
| `chemart.kinetics.k_to_c` | macroscopic → mesoscopic rate constant |

## Testing

```bash
uv run pytest                                # fast tests only
uv run pytest -q -m "slow or not slow"       # everything
CHEMART_ONLY=brusselator uv run pytest -q tests/test_contract.py
```

`CHEMART_ONLY` narrows the parametrized contract tests to one or more
comma-separated ids. The coverage gate that checks every catalogued chemistry is implemented
deliberately ignores it.
