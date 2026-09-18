# API and CLI reference

## Python API

### `chemart.list_chemistries()`

Every chemistry in the catalog.

```python
[{"id": "acgp", "name": "...", "summary": "...",
  "implemented": True, "fidelity": "reconstructed"}, ...]
```

### `chemart.describe_chemistry(chemistry)`

Everything known about one chemistry, JSON-ready.

| key | contents |
|---|---|
| `id`, `name`, `aliases`, `origin` | identity |
| `intuition` | plain-language explanation of how the chemistry works and what the idea behind it is — the place to start when you meet an entry for the first time |
| `family`, `kind`, `constructive` | classification |
| `implemented`, `fidelity` | status |
| `book`, `refs`, `sources` | provenance |
| `decisions` | every gap in the sources and how it was resolved |
| `molecules`, `reactions`, `reactor` | the (S, R, A) triple |
| `provides` | capability tags the chemistry claims |
| `phenomena` | what the model is known to produce |
| `notes` | anything else |
| `params` | JSON Schema of the parameters |

`phenomena` is the most underused field: it is the list of behaviours the
literature reports, which is usually what someone actually wants when they ask
for "a chemistry that does X".

### `chemart.generate_network(chemistry, seed=None, **params)`

Runs the generator and returns a [`Network`](record.md). Raises `ValueError`
for an unknown chemistry (suggesting close ids), an unknown parameter (listing
the valid ones) or an out-of-range value (naming the parameter and its
meaning).

### `chemart.tool_definitions()` and `chemart.call_tool(name, arguments)`

The three functions as JSON-Schema tool specs, and an executor for them.
`call_tool` returns JSON-ready data — `generate_network` comes back as the plain
record dict rather than a `Network` object.

```python
tools = chemart.tool_definitions()
[t["name"] for t in tools]
# ['list_chemistries', 'describe_chemistry', 'generate_network']

chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
```

The `chemistry` argument of `generate_network` is an **enum** of every
implemented id, so a model cannot ask for one that does not exist.

!!! tip "Keep records small when feeding them back to a model"
    Chemistries that carry genomes or sequences in `Species.structure` can
    produce large records — `aevol`'s default is about 1.2 MB of JSON. Drop
    structures, or summarise, before returning one to a model context.

### Also exported

`chemart.Network`, `chemart.Species`, `chemart.Reaction` — see
[The network record](record.md).

## Command line

```bash
uv run chemart list
uv run chemart describe <chemistry>
uv run chemart generate <chemistry> [-p NAME=VALUE ...] [--seed N]
                                    [--format summary|text|json]
```

`-p` may be repeated; values are parsed as JSON when possible and taken as
strings otherwise. Errors go to stderr and exit with status 2.

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

`Chemistry` carries `id`, `name`, `family`, `kind`, `constructive`, `provides`,
`fidelity`, `sources`, `decisions`, `phenomena`, `book`, `refs`, `params`, plus
the properties `implemented`, `module`, `tiers` and the method
`params_by_role(role)`.

Parameter objects carry `name`, `type`, `default`, `min`, `max`, `choices`,
`meaning`, `role` and `range`.

## Helpers

| import | for |
|---|---|
| `chemart.helpers.explicit.network` | building a network from written-down reactions |
| `chemart.helpers.params` | validating list/dict parameters (`vector`, `square_matrix`, `edges`, `apportion`) |
| `chemart.expand.expand` | closure of a constructive rule |
| `chemart.soup.soup` | well-stirred run that records which reactions fired |
| `chemart.kinetics.RATE_LAWS` | the rate-law vocabulary |
| `chemart.kinetics.k_to_c` | macroscopic → mesoscopic rate constant |

## Testing

```bash
uv run pytest                                # fast tests only
uv run pytest -q -m "slow or not slow"       # everything
CHEMART_ONLY=brusselator uv run pytest -q tests/test_contract.py
```

`CHEMART_ONLY` narrows the parametrized contract tests to one or more
comma-separated ids. The coverage gate that checks all 98 are implemented
deliberately ignores it.
