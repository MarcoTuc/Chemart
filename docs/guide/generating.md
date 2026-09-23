# Generating networks

## The three functions

```python
import chemart

chemart.list_chemistries()                  # [{id, name, summary, implemented, fidelity, type}]
chemart.describe_chemistry("brusselator")   # metadata + JSON Schema of the parameters
net = chemart.generate_network("brusselator", seed=1)
```

These three find a chemistry and give you its network. Anything else about the
network is a property or method of the returned
[network record](../reference/record.md). What you do next depends on the
chemistry's [type](../concepts.md#three-types-of-chemistry):
[simulate](simulating.md) a given or generated network,
[evolve](evolving.md) a gas, and [measure](measures.md) either.

`generate_network` works for every type. For a gas it returns the closure of
its rule if the gas has a generate face, and otherwise the network observed in
one run of its process.

## Parameters

Every parameter is optional. Pass them as keywords:

```python
net = chemart.generate_network("matrix-chemistry", seed=0, N=4)
```

To find out what a chemistry takes, ask it — the schema comes from the catalog,
so it is always current:

```python
d = chemart.describe_chemistry("matrix-chemistry")
d["params"]["properties"]["N"]
# {'type': 'integer', 'default': 4, 'minimum': 2,
#  'description': 'string length; also fixes |S| = 2^N - 1 (range: ...)'}
```

Each parameter's description carries its `meaning` and, where the paper used
larger values, a `range:` note. That note matters:

!!! tip "Defaults are small on purpose"
    Defaults are sized so a zero-argument call finishes in seconds — not to
    reproduce the paper. Before concluding a chemistry "doesn't do" something,
    check the `range` note and try the published setting. Several entries state
    outright which phenomena their defaults cannot reach.

## Seeds and reproducibility

`seed` is an argument of `generate_network`, never a catalog parameter. The same
seed always produces the same network:

```python
a = chemart.generate_network("gard", seed=7)
b = chemart.generate_network("gard", seed=7)
assert a.to_dict() == b.to_dict()
```

Omitting the seed draws from entropy, so a bare `chemart generate <id>` on the
command line will differ run to run for stochastic chemistries. Pass `--seed`
when you want to compare.

## Errors name what went wrong

Because the catalog knows each parameter's type, bounds and meaning, errors are
actionable rather than tracebacks:

```python
chemart.generate_network("brusselator", nonsense=1)
# ValueError: brusselator: unknown parameter 'nonsense'.
#   Valid parameters: a, b, k1, k2, k3, k4.

chemart.describe_chemistry("brusselater")
# ValueError: unknown chemistry 'brusselater'. Did you mean: brusselator?
```

## The command line

The CLI mirrors the API exactly, and is usually the fastest way to look at
something:

```bash
uv run chemart list                     # every chemistry, JSON
uv run chemart describe gamma           # one chemistry, JSON
uv run chemart generate brusselator --format summary
uv run chemart generate brusselator --format text
uv run chemart generate brusselator --format json
uv run chemart generate matrix-chemistry -p N=4 --seed 0 --format summary
```

`-p NAME=VALUE` may be repeated. Values are parsed as JSON when possible, so
lists and numbers work naturally and bare words are taken as strings:

```bash
uv run chemart generate kauffman-autocatalytic-sets -p N=6 -p p=0.01
uv run chemart generate gamma -p 'rules=["x, y -> x + y"]'
```

### Formats

| `--format` | gives you |
|---|---|
| `summary` | counts, status, capability tags, extras keys — the default |
| `text` | one reaction per line, with rate and firing count |
| `json` | the complete record |

## Looking at what you got

```python
net.summary()
# brusselator: 6 species, 4 reactions, status=complete
# provides: catalysts, initial-state, rate-constants, stoichiometry, topology
# seed: 1
# extras: buffered

net.to_text()
# A -> X  [mass-action k=1.0]
# B + X -> Y + D  [mass-action k=1.0]
# 2 X + Y -> 3 X  [mass-action k=1.0]
# X -> E  [mass-action k=1.0]
```

Always read `status` before interpreting the reaction list — see
[How it works](../concepts.md#status-changes-what-the-network-means).

## Finding the right chemistry

If you know the name, use it. If you know the *behaviour*, these are reliable
starting points:

| you want | try |
|---|---|
| oscillation, pattern formation | `brusselator`, `oregonator`, `repressilator` |
| enzyme kinetics | `michaelis-menten`, `hill-kinetics` |
| evolution, selection | `quasispecies`, `replicator-equation`, `ecolab` |
| autocatalysis, origin of life | `kauffman-autocatalytic-sets`, `raf`, `bagley-farmer`, `chemoton`, `gard` |
| constructive / open-ended | `matrix-chemistry`, `alchemy`, `combinator-chemistry`, `prime-number-chemistry` |
| digital organisms | `tierra`, `avida`, `corewar`, `coreworld` |
| rewriting formalisms | `gamma`, `p-systems`, `kappa-calculus`, `mgs` |
| spatial, agent-based | `squirm3`, `swarm-chemistry`, `sr-loops` |
| wet chemistry, in silico | `dna-hpp`, `dna-automaton`, `self-propelled-droplets` |

For anything more precise — "which chemistries give me both energies and a
conservation law?" — filter on capability rather than reputation:

```python
from chemart.catalog import load

for c in load():
    if {"energies", "mass-conservation"} <= set(c.provides):
        print(c.id, c.fidelity)
```

```bash
uv run python -m chemart.catalog query --provides rate-constants
uv run python -m chemart.catalog query --family origin-of-life
uv run python -m chemart.catalog show gard
uv run python -m chemart.catalog status
```

Or browse the [Catalog](../catalog/index.md), which has a page per chemistry
generated from the library itself.
