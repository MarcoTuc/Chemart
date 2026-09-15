# Chemart — a mart of artificial chemistries

Chemart aims to collect every artificial chemistry (AC) model in the
literature, expose them behind one library that *generates* them, and
provide interfacing tools that plug the generated chemical reaction
networks into the software people actually use to analyse CRNs.

## Status

Phase 1 (**collection**) is done for the field's reference survey:
Wolfgang Banzhaf & Lidia Yamamoto, *Artificial Chemistries*
(MIT Press, 2015).

**98 chemistries** are catalogued with their algorithm, their
hyperparameters (meaning, role, defaults and ranges), and what each one
can actually hand back — topology only, or kinetics, or energetics.

- [`docs/CATALOG.md`](docs/CATALOG.md) — generated index of all 98
- [`docs/OUTPUT-SCHEMA.md`](docs/OUTPUT-SCHEMA.md) — why
  *(stoichiometric matrix + species names)* is not enough, and what the
  standardised output record should be instead
- [`catalog/SCHEMA.md`](catalog/SCHEMA.md) — field definitions
- [`catalog/chemistries/`](catalog/chemistries/) — the catalog itself,
  one YAML file per chemistry

At a glance:

| | |
|---|---|
| chemistries catalogued | 98 |
| fully specified in the book (implementable as-is) | 53 |
| constructive (open, growing species set) | 49 |
| carry their own rate constants or rate law | 38 |
| carry energetics / thermodynamic consistency | 15 |
| define compartments | 16 |
| define space | 24 |

## Using the catalog

```bash
python -m chemart.catalog validate            # check invariants
python -m chemart.catalog index               # regenerate docs/CATALOG.md
python -m chemart.catalog show matrix-chemistry
python -m chemart.catalog query --provides rate-constants --ready yes
python -m chemart.catalog query --family origin-of-life
```

```python
from chemart.catalog import load

for c in load():
    if c.constructive and c.generator_ready == "yes":
        print(c.id, [p.name for p in c.params_by_role("structural")])
```

## Using the library

Environment is managed with [uv](https://docs.astral.sh/uv/): `uv sync`,
then run everything through `uv run`.

```python
import chemart

chemart.list_chemistries()                     # all 98, with implementation status
chemart.describe_chemistry("matrix-chemistry") # metadata + JSON Schema of the parameters
net = chemart.generate_network("matrix-chemistry", seed=0, N=4)

print(net.summary())
print(net.to_text())            # one reaction per line
species, R, P = net.matrices()  # sparse reactant/product matrices; S = P - R
net.to_dict()                   # plain JSON data; Network.from_dict inverts it
```

```bash
uv run chemart list
uv run chemart describe matrix-chemistry
uv run chemart generate matrix-chemistry -p N=4 --seed 0 --format json
```

For LLM agents, `chemart.tool_definitions()` returns the three functions as
JSON-Schema tool specs and `chemart.call_tool(name, arguments)` executes them.

## Next

Implementation proceeds in waves; see [`docs/PLAN.md`](docs/PLAN.md) and
`uv run python -m chemart.catalog status` for progress.

## Source

All page and section references are to Banzhaf & Yamamoto (2015).
Bracketed numbers in `refs:` are that book's bibliography entries.
