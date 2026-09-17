# Chemart — a mart of artificial chemistries

Chemart aims to collect every artificial chemistry (AC) model in the
literature, expose them behind one library that *generates* them, and
provide interfacing tools that plug the generated chemical reaction
networks into the software people actually use to analyse CRNs.

## Status

**All 98 chemistries are catalogued and implemented.** The source is the
field's reference survey: Wolfgang Banzhaf & Lidia Yamamoto, *Artificial
Chemistries* (MIT Press, 2015).

Every one is reachable through the same three functions, generates in
seconds from zero arguments, round-trips through JSON, and is covered by
tests that reproduce published results wherever any exist.

| fidelity | count | meaning |
|---|---|---|
| `book` | 4 | implemented exactly as the book specifies |
| `book+decisions` | 33 | the book left gaps; each filled choice is listed in `decisions` |
| `reconstructed` | 61 | built from the original papers, cited in `sources` |

Where a source could not be obtained, the entry says so in its own
`decisions` rather than inventing a value, and
[`to_decide.md`](to_decide.md) collects those cases in one place — the
low-confidence reconstructions, and the entries implementing only part of
their scope because half the specification is paywalled.

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
| chemistries catalogued and implemented | 98 |
| constructive (open, growing species set) | 51 |
| carry their own rate constants or rate law | 44 |
| carry energetics / thermodynamic consistency | 11 |
| define compartments | 10 |
| define space | 21 |

By kind: 80 generators, 7 formalisms, 4 analyses, 4 wet, 3 frameworks.
These counts are what the generated networks actually compute — several
were revised down during implementation, when a chemistry turned out not
to supply what the phase-1 catalog claimed for it.

## Using the catalog

```bash
python -m chemart.catalog validate            # check invariants
python -m chemart.catalog index               # regenerate docs/CATALOG.md
python -m chemart.catalog show matrix-chemistry
python -m chemart.catalog query --provides rate-constants
python -m chemart.catalog query --family origin-of-life
python -m chemart.catalog status              # implementation and fidelity counts
```

```python
from chemart.catalog import load

for c in load():
    if c.constructive and "rate-constants" in c.provides:
        print(c.id, c.fidelity, [p.name for p in c.params_by_role("structural")])
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

## Open questions

The mart is complete, so what is left is judgement rather than
implementation. [`to_decide.md`](to_decide.md) holds it:

- the generative vs algorithmic classification, deliberately parked;
- **licensing** — the Core War tests vendor two GPL-2 Redcode warriors
  from pMARS, and Chemart itself declares no license yet;
- two record conventions the genome-carrying chemistries disagree on
  (how much genome to store, and which non-replication events are
  reactions);
- a proposed catalysed Michaelis-Menten rate law, and the `arrhenius`
  gap in the shared ODE helper;
- the entries whose primary sources could not be reached, and those
  implementing only part of their scope.

[`docs/PLAN.md`](docs/PLAN.md) records how it was built.

## Source

All page and section references are to Banzhaf & Yamamoto (2015).
Bracketed numbers in `refs:` are that book's bibliography entries.
