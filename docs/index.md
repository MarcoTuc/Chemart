# Chemart

**A mart of artificial chemistries.** Every artificial chemistry in the field's
reference survey, available behind one interface — generate its reaction
network, inspect it, simulate it, compare it with any other.

All **98** chemistries from Wolfgang Banzhaf & Lidia Yamamoto,
*Artificial Chemistries* (MIT Press, 2015) are catalogued **and implemented**,
each with tests reproducing published results wherever any exist.

```python
import chemart

net = chemart.generate_network("brusselator", seed=1)
print(net.to_text())
# A -> X  [mass-action k=1.0]
# B + X -> Y + D  [mass-action k=1.0]
# 2 X + Y -> 3 X  [mass-action k=1.0]
# X -> E  [mass-action k=1.0]
```

<div class="grid cards" markdown>

-   **New here?**

    Start with [How it works](concepts.md) — what an artificial chemistry is,
    and why 98 unrelated models can share one output format.

-   **Want to use it?**

    [Generating networks](guide/generating.md) covers the three functions, the
    CLI, parameters and seeds.

-   **Looking for a specific chemistry?**

    Browse the [Catalog](catalog/index.md) — one page per chemistry, generated
    from the library itself.

-   **About to quote a number?**

    Read [Fidelity and trust](trust.md) first. Some entries rest on sources
    that could not be obtained, and they say so.

</div>

## The idea in one paragraph

An artificial chemistry is a triple **(S, R, A)**: a set of possible molecules,
the rules by which they react, and the algorithm that runs them. That definition
is broad enough to cover both the Brusselator — three coupled differential
equations — and Tierra, in which self-replicating assembly programs compete for
CPU time. Operationally they share nothing. But both can hand back a **chemical
reaction network**: which species exist, which reactions consume and produce
them, and the rates where the model defines any. Chemart's bet is that this
common output is enough to make incomparable models comparable.

## Install and first run

The environment is managed with [uv](https://docs.astral.sh/uv/). Run everything
through `uv run` so you get the project's own `./.venv`.

```bash
uv sync
uv run chemart list                          # every chemistry, as JSON
uv run chemart describe brusselator          # metadata + parameter JSON Schema
uv run chemart generate brusselator --format text
```

Every parameter is optional, and defaults are sized so a zero-argument call
finishes in seconds:

```python
import chemart

chemart.list_chemistries()                   # id, name, summary, fidelity
chemart.describe_chemistry("gard")           # everything known about one
net = chemart.generate_network("gard", seed=0)

# The same seed always gives the same network.
assert chemart.generate_network("gard", seed=0).to_dict() == net.to_dict()
```

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

## For LLM agents

The whole interface is three functions, and they ship as JSON-Schema tool specs:

```python
chemart.tool_definitions()   # list/describe/generate, ready to register
chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
```

The `generate_network` tool's `chemistry` argument is an enum of every
implemented id, so a model cannot ask for one that does not exist. See
[API and CLI](reference/api.md).
