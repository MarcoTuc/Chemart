# Chemart

**A mart of artificial chemistries.** Every artificial chemistry in the field's
reference survey, available behind one interface — generate its reaction
network, simulate it, evolve it, measure it, compare it with any other.

The chemistries come from Wolfgang Banzhaf & Lidia Yamamoto, *Artificial
Chemistries* (MIT Press, 2015). Every one is catalogued **and implemented**,
with tests reproducing published results wherever any exist. Chemistries
published after the book are added from their papers, such as
[Combinatory Chemistry](catalog/combinatory-chemistry.md) (2020) and
[BFF](catalog/bff.md) (2024).

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
    and why unrelated models can share one output format.

-   **Want to use it?**

    [Generating networks](guide/generating.md) covers the functions, the CLI,
    parameters and seeds; then [simulate](guide/simulating.md) a network,
    [evolve](guide/evolving.md) a Turing gas, and [measure](guide/measures.md)
    either.

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
equations — and AlChemy, in which λ-terms react by applying one to the other.
Operationally they share nothing. But both can hand back a **chemical
reaction network**: which species exist, which reactions consume and produce
them, and the rates where the model defines any. Chemart's bet is that this
common output is enough to make incomparable models comparable.

How you get that network depends on the chemistry's type. A **given**
chemistry writes it down; a **generator** computes it from arguments; a
**gas** (a Turing gas) never has one in advance, because its reactions exist
only as its molecules meet. Chemart simulates the first two and evolves the
third, and measures all of them the same way. [How it works](concepts.md#three-types-of-chemistry)
has the details.

## Install and first run

The environment is managed with [uv](https://docs.astral.sh/uv/). Run everything
through `uv run` so you get the project's own `./.venv`.

```bash
uv sync
uv run chemart list                          # every chemistry, as JSON
uv run chemart describe brusselator          # metadata + parameter JSON Schema
uv run chemart generate brusselator --format text
uv run chemart simulate brusselator --t-end 40
uv run chemart evolve alchemy --seed 1 --track shannon
uv run chemart measure raf --seed 1
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

The [Catalog](catalog/index.md) lists every chemistry and counts what they
supply: kinetics, energetics, conservation laws, space, compartments. Its
counts are computed from the library each time the site is built.

Entries that are not chemistry, such as the artificial-life systems (Tierra,
Avida, Core War), are kept in the catalog's [archive](catalog/index.md#archive):
still runnable by id, but left out of listings and the LLM tools.

## For LLM agents

The interface ships as six JSON-Schema tool specs: list, describe, generate,
simulate, evolve and measure.

```python
chemart.tool_definitions()   # ready to register
chemart.call_tool("generate_network", {"chemistry": "brusselator", "seed": 1})
chemart.call_tool("measure_network", {"chemistry": "raf", "seed": 1})
```

The `chemistry` argument is an enum of every implemented id in the catalog,
so a model cannot ask for one that does not exist, and the simulate and
evolve tools cap their answers to fit a model's context. See
[API and CLI](reference/api.md).
