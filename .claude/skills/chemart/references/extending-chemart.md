# Adding or changing a chemistry

The contract is deliberately small: one or two functions, one catalog entry,
one test file. There are no base classes and no registration step.

## Contents
- [The three files you own](#the-three-files-you-own)
- [The generate contract](#the-generate-contract)
- [The evolve contract](#the-evolve-contract)
- [The catalog entry](#the-catalog-entry)
- [Helpers](#helpers)
- [Tests](#tests)
- [Definition of done](#definition-of-done)
- [Working in parallel](#working-in-parallel)

## The three files you own

For chemistry id `<id>` (module name `<mod>` = id with `-` → `_`):

| file | content |
|---|---|
| `chemart/chemistries/<mod>.py` | `def generate(p, rng) -> Network`, and/or `def evolve(p, rng)` |
| `catalog/chemistries/<id>.yaml` | the catalog entry |
| `tests/chemistries/test_<mod>.py` | tests reproducing published results |

Registration is by naming convention — id `matrix-chemistry` resolves to
`chemart.chemistries.matrix_chemistry`, imported lazily. Nothing else needs
editing. The repository's own brief for implementers is
`docs/IMPLEMENTING.md`; read it if you are contributing back.

## The generate contract

```python
def generate(p, rng) -> Network:
    ...
```

- `p` holds validated parameters as attributes (`p.N`). Names, types and
  defaults come **only** from the catalog entry — there is no second source of
  truth in Python.
- `rng` is a seeded `numpy.random.Generator`. Use it for all randomness; never
  `random` or `np.random.*` globals, or the same seed will stop reproducing the
  same network (the contract test checks this for every entry).
- Return a `Network`. `generate_network` fills in `chemistry`, `params` and
  `seed`; you do not.
- Raise `ValueError` with an actionable message for parameter combinations the
  catalog cannot express (`K < N`, an empty rule set).

Everything in the returned network must be plain JSON: convert numpy scalars
with `int()` / `float()`, and keep dict keys strings.

## The evolve contract

A chemistry whose process is worth running — a Turing gas, a lattice — defines
an evolve face, besides `generate` or instead of it:

```python
from chemart.soup import Tally
from chemart.trajectory import Frame

def evolve(p, rng):
    tally = Tally()
    yield Frame(t=0, state=...)                                  # the initial population
    ...
    yield Frame(t=..., state=..., fired=tally.flush(), observables={...})
    return network                                               # what fired, with counts
```

- The first frame is the initial state and fires nothing; `t` never decreases
  and is counted in the entry's `clock` (required with an evolve face).
- `fired` lists the reactions since the previous frame as
  `[[reactants], [products], count]`; summed over frames, the counts must equal
  each reaction's `count` in the returned network (`Tally` keeps them in step).
- `observables` holds what only this chemistry can report (a fold's energy, a
  membrane count). Generic numbers such as richness come from
  `chemart.measures`; do not repeat them.
- For a well-stirred soup, `chemart.soup.stir(react, population, steps, rng,
  every=..., tally=...)` runs the loop and yields `(step, population, tally)`
  at each frame.
- A parameter only one face uses gets `face: evolve` (or `face: generate`) in
  the catalog. `every` is reserved: `chemart.evolve` uses it to merge frames.

## The catalog entry

One YAML file, `chemistries: [ <entry> ]`. Full field reference:
`catalog/SCHEMA.md`. The rules that bite:

- **Parameter types are JSON only**: `int | float | bool | str | enum | list |
  dict`. Callables become named choices (`distribution: enum [uniform,
  lognormal]`); rule sets become text in a documented mini-language, or JSON.
- **Every parameter needs a `default`, a `meaning` and a `role`.** Roles are
  `structural | kinetic | thermodynamic | population | spatial | stochastic |
  selection`, and they are what lets the library scale knobs coherently across
  chemistries.
- **`type` is required**: `given` (a network written down), `generator` (an
  algorithm computes it) or `gas` (a Turing gas), with the reason recorded in
  `catalog/TYPES.md`. `clock` is required when the module has `evolve`.
- **No `seed` parameter** — it is an argument of `generate_network`.
- **Small defaults.** A zero-argument call must finish in well under 5 seconds
  and give a readable network. Paper-scale values go in the parameter's `range`
  note, which is prose and shows up in `describe_chemistry`.
- **`fidelity` is required**, with `sources` for `reconstructed` and
  `decisions` for `book+decisions`.
- **`provides` must cover what your network computes.** The contract test
  enforces computed ⊆ claimed.

## Helpers

```python
from chemart.helpers.explicit import network      # written-down reactions
from chemart.helpers import params                # list/dict validation
from chemart.expand import expand                 # closure of a constructive rule
from chemart.soup import soup, stir, Tally        # well-stirred run; as frames; fired counts
```

`chemart.helpers.params` has `vector`, `square_matrix`, `edges` and
`apportion` for checking the inner shape of list/dict parameters, which the
catalog schema cannot express. They raise messages that name the parameter.

For `expand` and `soup` signatures and the `alternatives` flag, see
`comparing-and-analysing.md`.

**Rates**: use a law from `chemart.kinetics.RATE_LAWS` or `None`. Never invent
a law or a constant. Facts that are not a propensity — a threshold, a
temperature, a crowding capacity — go as extra scalar keys on the rate dict or
in `extras`; they are data, and nothing applies them automatically.

**Non-CRN mechanisms** (agents, CA, VMs, force laws) are exposed through the
reaction events they produce: an evolve face whose frames record what fired,
returning the observed reactions with counts, or species plus
`extras["interaction_law"]` when there genuinely are no reactions.

## Tests

At least one test must reproduce something **published** — a table, a count, a
closure, a steady state, a published figure's qualitative claim. If nothing
numeric exists, turn the catalog's `phenomena` into property tests.

Mark anything over ~2 seconds with `@pytest.mark.slow`; slow tests are excluded
by default and run with `-m "slow or not slow"`.

```bash
CHEMART_ONLY=<id> uv run pytest -q -m "slow or not slow" \
    tests/test_contract.py tests/chemistries/test_<mod>.py
```

`CHEMART_ONLY` narrows the parametrized contract tests to your chemistry so you
are not blocked by, or blamed for, someone else's work in progress. The
coverage gate deliberately ignores it.

The contract tests every implemented entry for: zero-argument generation under
5 seconds, exact JSON round trip, same seed → same network, computed
capabilities ⊆ catalog claims, and that unknown or out-of-range parameters
raise `ValueError` naming the parameter. A `given` entry must keep its topology
across seeds. With an evolve face it also checks the frames: the clock, frame
0, time order, JSON, same seed → same frames, and fired counts that add up to
the network's.

## Definition of done

```bash
uv run python -m chemart.catalog validate --only <id>       # 0 problems
CHEMART_ONLY=<id> uv run pytest -q -m "slow or not slow" \
    tests/test_contract.py tests/chemistries/test_<mod>.py
uv run chemart generate <id> --format summary
```

Then regenerate the index: `uv run python -m chemart.catalog index`.

## Working in parallel

The library was built by many agents implementing one chemistry each in a
shared tree. What made it work, if you repeat it:

- **File ownership.** Each implementer touches only their three files.
  Everything else — core modules, shared tests, docs, `pyproject.toml` — is
  read-only, and requested changes go to an integrator who applies them once,
  between batches, never while others are running.
- **Dependency changes happen between batches**, never with agents mid-flight.
- **Verify sources against the bibliography, including the brief.** Over the
  full build this caught errors in the catalog, in the book, in a published
  paper, and in six task briefs — one of which sent an implementer looking for
  a source paper for an algorithm the book proposes itself. Tell implementers
  to follow the bibliography and say so when it disagrees with their
  instructions.
