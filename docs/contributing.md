# Contributing a chemistry

The contract is deliberately small: one function, one catalog entry, one test
file. No base classes, no registration step.

## The three files you own

For chemistry id `<id>` (module name `<mod>` = id with `-` → `_`):

| file | content |
|---|---|
| `chemart/chemistries/<mod>.py` | `def generate(p, rng) -> Network` |
| `catalog/chemistries/<id>.yaml` | the catalog entry |
| `tests/chemistries/test_<mod>.py` | tests reproducing published results |

Registration is by naming convention — `matrix-chemistry` resolves to
`chemart.chemistries.matrix_chemistry`, imported lazily. Nothing else needs
editing.

## The generate contract

```python
def generate(p, rng) -> Network:
    ...
```

- `p` holds the validated parameters as attributes (`p.N`). Names, types and
  defaults come **only** from the catalog entry — there is no second source of
  truth in Python.
- `rng` is a seeded `numpy.random.Generator`. Use it for all randomness; never
  `random` or `np.random.*` globals, or the same seed will stop reproducing the
  same network. The contract test checks this for every entry.
- Return a `Network`. `generate_network` fills in `chemistry`, `params` and
  `seed`; you do not.
- Raise `ValueError` with an actionable message for parameter combinations the
  catalog cannot express (`K < N`, an empty rule set).

Everything in the returned network must be plain JSON: convert numpy scalars
with `int()` / `float()`, and keep dict keys strings.

## The catalog entry

One YAML file, `chemistries: [ <entry> ]`. Field definitions are in
`catalog/SCHEMA.md`. The rules that bite:

- **Parameter types are JSON only** — `int | float | bool | str | enum | list |
  dict`. Callables become named choices (`distribution: enum [uniform,
  lognormal]`); rule sets become text in a documented mini-language, or JSON.
- **Every parameter needs a `default`, a `meaning` and a `role`.** Roles are
  `structural | kinetic | thermodynamic | population | spatial | stochastic |
  selection`, and they are what lets the library scale knobs coherently across
  chemistries.
- **No `seed` parameter** — it is an argument of `generate_network`.
- **Small defaults.** A zero-argument call must finish in well under five
  seconds and give a readable network. Paper-scale values go in the parameter's
  `range` note, which is prose and shows up in `describe_chemistry`.
- **`fidelity` is required**, with `sources` for `reconstructed` and `decisions`
  for `book+decisions`.
- **`provides` must cover what your network computes** — the contract enforces
  computed ⊆ claimed.

## Sources, and what to do when they run out

Read the book section first; go to the original papers when the book is thin.
If a paper cannot be obtained, **do not stall and do not invent**: implement
from the book, set `fidelity: book+decisions`, and record the missing source in
`decisions`. An honest reconstruction that names its gaps is far more useful
than a fabricated number.

!!! tip "Verify against the bibliography — including the brief you were given"
    Over the full build, this caught errors in the catalog, in the book, in a
    published paper, and in six task briefs — one of which sent an implementer
    looking for a source paper for an algorithm the book proposes itself. When
    the bibliography disagrees with your instructions, follow the bibliography
    and say so.

## Rates

Use a law from `chemart.kinetics.RATE_LAWS`, or `None`. Never invent a law or a
constant. Facts that are not a propensity — a threshold, a temperature, a
crowding capacity — go as extra scalar keys on the rate dict or in `extras`;
they are data, and nothing applies them automatically.

## Non-CRN mechanisms

Agents, cellular automata, virtual machines and force laws are exposed through
the reaction events they produce: a `soup`-style run returning observed
reactions with counts, or species plus `extras["interaction_law"]` when there
genuinely are no reactions.

## Tests

At least one test must reproduce something **published** — a table, a count, a
closure, a steady state, a figure's qualitative claim. If nothing numeric
exists, turn the catalog's `phenomena` into property tests.

Mark anything over about two seconds with `@pytest.mark.slow`; slow tests are
excluded by default and run with `-m "slow or not slow"`.

```bash
CHEMART_ONLY=<id> uv run pytest -q -m "slow or not slow" \
    tests/test_contract.py tests/chemistries/test_<mod>.py
```

The contract tests every implemented entry for: zero-argument generation under
five seconds, an exact JSON round trip, same seed → same network, computed
capabilities ⊆ catalog claims, and that unknown or out-of-range parameters raise
`ValueError` naming the parameter.

## Definition of done

```bash
uv run python -m chemart.catalog validate --only <id>      # 0 problems
CHEMART_ONLY=<id> uv run pytest -q -m "slow or not slow" \
    tests/test_contract.py tests/chemistries/test_<mod>.py
uv run chemart generate <id> --format summary
```

Then write the chemistry's documentation page. Its prose lives in
`catalog/explainers/<id>.md`: what the chemistry is and why it exists, how it
works, how to use it, and what has been done with it. The generator weaves it
around the parts it builds from the catalog: the formal specification, the
default network, the parameters, the decisions and the sources.
`catalog/explainers/README.md` gives the format and the rules, and
`catalog/explainers/ccm.md` is a worked example.

```bash
uv run python tools/gen_catalog_pages.py --only <id>       # rebuild just this page
uv run pytest -q tests/test_explainers.py
```

Finally regenerate the index and all the docs pages:

```bash
uv run python -m chemart.catalog index
uv run python tools/gen_catalog_pages.py
```
