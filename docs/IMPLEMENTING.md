# Implementing one chemistry

The brief for whoever (person or agent) implements a single catalogued
chemistry. Follow it exactly: many chemistries are implemented in parallel
in the same working tree, and the library only stays coherent if everyone
respects the same contract and touches only their own files.

## 1. What you deliver

For chemistry id `<id>` (module name `<mod>` = id with `-` → `_`), you own
exactly three files:

| file | content |
|---|---|
| `chemart/chemistries/<mod>.py` | `def generate(p, rng) -> Network` |
| `catalog/chemistries/<id>.yaml` | the catalog entry, migrated to schema v2 |
| `tests/chemistries/test_<mod>.py` | tests reproducing published results |

**Everything else is read-only for you.** Do not edit `chemart/*.py`,
`chemart/helpers/*`, `tests/chemistries/odes.py`, `tests/test_*.py`, other
chemistries, `pyproject.toml`, `uv.lock` or `docs/`. Do not run `uv add`,
`uv sync` or `python -m chemart.catalog index`. Do not use git. If you
believe a shared file needs a change (a new helper, a new rate law, a
dependency, a core bug), **do not make it**: describe it in your report,
and the integrator will apply it.

Run everything with `uv run …` from the repository root (it uses the
project's `./.venv`). Never use pip or a system python.

## 2. The interface you implement

```python
def generate(p, rng) -> Network:
    ...
```

- `p` holds the validated parameters as attributes (`p.N`); names, types
  and defaults come **only** from your catalog entry. `rng` is a
  `numpy.random.Generator` seeded by the caller: use it for all randomness
  (never `random` or `np.random.*` globals), so the same seed gives the same
  network.
- Return a `chemart.network.Network`. `generate_network` fills in
  `chemistry`, `params` and `seed`; you don't.
- Build explicit reactions with `chemart.helpers.explicit.network`
  (`"2 X + Y -> 3 X"`, with spaces around `+` and a space after
  coefficients). Reuse `chemart.helpers.params` for list/dict parameter checks,
  `chemart.expand.expand` for closures of constructive chemistries
  (pass `alternatives=True` when one set of reactants can react in several
  ways, returning a list of outcomes), and
  `chemart.soup.soup` for well-stirred multiset runs that record the reactions
  that fired.
- Raise `ValueError` with an actionable message for invalid combinations
  of parameters that the catalog cannot express (e.g. `K < N`).

## 3. Network conventions

- **Species** ids are non-empty strings without spaces. If a species has
  internal structure (bitstring, λ-term, program, graph), put it in
  `Species(id, structure=...)` as a string.
- **Stoichiometry** is positive integers. Catalysts appear on both sides.
- **Rates** are `None` (topology only) or a dict whose `law` is in
  `chemart.kinetics.RATE_LAWS`: `mass-action`, `power`, `michaelis-menten`,
  `hill`, `saturating`, `arrhenius`. Read the definitions in
  `chemart/kinetics.py`. Never invent a law. Extra facts that are not a law
  (a threshold gate, a crowding capacity) go in the dict as extra scalar
  keys. Never put a made-up number in a rate: if the source gives no rate,
  use `None`.
- **status**: `complete` (the whole defined network or a finished closure),
  `truncated` (a closure cut off by a size parameter), or `observed` (reactions
  that fired in a simulation, each carrying `count`).
- **Flow**: the book's flow reactor with non-selective dilution is
  `outflow="constant-total"`. `inflow` maps species to a constant influx
  (amount per volume per time); a numeric `outflow` is a first-order
  removal rate, given per species (dict) or for all species (number).
  Species held at a constant concentration are listed in `extras["buffered"]`.
  `tests/chemistries/odes.py` applies all of these, so don't add flow terms
  in your tests yourself.
- **Facts without a formula** (an inhibitor with no published rate law, a
  threshold that gates an outflow) are recorded as extra scalar keys on the
  rate dict or in `extras`. They are data only: nothing applies them, and
  you must not invent a formula for them.
- **extras** reserved keys: `space`, `compartments`, `energies`,
  `conservation` (list of `{"name", "vector", optional "modulus"}`),
  `analysis`, `interaction_law`. Other keys are free-form but must be JSON.
- Everything must be JSON data: plain `int`/`float`/`str`/`bool`/`list`/
  `dict` with string keys (convert numpy scalars with `int()`/`float()`).
- Non-CRN mechanisms (agents, CA, VMs, force laws) are exposed through the
  reaction events they produce: typically a `soup`-style run returning the
  observed reactions with counts, or species plus `extras["interaction_law"]`
  when there are genuinely no reactions.

## 4. The catalog entry (schema v2)

Your file `catalog/chemistries/<id>.yaml` currently holds the v1 entry.
Migrate it in place, keeping `id`, `name`, `family`, `kind`, `S`, `R`, `A`
and `refs`, and correcting anything the sources contradict:

- `params`: every parameter has `type` in `int | float | bool | str | enum |
  list | dict`, a `default` that satisfies the spec, `min`/`max` for
  numbers where meaningful, `choices` for enums, a `meaning`, and a `role`.
  No `callable`, `matrix` or `seed` types. Callables become named choices;
  rule sets become text or JSON. Drop parameters that only size a simulation
  and do not change the network (and say so in `decisions`).
- **Small defaults**: `generate_network("<id>")` with no arguments must
  finish in well under 5 s and give a network that is readable (ideally
  ≤ ~20k reactions). Paper-scale values go in `range`.
- Remove `generator_ready`. Add `fidelity`:
  - `book`: implemented exactly as the book specifies;
  - `book+decisions`: the book leaves gaps you filled; list each in `decisions`;
  - `reconstructed`: built from the original papers; list them in `sources`
    (full citation + URL) and list remaining gaps in `decisions`.
- `provides` must include every capability tag your default network
  computes (`Network.provides`); the contract test checks this.
- Fix wrong facts in `phenomena`/`notes`.

See `catalog/SCHEMA.md` for the field definitions and any implemented
chemistry (e.g. `catalog/chemistries/chemoton.yaml`) for the style.

## 5. The spec: book first, papers when the book is not enough

1. Read your section(s) of Banzhaf & Yamamoto, *Artificial Chemistries*
   (2015). Extracted text: `/home/marco/.claude/jobs/a14e4701/tmp/book.txt`
   (use the line ranges in your assignment; bibliography entries `[n]` are
   near the end of that file).
2. If the book does not specify the algorithm, the parameters or the
   numbers you need, **reconstruct from the original papers** cited in
   `refs`: find open versions (WebSearch/WebFetch, arXiv, authors' pages,
   the Wayback Machine), download PDFs with `curl` and extract them with
   `pdftotext -layout` into your own scratch folder
   `/home/marco/.claude/jobs/a14e4701/tmp/agents/<id>/`.
3. **Verify every number you use against the text.** When a table and the
   prose disagree, trust the data and record the discrepancy in
   `decisions`. When the extracted text is garbled (lost primes, square
   roots, subscripts), resolve it from context and say how.
4. If a paper is inaccessible, do not stall: implement from the book, use
   `fidelity: book+decisions`, and record the missing source in `decisions`.
   Never present an invented value as published.

## 6. Tests

`tests/chemistries/test_<mod>.py` must contain at least one test that
reproduces something **published**: a table, a count, a closure, the rate
equations of the paper, a steady state, a qualitative dynamic behaviour.
If nothing numeric is published, turn the catalog's `phenomena` into
property tests. `tests/chemistries/odes.py` integrates mass-action (and
power, Michaelis-Menten, Hill, saturating) networks and supports
`constant-total` and `buffered`; import it with `from odes import integrate, rhs`.
Keep the file under ~10 s.

## 7. Definition of done

All three must pass, run from the repository root:

```bash
uv run python -m chemart.catalog validate --only <id>
CHEMART_ONLY=<id> uv run pytest -q tests/test_contract.py tests/chemistries/test_<mod>.py
uv run chemart generate <id> --format summary
```

Do not run the whole test suite: other chemistries are being written at the
same time. If a failure mentions another chemistry's files, wait a moment and
rerun.

## 8. Report

End with a short report:

- **Files**: the three paths.
- **Fidelity** and **sources** (citations with URLs actually used).
- **Decisions**: one line each.
- **Tests**: what each test reproduces and from where (book section, table,
  equation, paper figure).
- **Definition of done**: the output lines of the three commands.
- **Requests for shared changes**, if any, with the exact change proposed.
- **Open doubts** you could not resolve.

## 9. Examples to imitate

| pattern | module |
|---|---|
| written-down network with published dynamics | `brusselator.py`, `repressilator.py` |
| reconstructed from a paper, equation-level test | `chemoton.py`, `flow_ac.py`, `gard.py` |
| sampled topology with a model-level algorithm | `bigan_conservative_crn.py`, `jain_krishna.py` |
| structure → rates | `farmer_immune.py`, `quasispecies.py` |
| constructive closure (`expand`) | see `chemart/expand.py` and its tests in `tests/test_core.py` |
