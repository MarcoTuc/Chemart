# Finding the right chemistry

This is how to narrow the catalog down without reading every entry. Counts are
never written down here: compute them (below), or read `docs/CATALOG.md`,
which is generated from the catalog.

## Contents
- [The four axes](#the-four-axes)
- [The archive](#the-archive)
- [Querying from Python](#querying-from-python)
- [Querying from the CLI](#querying-from-the-cli)
- [What each capability tag means](#what-each-capability-tag-means)
- [Families](#families)
- [Reading a describe_chemistry result](#reading-a-describe_chemistry-result)

## The four axes

**`provides`** — what the generator can actually hand you. This is the axis
that matters most for analysis: if you need rate constants to simulate, filter
on `rate-constants`, not on the chemistry's reputation.

**`kind`** — what the thing *is*: `generator` (most entries), `formalism`
(e.g. Gamma and P systems: a rule language, where you supply the rules),
`analysis` (e.g. RAF: builds a model *and* runs an algorithm on it), `wet` (the
in-silico network of a published laboratory experiment), `framework` (larger
simulation platforms, minimally re-implemented).

**`family`** — subject area (see below).

**`constructive`** — true means the species set is open and grows at run time,
so the interesting object is a *closure* and the network may come back
`truncated`. False means the species set is fixed and enumerable.

## The archive

Entries whose YAML has an `archived:` field are not part of the chemistry
catalog: `archived: artificial-life` for the artificial-life systems (Tierra,
Avida, Core War, Coreworld, Swarm Chemistry), `archived: pruned` for entries
set aside on review. `list_chemistries()`, `chemart list` and the LLM tools
leave them out; `generate_network(id)` still runs them. Don't offer an archived
entry as "a chemistry" unless the user asks for it by name or asks for the
archive (`list_chemistries(include_archived=True)`, `chemart list --all`).

## Querying from Python

```python
from collections import Counter
from chemart.catalog import active, load

catalog = active(load())        # the chemistry catalog; load() alone includes the archive

# everything that supplies kinetics AND a conservation law
for c in catalog:
    if "rate-constants" in c.provides and "mass-conservation" in c.provides:
        print(c.id, c.fidelity)

# constructive chemistries, with whatever knob sizes each one
for c in catalog:
    if c.constructive:
        knobs = [p.name for p in c.params_by_role("structural")]
        print(f"{c.id:28s} {knobs}")

# counts, computed rather than remembered
Counter(tag for c in catalog for tag in c.provides)
Counter(c.family for c in catalog)
```

`load()` parses the YAML once and caches it; each call gets its own list, so it
is cheap to call repeatedly.

Useful `Chemistry` fields: `id`, `name`, `family`, `archived`, `kind`,
`constructive`, `provides`, `fidelity`, `sources`, `decisions`, `phenomena`,
`book`, `refs`, `params`, and the properties `implemented`, `module`, `tiers`,
`params_by_role(role)`.

Parameter objects carry `name`, `type`, `default`, `min`, `max`, `choices`,
`meaning`, `role`, `range`.

## Querying from the CLI

```bash
uv run python -m chemart.catalog query --provides rate-constants
uv run python -m chemart.catalog query --family origin-of-life   # add --all for the archive
uv run python -m chemart.catalog show gard          # human-readable entry
uv run python -m chemart.catalog status             # implementation + fidelity counts
uv run python -m chemart.catalog validate           # invariants; should print 0 problems
```

Or use the bundled `scripts/survey.py`, which does the common intersections
and prints a table.

## What each capability tag means

| tag | meaning |
|---|---|
| `topology` | species + who reacts with whom (every entry) |
| `stoichiometry` | separate reactant/product multiplicities (every entry) |
| `catalysts` | some species appears on both sides and must be kept |
| `initial-state` | the chemistry prescribes a starting multiset |
| `sequence-structure-function` | molecules carry structure the rule reads |
| `rate-constants` | the chemistry itself prescribes k's |
| `mass-conservation` | an atom/mass vector m with Sᵀm = 0 |
| `flow` | inflow/outflow is part of the definition |
| `space` | positions, lattice or diffusion |
| `compartments` | nested membranes/cells |
| `energies` | per-species free or bond energies |
| `rate-law` | a non-mass-action propensity is part of the model |
| `thermodynamic-consistency` | reverse rates constrained by ΔG |

The three-tag rule of thumb: `topology` alone means you can study structure;
add `rate-constants` and you can simulate; add `energies` or
`thermodynamic-consistency` and you can do thermodynamics.

## Families

`core`, `rewriting`, `automata`, `bio-inspired`, `origin-of-life`,
`evolutionary-dynamics`, `network`, `spatial`, `application`,
`systems-biology`, `wet`, `non-chemical`. List a family's members with
`python -m chemart.catalog query --family <family>`, or see them grouped in
`docs/CATALOG.md`.

Family is subject area, not capability, and it is the least reliable axis —
`oregonator` and `repressilator` sit in `wet` because they model real
chemistry, not because they need a laboratory. Filter on `provides` when the
question is technical. Some family assignments are still open; see
`to_decide.md` in the repository.

## Reading a describe_chemistry result

```python
d = chemart.describe_chemistry("gard")
d["molecules"]    # S: how a molecule is represented, explicit or implicit
d["reactions"]    # R: arity and the canonical reaction scheme
d["reactor"]      # A: reactor type and how population size is bounded
d["params"]       # JSON Schema: properties, defaults, min/max, enum, descriptions
d["phenomena"]    # what the model is known to produce — the informal spec
d["decisions"]    # where the sources were thin and what was chosen
d["sources"]      # citations actually used
d["book"], d["refs"]   # section of Banzhaf & Yamamoto, and its bibliography keys
d["archived"]     # None, or the archive group
```

`phenomena` is the most underused field: it is the list of behaviours the
literature reports, which is usually what someone actually wants when they ask
for "a chemistry that does X".
