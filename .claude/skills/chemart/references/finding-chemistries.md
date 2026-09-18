# Finding the right chemistry

98 entries. This is how to narrow them down without reading all of them.

## Contents
- [The four axes](#the-four-axes)
- [Querying from Python](#querying-from-python)
- [Querying from the CLI](#querying-from-the-cli)
- [What each capability tag means](#what-each-capability-tag-means)
- [Families, with counts](#families-with-counts)
- [Reading a describe_chemistry result](#reading-a-describe_chemistry-result)

## The four axes

**`provides`** — what the generator can actually hand you. This is the axis
that matters most for analysis: if you need rate constants to simulate, filter
on `rate-constants`, not on the chemistry's reputation.

**`kind`** — what the thing *is*: `generator` (80), `formalism` (7, e.g. Gamma
and P systems: a rule language, where you supply the rules), `analysis` (4,
e.g. RAF: builds a model *and* runs an algorithm on it), `wet` (4: the
in-silico network of a published laboratory experiment), `framework` (3:
larger simulation platforms, minimally re-implemented).

**`family`** — subject area (see below).

**`constructive`** — 51 of 98. True means the species set is open and grows at
run time, so the interesting object is a *closure* and the network may come
back `truncated`. False means the species set is fixed and enumerable.

## Querying from Python

```python
from chemart.catalog import load

# load() parses the YAML once and caches it; each call gets its own list,
# so it is cheap to call repeatedly.

# everything that supplies kinetics AND a conservation law
for c in load():
    if "rate-constants" in c.provides and "mass-conservation" in c.provides:
        print(c.id, c.fidelity)

# constructive chemistries, with whatever knob sizes each one
for c in load():
    if c.constructive:
        knobs = [p.name for p in c.params_by_role("structural")]
        print(f"{c.id:28s} {knobs}")
```

Useful `Chemistry` fields: `id`, `name`, `family`, `kind`, `constructive`,
`provides`, `fidelity`, `sources`, `decisions`, `phenomena`, `book`, `refs`,
`params`, and the properties `implemented`, `module`, `tiers`,
`params_by_role(role)`.

Parameter objects carry `name`, `type`, `default`, `min`, `max`, `choices`,
`meaning`, `role`, `range`.

## Querying from the CLI

```bash
uv run python -m chemart.catalog query --provides rate-constants
uv run python -m chemart.catalog query --family origin-of-life
uv run python -m chemart.catalog show gard          # human-readable entry
uv run python -m chemart.catalog status             # implementation + fidelity counts
uv run python -m chemart.catalog validate           # invariants; should print 0 problems
```

Or use the bundled `scripts/survey.py`, which does the common intersections
and prints a table.

## What each capability tag means

| tag | meaning | count |
|---|---|---|
| `topology` | species + who reacts with whom | 98 |
| `stoichiometry` | separate reactant/product multiplicities | 98 |
| `catalysts` | some species appears on both sides and must be kept | 74 |
| `initial-state` | the chemistry prescribes a starting multiset | 84 |
| `sequence-structure-function` | molecules carry structure the rule reads | 49 |
| `rate-constants` | the chemistry itself prescribes k's | 43 |
| `mass-conservation` | an atom/mass vector m with Sᵀm = 0 | 31 |
| `flow` | inflow/outflow is part of the definition | 23 |
| `space` | positions, lattice or diffusion | 21 |
| `compartments` | nested membranes/cells | 10 |
| `energies` | per-species free or bond energies | 10 |
| `rate-law` | a non-mass-action propensity is part of the model | 7 |
| `thermodynamic-consistency` | reverse rates constrained by ΔG | 3 |

The three-tag rule of thumb: `topology` alone means you can study structure;
add `rate-constants` and you can simulate; add `energies` or
`thermodynamic-consistency` and you can do thermodynamics.

## Families, with counts

- **systems-biology** (14) — aevol, bnc-cell, cpm-grn-evodevo,
  energy-gated-collision, french-flag, hbcb-psd, hill-kinetics,
  isologous-diversification, michaelis-menten, rna-folding-ac, smn, srsim,
  synthon, tominaga-stacked-strings
- **application** (13) — acgp, analog-function-crn, brusselator, ccm,
  disperser, fraglets, metabolic-robot-controller, molecular-tsp, music-ac,
  naming-game-ac, okamoto-switch, organization-computing, proof-ac
- **evolutionary-dynamics** (13) — ecolab, evolve-series, jain-krishna,
  logistic-chemistry, lotka-volterra, nk-landscape, quasispecies,
  random-catalytic-networks, rbn, replication-death, replicator-equation,
  selection-equation, urdar
- **automata** (12) — automata-reaction, avida, bondable-ca,
  ca-embedded-particles, corewar, coreworld, ikegami-hashimoto,
  laing-molecular-machines, mccaskill-polymer-tm, sr-loops, tierra,
  typogenetics
- **rewriting** (11) — alchemy, arms, brane-calculi, cham,
  combinator-chemistry, gamma, kappa-calculus, l-systems, mgs, p-systems,
  reflexive-ac
- **origin-of-life** (7) — autopoiesis-vmu, bagley-farmer, chemoton, gard,
  kauffman-autocatalytic-sets, ono-ikegami-protocell, raf
- **core** (5) — chameleon, dimerization, high-order-chem, matrix-chemistry,
  prime-number-chemistry
- **bio-inspired** (5) — conrad-enzymatic, farmer-immune, mcs-bl, sac, stringmol
- **non-chemical** (5) — mechanical-self-assembly, n-economy,
  nuclear-reaction-networks, soas, social-communication-ac
- **wet** (5) — dna-automaton, dna-hpp, oregonator, repressilator,
  self-propelled-droplets
- **network** (4) — arn, bigan-conservative-crn, nac, toychem
- **spatial** (4) — dorin-korb-ecosystem, flow-ac, squirm3, swarm-chemistry

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
```

`phenomena` is the most underused field: it is the list of behaviours the
literature reports, which is usually what someone actually wants when they ask
for "a chemistry that does X".
