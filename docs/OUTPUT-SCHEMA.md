# The output format

Every chemistry in Chemart returns the same record. This page states what that
record is, then explains why it is shaped that way — because the obvious choice,
a stoichiometric matrix plus a list of species names, turns out to lose
essential structure on most of the catalog.

## What a chemistry hands back

```python
Species(id: str, structure: str | None = None)

Reaction(reactants: dict[str, int],       # the R column
         products:  dict[str, int],       # the P column; catalysts appear on both sides
         rate: dict | None = None,        # {"law": "mass-action", "k": 1.0}, or None
         count: int | None = None)        # firings, when the network was observed

Network(species, reactions,
        status: "complete" | "truncated" | "observed",
        initial_state, inflow, outflow,
        extras,                           # space, compartments, energies,
                                          # conservation, analysis, interaction_law
        chemistry, params, seed)          # provenance: enough to reproduce it
```

Derived views, computed rather than stored: `provides`, `matrices()` returning
`(ids, R, P)`, `summary()`, `to_text()`, and `to_dict()` / `from_dict()`.

Everything is plain JSON — no numpy scalars, no objects — so
`Network.from_dict(net.to_dict())` round-trips exactly. That round trip is the
format's contract and is tested for all 98 entries.

Full field-by-field detail is in [The network record](reference/record.md).

## Why not just a stoichiometric matrix?

A single net matrix `S` with species names is the right *core*, and it is what
most CRN tooling consumes. But four of its implicit assumptions break on this
catalog, and each break is load-bearing for a large group of entries.

### 1. A net matrix erases catalysis, and catalysis is usually the point

The most common reaction scheme in the whole field is

```
s1 + s2  ->  s1 + s2 + s3
```

Both reactants survive. In `S = P − R` their entries are zero, so they vanish
from the record — yet they are exactly what sets the propensity
`k·[s1]·[s2]`, and *which* species catalyses *which* reaction is frequently the
entire content of the model.

**74 of 98** entries produce catalytic reactions. Whole analyses are defined in
terms of that relation: autocatalytic-set detection and the RAF algorithm ask
"which species catalyses which reaction", and the answer is unrecoverable from
`S` alone.

> **Consequence.** Chemart stores the reactant and product multiplicities
> **separately**. `S = P − R` is a view, never the stored form, and catalysts
> are recoverable as the species present on both sides.

### 2. For constructive chemistries there is no finite matrix

**51 of 98** entries are *constructive*: the species set is open and grows as
reactions produce new molecules. `prime-number-chemistry` ranges over the
naturals; `alchemy` over lambda normal forms; `stringmol`, `typogenetics` and
`squirm3` over arbitrary-length strings.

A birth event does not merely add a row — it changes the dimensionality of the
state space, which is why a differential-equation treatment cannot cross those
discontinuities.

So a chemistry's primitive cannot be a matrix. In Chemart the primitive is a
**generator function**, and the network is a product of running it:

```python
def generate(p, rng) -> Network      # what every chemistry implements
```

For constructive entries, the interesting object is a *closure* — react every
combination until nothing novel appears — computed by
`chemart.expand.expand(react, seed, max_species=…)`. A closure that hits its
budget is scientifically different from one that finished, so the network says
which it is.

> **Consequence.** `status` is part of the record: `complete`, `truncated`, or
> `observed`. A missing reaction in a `truncated` network means the budget
> stopped early; in an `observed` one it means that run didn't fire it.

### 3. Mass action is not the default, and sometimes there is no stoichiometry

Chemistries whose defining dynamics is not mass action on `S`:

| chemistry | actual rate law |
|---|---|
| `hill-kinetics`, `repressilator` | Hill function |
| `michaelis-menten` | saturating enzyme kinetics |
| `bigan-conservative-crn` | saturating kinetics modelling crowding |
| `farmer-immune` | rate set by a sequence alignment score |
| `arn` | response exponential in bitstring complementarity |
| `energy-gated-collision` | an Arrhenius gate on collision energy |
| `swarm-chemistry` | a **force law** — nothing is created or destroyed |

`swarm-chemistry` is the clean reductio: it is in the book as an artificial
chemistry, it is fully specified, and its stoichiometric matrix carries no
information at all. It is the one entry whose default network has **zero**
reactions, with the law in `extras["interaction_law"]`.

> **Consequence.** Rates are per-reaction and drawn from a fixed vocabulary
> (`mass-action`, `power`, `michaelis-menten`, `hill`, `saturating`,
> `arrhenius`), or `None`. An empty reaction list with a non-empty dynamics is
> legal.

A related rule: where a source publishes no rate constant, the entry stores
`None` rather than a plausible-looking number. That absence is information.

### 4. Flow, compartments and space are part of the model

- **Flow** — **23** entries define an open system. Catalysis cannot shift an
  equilibrium, so `gard`, `bagley-farmer` and `kauffman-autocatalytic-sets` only
  do anything interesting under a food-set inflow. Organisation theory is
  stricter still: with a leak, a set can be closed and self-maintaining and yet
  not be an organisation, because that also requires a flux `v > 0` with
  `S·v ≥ 0`. Drop the flow and the answer is silently wrong.
- **Compartments** — **10** entries are hierarchical (`p-systems`, `cham`,
  `fraglets`, `chemoton`, `gard`, …). A flat species list cannot express that
  species `a` inside membrane 7 is a different pool from `a` in membrane 2.
- **Space** — **21** entries place molecules on a lattice or in continuous
  space.

> **Consequence.** These are first-class fields (`inflow`, `outflow`) or
> reserved `extras` keys (`space`, `compartments`), not annotations.

## Knowing what you actually got

Consumers need to ask, before running an analysis, whether a network has real
rate constants or whether Chemart picked 1.0 for everything. `net.provides`
answers that, computed from content on every access. Rolled up, the catalog
falls into three tiers:

| tier | meaning | entries |
|---|---|---|
| **topology** | who reacts with whom; rates are yours to choose | all 98 |
| **kinetics** | the chemistry prescribes rate constants or a rate law | 44 |
| **thermodynamics** | per-species energies, or reverse rates constrained by ΔG | 11 |

Two details about the thermodynamic tier change how it must be stored:

- **Consistency is a constraint between forward and reverse rates, not extra
  data.** In `bigan-conservative-crn` the forward constants are drawn
  log-uniformly and the backward ones computed so detailed balance holds. Only
  3 entries claim `thermodynamic-consistency`.
- **Mass conservation deserves its own field.** **31** entries carry an exact
  atom or mass vector `m` with `Sᵀm = 0`, declared in
  `extras["conservation"]`. Some laws are **modular** rather than zero —
  `chameleon`'s hold mod 3 — so a real-valued nullspace computation will not
  find them, and a plain zero test reports them as violated.

## What the record deliberately does not do

The format was simplified during implementation, and the simplifications are
worth knowing because they are choices, not omissions:

- Rates are plain dicts from a fixed vocabulary, not a class hierarchy.
- Reversible reactions are two reactions, not one flagged pair.
- Macroscopic and mesoscopic rate constants are not stored as separate fields;
  `chemart.kinetics.k_to_c` converts between them on demand. The conversion
  carries a combinatorial factor for identical reactants, which bites on
  `brusselator`'s trimolecular step `2X + Y → 3X`.
- Space, compartments and energies live under `extras` rather than as typed
  blocks, because their shape differs too much between entries to fix in a
  schema.
- `provides` is computed, never stored, so it cannot drift from the content.
- There is no separate generator class and no materialised-network class. A
  chemistry is a module with one function; the network is a dataclass.

## What this buys downstream

The extra structure is exactly what interoperability targets need beyond
`(S, names)`:

| target | needs |
|---|---|
| SBML, COPASI, Tellurium | separate R and P, rate laws, compartments, initial amounts |
| Gillespie SSA | mesoscopic constants, volume, integer counts |
| FBA / COBRApy | `S`, plus bounds derived from inflow and outflow |
| chemical organisation theory | catalysts, and the flux condition `S·v ≥ 0, v > 0` |
| RAF detection | the catalysis relation, plus the food set |
| CRNT / deficiency analysis | reversibility pairing, which needs R and P |
| rule-based tools (BioNetGen, KaSim) | the rules *unflattened* — `kappa-calculus` and `srsim` must not be expanded |
| graph analysis | the bipartite species/reaction graph, which is R and P |

Exporters for these targets are not implemented. The record was designed so
that writing them is mechanical, and an exporter should refuse — or clearly
stamp its assumptions — rather than invent rate constants for a topology-only
network.

## Open questions

- **Higher-order chemistries.** `high-order-chem`, γ-calculus and
  `kappa-calculus` allow rules to be molecules, so a reaction can create a
  reaction. The record has no slot for that; the current approach lets a species
  carry a rule as its `structure`, which makes any exported network a snapshot.
- **Emergent species.** In `sr-loops` and `ca-embedded-particles`, "what counts
  as a molecule" is an observer's choice over cellular-automaton patterns.
  Chemart cannot derive the species set; it applies a documented
  coarse-graining and says so.
- **Non-CRN entries.** `swarm-chemistry` and the analysis-kind entries produce
  no meaningful CRN on their own. They stay in the catalog for completeness and
  expose their content through `extras` rather than exporting an empty matrix.
- **The rate-law vocabulary versus the ODE helper.** `arrhenius` is in the
  vocabulary but needs a temperature and gas constant the rate dict does not
  carry, so generic integrators cannot apply it. Whether to add a
  catalysed-Michaelis-Menten law, and how far the vocabulary and the integrator
  should track each other, is recorded as an open decision in `to_decide.md`.
