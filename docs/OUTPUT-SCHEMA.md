# What a Chemart chemistry must output

**Short answer to the question you asked: a stoichiometric matrix plus a
species-name list is too restrictive.** It is the right *core*, but it
loses essential structure on **59 of the 98** catalogued chemistries (49
are constructive, so no finite matrix exists; 27 are catalytic, so the
net matrix erases the species that drive them), and it cannot express
the dynamics of a further handful at all. It has to become the core of a
layered record, and four of its implicit assumptions have to be broken
explicitly.

This document says which assumptions break, on which chemistries, and
proposes the record we should emit instead.

---

## 1. The four ways `(S, names)` fails

### 1.1 A single net matrix destroys catalysis — and catalysis is the point

The dominant reaction scheme in the whole field is

```
s1 + s2  ->  s1 + s2 + s3
```

Both reactants survive. In the net matrix `S = P - R` their columns are
zero, so they vanish — yet they are exactly what sets the propensity
`a = k·[s1]·[s2]`, and their identity is the entire content of the model.

This scheme is used by: `matrix-chemistry`, `alchemy`,
`automata-reaction`, `random-catalytic-networks`, `mcs-bl`,
`ikegami-hashimoto`, `combinator-chemistry`, `replicator-equation`,
`prime-number-chemistry`, `disperser`, `fraglets` (`matchp`),
`kauffman-autocatalytic-sets`, `bagley-farmer`, `gard`, `raf`,
`jain-krishna`, the CNO cycle in `nuclear-reaction-networks`, …

> **Requirement 1.** Emit the reactant matrix `R` and the product matrix
> `P` separately. `S = P − R` is a derived view, never the stored one.
> Catalysts are recoverable as `min(R, P) > 0` per reaction.

This is not academic: `raf` and every autocatalytic-set analysis in
chapter 6 is *defined* in terms of "which species catalyses which
reaction". That relation is unrecoverable from `S`.

### 1.2 For constructive chemistries there is no finite matrix at all

**49 of 98** catalog entries are constructive: `S` is unbounded and grows
at run time. `prime-number-chemistry` ranges over all naturals;
`alchemy` over all lambda normal forms; `stringmol`, `typogenetics`,
`squirm3`, `bnc-cell` over arbitrary-length strings.

Chapter 15 is explicit about why this matters: a birth event does not
just add a row, it changes the *dimensionality of state space*, and "a
treatment using differential equations is not possible across these
discontinuities."

So a chemistry's primary interface cannot be a matrix. It has to be:

```python
class Chemistry:
    def react(self, *reactants) -> list[Molecule] | None   # None = elastic
    def species_space(self) -> SpeciesSpace                # possibly infinite
    def expand(self, seed, *, max_species, max_rounds) -> ReactionNetwork
```

`expand()` is the closure operator `G_C` from §12.2 — repeatedly react
all combinations until no novel species appears — with explicit budgets,
and it must report *why* it stopped (`closed` vs `truncated`). A
truncated network is scientifically different from a closed one and the
downstream tools must be told which they got. The matrix chemistry is
the test case: from `s(1)..s(15)` at `N=9` it provably closes at
`s(1)..s(27)` minus `s(20..23)`, and we can assert that.

> **Requirement 2.** The chemistry object is the primitive; the network
> is a *product* of it, carrying `closure_status` and the seed and
> budgets that produced it.

### 1.3 Mass action is not the default; sometimes there is no stoichiometry at all

Chemistries whose defining dynamics is *not* mass action on `S`:

| chemistry | actual rate law |
|---|---|
| `arn` | normalised exponential response `ċ_i = δ(a_i−h_i)c_i / Σc_j`, with `a_i, h_i` exponential in bitstring complementarity |
| `hill-kinetics`, `repressilator` | Hill function `P^n/(K^n+P^n)` |
| `michaelis-menten` (abridged) | `v_m[S]/(k_m+[S])` — fewer species than the elementary form |
| `bigan-conservative-crn` | saturating kinetics modelling crowding |
| `farmer-immune` | rate set by the alignment score `m_ij` |
| `energy-gated-collision` | Arrhenius gate on kinetic vs activation energy |
| `swarm-chemistry` | a **force law**. No species is ever created or destroyed. `S` is the zero matrix. |
| `mechanical-self-assembly` | rate factorises into a mass-action collision term `P_c` and a *geometric* orientation term `P_b` |

`swarm-chemistry` is the clean reductio: it is in the book as an
artificial chemistry, it is fully specified, and its stoichiometric
matrix carries literally no information.

> **Requirement 3.** Rate laws are per-reaction, typed, and open:
> `mass-action | michaelis-menten | hill | saturating | arrhenius |
> custom(expr)`. A chemistry may legitimately have an empty stoichiometry
> and a non-empty dynamics.

### 1.4 Space, compartments and flow are part of the model, not decoration

- **Flow.** Almost every origin-of-life chemistry is *defined* to be open.
  Catalysis provably cannot shift an equilibrium, so `gard`,
  `bagley-farmer` and `kauffman-autocatalytic-sets` only do anything
  interesting under a food-set inflow. Chapter 13 goes further: with a
  leak, `{A,B}` is *not* an organisation even though it is closed and
  logically self-maintaining — you have to check `S·v ≥ 0` for some flux
  `v > 0`. Drop the flow and you get the wrong answer, silently.
- **Compartments.** `p-systems`, `cham`, `fraglets`, `mcs-bl`, `sac`,
  `chemoton`, `bnc-cell`, `gard` are hierarchical. A flat species list
  cannot say that species `a` in membrane 7 is a different pool from `a`
  in membrane 2.
- **Space.** 24 entries place molecules: `squirm3`,
  `ono-ikegami-protocell`, `flow-ac`, `oregonator`, `srsim`,
  `isologous-diversification`, `avida`, `sr-loops`, …

Counts across the catalog: **29** entries define a flow, **16** define
compartments, **24** define space.

> **Requirement 4.** Flow terms, compartment tree and spatial/diffusion
> data are first-class fields, not annotations.

---

## 2. Your kinetics/thermodynamics flag: yes, and it should be three tiers

You suggested flagging whether a chemistry only gives topology or also
provides parameter-allocation primitives. That is right, and the catalog
already carries it as the `provides` list. Rolled up, it gives three
tiers:

| tier | meaning | count | examples |
|---|---|---|---|
| **T — topology** | who reacts with whom; rates are yours to pick | all | `matrix-chemistry`, `alchemy`, `automata-reaction`, `kauffman-autocatalytic-sets`, `jain-krishna` |
| **K — kinetics** | the chemistry prescribes rate constants or a rate law | 38 | `farmer-immune` (rates from the match score), `arn` (from bitstring complementarity), `gard` (β from the receptor affinity distribution), `repressilator`, `brusselator`, `oregonator`, `disperser` |
| **H — thermodynamics** | per-species energies and/or reverse rates constrained by ΔG | 15 | `toychem` (extended-Hückel orbital energies → activation energies → Arrhenius), `bigan-conservative-crn` (random ΔG_f, backward rates by detailed balance, conservativity enforced via `Sᵀm = 0`), `dorin-korb-ecosystem` (bond energies; energy must be spent in the same timestep or is lost), `bagley-farmer`, `gard` (catalysis must accelerate both directions equally), `rna-folding-ac` (Vienna free-energy folding), `nuclear-reaction-networks`, `bondable-ca` (mean polarity as an energy-like observable), `urdar` (string entropy), `hbcb-psd`, `energy-gated-collision` |

Two things worth knowing about tier H, because they change the API:

1. **Thermodynamic consistency is a constraint between forward and
   reverse rates, not extra data.** In `bigan-conservative-crn` the
   forward constants are drawn log-uniformly and the *backward* ones are
   then computed so that detailed balance holds. So the schema needs
   reversible reaction pairs to be linkable, and a validator that can
   check `∏ k_f / ∏ k_r = exp(−ΔG/RT)` around every cycle.

2. **Mass conservation deserves its own field.** `toychem`, `squirm3`,
   `bigan-conservative-crn`, `n-economy`, `smn` and
   `nuclear-reaction-networks` all have an exact atom/mass vector `m`
   with `Sᵀm = 0`. Some hand it to you (`squirm3`: atoms are never
   created or destroyed; `n-economy`: the prime-exponent vector *is* the
   composition); for others we must compute the left nullspace. And note
   `chameleon`: its conservation law is a **mod-3 integer invariant**
   that a real-valued nullspace computation will not find — so
   conservation detection has to be done over the integers.

---

## 3. The implemented record

> This section originally proposed a larger record. What is implemented is
> the slimmer `chemart.network.Network` below (see `docs/PLAN.md`, revision
> notes 7); the requirements above are all still met.

```python
Species(id: str, structure: str | None)          # structure: bitstring, λ-term, SMILES, …
Reaction(reactants: dict[str, int],              # R column
         products:  dict[str, int],              # P column (catalysts appear on both sides)
         rate: dict | None,                      # {"law": "mass-action" | "michaelis-menten" |
                                                 #  "hill" | "saturating" | "arrhenius", ...}
         count: int | None)                      # firings, for simulation-observed networks
Network(species, reactions,
        status: "complete" | "truncated" | "observed",
        initial_state, inflow, outflow,          # dict species -> number (outflow may be a scalar)
        extras: dict,                            # reserved keys: space, compartments, energies,
                                                 # conservation, analysis, interaction_law
        chemistry: str, params: dict, seed: int | None)
# computed: provides (derived from content), matrices() -> (ids, R, P),
#           summary(), to_text(), to_dict() / from_dict()
```

Deliberate simplifications relative to the original proposal: rates are plain
dicts instead of a class hierarchy; reversible pairs are two reactions;
macroscopic vs mesoscopic constants are converted with
`chemart.kinetics.k_to_c` rather than stored as a flag; compartments, space
and energies live under `extras`; `provides` is computed, never stored.

Three notes on fields that are easy to get wrong:

- **`rate_scale` and `volume`.** The book devotes an appendix section to
  this: macroscopic `k` and mesoscopic `c` differ by the Wolkenhauer
  relation `c = k / (N_A V)^(m−1) · ∏ l_i!`, where `m` is the total
  educt count and `l_i` the educt multiplicities. Exporting a number
  without saying which one it is is a bug waiting to happen, and the
  `l_i!` term bites on `brusselator`'s trimolecular step `2X + Y → 3X`.
- **`species[].structure`.** For constructive chemistries the species
  *is* a lambda term / bitstring / graph, and it is the only thing that
  makes the network reproducible or extensible. Keep it. Note `gard`
  needs a *compositional* identity (a multiset of lipid types), not a
  sequence.
- **`provides`.** Consumers must be able to ask "does this network have
  real rate constants, or did Chemart pick 1.0 for everything?" before
  they run an analysis that depends on the answer.

---

## 4. What this buys us downstream

The extra fields are exactly what the interop targets need:

| target | needs beyond `(S, names)` |
|---|---|
| SBML / COPASI / Tellurium / libRoadRunner | separate R and P, rate laws, compartments, initial amounts, units |
| Gillespie SSA (StochPy, `GillespieVessel`) | mesoscopic `c`, volume, integer initial counts |
| COBRApy / FBA | `S`, and *bounds* derived from inflow/outflow |
| Chemical Organisation Theory | catalysts, and the flux condition `S·v ≥ 0, v > 0` — §13.1 shows closure alone gives wrong organisations |
| RAF detection | the catalysis relation, plus the food set |
| CRNT (Feinberg) / ERNEST | reversibility pairing and deficiency, which need R and P |
| BioNetGen / KaSim | the rule-based form *un*flattened — `kappa-calculus` and `srsim` must not be expanded |
| network analysis (NetworkX) | the bipartite species/reaction graph, which is R and P |

---

## 5. Consequences for the library

1. `Chemistry` (a generator with `react`/`expand`) and `ChemistryNetwork`
   (a materialised, possibly truncated instance) are separate types.
2. Reactions store `R` and `P`. `S` is a property.
3. Every network carries `provides`, and exporters **refuse or warn**
   rather than invent: exporting a T-tier network to SBML must either
   demand rate constants from the caller or stamp them as defaults.
4. A `params` role vocabulary (`structural`, `kinetic`, `thermodynamic`,
   `population`, `spatial`, `stochastic`, `selection`) is already in the
   catalog. It is what lets us offer cross-chemistry scaling: *"shrink
   every structural knob until |S| < 500"* works uniformly over
   `matrix-chemistry`'s `N`, `kauffman-autocatalytic-sets`'s `L`, and
   `quasispecies`'s `L`.
5. Analyses (organisations, RAFs, conservation laws, deficiency) are
   functions on `ChemistryNetwork`, and each declares which `provides`
   tags it requires. Organisation theory needs flow; RAF needs catalysts;
   thermodynamic validation needs tier H.

## 6. Known open questions

- **Higher-order chemistries.** `high-order-chem`, `gamma`'s
  γ-calculus and `kappa-calculus` let rules be molecules. A reaction can
  then *create a reaction*. The record above has no slot for that; the
  cheap fix is to allow a species to carry a rule as its structure and to
  let `expand()` grow the reaction list, but the exported network is then
  only a snapshot.
- **Emergent species.** In `sr-loops` and `ca-embedded-particles` "what
  is a molecule" is an observer's choice over CA patterns. Chemart cannot
  derive the species set; it must accept a user-supplied coarse-graining.
- **Non-CRN chemistries.** `swarm-chemistry`, `l-systems`,
  `nk-landscape` produce no meaningful CRN. They stay in the catalog for
  completeness, but should be reachable through a different interface and
  must not silently export an empty matrix.
