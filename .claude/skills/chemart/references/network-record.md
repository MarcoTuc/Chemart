# The network record

Every chemistry returns the same object. It is plain JSON data throughout, so
`Network.from_dict(net.to_dict())` reproduces it exactly — that round trip is
the format's contract and is tested for all 98 entries.

## Contents
- [Fields](#fields)
- [Species and Reaction](#species-and-reaction)
- [status: what the network means](#status-what-the-network-means)
- [Flows](#flows)
- [extras](#extras)
- [provides: claimed vs computed](#provides-claimed-vs-computed)
- [Matrices](#matrices)
- [Rate laws](#rate-laws)
- [Building a network by hand](#building-a-network-by-hand)

## Fields

```python
Network(
    species: list[Species],
    reactions: list[Reaction],
    status: "complete" | "truncated" | "observed" = "complete",
    initial_state: dict[str, float] | None = None,
    inflow: dict[str, float] | None = None,
    outflow: dict[str, float] | float | str | None = None,
    extras: dict = {},
    # provenance, filled in by generate_network:
    chemistry: str = "", params: dict = {}, seed: int | None = None,
)
```

`to_dict()` keys, in order: `species`, `reactions`, `status`, `initial_state`,
`inflow`, `outflow`, `extras`, `chemistry`, `params`, `seed`.

The constructor validates and raises `ValueError` listing what is wrong —
unknown species in a reaction, non-positive stoichiometry, a rate law outside
the vocabulary, a bad status. You cannot build an inconsistent network by
accident.

## Species and Reaction

```python
Species(id: str, structure: str | None = None)
```

`structure` carries internal structure as a string when the chemistry has any —
a bitstring, a λ-term, a genome, a dot-bracket RNA fold, a graph encoding. For
the 49 entries tagged `sequence-structure-function` this is where the actual
molecule lives, and it is what makes those chemistries more than graphs.

```python
Reaction(reactants: dict[str, int], products: dict[str, int],
         rate: dict | None = None, count: int | None = None)

Reaction.of(["A", "A"], ["B"])        # build from id sequences
r.catalysts                            # species on both sides, surviving multiplicity
r.to_text()                            # "2 A -> B  [mass-action k=1.0]  (x17)"
```

`count` is the number of firings, present only when `status == "observed"`.

Build reaction dicts from **lists**, not dict literals: `Counter({a: 1, b: 1})`
silently collapses to one entry when `a == b`, which loses clonal offspring and
identical reactants. `Reaction.of([a, b], [c])` is safe.

## status: what the network means

| status | meaning | what absence of a reaction tells you |
|---|---|---|
| `complete` | the whole defined network, or a closure that finished | it is not in the chemistry |
| `truncated` | a closure cut off by a size budget | nothing — raise the budget parameter and look again |
| `observed` | the reactions that fired in one simulated run | nothing about the chemistry; only about that run and seed |

This distinction matters more than it looks. Counting reactions across a
`complete` and an `observed` network and comparing the numbers is meaningless:
one is a definition, the other is a sample.

## Flows

- `inflow` maps species to a constant influx (amount per volume per time).
- `outflow` is a first-order removal rate: a dict per species, a single number
  for all species, or the string `"constant-total"`.
- `"constant-total"` is the book's standard flow reactor: a non-selective
  dilution flux Φ(t) that holds total concentration constant,
  `x_i' = production_i − x_i·Φ/Σx_k`. Import the constant as
  `from chemart.network import CONSTANT_TOTAL` rather than typing the string.
- Species held at fixed concentration are listed in `extras["buffered"]`, not
  expressed as flows.

`scripts/simulate.py` applies all of these, so do not add flow terms yourself
on top of them.

## extras

Free-form JSON, with six reserved keys that have fixed meaning and feed
`provides`:

| key | holds |
|---|---|
| `space` | lattice, positions, diffusion, geometry |
| `compartments` | nested membranes / cells |
| `energies` | per-species or per-bond energies |
| `conservation` | list of `{"name", "vector", optional "modulus"}` |
| `analysis` | measurements the generator made during its run |
| `interaction_law` | the law, when the model has no transformational reactions |

Anything else is chemistry-specific and documented in that entry. Common ones
you will meet: `buffered`, `final_state`, `events`, `rules`, `reaction_rules`,
`species_encoding`, `published`.

## provides: claimed vs computed

`net.provides` is derived from content on every access:

- `topology` always; `stoichiometry` if there are reactions
- `catalysts` if any reaction has a species on both sides
- `rate-constants` if any reaction has a rate; `rate-law` if any rate is not
  mass-action
- `flow` if inflow or outflow; `initial-state` if initial_state
- `space` / `compartments` / `energies` / `mass-conservation` from the matching
  reserved extras keys

The catalog entry's `provides` is a *claim* about the chemistry. The contract
test enforces computed ⊆ claimed, so the claim may be a superset — typically
because a capability appears only at non-default parameters, or because the
chemistry genuinely carries structure that the tag vocabulary can't compute
(`sequence-structure-function` is never computed).

When you need to know what you have, read `net.provides`. When you need to know
what the chemistry can do, read the catalog and the parameter `range` notes.

## Matrices

```python
ids, R, P = net.matrices()     # scipy.sparse CSC, species x reactions, int64
S = (P - R)                    # net stoichiometry
```

Conservation laws are vectors m with `Sᵀm = 0`. 31 entries claim
`mass-conservation`; 30 of them declare laws in `extras["conservation"]` at
default parameters. Verifying one is cheap and worth doing whenever you are
unsure a network is well formed — but read the payload defensively:

```python
import numpy as np
ids, R, P = net.matrices()
S = (P - R).toarray().astype(int)

for i, law in enumerate(net.extras.get("conservation", [])):
    v = law["vector"]
    m = np.array([v.get(s, 0) for s in ids]) if isinstance(v, dict) else np.array(v)
    residual = S.T @ m
    modulus = law.get("modulus")                    # some laws are modular
    holds = (residual % modulus == 0).all() if modulus else not residual.any()
    assert holds, law.get("name", f"law {i}")
```

Measured across all 115 declared laws in the library:

- `vector` is always present, as a **dict** keyed by species id (87 laws) or a
  **list** aligned with `matrices()`'s `ids` (28). Handle both.
- `name` is *usually* present but not guaranteed — `chameleon`'s three laws
  have none. Default it rather than indexing.
- `modulus` makes the law **modular**: it holds when `Sᵀm ≡ 0 (mod modulus)`,
  not when it is zero. `chameleon` uses modulus 3, so a plain zero test
  reports its correct laws as violated.
- `unique` and `notes` also occur. Treat unknown keys as annotation.

With that reader, every declared law in the library verifies.

## Rate laws

`rate` is `None` (topology only) or a dict whose `law` is one of six. Extra
keys beyond the required ones are allowed and carry facts that are not part of
the propensity (a threshold, a temperature, units).

| law | required keys | propensity |
|---|---|---|
| `mass-action` | `k` | k·Πxᵢ^nᵢ |
| `power` | `k`, `order` | k·x^order |
| `michaelis-menten` | `vmax`, `km` | vmax·x/(km+x) |
| `hill` | `vmax`, `K`, `n` | vmax·h or vmax·(1−h), h = pⁿ/(Kⁿ+pⁿ) |
| `saturating` | `k`, `K` | k·Π(xᵢ/(1+xᵢ/K))^nᵢ |
| `arrhenius` | `A`, `Ea` | A·exp(−Ea/RT) |

Two traps worth knowing: Michaelis-Menten uses lowercase **`km`** while Hill
and saturating use uppercase **`K`**; and Hill rates carry extra keys
`regulator` (the species id whose concentration drives it) and `mode`
(`activation` or `repression`), because the regulator is generally not a
reactant.

Measured across all 98 default networks, only two of these six laws actually
appear: mass-action (39 entries) and arrhenius (2). Nothing emits
`michaelis-menten`, `hill`, `saturating` or `power` — `michaelis-menten` and
`hill-kinetics` give the elementary mechanism instead of the abridged law. The
unused laws stay in the vocabulary because a chemistry may legitimately
prescribe one.

`arrhenius` is in the vocabulary but the bundled simulator cannot integrate it:
temperature and the gas constant are not part of the law, and `chemart.kinetics`
is deliberately unit-agnostic. Entries using it validate their kinetics in
closed form instead. This is an open design question in the repository's
`to_decide.md`.

Never invent a rate. If a source gives no rate constant, the entry uses `None`,
and that is information — not an omission to be filled in.

## Building a network by hand

For your own chemistry, or to compare against a catalogued one:

```python
from chemart.helpers.explicit import network

net = network([
    ("A -> X",           1.0),                              # mass action, k=1.0
    ("2 X + Y -> 3 X",   {"law": "mass-action", "k": 2.0}),
    ("X -> ",            None),                             # empty side = ∅, no rate
], initial_state={"A": 1.0, "X": 1.0}, extras={"buffered": ["A"]})
```

Terms are separated by `" + "` with spaces, so species ids may contain `+`
(`e+`). A coefficient is separated by a space (`2 X`), so ids may start with
digits (`12C`). Species appear in the order given, then in order of appearance.
