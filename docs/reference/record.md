# The network record

Every chemistry returns the same object. It is plain JSON throughout, so
`Network.from_dict(net.to_dict())` reproduces it exactly — that round trip is
the format's contract and is tested for every entry.

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
    chemistry: str = "",
    params: dict = {},
    seed: int | None = None,
)
```

`to_dict()` emits those keys in that order.

The constructor validates and raises `ValueError` listing what is wrong —
unknown species in a reaction, non-positive stoichiometry, a rate law outside
the vocabulary, a bad status. You cannot build an inconsistent network by
accident.

## Species

```python
Species(id: str, structure: str | None = None)
```

`structure` carries internal structure as a string when the molecule has any: a
bitstring, a λ-term, a program, a genome, a dot-bracket RNA fold, a graph
encoding. **49** entries carry it, and for those it is where the molecule
actually lives — without it they would be bare graphs.

## Reaction

```python
Reaction(reactants: dict[str, int], products: dict[str, int],
         rate: dict | None = None, count: int | None = None)

Reaction.of(["A", "A"], ["B"])     # build from id sequences
r.catalysts                         # species on both sides, surviving multiplicity
r.to_text()                         # "2 A -> B  [mass-action k=1.0]  (x17)"
```

`count` is the number of firings, present only when `status == "observed"`.

!!! warning "Build sides from lists, not dict literals"
    `Counter({a: 1, b: 1})` silently collapses to a single entry when
    `a == b`, losing identical reactants and clonal offspring.
    `Reaction.of([a, b], [c])` and `Counter([a, b])` are safe.

## status

| status | meaning | what a missing reaction tells you |
|---|---|---|
| `complete` | the whole defined network, or a closure that finished | it is not in the chemistry |
| `truncated` | a closure cut off by a size budget | nothing — raise the budget parameter |
| `observed` | the reactions that fired in one simulated run | nothing about the chemistry; only about that run and seed |

## Flows

- **`inflow`** maps species to a constant influx (amount per volume per time).
- **`outflow`** is a first-order removal rate: a per-species dict, a single
  number for every species, or the string `"constant-total"`.
- **`"constant-total"`** is the book's standard flow reactor — a non-selective
  dilution flux Φ that keeps total concentration constant. Import it as
  `from chemart.network import CONSTANT_TOTAL` rather than typing the literal.
- Species held at a fixed concentration are listed in `extras["buffered"]`,
  not expressed as flows.

## extras

Free-form JSON with seven reserved keys that have a fixed meaning; all but
`food` feed `provides`:

| key | holds |
|---|---|
| `space` | lattice, positions, diffusion, geometry |
| `compartments` | nested membranes or cells |
| `energies` | per-species or per-bond energies |
| `conservation` | list of laws; see [Comparing and analysing](../guide/analysing.md#conservation-laws) |
| `analysis` | measurements the generator made during its run |
| `interaction_law` | the law, when the model has no transformational reactions |
| `food` | the food set: species supplied from outside, which organisation measures (RAFs, scope, flux modes) start from. Without it, `chemart.measures` takes the species with inflow, then the buffered ones, then those in the initial state |

Anything else is chemistry-specific and documented in that entry. Ones you will
meet often: `buffered`, `final_state`, `events`, `rules`, `reaction_rules`,
`species_encoding`, `published`.

## provides

`net.provides` is derived from content on every access:

- `topology` always; `stoichiometry` if there are reactions
- `catalysts` if any reaction has a species on both sides
- `rate-constants` if any reaction has a rate; `rate-law` if any rate is not
  mass-action
- `flow` if `inflow` or `outflow`; `initial-state` if `initial_state`
- `space` / `compartments` / `energies` / `mass-conservation` from the matching
  reserved `extras` keys

The catalog entry's `provides` is a **claim about the chemistry**. The contract
enforces computed ⊆ claimed, so the claim may be a superset — typically because
a capability only appears at non-default parameters, or because the tag is not
computable at all (`sequence-structure-function` never is).

## matrices()

```python
ids, R, P = net.matrices()    # scipy.sparse CSC, species x reactions, int64
S = P - R                     # net stoichiometry
```

## Rate laws

`rate` is `None` (topology only) or a dict whose `law` is one of six. Extra keys
beyond the required ones are allowed, and carry facts that are not part of the
propensity — a threshold, a temperature, units.

| law | required keys | propensity |
|---|---|---|
| `mass-action` | `k` | `k * prod(x_i ** n_i)` |
| `power` | `k`, `order` | `k * x ** order` |
| `michaelis-menten` | `vmax`, `km` | `vmax * x / (km + x)` |
| `hill` | `vmax`, `K`, `n` | `vmax * h` (activation) or `vmax * (1 - h)`, where `h = p**n / (K**n + p**n)` and `p` is the regulator |
| `saturating` | `k`, `K` | `k * prod((x_i / (1 + x_i / K)) ** n_i)` |
| `arrhenius` | `A`, `Ea` | `A * exp(-Ea / (R * T))` |

!!! warning "Two traps"
    Michaelis–Menten uses lowercase **`km`**, while Hill and saturating use
    uppercase **`K`**. And Hill rates carry extra keys `regulator` (the species
    id driving them) and `mode` (`activation` or `repression`), because the
    regulator is usually not a reactant.

!!! info "Only two of the six laws are currently emitted"
    Measured across the default networks: **mass-action**, and **arrhenius**
    (`toychem`, and the archived `energy-gated-collision`). Nothing emits
    `michaelis-menten`, `hill`, `saturating` or `power` — notably
    `michaelis-menten` and `hill-kinetics` give the *elementary* mechanism
    rather than the abridged law. The four unused laws remain in the vocabulary
    because a chemistry may legitimately prescribe one, and the validator
    accepts them; an integrator only needs mass-action to cover everything
    generated today.

`arrhenius` is in the vocabulary but needs a temperature and gas constant that
are not part of the law, because `chemart.kinetics` is deliberately
unit-agnostic. `chemart.simulate` applies it when they are known: on the rate
dict as `T` and `R`, or passed as `temperature=` and `gas_constant=`.

!!! note "Never invent a rate"
    If a source gives no rate constant, the entry uses `None`. That is
    information, not an omission waiting to be filled in.

## Units of amounts

Amounts in `initial_state`, and in the frames of a trajectory, are per unit
volume. The stochastic simulator converts them to molecule counts as
`count = amount × volume × avogadro`, so with its defaults (`volume=1`,
`avogadro=1`) amounts are counts, which is also how the Turing gases report
their soups.

## The trajectory record

Every run returns a `Trajectory`: integrating a network
(`chemart.simulate.ode`), sampling it stochastically (`chemart.simulate.ssa`), or
evolving a Turing gas. Like the network, it is plain JSON and
`Trajectory.from_dict(traj.to_dict())` reproduces it exactly.

```python
Trajectory(
    network: Network,          # the simulated network, or the gas's observed network
    frames: list[Frame],
    method: "ode" | "ssa" | "evolve",
    clock: str = "time",       # the unit of frame time: time, collisions, epochs, ...
    settings: dict = {},       # provenance: seed, solver, rate and state specs, ...
)

Frame(
    t: float,
    state: dict[str, float],   # amount per species; zeros are left out
    fired: list = [],          # [[reactant ids], [product ids], count] since the previous frame
    observables: dict = {},    # quantities the chemistry itself reports
)
```

The state is sparse, one dict per frame, so a soup with thousands of species
costs no more than the species actually present. For plotting,
`traj.array(species)` gives dense arrays `(ids, t, X)`. `traj.series(name)` is an
observable over time, `traj.window(i, width)` the network of the reactions fired
in a range of frames, and `traj.turnover()` a clock shared by every gas: the
cumulative number of reactions per molecule.

## Building a network by hand

For your own chemistry, or to compare against a catalogued one:

```python
from chemart.helpers.explicit import network

net = network([
    ("A -> X",         1.0),                              # mass action, k = 1.0
    ("2 X + Y -> 3 X", {"law": "mass-action", "k": 2.0}),
    ("X -> ",          None),                             # empty side = ∅, no rate
], initial_state={"A": 1.0, "X": 1.0}, extras={"buffered": ["A"]})
```

Terms are separated by `" + "` **with spaces**, so species ids may contain `+`
(`e+`). A coefficient is separated from the species by a space (`2 X`), so ids
may start with digits (`12C`). Species appear in the order given, then in order
of appearance.
