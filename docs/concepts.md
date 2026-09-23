# How it works

This page is the one to read if you want to *understand* Chemart rather than
just call it.

## What an artificial chemistry is

Banzhaf & Yamamoto define an artificial chemistry as a triple **(S, R, A)**:

| | |
|---|---|
| **S** | the set of possible **molecules** |
| **R** | the **reaction rules** — what happens when molecules meet |
| **A** | the **algorithm**, or reactor — how the rules are actually applied |

The breadth of that definition is the whole problem. In this catalog, a molecule
can be an integer, a bitstring, a λ-term, a lambda-like combinator, a graph, an
RNA sequence, a finite-state machine, or an oil droplet. A reactor can be a
well-stirred multiset, a system of differential equations, a 2D lattice, a
virtual machine with a scheduler, or a laboratory bench.

So the Brusselator and AlChemy are both artificial chemistries, and they have
essentially nothing in common as software.

## Why they can still share an output

Whatever the internals, a chemistry can be asked: *which species exist, and
which reactions relate them?* That is a **chemical reaction network** (CRN), and
it is the thing Chemart standardises.

This is not a lossless projection, and it is not meant to be. AlChemy's CRN
does not contain the λ-calculus machinery that reduces one term applied to
another. What it does contain is the reactions that actually fired in a run,
with counts — which is enough to compare AlChemy with any other chemistry by
the same measures.

Two consequences worth internalising:

- **A network is a view, not the model.** The generator is the model; the
  network is what it hands you.
- **A network can have no reactions at all.** Swarm Chemistry, kept in the
  [archive](catalog/index.md#archive), defines motion, not transformation. Its
  network has species and an empty reaction list, with the law in
  `extras["interaction_law"]`. That is a correct answer, not a bug.

## Three types of chemistry

Chemistries differ in where their network comes from, and that decides how
you work with them. Every catalog entry has a `type`:

| type | where the network comes from | what you do with it | examples |
|---|---|---|---|
| **given** | written down: the topology is fixed | pick rates and a starting state, then [simulate](guide/simulating.md) by rate equations or stochastically | brusselator, michaelis-menten, repressilator |
| **generator** | computed by an algorithm from its arguments; once built, it is treated as given | vary the arguments and see how the network's [measures](guide/measures.md) change (`measures.sweep`); simulate each network | kauffman-autocatalytic-sets, raf, random-catalytic-networks |
| **gas** | a Turing gas: structured molecules react by a procedure, so reactions exist only as they happen | [evolve](guide/evolving.md) it and follow its measures in chemical-evolutionary time | alchemy, bff, combinatory-chemistry |

The type is a judgement about the chemistry, recorded with a reason per entry
in `catalog/TYPES.md`. What you can *call* on a chemistry follows from its
code instead: its **faces**.

- The **generate** face (`chemart.generate_network`) returns one network.
- The **evolve** face (`chemart.evolve`) runs a process and returns a
  trajectory: a list of frames, each with the population, the reactions fired
  since the last frame and the chemistry's own observables.

Most gases have both faces. Their generate face returns the closure of their
rule from a seed set, and their evolve face runs the soup. A few given
chemistries also define a reactor worth running, such as the lattice of
autopoiesis-vmu, so they have an evolve face too.
`describe_chemistry(id)["faces"]` lists the faces of a chemistry.

## What you can do with a chemistry

| | given and generator | gas |
|---|---|---|
| **get a network** | `generate_network` | `generate_network` (the closure, or the network a run observed) |
| **dynamics** | `chemart.simulate.ode` or `.ssa` | `chemart.evolve` |
| **measures** | `chemart.measure(net)`; `measures.sweep` over arguments | `measures.over(traj)`, frame by frame |

All of it is also available from the command line (`chemart simulate`,
`chemart evolve`, `chemart measure`), as LLM tools, and in the
[simulation pit](hub.md#the-simulation-pit), a local app that plots runs as
they happen.


## The three moving parts

### 1. The catalog is the specification

Each chemistry is one YAML file in `catalog/chemistries/<id>.yaml`. It records
the (S, R, A) triple in prose, the parameters — type, default, bounds, meaning,
role, and the published range — what the generator can supply, the book section
and papers, and every decision taken where the sources were incomplete.

Crucially, the catalog is the **only** parameter specification. There is no
second copy in Python to drift from it:

- `describe_chemistry(id)` builds a JSON Schema from the YAML;
- `generate_network(id, **params)` validates your arguments against that YAML;
- the chemistry's own code receives parameters already validated.

A consequence you will feel: parameter errors are specific, because the catalog
knows what each knob means.

```python
chemart.generate_network("brusselator", b="not a number")
# ValueError: brusselator: b='not a number' is invalid: expected a number.
#   b: [B], held constant; sustained oscillations require b > 1 + a^2 ...
```

### 2. One or two functions per chemistry

```python
def generate(p, rng) -> Network
def evolve(p, rng):             # yields Frames, returns the observed Network
```

That is the entire contract. `p` holds the validated parameters as attributes;
`rng` is a seeded NumPy generator, used for *all* randomness so that a seed
reproduces a network, or a run, exactly. A module defines one face or both. A
parameter only one face uses says so in the catalog (`face: evolve`), and an
evolve face needs a `clock`, the unit its time is counted in.

There are no base classes, no registration and no plugin system. The id
`matrix-chemistry` resolves to `chemart.chemistries.matrix_chemistry` by naming
convention, imported lazily. Adding a chemistry means adding three files and
editing nothing else — see [Contributing a chemistry](contributing.md).

### 3. Everything returns the same record

```python
Network(
    species,        # [Species(id, structure)]
    reactions,      # [Reaction(reactants, products, rate, count)]
    status,         # "complete" | "truncated" | "observed"
    initial_state, inflow, outflow,
    extras,         # space, compartments, energies, conservation, analysis, ...
    chemistry, params, seed,
)
```

It is plain JSON throughout, so `Network.from_dict(net.to_dict())` round-trips
exactly — tested for every chemistry. A simulation or a run returns a
`Trajectory`: that network, plus its frames and the settings that produced
them, also plain JSON. Full detail in [The network record](reference/record.md).

`Species.structure` is where the molecule actually lives when it has internal
structure: the bitstring, the λ-term, the genome, the dot-bracket fold. Much of
the catalog carries it, and it is what makes these chemistries more than graphs.

## Three things that surprise people

### `status` changes what the network means

!!! warning "Comparing across statuses produces nonsense"
    Counting reactions in a `complete` network and an `observed` one and
    comparing the numbers is comparing a definition with a sample.

| status | meaning | what a missing reaction tells you |
|---|---|---|
| `complete` | the whole defined network, or a finished closure | it is not in the chemistry |
| `truncated` | a closure stopped by a size budget | nothing — raise the budget |
| `observed` | the reactions that fired in one run, each with a `count` | nothing about the chemistry; only about that run |

### Constructive chemistries grow their own species set

Many chemistries are **constructive**: S is open, and new molecules appear as
reactions produce them. For these the object of interest is a *closure* — apply
the rule until nothing new appears — and `truncated` is a routine answer rather
than a failure. `chemart.expand.expand` computes closures; see
[Comparing and analysing](guide/analysing.md). A closure is everything the
rule can make; a run of the gas is what a finite population does make.

The others have a fixed, enumerable species set, and their networks are simply
written down.

### Claimed capability is not computed capability

Every network reports what it contains:

```python
net.provides     # e.g. ['catalysts', 'initial-state', 'rate-constants', ...]
```

These tags are *derived from content* on every access. The catalog entry also
has a `provides` list, but that is a **claim about the chemistry**. The contract
guarantees computed ⊆ claimed, so an entry may advertise something its default
parameters don't exercise — `rbn` claims `mass-conservation`, but its default
network declares no conservation laws.

Rule of thumb: `net.provides` for what you hold; the catalog for what the
chemistry can do at other settings.

## Provenance is part of the data

Reproducing published models honestly means admitting where the sources ran
out. Every entry carries a `fidelity`:

| | meaning |
|---|---|
| `book` | implemented exactly as the book specifies |
| `book+decisions` | the book left gaps; each filled choice is in `decisions` |
| `reconstructed` | built from the original papers, cited in `sources` |

`reconstructed` is usually the *stronger* label, because a primary paper is more
precise than a survey chapter. What matters is the `decisions` list underneath
it, which records which sources could not be obtained, which rules were
inferred, and which published numbers are reproduced — and which are not.

Where a published result could not be reproduced, entries report that rather
than tuning parameters until it appeared. See
[Fidelity and trust](trust.md).

## Where the pieces live

```
chemart/
  api.py         list, describe, generate, evolve; the LLM tools
  catalog.py     loader, validator, index generator
  network.py     the Network / Species / Reaction record
  trajectory.py  the Trajectory / Frame record
  simulate.py    rate equations (ode) and Gillespie (ssa); rates and states (assign)
  measures/      the registry of measures, one module per section of the measures page
  kinetics.py    rate-law vocabulary, and k -> c conversion
  expand.py      closure of a constructive rule
  soup.py        well-stirred soup (stir) and the tally of what fired
  contract.py    the checks every chemistry passes
  cli.py         the `chemart` command
  chemistries/   one module per chemistry: generate(p, rng) and/or evolve(p, rng)
catalog/
  chemistries/   one YAML entry per chemistry — the specification
  explainers/    the prose of each chemistry's documentation page
  TYPES.md       the type of each entry, with the reason
hub/             the Chemart Hub and the simulation pit (chemart-hub)
```
