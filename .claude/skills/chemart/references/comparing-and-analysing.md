# Comparing and analysing networks

The reason different models share one record is so you can ask questions
across them. `chemart.measures` implements the questions: dozens of registered
measures, from stoichiometric rank to autocatalytic cores. This file covers
using them, and the analyses the record supports directly.

## Contents
- [Surveying the catalog](#surveying-the-catalog)
- [Measures](#measures)
- [Measures across arguments, seeds and null models](#measures-across-arguments-seeds-and-null-models)
- [Conservation laws](#conservation-laws)
- [Closures and constructive chemistries](#closures-and-constructive-chemistries)
- [Observed runs](#observed-runs)
- [Catalysis and autocatalysis](#catalysis-and-autocatalysis)
- [Exporting to graph tools](#exporting-to-graph-tools)

## Surveying the catalog

```bash
uv run python .claude/skills/chemart/scripts/survey.py --provides rate-constants
uv run python .claude/skills/chemart/scripts/survey.py --constructive --verbose
uv run python .claude/skills/chemart/scripts/survey.py --family origin-of-life
uv run python .claude/skills/chemart/scripts/survey.py --generate --limit 12
```

`--generate` actually builds each network with default parameters and reports
real sizes and computed capabilities rather than catalog claims. It is slower
(seconds per chemistry) but it is the honest version, and it is the quickest
way to find, say, every chemistry that yields under 200 reactions by default.

## Measures

```python
import chemart
from chemart import measures

net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
chemart.measure(net)                                   # every cheap measure that applies
chemart.measure(net, ["deficiency", "max_raf_fraction", "conservation_laws"])
chemart.measure(net, cost="moderate")                  # add the dearer tiers
measures.applicable(net)                               # {name: None, or why it does not apply}
measures.describe()                                    # name, section, input, needs, cost, meaning
```

```bash
uv run chemart measure kauffman-autocatalytic-sets --seed 1 --cost moderate --why
uv run chemart measure run.json          # a saved network, or a trajectory from `chemart evolve`
```

What to know before reading the numbers:

- **Sections.** A size, B stoichiometry (rank, conservation laws, deficiency,
  P-invariants, flux modes), C graph (degree statistics, bow-tie, spectra,
  modularity, motifs), D organisation (maxRAF, scope, autocatalytic cores,
  chemical organisations), E growth, F kinetics (stability, steady states,
  oscillation), G dynamics (diversity, turnover, attractor type), H
  robustness (knockouts, degeneracy), I information (compressibility).
- **Inputs.** A measure takes a network, a population state (one frame's
  `{species: amount}`) or a whole trajectory. `chemart.measure(traj)` gives the
  trajectory measures plus the network measures of what fired.
- **Needs.** Kinetic measures need a rate on every reaction; organisation
  measures need a food set (`extras["food"]`, else the inflow, buffered or
  initial species). A measure that does not apply is **left out**, never
  filled with a default: absence is information, and `applicable` says why.
- **Cost.** `cheap` runs by default; `moderate` (LPs, many shortest paths,
  steady-state searches) and `exponential` (extreme rays, organisations,
  cores) on request. Dear measures have a node limit and are skipped above it
  unless `force=True`.
- **Status still matters.** A measure of an `observed` network describes one
  run, not the chemistry, and a `truncated` closure is a lower bound.

## Measures across arguments, seeds and null models

```python
rows = measures.sweep("raf", {"n": [4, 5, 6]}, seeds=range(3),
                      names=["n_species", "n_reactions", "max_raf_fraction"])
measures.scaling(rows, "n_reactions")          # log-log exponent against n_species
measures.zscores(net, ["clustering", "reciprocity"], samples=20)   # against the null model
```

A **generator** is studied across its arguments (`sweep`, one row per grid
point and seed). A graph measure means little on its own: compare it with its
null model (`zscores`: reactant and product slots swapped between reactions,
keeping every species' degree and every reaction's arity) before calling it a
property of the chemistry. A **gas** is studied in evolutionary time with
`measures.over(traj, names, window=w)`; see `simulating-dynamics.md`.

The same questions by hand, from `S`, are a few lines if you need a variant:
`ids, R, P = net.matrices()`, `S = (P - R)`, and `len(ids) − rank(S)` bounds the
number of independent conservation laws — a quick sanity check against the
laws an entry declares in `extras["conservation"]`.

## Conservation laws

31 entries claim `mass-conservation` and 30 actually declare laws at default
parameters — `rbn` claims it without declaring any, which is the claim-versus-
computed gap in miniature. 115 laws in total. A law always carries `vector`,
and usually `name`; `modulus`, `unique` and `notes` also occur.

Copy the canonical reader from `network-record.md` rather than writing your
own, because two details bite: `vector` is a dict keyed by species for 87 of
the laws but a **list** aligned with `ids` for the other 28, and a law carrying
a `modulus` is **modular** — `chameleon`'s hold mod 3, so a plain zero test
flags correct laws as violated. With that reader, all 115 verify.

If you are writing a new chemistry that conserves atoms, declaring the law and
checking `Sᵀm = 0` in a test is the single most effective way to catch
stoichiometry bugs — several entries found real errors that way.

## Closures and constructive chemistries

51 entries are constructive: the species set grows as reactions produce new
molecules. The object of interest is the closure under the reaction rule,
computed by `chemart.expand.expand`:

```python
from chemart.expand import expand

def divide(a, b):
    """Prime-number chemistry: a + b -> a + b/a when a divides b."""
    return (a, b // a) if a < b and b % a == 0 else None

species, reactions, status = expand(divide, seed=[12, 2, 3])
# status "complete" = closed; "truncated" = hit max_species
```

Signature: `expand(react, seed, arity=2, max_species=1000, ordered=True,
alternatives=False)`.

- `react(*molecules)` returns the **complete** right-hand side (surviving
  reactants included, so catalysts must be listed) or `None` for an elastic
  collision.
- `ordered=True` when `react(a, b)` and `react(b, a)` can differ — the usual
  case for operator/operand chemistries.
- `alternatives=True` when one set of reactants can react several ways; `react`
  then returns an iterable of right-hand sides, each becoming its own reaction.

A `truncated` status means the budget stopped it, not that the chemistry is
finite. Raise `max_species` (or the entry's own size parameter) to see more.

**Organisations** (closed and self-maintaining sets, in the sense of chemical
organisation theory) are a natural next step from a closure:
`chemart.measure(net, ["organisations"])` counts them (an exponential measure,
limited to small networks), and `matrix-chemistry` also reports its own in
`extras["analysis"]`.

## Observed runs

To run a catalogued gas, use `chemart.evolve` (see `simulating-dynamics.md`).
For a rule of your own, `chemart.soup.soup` runs a well-stirred population and
records what fired:

```python
import numpy as np
from chemart.soup import soup

def divide(a, b):
    """Same rule as above, now run as a well-stirred population."""
    return (a, b // a) if a < b and b % a == 0 else None

fired, final_population = soup(
    divide, list(range(2, 40)), steps=2000, rng=np.random.default_rng(0),
    arity=2, dilution="constant",
)
# fired: [(reactants, products, count), ...] in order of first firing
busiest = max(fired, key=lambda f: f[2])
```

`dilution="constant"` removes random molecules after each reaction to hold the
population at its initial size — the book's flow reactor. `alternatives=True`
draws one outcome per collision (where `expand` would record them all).
`chemart.soup.stir` is the same loop as a generator of frames, which is what an
evolve face is built from (`extending-chemart.md`).

Set `alternatives` deliberately: molecules only have to be hashable, so tuples
are legal molecules, and a list of alternatives passed without the flag is
indistinguishable from a right-hand side of tuple-valued molecules and gets
injected into the population as-is.

## Catalysis and autocatalysis

```python
catalysed = [r for r in net.reactions if r.catalysts]
autocatalytic = [r for r in net.reactions
                 if any(r.products.get(s, 0) > n for s, n in r.reactants.items())]
```

`r.catalysts` gives species on both sides with the multiplicity that survives.
A reaction is autocatalytic in the loose sense when a reactant comes out with
higher multiplicity than it went in. For the rigorous notions use the
measures: `max_raf_fraction` and `irreducible_rafs` (Hordijk–Steel RAF sets,
on any network with catalysts and a food set) and `autocatalytic_cores`
(Blokhuis et al.'s stoichiometric cores). The `raf` entry generates the
Hordijk–Steel model itself.

## Exporting to graph tools

The record is deliberately plain, so conversion is mechanical:

```python
import networkx as nx
g = nx.DiGraph()
for i, r in enumerate(net.reactions):
    rid = f"r{i}"
    for s, n in r.reactants.items(): g.add_edge(s, rid, stoich=n)
    for s, n in r.products.items():  g.add_edge(rid, s, stoich=n)
```

That is the standard bipartite species/reaction graph. For the substrate graph
(species connected when they appear in a common reaction), project it. The
graph measures (section C) already work on this graph, so reach for them
before recomputing statistics yourself.
