# Comparing and analysing

The reason unrelated models share one record is so you can ask questions
across them. `chemart.measure` answers most of them; this page is about the
record underneath it, and about the questions it does not answer.
[Measures for comparison](measures.md) is the full list.

## Structural profiles

```python
import chemart
from chemart import measures

for cid in ["michaelis-menten", "brusselator", "oregonator", "matrix-chemistry"]:
    net = chemart.generate_network(cid, seed=1)
    print(cid, net.status, chemart.measure(net, ["n_species", "n_reactions",
                                                 "stoichiometric_rank", "conservation_laws",
                                                 "catalysed_fraction"]))
```

Leave `names` out for every cheap measure that applies, and pass
`cost="moderate"` or `"exponential"` for the dearer ones. For a generator, a
profile is a curve rather than a number:

```python
rows = measures.sweep("random-catalytic-networks", {"n": [10, 20, 40]}, seeds=range(5))
measures.scaling(rows, "n_reactions")            # how it grows with size
```

The stoichiometric matrix itself is one call away when you want to compute
something the registry does not have:

```python
ids, R, P = net.matrices()                       # sparse; S = P - R
```

!!! warning "Only compare within a status"
    A `complete` network is a definition; an `observed` one is a sample from a
    single run. Their reaction counts are not the same kind of number.

## Conservation laws

31 entries claim `mass-conservation`, and 30 declare explicit laws at default
parameters — `rbn` claims the capability without declaring any, which is the
claimed-versus-computed gap in miniature.

A law lives in `extras["conservation"]` and always carries a `vector`. Read the
payload defensively, because two details bite:

```python
import numpy as np

ids, R, P = net.matrices()
S = (P - R).toarray().astype(int)

for i, law in enumerate(net.extras.get("conservation", [])):
    v = law["vector"]
    # `vector` is a dict keyed by species for most laws, a list for others
    m = np.array([v.get(s, 0) for s in ids]) if isinstance(v, dict) else np.array(v)
    residual = S.T @ m
    # a law with a `modulus` is MODULAR: chameleon's hold mod 3, so a plain
    # zero test would flag correct laws as violated
    modulus = law.get("modulus")
    holds = (residual % modulus == 0).all() if modulus else not residual.any()
    assert holds, law.get("name", f"law {i}")       # `name` is not guaranteed
```

Across the library that reader verifies all 115 declared laws. If you write a
new chemistry that conserves atoms, declaring the law and checking `Sᵀm = 0` in
a test is the single most effective way to catch stoichiometry bugs — several
entries found real errors that way.

## Closures of constructive chemistries

51 entries are constructive: the species set grows as reactions produce new
molecules, so the object of interest is the closure under the reaction rule.

```python
from chemart.expand import expand

def divide(a, b):
    """Prime-number chemistry: a + b -> a + b/a when a divides b."""
    return (a, b // a) if a < b and b % a == 0 else None

species, reactions, status = expand(divide, seed=[12, 2, 3])
# status "complete" = closed; "truncated" = hit max_species
```

`expand(react, seed, arity=2, max_species=1000, ordered=True, alternatives=False)`

- `react(*molecules)` returns the **complete** right-hand side — surviving
  reactants included, so catalysts must be listed — or `None` for an elastic
  collision.
- `ordered=True` when `react(a, b)` and `react(b, a)` can differ, which is the
  usual case for operator/operand chemistries.
- `alternatives=True` when one set of reactants can react several ways; `react`
  then returns an iterable of right-hand sides, each becoming its own reaction.

**Organisations** — closed and self-maintaining sets, in the sense of chemical
organisation theory — are the natural next step from a closure:
`chemart.measure(net, ["organisations"], cost="exponential")` counts them and
gives the largest. Several entries also compute their own into
`extras["analysis"]`; `matrix-chemistry` is the worked example.

## Observed runs

A chemistry whose process is the point — a Turing gas, or a lattice — is run
rather than enumerated, and what you measure changes as it runs:

```python
traj = chemart.evolve("prime-number-chemistry", seed=1)
traj.network                                     # what fired, with counts
measures.over(traj, ["richness", "shannon", "n_reactions"], window=5)
```

[Evolving a chemistry](evolving.md) covers the frames, the clocks and the
sampling caveat: an observed network is what one run happened to visit, so it
must be compared with another observed network of the same window size, not
with a closure.

To run a soup of your own molecules, `chemart.soup.stir` is the same
well-stirred loop the chemistries use, as a generator of frames; `soup()` runs
it to the end and returns what fired.

## Catalysis and autocatalysis

```python
catalysed = [r for r in net.reactions if r.catalysts]

autocatalytic = [r for r in net.reactions
                 if any(r.products.get(s, 0) > n for s, n in r.reactants.items())]
```

`r.catalysts` gives the species on both sides with the multiplicity that
survives. A reaction is autocatalytic in the loose sense when a reactant comes
out amplified. The rigorous notions are measures of any network:

```python
net = chemart.generate_network("raf", seed=1)
chemart.measure(net, ["max_raf_fraction", "irreducible_rafs"], cost="exponential")
# {'max_raf_fraction': 0.764…, 'irreducible_rafs': 20}
```

`max_raf_fraction` is the Hordijk–Steel maxRAF as a share of the reactions and
`irreducible_rafs` counts the smallest RAFs inside it (up to 20).
`autocatalytic_cores` finds Blokhuis cores, on networks of at most 400 species
and reactions — `applicable` says so when a network is too big for it. The RAF
measures start from the food set:
`extras["food"]`, or, without it, the species with inflow, the buffered ones,
or those in the initial state. Pass `food=` to choose it yourself.

## As a graph

```python
import networkx as nx

g = nx.DiGraph()
for i, r in enumerate(net.reactions):
    rid = f"r{i}"
    for s, n in r.reactants.items():
        g.add_edge(s, rid, stoich=n)
    for s, n in r.products.items():
        g.add_edge(rid, s, stoich=n)
```

That is the standard bipartite species/reaction graph; project it for the
substrate graph. `chemart.measures.Context(net).graph` builds the same graph,
cached, and `species_graph` the projection — the graph measures (degrees,
assortativity, modularity, motifs, spectra) read them. Several entries also
report graph statistics of their own networks in `extras["analysis"]`, which is
worth checking before recomputing them.
