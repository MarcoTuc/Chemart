# Comparing and analysing

The reason 98 unrelated models share one record is so you can ask questions
across them.

## Structural profiles

```python
import chemart, numpy as np

def profile(cid, seed=1):
    net = chemart.generate_network(cid, seed=seed)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    rank = int(np.linalg.matrix_rank(S)) if S.size else 0
    return dict(
        id=cid,
        species=len(ids),
        reactions=len(net.reactions),
        status=net.status,
        rank=rank,
        conservation_bound=len(ids) - rank,      # upper bound on independent laws
        catalysed=sum(1 for r in net.reactions if r.catalysts),
    )

for cid in ["dimerization", "brusselator", "oregonator", "matrix-chemistry"]:
    print(profile(cid))
```

`rank(S)` and the deficiency of a CRN come straight off the stoichiometric
matrix, and `len(species) − rank(S)` bounds the number of independent
conservation laws — a quick sanity check against the laws an entry declares.

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
organisation theory — are the natural next step from a closure. Several entries
compute them into `extras["analysis"]`; `matrix-chemistry` is the worked
example.

## Observed runs

For chemistries that are simulated rather than enumerated:

```python
import numpy as np
from chemart.soup import soup

fired, final_population = soup(
    divide, list(range(2, 40)), steps=2000, rng=np.random.default_rng(0),
    arity=2, dilution="constant",
)
# fired: [(reactants, products, count), ...] in order of first firing
```

`dilution="constant"` removes random molecules after each reaction to hold the
population at its starting size — the book's flow reactor. `alternatives=True`
draws one outcome per collision, where `expand` would record them all.

!!! note "Set `alternatives` deliberately"
    Molecules only have to be hashable, so tuples are legal molecules. A list of
    alternatives passed without the flag is indistinguishable from a right-hand
    side of tuple-valued molecules, and gets injected into the population as-is.

## Catalysis and autocatalysis

```python
catalysed = [r for r in net.reactions if r.catalysts]

autocatalytic = [r for r in net.reactions
                 if any(r.products.get(s, 0) > n for s, n in r.reactants.items())]
```

`r.catalysts` gives the species on both sides with the multiplicity that
survives. A reaction is autocatalytic in the loose sense when a reactant comes
out amplified. For the rigorous notions — reflexively autocatalytic and
food-generated (RAF) sets — use the `raf` entry, which implements the
Hordijk–Steel algorithm and can run on a system you supply rather than only on
its own.

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
substrate graph. Several entries already report graph statistics of their own
networks in `extras["analysis"]`, which is worth checking before recomputing
them.
