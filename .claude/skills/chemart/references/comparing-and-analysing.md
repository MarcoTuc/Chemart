# Comparing and analysing networks

The reason different models share one record is so you can ask questions
across them. This file covers the analyses that the record supports directly.

## Contents
- [Surveying the catalog](#surveying-the-catalog)
- [Structural comparison](#structural-comparison)
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

## Structural comparison

```python
import chemart, numpy as np

def profile(cid, seed=1):
    net = chemart.generate_network(cid, seed=seed)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    return {
        "id": cid,
        "species": len(ids),
        "reactions": len(net.reactions),
        "status": net.status,
        "rank": int(np.linalg.matrix_rank(S)) if S.size else 0,
        "conservation_laws": len(ids) - (int(np.linalg.matrix_rank(S)) if S.size else 0),
        "catalysed": sum(1 for r in net.reactions if r.catalysts),
        "provides": net.provides,
    }
```

The **deficiency** of a CRN and its rank are computable from `S` directly, and
`len(species) − rank(S)` bounds the number of independent conservation laws —
a quick sanity check against the laws an entry declares in
`extras["conservation"]`.

Beware comparing raw counts across `status` values; see `network-record.md`.

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
organisation theory) are a natural next step from a closure, and several
entries compute them in `extras["analysis"]` — `matrix-chemistry` is the
worked example.

## Observed runs

For chemistries that are simulated rather than enumerated, use
`chemart.soup.soup`:

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
higher multiplicity than it went in. For the rigorous notions — reflexively
autocatalytic and food-generated sets — use the `raf` entry, which implements
the Hordijk-Steel RAF algorithm and can run on a supplied system, not only on
its own generated one.

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
(species connected when they appear in a common reaction), project it. Several
entries report graph statistics of their own networks in `extras["analysis"]`,
which is worth reading before recomputing them.
