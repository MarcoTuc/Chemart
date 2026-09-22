"""I. Information and algorithmic complexity."""

from __future__ import annotations

import zlib
from collections import Counter

import numpy as np

from chemart.measures import register


def canonical_text(ctx) -> str:
    """The reactions as text with species renamed by their degree profile
    (reactant uses, product uses; ties in order of first appearance), so the
    chemistry's naming scheme does not change the result. One sorted line per
    reaction."""
    uses = Counter()
    for r in ctx.net.reactions:
        for s, n in r.reactants.items():
            uses[s, "in"] += n
        for s, n in r.products.items():
            uses[s, "out"] += n
    order = sorted(ctx.ids, key=lambda s: (-uses[s, "in"], -uses[s, "out"], ctx.index[s]))
    name = {s: str(i) for i, s in enumerate(order)}

    def side(d):
        return "+".join(sorted(f"{n}{name[s]}" if n > 1 else name[s] for s, n in d.items()))

    return "\n".join(sorted(f"{side(r.reactants)}>{side(r.products)}" for r in ctx.net.reactions))


@register("compressibility", "I")
def compressibility(ctx) -> float:
    """Compressed size over raw size of the canonical reaction list (zlib,
    level 9): lower means more regular."""
    raw = canonical_text(ctx).encode()
    return len(zlib.compress(raw, 9)) / len(raw)


@register("degree_entropy", "I")
def degree_entropy(ctx) -> float:
    """Shannon entropy (nats) of the species degree distribution."""
    from chemart.measures.graph import species_degrees

    counts = np.array(list(Counter(species_degrees(ctx)).values()), dtype=float)
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum())


@register("structure_function_mi", "I", needs="str", cost="moderate", limit=20000)
def structure_function_mi(ctx, bins: int = 4) -> float | None:
    """Mutual information (nats) between the total structure length of a
    reaction's reactants and that of its products, each cut into `bins` (4)
    quantile classes: whether what comes out is predictable from what goes in.
    A crude estimate from lengths alone; None with fewer than 20 reactions."""
    size = {s.id: len(s.structure) if s.structure else len(s.id) for s in ctx.net.species}
    pairs = [(sum(size[s] * n for s, n in r.reactants.items()), sum(size[s] * n for s, n in r.products.items()))
             for r in ctx.net.reactions]
    if len(pairs) < 20:
        return None
    x, y = (np.array(v, dtype=float) for v in zip(*pairs))

    def classes(v):
        edges = np.unique(np.quantile(v, np.linspace(0, 1, bins + 1)[1:-1]))
        return np.searchsorted(edges, v, side="right")

    cx, cy = classes(x), classes(y)
    joint = Counter(zip(cx, cy))
    px, py, n = Counter(cx), Counter(cy), len(cx)
    return float(sum(c / n * np.log((c / n) / ((px[a] / n) * (py[b] / n))) for (a, b), c in joint.items()))
