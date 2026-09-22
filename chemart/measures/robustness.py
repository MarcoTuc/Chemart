"""H. Robustness and redundancy: how many ways a network has of doing the same thing."""

from __future__ import annotations

import numpy as np

from chemart.measures import register


@register("production_multiplicity", "H", needs="S")
def production_multiplicity(ctx) -> float:
    """Mean number of reactions with a net production of each species."""
    return float((ctx.S > 0).sum(axis=1).mean())


def _giant_share(g, species: list) -> float:
    import networkx as nx

    if g.number_of_nodes() == 0:
        return 0.0
    best = max((c for c in nx.connected_components(g)), key=lambda c: sum(1 for n in c if n[0] == "s"))
    return sum(1 for n in best if n[0] == "s") / len(species)


def _percolation_curve(g, order, species) -> float:
    g = g.copy()
    shares = [_giant_share(g, species)]
    for node in order:
        g.remove_node(node)
        shares.append(_giant_share(g, species))
    # area under the curve of giant-component share against fraction removed
    return float(np.trapezoid(shares, dx=1.0 / len(order)))


@register("percolation", "H", cost="moderate", limit=3000)
def percolation(ctx) -> dict[str, float]:
    """Area under the curve of the share of species in the largest connected
    piece as species are removed, at random and highest-degree first (Albert,
    Jeong & Barabási 2000). A robust network keeps a large area under both."""
    g = ctx.graph.to_undirected()
    species = [n for n in g.nodes if n[0] == "s"]
    rng = ctx.rng()
    random_order = [species[i] for i in rng.permutation(len(species))]
    targeted = sorted(species, key=lambda n: (-g.degree(n), str(n)))
    return {"random": _percolation_curve(g, random_order, species),
            "targeted": _percolation_curve(g, targeted, species)}
