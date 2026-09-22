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


def _scope_size(net, food, without: set[int]) -> int:
    from chemart.measures.organisation import expansion
    from chemart.network import Network

    kept = Network(species=net.species, reactions=[r for j, r in enumerate(net.reactions) if j not in without])
    return len(expansion(kept, food)[-1])


@register("knockout_tolerance", "H", needs="S F", cost="moderate", limit=1500)
def knockout_tolerance(ctx) -> float:
    """Share of reactions whose removal leaves the scope of the food set
    unchanged: how much of the network is backed up by alternatives."""
    full = _scope_size(ctx.net, ctx.food, set())
    return sum(_scope_size(ctx.net, ctx.food, {j}) == full for j in range(len(ctx.net.reactions))) / len(ctx.net.reactions)


@register("synthetic_lethal_pairs", "H", needs="S F", cost="moderate", limit=1500)
def synthetic_lethal_pairs(ctx, samples: int = 400) -> float | None:
    """Among pairs of reactions that are each dispensable on their own, the
    share whose joint removal shrinks the scope: reactions that back each other
    up. Sampled (400 pairs) on large networks; None with fewer than two
    dispensable reactions."""
    full = _scope_size(ctx.net, ctx.food, set())
    spare = [j for j in range(len(ctx.net.reactions)) if _scope_size(ctx.net, ctx.food, {j}) == full]
    if len(spare) < 2:
        return None
    pairs = [(a, b) for i, a in enumerate(spare) for b in spare[i + 1:]]
    if len(pairs) > samples:
        rng = ctx.rng()
        pairs = [pairs[k] for k in rng.choice(len(pairs), size=samples, replace=False)]
    return sum(_scope_size(ctx.net, ctx.food, {a, b}) < full for a, b in pairs) / len(pairs)


@register("degeneracy", "H", needs="S F", cost="moderate", limit=3000)
def degeneracy(ctx, samples: int = 20) -> float | None:
    """Mean number of reaction-disjoint routes from the food set to a species
    of the scope (up to 20 species sampled): structurally different ways of
    making the same thing (Edelman & Gally 2001). None when the food makes
    nothing new."""
    import networkx as nx

    from chemart.measures.organisation import expansion

    scope = expansion(ctx.net, ctx.food)[-1]
    targets = sorted(scope - set(ctx.food))
    if not targets:
        return None
    rng = ctx.rng()
    if len(targets) > samples:
        targets = [targets[k] for k in rng.choice(len(targets), size=samples, replace=False)]
    g = nx.DiGraph()
    for j, r in enumerate(ctx.net.reactions):
        g.add_edge(("r_in", j), ("r_out", j), capacity=1)         # a reaction carries one route
        for s in r.reactants:
            g.add_edge(("s", s), ("r_in", j), capacity=len(ctx.net.reactions))
        for s in r.products:
            g.add_edge(("r_out", j), ("s", s), capacity=len(ctx.net.reactions))
    for s in ctx.food:
        g.add_edge("food", ("s", s), capacity=len(ctx.net.reactions))
    routes = [nx.maximum_flow_value(g, "food", ("s", t)) for t in targets]
    return float(np.mean(routes))
