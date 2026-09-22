"""B. Stoichiometric structure: linear algebra on S = P - R and the complex graph.

The complex-graph measures follow chemical reaction network theory (Feinberg,
Foundations of Chemical Reaction Network Theory, 2019): a complex is one side
of a reaction as a multiset, and each reaction is an arrow between complexes.
"""

from __future__ import annotations

import numpy as np

from chemart.measures import register


def complex_graph(net):
    """The digraph of complexes (frozen multisets of species), one arrow per reaction."""
    import networkx as nx

    g = nx.DiGraph()
    for r in net.reactions:
        g.add_edge(frozenset(r.reactants.items()), frozenset(r.products.items()))
    return g


@register("stoichiometric_rank", "B", needs="S")
def stoichiometric_rank(ctx) -> int:
    """Rank of the stoichiometric matrix: the dimension of the space the
    concentrations can move in."""
    return ctx.rank


@register("rank_ratio", "B", needs="S")
def rank_ratio(ctx) -> float:
    """rank(S) / n: the share of independent directions of change."""
    return ctx.rank / len(ctx.ids)


@register("conservation_laws", "B", needs="S")
def conservation_laws(ctx) -> int:
    """n - rank(S): the number of independent conserved linear combinations
    of amounts (dimension of the left null space of S)."""
    return len(ctx.ids) - ctx.rank


@register("conservative", "B", needs="S")
def conservative(ctx) -> bool:
    """Whether a strictly positive mass is conserved: some m > 0 with mᵀS = 0
    (a linear programme)."""
    from scipy.optimize import linprog

    n = len(ctx.ids)
    res = linprog(np.zeros(n), A_eq=ctx.S.T.astype(float), b_eq=np.zeros(ctx.S.shape[1]),
                  bounds=[(1.0, None)] * n, method="highs")
    return bool(res.status == 0)


@register("deficiency", "B", needs="S")
def deficiency(ctx) -> int:
    """Feinberg's deficiency δ = n_c - ℓ - rank(S): complexes minus linkage
    classes minus the rank. With δ = 0 and weak reversibility, mass action has
    exactly one positive steady state in each stoichiometric compatibility
    class, and it is stable, whatever the rates."""
    import networkx as nx

    g = complex_graph(ctx.net)
    return g.number_of_nodes() - nx.number_weakly_connected_components(g) - ctx.rank


@register("weakly_reversible", "B", needs="S")
def weakly_reversible(ctx) -> bool:
    """Whether every reaction lies on a cycle of the complex graph (each
    linkage class is strongly connected)."""
    import networkx as nx

    g = complex_graph(ctx.net)
    return nx.number_strongly_connected_components(g) == nx.number_weakly_connected_components(g)


@register("flux_dimension", "B", needs="S")
def flux_dimension(ctx) -> int:
    """r - rank(S): the number of independent steady-state flux patterns of
    the closed network (dimension of the right null space of S)."""
    return len(ctx.net.reactions) - ctx.rank
