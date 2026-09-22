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


# --------------------------------------------------------------------------
# Moderate and exponential: extreme rays and flux analysis
# --------------------------------------------------------------------------

def extreme_rays(A: np.ndarray) -> np.ndarray:
    """The extreme rays of the pointed cone {x >= 0, A x = 0}, one per row
    (double description method, pycddlib)."""
    import cdd

    m, n = A.shape
    rows = [[0.0, *map(float, A[i])] for i in range(m)] + \
           [[0.0, *[1.0 if k == j else 0.0 for k in range(n)]] for j in range(n)]
    mat = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY, lin_set=set(range(m)))
    gen = cdd.copy_generators(cdd.polyhedron_from_matrix(mat))
    rays = [row[1:] for row in gen.array if row[0] == 0 and any(abs(v) > 1e-9 for v in row[1:])]
    return np.array(rays, dtype=float).reshape(len(rays), n)


@register("p_invariants", "B", needs="S", cost="exponential", limit=300)
def p_invariants(ctx) -> int:
    """Number of minimal semi-positive conservation laws (P-invariants,
    conserved moieties): the extreme rays of {m >= 0, mᵀS = 0}."""
    return len(extreme_rays(ctx.S.T.astype(float)))


def exchange_matrix(ctx) -> tuple[np.ndarray, int]:
    """S with boundary reactions: an input for each food species and an output
    for each species no reaction consumes. Returns (S_ext, number of internal reactions)."""
    consumed = {s for r in ctx.net.reactions for s, n in r.reactants.items() if n > r.products.get(s, 0)}
    columns = [ctx.S.astype(float)]
    for s in ctx.food:
        e = np.zeros((len(ctx.ids), 1))
        e[ctx.index[s]] = 1.0
        columns.append(e)
    for s in ctx.ids:
        if s not in consumed:
            e = np.zeros((len(ctx.ids), 1))
            e[ctx.index[s]] = -1.0
            columns.append(e)
    return np.hstack(columns), ctx.S.shape[1]


@register("elementary_flux_modes", "B", needs="S F", cost="exponential", limit=150)
def elementary_flux_modes(ctx) -> dict[str, float]:
    """Number of elementary flux modes (minimal steady-state pathways) of the
    network fed from its food set and drained where nothing consumes, and their
    mean number of internal reactions (Schuster, Fell & Dandekar 2000)."""
    S_ext, r = exchange_matrix(ctx)
    rays = extreme_rays(S_ext)
    lengths = [int((np.abs(ray[:r]) > 1e-9).sum()) for ray in rays]
    return {"count": len(rays), "mean_length": float(np.mean(lengths)) if lengths else 0.0}


@register("blocked_fraction", "B", needs="S F", cost="moderate", limit=3000)
def blocked_fraction(ctx) -> float:
    """Share of reactions that can carry no flux at steady state, with the food
    boundary of `elementary_flux_modes` (flux variability analysis, done as one
    linear programme: maximise the reactions that carry any flux)."""
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix, hstack, identity

    S_ext, r = exchange_matrix(ctx)
    m = S_ext.shape[1]
    # variables v (m fluxes, unbounded above) and z (r indicators, z_j <= min(v_j, 1));
    # steady states form a cone, so any flux that is live somewhere can be scaled to 1
    A_eq = hstack([csr_matrix(S_ext), csr_matrix((S_ext.shape[0], r))])
    A_ub = hstack([-identity(m, format="csr")[:r], identity(r, format="csr")])
    c = np.concatenate([np.zeros(m), -np.ones(r)])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(r), A_eq=A_eq, b_eq=np.zeros(S_ext.shape[0]),
                  bounds=[(0, None)] * m + [(0, 1)] * r, method="highs")
    if res.status != 0:
        return 1.0
    return float((res.x[m:] < 0.5).sum()) / r
