"""C. Graph topology, on the bipartite species-reaction graph.

Rule 1 of the measures page: one graph, kept for every measure. Species and
reactions are both nodes; an edge runs from each reactant to its reaction and
from each reaction to its products (a catalyst gets both). Projecting onto
species would turn every reaction into a clique that is not in the chemistry.
Two measures need the species view itself, and say so: reciprocity and the
bow-tie use the substrate -> product digraph.
"""

from __future__ import annotations

import numpy as np

from chemart.measures import register

#: Dense eigen-decompositions above this many nodes are left to force=True.
SPECTRAL_LIMIT = 3000


def species_degrees(ctx) -> np.ndarray:
    """In how many reactions each species takes part (either side, counted once)."""
    deg = np.zeros(len(ctx.ids))
    for r in ctx.net.reactions:
        for s in set(r.reactants) | set(r.products):
            deg[ctx.index[s]] += 1
    return deg


def gini(x: np.ndarray) -> float:
    """Gini coefficient: 0 when all values are equal, towards 1 when one holds everything."""
    x = np.sort(np.asarray(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = len(x)
    return float((2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum()))


@register("degree_cv", "C")
def degree_cv(ctx) -> float:
    """Coefficient of variation (sd / mean) of species degrees: how unevenly
    species take part in reactions."""
    deg = species_degrees(ctx)
    return float(deg.std() / deg.mean()) if deg.mean() else 0.0


@register("degree_gini", "C")
def degree_gini(ctx) -> float:
    """Gini coefficient of species degrees."""
    return gini(species_degrees(ctx))


@register("assortativity", "C")
def assortativity(ctx) -> float | None:
    """Degree assortativity of the undirected bipartite graph: whether
    high-degree species meet high-degree reactions. None when undefined
    (all degrees equal)."""
    import warnings

    import networkx as nx

    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        r = nx.degree_assortativity_coefficient(ctx.graph.to_undirected())
    return None if r is None or np.isnan(r) else float(r)


@register("clustering", "C")
def clustering(ctx) -> float:
    """Mean bipartite clustering of the species (Latapy et al. 2008): how much
    the reaction neighbourhoods of species that share a reaction overlap."""
    from networkx.algorithms import bipartite

    g = ctx.graph.to_undirected()
    return float(bipartite.average_clustering(g, [("s", s) for s in ctx.ids]))


@register("reciprocity", "C")
def reciprocity(ctx) -> float | None:
    """Share of substrate -> product links between species that also run back.
    None when there are no such links."""
    import networkx as nx

    g = ctx.species_graph
    return float(nx.reciprocity(g)) if g.number_of_edges() else None


@register("cycle_rank", "C")
def cycle_rank(ctx) -> int:
    """Independent loops of the undirected bipartite graph: edges - nodes + components."""
    import networkx as nx

    g = ctx.graph.to_undirected()
    return g.number_of_edges() - g.number_of_nodes() + nx.number_connected_components(g)


@register("flow_hierarchy", "C")
def flow_hierarchy(ctx) -> float:
    """Share of edges on no cycle (Luo & Magee 2011): 1 is purely
    feed-forward, 0 means everything feeds back."""
    import networkx as nx

    return float(nx.flow_hierarchy(ctx.graph))


@register("bow_tie", "C")
def bow_tie(ctx) -> dict[str, float]:
    """Shares of species in the core (the largest strongly connected set of
    the substrate -> product graph), upstream of it (in), downstream of it
    (out), and elsewhere (other)."""
    import networkx as nx

    g = ctx.species_graph
    core = max(nx.strongly_connected_components(g), key=len)
    probe = next(iter(core))
    upstream = nx.ancestors(g, probe) - core
    downstream = nx.descendants(g, probe) - core
    n = len(ctx.ids)
    return {"core": len(core) / n, "in": len(upstream) / n, "out": len(downstream) / n,
            "other": (n - len(core) - len(upstream) - len(downstream)) / n}


@register("efficiency", "C", cost="moderate", limit=5000)
def efficiency(ctx) -> float:
    """Global efficiency of the undirected bipartite graph: the mean inverse
    shortest-path length over all pairs of nodes."""
    import networkx as nx

    return float(nx.global_efficiency(ctx.graph.to_undirected()))


def _laplacian_spectrum(ctx) -> np.ndarray:
    """Eigenvalues of the normalised Laplacian, computed once per network."""
    cached = getattr(ctx, "_laplacian", None)
    if cached is None:
        import networkx as nx

        L = nx.normalized_laplacian_matrix(ctx.graph.to_undirected()).toarray()
        cached = ctx._laplacian = np.clip(np.linalg.eigvalsh(L), 0.0, None)
    return cached


@register("spectral_gap", "C", limit=20000)
def spectral_gap(ctx) -> float:
    """Second-smallest eigenvalue of the normalised Laplacian of the
    undirected bipartite graph (its algebraic connectivity): 0 when it falls
    apart, larger when it is well connected."""
    import networkx as nx

    g = ctx.graph.to_undirected()
    if g.number_of_nodes() < 3 or not nx.is_connected(g):
        return 0.0
    if g.number_of_nodes() <= 500:
        return float(_laplacian_spectrum(ctx)[1])
    return float(nx.algebraic_connectivity(g, normalized=True, method="tracemin_lu", seed=0))


@register("spectral_entropy", "C", cost="moderate", limit=SPECTRAL_LIMIT)
def spectral_entropy(ctx) -> float:
    """Shannon entropy (nats) of the normalised Laplacian eigenvalues, read as
    a distribution."""
    ev = _laplacian_spectrum(ctx)
    p = ev[ev > 0] / ev.sum()
    return float(-(p * np.log(p)).sum())


@register("catalytic_spectral_radius", "C", needs="C", limit=SPECTRAL_LIMIT)
def catalytic_spectral_radius(ctx) -> float:
    """Largest eigenvalue modulus of the catalytic graph, an arrow from each
    catalyst to each net product of the reactions it catalyses. At least 1
    exactly when the graph has a cycle: an autocatalytic set (Jain & Krishna 1998)."""
    n = len(ctx.ids)
    A = np.zeros((n, n))
    for r in ctx.net.reactions:
        for c in r.catalysts:
            for p, k in r.products.items():
                if k > r.reactants.get(p, 0):
                    A[ctx.index[p], ctx.index[c]] = 1.0
    return float(np.abs(np.linalg.eigvals(A)).max()) if n else 0.0


def _fold(X):
    """Identical rows of the incidence matrix, kept once with how many there are.

    Two rows of equal degree contribute nothing to NODF (neither is the poorer
    of the pair), so duplicates only have to be counted — and reaction sets
    repeat a lot, which is what makes the pair sum affordable on big networks.
    """
    groups: dict[bytes, list[int]] = {}
    for k, row in enumerate(X):
        groups.setdefault(row.tobytes(), []).append(k)
    keep = [rows[0] for rows in groups.values()]
    return X[keep], np.array([len(rows) for rows in groups.values()], dtype=float)


@register("nodf", "C")
def nodf(ctx) -> float:
    """Nestedness (NODF, Almeida-Neto et al. 2008) of the species x reaction
    incidence matrix, 0 to 100: whether rarer species take part only in
    reactions that commoner species also take part in."""
    M = np.zeros((len(ctx.ids), len(ctx.net.reactions)))
    for j, r in enumerate(ctx.net.reactions):
        for s in set(r.reactants) | set(r.products):
            M[ctx.index[s], j] = 1.0

    def paired(X, block: int = 256):
        # the pair sum in blocks of rows: the whole overlap matrix would be
        # 10^8 entries on the larger networks, and only the poorer row of each
        # pair contributes
        pairs = len(X) * (len(X) - 1) // 2
        X, counts = _fold(X)
        deg = X.sum(axis=1)
        order = np.argsort(-deg)
        X, counts, deg = X[order], counts[order], deg[order]
        share = np.where(deg > 0, 1.0 / np.where(deg > 0, deg, 1.0), 0.0)
        total = 0.0
        for start in range(0, len(X), block):
            stop = min(start + block, len(X))
            overlap = X[start:stop] @ X.T
            poorer = deg[None, :] < deg[start:stop, None]
            total += float((np.where(poorer, overlap * share[None, :], 0.0)
                            * counts[start:stop, None] * counts[None, :]).sum())
        return 100.0 * total, pairs

    rows, nr = paired(M)
    cols, nc = paired(M.T)
    return float((rows + cols) / (nr + nc)) if nr + nc else 0.0


@register("mean_hyperedge_size", "C", needs="S")
def mean_hyperedge_size(ctx) -> float:
    """Mean number of distinct species per reaction, reactions read as hyperedges."""
    return float(np.mean([len(set(r.reactants) | set(r.products)) for r in ctx.net.reactions]))


@register("modularity", "C", cost="moderate", limit=10000)
def modularity(ctx) -> float:
    """Modularity Q of the Louvain communities of the undirected bipartite
    graph: how well the network splits into semi-independent parts. Compare it
    with its null model (zscores) before reading it."""
    import networkx as nx

    g = ctx.graph.to_undirected()
    if g.number_of_edges() == 0:
        return 0.0
    parts = nx.community.louvain_communities(g, seed=ctx.seed)
    return float(nx.community.modularity(g, parts))


#: The 13 connected three-node patterns of a directed graph (triad census codes).
CONNECTED_TRIADS = ("021D", "021U", "021C", "111D", "111U", "030T", "030C",
                    "201", "120D", "120U", "120C", "210", "300")


@register("motif_profile", "C", cost="moderate", limit=3000)
def motif_profile(ctx, samples: int = 10) -> dict[str, float]:
    """Significance profile of the connected triads of the substrate -> product
    graph (Milo et al. 2004): each triad's z-score against networks randomised
    by `measures.randomize`, the vector normalised to length 1."""
    import networkx as nx

    from chemart.measures import Context, randomize

    def census(net):
        g = Context(net).species_graph
        g.remove_edges_from(nx.selfloop_edges(g))
        return nx.triadic_census(g)

    real = census(ctx.net)
    rng = ctx.rng()
    null = [census(randomize(ctx.net, rng)) for _ in range(samples)]
    z = {}
    for t in CONNECTED_TRIADS:
        values = np.array([c[t] for c in null], dtype=float)
        sd = values.std(ddof=1) if len(values) > 1 else 0.0
        z[t] = 0.0 if sd == 0 else (real[t] - values.mean()) / sd
    norm = np.sqrt(sum(v * v for v in z.values()))
    return {t: float(v / norm) if norm else 0.0 for t, v in z.items()}
