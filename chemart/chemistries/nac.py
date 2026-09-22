"""NAC: Suzuki's Network Artificial Chemistry (book 11.3.3). Catalog id: nac.

A graph stands in for the solvent. Every node is a molecule, classified as
hydrophilic or hydrophobic; an edge between two nodes is a weak bond (the
hydrophilic or hydrophobic interaction that water mediates), so the graph is
at once the space and the set of non-covalent contacts. Molecules therefore
come in clusters: a species here is a connected cluster of nodes, taken up to
isomorphism, and a reaction is what one rewiring event does to a cluster.

The dynamics is Suzuki's local rewiring rule (JSAI 2004 slides 9-13, after
Davidsen, Ebel & Bornholdt's acquaintance networks), at every time step:

    (1) a starting node A is chosen at random;
    (2) a stopping node B is chosen at random from the nodes adjacent to A;
    (3) a new stopping node C is chosen at random from the nodes at distance
        two from A;
    (4) edge A-B is eliminated and edge A-C is created.

Repeating it raises the clustering coefficient while leaving the average path
length alone, i.e. it turns a random network into a small-world network
(slide 14: N = 200, K = 10).

On top of that rule NAC adds the polarity constraint (slides 17-19): a
hydrophilic node and a hydrophobic node cannot be wired by the passive
rewiring, so if C has the other polarity the rewiring is cancelled. Mixed
edges are therefore destroyed and never created: the network demixes, the
hydrophilic nodes gather into a cluster (the "network cell" of Suzuki 2008)
and the hydrophobic nodes into their own.

A rewiring conserves the nodes and the number of edges of the cluster it
acts on, so every reaction is either an isomerisation

    cluster -> the same cluster rewired

or a fragmentation

    cluster -> cluster + cluster

with the node counts of both polarities and the edge count conserved.

Two faces: ``generate`` is the set of clusters reachable from the initial
ones by single rewirings (chemart.expand.expand, alternatives=True);
``evolve`` runs the rule on one graph, yields a frame every n_nodes steps and
returns the reactions that fired with their counts.

The active layer of NAC (nodes holding von Neumann programs that rewire the
covalent and hydrogen edges: polymerase, helicase, splitase, centrosome) is
not implemented; the papers that specify those programs were not obtainable.
See the catalog decisions.
"""

from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from itertools import combinations

from chemart.expand import expand
from chemart.helpers import params as check
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame

HYDROPHILIC, HYDROPHOBIC = "i", "o"
POLARITIES = HYDROPHILIC + HYDROPHOBIC

# Outcomes of one attempted rewiring.
MOVED = "moved"
NO_NEIGHBOUR = "no-neighbour"
NO_DISTANCE_TWO = "no-distance-two"
CANCELLED = "cancelled-polarity"


# ---------------------------------------------------------------------------
# Clusters: canonical form up to isomorphism
# ---------------------------------------------------------------------------
def _refine(labels: tuple[str, ...], adj: tuple[frozenset[int], ...]) -> list[int]:
    """1-dimensional Weisfeiler-Leman colours: an isomorphism-invariant partition."""
    colours = [POLARITIES.index(x) for x in labels]
    n = len(labels)
    while True:
        signature = [(colours[v], tuple(sorted(colours[u] for u in adj[v]))) for v in range(n)]
        rank = {s: i for i, s in enumerate(sorted(set(signature)))}
        new = [rank[s] for s in signature]
        if new == colours:
            return colours
        colours = new


def _twin_keys(labels, adj, v: int) -> tuple:
    """Two nodes sharing either key are interchangeable: swapping them is an
    automorphism, because they have the same label and the same neighbours
    (the open key) or the same neighbours plus each other (the closed key)."""
    return (("open", labels[v], frozenset(adj[v] - {v})),
            ("closed", labels[v], frozenset(adj[v] | {v})))


@lru_cache(maxsize=65536)
def canonical(labels: tuple[str, ...], edges: tuple[tuple[int, int], ...]) -> tuple[tuple[str, ...], tuple[tuple[int, int], ...]]:
    """Canonical labels and edges of a labelled graph: isomorphic graphs give the same pair.

    The nodes are ordered by maximising, position after position, the pair
    (label, adjacency to the nodes already placed). The candidates at each
    position are the unplaced nodes of the least refined colour, minus
    duplicates of a twin already tried, both of which are isomorphism
    invariant, so the order that wins is a canonical one.
    """
    n = len(labels)
    adj = tuple(frozenset(v for u, v in edges if u == w) | frozenset(u for u, v in edges if v == w)
                for w in range(n))
    colours = _refine(labels, adj)
    twins = [_twin_keys(labels, adj, v) for v in range(n)]

    placed = [False] * n
    order: list[int] = []
    best: list[int] | None = None
    best_code: tuple | None = None

    def extend(code: tuple) -> None:
        nonlocal best, best_code
        if len(order) == n:
            if best_code is None or code > best_code:
                best, best_code = list(order), code
            return
        least = min(colours[v] for v in range(n) if not placed[v])
        entries, seen = [], set()
        for v in range(n):
            if placed[v] or colours[v] != least or seen.intersection(twins[v]):
                continue
            seen.update(twins[v])
            bits = 0
            for j, u in enumerate(order):
                if u in adj[v]:
                    bits |= 1 << (n - 1 - j)
            entries.append(((labels[v], bits), v))
        top = max(entry for entry, _ in entries)
        for entry, v in entries:
            if entry != top:
                continue
            placed[v] = True
            order.append(v)
            extend(code + (entry,))
            order.pop()
            placed[v] = False

    extend(())
    position = {v: i for i, v in enumerate(best or [])}
    canon_edges = tuple(sorted(
        (min(position[u], position[v]), max(position[u], position[v])) for u, v in edges))
    return tuple(labels[v] for v in (best or [])), canon_edges


def cluster_id(labels: tuple[str, ...], edges: tuple[tuple[int, int], ...]) -> str:
    """Species id: the canonical labels, then the canonical edges, e.g. ``iio:0-1.1-2``."""
    canon_labels, canon_edges = canonical(labels, edges)
    body = ".".join(f"{u}-{v}" for u, v in canon_edges)
    return "".join(canon_labels) + (f":{body}" if body else "")


def cluster_text(labels: tuple[str, ...], edges: tuple[tuple[int, int], ...]) -> str:
    """Readable structure: ``nodes=iio edges=0-1,1-2``, in canonical order."""
    canon_labels, canon_edges = canonical(labels, edges)
    return (f"nodes={''.join(canon_labels)} "
            f"edges={','.join(f'{u}-{v}' for u, v in canon_edges) or 'none'}")


def subgraph(labels, adj, nodes: list[int]) -> tuple[tuple[str, ...], tuple[tuple[int, int], ...]]:
    """The labelled graph induced on `nodes`, renumbered 0..k-1."""
    index = {v: i for i, v in enumerate(sorted(nodes))}
    edges = tuple(sorted(
        (index[u], index[v]) for u in index for v in adj[u] if u < v and v in index))
    return tuple(labels[v] for v in sorted(nodes)), edges


# ---------------------------------------------------------------------------
# The graph
# ---------------------------------------------------------------------------
def random_polarities(n: int, hydrophilic_fraction: float, rng) -> tuple[str, ...]:
    """n node labels, round(n * hydrophilic_fraction) of them hydrophilic, at random positions."""
    hydrophilic = {int(v) for v in rng.permutation(n)[:round(n * hydrophilic_fraction)]}
    return tuple(HYDROPHILIC if v in hydrophilic else HYDROPHOBIC for v in range(n))


def random_edges(n: int, mean_degree: float, rng) -> list[set[int]]:
    """Adjacency of round(n * mean_degree / 2) distinct undirected edges drawn uniformly."""
    adj: list[set[int]] = [set() for _ in range(n)]
    pairs = list(combinations(range(n), 2))
    m = min(round(n * mean_degree / 2), len(pairs))
    if m:
        for i in rng.choice(len(pairs), size=int(m), replace=False):
            u, v = pairs[int(i)]
            adj[u].add(v)
            adj[v].add(u)
    return adj


def random_graph(n: int, hydrophilic_fraction: float, mean_degree: float, rng):
    """A random NAC graph: polarities and wiring (the starting point of every run)."""
    return random_polarities(n, hydrophilic_fraction, rng), random_edges(n, mean_degree, rng)


def distance_two(adj, a: int) -> list[int]:
    """The nodes at distance exactly two from a."""
    return sorted({c for b in adj[a] for c in adj[b]} - adj[a] - {a})


def rewire(labels, adj, rng, polarity_constraint: bool = True):
    """One local rewiring step; returns (outcome, A, B, C) with B, C possibly None.

    Steps (1)-(4) of the rule, with the NAC constraint of slides 18-19: if the
    new stopping node C has the other polarity, the whole step is cancelled.
    """
    n = len(labels)
    a = int(rng.integers(n))
    if not adj[a]:
        return NO_NEIGHBOUR, a, None, None
    neighbours = sorted(adj[a])
    b = neighbours[int(rng.integers(len(neighbours)))]
    far = distance_two(adj, a)
    if not far:
        return NO_DISTANCE_TWO, a, b, None
    c = far[int(rng.integers(len(far)))]
    if polarity_constraint and labels[c] != labels[a]:
        return CANCELLED, a, b, c
    adj[a].discard(b)
    adj[b].discard(a)
    adj[a].add(c)
    adj[c].add(a)
    return MOVED, a, b, c


def component(adj, start: int) -> list[int]:
    seen, queue = {start}, deque([start])
    while queue:
        v = queue.popleft()
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                queue.append(u)
    return sorted(seen)


def components(adj) -> list[list[int]]:
    seen: set[int] = set()
    out = []
    for v in range(len(adj)):
        if v not in seen:
            part = component(adj, v)
            seen.update(part)
            out.append(part)
    return out


# ---------------------------------------------------------------------------
# Network statistics (the published small-world measurement)
# ---------------------------------------------------------------------------
def clustering(adj) -> float:
    """Watts-Strogatz clustering coefficient: the average over all nodes of the
    fraction of a node's neighbour pairs that are themselves linked (a node of
    degree < 2 counts as 0)."""
    total = 0.0
    for v in range(len(adj)):
        k = len(adj[v])
        if k < 2:
            continue
        links = sum(1 for x, y in combinations(sorted(adj[v]), 2) if y in adj[x])
        total += 2 * links / (k * (k - 1))
    return total / len(adj) if adj else 0.0


def mean_path_length(adj) -> float:
    """Average shortest-path length over the ordered pairs that are connected."""
    total, pairs = 0, 0
    for s in range(len(adj)):
        dist = {s: 0}
        queue = deque([s])
        while queue:
            v = queue.popleft()
            for u in adj[v]:
                if u not in dist:
                    dist[u] = dist[v] + 1
                    queue.append(u)
        total += sum(dist.values())
        pairs += len(dist) - 1
    return total / pairs if pairs else 0.0


def mixed_edges(labels, adj) -> int:
    """Edges between a hydrophilic and a hydrophobic node."""
    return sum(1 for u in range(len(adj)) for v in adj[u] if u < v and labels[u] != labels[v])


# ---------------------------------------------------------------------------
# Species bookkeeping
# ---------------------------------------------------------------------------
class Clusters:
    """The clusters seen so far, by canonical id."""

    def __init__(self) -> None:
        self.graphs: dict[str, tuple[tuple[str, ...], tuple[tuple[int, int], ...]]] = {}

    def add(self, labels, edges) -> str:
        sid = cluster_id(labels, edges)
        self.graphs.setdefault(sid, (labels, edges))
        return sid

    def of(self, adj, labels, nodes: list[int]) -> str:
        return self.add(*subgraph(labels, adj, nodes))

    def species(self, ids) -> list[Species]:
        return [Species(s, structure=cluster_text(*self.graphs[s])) for s in ids]

    def nodes_of(self, sid: str, polarity: str) -> int:
        return sum(1 for x in self.graphs[sid][0] if x == polarity)

    def edges_of(self, sid: str) -> int:
        return len(self.graphs[sid][1])

    # -- the reaction, for the closure --------------------------------------
    def outcomes(self, sid: str, polarity_constraint: bool) -> list[tuple[str, ...]] | None:
        """Every cluster (or pair of clusters) one rewiring can turn `sid` into."""
        labels, edges = self.graphs[sid]
        n = len(labels)
        adj = [set() for _ in range(n)]
        for u, v in edges:
            adj[u].add(v)
            adj[v].add(u)
        found: dict[tuple[str, ...], None] = {}
        for a in range(n):
            for b in sorted(adj[a]):
                for c in distance_two(adj, a):
                    if polarity_constraint and labels[c] != labels[a]:
                        continue
                    adj[a].discard(b)
                    adj[b].discard(a)
                    adj[a].add(c)
                    adj[c].add(a)
                    products = tuple(sorted(self.of(adj, labels, part) for part in components(adj)))
                    adj[a].discard(c)
                    adj[c].discard(a)
                    adj[a].add(b)
                    adj[b].add(a)
                    found[products] = None
        return list(found) or None


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------
def _initial_graph(p, rng):
    if p.polarities:
        if not isinstance(p.polarities, str) or set(p.polarities) - set(POLARITIES):
            raise ValueError(
                f"polarities must be a string over {POLARITIES!r} (i hydrophilic, o hydrophobic), "
                f"got {p.polarities!r}")
        labels = tuple(p.polarities)
        adj = random_edges(len(labels), p.mean_degree, rng)
    else:
        labels, adj = random_graph(p.n_nodes, p.hydrophilic_fraction, p.mean_degree, rng)
    if p.edges:
        adj = [set() for _ in labels]
        for u, v in check.edges("edges", p.edges):
            if not (0 <= u < len(labels) and 0 <= v < len(labels)):
                raise ValueError(f"edges: node ids must be in 0..{len(labels) - 1}, got [{u}, {v}]")
            adj[u].add(v)
            adj[v].add(u)
    return labels, adj


def _space(labels, adj, initial_edges, p) -> dict:
    return {
        "kind": "solvent graph: every node is a molecule, every edge a weak "
                "(hydrophilic or hydrophobic) bond",
        "nodes": list(labels),
        "initial_edges": [[u, v] for u, v in initial_edges],
        "edges": [[u, v] for u in range(len(adj)) for v in sorted(adj[u]) if u < v],
        "rewiring_rule": "(1) starting node A at random; (2) stopping node B at random among "
                         "the nodes adjacent to A; (3) new stopping node C at random among the "
                         "nodes at distance two from A; (4) edge A-B is eliminated and edge A-C "
                         "is created (Suzuki, JSAI 2004, slides 9-13)",
        "polarity_constraint": bool(p.polarity_constraint),
    }


def _conservation(clusters: Clusters, ids) -> list[dict]:
    return [
        {"name": "hydrophilic nodes",
         "vector": {s: clusters.nodes_of(s, HYDROPHILIC) for s in ids}},
        {"name": "hydrophobic nodes",
         "vector": {s: clusters.nodes_of(s, HYDROPHOBIC) for s in ids}},
        {"name": "weak edges",
         "vector": {s: clusters.edges_of(s) for s in ids}},
    ]


def _start(p, rng):
    labels, adj = _initial_graph(p, rng)
    initial_edges = [(u, v) for u in range(len(adj)) for v in sorted(adj[u]) if u < v]
    return labels, adj, initial_edges


def generate(p, rng):
    """The closure: every cluster reachable from the initial ones by single rewirings."""
    labels, adj, initial_edges = _start(p, rng)
    return _closure(p, labels, adj, initial_edges)


def _observe(labels, adj) -> dict:
    """The demixing and small-world measures of the whole graph."""
    parts = components(adj)
    hydrophilic = [len(part) for part in parts if all(labels[v] == HYDROPHILIC for v in part)]
    return {
        "mixed_edges": mixed_edges(labels, adj),
        "largest_hydrophilic_cluster": max(hydrophilic, default=0),
        "clustering": round(clustering(adj), 6),
        "path_length": round(mean_path_length(adj), 6),
    }


def _frame(step, labels, adj, clusters, tally) -> Frame:
    state = Counter(clusters.of(adj, labels, part) for part in components(adj))
    return Frame(t=float(step), state={s: float(n) for s, n in state.items()},
                 fired=tally.flush(), observables=_observe(labels, adj))


def evolve(p, rng):
    """The rewiring run on one graph: a frame every n_nodes attempted steps."""
    labels, adj, initial_edges = _start(p, rng)
    clusters = Clusters()
    start = Counter(clusters.of(adj, labels, part) for part in components(adj))
    tally = Tally()
    outcome = Counter()
    every = max(1, len(labels))
    yield _frame(0, labels, adj, clusters, tally)

    for step in range(p.steps):
        kind, a, b, c = rewire(labels, adj, rng, p.polarity_constraint)
        outcome[kind] += 1
        if kind == MOVED:
            # A rewiring touches one cluster: A, B and C were all connected before
            # the move, and afterwards the cluster is whole or split in two.
            after = [component(adj, a)]
            if b not in after[0]:
                after.append(component(adj, b))
            before = sorted(v for part in after for v in part)
            lhs = (clusters.add(*_before_subgraph(labels, adj, before, a, b, c)),)
            rhs = tuple(sorted(clusters.of(adj, labels, part) for part in after))
            if Counter(lhs) != Counter(rhs):
                tally.add(lhs, rhs)
        if (step + 1) % every == 0:
            yield _frame(step + 1, labels, adj, clusters, tally)

    if p.steps % every:
        yield _frame(p.steps, labels, adj, clusters, tally)
    fired = tally.reactions()
    final = Counter(clusters.of(adj, labels, part) for part in components(adj))
    ids = list(dict.fromkeys([*start, *(s for lhs, rhs, _ in fired for s in (*lhs, *rhs)), *final]))
    reactions = [Reaction.of(lhs, rhs, count=n) for lhs, rhs, n in fired]
    extras = {
        "space": _space(labels, adj, initial_edges, p),
        "conservation": _conservation(clusters, ids),
        "analysis": {"attempts": dict(outcome)},
        "final_state": {s: n for s, n in final.most_common()},
    }
    return Network(
        species=clusters.species(ids),
        reactions=reactions,
        status="observed",
        initial_state={s: float(n) for s, n in start.items()},
        extras=extras,
    )


def _before_subgraph(labels, adj, nodes, a, b, c):
    """The cluster `nodes` as it was before the move a-b -> a-c."""
    index = {v: i for i, v in enumerate(sorted(nodes))}
    edges = {(index[u], index[v]) for u in index for v in adj[u] if u < v and v in index}
    edges.discard((min(index[a], index[c]), max(index[a], index[c])))
    edges.add((min(index[a], index[b]), max(index[a], index[b])))
    return tuple(labels[v] for v in sorted(nodes)), tuple(sorted(edges))


def _closure(p, labels, adj, initial_edges):
    clusters = Clusters()
    seed = [clusters.of(adj, labels, part) for part in components(adj)]
    start = Counter(seed)

    def react(sid):
        return clusters.outcomes(sid, p.polarity_constraint)

    found, pairs, status = expand(react, list(dict.fromkeys(seed)), arity=1,
                                  max_species=p.max_species, ordered=True, alternatives=True)
    reactions = [Reaction.of(lhs, rhs) for lhs, rhs in pairs]
    extras = {
        "space": _space(labels, adj, initial_edges, p),
        "conservation": _conservation(clusters, found),
        "analysis": {
            "seed": list(dict.fromkeys(seed)),
            "cluster_sizes": {s: len(clusters.graphs[s][0]) for s in found},
        },
    }
    return Network(
        species=clusters.species(found),
        reactions=reactions,
        status=status,
        initial_state={s: float(n) for s, n in start.items()},
        extras=extras,
    )
