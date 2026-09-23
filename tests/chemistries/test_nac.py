"""NAC reproduces book 11.3.3 and Suzuki's JSAI 2004 slides (2004.06_JSAI.presentation.pdf).

Slides: the local rewiring rule step by step (pp. 9-13), its effect on the clustering
coefficient and the average path length at N = 200, K = 10 (p. 14), the four edge types
(p. 17) and the hydrophilic/hydrophobic cancellation (pp. 18-19). Book 11.3.3 and the
abstract of Suzuki (2008): hydrophilic nodes cluster into a cell, hydrophobic nodes are
repelled (figure 11.14(a)-(b)).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import evolve, evolve_frames, generate_network
from chemart.chemistries.nac import (
    CANCELLED, HYDROPHILIC, MOVED, NO_DISTANCE_TWO, canonical, cluster_id, clustering,
    components, distance_two, mean_path_length, mixed_edges, random_graph, rewire,
)

ID = "nac"


def edge_set(adj) -> set[tuple[int, int]]:
    return {(u, v) for u in range(len(adj)) for v in adj[u] if u < v}


def path(labels: str):
    """The labelled path 0-1-2-... as (labels, adjacency)."""
    n = len(labels)
    adj = [set() for _ in range(n)]
    for v in range(n - 1):
        adj[v].add(v + 1)
        adj[v + 1].add(v)
    return tuple(labels), adj


# --- the rewiring rule (slides 9-13) ------------------------------------------------
def test_distance_two_is_measured_from_the_starting_node():
    _, adj = path("iiii")
    assert distance_two(adj, 0) == [2] and distance_two(adj, 1) == [3]
    assert distance_two(adj, 3) == [1]


def test_rule_moves_one_edge_from_b_to_c():
    # On the path 0-1-2 the only move from A = 0 is 0-1 -> 0-2 (and mirrored from A = 2);
    # node 1 has nothing at distance two, so its step is inert.
    kinds = set()
    for seed in range(40):
        labels, adj = path("iii")
        kind, a, b, c = rewire(labels, adj, np.random.default_rng(seed))
        kinds.add(kind)
        if kind == MOVED:
            assert (a, b, c) in {(0, 1, 2), (2, 1, 0)}
            assert edge_set(adj) == ({(0, 2), (1, 2)} if a == 0 else {(0, 1), (0, 2)})
        else:
            assert kind == NO_DISTANCE_TWO and a == 1
            assert edge_set(adj) == {(0, 1), (1, 2)}
        assert len(edge_set(adj)) == 2, "a rewiring deletes one edge and creates one"
    assert kinds == {MOVED, NO_DISTANCE_TWO}


def test_rewiring_conserves_the_edges_and_the_degree_of_the_starting_node():
    rng = np.random.default_rng(3)
    labels, adj = random_graph(30, 0.5, 4, rng)
    for _ in range(500):
        before = len(edge_set(adj))
        degrees = [len(a) for a in adj]
        kind, a, b, c = rewire(labels, adj, rng)
        assert len(edge_set(adj)) == before
        if kind == MOVED:
            assert len(adj[a]) == degrees[a]                  # A loses B and gains C
            assert len(adj[b]) == degrees[b] - 1 and len(adj[c]) == degrees[c] + 1
            assert labels[c] == labels[a]                     # slides 18-19


# --- the hydrophilic/hydrophobic constraint (slides 18-19) ---------------------------
def test_a_mixed_pair_cannot_be_wired():
    # i-i-o: from A = 0 the only node at distance two is the hydrophobic 2, and from
    # A = 2 it is the hydrophilic 0, so with the constraint nothing ever moves.
    kinds = set()
    for seed in range(40):
        labels, adj = path("iio")
        kind, *_ = rewire(labels, adj, np.random.default_rng(seed), polarity_constraint=True)
        kinds.add(kind)
        assert edge_set(adj) == {(0, 1), (1, 2)}, "the rewiring is cancelled"
    assert CANCELLED in kinds

    moved = 0
    for seed in range(40):
        labels, adj = path("iio")
        kind, *_ = rewire(labels, adj, np.random.default_rng(seed), polarity_constraint=False)
        moved += kind == MOVED
    assert moved > 0, "without the constraint the same steps are carried out"


# --- the published small-world measurement (slide 14) --------------------------------
def test_local_rewiring_makes_a_small_world():
    """Slide 14, N = 200, K = 10: the local rewiring increases C while keeping L constant."""
    rng = np.random.default_rng(0)
    labels, adj = random_graph(200, 1.0, 10, rng)              # one polarity: constraint inert
    c_random, l_random = clustering(adj), mean_path_length(adj)
    assert c_random == pytest.approx(10 / 200, abs=0.02)       # a random network has C ~ K/N
    for _ in range(2000):
        rewire(labels, adj, rng)
    assert len(edge_set(adj)) == 1000, "the rewiring conserves the number of edges"
    assert clustering(adj) > 2.5 * c_random                    # C rises well above random
    assert mean_path_length(adj) == pytest.approx(l_random, rel=0.15)   # L stays put


# --- species identity ----------------------------------------------------------------
def test_cluster_ids_are_invariant_under_relabelling():
    rng = np.random.default_rng(5)
    for _ in range(60):
        n = int(rng.integers(3, 9))
        labels = tuple("io"[int(rng.integers(2))] for _ in range(n))
        pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
        chosen = rng.choice(len(pairs), size=int(rng.integers(1, len(pairs) + 1)), replace=False)
        edges = tuple(sorted(pairs[int(i)] for i in chosen))
        perm = [int(x) for x in rng.permutation(n)]
        other_labels = tuple(labels[perm.index(i)] for i in range(n))
        other_edges = tuple(sorted((min(perm[u], perm[v]), max(perm[u], perm[v])) for u, v in edges))
        assert cluster_id(labels, edges) == cluster_id(other_labels, other_edges)
        assert canonical(labels, edges) == canonical(other_labels, other_edges)


def test_different_clusters_get_different_ids():
    assert cluster_id(("i", "o"), ((0, 1),)) != cluster_id(("i", "i"), ((0, 1),))
    assert cluster_id(("i",), ()) == "i"
    # a path and a star of three nodes are the same graph; a triangle is not
    assert cluster_id(("i",) * 3, ((0, 1), (1, 2))) == cluster_id(("i",) * 3, ((0, 2), (1, 2)))
    assert cluster_id(("i",) * 3, ((0, 1), (1, 2))) != cluster_id(("i",) * 3, ((0, 1), (0, 2), (1, 2)))


# --- the default network --------------------------------------------------------------
def nodes_per_species(net) -> dict[str, int]:
    hydrophilic, hydrophobic, _ = net.extras["conservation"]
    return {s: hydrophilic["vector"][s] + hydrophobic["vector"][s] for s in hydrophilic["vector"]}


def test_default_network_is_an_observed_rewiring_run():
    traj = evolve(ID, seed=1)
    net = traj.network
    assert traj.clock == "steps" and traj.times()[:3] == [0.0, 12.0, 24.0] and traj.times()[-1] == 600.0
    assert net.status == "observed" and net.reactions
    assert all(r.rate is None and r.count >= 1 for r in net.reactions)
    for r in net.reactions:
        assert sum(r.reactants.values()) == 1, "one rewiring acts on one cluster"
        assert 1 <= sum(r.products.values()) <= 2, "it rewires that cluster or splits it in two"
    size = nodes_per_species(net)
    n_nodes = len(net.extras["space"]["nodes"])
    assert sum(size[s] * n for s, n in net.initial_state.items()) == n_nodes
    assert sum(size[s] * n for s, n in net.extras["final_state"].items()) == n_nodes
    assert net.extras["space"]["polarity_constraint"] is True
    assert sum(net.extras["analysis"]["attempts"].values()) == net.params["steps"]


def test_every_reaction_conserves_nodes_and_edges():
    net = evolve(ID, seed=2).network
    ids, R, P = net.matrices()
    stoichiometry = (P - R).toarray()
    for law in net.extras["conservation"]:
        vector = np.array([law["vector"][s] for s in ids])
        assert not (vector @ stoichiometry).any(), f"{law['name']} is not conserved"


def test_hydrophilic_and_hydrophobic_nodes_demix():
    """Book fig. 11.14(a)-(b) and Suzuki (2008): mixed edges vanish, hydrophilic nodes cluster."""
    traj = evolve(ID, seed=1)
    mixed = traj.series("mixed_edges")
    assert mixed[0] > 0 and mixed[-1] == 0
    assert all(a >= b for a, b in zip(mixed, mixed[1:])), "a mixed edge can be deleted but never created"
    hydrophilic = traj.network.extras["space"]["nodes"].count(HYDROPHILIC)
    assert traj.series("largest_hydrophilic_cluster")[-1] >= hydrophilic - 1
    clustering_ = traj.series("clustering")
    assert clustering_[-1] > clustering_[0]

    loose = evolve(ID, seed=1, polarity_constraint=False)
    assert loose.series("mixed_edges")[-1] > 0, "the constraint is what demixes"


def test_given_graph_is_used_as_the_initial_state():
    graph = dict(polarities="iioo", edges=[[0, 1], [1, 2], [2, 3]])
    whole = cluster_id(tuple("iioo"), ((0, 1), (1, 2), (2, 3)))

    # steps=0 rewires without end, so the graph is read off its first frame
    run = evolve_frames(ID, seed=1, steps=0, **graph)
    first = next(run)
    run.close()
    assert first.state == {whole: 1.0} and not first.fired

    net = evolve(ID, seed=1, steps=1, **graph).network
    assert net.extras["space"]["nodes"] == list("iioo")
    assert net.extras["space"]["initial_edges"] == [[0, 1], [1, 2], [2, 3]]
    assert net.initial_state == {whole: 1.0}
    assert net.status == "observed"


# --- the closure -----------------------------------------------------------------------
def test_closure_terminates_and_conserves():
    net = generate_network(ID, seed=1, n_nodes=8, mean_degree=1.5)
    assert net.status == "complete" and net.species
    ids, R, P = net.matrices()
    stoichiometry = (P - R).toarray()
    for law in net.extras["conservation"]:
        vector = np.array([law["vector"][s] for s in ids])
        assert not (vector @ stoichiometry).any()
    assert all(sum(r.reactants.values()) == 1 for r in net.reactions)
    # clusters only fragment, so no product ever has more nodes than its reactant
    size = nodes_per_species(net)
    for r in net.reactions:
        (reactant,) = r.reactants
        assert sum(size[s] * n for s, n in r.products.items()) == size[reactant]


def test_closure_truncates_on_its_budget():
    net = generate_network(ID, seed=1, n_nodes=8, mean_degree=2.0, max_species=40)
    assert net.status == "truncated" and len(net.species) <= 40


def test_same_seed_same_graph():
    a = evolve(ID, seed=11)
    b = evolve(ID, seed=11)
    assert a.to_dict() == b.to_dict()
    assert evolve(ID, seed=12).network.extras["space"] != a.network.extras["space"]


# --- parameters -------------------------------------------------------------------------
def test_bad_parameters():
    with pytest.raises(ValueError, match="polarities"):
        generate_network(ID, polarities="ixo")
    with pytest.raises(ValueError, match="node ids"):
        generate_network(ID, polarities="iio", edges=[[0, 5]])
    with pytest.raises(ValueError, match="self-loop"):
        generate_network(ID, edges=[[1, 1]])
    with pytest.raises(ValueError, match="n_nodes"):
        generate_network(ID, n_nodes=1)
