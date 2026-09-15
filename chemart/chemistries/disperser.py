"""Chemical disperser load balancing, C_ij + X_i -> C_ij + X_j (book eq. 17.3). Catalog id: disperser."""

from chemart.helpers.explicit import network
from chemart.helpers.params import edges


def generate(p, rng):
    E = edges("graph", p.graph)
    nodes = sorted({v for e in E for v in e})
    links = sorted({(a, b) for a, b in E} | {(b, a) for a, b in E})
    for key, value in p.initial_jobs.items():
        if key not in {str(v) for v in nodes}:
            raise ValueError(f"initial_jobs key {key!r} is not a node of the graph {nodes}")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f"initial_jobs[{key!r}] must be a non-negative number, got {value!r}")

    X = {v: f"X{v}" for v in nodes}
    C = {link: f"C{link[0]}_{link[1]}" for link in links}
    reactions = [(f"{C[(a, b)]} + {X[a]} -> {C[(a, b)]} + {X[b]}", p.k) for a, b in links]
    initial = {X[v]: float(p.initial_jobs.get(str(v), 0.0)) for v in nodes}
    initial.update({c: p.catalyst_concentration for c in C.values()})
    compartments = {str(v): [X[v]] + [C[(a, b)] for a, b in links if a == v] for v in nodes}
    return network(
        reactions,
        species=list(X.values()) + list(C.values()),
        initial_state=initial,
        extras={"compartments": compartments},
    )
