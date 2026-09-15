"""Jain-Krishna autocatalytic set model (book 15.2.2). Catalog id: jain-krishna.

Fast dynamics (eq. 15.7): x_i' = sum_j c_ij x_j - x_i sum_kj c_kj x_j, i.e. each
link c_ij = 1 is the catalytic reaction X_j -> X_j + X_i under a constant-total
dilution. Slow dynamics: repeatedly remove a least-populated species at the
fast attractor and rewire the replacement node with link probability p.
"""

import numpy as np

from chemart.helpers.explicit import network
from chemart.network import CONSTANT_TOTAL


def attractor(C: np.ndarray) -> tuple[float, np.ndarray]:
    """Perron-Frobenius eigenvalue and the normalised attractor of eq. 15.7."""
    m = len(C)
    values, vectors = np.linalg.eig(C)
    k = int(np.argmax(values.real))
    lam = float(values.real[k])
    if lam > 1e-9:
        x = np.abs(vectors[:, k].real)
    else:
        # No autocatalytic set: iterate (C + I), which concentrates on the ends of the longest chains.
        x, A = np.ones(m), C + np.eye(m)
        for _ in range(m):
            x = A @ x
            x /= x.sum()
    return lam, x / x.sum()


def _rewire(C, r, rng, p, self_loops):
    C[r, :] = rng.random(len(C)) < p
    C[:, r] = rng.random(len(C)) < p
    if not self_loops:
        C[r, r] = 0


def generate(p, rng):
    m = p.m
    C = (rng.random((m, m)) < p.p).astype(int)
    if not p.self_loops:
        np.fill_diagonal(C, 0)
    for _ in range(p.graph_updates):
        _, x = attractor(C)
        least = np.flatnonzero(x <= x.min() + 1e-12)
        _rewire(C, int(rng.choice(least)), rng, p.p, p.self_loops)
    lam, x = attractor(C)

    X = [f"X{i + 1}" for i in range(m)]
    reactions = [(f"{X[j]} -> {X[j]} + {X[i]}", 1.0) for i, j in zip(*np.nonzero(C))]
    return network(
        reactions,
        species=X,
        initial_state={s: float(v) for s, v in zip(X, x)},
        outflow=CONSTANT_TOTAL,
        extras={"analysis": {
            "graph_updates": p.graph_updates,
            "perron_frobenius_eigenvalue": lam,
            "attractor_support": [s for s, v in zip(X, x) if v > 1e-9],
        }},
    )
