"""Replicator equation as bimolecular reactions (book eqs. 7.22-7.24). Catalog id: replicator-equation."""

from chemart.helpers.explicit import network
from chemart.helpers.params import square_matrix, vector
from chemart.network import CONSTANT_TOTAL


def generate(p, rng):
    A = square_matrix("A", p.A)
    n = len(A)
    x0 = vector("x0", p.x0, n, minimum=0.0) if p.x0 else [1.0 / n] * n
    X = [f"X{i + 1}" for i in range(n)]
    reactions = []
    for i in range(n):
        for j in range(n):
            a = A[i][j]
            if a == 0:
                continue  # elastic collision
            if i == j:
                text = f"2 {X[i]} -> 3 {X[i]}" if a > 0 else f"2 {X[i]} -> {X[i]}"
            else:
                text = f"{X[i]} + {X[j]} -> 2 {X[i]} + {X[j]}" if a > 0 else f"{X[i]} + {X[j]} -> {X[j]}"
            reactions.append((text, abs(a)))
    return network(reactions, species=X, initial_state=dict(zip(X, x0)), outflow=CONSTANT_TOTAL)
