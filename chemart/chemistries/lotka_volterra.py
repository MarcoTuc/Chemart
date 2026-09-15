"""Generalised Lotka-Volterra, x_i' = x_i (r_i + sum_j b_ij x_j). Catalog id: lotka-volterra.

Each term becomes the mass-action reaction with the same ODE contribution.
A predator-prey pair (b_ij = -b_ji > 0) becomes the single reaction
X_i + X_j -> 2 X_i of book eq. 7.26 rather than two separate reactions.
"""

from chemart.helpers.explicit import network
from chemart.helpers.params import square_matrix, vector


def _interaction(X, i, j, b):
    if b > 0:
        return [(f"{X[i]} + {X[j]} -> 2 {X[i]} + {X[j]}", b)]
    if b < 0:
        return [(f"{X[i]} + {X[j]} -> {X[j]}", -b)]
    return []


def generate(p, rng):
    r = vector("r", p.r)
    n = len(r)
    if n == 0:
        raise ValueError("r must list at least one growth rate")
    B = square_matrix("B", p.B, n)
    x0 = vector("x0", p.x0, n, minimum=0.0)
    X = [f"X{i + 1}" for i in range(n)]

    reactions = []
    for i in range(n):
        if r[i] > 0:
            reactions.append((f"{X[i]} -> 2 {X[i]}", r[i]))
        elif r[i] < 0:
            reactions.append((f"{X[i]} -> ", -r[i]))
    for i in range(n):
        if B[i][i] > 0:
            reactions.append((f"2 {X[i]} -> 3 {X[i]}", B[i][i]))
        elif B[i][i] < 0:
            reactions.append((f"2 {X[i]} -> {X[i]}", -B[i][i]))
    for i in range(n):
        for j in range(i + 1, n):
            bij, bji = B[i][j], B[j][i]
            if bij > 0 and bji == -bij:
                reactions.append((f"{X[i]} + {X[j]} -> 2 {X[i]}", bij))
            elif bji > 0 and bij == -bji:
                reactions.append((f"{X[j]} + {X[i]} -> 2 {X[j]}", bji))
            else:
                reactions += _interaction(X, i, j, bij) + _interaction(X, j, i, bji)
    return network(reactions, species=X, initial_state=dict(zip(X, x0)))
