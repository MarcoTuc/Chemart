"""Random catalytic reaction networks (book eqs. 7.31-7.33). Catalog id: random-catalytic-networks.

For every unordered reactant pair {X_i, X_j} and product X_k, the reaction
X_i + X_j -> X_i + X_j + X_k exists with probability `density`, with a random
rate alpha. The dilution flux phi of eq. 7.24 is outflow = constant-total.
"""

from itertools import combinations_with_replacement

from chemart.helpers.explicit import network, term
from chemart.network import CONSTANT_TOTAL


def generate(p, rng):
    X = [f"X{i + 1}" for i in range(p.n)]
    reactions = []
    for i, j in combinations_with_replacement(range(p.n), 2):
        lhs = term(2, X[i]) if i == j else f"{X[i]} + {X[j]}"
        for k in range(p.n):
            if not p.allow_direct_replication and k in (i, j):
                continue
            if rng.random() >= p.density:
                continue
            alpha = rng.random() if p.alpha_distribution == "uniform" else rng.exponential()
            reactions.append((f"{lhs} -> {lhs} + {X[k]}", float(alpha)))
    return network(reactions, species=X, outflow=CONSTANT_TOTAL)
