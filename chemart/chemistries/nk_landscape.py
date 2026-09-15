"""Kauffman's NK fitness landscape (book 18.4.1), as a replicator-mutator network. Catalog id: nk-landscape.

Each gene i reads its own state and those of K epistatic partners; the
concatenated bits (partners first, own state last) index a random table
M(i, j) in [0, 1); genome fitness is the mean over genes (eq. 18.22). The
network replicates every genotype at its fitness with single point
mutations, under constant-total dilution; the landscape itself is in
extras.analysis.
"""

from itertools import product

import numpy as np

from chemart.helpers.explicit import network


def generate(p, rng):
    N, K = p.N, p.K
    if K >= N:
        raise ValueError(f"K must satisfy 0 <= K <= N - 1, got K = {K} with N = {N}")
    if p.topology == "adjacent":
        partners = [[(i + d) % N for d in range(1, K + 1)] for i in range(N)]
    else:
        partners = [[int(j) for j in rng.choice([j for j in range(N) if j != i], size=K, replace=False)] for i in range(N)]
    table = rng.random((N, 2 ** (K + 1)))

    genomes = ["".join(bits) for bits in product("01", repeat=N)]

    def fitness(g: str) -> float:
        return float(np.mean([table[i, int("".join(g[k] for k in partners[i]) + g[i], 2)] for i in range(N)]))

    W = {g: fitness(g) for g in genomes}

    def flips(g: str):
        return [g[:i] + ("1" if g[i] == "0" else "0") + g[i + 1:] for i in range(N)]

    reactions = []
    for g in genomes:
        reactions.append((f"{g} -> 2 {g}", W[g] * (1 - p.mu) ** N))
        if p.mu > 0:
            reactions += [(f"{g} -> {g} + {h}", W[g] * p.mu * (1 - p.mu) ** (N - 1)) for h in flips(g)]
    return network(
        reactions,
        species=genomes,
        outflow="constant-total",
        extras={"analysis": {
            "fitness": W,
            "epistatic_partners": partners,
            "local_optima": [g for g in genomes if all(W[g] >= W[h] for h in flips(g))],
            "global_optimum": max(genomes, key=W.get),
        }},
    )
