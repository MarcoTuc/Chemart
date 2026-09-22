"""Random Boolean networks (book 18.4.2). Catalog id: rbn.

Kauffman's classic RBN. N nodes, node i reads k_i distinct inputs (k_i = K,
or drawn with mean K) through a random truth table whose outputs are 1 with
probability function_bias. The network is the Boolean network written as
reactions: node i has the species x{i}_0 and x{i}_1, and truth-table row r of
node i is the reaction in which the input species of row r catalyse the
switch of x{i} to the row's output. With one molecule per node this is exactly
the node-by-node (asynchronous) semantics; the synchronous map and its
attractors are in extras.analysis.

The RBN machinery here (random_rbn, successor_map, attractors) is also what
the atoms of RBN World are built from: see chemart.chemistries.rbn_world.
"""

from __future__ import annotations

import numpy as np

from chemart.network import Network, Reaction, Species


def generate(p, rng):
    if p.K > p.N:
        raise ValueError(f"K must satisfy 1 <= K <= N, got K = {p.K} with N = {p.N}")
    return _classic(p, rng)


# ============================================================================
# Classic RBN
# ============================================================================
def power_law_exponent(K: float, N: int) -> float:
    """gamma such that P(k) ~ k^-gamma on k = 1..N has mean K (bisection)."""
    ks = np.arange(1, N + 1, dtype=float)

    def mean(g):
        w = ks ** -g
        return float((ks * w).sum() / w.sum())

    lo, hi = 0.0, 1.0
    while mean(hi) > K:
        hi *= 2
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mean(mid) > K else (lo, mid)
    return (lo + hi) / 2


def in_degrees(N: int, K: int, distribution: str, rng) -> tuple[list[int], float | None]:
    if distribution == "constant":
        return [K] * N, None
    if distribution == "poisson":
        return [int(min(N, 1 + rng.poisson(K - 1))) for _ in range(N)], None
    if K == 1:
        return [1] * N, None              # mean 1 on k >= 1 forces every k = 1 (gamma -> infinity)
    if K >= (N + 1) / 2:
        raise ValueError(f"a power law on 1..N has mean below (N + 1) / 2 = {(N + 1) / 2}; got K = {K}")
    gamma = power_law_exponent(K, N)
    ks = np.arange(1, N + 1)
    w = ks ** -gamma
    return [int(k) for k in rng.choice(ks, size=N, p=w / w.sum())], gamma


def random_rbn(N: int, degrees: list[int], bias: float, rng) -> tuple[list[list[int]], list[list[int]]]:
    """Distinct inputs per node (possibly itself) and truth tables; row r reads the inputs as bits, first input highest."""
    inputs = [[int(j) for j in rng.choice(N, size=k, replace=False)] for k in degrees]
    tables = [[int(b) for b in rng.random(2 ** k) < bias] for k in degrees]
    return inputs, tables


def successor_map(N: int, inputs, tables) -> np.ndarray:
    """Synchronous update of all 2^N states; bit i of a state is node i."""
    states = np.arange(2 ** N, dtype=np.int64)
    nxt = np.zeros_like(states)
    for i, (ins, table) in enumerate(zip(inputs, tables)):
        row = np.zeros_like(states)
        for j in ins:
            row = (row << 1) | ((states >> j) & 1)
        nxt |= np.asarray(table, dtype=np.int64)[row] << i
    return nxt


def attractors(f: np.ndarray) -> list[dict]:
    """Cycles of a functional graph with their basin sizes, largest basin first."""
    g = f.copy()
    n = max(1, int(len(f)).bit_length())
    for _ in range(n):                      # g = f^(2^n) maps every state onto its attractor
        g = g[g]
    label = np.full(len(f), -1, dtype=np.int64)
    cycles = []
    for s in np.unique(g):
        s = int(s)
        if label[s] >= 0:
            continue
        cycle = [s]
        t = int(f[s])
        while t != s:
            cycle.append(t)
            t = int(f[t])
        label[cycle] = len(cycles)
        cycles.append(cycle)
    basins = np.bincount(label[g], minlength=len(cycles))
    out = [{"length": len(c), "basin_size": int(b), "cycle": c} for c, b in zip(cycles, basins)]
    return sorted(out, key=lambda a: (-a["basin_size"], a["cycle"][0]))


def average_sensitivity(inputs, tables) -> float:
    """Mean over nodes of the summed probability that flipping one input flips the output (Derrida's lambda)."""
    total = 0.0
    for ins, table in zip(inputs, tables):
        k = len(ins)
        for t in range(k):
            mask = 1 << (k - 1 - t)
            total += sum(table[r] != table[r ^ mask] for r in range(2 ** k)) / 2 ** k
    return total / len(inputs)


def bits(state: int, n: int) -> str:
    """A state as a bitstring, node 0 first."""
    return "".join(str((state >> i) & 1) for i in range(n))


def _classic(p, rng):
    N = p.N
    degrees, gamma = in_degrees(N, p.K, p.K_distribution, rng)
    inputs, tables = random_rbn(N, degrees, p.function_bias, rng)
    initial = [int(b) for b in rng.integers(0, 2, size=N)]

    X = [[f"x{i}_0", f"x{i}_1"] for i in range(N)]
    reactions = []
    for i, (ins, table) in enumerate(zip(inputs, tables)):
        k = len(ins)
        for r in range(2 ** k):
            row = {j: (r >> (k - 1 - t)) & 1 for t, j in enumerate(ins)}
            out = table[r]
            if row.get(i, 1 - out) == out:        # node i already shows the output in this row
                continue
            lhs = {X[j][b]: 1 for j, b in row.items()}
            lhs[X[i][1 - out]] = 1
            rhs = {X[j][b]: 1 for j, b in row.items() if j != i}
            rhs[X[i][out]] = 1
            reactions.append(Reaction(lhs, rhs))

    f = successor_map(N, inputs, tables)
    found = attractors(f)
    x0 = sum(b << i for i, b in enumerate(initial))
    listed = [{"length": a["length"], "basin_size": a["basin_size"],
               "cycle": [bits(s, N) for s in a["cycle"][:64]]} for a in found]
    trajectory, seen, s = [], set(), x0
    while s not in seen:
        seen.add(s)
        trajectory.append(s)
        s = int(f[s])
    reached = next(k for k, a in enumerate(found) if s in a["cycle"])
    analysis = {
        "inputs": inputs,
        "functions": ["".join(map(str, t)) for t in tables],
        "in_degrees": degrees,
        "state_encoding": "bitstrings list node 0 first; function row r reads the node's inputs as bits, first input highest",
        "average_sensitivity": average_sensitivity(inputs, tables),
        "annealed_sensitivity": 2 * float(np.mean(degrees)) * p.function_bias * (1 - p.function_bias),
        "attractors": listed,
        "initial_attractor": reached,
        "transient_length": trajectory.index(s),
    }
    if gamma is not None:
        analysis["power_law_exponent"] = gamma
    return Network(
        species=[Species(x) for pair in X for x in pair],
        reactions=reactions,
        initial_state={X[i][b]: 1.0 for i, b in enumerate(initial)},
        extras={"analysis": analysis},
    )
