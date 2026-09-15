"""Bitstring immune system model of Farmer, Packard & Perelson (book 11.2.2). Catalog id: farmer-immune.

Antibodies X_i carry an epitope e and a paratope p; antigens Y_j carry an
epitope. Matching strengths (eq. 11.4) fix every rate constant of eqs.
11.5-11.6, each ODE term becoming one mass-action reaction.
"""

from chemart.helpers.explicit import network, term
from chemart.network import Species


def match(epitope: str, paratope: str, s: int) -> int:
    """m = sum over alignments k of G(complementary bits - s + 1), G(x) = max(x, 0)."""
    le, lp = len(epitope), len(paratope)
    total = 0
    for k in range(-(lp - 1), le):
        score = sum(epitope[n + k] != paratope[n] for n in range(lp) if 0 <= n + k < le)
        total += max(score - s + 1, 0)
    return total


def _bits(rng, n: int) -> str:
    return "".join(str(int(b)) for b in rng.integers(0, 2, n))


def generate(p, rng):
    antibodies = [(_bits(rng, p.l_e), _bits(rng, p.l_p)) for _ in range(p.N)]
    antigens = [_bits(rng, p.l_e) for _ in range(p.M)]
    X = [f"X{i + 1}" for i in range(p.N)]
    Y = [f"Y{j + 1}" for j in range(p.M)]
    # m[i][j]: epitope of X_i recognised by the paratope of X_j
    m = [[match(antibodies[i][0], antibodies[j][1], p.s) for j in range(p.N)] for i in range(p.N)]

    reactions = []
    for i in range(p.N):
        for j in range(p.N):
            if m[j][i]:   # X_i recognises X_j: stimulation of X_i
                text = f"{term(2, X[i])} -> {term(3, X[i])}" if i == j else f"{X[i]} + {X[j]} -> {term(2, X[i])} + {X[j]}"
                reactions.append((text, p.c * m[j][i]))
            if m[i][j]:   # X_j recognises X_i: suppression of X_i
                text = f"{term(2, X[i])} -> {X[i]}" if i == j else f"{X[i]} + {X[j]} -> {X[j]}"
                reactions.append((text, p.c * p.k1 * m[i][j]))
    for j in range(p.M):
        for i in range(p.N):
            strength = match(antigens[j], antibodies[i][1], p.s)
            if strength:
                reactions.append((f"{X[i]} + {Y[j]} -> {term(2, X[i])} + {Y[j]}", p.c * strength))
                reactions.append((f"{Y[j]} + {X[i]} -> {X[i]}", p.k3 * strength))
    reactions += [(f"{x} -> ", p.k2) for x in X]

    species = [Species(x, f"e={e} p={pa}") for x, (e, pa) in zip(X, antibodies)]
    species += [Species(y, f"e={e}") for y, e in zip(Y, antigens)]
    return network(reactions, species=species)
