"""Matrix chemistry (Banzhaf 1993; book chapter 3, 12.5.2, 13.2). Catalog id: matrix-chemistry.

Molecules are binary strings of length N (N a perfect square), named by the
integer they encode with s_1 as the least significant bit (s(5) = 1010). A
string s1 folds into a sqrt(N) x sqrt(N) operator P and acts on s2 in
sqrt(N)-sized chunks (eq. 3.4): output bit i + k sqrt(N) is 1 iff
sum_j P_ij s2_(j + k sqrt(N)) > Theta. The reaction s1 + s2 -> s1 + s2 + s3
keeps both reactants; a product equal to the all-zero destructor s(0) is an
elastic collision.

method="closure" (default) returns the closure from the seed strings (book
3.3), with mass-action constants that make the ODE with a constant-total
outflow the book's eq. 3.24. method="soup" runs the book's algorithm (3.1)
and returns the reactions that fired.
"""

from collections import Counter
from math import isqrt

import numpy as np

from chemart.expand import expand
from chemart.helpers.params import apportion
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import soup

FOLDINGS = (1, 2, 3, 4)
TABLE_MAX_N = 9


def side(N: int) -> int:
    n = isqrt(N)
    if N < 1 or n * n != N:
        raise ValueError(f"N must be a perfect square (1, 4, 9, 16, 25, ...), got N={N}")
    return n


def layout(n: int, folding: int) -> list[list[int]]:
    """layout[i][j] = index (0-based) of the string component placed at P_ij.

    1 row-major (eq. 3.16), 2 transposed (eq. 3.17), 3 row snake (odd rows
    reversed), 4 column snake (odd columns reversed).
    """
    if folding == 1:
        return [[i * n + j for j in range(n)] for i in range(n)]
    if folding == 2:
        return [[j * n + i for j in range(n)] for i in range(n)]
    if folding == 3:
        return [[i * n + (n - 1 - j if i % 2 else j) for j in range(n)] for i in range(n)]
    if folding == 4:
        return [[j * n + (n - 1 - i if j % 2 else i) for j in range(n)] for i in range(n)]
    raise ValueError(f"folding must be one of {FOLDINGS}, got {folding!r}")


def operation(N: int, folding: int = 1, Theta: float = 0.0):
    """Return product(operator, string) -> the integer name of the new string."""
    n = side(N)
    lay = layout(n, folding)
    chunk = (1 << n) - 1
    rows_cache: dict[int, list[int]] = {}

    def rows(op: int) -> list[int]:
        masks = rows_cache.get(op)
        if masks is None:
            masks = [sum(((op >> lay[i][j]) & 1) << j for j in range(n)) for i in range(n)]
            rows_cache[op] = masks
        return masks

    def product(op: int, s: int) -> int:
        out = 0
        masks = rows(op)
        for k in range(n):
            c = (s >> (k * n)) & chunk
            for i, m in enumerate(masks):
                if (m & c).bit_count() > Theta:
                    out |= 1 << (i + k * n)
        return out

    return product


def reaction_table(N: int, folding: int = 1, Theta: float = 0.0) -> np.ndarray:
    """table[operator, string] = product, for all 2^N strings including s(0) (N <= 9)."""
    n = side(N)
    if N > TABLE_MAX_N:
        raise ValueError(f"reaction_table needs N <= {TABLE_MAX_N} (2^N x 2^N entries), got N={N}")
    strings = np.arange(1 << N)
    bits = (strings[:, None] >> np.arange(N)) & 1
    P = bits[:, np.array(layout(n, folding))]
    chunks = bits.reshape(-1, n, n)
    out = np.einsum("aij,bkj->abki", P, chunks) > Theta
    return (out.reshape(1 << N, 1 << N, N).astype(np.int64) << np.arange(N)).sum(-1)


def bitstring(s: int, N: int) -> str:
    """Components s_1 ... s_N in the book's order: bitstring(5, 4) == '1010'."""
    return "".join(str((s >> i) & 1) for i in range(N))


def _seeds(p) -> list[int]:
    seeds = p.seed_species
    if not seeds or not all(isinstance(s, int) and not isinstance(s, bool) for s in seeds):
        raise ValueError(f"seed_species must be a non-empty list of string numbers, got {seeds!r}")
    if len(set(seeds)) != len(seeds):
        raise ValueError(f"seed_species has duplicates: {seeds!r}")
    top = (1 << p.N) - 1
    bad = [s for s in seeds if not 0 <= s <= top]
    if bad:
        raise ValueError(f"seed_species entries must be in 0..{top} for N={p.N}, got {bad}")
    if p.destructor_elastic and 0 in seeds:
        raise ValueError("seed_species contains 0, the destructor, which cannot react while "
                         "destructor_elastic is true; remove it or set destructor_elastic=false")
    return list(seeds)


def generate(p, rng):
    side(p.N)
    seeds = _seeds(p)
    product = operation(p.N, p.folding, p.Theta)

    def react(a, b):
        c = product(a, b)
        if p.destructor_elastic and c == 0:
            return None
        return (a, b, c)

    def multiplicity(a, b, c):
        """Ordered (operator, operand) pairs that give the multiset reaction a + b -> a + b + c."""
        if a == b:
            return 1
        return sum(react(x, y) == (x, y, c) for x, y in ((a, b), (b, a)))

    def species_id(s):
        return f"s{s}"

    def reaction(a, b, c, count=None):
        k = multiplicity(a, b, c)
        return Reaction.of([species_id(a), species_id(b)], [species_id(a), species_id(b), species_id(c)],
                           rate={"law": "mass-action", "k": float(k)}, count=count)

    extras = {"encoding": "species sK is the string whose integer name is K; structure lists "
                          "s_1..s_N, with s_1 the least significant bit"}

    if p.method == "closure":
        found, _, status = expand(react, seeds, arity=2, max_species=p.max_species, ordered=True)
        strings = sorted(found)
        known = set(strings)
        reactions = {}
        for a in strings:
            for b in strings:
                out = react(a, b)
                if out is None or out[2] not in known:
                    continue
                key = (min(a, b), max(a, b), out[2])
                if key not in reactions:
                    reactions[key] = reaction(*key)
        initial = {species_id(s): 1.0 / len(seeds) for s in seeds if s in known}
        extras["analysis"] = {
            "self_replicators": [species_id(s) for s in strings if react(s, s) == (s, s, s)],
        }
        reaction_list = list(reactions.values())
    else:
        counts = apportion(p.M, [1.0] * len(seeds))
        population = [s for s, m in zip(seeds, counts) for _ in range(m)]
        fired, final = soup(react, population, p.steps, rng, arity=2, dilution="constant")
        status = "observed"
        reaction_list = []
        for lhs, rhs, count in fired:
            a, b = sorted(lhs)
            (c,) = (Counter(rhs) - Counter(lhs)).elements()
            reaction_list.append(reaction(a, b, c, count))
        strings = sorted({*seeds, *(s for lhs, rhs, _ in fired for s in rhs)})
        initial = {species_id(s): int(m) for s, m in zip(seeds, counts) if m > 0}
        extras["analysis"] = {
            "steps": p.steps,
            "final_population": {species_id(s): n for s, n in sorted(Counter(final).items())},
        }

    return Network(
        species=[Species(species_id(s), structure=bitstring(s, p.N)) for s in strings],
        reactions=reaction_list,
        status=status,
        initial_state=initial or None,
        outflow=CONSTANT_TOTAL,
        extras=extras,
    )
