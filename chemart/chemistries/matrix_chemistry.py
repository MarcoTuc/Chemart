"""Matrix chemistry (Banzhaf 1993; book chapter 3, 12.5.2, 13.2). Catalog id: matrix-chemistry.

Molecules are binary strings of length N (N a perfect square), named by the
integer they encode with s_1 as the least significant bit (s(5) = 1010). A
string s1 folds into a sqrt(N) x sqrt(N) operator P and acts on s2 in
sqrt(N)-sized chunks (eq. 3.4): output bit i + k sqrt(N) is 1 iff
sum_j P_ij s2_(j + k sqrt(N)) > Theta. The reaction s1 + s2 -> s1 + s2 + s3
keeps both reactants; a product equal to the all-zero destructor s(0) is an
elastic collision.

generate returns the closure from the seed strings (book 3.3), with
mass-action constants that make the ODE with a constant-total outflow the
book's eq. 3.24. evolve runs the book's algorithm (3.1), a frame every M
collisions, and returns the reactions that fired.
"""

from collections import Counter
from math import isqrt

import numpy as np

from chemart.expand import expand
from chemart.helpers.params import apportion
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

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


class _Chemistry:
    """The reaction rule and reaction records shared by both faces."""

    def __init__(self, p):
        side(p.N)
        self.p = p
        self.product = operation(p.N, p.folding, p.Theta)

    def react(self, a, b):
        c = self.product(a, b)
        if self.p.destructor_elastic and c == 0:
            return None
        return (a, b, c)

    def multiplicity(self, a, b, c):
        """Ordered (operator, operand) pairs that give the multiset reaction a + b -> a + b + c."""
        if a == b:
            return 1
        return sum(self.react(x, y) == (x, y, c) for x, y in ((a, b), (b, a)))

    def reaction(self, a, b, c, count=None):
        k = self.multiplicity(a, b, c)
        return Reaction.of([species_id(a), species_id(b)], [species_id(a), species_id(b), species_id(c)],
                           rate={"law": "mass-action", "k": float(k)}, count=count)

    def network(self, strings, reactions, status, initial, analysis, **extras):
        return Network(
            species=[Species(species_id(s), structure=bitstring(s, self.p.N)) for s in strings],
            reactions=reactions,
            status=status,
            initial_state=initial or None,
            outflow=CONSTANT_TOTAL,
            extras={"encoding": ENCODING, "analysis": analysis, **extras},
        )


ENCODING = ("species sK is the string whose integer name is K; structure lists "
            "s_1..s_N, with s_1 the least significant bit")


def species_id(s: int) -> str:
    return f"s{s}"


def generate(p, rng):
    """The closure from seed_species (book 3.3), cut off by max_species."""
    chem = _Chemistry(p)
    seeds = _seeds(p)
    react = chem.react
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
                reactions[key] = chem.reaction(*key)
    initial = {species_id(s): 1.0 / len(seeds) for s in seeds if s in known}
    analysis = {"self_replicators": [species_id(s) for s in strings if react(s, s) == (s, s, s)]}
    return chem.network(strings, list(reactions.values()), status, initial, analysis)


def _ids(fired):
    return [[[species_id(s) for s in lhs], [species_id(s) for s in rhs], n] for lhs, rhs, n in fired]


def evolve(p, rng):
    """The book's algorithm (3.1): M strings, `steps` collisions, a frame every M collisions."""
    chem = _Chemistry(p)
    seeds = _seeds(p)
    counts = apportion(p.M, [1.0] * len(seeds))
    population = [s for s, m in zip(seeds, counts) for _ in range(m)]
    tally = Tally()
    final = population
    for step, final, tally in stir(chem.react, population, p.steps, rng, arity=2, dilution="constant",
                                   tally=tally):
        yield Frame(t=float(step), state={species_id(s): float(n) for s, n in sorted(Counter(final).items())},
                    fired=_ids(tally.flush()))
    fired = tally.reactions()
    reactions = []
    for lhs, rhs, count in fired:
        a, b = sorted(lhs)
        (c,) = (Counter(rhs) - Counter(lhs)).elements()
        reactions.append(chem.reaction(a, b, c, count))
    strings = sorted({*seeds, *(s for lhs, rhs, _ in fired for s in rhs)})
    initial = {species_id(s): int(m) for s, m in zip(seeds, counts) if m > 0}
    analysis = {"steps": p.steps,
                "final_population": {species_id(s): n for s, n in sorted(Counter(final).items())}}
    return chem.network(strings, reactions, "observed", initial, analysis)
