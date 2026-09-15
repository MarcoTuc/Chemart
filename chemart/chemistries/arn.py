"""Artificial Regulatory Network of Banzhaf (2003) (book 11.3.4). Catalog id: arn.

A bitstring genome is scanned for promoters. Each gene has, immediately
upstream of its promoter, a 32-bit inhibitor site followed by a 32-bit
enhancer site, and immediately downstream five 32-bit segments that encode a
32-bit protein by bitwise majority voting (Kuo, Leier & Banzhaf 2004, Fig. 1).
Protein j regulates gene i with strength exp(beta (d_ij - d_max)), where
d_ij is the number of complementary bits (XOR popcount) between protein j and
a regulatory site of gene i.

On the simplex sum_j c_j = 1 that the model imposes, eqs. 11.12-11.14 are
exactly mass-action kinetics under the non-selective constant-total dilution:

    P_j + P_i -> P_j + 2 P_i   k = delta/N exp(beta (d^a_ij - d_max))   (enhancement)
    P_j + P_i -> P_j           k = delta/N exp(beta (d^h_ij - d_max))   (inhibition)

(with 2 P_i on the left when j = i), since their production term is
sum_j k c_j c_i = delta (a_i - h_i) c_i.
"""

import numpy as np

from chemart.helpers.explicit import network, term
from chemart.network import CONSTANT_TOTAL, Species

SEGMENTS = 5  # l_gene: 32-bit integers per gene (all sources)


def bits(text: str) -> np.ndarray:
    return np.frombuffer(text.encode(), dtype=np.uint8) - ord("0")


def text(b) -> str:
    return "".join("1" if v else "0" for v in b)


def duplication_divergence(rng, duplications: int, mutation_rate: float, word: int = 32) -> np.ndarray:
    """Random `word`-bit string, then whole-length duplications, each followed by
    independent bit flips with probability `mutation_rate` over the whole genome."""
    g = rng.integers(0, 2, word, dtype=np.uint8)
    for _ in range(duplications):
        g = np.concatenate([g, g])
        g ^= (rng.random(g.size) < mutation_rate).astype(np.uint8)
    return g


def promoters(genome: np.ndarray, promoter: np.ndarray) -> np.ndarray:
    """Start positions of promoters. A match that overlaps an earlier match (a
    periodic extension such as 0101010101) is not a new promoter: the whole run
    counts once, at its first bit."""
    n = promoter.size
    if genome.size < n:
        return np.zeros(0, dtype=np.int64)
    windows = np.lib.stride_tricks.sliding_window_view(genome, n)
    raw = np.flatnonzero((windows == promoter).all(axis=1))
    if raw.size == 0:
        return raw
    keep = np.concatenate([[True], np.diff(raw) >= n])
    return raw[keep]


def genes(genome: np.ndarray, promoter: np.ndarray, word: int = 32):
    """(position, inhibitor site, enhancer site, coding region) for every promoter
    with both regulatory sites and the whole coding region inside the genome."""
    out = []
    n = promoter.size
    for p in promoters(genome, promoter):
        start, end = p + n, p + n + SEGMENTS * word
        if p >= 2 * word and end <= genome.size:
            out.append((int(p), genome[p - 2 * word:p - word], genome[p - word:p], genome[start:end]))
    return out


def majority(coding: np.ndarray, word: int, rng) -> np.ndarray:
    """Bit k of the protein is the majority of bit k over the gene's segments;
    ties (only possible with an even number of segments) are resolved by chance."""
    votes = coding.reshape(-1, word).sum(axis=0).astype(int) * 2
    segments = coding.size // word
    protein = (votes > segments).astype(np.uint8)
    ties = votes == segments
    if ties.any():
        protein[ties] = rng.integers(0, 2, int(ties.sum()), dtype=np.uint8)
    return protein


def generate(p, rng):
    word = p.word_bits
    if set(p.promoter) - {"0", "1"} or not p.promoter:
        raise ValueError(f"promoter must be a non-empty string of 0s and 1s, got {p.promoter!r}")
    promoter = bits(p.promoter)
    if p.genome_source == "random":
        genome = rng.integers(0, 2, p.genome_length, dtype=np.uint8)
    else:
        doublings = p.genome_length // word
        if p.genome_length % word or doublings & (doublings - 1):
            raise ValueError(
                f"genome_length {p.genome_length} is not word_bits * 2**k: a duplication-divergence genome "
                f"doubles a {word}-bit string k times (e.g. {word * 2 ** 12} for the papers' 12 duplications)")
        genome = duplication_divergence(rng, doublings.bit_length() - 1, p.mutation_rate, word)

    found = genes(genome, promoter, word)
    N = len(found)
    P = [f"P{i + 1}" for i in range(N)]
    species, reactions = [], []
    if N:
        inhibitor = np.array([g[1] for g in found])
        enhancer = np.array([g[2] for g in found])
        protein = np.array([majority(g[3], word, rng) for g in found])
        # d[i, j]: complementary bits between protein j and a site of gene i
        d_a = (enhancer[:, None, :] != protein[None, :, :]).sum(axis=2)
        d_h = (inhibitor[:, None, :] != protein[None, :, :]).sum(axis=2)
        species = [
            Species(P[i], f"position={pos} inhibitor={text(inhibitor[i])} "
                          f"enhancer={text(enhancer[i])} protein={text(protein[i])}")
            for i, (pos, *_) in enumerate(found)
        ]
        for i in range(N):
            for j in range(N):
                for site, d in (("enhancer", int(d_a[i, j])), ("inhibitor", int(d_h[i, j]))):
                    if d < p.link_threshold:
                        continue
                    k = p.delta / N * float(np.exp(p.beta * (d - word)))
                    if i == j:
                        rhs = term(3, P[i]) if site == "enhancer" else P[i]
                        txt = f"{term(2, P[i])} -> {rhs}"
                    else:
                        rhs = f"{P[j]} + {term(2, P[i])}" if site == "enhancer" else P[j]
                        txt = f"{P[j]} + {P[i]} -> {rhs}"
                    reactions.append((txt, {"law": "mass-action", "k": k, "site": site, "match": d}))

    return network(
        reactions,
        species=species,
        initial_state={s: 1.0 / N for s in P} if N else None,
        outflow=CONSTANT_TOTAL if N else None,
        extras={
            "genome": text(genome),
            "analysis": {"genes": N, "genome_length": int(genome.size), "d_max": word},
        },
    )
