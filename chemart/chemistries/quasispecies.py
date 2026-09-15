"""Quasispecies equation over all B^L genomes (book eqs. 7.17-7.19). Catalog id: quasispecies."""

from itertools import product

from chemart.helpers.explicit import network
from chemart.network import CONSTANT_TOTAL

MAX_GENOMES = 256


def generate(p, rng):
    genomes_total = p.B ** p.L
    if genomes_total > MAX_GENOMES:
        raise ValueError(
            f"B^L = {genomes_total} genomes exceeds the limit of {MAX_GENOMES} "
            "(the network has up to (B^L)^2 reactions); lower L or B"
        )
    if p.m > p.L:
        raise ValueError(f"m (mutations per genome) cannot exceed L = {p.L}, got {p.m}")
    genomes = ["".join(g) for g in product("0123"[: p.B], repeat=p.L)]
    per_site = p.m / p.L
    top = str(p.B - 1)

    def fitness(g: str) -> float:
        if p.fitness == "ones-fraction":
            return g.count(top) / p.L
        return 2.0 if g == top * p.L else 1.0

    reactions = []
    for gj in genomes:
        fj = fitness(gj)
        if fj == 0:
            continue
        for gi in genomes:
            d = sum(a != b for a, b in zip(gj, gi))
            q = (per_site / (p.B - 1)) ** d * (1 - per_site) ** (p.L - d)
            if q > 0:
                text = f"{gj} -> 2 {gj}" if gi == gj else f"{gj} -> {gj} + {gi}"
                reactions.append((text, fj * q))
    return network(reactions, species=genomes, outflow=CONSTANT_TOTAL)
