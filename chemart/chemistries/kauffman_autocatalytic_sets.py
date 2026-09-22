"""Kauffman autocatalytic sets: the binary polymer model (book 6.3.1, eqs. 6.2-6.3).

Catalog id: kauffman-autocatalytic-sets. Every polymer up to max_length is a
species; every way to split a polymer gives a condensation/cleavage pair
A + B <-> AB; each polymer catalyses each pair with probability P.
Topology only.
"""

from itertools import product

import numpy as np

from chemart.helpers.explicit import network, term

ALPHABET = "abcd"
MAX_SPECIES = 512


def _join(a: str, b: str) -> str:
    return term(2, a) if a == b else f"{a} + {b}"


def generate(p, rng):
    letters = ALPHABET[: p.B]
    species = ["".join(s) for n in range(1, p.max_length + 1) for s in product(letters, repeat=n)]
    if len(species) > MAX_SPECIES:
        raise ValueError(f"B = {p.B}, max_length = {p.max_length} gives {len(species)} polymers, above the limit of {MAX_SPECIES}")
    unknown = [f for f in p.food_set if f not in set(species)]
    if unknown:
        raise ValueError(f"food_set entries {unknown} are not polymers over {letters!r} of length <= {p.max_length}")

    pairs = [(s[:cut], s[cut:], s) for s in species for cut in range(1, len(s))]
    reactions = []
    for a, b, c in pairs:
        reactions += [(f"{_join(a, b)} -> {c}", None), (f"{c} -> {_join(a, b)}", None)]
    catalysed = rng.random((len(pairs), len(species))) < p.P
    for pair, catalyst in zip(*np.nonzero(catalysed)):
        a, b, c = pairs[pair]
        e = species[catalyst]
        reactions += [(f"{_join(a, b)} + {e} -> {c} + {e}", None), (f"{c} + {e} -> {_join(a, b)} + {e}", None)]

    return network(
        reactions,
        species=species,
        extras={
            "food": list(p.food_set),
            "conservation": [
                {"name": f"monomer {x}", "vector": {s: s.count(x) for s in species if x in s}}
                for x in letters
            ],
        },
    )
