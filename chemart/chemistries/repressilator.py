"""Repressilator. Catalog id: repressilator."""

from chemart.helpers.explicit import network, term

GENES = (1, 2, 3)


def generate(p, rng):
    reactions = []
    for i in GENES:
        repressor = term(p.n, f"P{3 if i == 1 else i - 1}")
        reactions += [
            (f"G{i} + {repressor} -> C{i}", p.ke),
            (f"C{i} -> G{i} + {repressor}", p.kr),
        ]
    for i in GENES:
        reactions += [
            (f"G{i} -> G{i} + M{i}", p.km),
            (f"M{i} -> M{i} + P{i}", p.kp),
            (f"M{i} -> ", p.mu_m),
            (f"P{i} -> ", p.mu_p),
        ]
    species = [f"{kind}{i}" for kind in "GCMP" for i in GENES]
    initial = {s: 0.0 for s in species}
    initial.update(G1=1.0, C2=1.0, C3=1.0)
    return network(reactions, species=species, initial_state=initial)
