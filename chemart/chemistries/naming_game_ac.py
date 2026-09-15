"""Naming game as an artificial chemistry (book eqs. 16.2-16.6). Catalog id: naming-game-ac."""

from chemart.helpers.explicit import network


def generate(p, rng):
    names = range(1, p.N_s + 1)
    reactions = []
    for j in names:
        reactions += [
            (f"M + C{j} -> S{j}", p.ks),
            (f"M + S{j} -> C{j}", p.kc),
            (f"M + S{j} + C{j} -> 2 S{j}", p.ka),
        ]
    for j in names:
        for k in names:
            if k != j:
                reactions += [
                    (f"M + S{j} + C{k} -> M + S{j} + C{j}", p.kappa1),  # replace mismatching adaptor
                    (f"M + S{k} + C{j} -> M + S{j} + C{j}", p.kappa2),  # replace mismatching word
                ]
    species = ["M"] + [f"S{j}" for j in names] + [f"C{j}" for j in names]
    return network(reactions, species=species)
