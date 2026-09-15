"""Logistic growth as 'replicate and fight'. Catalog id: logistic-chemistry."""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network(
        [("X -> 2 X", p.r), ("2 X -> X", p.r / p.K)],
        initial_state={"X": p.x0},
    )
