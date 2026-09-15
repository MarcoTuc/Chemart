"""Reversible dimerization. Catalog id: dimerization."""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network(
        [("A + B -> C", p.k_f), ("C -> A + B", p.k_r)],
        initial_state={"A": p.A0, "B": p.B0, "C": 0.0},
    )
