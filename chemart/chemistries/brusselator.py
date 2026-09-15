"""Brusselator. Catalog id: brusselator."""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network(
        [
            ("A -> X", p.k1),
            ("B + X -> Y + D", p.k2),
            ("2 X + Y -> 3 X", p.k3),
            ("X -> E", p.k4),
        ],
        species=["A", "B", "X", "Y", "D", "E"],
        initial_state={"A": p.a, "B": p.b, "X": 0.0, "Y": 0.0, "D": 0.0, "E": 0.0},
        extras={"buffered": ["A", "B"]},
    )
