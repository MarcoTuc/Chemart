"""Colored chameleon chemistry. Catalog id: chameleon."""

from chemart.helpers.explicit import network
from chemart.helpers.params import apportion, vector

COLORS = ("r", "g", "b")


def generate(p, rng):
    x0 = vector("x0", p.x0, 3, minimum=0.0)
    if abs(sum(x0) - 1.0) > 1e-9:
        raise ValueError(f"x0 must be fractions of (r, g, b) summing to 1, got {p.x0!r}")
    return network(
        [("r + g -> 2 b", 1.0), ("r + b -> 2 g", 1.0), ("g + b -> 2 r", 1.0)],
        species=COLORS,
        initial_state=dict(zip(COLORS, apportion(p.M, x0))),
        extras={
            "conservation": [
                {"vector": {"r": 1, "g": 1, "b": 1}},
                {"vector": {"r": 1, "g": -1}, "modulus": 3},
                {"vector": {"g": 1, "b": -1}, "modulus": 3},
            ]
        },
    )
