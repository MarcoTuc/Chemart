"""Okamoto's biochemical switch (book eqs. 17.8-17.13). Catalog id: okamoto-switch."""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network(
        [
            ("I1 -> I1 + X1", p.k_in),
            ("I2 -> I2 + X3", p.k_in),
            ("A + X3 -> B + X4", p.k2),
            ("B + X1 -> A + X2", p.k1),
            ("X2 -> ", p.k3),
            ("X4 -> ", p.k4),
            ("I1 -> I2", p.k_conv),
        ],
        species=["I1", "I2", "X1", "X2", "X3", "X4", "A", "B"],
        initial_state={
            "I1": p.I1_0, "I2": p.I2_0, "X1": 0.0, "X2": 8.0, "X3": 0.0, "X4": 8.0,
            "A": 1.0, "B": 0.0,
        },
    )
