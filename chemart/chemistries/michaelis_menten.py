"""Michaelis-Menten enzyme kinetics. Catalog id: michaelis-menten."""

from chemart.helpers.explicit import network


def generate(p, rng):
    if p.form == "elementary":
        return network(
            [("E + S -> ES", p.ka), ("ES -> E + S", p.ka_rev), ("ES -> E + P", p.kb)],
            species=["E", "S", "ES", "P"],
            initial_state={"E": p.E0, "S": p.S0, "ES": 0.0, "P": 0.0},
        )
    if p.ka == 0:
        raise ValueError("ka must be > 0 for the abridged form, since k_m = (ka_rev + kb) / ka")
    rate = {"law": "michaelis-menten", "vmax": p.kb * p.E0, "km": (p.ka_rev + p.kb) / p.ka}
    return network([("S -> P", rate)], species=["S", "P"], initial_state={"S": p.S0, "P": 0.0})
