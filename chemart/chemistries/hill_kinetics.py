"""Hill kinetics (cooperative binding). Catalog id: hill-kinetics."""

from chemart.helpers.explicit import network, term


def generate(p, rng):
    if p.form == "elementary":
        bound = term(p.n, "P")
        expression = ("C -> C + X" if p.regulation == "activation" else "G -> G + X", p.k_expr)
        return network(
            [
                (f"G + {bound} -> C", p.kf),
                (f"C -> G + {bound}", p.kf * p.K ** p.n),  # K^n = K_d = kr / kf
                expression,
            ],
            species=["G", "P", "C", "X"],
            initial_state={"G": p.G0, "P": p.P0, "C": 0.0, "X": 0.0},
        )
    rate = {
        "law": "hill", "vmax": p.k_expr * p.G0, "K": p.K, "n": p.n,
        "regulator": "P", "mode": p.regulation,
    }
    return network([("-> X", rate)], species=["P", "X"], initial_state={"P": p.P0, "X": 0.0})
