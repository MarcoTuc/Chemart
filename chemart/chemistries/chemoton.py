"""Chemoton: model I of Fernando & Di Paolo (2004), after Csendes (1984). Catalog id: chemoton.

Metabolic A-wheel (eqs. 1-5), template polycondensation by replication stage
(eqs. 6-10) and membrane growth (eqs. 11-14), with the paper's appendix rate
constants and initial conditions.
"""

from chemart.helpers.explicit import network

RATES = {
    "k1": 2.0, "k1r": 0.1, "k2": 100.0, "k2r": 0.1, "k3": 100.0, "k3r": 0.1,
    "k4": 100.0, "k4r": 0.1, "k5": 10.0, "k5r": 0.1, "k6": 10.0, "k6r": 1.0,
    "k7": 10.0, "k8": 10.0, "k9": 10.0, "k9r": 0.1, "k10": 10.0,
}


def generate(p, rng):
    for name, value in p.rates.items():
        if name not in RATES:
            raise ValueError(f"rates: unknown constant {name!r}; valid names: {sorted(RATES)}")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f"rates[{name!r}] must be a non-negative number, got {value!r}")
    k = {**RATES, **p.rates}
    N = p.N
    initiation = {"law": "mass-action", "k": k["k6"], "threshold_species": "V", "threshold": p.V_threshold}

    reactions = [
        # metabolic subsystem (eqs. 1-5)
        ("A1 + X -> A2", k["k1"]), ("A2 -> A1 + X", k["k1r"]),
        ("A2 -> A3 + Y", k["k2"]), ("A3 + Y -> A2", k["k2r"]),
        ("A3 -> A4 + V", k["k3"]), ("A4 + V -> A3", k["k3r"]),
        ("A4 -> A5 + T0", k["k4"]), ("A5 + T0 -> A4", k["k4r"]),
        ("A5 -> 2 A1", k["k5"]), ("2 A1 -> A5", k["k5r"]),
        # template subsystem (eqs. 6-10): pV{r} carries r extra bound monomers
        ("pV0 + V -> pV1 + R", initiation), ("pV1 + R -> pV0 + V", k["k6r"]),
        *((f"pV{r} + V -> pV{r + 1} + R", k["k7"]) for r in range(1, N - 1)),
        (f"pV{N - 1} + V -> 2 pV0 + R", k["k7"]),
        # membrane subsystem (eqs. 11-14)
        ("T0 -> T*", k["k8"]), ("T* + R -> T", k["k9"]), ("T -> T* + R", k["k9r"]),
        ("T + S -> 2 S", k["k10"]),
    ]
    species = ["X", "Y", "A1", "A2", "A3", "A4", "A5", "V", "R", "T0", "T*", "T", "S"]
    species += [f"pV{r}" for r in range(N)]
    initial = {
        "X": p.X, "Y": 0.1, "A1": 1.0, "A2": 1.8, "A3": 1.9, "A4": 1.7, "A5": 10.0,
        "V": 26.0, "R": 0.0, "T0": 17.0, "T*": 14.0, "T": 0.0, "S": 1.0, "pV0": 0.01,
    }
    return network(
        reactions,
        species=species,
        initial_state=initial,
        extras={
            "buffered": ["X", "Y"],
            "compartments": {
                "cell": {
                    "volume": "Q = S^1.5, initially 1; concentrations are rescaled by Q(t)/Q(t+dt)",
                    "division": "when S has doubled, the cell splits into two with half of every amount",
                }
            },
        },
    )
