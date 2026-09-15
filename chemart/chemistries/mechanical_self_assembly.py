"""Mechanical self-assembly of magnetic triangles (book eq. 20.1). Catalog id: mechanical-self-assembly."""

from chemart.helpers.explicit import network

SIZES = {"x": 1, "x2": 2, "x3": 3, "x4": 4, "x5": 5, "x6": 6}
PAIRS = [
    ("x", "x", "x2"), ("x", "x2", "x3"), ("x", "x3", "x4"), ("x", "x4", "x5"),
    ("x", "x5", "x6"), ("x2", "x2", "x4"), ("x2", "x3", "x5"), ("x2", "x4", "x6"),
    ("x3", "x3", "x6"),
]


def generate(p, rng):
    valid = {f"{a}+{b}" for a, b, _ in PAIRS}
    for key, value in p.P_b.items():
        if key not in valid:
            raise ValueError(f"P_b key {key!r} is not a bonding pair; valid keys: {sorted(valid)}")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError(f"P_b[{key!r}] must be a probability in [0, 1], got {value!r}")
    reactions = []
    for a, b, c in PAIRS:
        lhs = f"2 {a}" if a == b else f"{a} + {b}"
        reactions.append((f"{lhs} -> {c}", p.agitation_rate * p.P_b.get(f"{a}+{b}", 1.0)))
    return network(
        reactions,
        species=list(SIZES),
        initial_state={"x": p.n_monomers},
        extras={"conservation": [{"name": "monomer count", "vector": dict(SIZES)}]},
    )
