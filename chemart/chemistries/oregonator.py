"""Oregonator (Belousov-Zhabotinsky). Catalog id: oregonator."""

from chemart.helpers.explicit import network, term

SPECIES = ["A", "B", "X", "Y", "Z", "P", "Q"]


def _regenerate_inhibitor(f: float, k5: float) -> list[tuple[str, float]]:
    """Z --k5--> f Y for real f, as integer-stoichiometry reactions with the same ODE.

    With lo = floor(f): Z -> lo Y at k5 (1 - frac) and Z -> (lo+1) Y at k5 frac.
    """
    lo = int(f)
    frac = f - lo
    out = []
    for n, weight in ((lo, 1.0 - frac), (lo + 1, frac)):
        if weight > 0:
            out.append((f"Z -> {term(n, 'Y') if n else ''}", k5 * weight))
    return out


def generate(p, rng):
    for name, value in p.D.items():
        if name not in SPECIES or isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f"D must map species in {SPECIES} to non-negative numbers, got {name!r}: {value!r}")
    return network(
        [
            ("A + Y -> X", p.k1),
            ("X + Y -> P", p.k2),
            ("B + X -> 2 X + Z", p.k3),
            ("2 X -> Q", p.k4),
            *_regenerate_inhibitor(p.f, p.k5),
        ],
        species=SPECIES,
        initial_state={"A": p.A, "B": p.B},
        extras={
            "buffered": ["A", "B"],
            "space": {"lattice": p.lattice, "diffusion": dict(p.D)},
        },
    )
