"""Selection equation under a dilution flow, x_i' = f_i x_i^c - Phi x_i. Catalog id: selection-equation."""

from chemart.helpers.explicit import network
from chemart.helpers.params import vector
from chemart.network import CONSTANT_TOTAL


def generate(p, rng):
    f = vector("f", p.f, minimum=0.0)
    n = len(f)
    if n == 0:
        raise ValueError("f must list at least one fitness")
    x0 = vector("x0", p.x0, n, minimum=0.0) if p.x0 else [1.0 / n] * n
    X = [f"X{i + 1}" for i in range(n)]
    reactions = []
    for s, fi in zip(X, f):
        if fi > 0:
            rate = fi if p.c == 1.0 else {"law": "power", "k": fi, "order": p.c}
            reactions.append((f"{s} -> 2 {s}", rate))
    return network(reactions, species=X, initial_state=dict(zip(X, x0)), outflow=CONSTANT_TOTAL)
