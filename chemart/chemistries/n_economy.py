"""N-economy production system: natural numbers as goods (book eqs. 20.6-20.7). Catalog id: n-economy."""

from chemart.helpers.explicit import network, term
from chemart.network import Species


def _factorisation(n: int) -> str:
    factors, d = [], 2
    while d * d <= n:
        k = 0
        while n % d == 0:
            n //= d
            k += 1
        if k:
            factors.append(f"{d}^{k}" if k > 1 else str(d))
        d += 1
    if n > 1:
        factors.append(str(n))
    return "*".join(factors)


def _side(goods: set[str], side: dict, where: str) -> str:
    for good, amount in side.items():
        if good not in goods:
            raise ValueError(f"{where}: good {good!r} is not in product_set")
        if isinstance(amount, bool) or not isinstance(amount, int) or amount < 1:
            raise ValueError(f"{where}: amount of {good!r} must be a positive integer, got {amount!r}")
    return " + ".join(term(amount, good) for good, amount in side.items())


def generate(p, rng):
    if not all(isinstance(g, int) and not isinstance(g, bool) and g >= 2 for g in p.product_set):
        raise ValueError(f"product_set must list integers >= 2, got {p.product_set!r}")
    goods = [str(g) for g in dict.fromkeys(p.product_set)]
    reactions = []
    for i, step in enumerate(p.technology):
        if not isinstance(step, dict) or set(step) != {"in", "out"}:
            raise ValueError(f"technology[{i}] must be an object with keys 'in' and 'out', got {step!r}")
        lhs = _side(set(goods), step["in"], f"technology[{i}].in")
        rhs = _side(set(goods), step["out"], f"technology[{i}].out")
        reactions.append((f"{lhs} -> {rhs}", None))
    species = [Species(g, _factorisation(int(g))) for g in goods]
    return network(reactions, species=species)
