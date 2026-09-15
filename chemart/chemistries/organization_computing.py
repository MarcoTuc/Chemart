"""Computing with chemical organizations: XOR gate and maximal independent set (book 17.3.3).

Catalog id: organization-computing.
"""

from chemart.helpers.explicit import network, term
from chemart.helpers.params import edges


def _xor(p):
    reactions = [(f"a{a} + b{b} -> c{a ^ b}", None) for a in (0, 1) for b in (0, 1)]
    reactions += [(f"{v}0 + {v}1 -> ", None) for v in "abc"]
    return network(
        reactions,
        species=[f"{v}{s}" for v in "abc" for s in (0, 1)],
        initial_state={f"a{int(p.a)}": 1, f"b{int(p.b)}": 1},
    )


def _maximal_independent_set(p):
    E = edges("graph", p.graph)
    vertices = sorted({v for e in E for v in e})
    neighbours = {v: sorted({b for a, b in E if a == v} | {a for a, b in E if b == v}) for v in vertices}
    reactions = [(f"s0_{v} + s1_{v} -> ", None) for v in vertices]                      # eq. 17.5
    for v in vertices:                                                                 # eq. 17.6
        lhs = " + ".join(f"s0_{u}" for u in neighbours[v])
        reactions.append((f"{lhs} -> {term(len(neighbours[v]), f's1_{v}')}", None))
    for v in vertices:                                                                 # eq. 17.7
        reactions += [(f"s1_{u} -> s0_{v}", None) for u in neighbours[v]]
    species = [f"s{state}_{v}" for v in vertices for state in (0, 1)]
    return network(reactions, species=species)


def generate(p, rng):
    return _xor(p) if p.problem == "xor" else _maximal_independent_set(p)
