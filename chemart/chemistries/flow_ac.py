"""Reaction flow artificial chemistries (Kreyssig & Dittrich, ECAL 2011). Catalog id: flow-ac.

The reaction networks of the paper's examples, with the vector field and
stochastic-simulation settings of each example recorded in extras.space.
"""

from itertools import combinations, product

from chemart.helpers.explicit import network, term

SWIRL = "V1(x, y) = (1/sqrt(2)) ((1 - sqrt(2)) x - y, x + (1 - sqrt(2)) y)"

R1 = [
    "a + b -> a + 2 b", "a + d -> a + 2 d", "b + c -> 2 c", "c -> b",
    "b + d -> c", "b -> ", "c -> ", "d -> ",
]


def _membrane():
    producers = ("p1", "p2", "p3")
    build = []
    for counts in product((0, 1), repeat=3):
        if any(counts):  # the all-zero combination would make {} not closed
            build.append(" + ".join(p for p, c in zip(producers, counts) if c) + " -> m")
    regenerate = [
        "p1 + p2 + p3 -> 4 p1 + p2", "p1 + p2 + p3 -> 2 p1 + p3",
        "p1 + p2 + p3 -> 4 p2 + p1", "p1 + p2 + p3 -> 2 p2 + p3",
        "p1 + p2 + p3 -> 4 p3 + p1", "p1 + p2 + p3 -> 2 p3 + p2",
    ]
    return build + regenerate


def _rock_scissors_paper():
    s = ("s1", "s2", "s3")
    out = []
    for x in s:
        out += [f"{x} -> {term(2, x)}", f"{term(2, x)} -> {x}"]
    for x, y in combinations(s, 2):
        out += [f"{x} + {y} -> {x}", f"{x} + {y} -> {y}"]
    return out


EXAMPLES = {
    "swirl": (R1, {
        "field": SWIRL, "molecules": 2500, "initial_radius": 1.0,
        "reaction_radius": 0.1, "fraction_reacting": 1.0, "iterations": 400,
    }),
    "radial": (R1, {
        "field": "V2(x, y) = cos(10 r) / (10 r) (x, y), r = sqrt(x^2 + y^2)", "molecules": 10000,
        "initial_radius": 1.0, "reaction_radius": 0.01, "fraction_reacting": 1.0, "iterations": 3,
    }),
    "membrane": (_membrane(), {
        "field": "V(x, y) = 0.0005 exp(5 r) (x, y) for 0.2 < r < 0.8; V1(x, y) otherwise",
        "molecules": 2500, "initial_radius": 0.2, "reaction_radius": 0.15,
        "fraction_reacting": 1.0, "iterations": 100,
    }),
    "compartments": (_rock_scissors_paper(), {
        "field": "V3(x, y) = -0.01 exp(5 x) (x, 0) + V1(x + 0.4, y) for x < 0; 0.01 exp(-5 x) (x, 0) + V1(x - 0.4, y) for x > 0",
        "molecules": 2500, "initial_radius": 0.8, "reaction_radius": 0.1,
        "fraction_reacting": 1.0, "iterations": 150,
    }),
}


def generate(p, rng):
    reactions, simulation = EXAMPLES[p.example]
    return network(
        [(text, p.k) for text in reactions],
        extras={"space": {"dimensions": 2, "domain": [[-1.0, 1.0], [-1.0, 1.0]], "swirl": SWIRL, **simulation}},
    )
