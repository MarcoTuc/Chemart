"""Analog computation of algebraic functions at steady state. Catalog id: analog-function-crn.

The answer is the steady-state concentration of Y. `sqrt` is the book's
example (eqs. 17.18-17.24, k1 = 2 k2); the others follow the same pattern
and are checked against their steady state in the tests.
"""

import math

from chemart.helpers.explicit import network

FUNCTIONS = {
    #          reactions                                                        y*
    "sqrt": ([("X1 -> X1 + Y", 2.0), ("2 Y -> ", 1.0)], lambda x1, x2: math.sqrt(x1)),
    "square": ([("2 X1 -> 2 X1 + Y", 1.0), ("Y -> ", 1.0)], lambda x1, x2: x1 * x1),
    "add": ([("X1 -> X1 + Y", 1.0), ("X2 -> X2 + Y", 1.0), ("Y -> ", 1.0)], lambda x1, x2: x1 + x2),
    "multiply": ([("X1 + X2 -> X1 + X2 + Y", 1.0), ("Y -> ", 1.0)], lambda x1, x2: x1 * x2),
    "divide": ([("X1 -> X1 + Y", 1.0), ("X2 + Y -> X2", 1.0)], lambda x1, x2: x1 / x2),
}


def generate(p, rng):
    if p.function == "divide" and p.x2 == 0:
        raise ValueError("divide needs x2 > 0")
    reactions, expected = FUNCTIONS[p.function]
    net = network(reactions, species=["X1", "X2", "Y"] if p.function in ("add", "multiply", "divide") else ["X1", "Y"])
    net.initial_state = {s.id: (p.x1 if s.id == "X1" else p.x2 if s.id == "X2" else 0.0) for s in net.species}
    net.extras = {"readout": {"species": "Y", "steady_state": expected(p.x1, p.x2)}}
    return net
