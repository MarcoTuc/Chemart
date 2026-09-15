"""Replication and death. Catalog id: replication-death."""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network([("X -> 2 X", p.b), ("X -> ", p.d)], initial_state={"X": p.x0})
