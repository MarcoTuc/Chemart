"""Flow artificial chemistry networks of Kreyssig & Dittrich (2011)."""

import numpy as np
import pytest
from odes import rhs

from chemart import generate_network


def test_r1_rate_equations_match_the_papers_pde_reaction_terms():
    net = generate_network("flow-ac")
    ids, f = rhs(net)
    c = dict(zip(ids, np.random.default_rng(1).uniform(0.1, 2.0, len(ids))))
    ours = dict(zip(ids, f(0.0, np.array([c[s] for s in ids]))))
    a, b, cc, d = c["a"], c["b"], c["c"], c["d"]
    paper = {
        "a": 0.0,
        "b": a * b - b * cc + cc - b * d - b,
        "c": b * cc - cc + b * d - cc,
        "d": a * d - b * d - d,
    }
    for species, value in paper.items():
        assert ours[species] == pytest.approx(value, abs=1e-12)


def test_membrane_network_has_no_spontaneous_membrane():
    net = generate_network("flow-ac", example="membrane")
    builds = [r for r in net.reactions if r.products == {"m": 1}]
    assert len(builds) == 7 and all(r.reactants for r in builds)
    assert len(net.reactions) == 13


def test_compartments_is_rock_scissors_paper():
    lines = generate_network("flow-ac", example="compartments").to_text().splitlines()
    assert len(lines) == 12
    assert "s1 + s2 -> s1  [mass-action k=1.0]" in lines and "2 s3 -> s3  [mass-action k=1.0]" in lines


def test_space_records_the_field():
    space = generate_network("flow-ac", example="membrane").extras["space"]
    assert space["reaction_radius"] == 0.15 and "exp(5 r)" in space["field"]
