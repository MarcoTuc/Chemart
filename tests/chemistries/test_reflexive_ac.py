"""Reflexive AC: machines composing machines.

Published results from Salzberg & Sayama (2025), 'Reflexive composition of elementary
state machines', arXiv:2505.07186 (Table 1, Table 2, eqs. 1-7, figs. 2-17), which restates
the formulation of Salzberg (2007) [737]; tapes as graphs from Salzberg, Sayama & Ikegami
(2004), ALife IX, fig. 2.
"""

from collections import Counter

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries.reflexive_ac import canonical, compose, elementary, number, parse, to_id

ID = "reflexive-ac"


# --- helpers: the paper's transformations and lattice ----------------------------
def complement(m):
    """Swap the message values 0 and 1."""
    return tuple(tuple((d, 1 - o) for d, o in (node[1], node[0])) for node in m)


def mirror(m):
    """Swap the states S0 and S1."""
    return tuple(tuple((1 - d, o) for d, o in node) for node in (m[1], m[0]))


def transpose(m):
    """Fig. 17: states become messages and messages states."""
    return tuple(tuple((m[x][s][1], m[x][s][0]) for x in range(2)) for s in range(2))


def lattice_step(m, q, i):
    """Eqs. 3-4: pipe the input through the lattice; states land in reversed positions."""
    new, msg = [None] * len(q), i
    for k, state in enumerate(q):
        d, msg = m[state][msg]
        new[len(q) - 1 - k] = d
    return new, msg


def rule90(c):
    return [(c[j - 1] if j else 0) ^ (c[j + 1] if j < len(c) - 1 else 0) for j in range(len(c))]


# --- naming scheme -----------------------------------------------------------------
def test_table1_naming_scheme():
    # fig. 2: M61 = 00111101; S0: 0 -> S0/0, 1 -> S1/1; S1: 0 -> S1/1, 1 -> S0/1
    assert elementary(61) == (((0, 0), (1, 1)), ((1, 1), (0, 1)))
    # fig. 4: M45 = 00101101, the adding machine modulo 2 that reports its state
    m45 = elementary(45)
    assert all(m45[q][i] == (q ^ i, q) for q in (0, 1) for i in (0, 1))
    assert all(number(elementary(n)) == n for n in range(256))
    assert parse("M45") == parse("0/0,1/0;1/1,0/1") and to_id(parse("M45")) == "0/0,1/0;1/1,0/1"


def test_table2_76_unique_machines_and_equivalents():
    orbits = {frozenset(number(t) for t in (m, complement(m), mirror(m), mirror(complement(m))))
              for m in map(elementary, range(256))}
    assert len(orbits) == 76
    table2 = {7: {47, 88, 218}, 44: {104, 199, 214}, 45: {120, 135, 210},
              54: {99, 156, 201}, 60: {105, 150, 195}, 61: {121, 131, 146}}
    for n, eq in table2.items():
        assert next(o for o in orbits if n in o) == {n} | eq


def test_fig17_transposition():
    assert number(transpose(elementary(45))) == 54
    assert number(transpose(elementary(54))) == 45
    assert number(transpose(elementary(60))) == 60          # footnote 3: M60 transposes into itself


# --- the reaction: reflexive composition ----------------------------------------------
def test_product_follows_eqs_1_2():
    """The product graph, walked from its pointer, outputs phi_r(phi_s(i)) with sender and receiver swapping."""
    rng = np.random.default_rng(0)
    for _ in range(300):
        a, b = (canonical(elementary(int(rng.integers(256))), int(rng.integers(2))) for _ in range(2))
        p = compose(a, b)
        same = a == b
        assert len(p) <= (len(a) ** 2 if same else 2 * len(a) * len(b))
        machines, x, y, node = [a, b], 0, 0, 0
        for i in rng.integers(2, size=30):
            s, r = machines
            dx, msg = s[x][int(i)]
            dy, o = r[y][msg]
            node, out = p[node][int(i)]
            assert out == o
            machines, x, y = [r, s], dy, dx


def test_nested_products_are_the_lattice_of_eqs_3_4():
    cells = [0, 1, 0, 0, 1]                                  # fig. 4(b) initial states
    product = canonical(elementary(45), cells[0])
    for c in cells[1:]:
        product = compose(product, canonical(elementary(45), c))
    assert len(product) == 2 * 2 ** len(cells)               # every configuration, in both orientations
    rng = np.random.default_rng(1)
    for _ in range(20):
        q, node = list(cells), 0
        for i in map(int, rng.integers(2, size=40)):
            q, o = lattice_step(elementary(45), q, i)
            node, out = product[node][i]
            assert out == o


def test_m45_is_rule_90_and_draws_the_sierpinski_triangle():
    # figs. 5-7: 19 cells, one centred S1; inputs i(t_even) = 0, i(t_odd) = o(t - 1)
    q = [0] * 19
    q[9] = 1
    for t in range(9):
        assert sum(q) == 2 ** bin(t).count("1")              # Sierpinski rows
        half, o = lattice_step(elementary(45), q, 0)
        nxt, _ = lattice_step(elementary(45), half, o)
        assert nxt == rule90(q)
        q = nxt


def test_m54_runs_rule_90_backwards():
    # figs. 12 and 14: 8 cells, two centred 1 cells; eqs. 6-7: i(t_odd) = o(t - 1), o(t_odd) = 0, try i = 0 first
    m54 = elementary(54)
    q = [0, 0, 0, 1, 1, 0, 0, 0]
    for _ in range(8):
        for i in (0, 1):
            half, o = lattice_step(m54, q, i)
            nxt, o_odd = lattice_step(m54, half, o)
            if o_odd == 0:
                break
        assert o_odd == 0
        assert rule90(nxt) == q                              # the new row is a rule-90 preimage
        q = nxt


def test_m44_is_inverted_by_flipping_the_lattice():
    # fig. 8: 'M44 is reversible and can be inverted simply by flipping the lattice of states'
    rng = np.random.default_rng(2)
    for _ in range(50):
        c = [int(x) for x in rng.integers(2, size=20)]
        after, _ = lattice_step(elementary(44), c, 0)
        back, _ = lattice_step(elementary(44), after[::-1], 0)
        assert back[::-1] == c


def test_tapes_are_graphs_too():
    # ALife IX fig. 2: a tape is a chain of xi/a links (xi = symbol 2 here); reading it with M45 composes
    # the first tape link with the machine; the machine has no xi link, so the product stops there.
    tape = parse("-,-,1/1;-,-,2/0;-,-,-")
    m45 = parse("0/0,1/0,-;1/1,0/1,-")
    assert compose(tape, m45) == ((None, None, (1, 0)), (None, None, None))
    assert compose(m45, tape) == ((None, None, None),)


# --- species identity --------------------------------------------------------------
def test_ids_are_invariant_under_relabelling():
    rng = np.random.default_rng(3)
    for _ in range(50):
        n = 6
        nodes = [tuple((int(rng.integers(n)), int(rng.integers(2))) for _ in range(2)) for _ in range(n)]
        perm = [int(x) for x in rng.permutation(n)]
        relabelled = [None] * n
        for j, node in enumerate(nodes):
            relabelled[perm[j]] = tuple((perm[d], o) for d, o in node)
        a, b = canonical(tuple(nodes), 0), canonical(tuple(relabelled), perm[0])
        assert a == b and parse(to_id(a)) == a
    assert to_id(parse("M45@S1")) == "0/1,1/1;1/0,0/0"


# --- networks ------------------------------------------------------------------------
def states(net, sid):
    return net.extras["analysis"]["states"][sid]


def test_default_closure_grows_exponentially():
    net = generate_network(ID)
    assert net.status == "complete"
    m45 = to_id(parse("M45"))
    assert net.species[0].id == m45 and net.species[0].structure.startswith("M45@S0")
    assert sorted(states(net, s.id) for s in net.species) == [2, 4, 16, 16, 32, 64]
    for r in net.reactions:
        assert r.catalysts == r.reactants and r.rate is None
        (new,) = (Counter(r.products) - Counter(r.reactants)).elements()
        a, b = Counter(r.reactants).elements()
        bound = states(net, a) ** 2 if a == b else 2 * states(net, a) * states(net, b)
        assert states(net, new) <= bound                     # only the part reachable from the pointer is kept
    p4 = next(r for r in net.reactions if Counter(r.reactants) == Counter({m45: 2}))
    assert states(net, (Counter(p4.products) - Counter(p4.reactants)).most_common(1)[0][0]) == 4
    # s + s -> 2 s + s o s: the machine reads itself
    assert Counter({m45: 2}) in [Counter(r.reactants) for r in net.reactions]
    assert net.extras["analysis"]["elastic"]["too_large"] > 0


def test_max_states_makes_large_products_elastic():
    net = generate_network(ID, max_states=4)
    assert len(net.species) == 2 and len(net.reactions) == 1
    assert net.extras["analysis"]["elastic"]["too_large"] == 3


def test_closure_budget_truncates():
    net = generate_network(ID, machines=[], M=5, max_states=16, max_species=30, seed=4)
    assert net.status == "truncated" and len(net.species) == 30


def test_soup():
    traj = evolve(ID, machines=[], M=50, seed=3)
    net = traj.network
    assert net.status == "observed" and net.outflow == "constant-total"
    assert traj.times()[:2] == [0.0, 50.0] and traj.times()[-1] == 2000.0
    assert sum(net.initial_state.values()) == 50 == sum(net.extras["final_state"].values())
    assert all(r.count >= 1 and r.catalysts == r.reactants for r in net.reactions)
    given = evolve(ID, machines=["M45", "M61"], M=10, collisions=50, seed=0).network
    assert given.initial_state == {to_id(parse("M45")): 5.0, to_id(parse("M61")): 5.0}


def test_bad_parameters():
    with pytest.raises(ValueError, match="machine"):
        generate_network(ID, machines=["0/0,1/0;7/1,0/1"])
    with pytest.raises(ValueError, match="link"):
        generate_network(ID, machines=["0-0,1/0"])
    with pytest.raises(ValueError, match="symbols"):
        generate_network(ID, machines=["M45", "0/0,1/0,0/1"])
    with pytest.raises(ValueError, match="smaller"):
        evolve(ID, machines=["M45", "M61", "M60"], M=2)
