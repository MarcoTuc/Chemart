"""Gamma reproduces the published programs and their order-independent results.

Book 9.2 (eqs. 9.7-9.9); Gamma15 = Banatre, Fradet & Le Metayer, "Gamma and the
chemical reaction model: fifteen years after", LNCS 2235 (2001), secs. 1, 2.1, 3.1;
RULE04 = Banatre, Fradet & Radenac, "Principles of chemical programming" (2004),
fig. 1; York = the same authors' "Higher-order chemical model of computation".
"""

from collections import Counter
from itertools import combinations

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.gamma import Rule, parse_program, run

ID = "gamma"


def analysis(net):
    return net.extras["analysis"]


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


# --- the reaction language ------------------------------------------------------------
def test_book_and_chemart_notation_agree():
    book = Rule("max : x, y → y ⇐ x ≤ y")                          # book eq. 9.7
    text = Rule("max: x, y -> y if x <= y")
    for pair in [(3, 8), (8, 3), (5, 5)]:
        assert book.apply(*pair) == text.apply(*pair)
    assert text.apply(3, 8) == (8,) and text.apply(8, 3) is None
    sort = Rule("sort: (i, x), (j, y) -> (i, y), (j, x) if i > j and x < y")
    assert sort.apply((2, 1), (1, 5)) == ((2, 5), (1, 1))
    assert sort.apply((1, 1), (2, 5)) is None
    assert Rule("rem: x, y -> y if multiple(x, y)").apply(12, 4) == (4,)
    assert Rule("maj: x, y -> if x != y").apply(1, 2) == ()
    iota = Rule("(x, y) -> (x, [(x + y) / 2]), ([(x + y) / 2] + 1, y) if x != y")
    assert iota.apply((2, 9)) == ((2, 5), (6, 9))


def test_ill_typed_tuples_do_not_react():
    rule = Rule("x, y -> y if x <= y")
    assert rule.apply((1, 2), 3) is None                               # tuple compared with <
    assert Rule("x, y -> x / y").apply(4, 0) is None                   # division by zero
    assert Rule("x, x -> x").apply(2, 2) == (2,) and Rule("x, x -> x").apply(2, 3) is None


@pytest.mark.parametrize("rules, message", [
    ("x -> __import__('os')", "cannot read"),
    ("x -> y", "variable 'y'"),
    ("x -> exec(x)", "unknown function"),
    ("x, y", "expected '->'"),
    ("then\nx -> x", "then"),
])
def test_rules_are_parsed_never_evaluated(rules, message):
    with pytest.raises(ValueError, match=message):
        parse_program(rules)


# --- max (book eq. 9.7) -----------------------------------------------------------------
def test_default_max_network():
    net = generate_network(ID)
    assert net.status == "complete"
    assert sorted(s.id for s in net.species) == ["n1", "n3", "n4", "n8"]
    assert net.initial_state == {"n1": 1.0, "n3": 1.0, "n4": 1.0, "n8": 2.0}
    # every ground instance x + y -> y with x <= y, including 2 x -> x
    assert len(net.reactions) == 6 + 4
    for r in net.reactions:
        values = sorted(int(s[1:]) for s, n in r.reactants.items() for _ in range(n))
        assert r.products == {f"n{values[1]}": 1}
    a = analysis(net)
    assert a["final_values"] == [8] and a["stable"] and a["steps"] == 4
    assert a["results"] == [{"n8": 1}] and a["deterministic"] is True


def test_max_is_independent_of_reaction_order():
    rng = np.random.default_rng(0)
    for _ in range(5):
        values = [int(v) for v in rng.integers(-20, 20, size=7)]
        finals = {tuple(analysis(generate_network(ID, multiset=values, seed=s))["final_values"])
                  for s in range(4)}
        assert finals == {(max(values),)}
        assert analysis(generate_network(ID, multiset=values))["results"] == [{f"n{max(values)}": 1}]


# --- sort (Gamma15 sec. 1 and 3.1; book eq. 9.8) ----------------------------------------------
def test_sort_orders_values_by_index():
    net = generate_network(ID, program="sort")
    a = analysis(net)
    assert a["final_values"] == [[1, 1], [2, 3], [3, 6], [4, 8], [5, 9]]
    assert a["deterministic"] is True and len(a["results"]) == 1
    # invariant (Gamma15 sec. 3.1): the index set and the multiset of values are preserved
    for r in net.reactions:
        def parts(side):
            pairs = [tuple(map(int, s.strip("()").split(","))) for s, n in side.items() for _ in range(n)]
            return sorted(i for i, _ in pairs), sorted(v for _, v in pairs)
        assert parts(r.reactants) == parts(r.products)
    assert net.status == "complete"
    assert {s.structure for s in net.species} <= {f"({i},{v})" for i in range(1, 6) for v in (8, 3, 6, 1, 9)}


def test_sort_random_sequences_from_any_order():
    rng = np.random.default_rng(1)
    for _ in range(3):
        values = [int(v) for v in rng.integers(0, 50, size=6)]
        pairs = [[i + 1, v] for i, v in enumerate(values)]
        for seed in range(3):
            net = generate_network(ID, program="sort", multiset=pairs, seed=seed, max_states=0)
            assert analysis(net)["final_values"] == [[i + 1, v] for i, v in enumerate(sorted(values))]


def test_book_printed_sort_condition_sorts_decreasingly():
    """Book eq. 9.8 as printed, (i >= j) and (x >= y), reverses the order (see decisions)."""
    net = generate_network(ID, program="custom", multiset=[[1, 8], [2, 3], [3, 6], [4, 1], [5, 9]],
                           rules="sort: (i, x), (j, y) -> (i, y), (j, x) if i >= j and x >= y")
    assert analysis(net)["final_values"] == [[1, 9], [2, 8], [3, 6], [4, 3], [5, 1]]


# --- primes (Gamma15 sec. 2.1; RULE04 fig. 1) ----------------------------------------------------
def test_primes_rem_of_iota():
    expected = [n for n in range(2, 31) if is_prime(n)]
    for seed in range(3):
        net = generate_network(ID, program="primes", multiset=[[2, 30]], seed=seed, max_states=0)
        a = analysis(net)
        assert a["stage_multisets"][0] == {f"n{k}": 1 for k in range(2, 31)}     # iota
        assert a["final_values"] == expected
    # between distinct numbers, rem consumes exactly the composites and never a prime
    rules = net.extras["reaction_rules"]
    consumed = {s for r, rule in zip(net.reactions, rules) if rule == "rem" and len(r.reactants) == 2
                for s in r.reactants if s not in r.products}
    assert {int(s[1:]) for s in consumed} == {n for n in range(4, 31) if not is_prime(n)}
    assert net.status == "complete"


def test_rem_on_2_to_n_is_deterministic():
    """RULE04 / York: primes applied to the multiset 2..N leaves the primes <= N, in any order."""
    net = generate_network(ID, program="custom", rules="rem: x, y -> y if multiple(x, y)",
                           multiset=list(range(2, 13)))
    assert analysis(net)["results"] == [{f"n{p}": 1 for p in (2, 3, 5, 7, 11)}]


def test_largest_prime_below_10():
    """York abstract: the gamma-cn molecule computes the largest prime lower than 10."""
    a = analysis(generate_network(ID, program="largest-prime"))
    assert a["stage_multisets"][0] == {"n2": 1, "n3": 1, "n5": 1, "n7": 1}
    assert a["final_values"] == [7] and a["results"] == [{"n7": 1}]


# --- fibonacci and maximum segment sum (Gamma15 sec. 2.1) ---------------------------------------
@pytest.mark.parametrize("n, fib", [(0, 1), (1, 1), (2, 2), (3, 3), (4, 5), (5, 8), (6, 13), (7, 21)])
def test_fibonacci(n, fib):
    net = generate_network(ID, program="fibonacci", multiset=[n], max_species=30, max_states=0, seed=n)
    a = analysis(net)
    assert a["stage_multisets"][0] == {"n1": fib}                  # dec1 expands into ones
    assert a["final_values"] == [fib]


def test_fibonacci_closure_is_truncated_but_result_unique():
    net = generate_network(ID, program="fibonacci", multiset=[4], max_species=40)
    assert net.status == "truncated" and len(net.species) == 40
    assert analysis(net)["results"] == [{"n5": 1}]


def brute_force_mss(xs):
    return max(sum(xs[i:j]) for i, j in combinations(range(len(xs) + 1), 2))


def test_max_segment_sum():
    a = analysis(generate_network(ID, program="max-segment-sum"))
    assert {s for _, _, s in a["final_values"]} == {6} and a["deterministic"] is True
    rng = np.random.default_rng(2)
    for _ in range(4):
        xs = [int(v) for v in rng.integers(-6, 7, size=6)]
        net = generate_network(ID, program="max-segment-sum", multiset=[[i + 1, x, x] for i, x in enumerate(xs)],
                               max_states=0, seed=3)
        assert {s for _, _, s in analysis(net)["final_values"]} == {brute_force_mss(xs)}


# --- majority (RULE04 fig. 1): the result is not unique ----------------------------------------
def test_majority_leaves_only_the_majority_element():
    a = analysis(generate_network(ID, program="majority"))
    # 1 occurs 4 times among 7; pairs of distinct elements are removed in any order
    assert a["results"] == [{"n1": 1}, {"n1": 3}] and a["deterministic"] is False
    assert set(a["final_multiset"]) == {"n1"}


# --- re-injection (book 9.2) ----------------------------------------------------------------------
def test_reinjection_restarts_the_computation():
    stages = parse_program("max: x, y -> y if x <= y")
    rng = np.random.default_rng(5)
    stable, _, done, _ = run(stages, Counter([4, 9, 2]), rng, 1000)
    assert done and stable == Counter([9])
    again, _, done, _ = run(stages, stable + Counter([11, 3]), rng, 1000)
    assert done and again == Counter([11])


def test_non_terminating_program_is_reported():
    net = generate_network(ID, program="custom", rules="x -> x + 1", multiset=[0], max_steps=50,
                           max_species=10, max_states=100)
    a = analysis(net)
    assert a["stable"] is False and a["steps"] == 50 and a["results"] is None
    assert net.status == "truncated"


# --- parameters -----------------------------------------------------------------------------------
@pytest.mark.parametrize("given, message", [
    (dict(rules="x -> x"), "only used with program='custom'"),
    (dict(program="custom", multiset=[1]), "needs rules"),
    (dict(program="custom", rules="x -> x"), "non-empty multiset"),
    (dict(multiset=[1.5]), "integers or lists"),
    (dict(multiset=[[1]]), "at least two"),
])
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
