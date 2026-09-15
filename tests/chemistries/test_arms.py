"""ARMS reproduces the book's example and the published behaviour of [831].

Book: Banzhaf & Yamamoto (2015), section 9.4, eqs. 9.13-9.18.
[831]: Suzuki & Tanaka, "Order parameter for a symbolic chemical system",
Artificial Life VI (1998) 130-139: fig. 3 (Ru1 walk-through), fig. 8
(Brusselator rules), eq. 3 (lambda_e), 'rule set' and figs. 15-17 (termination
and kinds of periods against the heating probability p).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.arms import BOOK_RULES, RU1_RULES, Rule, run, two_symbol_rules

ID = "arms"


def rules(texts):
    return [Rule(t, f"r{i + 1}") for i, t in enumerate(texts)]


def stoich(net):
    return [(r.reactants, r.products) for r in net.reactions]


# --- book 9.4 --------------------------------------------------------------------------
def test_book_network_is_eq_9_15():
    net = generate_network(ID)
    assert net.status == "complete"
    assert {s.id for s in net.species} == set("abcdefh")                 # eq. 9.14
    assert stoich(net) == [                                                # eq. 9.15
        ({"a": 3}, {"c": 1}), ({"b": 1}, {"d": 1}), ({"c": 1}, {"e": 1}),
        ({"d": 1}, {"f": 2}), ({"a": 1}, {"a": 2, "b": 2}), ({"f": 1}, {"h": 1}),
    ]
    assert all(r.rate is None for r in net.reactions)
    assert net.initial_state == {"a": 3.0, "b": 1.0}                     # eq. 9.16
    assert net.extras["rule_kinds"] == ["cooling", "neutral", "neutral", "heating", "heating", "neutral"]


def test_index_order_halts_as_eq_9_17():
    out = run(rules(BOOK_RULES), {"a": 3, "b": 1}, selection="order", steps=50,
              rng=np.random.default_rng(0), trace=True)
    assert out["trace"] == [{"a": 3, "b": 1}, {"c": 1, "d": 1}, {"e": 1, "f": 2}, {"e": 1, "h": 2}]
    analysis = generate_network(ID).extras["analysis"]
    assert analysis["class"] == "halting" and analysis["steps"] == 3
    assert analysis["final_state"] == {"e": 1, "h": 2}


def test_reverse_order_expands_as_eq_9_18():
    out = run(rules(BOOK_RULES), {"a": 3, "b": 1}, selection="reverse-order", steps=12,
              rng=np.random.default_rng(0), trace=True)
    assert out["trace"][1] == {"a": 6, "b": 6, "d": 1}                    # {a,b,b,a,d,a,b,b,a,a,b,b,a}
    sizes = [sum(s.values()) for s in out["trace"]]
    assert all(later > earlier for earlier, later in zip(sizes, sizes[1:]))
    assert out["class"] == "non-recurrent" and not out["halted"]
    net = generate_network(ID, selection="reverse-order", steps=30)
    assert net.extras["analysis"]["final_size"] > 2 ** 30


# --- [831] -----------------------------------------------------------------------------
def test_ru1_walk_through_of_figure_3():
    # maximal size 4, input a, rule order r4 => r1 => r3 => r2
    out = run(rules(RU1_RULES), {"a": 4}, selection="random", steps=10, rng=np.random.default_rng(0),
              inputs=["a"], max_size=4, sequence=[3, 0, 2, 1], trace=True)
    assert out["trace"][:5] == [
        {"a": 4},              # initial
        {"a": 4},              # step 1: cannot input a, cannot apply r4 (size would be 5)
        {"a": 1, "b": 1},      # step 2: cannot input a, r1 gives {b, a}
        {"a": 2, "c": 1},      # step 3: input a, r3 gives {c, a, a}
        {"a": 3, "c": 1},      # step 4: input a, r2 cannot apply
    ]


def test_ru1_and_brusselator_rule_sets():
    ru1 = generate_network(ID, system="ru1", seed=1)
    assert stoich(ru1)[:4] == [({"a": 3}, {"b": 1}), ({"b": 1}, {"a": 1}), ({"b": 1}, {"c": 1}),
                               ({"a": 1}, {"b": 2})]
    assert ru1.extras["reaction_rules"] == ["r1", "r2", "r3", "r4", "input"]
    assert ru1.extras["max_size"] == 4
    assert ru1.extras["analysis"]["final_size"] <= 4

    bru = generate_network(ID, system="brusselator", steps=200, seed=1)   # fig. 8
    assert stoich(bru) == [({"A": 1}, {"X": 1}), ({"B": 1, "X": 1}, {"Y": 1, "D": 1}),
                           ({"X": 2, "Y": 1}, {"X": 3}), ({"X": 1}, {"E": 1}),
                           ({}, {"A": 1}), ({}, {"B": 1})]
    assert bru.extras["max_size"] == 5000 and not bru.initial_state


def test_two_symbol_rule_set():
    texts = two_symbol_rules()
    assert len(texts) == len(set(texts)) == 380
    kinds = Counter(r.kind for r in rules(texts))
    assert kinds == {"heating": 155, "cooling": 155, "neutral": 70}
    for r in rules(texts):
        assert set(r.lhs) | set(r.rhs) <= {"a", "b"} and r.lhs != r.rhs
        assert 1 <= sum(r.lhs.values()) <= 5 and 1 <= sum(r.rhs.values()) <= 5
    net = generate_network(ID, system="two-symbol", rule_count=12, seed=4)
    assert len(net.reactions) == 12 and 1 <= sum(net.initial_state.values()) <= 10
    assert net.to_dict() == generate_network(ID, system="two-symbol", rule_count=12, seed=4).to_dict()


def test_termination_at_extreme_heating_probability():
    """[831] fig. 15 and eq. 3: p = 0 and p = 1 terminate; p = 0.5 does not and lambda_e is near 1."""
    for seed in range(5):
        cold = generate_network(ID, system="two-symbol", p=0.0, seed=seed).extras["analysis"]
        assert cold["halted"] and cold["lambda_e"] == 0.0 and cold["final_size"] == 1
        hot = generate_network(ID, system="two-symbol", p=1.0, seed=seed).extras["analysis"]
        assert hot["halted"] and hot["cooling_uses"] == 0 and hot["lambda_e"] is None
        assert hot["final_size"] == 10
        mid = generate_network(ID, system="two-symbol", p=0.5, seed=seed).extras["analysis"]
        assert not mid["halted"] and mid["class"] == "recurrent" and mid["steps"] == 1000
        assert 0.8 < mid["lambda_e"] < 1.25


@pytest.mark.slow
def test_kinds_of_periods_peak_near_equal_heating_and_cooling():
    """[831] figs. 16-17: cycles for intermediate p, most kinds of periods near p = 0.5.

    Statistic: kinds_of_periods, the number of distinct return times to a multiset
    in a run of 1000 rewriting steps, averaged over 20 random initial states.
    """
    def mean_kinds(p):
        runs = [generate_network(ID, system="two-symbol", p=p, seed=s).extras["analysis"] for s in range(20)]
        assert not any(a["halted"] for a in runs)
        return np.mean([a["kinds_of_periods"] for a in runs])

    low, mid, high = mean_kinds(0.1), mean_kinds(0.5), mean_kinds(0.9)
    assert mid > low and mid > high


# --- properties and parameters -----------------------------------------------------------
def test_custom_rules_equal_the_book_example():
    custom = generate_network(ID, system="custom", rules=BOOK_RULES, initial={"a": 3, "b": 1},
                              selection="order")
    book = generate_network(ID)
    assert stoich(custom) == stoich(book)
    assert custom.extras["analysis"] == book.extras["analysis"]


def test_deterministic_ordered_run_detects_its_period():
    net = generate_network(ID, system="custom", rules=["a -> b", "b -> a + c", "c -> "],
                           initial={"a": 1}, selection="order")
    analysis = net.extras["analysis"]
    # {a} -> {b} -> {a, c} -> {b} : back to {b} after 2 steps
    assert analysis["class"] == "periodic" and analysis["period"] == 2


def test_bad_parameters():
    with pytest.raises(ValueError, match="only used with system='custom'"):
        generate_network(ID, rules=["a -> b"])
    with pytest.raises(ValueError, match="rule_count"):
        generate_network(ID, rule_count=3)
    with pytest.raises(ValueError, match="needs rules"):
        generate_network(ID, system="custom")
    with pytest.raises(ValueError, match="exactly one '->'"):
        generate_network(ID, system="custom", rules=["a b"], initial={"a": 1})
    with pytest.raises(ValueError, match="max_size"):
        generate_network(ID, system="custom", rules=["a -> b"], initial={"a": 5}, max_size=2)
