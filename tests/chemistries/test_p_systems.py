"""P systems: book 9.5 (divisibility, eq. 9.21 / fig. 9.3) and Paun (2006), Introduction to Membrane Computing."""

from math import isqrt

import pytest

from chemart import generate_network


def custom(**kw):
    return generate_network("p-systems", system="custom", **kw)


def texts(net):
    return sorted(r.to_text() for r in net.reactions)


# --- book 9.5: is n a multiple of k? ----------------------------------------------------------
def test_book_divisibility_network():
    """Eq. 9.21: R1 = {r4: dcc' -> a_in3}, R2 = {r1: ac -> c', r2: ac' -> c, r3: d -> d delta}."""
    net = generate_network("p-systems")
    assert net.status == "complete"
    assert texts(net) == sorted([
        "d@1 + c@1 + c'@1 -> a@3",
        "a@2 + c@2 -> c'@2",
        "a@2 + c'@2 -> c@2",
        "d@2 -> d@2",
    ])
    assert net.initial_state == {"a@2": 7.0, "c@2": 3.0, "d@2": 1.0}
    assert net.extras["priorities"] == {"2": [["r1", "r3"], ["r2", "r3"]]}
    assert net.extras["dissolution"]["2"]["rules"] == ["2/r3"]
    assert net.extras["dissolution"]["2"]["contents_move_to"] == ["1"]
    assert net.extras["compartments"]["2"]["parent"] == "1"
    assert net.extras["compartments"]["1"]["children"] == ["2", "3"]
    assert set(net.provides) >= {"compartments", "initial-state"}


def test_book_divisibility_decides_multiples():
    """Book 9.5: membrane 3 stays empty exactly when n is a multiple of k (also for n < k)."""
    for n in range(1, 13):
        for k in range(1, 6):
            a = generate_network("p-systems", n=n, k=k).extras["analysis"]
            assert a["halted"] and not a["output_dissolved"]
            assert a["n_multiple_of_k"] == (n % k == 0), (n, k)
            assert a["output"] == ({} if n % k == 0 else {"a": 1})


def test_book_divisibility_computation_steps():
    """Book 9.5: r1 fires k times, r2 k times, ...; r3 only once a is exhausted; then r4."""
    a = generate_network("p-systems", n=7, k=3).extras["analysis"]
    assert [s["applied"] for s in a["trace"]] == [
        {"2/r1": 3}, {"2/r2": 3}, {"2/r1": 1}, {"2/r3": 1}, {"1/r4": 1},
    ]
    assert a["trace"][3]["dissolved"] == ["2"]
    assert a["steps"] == 5 and a["final_structure"] == "[1 [3 ]3 ]1"
    assert a["final_configuration"] == {"1": {"c": 1}, "3": {"a": 1}, "env": {}}
    small = generate_network("p-systems", n=2, k=5).extras["analysis"]     # n < k: a mix of c and c'
    assert [s["applied"] for s in small["trace"]] == [{"2/r1": 2}, {"2/r3": 1}, {"1/r4": 1}]


def test_book_system_needs_strong_priorities():
    """Under weak priorities r3 fires as soon as r1 has used up the c, and the answer is wrong."""
    a = generate_network("p-systems", n=7, k=3, priority="weak").extras["analysis"]
    assert a["trace"][0]["applied"] == {"2/r1": 3, "2/r3": 1}
    assert a["n_multiple_of_k"] is True                       # 7 is not a multiple of 3


# --- Paun 2006 ----------------------------------------------------------------------------------
def test_paun_maximal_parallelism_example():
    """Sec. 4: a^5 b^2 c^6 with aab -> abcc and bb -> aac gives a^3 b^2 c^10 or a^7 c^7."""
    outcomes = set()
    for seed in range(40):
        net = custom(membranes="[1 ]1", objects={"1": "5 a + 2 b + 6 c"},
                     rules={"1": ["2 a + b -> a + b + 2 c", "2 b -> 2 a + c"]},
                     output="1", max_steps=1, seed=seed)
        outcomes.add(tuple(sorted(net.extras["analysis"]["final_configuration"]["1"].items())))
    assert outcomes == {(("a", 3), ("b", 2), ("c", 10)), (("a", 7), ("c", 7))}


def test_paun_target_indications():
    """Sec. 7: aab -> (a,here)(b,out)(c,here)(c,in) is used twice on a^5 b^2 c^6 in the skin."""
    net = custom(membranes="[1 [2 ]2 ]1", objects={"1": {"a": 5, "b": 2, "c": 6}},
                 rules={"1": ["2 a + b -> a + b@out + c + c@in"]}, output="0")
    assert texts(net) == ["2 a@1 + b@1 -> a@1 + b@env + c@1 + c@2"]
    a = net.extras["analysis"]
    assert a["halted"] and a["steps"] == 1 and a["trace"][0]["applied"] == {"1/r1": 2}
    assert a["final_configuration"] == {"1": {"a": 3, "c": 8}, "2": {"c": 2}, "env": {"b": 2}}
    assert a["result"] == 2


def test_paun_strong_and_weak_priorities():
    """Sec. 11: r1: ff -> f > r2: cf -> cd delta on fffc."""
    spec = dict(membranes="[1 [2 ]2 ]1", objects={"2": "3 f + c"},
                rules={"2": ["r1: 2 f -> f", "r2: c + f -> c + d + delta"]},
                priorities={"2": ["r1 > r2"]}, output="1")
    strong = custom(**spec).extras["analysis"]["trace"][0]
    assert strong["applied"] == {"2/r1": 1} and strong["dissolved"] == []
    weak = custom(priority="weak", **spec).extras["analysis"]["trace"][0]
    assert weak["applied"] == {"2/r1": 1, "2/r2": 1} and weak["dissolved"] == ["2"]


def test_paun_n_squared_system():
    """Fig. 3: every halting computation sends n^2 (n >= 1) copies of e to the environment."""
    net = generate_network("p-systems", system="n-squared", seed=0)
    assert {"a@3 -> a@3 + b@3", "a@3 -> b@3", "f@3 -> 2 f@3", "b@2 -> d@2", "d@2 -> d@2 + e@2",
            "2 f@2 -> f@2", "c@2 + f@2 -> c@2 + d@2", "e@1 -> e@env", "f@1 -> f@1"} == set(texts(net))
    assert net.extras["output_region"] == "env"
    results, stuck = set(), 0
    for seed in range(150):
        a = generate_network("p-systems", system="n-squared", max_steps=60, seed=seed).extras["analysis"]
        if a["halted"]:
            assert a["result_is_square"] and a["result"] >= 1
            assert a["output"] == {"e": a["result"]}
            n = isqrt(a["result"])
            # n-1 growth steps, dissolve 3, b -> d, n-1 halvings, cf -> cd delta, e out
            assert a["steps"] == 2 * n + 2
            results.add(a["result"])
        else:
            assert a["result"] is None
            stuck += 1
    assert {1, 4} <= results and stuck > 0


# --- semantics details ---------------------------------------------------------------------------
def test_simultaneous_dissolution_goes_to_first_surviving_membrane():
    """Paun 2006 sec. 7: contents stop in the first upper membrane that is not dissolved."""
    net = custom(membranes="[1 [2 [3 ]3 ]2 ]1", objects={"2": "x", "3": "y"},
                 rules={"2": ["x -> x + delta"], "3": ["y -> z + delta"]}, output="1")
    a = net.extras["analysis"]
    assert a["trace"][0]["dissolved"] == ["2", "3"]
    assert a["final_structure"] == "[1 ]1"
    assert a["final_configuration"] == {"1": {"x": 1, "z": 1}, "env": {}}
    assert a["halted"] and a["result"] == 2


def test_network_covers_structures_reachable_by_dissolution():
    net = custom(membranes="[1 [2 [3 ]3 ]2 ]1",
                 rules={"3": ["y -> y@out"], "2": ["x -> delta"], "1": ["q -> q@in"]}, output="1")
    assert texts(net) == sorted(["y@3 -> y@2", "y@3 -> y@1", "x@2 -> ∅", "q@1 -> q@2", "q@1 -> q@3"])
    both = custom(membranes="[1 [2 ]2 [3 ]3 ]1", rules={"1": ["a -> 2 b@in"]}, output="1")
    assert texts(both) == sorted(["a@1 -> 2 b@2", "a@1 -> b@2 + b@3", "a@1 -> 2 b@3"])
    lost = custom(membranes="[1 [2 ]2 ]1", rules={"2": ["a -> b@in"]}, output="1")
    assert lost.reactions == [] and lost.extras["inapplicable_rules"] == ["2/r1"]


def test_dissolved_output_membrane_gives_no_result():
    a = custom(membranes="[1 [2 ]2 ]1", objects={"2": "a"}, rules={"2": ["a -> b + delta"]},
               output="2").extras["analysis"]
    assert a["halted"] and a["output_dissolved"] and a["result"] is None


@pytest.mark.parametrize("kw, message", [
    (dict(membranes="[1 ]1", rules={"1": ["a -> b + delta"]}, output="1"), "skin"),
    (dict(membranes="[1 [2 ]1", output="1"), "membranes"),
    (dict(membranes="[1 ]1", output="7"), "output"),
    (dict(membranes="[1 ]1", rules={"1": ["r1: a -> b", "r2: b -> a"]},
          priorities={"1": ["r1 > r2", "r2 > r1"]}, output="1"), "cycle"),
    (dict(membranes="[1 ]1", rules={"1": ["a -> b@in9"]}, output="1"), "target"),
])
def test_rejects_invalid_custom_systems(kw, message):
    with pytest.raises(ValueError, match=message):
        custom(**kw)


def test_named_systems_reject_custom_fields():
    with pytest.raises(ValueError, match="system='custom'"):
        generate_network("p-systems", rules={"1": ["a -> b"]})
