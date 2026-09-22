"""High-order chemistry: book appendix 'A High-Order Chemistry' (figure 3, divrule),
PyCellChemistry HighOrderChem.py, NumberChemHO.py and MolecularTSP.py."""

import ast
import inspect
from collections import Counter
from statistics import mean

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries import high_order_chem as hoc

ID = "high-order-chem"
RULE = "rule:divrule"


def value(species_id):
    return int(species_id[1:])


def without_rules(reaction):
    strip = lambda side: frozenset((s, n) for s, n in side.items() if not s.startswith("rule:"))
    return strip(reaction.reactants), strip(reaction.products)


def consumed(reaction):
    return {s for s, n in reaction.reactants.items() if reaction.products.get(s, 0) < n}


def assert_division(reaction):
    lhs = sorted(value(s) for s, n in reaction.reactants.items() if s != RULE for _ in range(n))
    rhs = sorted(value(s) for s, n in reaction.products.items() if s != RULE for _ in range(n))
    s1, s2 = lhs
    assert s1 < s2 and s2 % s1 == 0
    assert rhs == sorted([s1, s2 // s1])


# --- the book's example: divrule re-implements the prime number chemistry ------------
def test_divrule_closure_equals_prime_number_chemistry():
    """Appendix: 'The results should be the same as the original NumberChem implementation'."""
    mine = generate_network(ID, data=[12, 2, 3])
    assert mine.status == "complete"
    as_sets = {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in mine.reactions}
    assert (frozenset({RULE: 1, "n2": 1, "n12": 1}.items()),
            frozenset({RULE: 1, "n2": 1, "n6": 1}.items())) in as_sets
    assert (frozenset({RULE: 1, "n3": 1, "n12": 1}.items()),
            frozenset({RULE: 1, "n3": 1, "n4": 1}.items())) in as_sets
    for numbers, seed in (([12, 2, 3], None), (None, 3)):
        kwargs = dict(numbers=numbers) if numbers else {}
        ref = generate_network("prime-number-chemistry", seed=seed, **kwargs)
        kwargs = dict(data=numbers) if numbers else {}
        net = generate_network(ID, seed=seed, **kwargs)
        # same rng stream: M = 100 draws from [2, 1000] in both chemistries
        assert {s.id for s in net.species} - {RULE} == {s.id for s in ref.species}
        assert {without_rules(r) for r in net.reactions} == {without_rules(r) for r in ref.reactions}


def test_rules_are_catalysts():
    traj = evolve(ID, seed=1)
    net = traj.network
    assert net.status == "observed" and net.reactions
    assert net.species[0].id == RULE and net.species[0].structure == "divrule(m1, m2)"
    for r in net.reactions:
        assert r.catalysts.get(RULE) == 1 and r.reactants[RULE] == r.products[RULE] == 1
        assert_division(r)
    assert net.initial_state[RULE] == net.extras["final_state"][RULE] == 4
    assert net.extras["analysis"]["rule_draws"] == {RULE: 10000}
    assert all(f.state[RULE] == 4 for f in traj.frames)


@pytest.mark.slow
def test_divrule_soup_behaves_like_numberchem():
    """NumberChemHO.run (M = 100 from [2, 1000], 10000 iterations) against NumberChem:
    only composites are consumed, primes accumulate, and the final prime fraction matches."""
    mine, ref = [], []
    for seed in range(6):
        traj = evolve(ID, seed=seed)
        net = traj.network
        a = net.extras["analysis"]
        assert sum(net.initial_state.values()) - 4 == 100
        assert sum(net.extras["final_state"].values()) - 4 == 100     # divrule returns 2 for 2
        assert a["effective_collisions"] == sum(r.count for r in net.reactions) > 0
        assert a["idle_draws"] == 0
        for r in net.reactions:
            assert_division(r)
            assert not any(hoc.is_prime(value(s)) for s in consumed(r))
        fraction = traj.series("prime_fraction")
        assert len(fraction) == 10000 // 100 + 1
        assert all(b >= a_ for a_, b in zip(fraction, fraction[1:]))   # primes are never consumed
        assert fraction[-1] > fraction[0]
        mine.append(fraction[-1])
        ref.append(evolve("prime-number-chemistry", seed=seed).series("prime_fraction")[-1])
    assert mean(mine) > 0.95
    assert abs(mean(mine) - mean(ref)) < 0.05


# --- the algorithm of figure 3 -------------------------------------------------------
def test_iterate_draws_exactly_the_rule_arity():
    rules = hoc.parse_rules({"divrule(m1, m2)": 2, "one: x -> x": 1, "three: x, y, z -> x + y + z": 1})
    arity = {r.id: r.arity for r, _ in rules}
    assert arity == {RULE: 2, "rule:one": 1, "rule:three": 3}
    chem = hoc.HighOrderChem(rules, list(range(2, 60)), np.random.default_rng(0))
    rset = Counter(r.id for r in chem.rset)
    seen = Counter()
    for _ in range(300):
        before = len(chem.mset)
        rule, educts, products = chem.iterate()
        seen[rule.id] += 1
        if before >= rule.arity:
            assert len(educts) == rule.arity
            assert len(chem.mset) == before - rule.arity + len(products)
        else:
            assert educts == [] and products == [] and len(chem.mset) == before
        assert Counter(r.id for r in chem.rset) == rset      # the rule is reinjected
    assert set(seen) == set(arity)
    # the three-site rule shrinks the soup until it idles: too few molecules to bind
    idle = hoc.HighOrderChem(hoc.parse_rules({"x, y, z -> x": 1}), [5, 7], np.random.default_rng(0))
    assert [len(e) for _, e, _ in (idle.iterate() for _ in range(5))] == [0] * 5
    assert sorted(idle.mset) == [5, 7]


def test_binding_sites_are_counted_like_the_reference():
    assert hoc.binding_sites("self.divrule(%d,%d)") == 2           # NumberChemHO.py
    assert hoc.binding_sites("self.recombinationMachine(%s,%s)") == 2   # MolecularTSP.py
    assert hoc.binding_sites("bimolecular(m1, m2)") == 2           # book example
    assert hoc.binding_sites("self.fold(m)") == 1
    net = evolve(ID, rules={"self.divrule(%d,%d)": 1, "divrule": 1}, iterations=10, seed=0).network
    assert net.initial_state[RULE] == 2                            # both spellings are one rule


def test_rules_are_drawn_in_proportion_to_multiplicity():
    net = evolve(ID, rules={"a: x -> x": 3, "b: x -> x": 1}, iterations=4000, seed=2).network
    draws = net.extras["analysis"]["rule_draws"]
    assert draws["rule:a"] + draws["rule:b"] == 4000
    assert abs(draws["rule:a"] / 4000 - 0.75) < 0.03
    assert net.reactions == []     # identity rules are elastic


def test_rules_compete_for_substrate():
    """Destroying substrate with one rule starves the other (book: 'look at how rules compete')."""
    kept = evolve(ID, rules={"divrule": 1}, seed=5, iterations=3000).network
    eaten = evolve(ID, rules={"divrule": 1, "eat: x ->": 1}, seed=5, iterations=3000).network
    assert sum(eaten.extras["final_state"].values()) - 2 == 0
    assert eaten.extras["analysis"]["idle_draws"] > 0
    div = lambda n: sum(r.count for r in n.reactions if RULE in r.reactants)
    assert div(eaten) < div(kept)


# --- expression rules ----------------------------------------------------------------
def test_expression_rule_reproduces_divrule():
    expr = "div: x, y -> max(x, y) / min(x, y), min(x, y) if x != y and max(x, y) % min(x, y) == 0"
    a = generate_network(ID, seed=4)
    b = generate_network(ID, seed=4, rules={expr: 1})
    assert b.species[0].id == "rule:div"
    assert {s.id for s in a.species[1:]} == {s.id for s in b.species[1:]}
    assert {without_rules(r) for r in a.reactions} == {without_rules(r) for r in b.reactions}


def test_expression_semantics():
    (rule, _), = hoc.parse_rules({"x, y -> x / y, y % x if y != 0": 1})
    assert rule.id == "rule:expr1" and rule.arity == 2
    assert rule.apply(None, None, 7, 3) == [2, 3]
    assert rule.apply(None, None, 7, 0) == [7, 0]              # condition false: elastic
    assert rule.apply(None, None, 0, 3) == [0, 3]              # 3 % 0: ill-typed, elastic
    assert rule.apply(None, None, (1, 2), 3) == [(1, 2), 3]    # list in arithmetic: elastic
    (swap, _), = hoc.parse_rules({"x, y -> y, x, -x if not (x == y) or false": 1})
    assert swap.apply(None, None, 3, 4) == [4, 3, -3]
    assert swap.apply(None, None, (1, 2), 4) == [(1, 2), 4]    # -x on a list: elastic
    (make, _), = hoc.parse_rules({"x -> x, abs(0 - x), min(x, 5, 9) * 2": 1})
    assert make.apply(None, None, -7) == [-7, 7, -14]


def test_closure_truncates_on_budget():
    net = generate_network(ID, rules={"x, y -> x + y": 1}, data=[1], max_species=5)
    assert net.status == "truncated"
    assert [s.id for s in net.species] == ["rule:expr1", "n1", "n2", "n3", "n4", "n5"]


# --- no string is ever executed ------------------------------------------------------
@pytest.mark.parametrize("rule, message", [
    ("__import__('os').system('touch pwned')", "unknown rule"),
    ("os.system(m)", "unknown rule"),
    ("exec(m1, m2)", "unknown rule"),
    ("self.nprimes()", "unknown rule"),
    ("divrule(m)", "binding sites"),
    ("x -> eval(x)", "unknown function"),
    ("x -> __import__(x)", "unknown function"),
    ("x -> y", "unknown variable"),
    ("x -> x.__class__", "cannot read"),
    ("x -> 'rm -rf'", "cannot read"),
])
def test_unknown_rules_are_rejected_not_executed(rule, message, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match=message):
        evolve(ID, rules={rule: 1}, iterations=10)
    assert list(tmp_path.iterdir()) == []


def test_module_never_calls_eval_or_exec():
    tree = ast.parse(inspect.getsource(hoc))
    called = {n.func.id for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert not called & {"eval", "exec", "compile", "__import__", "getattr"}
    imported = {a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom)) for a in n.names}
    assert "builtins" not in imported


# --- MolecularTSP.py rules -----------------------------------------------------------
def test_ring_graph_follows_moleculartsp():
    g = hoc.TourGraph(10)
    assert g.coord[0] == pytest.approx((20.0, 10.0))
    assert g.toofar == pytest.approx(20 * 2 ** 0.5)
    polygon = 2 * 10 * 10 * np.sin(np.pi / 10)
    assert g.fitness(tuple(range(10))) == pytest.approx(polygon)
    d01 = g.distance(g.coord[0], g.coord[1])
    # no road from a city to itself (2 * toofar), city 1 visited twice (+ toofar)
    assert g.fitness((0, 1, 1)) == pytest.approx(2 * d01 + 3 * g.toofar)
    assert hoc.split_tour([1, 2, 3, 4, 5, 6, 7], 1, 3) == ([2, 3], [4, 5, 6, 7, 1])
    assert hoc.split_tour([1, 2, 3, 4, 5, 6, 7], 3, 1) == ([4, 5, 6, 7, 1], [2, 3])


def test_tsp_machines_only_release_better_tours():
    rules = {"exchangeMachine": 100, "cutMachine": 100, "invertMachine": 100, "recombinationMachine": 1}
    net = evolve(ID, rules=rules, init="tours", M=9, cities=10, iterations=3000, seed=3).network
    g = hoc.TourGraph(10)
    tour = lambda s: tuple(int(c) for c in s.strip("[]").split(","))
    for r in net.reactions:
        machine = next(s for s in r.reactants if s.startswith("rule:"))
        educts = [tour(s) for s, n in r.reactants.items() if s != machine for _ in range(n)]
        products = [tour(s) for s, n in r.products.items() if s != machine for _ in range(n)]
        assert all(g.is_tour(t) for t in products)
        if machine == "rule:recombinationMachine":
            assert len(educts) == len(products) == 2
        else:
            assert len(educts) == len(products) == 1
            assert g.fitness(products[0]) < g.stored_fitness(educts[0])
    lengths = lambda state: [g.fitness(tour(s)) for s, c in state.items() if not s.startswith("rule:")
                             for _ in range(c)]
    start, end = lengths(net.initial_state), lengths(net.extras["final_state"])
    assert len(start) == len(end) == 9
    assert min(end) < min(start) and mean(end) < mean(start)


# --- parameter errors ----------------------------------------------------------------
@pytest.mark.parametrize("given, message", [
    (dict(minn=500, maxn=100), "maxn must be >= minn"),
    (dict(data=[4, True]), "integers or lists"),
    (dict(data=["12"]), "integers or lists"),
    (dict(rules={}), "non-empty"),
    (dict(rules={"divrule": 0}), "positive integer"),
    (dict(rules={"divrule": 1, "divrule: x, y -> x": 1}), "duplicate"),
    (dict(rules={"cutMachine": 1}), "deterministic"),
    (dict(rules={"cutMachine": 1}, data=[[0, 1, 1]], cities=3), "permutations"),
])
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
