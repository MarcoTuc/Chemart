"""Proof AC (RESAC): book 16.5; Busch 2004 thesis, Def. 1.2.12-1.2.13, ch. 3, ch. 6."""

import itertools

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.proof_ac import (
    EMPTY, PROBLEMS, Chemistry, binary_resolvents, canonical, factors, parse_clauses,
)

ID = "proof-ac"


def cid(text):
    return canonical(parse_clauses(text)[0])[0]


def new_products(reaction):
    return [s for s in reaction.products if reaction.products[s] > reaction.reactants.get(s, 0)]


# --- the reaction rule --------------------------------------------------------------
def test_resolution_example_1_2_1():
    """Thesis example 1.2.1: {~owns(john,X), ~cat(X)} and {cat(pussy)} give ~owns(john,pussy).
    (The extracted text lost the negation signs; the prose 'John does not own Pussy' fixes them.)"""
    [a, b] = parse_clauses("~owns(john,X), ~cat(X); cat(pussy)")
    assert [canonical(r)[0] for r in binary_resolvents(a, b)] == ["~owns(john,pussy)"]


def test_factoring_example():
    """Thesis Def. 1.2.13 example: {P(f(Y)), P(X), R(X,Y)} has the factor {P(f(Y)), R(f(Y),Y)}."""
    [c] = parse_clauses("p(f(Y)), p(X), r(X,Y)")
    assert [canonical(f)[0] for f in factors(c)] == [cid("p(f(Y)), r(f(Y),Y)")]


def test_occurs_check_and_standardising_apart():
    [a, b] = parse_clauses("p(X); ~p(f(X))")
    # variables of different clauses are distinct, so X unifies with f(X') ...
    assert [canonical(r)[0] for r in binary_resolvents(a, b)] == [EMPTY]
    # ... but a term never unifies with a term containing the same variable
    [c] = parse_clauses("~p(X), ~q(X,f(X))")
    [d] = parse_clauses("p(Y), q(Y,Y)")
    chem = Chemistry(factoring=True, max_length=0, strategy="unrestricted")
    outs = set(chem.outcomes(chem.add(c), chem.add(d)))
    assert EMPTY not in outs


def test_alphabetic_variants_are_one_species():
    assert cid("~has(X,Y), mice(Y)") == cid("mice(B), ~has(A,B)")
    assert cid("p(X,Y), p(Y,X)") == cid("p(B,A), p(A,B)")
    assert cid("p(X), p(X)") == cid("p(Z)")          # identical literals merge
    assert cid("p(X), q(Y)") != cid("p(X), q(X)")


def test_published_proof_of_table_6_1():
    """Thesis table 6.1 (Gini 1995): every step of the refutation found by RESAC is a
    resolvent of its two premises, and all its clauses belong to the closure."""
    axioms, goal = PROBLEMS["gini-1995"]
    s = {f"s{i + 1}": c for i, c in enumerate(parse_clauses(axioms + ";" + goal))}
    steps = [
        ("cat(s), howls(s)", "s5", "s1"),
        ("cat(s), ~has(X,s), ~lightsleep(X)", "1", "s3"),
        ("~has(X,s), ~has(Y,Z), ~has(Y,s), ~lightsleep(X), ~mice(Z)", "2", "s2"),
        ("~has(X,Y), ~has(X,s), ~has(john,s), ~mice(Y)", "3", "s6"),
        ("~has(X,m), ~has(X,s), ~has(john,s)", "4", "s8"),
        ("~has(john,s)", "5", "s7"),
        ("[]", "6", "s4"),
    ]
    chem = Chemistry(factoring=True, max_length=0, strategy="unrestricted")
    ids = {k: chem.add(v) for k, v in s.items()}
    for n, (text, a, b) in enumerate(steps, 1):
        ids[str(n)] = cid(text)
        assert ids[str(n)] in chem.outcomes(ids[a], ids[b]), f"step {n}"
    # the clauses first grow to 5 literals, then shrink to the empty clause
    assert [len(parse_clauses(t)[0]) for t, _, _ in steps] == [2, 3, 5, 4, 3, 1, 0]

    net = generate_network(ID)
    assert net.status == "complete"
    species = {sp.id for sp in net.species}
    assert {ids[str(n)] for n in range(1, 8)} <= species
    lines = set(net.to_text().splitlines())
    assert f"{ids['s1']} + {ids['s5']} -> {ids['s1']} + {ids['s5']} + {ids['1']}" in lines or \
        f"{ids['s5']} + {ids['s1']} -> {ids['s5']} + {ids['s1']} + {ids['1']}" in lines


def test_default_closure_proves_table_6_1():
    net = generate_network(ID)
    a = net.extras["analysis"]
    assert net.status == "complete" and a["proved"] and a["target"] == EMPTY
    start = set(net.extras["axioms"]) | set(net.extras["goal"])
    assert len(start) == 8
    # the proof is a derivation: every premise is a start clause or an earlier product
    have = set(start)
    for step in a["proof"]:
        assert set(step["reactants"]) <= have
        have.add(step["product"])
    assert a["proof"][-1]["product"] == EMPTY
    chem = Chemistry(factoring=True, max_length=0, strategy="unrestricted")
    for sp in net.species:
        chem.add(parse_clauses(sp.id)[0])
    for step in a["proof"]:
        assert step["product"] in chem.outcomes(*step["reactants"])


# --- soundness and completeness ------------------------------------------------------
def _satisfiable(clauses, atoms):
    for values in itertools.product([False, True], repeat=len(atoms)):
        env = dict(zip(atoms, values))
        if all(any(env[n] == pos for pos, n, _ in c) for c in clauses):
            return True
    return False


def test_sound_and_refutation_complete_on_random_propositional_sets():
    """Resolution is sound, and the closure contains [] exactly for unsatisfiable sets
    (thesis 1.2.2, properties 1 and 3); checked against truth tables."""
    rng = np.random.default_rng(0)
    atoms = ["p", "q", "r", "s"]
    unsat = 0
    for _ in range(60):
        text = "; ".join(
            ", ".join(("~" if rng.random() < 0.5 else "") + a
                      for a in rng.choice(atoms, size=int(rng.integers(1, 4)), replace=False))
            for _ in range(int(rng.integers(3, 9)))
        )
        net = generate_network(ID, problem="custom", clauses=text, max_species=5000)
        sat = _satisfiable(parse_clauses(text), atoms)
        unsat += not sat
        assert net.status == "complete"
        assert net.extras["analysis"]["proved"] == (not sat), text
        struct = {sp.id: parse_clauses(sp.id)[0] for sp in net.species}
        for r in net.reactions:
            premises = [struct[s] for s in r.reactants]
            for s in new_products(r):
                # every model of the premises satisfies the product
                for values in itertools.product([False, True], repeat=len(atoms)):
                    env = dict(zip(atoms, values))
                    if all(any(env[n] == pos for pos, n, _ in c) for c in premises):
                        assert any(env[n] == pos for pos, n, _ in struct[s]), r.to_text()
    assert unsat > 5


def test_classic_unsatisfiable_set_needs_merging():
    """{p,q} {~p,q} {p,~q} {~p,~q}: every refutation merges identical literals."""
    net = generate_network(ID, problem="custom", clauses="p, q; ~p, q; p, ~q; ~p, ~q")
    assert net.status == "complete" and net.extras["analysis"]["proved"]
    assert {"p", "q", "~p", "~q"} <= {sp.id for sp in net.species}


def test_no_proof_without_the_negated_goal():
    """The theory of table 6.1 alone is consistent: its closure is finite and has no []."""
    axioms, _ = PROBLEMS["gini-1995"]
    net = generate_network(ID, problem="custom", clauses=axioms)
    assert net.status == "complete"
    assert not net.extras["analysis"]["proved"] and net.extras["analysis"]["proof"] == []
    assert EMPTY not in {sp.id for sp in net.species}


def test_first_order_unsatisfiable_set_needs_factoring():
    """{p(X),p(Y)} and {~p(U),~p(V)}: without factoring only 2-literal resolvents appear
    (resolution alone is not refutation-complete, thesis 1.2.2)."""
    text = "p(X), p(Y); ~p(U), ~p(V)"
    assert generate_network(ID, problem="custom", clauses=text).extras["analysis"]["proved"]
    assert not generate_network(ID, problem="custom", clauses=text, factoring=False).extras["analysis"]["proved"]


def test_every_closure_reaction_keeps_its_premises():
    net = generate_network(ID)
    for r in net.reactions:
        assert all(r.products[s] >= n for s, n in r.reactants.items())
        assert len(new_products(r)) == 1 and sum(r.reactants.values()) == 2


def test_strategies_restrict_but_still_refute():
    """Thesis 1.3.3: negative resolution and set of support keep refutation-completeness."""
    full = generate_network(ID)
    sizes = {}
    for strategy in ("set-of-support", "negative", "negative-set-of-support"):
        net = generate_network(ID, strategy=strategy)
        assert net.extras["analysis"]["proved"], strategy
        assert {sp.id for sp in net.species} <= {sp.id for sp in full.species}
        sizes[strategy] = len(net.species)
    struct = {sp.id: parse_clauses(sp.id)[0] for sp in full.species}
    net = generate_network(ID, strategy="negative")
    for r in net.reactions:
        assert any(all(not p for p, _, _ in struct[s]) for s in r.reactants)
    # the published proof starts with s5 + s1, neither negative nor supported
    assert cid("cat(s), howls(s)") not in {sp.id for sp in generate_network(ID, strategy="set-of-support").species}
    assert max(sizes.values()) < len(full.species)


def test_max_length_makes_long_clauses_unstable():
    net = generate_network(ID, max_length=6)
    assert all(sum(1 for _ in parse_clauses(sp.id)[0]) <= 3 for sp in net.species if sp.id not in net.extras["axioms"])
    # the published proof needs a 5-literal clause; shorter refutations may remain
    assert cid("~has(X,s), ~has(Y,Z), ~has(Y,s), ~lightsleep(X), ~mice(Z)") not in {sp.id for sp in net.species}


def test_group_theory_closure_is_infinite():
    """Table 6.2 has function symbols: its resolution closure never ends, so it is truncated."""
    net = generate_network(ID, problem="group-right-inverse", max_species=100)
    assert len(net.extras["axioms"]) == 5 and len(net.extras["goal"]) == 1
    assert net.status == "truncated" and len(net.species) == 100
    assert cid("p(X,Y,f(X,Y))") in net.extras["axioms"]


def test_closure_truncates_on_budget():
    net = generate_network(ID, max_species=20)
    assert net.status == "truncated" and len(net.species) == 20
    assert not net.extras["analysis"]["proved"]


# --- the reactor ------------------------------------------------------------------------
def test_soup_follows_algorithms_3_1_and_3_3():
    net = generate_network(ID, method="soup", seed=3)
    a = net.extras["analysis"]
    assert net.status == "observed"
    assert net.initial_state == {c: 20.0 for c in net.extras["axioms"] + net.extras["goal"]}
    assert sum(net.extras["final_state"].values()) == 160           # educt replacement: size fixed
    assert a["productive_collisions"] == sum(r.count for r in net.reactions)
    assert a["inflows"] == a["collisions"] - a["productive_collisions"]   # elastic inflow
    for r in net.reactions:                                            # s1 + s2 -> s_i + s3
        assert sum(r.reactants.values()) == 2 and sum(r.products.values()) == 2
    if a["proved"]:
        assert a["collisions"] == a["collisions_to_proof"]             # stops at the proof
        have = set(net.initial_state)
        for step in a["proof"]:
            assert set(step["reactants"]) <= have
            have.add(step["product"])


def test_soup_proof_times_fig_6_3():
    """Thesis fig. 6.3: table 6.1, elastic inflow, multiplicity 20: most runs need at most
    8500 collisions (the first class; the class mean is 0.36)."""
    times = [generate_network(ID, method="soup", seed=s).extras["analysis"]["collisions_to_proof"]
             for s in range(20)]
    assert sum(t is not None for t in times) >= 18
    assert sum(t is not None and t <= 8500 for t in times) >= 12


def test_reactor_size_matters_figs_6_1_6_2():
    """Thesis figs. 6.1-6.2: a too-small reactor can prevent the proof."""
    tiny = [generate_network(ID, method="soup", seed=s, multiplicity=1, max_collisions=8500)
            .extras["analysis"]["proved"] for s in range(12)]
    normal = [generate_network(ID, method="soup", seed=s, max_collisions=8500)
              .extras["analysis"]["proved"] for s in range(12)]
    assert sum(tiny) < sum(normal)


def test_free_replacement_and_rate_inflow():
    net = generate_network(ID, method="soup", seed=1, replacement="free")
    assert sum(net.extras["final_state"].values()) == 160
    for r in net.reactions:
        assert sum(r.reactants.values()) == 2 and sum(r.products.values()) in (2, 3)
    net = generate_network(ID, method="soup", seed=1, elastic_inflow=False, inflow_rate=0.5, max_collisions=40)
    a = net.extras["analysis"]
    assert a["inflows"] == a["collisions"] // 2


@pytest.mark.parametrize(
    "given, message",
    [
        (dict(clauses="p"), "only read when problem is custom"),
        (dict(problem="custom"), "needs clauses"),
        (dict(problem="custom", clauses="p(X"), "cannot parse"),
        (dict(target="p; q"), "exactly one clause"),
        (dict(problem="custom", clauses="p; ~p", strategy="set-of-support"), "need goal clauses"),
        (dict(method="soup", inflow_rate=0.5), "elastic_inflow=false"),
        (dict(method="soup", problem="custom", clauses="p", multiplicity=1), "at least 2 molecules"),
    ],
)
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)
