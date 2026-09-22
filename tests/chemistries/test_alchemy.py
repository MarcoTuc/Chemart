"""AlChemy reproduces Fontana & Buss (1994), 'The arrival of the fittest', Bull. Math. Biol. 56:1-64.

Section numbers, equations and figures refer to the SFI working paper 93-09-055.
Reduction correctness is checked on standard combinators first.
"""

from collections import Counter

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import evolve, generate_network
from chemart.chemistries.alchemy import Diverged, app, lam, normal_form, parse, to_id, to_text, var

I, K = "λx.x", "λx.λy.x"
S = "λx.λy.λz.((x)z)(y)z"


def nf(text):
    return to_text(normal_form(parse(text))[0])


def same(a, b):
    return to_id(parse(a)) == to_id(parse(b))


# --- the calculus --------------------------------------------------------------
def test_combinators():
    assert nf(f"(({S}){K}){K}") == "λx1.x1"                   # S K K = I
    assert nf(f"(({K}){I}){S}") == "λx1.x1"                   # K I S = I
    assert nf(f"(({S}){K})λu.λv.v") == "λx1.x1"               # S K anything = I
    church = lambda n: "λf.λx." + "(f)" * n + "x"
    plus = "λm.λn.λf.λx.((m)f)((n)f)x"
    times = "λm.λn.λf.(m)(n)f"
    assert same(nf(f"(({plus}){church(2)}){church(3)}"), church(5))
    assert same(nf(f"(({times}){church(2)}){church(3)}"), church(6))
    assert same(nf(f"({church(2)}){church(3)}"), church(9))  # (m)n = n^m


def test_capture_avoiding_substitution():
    # (λx.λw.x) w under a binder w must not capture: λw.λw'.w, not the identity
    assert nf("λw.(λx.λw.x)w") == "λx1.λx2.x1"
    # Mathis et al. 2.1: (λx.λy.(x)y) applied to a term with a free y keeps it free
    t = app(parse("λx.λy.(x)y"), lam(app(var(0), var(1))))   # argument λx.(x)y, y free
    out, steps = normal_form(lam(t))                          # bind y outside
    assert to_text(out) == "λx1.λx2.(x2)x1" and steps == 2


def test_standardized_names_and_ids():
    # paper 4.3.3 (eq. 11): renaming gives λx1.λx2.((x1)x2)x1 (the scan shows the operand (y) in parentheses)
    t = parse("λx.λy.((x)y)x")
    assert to_text(t) == "λx1.λx2.((x1)x2)x1" and to_id(t) == "^^((2)1)2"
    assert parse("λu.λv.((u)v)u") == t == parse(to_id(t))
    # binders are numbered in order of occurrence (example 4, eq. 64)
    b = "λx1.λx2.((x2)λx3.((x1)λx4.x2)λx5.λx6.x5)x1"
    assert to_text(parse(b)) == b


def test_omega_has_no_normal_form():
    omega = "(λx.(x)x)λx.(x)x"                                # paper eq. 9
    with pytest.raises(Diverged, match="steps"):
        normal_form(parse(omega), max_steps=500)
    with pytest.raises(Diverged, match="characters"):
        normal_form(parse("(λx.((x)x)x)λx.((x)x)x"), max_size=300)
    net = generate_network("alchemy", terms=["λx.(x)x", I], max_steps=200)
    # (λx.(x)x) applied to itself is elastic; everything else closes on {ω, I}
    ids = {s.id for s in net.species}
    assert net.status == "complete" and ids == {"^(1)1", "^1"}
    assert not any(Counter(r.reactants) == Counter({"^(1)1": 2}) for r in net.reactions)
    assert net.extras["analysis"]["elastic"]["no_normal_form"] == 1


def test_parse_errors():
    with pytest.raises(ValueError, match="free variable"):
        parse("λx.y")
    with pytest.raises(ValueError, match="cannot parse"):
        parse("λx.(x")
    with pytest.raises(ValueError, match="terms"):
        generate_network("alchemy", terms=["(λx.(x)x)λx.(x)x"], max_steps=100)


# --- boundary conditions -------------------------------------------------------
def test_no_copy_filter_makes_copies_elastic():
    free = generate_network("alchemy", terms=[I, K], max_species=30)
    barred = generate_network("alchemy", terms=[I, K], filter="no-copy", max_species=30)
    copies = [r for r in free.reactions if not set(Counter(r.products) - Counter(r.reactants)) - set(r.reactants)]
    assert copies, "I o v = v is a copy action"
    for r in barred.reactions:
        (new,) = (Counter(r.products) - Counter(r.reactants)).elements()
        assert new not in r.reactants
    assert barred.extras["analysis"]["elastic"]["copy"] >= 2


def test_syntactic_filter_bars_three_abstractions():
    # paper 6.2.3: objects containing three consecutive abstractions are barred
    net = generate_network("alchemy", terms=[K], forbidden_patterns=[r"\^\^\^"])
    assert {s.id for s in net.species} == {"^^2"} and not net.reactions
    with pytest.raises(ValueError, match="regular expression"):
        generate_network("alchemy", forbidden_patterns=["("])


def test_mediator_is_the_generalized_collision_rule():
    # eq. 14: Phi = λx1.λx2.(x1)x2 reproduces plain application
    a = generate_network("alchemy", terms=[K, "λx.λy.λz.(z)x"], max_species=40)
    b = generate_network("alchemy", terms=[K, "λx.λy.λz.(z)x"], max_species=40,
                         mediator="λx1.λx2.(x1)x2")
    assert a.to_dict()["reactions"] == b.to_dict()["reactions"]
    # Phi = λa.λb.(a)(a)b lets the operator act twice: K o I = λy.λy'.I instead of λy.I
    c = generate_network("alchemy", terms=[K, I], mediator="λa.λb.(a)(a)b", max_species=20)
    made = {(Counter(r.products) - Counter(r.reactants)).most_common(1)[0][0]
            for r in c.reactions if Counter(r.reactants) == Counter(["^^2", "^1"])}
    assert "^^^1" in made and "^^1" not in made


# --- Level 0 (section 6.1) -----------------------------------------------------
FIG1 = {
    "A": "λx1.(x1)λx2.λx3.(x3)λx4.λx5.(x5)x4",
    "B": "λx1.(x1)λx2.λx3.(x3)x2",
    "C": "λx1.λx2.(x1)λx3.x1",
    "D": "λx1.λx2.λx3.(x2)λx4.x2",
    "E": "λx1.λx2.λx3.λx4.(x3)λx5.x3",
}


def product(a, b):
    """a o b in standardized notation, or None when pragmatic reduction fails (elastic)."""
    try:
        return nf(f"({a}){b}")
    except Diverged:
        return None


def copy_ecology(ids_terms):
    """eq. 33: every object is a left or right fixed point of some object (elastic pairs do not count)."""
    def fixed(x, f):
        return x is not None and same(x, f)

    def ok(f):
        return any(fixed(product(g, f), f) or fixed(product(f, g), f) for g in ids_terms)
    return all(ok(f) for f in ids_terms)


@pytest.mark.parametrize("members", ["AB", "CDE"])
def test_fig1_level0_ecologies_are_closed_copy_ecologies(members):
    terms = [FIG1[m] for m in members]
    net = generate_network("alchemy", terms=terms)
    assert net.status == "complete" and len(net.species) == len(members)
    assert copy_ecology(terms)


def test_fig1_caption_actions():
    A, B, C, D, E = (FIG1[k] for k in "ABCDE")
    assert same(product(E, C), D)                             # E o C = D
    assert same(product(B, A), A) and same(product(B, B), B)  # both are fixed points of B
    assert product(B, C) == "λx1.λx2.(x2)λx3.λx4.λx5.(x5)x4"  # B o C outside the ecology


@pytest.mark.parametrize("seed", [0, 1])
def test_level0_reactor_collapses_to_copiers(seed):
    # 6.1: 'in many instances the system reduces to just one object species that is a self-copier'.
    # Not every run does: e.g. seed 2 ends in a nearly inert pair whose collisions are almost all elastic.
    traj = evolve("alchemy", seed=seed, M=100, collisions=20000)
    net = traj.network
    final = net.extras["final_state"]
    assert len(traj.frames[0].state) == 100 and len(traj.frames[-1].state) == len(final)
    assert len(final) <= 2
    terms = [next(s.structure for s in net.species if s.id == sid) for sid in final]
    assert copy_ecology(terms)
    closed = generate_network("alchemy", terms=terms)
    assert closed.status == "complete" and len(closed.species) == len(final)


# --- Level 1 (section 6.2) -----------------------------------------------------
def projector(i, j):
    """eq. 34: A(i,j) = λx1. ... λxi.xj."""
    t = var(i - j)
    for _ in range(i):
        t = lam(t)
    return t


def projector_indices(sid):
    body = sid.lstrip("^")
    return (len(sid) - len(body), len(sid) - len(body) - int(body) + 1) if body.isdigit() else None


def test_example1_laws():
    ids = {(i, j): to_id(projector(i, j)) for i in range(1, 9) for j in range(1, i + 1)}
    for (i, j), a in ids.items():
        if i > 6:
            continue
        for (k, l), b in ids.items():
            if k > 6:
                continue
            c = to_id(normal_form(app(projector(i, j), projector(k, l)))[0])
            if j > 1:
                assert c == to_id(projector(i - 1, j - 1))                   # eq. 35
            else:
                assert c == to_id(projector(k + i - 1, l + i - 1))           # eq. 36
    assert nf("(λx1.λx2.λx3.x3)λx1.λx2.x1") == "λx1.λx2.x2"                  # the axiom-5 example


def test_example1_basic_cycle():
    # eq. 41: A10,1 o A10,1 -> A19,10, A19,10 o A19,10 -> A18,9, ..., A11,2 o A11,2 -> A10,1
    t = projector(10, 1)
    t = normal_form(app(t, t))[0]
    assert t == projector(19, 10)
    for i in range(19, 11, -1):
        t = normal_form(app(t, t))[0]
        assert t == projector(i - 1, i - 10)
    assert normal_form(app(t, t))[0] == projector(10, 1)


@pytest.mark.parametrize("family", [2, 3, 4])
def test_example1_center_is_self_maintaining(family):
    # 6.2.2 'Families': the center of family i is its i smallest objects; fewer do not seed it
    center = [to_text(projector(family + l, 1 + l)) for l in range(family)]
    net = generate_network("alchemy", terms=center, filter="no-copy", max_species=25)
    assert net.status == "truncated"
    assert all(i - j == family - 1 for i, j in map(projector_indices, (s.id for s in net.species)))
    whole = generate_network("alchemy", terms=center, filter="no-copy", max_species=family)
    assert whole.extras["analysis"]["self_maintaining"]
    part = generate_network("alchemy", terms=center[:-1], filter="no-copy",
                            max_species=family - 1)
    assert not part.extras["analysis"]["self_maintaining"]


def test_level1_reactor_settles_in_a_projector_family():
    # 6.2.2: random objects, copy actions barred -> the projector organisation O1 (one family)
    net = evolve("alchemy", seed=1, M=100, collisions=20000, filter="no-copy").network
    final = net.extras["final_state"]
    indices = [projector_indices(s) for s in final]
    assert len(final) >= 3 and None not in indices
    assert len({i - j for i, j in indices}) == 1
    for r in net.reactions:
        (new,) = (Counter(r.products) - Counter(r.reactants)).elements()
        assert new not in r.reactants


def test_level1_soup_from_the_center_stays_in_its_family():
    center = [to_text(projector(3 + l, 1 + l)) for l in range(3)]
    net = evolve("alchemy", seed=0, terms=center, M=60, collisions=3000, filter="no-copy").network
    assert net.initial_state == {to_id(parse(t)): 20 for t in center}
    assert all(projector_indices(s.id)[0] - projector_indices(s.id)[1] == 2 for s in net.species)
    assert len(net.extras["final_state"]) >= 3


# example 2 (6.2.3): a = λx.(x) prefix, A = λxj.aaxj, B = λxk.xk
def a_(v):
    return f"λa.(a){v}"


def numeral(family, i):
    base, n = ("λy.λu.(u)λw.(w)y", i + 2) if family == "A" else ("λy.y", i)
    for _ in range(n):
        base = a_(base)
    return base


def test_example2_laws():
    v1, v2 = numeral("A", 1), numeral("B", 2)
    assert same(product(a_(v1), v2), product(v2, v1))                        # eq. 50
    assert same(product("λy.λu.(u)λw.(w)y", v1), a_(a_(v1)))                 # eq. 51
    assert same(product("λy.y", v2), v2)                                     # eq. 52
    assert same(product(a_("λy.y"), numeral("A", 1)), a_(a_("λy.λu.(u)λw.(w)y")))  # eq. 53: aB o aaaA = aaA


def test_example2_numeral_arithmetic():
    # eqs. 59-61: i o j = j - i if i <= j (+0, -2, +2 across families), i - j - 1 otherwise
    lookup = {to_id(parse(numeral(f, k))): (f, k) for f in "AB" for k in range(-2 if f == "A" else 0, 14)}
    shift = {("A", "A"): 0, ("B", "B"): 0, ("A", "B"): -2, ("B", "A"): 2}
    for x in "AB":
        for y in "AB":
            for i in range(-2 if x == "A" else 0, 5):
                for j in range(-2 if y == "A" else 0, 5):
                    got = lookup[to_id(normal_form(app(parse(numeral(x, i)), parse(numeral(y, j))))[0])]
                    expected = (y, j - i) if i <= j + shift[(x, y)] else (x, i - j - 1)
                    assert got == expected, (x, i, y, j)


def test_example2_center_is_self_maintaining():
    center = [numeral("A", -2), numeral("A", 0), numeral("B", 0), numeral("B", 2)]
    net = generate_network("alchemy", terms=center, filter="no-copy",
                           forbidden_patterns=[r"\^\^\^"], max_species=4)
    assert net.extras["analysis"]["self_maintaining"]
    for sub in ([center[0], center[1], center[2]], [center[1], center[2], center[3]]):
        part = generate_network("alchemy", terms=sub, filter="no-copy",
                                forbidden_patterns=[r"\^\^\^"], max_species=3)
        assert not part.extras["analysis"]["self_maintaining"]


def test_example4_A_makes_cB():
    A = "λx1.((x1)λx2.λx3.x1)λx4.λx5.x4"                                     # eq. 63
    B = "λx1.λx2.((x2)λx3.((x1)λx4.x2)λx5.λx6.x5)x1"                        # eq. 64
    cB = "λc." + B
    assert same(product(A, B), cB) and same(product(cB, A), B)


# --- Level 2 (section 6.4.2) ---------------------------------------------------
def test_level2_organisation_A_is_built_by_T():
    T = "λa.λb.λc.λd.(d)λe.a"                                                # eq. 77
    one, zero = (lambda v: f"λq.(q){v}"), (lambda v: f"λq.{v}")
    cycle = [T, zero(zero(one(zero(T)))), zero(one(zero(T))), one(zero(T)), T]
    for x, y in zip(cycle, cycle[1:]):                                       # fig. 6: T o T -> 0010T -> ... -> T
        assert same(product(x, x), y)
    assert same(product(one(zero(T)), T), zero(zero(one(zero(zero(T))))))    # 10T o T -> 00100T
    v = one(zero(T))
    assert same(product(T, v), zero(zero(one(zero(v)))))                     # eq. 80
    assert same(product(one(v), zero(T)), T)                                 # eq. 82
    assert same(product(zero(v), T), v)                                      # eq. 84


# --- network ---------------------------------------------------------------------
def test_rates_give_the_flow_reactor_equation():
    # eq. 18 with unit rates on ordered collisions and a dilution flux keeping sum x = 1
    terms = [FIG1[k] for k in "CDE"]
    net = generate_network("alchemy", terms=terms)
    ids, f = rhs(net)
    x = np.random.default_rng(0).uniform(0.1, 1.0, len(ids))
    x /= x.sum()
    conc = dict(zip(ids, x))
    production = Counter()
    for a in ids:
        for b in ids:
            c = to_id(normal_form(app(parse(a), parse(b)))[0])
            production[c] += conc[a] * conc[b]
    flux = sum(production.values())
    for s, dx in zip(ids, f(0.0, x)):
        assert dx == pytest.approx(production[s] - conc[s] * flux, abs=1e-12)


def test_default_soup_network():
    traj = evolve("alchemy", seed=3)
    net = traj.network
    assert traj.clock == "collisions" and [f.t for f in traj.frames][:3] == [0.0, 100.0, 200.0]
    assert net.status == "observed" and net.outflow == "constant-total"
    assert sum(net.initial_state.values()) == 100 and len(net.initial_state) == 100
    assert sum(net.extras["final_state"].values()) == 100
    for s in net.species:
        assert parse(s.structure) == parse(s.id) and s.structure.startswith("λ")
    for r in net.reactions:
        assert r.catalysts == r.reactants and r.count >= 1 and r.rate["k"] in (1.0, 2.0)


def test_bad_parameters():
    with pytest.raises(ValueError, match="p_variable"):
        generate_network("alchemy", p_variable=0.7, p_abstraction=0.5)
    with pytest.raises(ValueError, match="smaller"):
        evolve("alchemy", terms=[I, K, S], M=2)
    with pytest.raises(ValueError, match="belongs to the evolve face"):
        generate_network("alchemy", collisions=10)
    with pytest.raises(ValueError, match="is gone"):
        generate_network("alchemy", method="soup")
