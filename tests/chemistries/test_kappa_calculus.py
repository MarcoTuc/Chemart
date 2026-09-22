"""Kappa calculus: Kappa rules flattened into reaction networks over complexes.

Reference guide: Fontana, Boutillier, Feret & Krivine, The Kappa Language v4
(sec. 1.3 ABC model and fig. 3; 2.3.1.1 dimerisation rules; 3.3.1 symmetry
stances; 3.5 footnote on homodimerisation activity). KaDE (Camporesi, Feret &
Ly, CMSB 2017) and Blinov et al. (2006): the EGFR network has 356 species and
3749 reactions. KappaTools examples/abc.ka: the ABC model in edit notation.
"""

import re

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.kappa_calculus import ABC, _split_top, canonical, kappa3_to_kappa4
from chemart.simulate import integrate

ID = "kappa-calculus"

A = "A(x[.],c[.])"
B = "B(x[.])"
C_UU, C_PU, C_PP = "C(x1{u}[.],x2{u}[.])", "C(x1{p}[.],x2{u}[.])", "C(x1{p}[.],x2{p}[.])"
AB = "A(x[1],c[.]),B(x[1])"
ABC_U = "A(x[1],c[2]),B(x[1]),C(x1{u}[2],x2{u}[.])"
ABC_P = "A(x[1],c[2]),B(x[1]),C(x1{p}[.],x2{u}[2])"
AC_X1 = "A(x[.],c[1]),C(x1{u}[1],x2{u}[.])"
AC_X2 = "A(x[.],c[1]),C(x1{p}[.],x2{u}[1])"


def reaction_set(net):
    return {(tuple(sorted(r.reactants.items())), tuple(sorted(r.products.items())), round(r.rate["k"], 12))
            for r in net.reactions}


def R(lhs, rhs, k):
    from collections import Counter
    return (tuple(sorted(Counter(lhs).items())), tuple(sorted(Counter(rhs).items())), k)


def assert_conserves_agents(net):
    ids, Rm, Pm = net.matrices()
    S = (Pm - Rm).toarray()
    assert net.extras["conservation"]
    for law in net.extras["conservation"]:
        m = np.array([law["vector"][s] for s in ids])
        assert np.all(m @ S == 0), law["name"]


# --- the ABC model (reference guide 1.3) -------------------------------------------
def test_abc_flattens_into_its_reactions():
    """The five ABC rules (rule 1 reversible) and their rates give these 11 reactions over 10 complexes."""
    net = generate_network(ID)
    assert net.status == "complete"
    assert {s.id for s in net.species} == {A, B, C_UU, C_PU, C_PP, AB, ABC_U, ABC_P, AC_X1, AC_X2}
    assert reaction_set(net) == {
        R([A, B], [AB], 1e-4), R([AB], [A, B], 0.1),                       # rule 1
        R([AC_X1, B], [ABC_U], 1e-4), R([ABC_U], [AC_X1, B], 0.1),
        R([AC_X2, B], [ABC_P], 1e-4), R([ABC_P], [AC_X2, B], 0.1),
        R([AB, C_UU], [ABC_U], 1e-4),                                        # rule 2
        R([ABC_U], [AB, C_PU], 1.0), R([AC_X1], [A, C_PU], 1.0),            # rule 3
        R([A, C_PU], [AC_X2], 1e-4),                                         # rule 4
        R([AC_X2], [A, C_PP], 1.0),                                          # rule 5
    }
    assert net.initial_state == {A: 1000.0, B: 1000.0, C_UU: 10000.0}
    labels = dict(zip([r.to_text().split("  [")[0] for r in net.reactions], net.extras["reaction_rules"]))
    assert labels[f"{AB} -> {A} + {B}"] == ["rule 1_op"]
    assert_conserves_agents(net)


def test_abc_dynamics_follow_figure_3():
    """Fig. 3: Cuu is replaced first by Cpu, then by Cpp; AB settles slightly below 400 after about 20 s."""
    net = generate_network(ID)
    t_eval = [5.0, 20.0, 50.0, 100.0, 250.0]
    y, _ = integrate(net, 250.0, t_eval=t_eval)
    ab = y[AB] + y[ABC_U] + y[ABC_P]                     # observable |A(x[x.B])|
    assert np.all((ab[1:] > 350) & (ab[1:] < 400))
    assert abs(ab[-1] - ab[1]) < 0.05 * ab[-1]
    assert y[C_UU][0] > y[C_PU][0] > y[C_PP][0]           # early: mostly unmodified
    assert y[C_PU][2] > y[C_PU][0] and y[C_PU][2] > y[C_PU][-1]   # Cpu rises, then falls
    assert y[C_PP][-1] > 0.99 * 10000                     # stationary: doubly modified


def test_edit_notation_of_abc_ka_gives_the_same_network():
    """KappaTools examples/abc.ka writes the same model in edit notation and initialises C()."""
    rules = [
        "'a.b' A(x[./1]),B(x[./1]) @ 1.0E-4",
        "'a..b' A(x[1/.]),B(x[1/.]) @ 0.1",
        "'ab.c' A(x[_],c[./2]),C(x1{u}[./2]) @ 1.0E-4",
        "'mod x1' C(x1{u/p}[1/.]),A(c[1/.]) @ 1",
        "'a.c' A(x[.],c[./1]),C(x1{p}[.],x2{u}[./1]) @ 1.0E-4",
        "'mod x2' A(x[.],c[1/.]),C(x1{p}[.],x2{u/p}[1/.]) @ 1",
    ]
    edit = generate_network(ID, model="custom", signatures=ABC["signatures"], rules=rules,
                            init={"A(),B()": 1000, "C()": 10000})
    arrow = generate_network(ID)
    assert reaction_set(edit) == reaction_set(arrow)
    assert edit.initial_state == arrow.initial_state


def test_rules_are_kept_unflattened():
    net = generate_network(ID)
    assert net.extras["rules"] == ABC["rules"]
    assert net.extras["signatures"] == ABC["signatures"]
    kappa = net.extras["kappa"]
    assert "%agent: C(x1{u p},x2{u p})" in kappa and "%init: 10000 C(x1{u},x2{u})" in kappa
    assert ABC["rules"][1] in kappa.splitlines()
    assert all(s.structure == s.id for s in net.species)


# --- canonical complexes -------------------------------------------------------------
def test_isomorphic_complexes_are_one_species():
    ring = ["A(x,y)"]
    assert canonical("A(x[2],y[1]),A(x[3],y[2]),A(x[1],y[3])", ring) == \
        canonical("A(y[7],x[9]),A(y[9],x[8]),A(x[7],y[8])", ring)
    assert canonical("A(x[1],y[2]),A(x[2],y[1])", ring) != canonical("A(x[1],y[.]),A(x[.],y[1])", ring)
    assert canonical("B(x[5]),A(c,x[5])", ABC["signatures"]) == [AB]

    rng = np.random.default_rng(0)
    net = generate_network(ID, model="egfr", max_complex_size=6, max_species=150)
    for s in net.species:
        agents = _split_top(s.id, ",")
        shuffled = [agents[i] for i in rng.permutation(len(agents))]
        relabelled = re.sub(r"\[(\d+)\]", lambda m: f"[{int(m.group(1)) + 40}]", ",".join(shuffled))
        assert canonical(relabelled, net.extras["signatures"]) == [s.id]


def test_isomorphic_products_of_different_rules_merge():
    """Two rules that build the same heterodimer with its bond written in opposite order give one reaction."""
    net = generate_network(ID, model="custom", signatures=["A(x)", "B(x)"],
                           rules=["'r' A(x[.]), B(x[.]) -> A(x[1]), B(x[1]) @ 1",
                                  "'s' B(x[.]), A(x[.]) -> B(x[2]), A(x[2]) @ 2"],
                           init={"A()": 1, "B()": 1})
    assert len(net.species) == 3 and len(net.reactions) == 1
    assert net.reactions[0].rate["k"] == 3.0 and net.extras["reaction_rules"] == [["r", "s"]]


# --- symmetries and molecularity (reference guide 2.3.1.1, 3.3.1, 3.5) -----------------------
def _k(**params):
    net = generate_network(ID, model="custom", **params)
    return {r.to_text().split("  [")[0]: r.rate["k"] for r in net.reactions}


@pytest.mark.parametrize("convention, k1p, k2p", [("lhs", 1.0, 1.0), ("rule", 2.0, 1.0), ("none", 2.0, 2.0)])
def test_symmetry_stances_of_rules_1p_and_2p(convention, k1p, k2p):
    """3.3.1: under stance ND '1p' and '2p' fire equally; under D '1p' fires twice as often as '2p' and as under ND."""
    common = dict(signatures=["A(x{u p})"], init={"A(x{u}[1]),A(x{u}[1])": 1}, symmetry=convention)
    one = _k(rules=["'1p' A(x{u}[1]), A(x{u}[1]) -> A(x{p}[1]), A(x{u}[1]) @ 1"], **common)
    two = _k(rules=["'2p' A(x{u}[1]), A(x{u}[1]) -> A(x{p}[1]), A(x{p}[1]) @ 1"], **common)
    assert one == {"A(x{u}[1]),A(x{u}[1]) -> A(x{p}[1]),A(x{u}[1])": k1p}
    assert two == {"A(x{u}[1]),A(x{u}[1]) -> A(x{p}[1]),A(x{p}[1])": k2p}


@pytest.mark.parametrize("convention, forward, backward", [("rule", 0.5, 1.0), ("lhs", 0.5, 1.0), ("none", 1.0, 2.0)])
def test_homodimerisation_activity(convention, forward, backward):
    """3.5 footnote: A(x[.]),A(x[.]) <-> A(x[1]),A(x[1]) @ g1, g-1 proceeds with activity g1 nA(nA-1)/2
    once symmetry is accounted for (flux k [A]^2 with k = g1/2); KaSim's uncorrected count gives g1 nA(nA-1)."""
    ks = _k(signatures=["A(x)"], rules=["'d' A(x[.]), A(x[.]) <-> A(x[1]), A(x[1]) @ 1, 1"],
            init={"A()": 10}, symmetry=convention)
    assert ks == {"2 A(x[.]) -> A(x[1]),A(x[1])": forward, "A(x[1]),A(x[1]) -> 2 A(x[.])": backward}


def test_symmetric_dimerization_gives_strict_dimers():
    """2.3.1.1: repeated application of 'symmetric dimerization' results in a population of strict dimers."""
    net = generate_network(ID, model="dimerization", max_complex_size=50)
    assert net.status == "complete"
    assert [s.id for s in net.species] == ["A(x[.])", "A(x[1]),A(x[1])"]
    assert net.initial_state is None


def test_polymerization_is_truncated_by_complex_size():
    """2.3.1.1: asymmetric binding polymerises; with the unimolecular rate {0.1} chain ends close into rings."""
    n = 6
    net = generate_network(ID, model="polymerization", max_complex_size=n)
    assert net.status == "truncated" and net.extras["analysis"]["truncated_by_max_complex_size"]
    sizes = {s.id: s.id.count("A(") for s in net.species}
    rings = {s for s in sizes if "[.]" not in s}
    assert sorted(sizes[s] for s in rings) == list(range(2, n + 1))
    assert sorted(sizes[s] for s in sizes if s not in rings) == list(range(1, n + 1))
    closures = [r for r in net.reactions if sum(r.reactants.values()) == 1]
    joins = [r for r in net.reactions if sum(r.reactants.values()) == 2]
    assert len(closures) == n - 1 and all(r.rate["k"] == 0.1 for r in closures)
    assert len(joins) == sum(1 for i in range(1, n) for j in range(i, n) if i + j <= n)
    for r in joins:     # two orientations for different chains, one ordered pair of copies for equal ones
        assert r.rate["k"] == pytest.approx(0.001 if len(r.reactants) == 1 else 0.002)
    assert_conserves_agents(net)

    small = generate_network(ID, model="polymerization", max_species=5)
    assert small.status == "truncated" and len(small.species) == 5


# --- EGFR (Blinov et al. 2006; KaDE benchmark egfr_net.ka) -----------------------------------
def test_kappa3_translation():
    assert kappa3_to_kappa4("egfr(l!_,r) , egf(r)") == "egfr(l[_],r[.]) , egf(r[.])"
    assert kappa3_to_kappa4("Shc(PTB!1,Y317~pY),x(a?,b!s.T) @'k'{0}", {"k": 2}) == \
        "Shc(PTB[1],Y317{pY}[.]),x(a,b[s.T]) @2.0{0}"


@pytest.mark.slow
def test_egfr_network_has_published_size():
    """Blinov et al. (2006): 356 molecular species connected by 3749 unidirectional reactions (KaDE: 356 ground species)."""
    net = generate_network(ID, model="egfr")
    assert net.status == "complete"
    assert (len(net.species), len(net.reactions)) == (356, 3749)
    assert net.extras["analysis"]["largest_complex"] == 14
    assert_conserves_agents(net)


# --- parameter checks ---------------------------------------------------------------------
def test_invalid_models_are_rejected():
    sig = ["A(x{u p})", "B(x)"]
    with pytest.raises(ValueError, match="only used with model='custom'"):
        generate_network(ID, rules=["A(x[.]) -> A(x[.]) @ 1"])
    with pytest.raises(ValueError, match="creation"):
        generate_network(ID, model="custom", signatures=sig, rules=["A(x[.]), . -> A(x[1]), B(x[1]) @ 1"],
                         init={"A()": 1})
    with pytest.raises(ValueError, match="no site"):
        generate_network(ID, model="custom", signatures=sig, rules=["A(y[.]) -> A(y[.]) @ 1"], init={"A()": 1})
    with pytest.raises(ValueError, match="unimolecular"):
        generate_network(ID, model="custom", signatures=sig, rules=["A(x{u}) -> A(x{p}) @ 1 {2}"], init={"A()": 1})
    with pytest.raises(ValueError, match="rates are numbers"):
        generate_network(ID, model="custom", signatures=sig, rules=["A(x{u}) -> A(x{p}) @ 'k'"], init={"A()": 1})
    with pytest.raises(ValueError, match="init"):
        generate_network(ID, model="custom", signatures=sig, rules=["A(x{u}) -> A(x{p}) @ 1"])
