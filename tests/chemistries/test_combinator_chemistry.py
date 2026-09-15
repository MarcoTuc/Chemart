"""The combinator chemistry reproduces Speroni di Fenizio (2000; ECAL 2001; thesis 2007).

Sources: "A less abstract artificial chemistry" (Artificial Life VII, 2000), "Stability of
metabolic and balanced organisations" (ECAL 2001, with Banzhaf), and the PhD thesis
"Chemical Organization Theory" (Jena, 2007), chapters 6 and 7.
"""

from collections import Counter

import pytest

from chemart import generate_network
from chemart.chemistries.combinator_chemistry import (
    Chemistry, Elastic, apply, atom_counts, normal_form, parse, reduce_terms, show,
)


# --- the atoms (thesis table 6.1; ECAL 2001 table 1) ------------------------------
@pytest.mark.parametrize("term, result", [
    ("Bxyz", ["x(yz)"]),
    ("Cxyz", ["xzy"]),
    ("Ix", ["x"]),
    ("Kxy", ["x"]),
    ("Sxyz", ["xz(yz)"]),
    ("Wxy", ["xyy"]),
    ("Rxy", ["x", "y"]),            # R releases y as a separate molecule
    ("Bxyzu", ["x(yz)u"]),          # s0, the rest of the string, is kept
    ("Bxy", ["Bxy"]),               # too few items: normal form
])
def test_atom_rules(term, result):
    assert normal_form(term) == result


def test_k_releases_in_the_2000_paper():
    assert normal_form("Kxy", k_release=True) == ["x", "y"]


def test_thesis_reduction_examples():
    # thesis 6.2.1
    assert normal_form("BBBBBB") == ["BB(BB)"]
    assert normal_form("B(BBBB)") == ["B(B(BB))"]
    assert normal_form("BBB(BBBB)") == ["B(B(B(BB)))"]
    for nf in ["B", "BB", "BBB", "B(BB)", "BB(BB)", "B(B(BB))"]:
        assert normal_form(nf) == [nf]
    # CKB recreates one of S, K, I (the reader's riddle): it is I
    assert normal_form("CKBx") == ["x"]


def test_parenthesis_simplification():
    # ECAL 2001 and thesis 6.2.1
    assert show(parse("(S(I)((K)K))")) == "SI(KK)"
    assert show(parse("(((B)B)(((B)B)B))")) == "BB(BBB)"
    assert show(parse("B((BB))")) == "B(BB)"
    assert show(parse("B()B")) == "BB"
    assert show(parse("(SI)KKS")) == "SIKKS"                           # 2000 paper
    assert show(apply(parse("SSK"), parse("KI"))) == "SSK(KI)"         # 2000 paper


def test_2000_paper_reduction_examples():
    # KK(SKI) reduces to K; K(SKKI) -> K(KI(KI)) -> KI; K(SKI) is in normal form
    assert normal_form("KK(SKI)") == ["K"]
    assert normal_form("K(SKKI)") == ["KI"]
    assert normal_form("K(SKI)") == ["K(SKI)"]


def test_reduction_order_of_table_6_3():
    # destructors first: K discards the non-terminating WWW before it is touched
    assert normal_form("KI(WWW)") == ["I"]
    # thesis 6.5: Kc(Rab) normalises to c; K acts before R, so b is never released
    assert normal_form("Kc(Rab)") == ["c"]


# --- limits ----------------------------------------------------------------------
def test_reduction_budget():
    assert normal_form("WWW") is None          # WWW reduces to itself (2000 paper)
    assert normal_form("SII(SII)") is None
    with pytest.raises(Elastic, match="reductions"):
        reduce_terms([parse("WWW")], max_reductions=50)
    done, steps, _ = reduce_terms([parse("BBBBBB")], max_reductions=2)
    assert steps == 2 and [show(t) for t in done] == ["BB(BB)"]
    with pytest.raises(Elastic, match="reductions"):
        reduce_terms([parse("BBBBBB")], max_reductions=1)


def test_size_limits_make_reactions_elastic():
    # WW(WW) keeps growing: it is cut off by the size during the reduction
    with pytest.raises(Elastic, match="longer"):
        reduce_terms([parse("WW(WW)")], max_reductions=10_000, max_react_size=30)
    assert Chemistry(max_size=100, max_depth=20).react("WR", "SKI") is not None
    assert Chemistry(max_size=2).react("WR", "SKI") is None             # SKI has 3 atoms
    assert Chemistry(max_depth=0).react("K", "BKK") is None             # K(BKK) has depth 1
    assert Chemistry(max_react_size=4).react("WR", "SKI") is None       # WR(SKI) has 5 atoms


# --- published reactions ---------------------------------------------------------
def test_ecal_2001_example_wr_on_ski():
    # eq. 2: WR * SKI -> WR(SKI) -> R(SKI)(SKI) -> SKI, SKI
    products, _ = Chemistry().react("WR", "SKI")
    assert products == ("SKI", "SKI")
    net = generate_network("combinator-chemistry", molecules=["WR", "SKI"])
    by_lhs = {frozenset(r.reactants.items()): r for r in net.reactions}
    r = by_lhs[frozenset({"WR": 1, "SKI": 1, "free:S": 1, "free:K": 1, "free:I": 1}.items())]
    assert r.products == {"SKI": 2, "free:W": 1, "free:R": 1}


def test_2000_paper_pool_example():
    # pool of 20 atoms per type; SKI applied to KW goes through three published pool states
    _, _, snapshots = reduce_terms([apply(parse("SKI"), parse("KW"))])
    pools = [{x: 20 - s[x] for x in "BCIKSW"} for s in snapshots]
    assert pools == [
        dict(B=20, C=20, I=19, K=18, S=19, W=19),     # SKI(KW) taken from the pool
        dict(B=20, C=20, I=19, K=17, S=20, W=18),     # K(KW)(IKW): S back, K and W copied
        dict(B=20, C=20, I=20, K=19, S=20, W=19),     # KW: 2 K, 1 I, 1 W back
    ]


ALPHA = "C(C(K(CKK))(WC))(C(K(CKK))(WC))"
GAMMA = "S(K(SSK))(K(K(SSK)))"
DELTA = "B(WW)(W(B(WK)W))"


def _products(a, b, k_release):
    done, _, _ = reduce_terms([apply(parse(a, variables=True), parse(b, variables=True))], k_release, 10_000, 1_000)
    return Counter(show(t) for t in done)


def test_2000_paper_organisation_reactions():
    # K releases its second argument; 'a' is an arbitrary element
    assert _products(ALPHA, "a", True) == Counter({ALPHA: 1, "K": 1, "a": 1})
    assert _products("K", "a", True) == Counter({"Ka": 1})
    assert _products("K(K(Kb))", "a", True) == Counter({"K(Kb)": 1, "a": 1})
    assert _products("Kb", "a", True) == Counter({"b": 1, "a": 1})
    assert _products(GAMMA, "a", True) == Counter({GAMMA: 1, "a": 2})
    assert _products(DELTA, "a", True) == Counter({"aaa(aaa)(aaa)": 1, "Wa": 3})


# Thesis ch. 7 (spatial run, reactive, atoms BCIKRSW, MaxLength 100, MaxDepth 20).
MOL = dict(g="B(WR)(WR)", b="SIR", v1="R", v2="RR", y1="W(SS)", w="S(R(W(SS)))(W(SS))",
           r="SII", y2="W(SI)", w2="S(B(WR)(WR))(B(WR)(WR))", w3="S(SIR)(SIR)")
NAME = {v: k for k, v in MOL.items()}


def _named(molecule):
    """The thesis lists R(R(...(x)...)) as x."""
    term = parse(molecule)
    while term[0] == "R" and len(term) == 2 and isinstance(term[1], tuple):
        term = term[1]
    text = show(term)
    return NAME.get(text, text)


CH7 = Chemistry(max_reductions=10_000, max_react_size=100, max_size=100, max_depth=20)


@pytest.mark.parametrize("a, b, published", [
    ("b", "g", "4 g"), ("b", "b", "2 b"), ("b", "v1", "v2"), ("b", "v2", "v1 v2"),
    ("b", "w", None), ("b", "r", "2 r"), ("b", "y2", "3 y2"),
    ("y1", "g", "w2 3 g"), ("y1", "b", "w3 b"), ("w", "g", None), ("w", "b", None),
    ("r", "g", "4 g"), ("r", "b", "2 b"), ("y2", "v1", "v2"), ("y2", "v2", "2 v1 v2"),
    ("b", "w3", None), ("w3", "w3", None), ("g", "w2", "4 w2"), ("w2", "w2", None),
])
def test_thesis_chapter_7_reaction_tables(a, b, published):
    out = CH7.react(MOL[a], MOL[b])
    if published is None:
        assert out is None
        return
    expected, n = Counter(), 1
    for token in published.split():
        if token.isdigit():
            n = int(token)
        else:
            expected[token] += n
            n = 1
    assert Counter(_named(m) for m in out[0]) == expected


@pytest.mark.parametrize("y", ["SIR", "W(SS)", "SII", "KB", "S(KS)(K(SI))"])
def test_thesis_chapter_7_general_rules(y):
    # g + y -> 4y ; R(x) + y -> x, y ; RR(x) + y -> R(y), x
    g, _ = CH7.react("B(WR)(WR)", y)
    assert Counter(g) == Counter({y: 4})
    x = "BC"
    assert Counter(CH7.react(f"R({x})", y)[0]) == Counter({x: 1, y: 1})
    assert Counter(CH7.react(f"RR({x})", y)[0]) == Counter({f"R({y})": 1, x: 1})


# --- organisations as closures --------------------------------------------------
LEVEL1 = dict(reaction="catalytic", atoms="BCIKSW", filter_reproduction=True,
              max_size=15, max_react_size=15, max_depth=7, max_reductions=100)   # thesis fig. 6.4


def _ladder(base, max_depth, base_depth):
    out, m = [], base
    for _ in range(max_depth - base_depth + 1):
        out.append(m)
        m = f"K({m})"
    return out


@pytest.mark.parametrize("seeds, bases", [
    (["BKK"], [("BKK", 0)]),                                                    # generation 605
    (["BKK", "BK(BKK)"], [("BKK", 0), ("BK(BKK)", 1)]),                         # generation 405
    (["BKK", "BK(BKK)", "BK(K(K(BKK)))"], [("BKK", 0), ("BK(BKK)", 1), ("BK(K(K(BKK)))", 3)]),  # 137
])
def test_thesis_fig_6_5_level_1_organisations(seeds, bases):
    # infinite ladders K^n(x) (Fontana & Buss' L1), cut at max_depth 7 by the limits
    net = generate_network("combinator-chemistry", molecules=seeds, **LEVEL1)
    assert net.status == "complete"
    expected = {m for base, d in bases for m in _ladder(base, 7, d)}
    assert {s.id for s in net.species} == expected
    for r in net.reactions:   # catalytic: both reactants survive
        assert all(r.products.get(s, 0) >= n for s, n in r.reactants.items())


def test_2000_paper_type_a_organisation():
    # O = {alpha, K, KK, K^n(alpha), K^n(KK)}, cut at max_depth 6 (alpha has depth 3)
    net = generate_network("combinator-chemistry", molecules=[ALPHA], k_action="release",
                           atoms="BCIKSW", max_size=100, max_react_size=100, max_depth=6)
    assert net.status == "complete"
    expected = set(_ladder(ALPHA, 6, 3)) | {"K"} | set(_ladder("KK", 6, 0))
    assert {s.id for s in net.species if not s.id.startswith("free:")} == expected


def _conserves_atoms(net):
    for law in net.extras["conservation"]:
        m = law["vector"]
        for r in net.reactions:
            assert sum(m[s] * n for s, n in r.reactants.items()) == sum(m[s] * n for s, n in r.products.items())


def test_default_closure_conserves_atoms():
    net = generate_network("combinator-chemistry")
    assert net.status == "truncated"
    assert {law["name"] for law in net.extras["conservation"]} == {f"atom {x}" for x in "BCIKRSW"}
    _conserves_atoms(net)


# --- soups ----------------------------------------------------------------------
def test_reactive_soup_conserves_atoms_and_varies_population():
    # 2000 paper / ECAL 2001: the number of atoms is constant, the number of molecules is not
    net = generate_network("combinator-chemistry", seed=3, method="soup", M=60, generations=15,
                           atoms_per_type=150, prob_destroy=0.01)
    assert net.status == "observed" and all(r.count for r in net.reactions)
    _conserves_atoms(net)
    total = Counter()
    for m, c in net.extras["final_state"].items():
        for x, k in atom_counts(m).items():
            total[x] += k * c
    for x, free in net.extras["final_free_atoms"].items():
        assert total[x] + free == 150
    assert len(set(net.extras["analysis"]["population"])) > 1
    # random molecules are assembled from the pool, and molecules decay back into it
    assert any(all(s.startswith("free:") for s in r.reactants) for r in net.reactions)
    assert any(all(s.startswith("free:") for s in r.products) for r in net.reactions)


def test_catalytic_soup_keeps_its_size():
    net = generate_network("combinator-chemistry", seed=0, method="soup", reaction="catalytic",
                           atoms="BCIKSW", M=50, generations=5)
    assert sum(net.extras["final_state"].values()) == 50
    assert net.outflow == "constant-total"
    for r in net.reactions:
        assert all(r.products.get(s, 0) >= n for s, n in r.reactants.items())


def test_invalid_molecules():
    with pytest.raises(ValueError, match="normal form"):
        generate_network("combinator-chemistry", molecules=["BBBB"])
    with pytest.raises(ValueError, match="outside"):
        generate_network("combinator-chemistry", molecules=["RK"], atoms="BCIKSW")
    with pytest.raises(ValueError, match="atoms"):
        generate_network("combinator-chemistry", atoms="BX")
