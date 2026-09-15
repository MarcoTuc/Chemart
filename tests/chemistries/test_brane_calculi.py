"""Brane calculi reproduce the reductions of Cardelli (2004), book 9.7 [157].

Cardelli, "Brane calculi: interactions of biological membranes", CMSB 2004,
LNCS 3082 (author's copy 2004-06-01): structural congruence (sec. 2.1),
Phago/Exo/Pino (sec. 3.1), Mate/Bud/Drip and their encodings (sec. 3.2),
viral infection figs. 7-8 (sec. 3.3), bind&release (4.1), plant vacuole (4.4),
Eat Me and Seek and Store (4.5), nucleocapsid replication (4.6).
"""

import pytest

from chemart import generate_network
from chemart.chemistries import brane_calculi as bc

ID = "brane-calculi"
BITONAL = {"phago": 1, "exo": -1, "pino": 1, "mate": -1, "bud": 1, "drip": 1}


def S(text, defs=None):
    return bc.parse_system(text, defs)


def once(start, defs=None):
    """The unique one-step successor of a configuration and the rule that fired."""
    nxt = bc.steps(start if isinstance(start, bc.System) else S(start, defs))
    assert len(nxt) == 1, [t.text for t, _ in nxt]
    return nxt[0]


def chain(texts, rules, defs=None):
    """Each configuration reduces in exactly one way, by `rule`, to the next one."""
    current = S(texts[0], defs)
    for text, rule in zip(texts[1:], rules):
        nxt, fired = once(current)
        assert (nxt, fired) == (S(text, defs), rule), f"{current.text} -> {nxt.text} by {fired}"
        current = nxt
    return current


def defs_of(system):
    return bc.parse_definitions(bc.PROGRAMS[system])


def check_wellformed(net):
    by_id = {sp.id: S(sp.id) for sp in net.species}
    for sid, term in by_id.items():
        assert " " not in sid and term.text == sid           # canonical and parseable
    for r, rule in zip(net.reactions, net.extras["reaction_rules"]):
        (lhs,), (rhs,) = r.reactants, r.products
        a, b = by_id[lhs], by_id[rhs]
        if rule in BITONAL:
            assert bc.molecule_parities(a) == bc.molecule_parities(b), rule
            assert bc.membranes(b) - bc.membranes(a) == BITONAL[rule], rule


# --- syntax and structural congruence (sec. 2.1) -------------------------------------
def test_structural_congruence_normal_form():
    assert S("a|b[P, Q]") == S("b|a[Q, P]")                          # commutative monoids
    assert S("a|0[P, 0]") == S("a[P]")                               # units
    assert S("[], m") == S("m") and S("0[0]") == S("0")              # 0⟨◇⟩ ≡ ◇
    assert S("!m, m") == S("!m") and S("!!m") == S("!m")             # !P ≡ P∘!P, !!P ≡ !P
    assert S("!(m, n)") == S("!m, !n") and S("!0") == S("0")
    assert S("!(a|b)[P]") == S("!a|!b[P]") and S("a|!a[P]") == S("!a[P]")
    assert S("!(a[P])") != S("!a[P]")                                # replicated membrane vs brane
    assert S("phago⟨P⟩ ∘ ⟨◇⟩ ∘ Q") == S("phago[P], Q")               # the paper's glyphs
    assert S("a.(b|c)|d[P]").text == "a.(b|c)|d[P]"
    assert S("m, a=>b[P]") != S("m, (a=>b)[P]")                       # molecule lists are greedy
    for text in ("x()=>y,z().phago_n'[P]", "(H+)=>(Cl-,H+)[]", "!(s[!m, t[]]), m", "cophago_n(r|s).!t[]",
                 "m, (a(b)=>c(d).s|s0)[P]", "!m, (!a=>b)[], n"):
        assert S(S(text).text) == S(text)


def test_program_errors():
    with pytest.raises(ValueError, match="main"):
        bc.parse_program("a := m")
    with pytest.raises(ValueError, match="recursive"):
        bc.parse_program("a := b[a]\nmain := a")
    with pytest.raises(ValueError, match="must be on a membrane"):
        S("phago, m")
    with pytest.raises(ValueError, match="cannot parse"):
        S("a[P")
    with pytest.raises(ValueError, match="custom"):
        generate_network(ID, system="custom")
    with pytest.raises(ValueError, match="only used"):
        generate_network(ID, program="main := m")


# --- reaction rules (figs. 3, 4, sec. 4.1) --------------------------------------------
def test_phago_exo_pino_rules():
    assert once("phago_n.s|s0[P], cophago_n(r).t|t0[Q]") == (S("t|t0[r[s|s0[P]], Q]"), "phago")
    assert once("coexo_n.t|t0[exo_n.s|s0[P], Q]") == (S("P, s|s0|t|t0[Q]"), "exo")
    assert once("pino(r).s|s0[P]") == (S("s|s0[r[], P]"), "pino")
    assert bc.steps(S("phago_n[P], cophago_m(r)[Q]")) == ()           # names pair actions
    assert bc.steps(S("coexo[exo_n[P]]")) == ()
    assert once("phago[P], cophago(r)[Q]")[1] == "phago"               # omitted names match


def test_mate_bud_drip_rules():
    assert once("mate_n.s|s0[P], comate_n.t|t0[Q]") == (S("s|s0|t|t0[P, Q]"), "mate")
    assert once("cobud_n(r).t|t0[bud_n.s|s0[P], Q]") == (S("r[s|s0[P]], t|t0[Q]"), "bud")
    assert once("drip(r).s|s0[P]") == (S("r[], s|s0[P]"), "drip")


def test_bind_and_release_rule():
    assert once("a, (a(b)=>c(d).s|s0)[b, P]") == (S("c, s|s0[d, P]"), "bind&release")
    assert bc.steps(S("a(b)=>c(d)[P]")) == ()                           # nothing to bind
    assert once("!a, (a=>b)[]") == (S("!a, b"), "bind&release")        # unfolds !a; 0[] ≡ 0 afterwards
    assert once("a, a, a, (a,a=>b)[]") == (S("a, b"), "bind&release")  # binds a multiset
    assert bc.steps(S("a, (a,a=>b)[]")) == ()


# --- sec. 3.2: Mate, Bud, Drip encoded by Phago, Exo, Pino (page 6 derivations) ---------
def test_mate_encoding_derivation():
    end = chain([
        "phago_n.exo_n'.s|s0[P], cophago_n(coexo_n'.exo_n'').coexo_n''.t|t0[Q]",
        "coexo_n''.t|t0[coexo_n'.exo_n''[exo_n'.s|s0[P]], Q]",
        "coexo_n''.t|t0[exo_n''|s|s0[], P, Q]",
        "s|s0|t|t0[P, Q]",
    ], ["phago", "exo", "exo"])
    assert end == once("mate_n.s|s0[P], comate_n.t|t0[Q]")[0]


def test_bud_encoding_derivation():
    end = chain([
        "pino(cophago_n(r).exo_n').coexo_n'.t|t0[phago_n.s|s0[P], Q]",
        "coexo_n'.t|t0[cophago_n(r).exo_n'[], phago_n.s|s0[P], Q]",
        "coexo_n'.t|t0[exo_n'[r[s|s0[P]]], Q]",
        "r[s|s0[P]], t|t0[Q]",
    ], ["pino", "phago", "exo"])
    assert end == once("cobud_n(r).t|t0[bud_n.s|s0[P], Q]")[0]


def test_drip_encoding_derivation():
    end = chain([
        "pino(pino(r).exo_n).coexo_n.s|s0[P]",
        "coexo_n.s|s0[pino(r).exo_n[], P]",
        "coexo_n.s|s0[exo_n[r[]], P]",
        "r[], s|s0[P]",
    ], ["pino", "pino", "exo"])
    assert end == once("drip(r).s|s0[P]")[0]


# --- sec. 3.3: viral infection and reproduction (figs. 7, 8) ------------------------------
def test_viral_infection_figure_7():
    defs = defs_of("viral-infection")
    chain([
        "virus, cell",
        "membrane[mate[exo[nucap]], endosome, Z]",
        "membrane[!comate|!coexo[exo[nucap]], Z]",
        "membrane[!comate|!coexo[], nucap, Z]",
    ], ["phago", "mate", "exo"], defs)
    net = generate_network(ID, seed=1)
    assert net.status == "complete" and len(net.species) == 4
    assert net.extras["reaction_rules"] == ["phago", "mate", "exo"]
    assert net.extras["analysis"]["infection (fig. 7)"]["reached"]
    assert net.extras["analysis"]["terminal"] == [S("membrane[nucap, cytosol]", defs).text]
    check_wellformed(net)


def test_viral_reproduction_figure_8():
    defs = defs_of("viral-reproduction")
    chain([
        "membrane[envelope-vesicle, nucap, Z']",
        "membrane|viral-envelope[nucap, Z']",
        "membrane[Z'], virus",
    ], ["exo", "bud"], defs)
    net = generate_network(ID, system="viral-reproduction")
    assert net.extras["analysis"]["reproduction (fig. 8)"]["reached"]
    assert net.extras["reaction_rules"] == ["exo", "bud", "phago"]      # the new virus can re-enter
    check_wellformed(net)


# --- sec. 4: molecules ------------------------------------------------------------------------
def test_eat_me_and_seek_and_store():
    defs = defs_of("eat-me")
    chain(["A, B", "phago[P], n, B", "phago[P], cophago(rho)[Q]", "[rho[[P]], Q]"],
          ["bind&release", "bind&release", "phago"], defs)
    defs = defs_of("seek-and-store")
    chain(["n, C", "pino(=>(n).mate_store)|seek_n[store[]]",
           "seek_n[=>(n).mate_store[], store[]]", "seek_n[mate_store[n], store[]]", "seek_n[store[n]]"],
          ["bind&release", "pino", "bind&release", "mate"], defs)
    for system in ("eat-me", "seek-and-store"):
        net = generate_network(ID, system=system)
        assert net.status == "complete"
        assert all(v["reached"] for k, v in net.extras["analysis"].items() if isinstance(v, dict))


def test_plant_vacuole_pumps_and_channels():
    defs = defs_of("plant-vacuole")
    vac = "ProtonPump|IonChannel|ProtonAntiporter"
    pumped = S(f"{vac}[H+, H+], ADP, Pi, Cl-, Na+", defs)
    assert once("ATP, Cl-, Na+, PlantVacuole", defs) == (pumped, "bind&release")   # nothing else fires first
    assert {t for t, _ in bc.steps(pumped)} == {
        S(f"{vac}[H+, H+, Cl-], ADP, Pi, Na+", defs),               # ion channel: Cl- in, H+ given back
        S(f"{vac}[H+, Na+], ADP, Pi, Cl-, H+", defs),               # antiporter: Na+ in, H+ out
    }
    net = generate_network(ID, system="plant-vacuole")
    assert net.status == "complete"
    assert net.extras["analysis"]["terminal"] == [S(f"{vac}[H+, Cl-, Na+], ADP, Pi, H+", defs).text]


def test_nucleocapsid_replication_section_4_6():
    defs = defs_of("viral-replication")
    assert (S("cytosol, vRNA, !bud[]", defs), "bind&release") in bc.steps(S("nucap, cytosol", defs))
    assert once("vRNA-repl, vRNA", defs) == (S("vRNA-repl, vRNA, vRNA", defs), "bind&release")       # (a)
    tran = "!vRNA=>vRNA.drip(capsomers)"
    s1 = S(f"drip(capsomers)|{tran}[], vRNA", defs)                                                    # (b)
    assert (s1, "bind&release") in bc.steps(S("capsomer-tran, vRNA", defs))
    s2 = S("capsomers[], capsomer-tran, vRNA", defs)
    assert (s2, "drip") in bc.steps(s1)
    assert (S("capsomer-tran, nucap", defs), "bind&release") in bc.steps(s2)
    er = "!vRNA=>vRNA.drip(exo.viral-envelope)"                                                        # (c)
    s1 = S(f"drip(exo.viral-envelope)|{er}[Nucleus], vRNA", defs)
    assert (s1, "bind&release") in bc.steps(S("ER, vRNA", defs))
    assert (S("ER, vRNA, envelope-vesicle", defs), "drip") in bc.steps(s1)


def test_nucleocapsid_replication_closure_makes_nucaps_and_vesicles():
    """nucap ∘ cytosol →* nucap^n ∘ envelope-vesicle^m ∘ cytosol ∘ vRNA^p ∘ !bud⟨⟩ (sec. 4.6)."""
    defs = defs_of("viral-replication")
    net = generate_network(ID, system="viral-replication")               # 1000 configurations, ~0.3 s
    assert net.status == "truncated"
    assert net.extras["analysis"]["disassembly (sec. 4.6)"]["reached"]
    assert net.extras["analysis"]["max_nucap_copies"] >= 2
    nucap, vesicle, residue = (next(iter(dict(S(t, defs).plain))) for t in ("nucap", "envelope-vesicle", "!bud[]"))
    both = [s for s in map(S, (sp.id for sp in net.species))
            if bc.copies(s, nucap) >= 1 and bc.copies(s, vesicle) >= 1 and bc.copies(s, residue) >= 1]
    assert both


def test_truncation_and_custom_program():
    net = generate_network(ID, system="custom", program="main := !pino(r)[P]", max_species=10)
    assert net.status == "truncated" and len(net.species) == 10
    assert set(net.extras["reaction_rules"]) == {"pino"}
    check_wellformed(net)
    net = generate_network(ID, system="custom", program="""
        A := =>n.phago[P]   # Eat Me, written as a custom program
        B := n=>
             .cophago(rho)[Q]
        main := A, B""")
    assert net.extras["analysis"]["terminal"] == ["[Q,rho[[P]]]"]
