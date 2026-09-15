"""Fraglets: published programs from Tschudin (2003), the book, PyCellChemistry and the 2007 tutorial."""

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import fraglets as fr


def run(fraglet, dialect="pycellchem", segments=None, node="n"):
    """Products of one transformation, as symbol strings."""
    chem = fr.Chemistry(dialect, segments or {})
    syms = tuple(fraglet.split())
    out = chem.transform((node, syms))
    return sorted(" ".join(s) for _, s in out)


def reactions(net):
    return {(tuple(sorted(r.reactants.items())), tuple(sorted(r.products.items()))) for r in net.reactions}


def ssa_final(program, dialect="fraglets-2007", seed=0, steps=10_000):
    net = generate_network("fraglets", program=program, dialect=dialect, method="ssa", steps=steps, seed=seed)
    return net, net.extras["final_state"]


# --- instruction set --------------------------------------------------------------
def test_book_table_16_1_transformations():
    """Book table 16.1 (from AINS 2003), in the pycellchem dialect."""
    assert run("dup A B t") == ["A B B t"]
    assert run("exch A B C t") == ["A C B t"]
    assert run("nop t u") == ["t u"]
    assert run("split s1 s2 * t u") == ["s1 s2", "t u"]
    chem = fr.Chemistry("pycellchem", {"net": ["A", "B"]})
    assert chem.transform(("A", ("send", "B", "t", "u"))) == [("B", ("t", "u"))]
    assert chem.transform(("A", ("send", "C", "t"))) == []          # no link: lost


def test_book_table_16_1_match_rules():
    chem = fr.Chemistry("pycellchem", {})
    m = ("n", ("match", "A", "x", "y"))
    mp = ("n", ("matchp", "A", "x", "y"))
    data = ("n", ("A", "z"))
    assert chem.match(m, data) == [("n", ("x", "y", "z"))]
    assert chem.match(mp, data) == [mp, ("n", ("x", "y", "z"))]
    assert chem.match(m, ("n", ("B", "z"))) is None
    assert chem.match(m, ("other", ("A", "z"))) is None               # different vessel
    assert chem.match(m, ("n", ("dup", "A", "z"))) is None            # instructions are never passive


def test_pycellchem_reference_edge_cases():
    """PyCellChemistry Fraglets.r_dup, r_exch, r_pop, r_fork, r_split on short fraglets."""
    assert run("dup a") == ["a"]
    assert run("exch a b") == ["a b"]
    assert run("pop a b") == ["a"]
    assert run("pop h a t") == ["h t"]
    assert run("fork a") == ["a"]
    assert run("fork a b t") == ["a t", "b t"]
    assert run("nul a b") == []
    assert run("split a b") == ["a b"]                 # no star: tail only
    assert run("split * a") == ["a"]
    assert not fr.keep(("dup",), "pycellchem")         # one-symbol instruction never stored
    assert fr.keep(("dup",), "fraglets-2007")


def test_fraglets_2007_instruction_set_examples():
    """Examples of the 2007 instruction set and of the 2007 tutorial (slides 16-18)."""
    d = "fraglets-2007"
    assert run("dup t a b", d) == ["t a a b"]
    assert run("exch t a b c", d) == ["t b a c"]
    assert run("fork a b x y", d) == ["a x y", "b x y"]
    assert run("pop2 h t a b c", d) == ["h a", "t b c"]
    assert run("length t a b c", d) == ["t 3 a b c"]
    assert run("sum total 3 4 rest", d) == ["total 7 rest"]
    assert run("lt y n 1 2 rest", d) == ["y 1 2 rest"]
    assert run("lt y n 9 7 rest", d) == ["n 9 7 rest"]
    assert run("eq y n x x", d) == ["y x x"]
    assert run("copy this is a fraglet", d) == ["this is a fraglet", "this is a fraglet"]
    assert run("empty finish continue 6 7 8", d) == ["continue 6 7 8"]
    assert run("empty finish continue", d) == ["finish"]
    assert run("sub t 3 5", d) == ["t -2"]
    assert run("abs t -2 x", d) == ["t 2 x"]
    assert run("div t -7 2", d) == ["t -3"] and run("mod t -7 2", d) == ["t -1"]   # C truncation
    assert run("div t 1 0", d) == []
    assert run("pow t 2 10", d) == ["t 1024"]


def test_frag_c_drops_short_fraglets():
    """frag.c: too-short transformations remove the fraglet, unlike PyCellChemistry."""
    d = "fraglets-2007"
    assert run("dup a", d) == [] and run("exch a b", d) == [] and run("fork a", d) == []
    assert run("pop a b", d) == ["a"]
    assert run("lt y n a b", d) == []                  # not numbers
    assert run("send stdout hello", d) == []
    seg = {"ch": ["a", "b"]}
    assert run("send ch b x y", d, seg, node="a") == ["x y"]
    assert run("send ch c x", d, seg, node="a") == []


# --- published programs -------------------------------------------------------------
def test_active_cdp_closure_is_the_ains_trace():
    """Tschudin (2003) sec. IV.B active CDP trace = book fig. 16.10 (default program)."""
    net = generate_network("fraglets")
    assert net.status == "complete"
    code = "a[matchp,cdp,send,b,split,send,a,ack,*]"
    assert reactions(net) == {
        ((("a[cdp,data]", 1), (code, 1)), ((code, 1), ("a[send,b,split,send,a,ack,*,data]", 1))),
        ((("a[send,b,split,send,a,ack,*,data]", 1),), (("b[split,send,a,ack,*,data]", 1),)),
        ((("b[split,send,a,ack,*,data]", 1),), (("b[data]", 1), ("b[send,a,ack]", 1))),
        ((("b[send,a,ack]", 1),), (("a[ack]", 1),)),
    }
    assert [r.catalysts for r in net.reactions if r.catalysts] == [{code: 1}]
    assert net.extras["compartments"]["b"] == ["b[split,send,a,ack,*,data]", "b[send,a,ack]", "b[data]"]
    assert net.initial_state == {code: 1.0, "a[cdp,data]": 1.0}


def test_passive_cdp_with_installed_code():
    """Tschudin (2003) sec. IV.B, fig. 4: static CDP with the code pre-installed at B (':' syntax)."""
    program = """
    a A net
    a B net
    f A[matchP : cdp : send : B : deliver]
    f B[matchP : deliver : split : send : A : ack : *]
    f A[cdp : data]
    """
    # The paper writes matchP; the interpreters spell it matchp.
    program = program.replace("matchP", "matchp")
    net, final = ssa_final(program, dialect="pycellchem")
    assert final == {
        "A[matchp,cdp,send,B,deliver]": 1, "B[matchp,deliver,split,send,A,ack,*]": 1,
        "A[ack]": 1, "B[data]": 1,
    }
    assert net.extras["bimolecular_events"] == 2 and net.extras["inert"]


def test_header_rewriting_and_code_mobility_2007_tutorial():
    """AINS 2003 fig. 3 header rewriting [matchP in out], and tutorial slide 13 code mobility.

    'out' is a reserved instruction in fraglets0.32, so header rewriting runs in the pycellchem dialect.
    """
    _, final = ssa_final("f [matchp in out]\nf [in 1 2]", dialect="pycellchem")
    assert final == {"[matchp,in,out]": 1, "[out,1,2]": 1}
    program = """
    a a ch
    a b ch
    f a[send ch b match temp send ch a tempis]
    f b[temp 30]
    """
    _, final = ssa_final(program)
    assert final == {"a[tempis,30]": 1}


@pytest.mark.parametrize("program, inputs, expected", [
    ("f [matchp incr exch sum 1]", "f [incr x 5]", {"[x,6]": 1}),                                  # slide 22
    ("f [matchp prepend match store store]", "f [store 7 8]\nf [prepend 4 5 6]", {"[store,4,5,6,7,8]": 1}),  # slide 23
    ("f [matchp append split match store match app1 store * app1]",
     "f [store 1 2]\nf [append 3 4 5]", {"[store,1,2,3,4,5]": 1}),                                  # slide 24
    ("f [matchp del exch tmp2]\nf [matchp tmp2 exch tmp1 * ]\nf [matchp tmp1 split nul]",
     "f [del tag x y z]", {"[tag,y,z]": 1}),                                                         # slide 25
    ("f [counter 0]\nf [matchp count empty stop cnt]\nf [matchp stop match counter total]\n"
     "f [matchp cnt pop cnt1]\nf [matchp cnt1 split match counter incr counter * count]\n"
     "f [matchp incr exch sum 1]", "f [count a b c]", {"[total,3]": 1}),                             # slides 26-27
    ("f [matchp getmin length len1]\nf [matchp len1 lt getmin2 min2 1]\nf [matchp min2 pop d1]\n"
     "f [matchp d1 pop min]\nf [matchp getmin2 pop d11]\nf [matchp d11 pop getmin3]\n"
     "f [matchp getmin3 lt islt nlt]\nf [matchp nlt pop getmin]\nf [matchp islt exch nlt]",
     "f [getmin 8 99 4 23]", {"[min,4]": 1}),                                                        # slides 42-43
])
def test_tutorial_programs_compute_their_goal(program, inputs, expected):
    """Fraglets tutorial (Yamamoto 2007): each program's published goal is its output."""
    net, final = ssa_final(program + "\n" + inputs)
    code = {s for s in net.initial_state if s.startswith("[matchp")}
    assert {s: c for s, c in final.items() if s not in code} == expected
    assert net.extras["inert"]


def test_pycellchem_quine_regenerates_itself():
    """PyCellChemistry Fraglets.quine(): 'FNbMbFNb' = [fork nop b match b fork nop b]."""
    net = generate_network("fraglets", program="f [fork nop b match b fork nop b]")
    q = "[fork,nop,b,match,b,fork,nop,b]"
    assert net.status == "complete" and len(net.species) == 4
    assert reactions(net) == {
        (((q, 1),), (("[b,match,b,fork,nop,b]", 1), ("[nop,match,b,fork,nop,b]", 1))),
        ((("[nop,match,b,fork,nop,b]", 1),), (("[match,b,fork,nop,b]", 1),)),
        ((("[b,match,b,fork,nop,b]", 1), ("[match,b,fork,nop,b]", 1)), ((q, 1),)),
    }
    # Transformations are instantaneous, so every settled state holds the two halves, never the quine.
    run_net, final = ssa_final("f [fork nop b match b fork nop b]3", dialect="pycellchem", steps=50)
    assert final == {"[b,match,b,fork,nop,b]": 3, "[match,b,fork,nop,b]": 3}
    forks = [r.count for r in run_net.reactions if r.reactants == {q: 1}]
    assert forks == [3 + 50] and not run_net.extras["inert"]


def test_pycellchem_codegrowth_elongates_without_end():
    """PyCellChemistry Fraglets.codegrowth(): 'SZaDa*aa' = [split matchp a dup a * a a] grows tails."""
    net = generate_network("fraglets", program="f [split matchp a dup a * a a]", max_species=30)
    assert net.status == "truncated"
    lengths = [len(s.id.split(",")) for s in net.species if set(s.id.strip("[]").split(",")) == {"a"}]
    assert max(lengths) >= 10 and sorted(lengths) == list(range(2, max(lengths) + 1))
    _, final = ssa_final("f [split matchp a dup a * a a]", dialect="pycellchem", steps=20)
    assert final == {"[matchp,a,dup,a]": 1, "[" + ",".join(["a"] * 22) + "]": 1}


def test_lossy_link_emulation_2007_tutorial():
    """Tutorial slides 14-15: b[msg]44 of 100 for 50% loss and b[msg]74 for 25% loss."""
    base = "a a ch\na b ch\nf a[transmit b msg]100\nf a[matchp transmit nul]\n"
    for copies, p, published in ((1, 0.5, 44), (3, 0.75, 74)):
        program = base + f"f a[matchp transmit send ch]{copies}\n"
        got = [ssa_final(program, seed=s, steps=200)[1].get("b[msg]", 0) for s in range(40)]
        assert abs(np.mean(got) - 100 * p) < 3
        sd = (100 * p * (1 - p)) ** 0.5
        assert abs(published - 100 * p) < 2.5 * sd       # the published single run is a typical draw


def test_ssa_observed_network_balances():
    program = fr.CDP.replace("f a[cdp data]", "f a[cdp data]5")
    net = generate_network("fraglets", program=program, method="ssa", seed=3)
    assert net.status == "observed"
    total = dict(net.initial_state)
    for r in net.reactions:
        for s, n in r.reactants.items():
            total[s] = total.get(s, 0) - n * r.count
        for s, n in r.products.items():
            total[s] = total.get(s, 0) + n * r.count
    assert {s: v for s, v in total.items() if v} == net.extras["final_state"]
    assert net.extras["final_state"]["a[ack]"] == 5 and net.extras["final_state"]["b[data]"] == 5


def test_bad_programs():
    with pytest.raises(ValueError, match="no fraglet"):
        generate_network("fraglets", program="a a net")
    with pytest.raises(ValueError, match="line 1"):
        generate_network("fraglets", program="x [a b]")
    with pytest.raises(ValueError, match="not supported"):
        generate_network("fraglets", program="f [wait a b]", dialect="fraglets-2007")
    with pytest.raises(ValueError, match="end marker"):
        generate_network("fraglets", program="f [a b]\ne\nf [c]")
