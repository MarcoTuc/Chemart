"""Typogenetics reproduces Hofstadter's GEB examples, book fig. 10.6 and Snare (1999).

GEB pp. 508-511 examples as reproduced by kortschak/typo and mkpankov/typogenetics;
book 10.5.2 table 10.1 and fig. 10.6 (from Morris [597]); Snare's thesis, chapters 2-3.
"""

from collections import Counter

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.typogenetics import (
    CODE_TABLES, Complex, apply_enzyme, binding_preference, enzyme, enzyme_text, outcomes, translate,
)


def names(enzymes, fold="hofstadter"):
    return [enzyme_text(e, fold) for e in enzymes]


# --- code table and folding ------------------------------------------------------
def test_code_table_hofstadter_table_10_1():
    code = CODE_TABLES["hofstadter"]
    row = lambda first: [code[first + b] for b in "ACGT"]
    assert row("A") == [None, ("cut", "s"), ("del", "s"), ("swi", "r")]
    assert row("C") == [("mvr", "s"), ("mvl", "s"), ("cop", "r"), ("off", "l")]
    assert row("G") == [("ina", "s"), ("inc", "r"), ("ing", "r"), ("int", "l")]
    assert row("T") == [("rpy", "r"), ("rpu", "l"), ("lpy", "l"), ("lpu", "l")]


def test_code_table_varetto_differs_only_in_row_g():
    # Snare table 2.1(b), after Varetto (1993) p. 187
    h, v = CODE_TABLES["hofstadter"], CODE_TABLES["varetto"]
    assert {d for d in h if h[d] != v[d]} == {"GA", "GC", "GG", "GT"}
    assert [v["G" + b] for b in "ACGT"] == [("inc", "s"), ("ing", "r"), ("int", "r"), ("ina", "l")]


def test_translation_examples_snare_2_3():
    assert names(translate("CGCTAATAAGT")) == ["cop-off:A", "rpy-del:A"]   # odd base ignored, AA off-frame
    assert names(translate("CGAAAAGCAAT")) == ["cop:A", "inc:A"]           # null enzymes dropped


def test_geb_p511_binding_preference():
    # rpy-ina-rpu-mvr-int-mvl-cut-swi-cop binds C (GEB p. 511); Morris/Varetto folding gives G (Snare 2.4)
    (e,) = translate("TAGATCCAGTCCACATCGA")
    assert names([e]) == ["rpy-ina-rpu-mvr-int-mvl-cut-swi-cop:C"]
    assert binding_preference(e, "morris") == "G"


# --- enzymes at work -------------------------------------------------------------
def test_geb_worked_example_p508():
    # rpu-inc-cop-mvr-mvl-swi-lpu-int bound to the middle G of TAGATCCAGTCCATCGA
    strand = "TAGATCCAGTCCATCGA"
    e = enzyme("rpu-inc-cop-mvr-mvl-swi-lpu-int")
    assert Complex(strand).sites("G") == [2, 8, 15]
    assert sorted(apply_enzyme(e, strand, 8)) == ["ATG", "TAGATCCAGTCCACATCGA"]
    # and the product codes the enzyme of GEB p. 511
    assert names(translate("TAGATCCAGTCCACATCGA")) == ["rpy-ina-rpu-mvr-int-mvl-cut-swi-cop:C"]


def test_book_fig_10_6():
    # two enzymes of CGACCCAACGATTTTTCAT; first attached at g (unit 9), second at a (unit 2)
    strand = "CGACCCAACGATTTTTCAT"
    first, second = translate(strand)
    assert names([first, second], "morris") == ["cop-cut-mvl:G", "cop-swi-lpu-lpu-mvr:A"]
    c = Complex(strand)
    c.apply(first, 9)
    assert c.pos == 8 and "".join(b or " " for b in c.rows[1]).strip() == "GC"
    assert c.sites("A") == [2, 6, 7]        # only the part the enzyme stayed on
    c.apply(second, 2)
    expected = ("ATTTTTCAT", "CGACCCAACG", "CGTTGGGT")
    assert tuple(sorted(c.strands())) == expected
    assert expected in outcomes([first, second], strand, "morris", "all")


@pytest.mark.parametrize("strand", ["CGATTCGAATCG", "CGATTAATTAATCG"])
def test_snare_self_replicators(strand):
    # Snare 3.3, fig. 3.2: the enzyme binds the rightmost G and copies the inverted repeat
    (e,) = translate(strand)
    assert binding_preference(e) == "G"
    assert outcomes([e], strand, "hofstadter", "rightmost") == [(strand, strand)]
    net = generate_network("typogenetics", strands=[strand])
    assert Counter({strand: 2}) in [Counter(r.products) for r in net.reactions if r.reactants == {strand: 1}]


def test_dud_strand_has_no_reactions():
    # Snare 3.1.1: CGGC codes cop-inc, which binds A
    net = generate_network("typogenetics", strands=["CGGC"])
    assert net.reactions == [] and net.status == "complete"
    assert net.species[0].structure == "cop-inc:A"


# --- networks ---------------------------------------------------------------------
def test_default_closure_lists_alternative_outcomes():
    net = generate_network("typogenetics")
    replicator = "CGATTCGAATCG"
    from_seed = [r for r in net.reactions if r.reactants == {replicator: 1}]
    assert len(from_seed) == 3                           # three G sites, three outcomes
    ids = {s.id for s in net.species}
    for r in net.reactions:
        assert set(r.reactants) | set(r.products) <= ids
        assert len(r.reactants) == 1 and sum(r.reactants.values()) == 1
        assert r.rate is None


def test_pair_reaction_keeps_the_gene():
    net = generate_network("typogenetics", reaction="pair", strands=["CGATTCGAATCG", "CGGC"],
                           binding_tiebreak="rightmost", max_species=30)
    assert net.reactions
    for r in net.reactions:
        reactants = list(Counter(r.reactants).elements())
        assert len(reactants) == 2 and r.catalysts, "the gene strand survives"
    # cop-swi-rpu-ina-swi-cop binds the rightmost G of CGGC (unit 2), copies units 2..0 on the
    # upper strand and runs off its end: the complement of CGG read backwards, CCG (hand-derived)
    r = [x.to_text() for x in net.reactions if x.reactants == {"CGATTCGAATCG": 1, "CGGC": 1}]
    assert "CGATTCGAATCG + CGGC -> CGATTCGAATCG + CCG + CGGC" in r


def test_max_length_truncates():
    net = generate_network("typogenetics", strands=["GTGTGTGT"], max_length=8)   # int-int-int-int inserts
    assert net.status == "truncated"
    assert all(len(s.id) <= 8 for s in net.species)


def test_soup_replicator_outgrows_dud():
    traj = evolve("typogenetics", seed=3, strands=["CGATTCGAATCG", "CGGC"],
                  binding_tiebreak="rightmost", copies=50, steps=3000)
    net = traj.network
    assert net.status == "observed" and net.outflow == "constant-total"
    assert net.initial_state == {"CGATTCGAATCG": 50, "CGGC": 50}
    assert all(r.count >= 1 for r in net.reactions)
    final = net.extras["final_state"]
    assert sum(final.values()) == 100
    assert final.get("CGATTCGAATCG", 0) > 90
    assert traj.times()[:2] == [0.0, 100.0] and traj.frames[-1].state == {s: float(n) for s, n in final.items()}


def test_evolve_reads_tiebreak_all_as_random():
    # a population draws one outcome per collision, so the closure's 'all' becomes 'random'
    default = evolve("typogenetics", seed=1).network.to_dict()
    random = evolve("typogenetics", seed=1, binding_tiebreak="random").network.to_dict()
    assert default["reactions"] == random["reactions"] and default["species"] == random["species"]


def test_bad_parameters():
    with pytest.raises(ValueError, match="A, C, G, T"):
        generate_network("typogenetics", strands=["ACGU"])
    with pytest.raises(ValueError, match="evolve"):
        generate_network("typogenetics", steps=10)
    with pytest.raises(ValueError, match="max_length"):
        generate_network("typogenetics", strands=["ACGT" * 10])
