"""Tominaga's stacked-string chemistry reproduces the published systems.

2007: Tominaga, Watanabe, Kobayashi, Nakamura, Kishi & Kazuno, Artificial Life
13(3):223-247 (book [856]): notation and matching examples (sec. 2.2), the
Adleman-Lipton model (sec. 6.1, fig. 6), the Benenson DNA automaton (sec. 6.2,
rules 3-11 and the published molecules).
2009: Tominaga, Suzuki, Kobayashi, Watanabe, Koizumi & Kishi, Artificial Life
15(1):115-129 (book [855]): notation (sec. 2, eqs. 1-3), transcription
(sec. 3.2, rules 12-18), fatty acid oxidation (sec. 4, rules 25-28, figs. 1, 3).
"""

from collections import Counter

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.tominaga_stacked_strings import (
    AB_RULES, ACETYL_COA, ADLEMAN_ANSWER, ADLEMAN_RULES, BENENSON_RULES, DETECTOR_S0, DETECTOR_S1,
    FATTY_ACID_RULES, FOK_I, TRANSCRIPTION_RULES, TRANSITIONS, Rule, acyl_coa, benenson_input,
    benenson_reporter, chromosome, layout, match, parse_molecule, parse_pattern,
)
from chemart.chemistries.tominaga_stacked_strings import molecule_text

ID = "tominaga-stacked-strings"


def rules(table):
    return {label: Rule(text, label) for label, text in table}


def once(rule, *mols):
    """The unique outcome of a rule on the given molecules, as texts."""
    out = rule.apply(*mols)
    assert len(out) == 1, f"rule {rule.label} gives {len(out)} outcomes"
    return [molecule_text(m) for m in out[0]]


# --- the notation (2007 sec. 2.2.3, 2009 sec. 2) --------------------------------------
def test_pattern_matching_examples():
    p = parse_pattern("0#1BC2/")                                   # 2009 sec. 2
    assert match(p, parse_molecule("0#ABCD/")) and match(p, parse_molecule("0#BBCF/"))
    assert not match(p, parse_molecule("0#ABC/"))
    assert match(parse_pattern("0#*1AB/1#CD2*/"), parse_molecule("0#XXAB/3#CD/")) == [
        {"*1": ("X", "X"), "2*": ()}]                               # *1 = XX, 2* = null
    q = parse_pattern("0#*1AB/1#CD/")                              # 2007 sec. 2.2.3
    assert match(q, parse_molecule("0#ABAB/3#CD/")) and not match(q, parse_molecule("0#ABAB/1#CD/"))
    assert not match(q, parse_molecule("0#ABAB/"))                 # line counts must agree
    assert match(parse_pattern("0#A1/1#*23E/0#F4*/"), parse_molecule("0#AB/0#CDE/0#FG/"))
    assert parse_molecule("0#OrcTATATT/0#CapRp/")[1] == (0, ("Cap", "Rp"))   # multi-letter v-elements
    assert layout(parse_molecule("0#ABAB/3#CD/")) == "ABAB\n   CD"


def test_ab_concatenation_with_catalyst_source_and_drain():
    r = rules(AB_RULES)
    assert once(r["(1)"], "0#ABAB/", "0#CD/") == ["0#ABAB/3#CD/"]
    assert once(r["(2)"], "0#ABAB/3#CD/", "0#AB/") == ["0#ABABAB/3#CD/"]   # 2007 sec. 2.2.5
    assert once(r["(3)"], "0#ABABAB/3#CD/") == ["0#ABABAB/", "0#CD/"]
    drain = parse_pattern("0#ABABAB1*/")                            # 2009 sec. 2
    assert match(drain, parse_molecule("0#ABABAB/")) and match(drain, parse_molecule("0#ABABABAB/"))
    assert not match(drain, parse_molecule("0#ABAB/"))

    net = generate_network(ID, system="ab-concatenation", max_species=60)
    assert net.status == "truncated"
    labels = net.extras["reaction_rules"]
    ids = {s.id for s in net.species}
    assert {"0#AB/", "0#CD/", "0#ABAB/", "0#ABABAB/"} <= ids
    sources = [r for r, lab in zip(net.reactions, labels) if lab == "source"]
    assert [(r.reactants, r.products) for r in sources] == [({}, {"0#AB/": 1})]
    drained = [next(iter(r.reactants)) for r, lab in zip(net.reactions, labels) if lab == "drain"]
    assert drained and all(len(parse_molecule(s)) == 1 and len(parse_molecule(s)[0][1]) >= 10
                           for s in drained)


def test_custom_rules_equal_the_builtin_system():
    custom = generate_network(ID, system="custom", rules=[t for _, t in AB_RULES], pool={"0#CD/": 1},
                              sources=["0#AB/"], drains=["0#ABABABABAB1*/"], max_species=40)
    builtin = generate_network(ID, system="ab-concatenation", max_species=40)
    assert [s.id for s in custom.species] == [s.id for s in builtin.species]
    assert custom.reactions == builtin.reactions


# --- Adleman-Lipton (2007 sec. 6.1) -----------------------------------------------------
def test_adleman_hybridisation_and_answer_molecule():
    r = rules(ADLEMAN_RULES)
    # fig. 6: U22U31 and L31L32 hybridise at 3' ends, then U32U41 at 5' ends
    first = once(r["(1) U_31"], "0#U_22U_31/1#/", "0#/-1#L_31L_32/")
    assert first == ["0#U_22U_31/1#L_31L_32/"]
    assert once(r["(2) L_32"], first[0], "0#U_32U_41/1#/") == ["0#U_22U_31U_32U_41/1#L_31L_32/"]

    # the answer: path 0 -> 1 -> 2 -> 3 -> 4 -> 5 -> 6, a fully double-stranded molecule
    [mol] = once(r["start"], "0#U_01U_02U_11/2#/", "0#/0#L_01L_02/")
    edges = {1: "0#U_12U_21/1#/", 2: "0#U_22U_31/1#/", 3: "0#U_32U_41/1#/", 4: "0#U_42U_51/1#/",
             5: "0#U_52U_61U_62/1#/"}
    for i in range(1, 6):
        [mol] = once(r[f"(1) U_{i}1"], mol, f"0#/-1#L_{i}1L_{i}2/")
        [mol] = once(r[f"(2) L_{i}2"], mol, edges[i])
    [mol] = once(r["end"], mol, "0#/-2#L_61L_62/")
    assert mol == ADLEMAN_ANSWER
    assert [(d, len(line)) for d, line in parse_molecule(mol)] == [(0, 14), (0, 14)]


def test_adleman_closure_is_truncated():
    net = generate_network(ID, system="adleman-hamiltonian-path", max_species=300)
    assert net.status == "truncated" and len(net.species) <= 300
    assert "0#U_22U_31/1#L_31L_32/" in {s.id for s in net.species}


# --- the Benenson DNA automaton (2007 sec. 6.2) -------------------------------------------
def test_published_input_molecules():
    assert benenson_input("abb", 1) == ("0#GGATGXXXXXXXCTGGCTCGCAGCCGCAGCTGTCGCX/"
                                        "0#CCTACXXXXXXXGACCGAGCGTCGGCGTCGACAGCGX/")
    assert benenson_input("aba", 2) == ("0#GGATGXXXXXXXCTGGCTCGCAGCCTGGCTTGTCGCXX/"
                                        "0#CCTACXXXXXXXGACCGAGCGTCGGACCGAACAGCGXX/")


def fok(r, mol):
    """Rules (3), (4), (5): Fok I binds, cuts 9/13 bases after GGATG and leaves; returns the remnant."""
    [bound] = once(r["(3)"], FOK_I, mol)
    enzyme, remnant = once(r["(4)"], bound)
    assert once(r["(5)"], enzyme)[0] == FOK_I
    return remnant


@pytest.mark.parametrize("word, trail, steps, state", [
    ("abb", 1, [("(6)", TRANSITIONS["S0", "a"], "GGCT"), ("(7)", TRANSITIONS["S0", "b"], "CAGC"),
                ("(9)", TRANSITIONS["S1", "b"], "CGCA"), ("(10)", DETECTOR_S0, "TCGC")], "S0"),
    ("aba", 2, [("(6)", TRANSITIONS["S0", "a"], "GGCT"), ("(7)", TRANSITIONS["S0", "b"], "CAGC"),
                ("(8)", TRANSITIONS["S1", "a"], "CTGG"), ("(11)", DETECTOR_S1, "TGTC")], "S1"),
])
def test_automaton_reaction_sequence(word, trail, steps, state):
    r = rules(BENENSON_RULES)
    partners = list(TRANSITIONS.values()) + [DETECTOR_S0, DETECTOR_S1]
    remnant = fok(r, benenson_input(word, trail))
    assert remnant == ("0#GGCTCGCAGCCGCAGCTGTCGCX/4#GCGTCGGCGTCGACAGCGX/" if word == "abb"
                       else "0#GGCTCGCAGCCTGGCTTGTCGCXX/4#GCGTCGGACCGAACAGCGXX/")
    for label, partner, sticky in steps:
        assert parse_molecule(remnant)[0][1][:4] == tuple(sticky)
        # the sticky end selects exactly one transition molecule (or detector)
        chosen = [(lab, p) for lab in ("(6)", "(7)", "(8)", "(9)", "(10)", "(11)") for p in partners
                  if r[lab].apply(p, remnant)]
        assert chosen == [(label, partner)]
        [mol] = once(r[label], partner, remnant)
        if label in ("(10)", "(11)"):
            break
        remnant = fok(r, mol)
    assert mol == benenson_reporter(state, trail)
    assert mol == ("0#XTCGCX/0#XAGCGX/" if word == "abb" else "0#XTGTCGCXX/0#XACAGCGXX/")


def test_default_network_accepts_abb_only():
    net = generate_network(ID)
    assert net.status == "complete"
    assert net.extras["analysis"] == {"reporters": {"abb": ["S0"], "aba": []}, "accepted": ["abb"]}
    assert net.initial_state[FOK_I] == 100 and net.initial_state[DETECTOR_S0] == 20
    assert all(net.initial_state[t] == 20 for t in TRANSITIONS.values())
    assert net.initial_state[benenson_input("abb", 1)] == 10
    ids = {s.id: s.structure for s in net.species}
    assert ids[TRANSITIONS["S0", "a"]] == "GGATGTAC\nCCTACATGCCGA"
    # recombination consumes its reactants; Fok I is released by rule (5)
    fok_uses = [lab for rx, lab in zip(net.reactions, net.extras["reaction_rules"]) if FOK_I in rx.reactants]
    assert set(fok_uses) == {"(3)"} and not any(rx.catalysts for rx in net.reactions)


def test_automaton_accepts_words_with_an_even_number_of_b():
    words = ["", "a", "b", "ab", "bb", "bab", "abba", "babb"]
    net = generate_network(ID, words=words, s1_detector=True)
    assert net.status == "complete"
    reporters = net.extras["analysis"]["reporters"]
    assert reporters == {w: ["S0" if w.count("b") % 2 == 0 else "S1"] for w in words}


def test_soup_run_of_the_published_pool():
    traj = evolve(ID, steps=4000, seed=3)
    net = traj.network
    assert net.status == "observed"
    final = Counter(net.extras["final_state"])
    fok_elements = sum(c * parse_molecule(m)[0][1].count("F") for m, c in final.items())
    assert fok_elements == 900, "the 100 Fok I objects are never destroyed"
    assert net.extras["analysis"]["accepted"] == ["abb"]
    assert final[benenson_reporter("S0", 1)] >= 1
    assert sum(rx.count for rx in net.reactions) <= 4000
    # a frame per generation: the published pool holds 220 objects
    assert traj.times()[:3] == [0.0, 220.0, 440.0] and traj.times()[-1] == 4000.0


# --- transcription (2009 sec. 3.2) -----------------------------------------------------------
def test_transcription_rules_12_to_18():
    r = rules(TRANSCRIPTION_RULES)
    chrom = chromosome("TATATTCGCAATGCTGAGCTAGTTTT")
    assert chrom == "0#OrcTATATTCGCAATGCTGAGCTAGTTTT/0#OrcATATAAGCGTTACGACTCGATCAAAA/"
    [mol] = once(r["(12)"], chrom, "0#Rp/")
    assert mol.endswith("/7#Rp/")                       # Rp sits at the first base after TATATT
    [mol] = once(r["(13)"], mol, "0#Cap/")
    pair = {"A": ("(14)", "0#U/"), "T": ("(15)", "0#A/"), "G": ("(16)", "0#C/"), "C": ("(17)", "0#G/")}
    for base in "GCGTTACGACTCGATCAAAA":                  # lower strand after the promoter
        label, nucleotide = pair[base]
        [mol] = once(r[label], mol, nucleotide)
    assert once(r["(18)"], mol) == [chrom, "0#CapCGCAAUGCUGAGCUAGUUUU/", "0#Rp/"]


def test_transcription_network():
    net = generate_network(ID, system="transcription")
    assert net.status == "complete"
    assert net.extras["analysis"] == {"mrna": "0#CapCGCAAUGCUGAGCUAGUUUU/", "mrna_found": True}
    labels = net.extras["reaction_rules"]
    assert set(labels) == {"(12)", "(13)", "(14)", "(15)", "(16)", "(17)", "(18)", "source"}
    assert labels.count("source") == 5


# --- fatty acid oxidation (2009 sec. 4, figs. 1 and 3) ---------------------------------------
def cycle(r, chain):
    """Rules (25)-(28) in order; returns the chains after each step and the byproducts."""
    chain1, fadh2 = once(r["(25)"], chain, "0#Fad/")
    [chain2] = once(r["(26)"], chain1, "0#HOH/")
    chain3, nadh, proton = once(r["(27)"], chain2, "0#NadPo/")
    acetyl, chain4 = once(r["(28)"], chain3, "0#CoaSH/")
    assert (fadh2, nadh, proton, acetyl) == ("0#HFadH/", "0#NadH/", "0#HPo/", ACETYL_COA)
    return [chain1, chain2, chain3, chain4]


def test_figure_3_reaction_path():
    r = rules(FATTY_ACID_RULES)
    # fig. 3 read bottom to top, starting from the printed molecule
    first = cycle(r, "0#HHHHHO/0#CCCCCC/0#HHHHHSCoa/")
    assert first == ["0#HHHHHO/0#CCCCCC/0#HHHXXSCoa/", "0#H/-3#HHHOHO/-3#CCCCCC/-3#HHHHHSCoa/",
                     "0#HHHOHO/0#CCCCCC/0#HHHXHSCoa/", "0#HHHO/0#CCCC/0#HHHSCoa/"]
    second = cycle(r, first[-1])
    assert second[:3] == ["0#HHHO/0#CCCC/0#HXXSCoa/", "0#H/-1#HOHO/-1#CCCC/-1#HHHSCoa/",
                          "0#HOHO/0#CCCC/0#HXHSCoa/"]
    # the printed start lacks the -1#H of section 4's form: its last residue is not acetyl CoA
    assert second[3] == "0#HO/0#CC/0#HSCoa/" != ACETYL_COA
    # with the section 4 form, C6 acyl CoA gives three acetyl CoA (the question of sec. 5)
    assert cycle(r, acyl_coa(6))[-1] == acyl_coa(4)
    assert cycle(r, acyl_coa(4))[-1] == ACETYL_COA


def test_fatty_acid_network_conserves_atoms_except_the_vacancy():
    assert acyl_coa(10) == "0#HHHHHHHHHO/-1#HCCCCCCCCCC/0#HHHHHHHHHSCoa/"   # 2009 sec. 4
    assert acyl_coa(2) == ACETYL_COA
    net = generate_network(ID, system="fatty-acid-oxidation")
    assert net.status == "complete"
    assert net.extras["analysis"] == {"acyl_coa_carbons": [2, 4, 6, 8, 10], "acetyl_coa": True}

    def atoms(side):
        total = Counter()
        for s, n in side.items():
            for _, line in parse_molecule(s):
                total.update(line * n)
        total.pop("X", None)
        return total

    for rx in net.reactions:
        assert atoms(rx.reactants) == atoms(rx.products), rx.to_text()
    ids = {s.id for s in net.species}
    assert {"0#HFadH/", "0#NadH/", "0#HPo/", ACETYL_COA} <= ids


# --- parameters ------------------------------------------------------------------------------
def test_bad_parameters():
    with pytest.raises(ValueError, match="words"):
        generate_network(ID, words=["abc"])
    with pytest.raises(ValueError, match="even"):
        generate_network(ID, system="fatty-acid-oxidation", carbons=5)
    with pytest.raises(ValueError, match="dna"):
        generate_network(ID, system="transcription", dna="TAXU")
    with pytest.raises(ValueError, match="only used with system='custom'"):
        generate_network(ID, rules=["0#A/ -> 0#B/"])
    with pytest.raises(ValueError, match="exactly once"):
        generate_network(ID, system="custom", rules=["0#A1/ -> 0#A/"], pool={"0#AB/": 1})
    with pytest.raises(ValueError, match="twice"):
        generate_network(ID, system="custom", rules=["0#A1/ + 0#B1/ -> 0#A1/ + 0#B/"], pool={"0#AB/": 1})
    with pytest.raises(ValueError, match="notation"):
        generate_network(ID, system="custom", rules=["0#A/ -> 0#B/"], pool={"AB": 1})
    with pytest.raises(ValueError, match="line"):
        generate_network(ID, system="custom", rules=["0#A1*B/ -> 0#AB1*/"], pool={"0#AB/": 1})
