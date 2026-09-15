"""SAC reproduces book 11.1.3 and Suzuki's ECAL 2003 slides for [823].

Book: table 11.4 and eqs. 11.2-11.3. Slides (2003.09_ECAL.presentation.pdf):
reaction examples pp. 10-13, genetic information p. 15, ancestral string sets
of models (i)-(iii) pp. 16-18, spindle and membrane p. 19.
"""

from collections import Counter

import pytest

from chemart import generate_network
from chemart.chemistries.sac import (
    MEMBRANE_BREEDER, MEMBRANE_SEED, MODEL_I_CONSTRUCTOR, MODEL_I_COPIER, MODEL_II_CONSTRUCTOR,
    MODEL_II_COPIER, MODEL_III_CONSTRUCTOR, MODEL_III_COPIER, ancestral_cell, decode, escape,
    membrane_m, react, rewrite, sid, split, translate,
)

P0 = r"""\0&\!$0'0"1*0'1"2"""


# --- the language ----------------------------------------------------------------
def test_book_eq_11_2_and_11_3():
    # [\0]&[\!]$[00*01 -> 01*02]: a 0 is found, no ! is found, so the third rule rewrites
    assert decode(P0) == r"[\0]&[\!]$[00*01 → 01*02]"
    assert rewrite(P0, "002301") == "012302"            # book eq. 11.3
    assert react(P0, "003301") == (P0, "013302")         # slides pp. 10-11
    # the guards: an operand containing ! stops before the rewrite; one without 0 too
    assert react(P0, "00!301") is None
    assert react(P0, "22311") is None


def test_book_erratum_string_as_printed():
    # book eq. 11.2 prints \0!$0'0"1*0'1"2 (no &, ! unsuppressed): rule [0!] matches, $ stops
    printed = r"""\0!$0'0"1*0'1"2"""
    assert decode(printed) == r"[\0!]$[00*01 → 01*02]"
    assert react(printed, "002301") is None


def test_slide_permutation_and_separation():
    assert decode(r'/2/%/"%') == "[2% → %2%]" and rewrite(r'/2/%/"%', "002300") == "0032300"   # p. 12
    assert decode(r'02"\.30') == "[0230 → 02.30]"                                              # p. 13
    assert react(r'02"\.30', "002300") == (r'02"\.30', "002", "300")


def test_wild_cards_star_and_suppression():
    # string rewriting grammar slide, rule 5?*4 -> ?5*?4: ".51234." -> ".152314."
    assert rewrite("\"?5'?*\"?4", "51234") == "152314"      # "? recreates the matched character
    assert rewrite("0*'1", "0221021") == "022021"          # * is the shortest match: first 1 deleted
    assert rewrite("!!'0", "12330") == "1233"               # repeated ! must match the same character
    assert split(r"12.3\.4.") == ["12", r"3\.4"]           # suppressed '.' does not cut
    assert rewrite(r"\*'1", "*1") == "*"


def test_genotype_is_suppressed_phenotype():
    # slide p. 15
    assert escape(r"""\0&\!$'0"1*0'1"2""") == r"""\\\0\&\\\!\$\'\0\"\1\*\0\'\1\"\2"""
    assert translate(r"""\\\0\&\\\!\$\'\0\"\1\*\0\'\1\"\2""") == r"""\0&\!$'0"1*0'1"2"""


# slide p. 16, transcribed independently of the enzymes
MODEL_I_GENES = [
    r"""0000\0\0\'\0\'\0\"\1\"\1\"\0\"\2\"\0\"\2\*\0\0\1\1\"\0\"\2\"\0\"\20011""",
    r"""0000\/\\\0\\\2\\\0\\\2\/\\\\\%\/\*\"\\\\\"\%\0\2\0\20011""",
    r"""0000\0\0\'\1\'\1\"\0\"\0\*\'\0\'\2\'\0\'\2\0\0\1\1\"\\\.\"\0\"\0\"\0\"\0\*\'\0\'\2\'\0\'\2\"\0\"\0\"\1\"\10011""",
    r"""0000\0\0\'\0\'\0\"\1\"\1\"\0\"\2\"\2\"\0\*\0\0\1\1\"\0\"\2\"\2\"\00011""",
    r"""0000\/\\\0\\\2\\\2\\\0\/\\\\\%\/\*\"\%\0\2\2\00011""",
    r"""0000\0\0\'\1\'\1\"\0\"\0\*\'\0\'\2\'\2\'\0\0\0\1\1\"\\\.\*\'\0\'\2\'\2\'\00011""",
]
# slide p. 17, the single chromosome as printed (eight lines)
MODEL_II_CHROMOSOME = "".join([
    r"""0000\\\0\\\0\'\0\'\0\"\1\"\1\"\0\"\2\"\0\"\2\*\3\3\0\0\"\0\"\2\"\0\"\201010110\/""",
    r"""\\\0\\\2\\\0\\\2\/\\\\\%\/\*\"\\\\\"\%\\\0\\\2\\\0\\\201010110\/\\\0\\\2\\\0\\\2\/""",
    r"""\!\/\*\"\!\\\0\\\2\\\0\\\201010110\\\0\\\0\'\1\'\1\"\0\"\0\*\'\0\'\2\'\0\'\2\3\3""",
    r"""\0\0\"\\\.\"\0\"\0\"\0\"\0\*\'\0\'\2\'\0\'\2\"\3\"\3\"\0\"\001010110\\\0\\\0\'\0\'""",
    r"""\0\"\1\"\1\"\0\"\2\"\2\"\0\*\3\3\0\0\"\0\"\2\"\2\"\001010110\/\\\0\\\2\\\2\\\0\/""",
    r"""\\\\\%\/\*\"\%\\\0\\\2\\\2\\\001010110\/\\\0\\\2\\\2\\\0\/\0\1\0\1\*\0\1\1\0\/\*\3""",
    r"""\3\0\0\/\"\\\.\*\/\\\0\\\2\\\2\\\0\/01010110\\\0\\\0\'\1\'\1\"\0\"\0\*\'\0\'\2\'""",
    r"""\2\'\0\3\3\0\0\"\\\.\*\'\0\'\2\'\2\'\03300""",
])


def test_published_genes_match_their_enzymes():
    assert ancestral_cell("independent-genes")["genes"] == MODEL_I_GENES
    assert ancestral_cell("single-chromosome")["genes"] == [MODEL_II_CHROMOSOME]
    assert len(MODEL_I_GENES) + len(MODEL_I_COPIER) + len(MODEL_I_CONSTRUCTOR) == 12   # book: 12 strings


# --- replication and translation ----------------------------------------------------
def drive(enzymes, gene, limit=3000):
    """Let the first applicable enzyme act on the growing gene until the gene reappears."""
    s, released = gene, []
    for _ in range(limit):
        for e in enzymes:
            out = react(e, s)
            if out:
                pieces = list(out[1:])
                break
        else:
            raise AssertionError(f"no enzyme applies to {s!r}")
        if any(gene in q for q in pieces):
            return pieces + released
        s, released = pieces[0], released + pieces[1:]
    raise AssertionError("did not finish")


@pytest.mark.parametrize("i", range(6))
def test_model_i_copier_and_constructor(i):
    gene = MODEL_I_GENES[i]
    assert drive(MODEL_I_COPIER, gene) == [gene, gene]
    assert drive(MODEL_I_CONSTRUCTOR, gene) == [gene, (MODEL_I_COPIER + MODEL_I_CONSTRUCTOR)[i]]


def test_constructor_is_universal():
    # slide p. 15: the constructor can create any operator string, e.g. the book's P0
    gene = "0000" + escape(P0) + "0011"
    assert drive(MODEL_I_CONSTRUCTOR, gene) == [gene, P0]


def test_model_ii_single_chromosome():
    gene = MODEL_II_CHROMOSOME
    c = MODEL_II_COPIER
    # the copy is cut off at the terminator when the last copier string acts before the
    # single-character mover (see decisions)
    assert drive([c[0], c[3], c[1], c[2]], gene) == [gene, gene]
    out = drive(MODEL_II_CONSTRUCTOR, gene)
    assert out[0] == gene and Counter(out[1:]) == Counter(MODEL_II_COPIER + MODEL_II_CONSTRUCTOR)


def test_model_iii_spindle_and_membrane():
    cell = ancestral_cell("spindle-membrane")
    copier = [MODEL_III_COPIER[i] for i in (0, 1, 3)]
    for gene, enzyme in zip(cell["genes"], MODEL_III_COPIER + MODEL_III_CONSTRUCTOR):
        out = drive(copier, gene)
        # slide p. 19: tagged copies for the two daughters, membrane seed and breeder
        assert out[:4] == ["L" + gene, "R" + gene, MEMBRANE_SEED, MEMBRANE_BREEDER]
        # replication is regulated: tagged genes are not copied again
        assert react(MODEL_III_COPIER[0], "L" + gene) is None
        assert react(MODEL_III_COPIER[0], "R" + gene) is None
        assert drive(MODEL_III_CONSTRUCTOR, gene) == [gene, enzyme]
    assert decode(MEMBRANE_BREEDER) == "[EM → MEM]"
    membrane = MEMBRANE_SEED
    for _ in range(10):
        membrane = react(MEMBRANE_BREEDER, membrane)[1]
    assert membrane == "MMMMMMMMMMEM" and membrane_m(membrane) == 10   # "ten Ms": division


# --- networks ---------------------------------------------------------------------
def test_default_closure_reproduces_the_12_strings():
    net = generate_network("sac")
    assert net.status == "complete"
    seed = net.extras["seed"]
    assert len(seed) == 12 and sorted(net.extras["analysis"]["reproduced_seed"]) == sorted(seed)
    assert all(r.catalysts for r in net.reactions), "the operator is always put back"
    ids = {s.id: s.structure for s in net.species}
    assert all(sid(st) == i for i, st in ids.items())
    g = sid(MODEL_I_GENES[0])
    doubling = [r for r in net.reactions if r.products.get(g) == 2]
    assert doubling, "some reaction releases two copies of a gene"


def test_closure_truncates_on_budget():
    net = generate_network("sac", max_species=50)
    assert net.status == "truncated" and len(net.species) <= 50


def test_soup_grows_the_cell():
    # a copy needs ~35 consecutive effective collisions among 12 strings; 20000 collisions
    # grow the ancestral cell of model (i) to about three times its size
    net = generate_network("sac", seed=3, method="soup", steps=20000)
    a = net.extras["analysis"]
    assert net.status == "observed" and a["population_size"][0] == 12
    assert a["population_size"][-1] >= 24
    assert all(b >= a_ for a_, b in zip(a["population_size"], a["population_size"][1:])), "catalytic: never shrinks"
    assert sum(r.count for r in net.reactions) <= 20000
    final = Counter(net.extras["final_state"])
    assert sum(final.values()) == a["population_size"][-1]
    assert any(c > 1 for c in a["seed_copies_final"].values())


def test_soup_constant_dilution_keeps_size():
    net = generate_network("sac", seed=1, method="soup", steps=500, dilution="constant")
    assert set(net.extras["analysis"]["population_size"]) == {12}


def test_custom_strings():
    net = generate_network("sac", strings=[P0, "002301"])
    assert {s.structure for s in net.species} >= {P0, "002301", "012302"}


def test_bad_parameters():
    with pytest.raises(ValueError, match="SAC strings"):
        generate_network("sac", strings=["00x1"])
    with pytest.raises(ValueError, match="unsuppressed"):
        generate_network("sac", strings=["00.11"])
