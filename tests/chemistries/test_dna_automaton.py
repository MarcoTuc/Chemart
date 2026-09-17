"""The DNA automaton reproduces Benenson et al., PNAS 100:2191-2196 (2003).

Everything checked here is published: the oligonucleotides and the transition
table of fig. 2, the inputs of fig. 1D and Materials and Methods, the expected
outputs of fig. 3B, the output band lengths, the software reuse of fig. 3C and
the throughput figures of the abstract. The book's automaton (fig. 19.13) is
checked as a property.
"""

from collections import Counter

import pytest

import chemart
from chemart.chemistries import dna_automaton as m
from chemart.kinetics import AVOGADRO

# --- published oligonucleotides (PNAS 2003, Materials and Methods) ------------
SENSE_OLIGOS = {
    "TN1368": "AAGAGCTAGAGTCGGATGC",
    "TN24": "AAGAGCTAGAGTCGGATGCC",
    "TN57": "AAGAGCTAGAGTCGGATG",
}
SOFTWARE = {   # name: (sense oligo, antisense oligo, rule of fig. 2C)
    "T1": ("TN1368", "AGCCGCATCCGACTCTAGCTCT", ("S0", "a", "S0")),
    "T2": ("TN24", "AGCCGGCATCCGACTCTAGCTCT", ("S0", "a", "S1")),
    "T3": ("TN1368", "CCTGGCATCCGACTCTAGCTCT", ("S0", "b", "S0")),
    "T4": ("TN24", "CCTGGGCATCCGACTCTAGCTCT", ("S0", "b", "S1")),
    "T5": ("TN57", "GCCACATCCGACTCTAGCTCT", ("S1", "a", "S0")),
    "T6": ("TN1368", "GCCAGCATCCGACTCTAGCTCT", ("S1", "a", "S1")),
    "T7": ("TN57", "CTGCCATCCGACTCTAGCTCT", ("S1", "b", "S0")),
    "T8": ("TN1368", "CTGCGCATCCGACTCTAGCTCT", ("S1", "b", "S1")),
}
INPUT_OLIGOS = {   # complete input molecules
    "abb": ("GGCTGCCGCAGGGCCGCAGGGCCGTCGGTACCGATTAAGTTGGA",
            "CCAACTTAATCGGTACCGACGGCCCTGCGGCCCTGCGGC"),
    "abba": ("GGCTGCCGCAGGGCCGCAGGGCCTGGCTGCCGTCGGTACCGATTAAGTTGGA",
             "CCAACTTAATCGGTACCGACGGCAGCCAGGCCCTGCGGCCCTGCGGC"),
}
INPUT_BLOCKS = {   # building blocks ligated into the longer inputs
    "babb2": ("CAGGGCCTGGCTGCCGCAGGGCCGCAGGGCCT", "AGCCAGGCCCTGCGGCCCTGCGGCAGCCAGGC"),
    "baaa2": ("CAGGGCCTGGCTGCCTGGCTGCCTGGCTGCCT", "AGCCAGGCAGCCAGGCAGCCAGGCAGCCAGGC"),
    "abbb3": ("GGCTGCCGCAGGGCCGCAGGGCCGCAGGGCCG", "CCTGCGGCCCTGCGGCCCTGCGGCCCTGCGGC"),
}
LIGATIONS = {
    "I3": ["babb2", "abb"], "I4": ["babb2", "abba"], "I5": ["baaa2", "abb"],
    "I6": ["baaa2", "abba"], "I7": ["abbb3", "babb2", "abb"],
    "I8": ["abbb3", "baaa2", "abba"],
}
# fig. 3B: expected output state of automata A1-A3 on inputs I1-I8 (0 = S0).
FIG_3B = {
    "A1": [1, 0, 0, 1, 0, 1, 1, 0],
    "A2": [1, 0, 1, 0, 1, 0, 1, 0],
    "A3": [1, 0, 1, 0, 1, 0, 1, 0],
}


@pytest.fixture(scope="module")
def net():
    return chemart.generate_network("dna-automaton", seed=1)


def reference(program: str, symbols: str) -> str:
    """The abstract automaton of fig. 1A/1C, computed independently."""
    if program == "A1":
        return "S0" if symbols.count("a") % 2 == 0 else "S1"
    if program == "A2":
        return "S0" if len(symbols) % 2 == 0 else "S1"
    if program == "A3":
        return "S1" if symbols.endswith("b") else "S0"
    return "S0" if symbols.count("b") % 2 == 0 else "S1"   # book fig. 19.13


# --- the molecules (fig. 2) ---------------------------------------------------
def test_symbol_and_state_encoding():
    # fig. 2A: 5-bp symbols; the leftmost 4-nt window means S1, the rightmost S0
    assert m.SYMBOL_DUPLEX == {"a": "TGGCT", "b": "GCAGG", "t": "GTCGG"}
    assert m.SYMBOL_SPACER == "GCC" and len(m.SYMBOL_SPACER) == 3
    assert (m.sticky_end("S1", "a"), m.sticky_end("S0", "a")) == ("TGGC", "GGCT")
    assert (m.sticky_end("S1", "b"), m.sticky_end("S0", "b")) == ("GCAG", "CAGG")
    assert (m.sticky_end("S1", "t"), m.sticky_end("S0", "t")) == ("GTCG", "TCGG")


@pytest.mark.parametrize("name", sorted(SOFTWARE))
def test_software_molecules_are_the_published_duplexes(name):
    sense_oligo, antisense, rule = SOFTWARE[name]
    sense = SENSE_OLIGOS[sense_oligo]
    assert m.software_duplex(name) == (sense, antisense)
    # a duplex with a 4-nt 5' antisense overhang (the detector) and a 1-nt 5'
    # A extrusion on the sense strand (fig. 3A)
    assert antisense[4:] == m.revcomp(sense)[:len(sense) - 1]
    assert sense[0] == "A"
    # the rule Chemart reads off the sequences is the one of fig. 2C
    derived = m.TRANSITIONS[name]
    assert (derived["state"], derived["symbol"], derived["next"]) == rule
    assert m.revcomp(derived["sticky_end"]) == m.sticky_end(rule[0], rule[1])
    # the spacer is what changes the state: 0 bp S1->S0, 1 bp keep, 2 bp S0->S1
    expected_spacer = {("S0", "S0"): 1, ("S0", "S1"): 2, ("S1", "S0"): 0, ("S1", "S1"): 1}
    assert derived["spacer"] == expected_spacer[(rule[0], rule[2])]
    assert sense.index(m.FOKI_SITE) + 5 + derived["spacer"] == len(sense)


def test_all_eight_two_state_two_symbol_transitions_exist():
    # book 19.3.1: "it is possible to implement any of the eight possible state
    # transitions for a two-state two-symbol FSM"
    rules = {(r["state"], r["symbol"], r["next"]) for r in m.TRANSITIONS.values()}
    assert rules == {(s, x, n) for s in ("S0", "S1") for x in ("a", "b") for n in ("S0", "S1")}


def test_published_input_molecules():
    # fig. 1D and Materials and Methods: abb and abba are synthesised whole ...
    for symbols, oligos in INPUT_OLIGOS.items():
        assert m.input_molecule("S0", symbols) == oligos
    # ... and I3-I8 are ligations of the published building blocks
    parts = dict(INPUT_BLOCKS, **INPUT_OLIGOS)
    for name, blocks in LIGATIONS.items():
        sense = "".join(parts[b][0] for b in blocks)
        antisense = "".join(parts[b][1] for b in reversed(blocks))
        assert m.input_molecule("S0", m.PUBLISHED_INPUTS[name]) == (sense, antisense)


def test_input_molecules_have_a_sticky_end_and_an_extrusion():
    for symbols in m.PUBLISHED_INPUTS.values():
        sense, antisense = m.input_molecule("S0", symbols)
        assert sense[:4] == m.sticky_end("S0", symbols[0])
        assert antisense == m.revcomp(sense[4:-1])      # 4-nt overhang, 1-nt A extrusion
        assert sense[-1] == "A"
        # 4-nt sticky end + (5 bp symbol + 3 bp spacer) per symbol + terminator + tail
        assert len(sense) == 4 + 8 * len(symbols) + len(m.TERMINATOR_TAIL)


# --- FokI (fig. 2B) -----------------------------------------------------------
def test_fok_i_cut_positions():
    assert (m.FOKI_SITE, m.FOKI_SENSE_CUT, m.FOKI_ANTISENSE_CUT) == ("GGATG", 9, 13)
    sense, _ = m.input_molecule("S0", "abb")
    for spacer in (0, 1, 2):
        fragment_sense, fragment_antisense, rest = m.cleave(sense, spacer)
        # In the hybrid the sense strand runs GGATG + spacer + input, so the cut
        # 9 nt downstream of the site falls 9 - spacer bases into the input ...
        assert spacer + len(fragment_sense) == 9
        # ... while on the antisense strand the software's 4-nt detector overhang
        # comes first (the two are only nicked, never ligated), so the cut 13 nt
        # downstream falls 9 - spacer bases into the input's antisense strand.
        assert spacer + 4 + len(fragment_antisense) == 13
        assert fragment_antisense == m.revcomp(sense[4:13 - spacer])
        # the scattered symbol is a 7-9 nt duplex with only 3-5 bp paired
        assert 3 <= len(fragment_sense) - 4 <= 5
        # and the remaining molecule keeps a 4-nt 5' sticky end
        assert rest == sense[9 - spacer:]
        assert m.revcomp(rest[4:-1]) == m.revcomp(sense[13 - spacer:-1])


def test_every_transition_keeps_the_molecule_well_formed():
    for name, rule in m.TRANSITIONS.items():
        symbols = rule["symbol"] + "ab"
        sense, _ = m.input_molecule(rule["state"], symbols)
        assert m.revcomp(sense[:4]) == rule["sticky_end"], name
        *_, rest = m.cleave(sense, rule["spacer"])
        assert rest == m.input_molecule(rule["next"], symbols[1:])[0], name


# --- the computations (fig. 3B) ----------------------------------------------
@pytest.mark.parametrize("program", ["A1", "A2", "A3"])
def test_fig_3b_expected_outputs(program):
    rules = {n: m.TRANSITIONS[n] for n in m.PROGRAMS[program]["transitions"]}
    states = []
    for i in range(1, 9):
        symbols = m.PUBLISHED_INPUTS[f"I{i}"]
        trace = m.run(rules, symbols)
        assert not trace["suspended"]
        assert len(trace["steps"]) == len(symbols)
        assert trace["final_state"] == reference(program, symbols)
        states.append(int(trace["final_state"][1]))
    assert states == FIG_3B[program]


@pytest.mark.parametrize("program", ["A1", "A2", "A3"])
def test_output_band_lengths(program):
    # Materials and Methods: "S0 and S1 outputs were represented by 15-nt and
    # 16-nt bands" on the labelled antisense strand
    rules = {n: m.TRANSITIONS[n] for n in m.PROGRAMS[program]["transitions"]}
    for symbols in m.PUBLISHED_INPUTS.values():
        trace = m.run(rules, symbols)
        assert trace["output_antisense_nt"] == {"S0": 15, "S1": 16}[trace["final_state"]]


def test_the_books_automaton_counts_b_symbols():
    # book fig. 19.13: two states, accept an even number of b's
    net = chemart.generate_network("dna-automaton", program="book",
                                   inputs=list(m.PUBLISHED_INPUTS))
    analysis = net.extras["analysis"]
    assert set(analysis["transition_table"]) == {"T1", "T4", "T6", "T7"}
    assert analysis["accepting_state"] == "S0"
    for run in analysis["runs"]:
        assert run["final_state"] == reference("book", run["symbols"])
        assert run["accepted"] == (run["symbols"].count("b") % 2 == 0)


def test_software_reuse_of_fig_3c():
    # fig. 3C: A1 on I8; T2, T5 and T8 performed 29, 21 and 54 transitions per
    # molecule, i.e. T8 about twice T2 and T5, and T3 is never used
    net = chemart.generate_network("dna-automaton", program="A1", inputs=["I8"],
                                   conditions="reuse")
    uses = net.extras["analysis"]["transitions_per_software_molecule"]
    assert uses == {"T2": 3, "T3": 0, "T5": 3, "T8": 6}
    published = net.extras["published"]["transitions_per_software_molecule_fig3C"]
    order = sorted(uses, key=lambda n: -uses[n])[:3]
    assert order[0] == "T8" and set(order) == set(published)
    assert sorted(published, key=lambda n: -published[n])[0] == "T8"
    assert net.initial_state["T8"] == 0.075 and net.initial_state["FokI"] == 1.0


def test_suspension_when_no_transition_matches():
    # a program with no rule for <S1, b> suspends on the first b read in S1
    net = chemart.generate_network("dna-automaton", transitions=["T2", "T5"],
                                   inputs=["ab"])
    (run,) = net.extras["analysis"]["runs"]
    assert run["suspended"] and run["unprocessed"] == "b"
    assert len(run["transitions"]) == 1 and run["final_state"] == "S1"
    assert not run["accepted"]


# --- the network --------------------------------------------------------------
def test_default_network_is_the_published_reaction(net):
    analysis = net.extras["analysis"]
    assert analysis["program"] == "A1"
    assert net.status == "complete"
    (run,) = analysis["runs"]
    assert (run["input"], run["symbols"]) == ("I8", "abbbbaaaabba")
    assert run["final_state"] == "S0" and run["accepted"]
    steps = len(run["transitions"])
    assert steps == 12
    # 2 assembly reactions per software molecule, 3 reactions per transition
    assert len(net.reactions) == 2 * 4 + 3 * steps
    assert net.extras["conditions"]["experiment"].startswith("fig. 3B")
    assert net.initial_state == {"FokI": 4.0, "T2": 1.0, "T3": 1.0, "T5": 1.0,
                                "T8": 1.0, "S0_abbbbaaaabbat": 1.0}
    assert {"FokI", "FokI-T2", "S0_abbbbaaaabbat", "S0_t"} <= {s.id for s in net.species}


def test_hardware_and_software_are_recycled(net):
    # fig. 2E: the complex is consumed by the hybridisation and released by the
    # cleavage, so the recycling is a cycle, not a species on both sides
    keys = {(frozenset(r.reactants.items()), frozenset(r.products.items()))
            for r in net.reactions}
    cleavages = [r for r in net.reactions if len(r.products) == 3]
    assert len(cleavages) == 12
    for reaction in cleavages:
        (hybrid,) = reaction.reactants
        (loaded,) = [s for s in reaction.products if s.count("-") == 1]
        molecule = hybrid[len(loaded) + 1:]
        assert loaded.startswith("FokI-T") and hybrid == f"{loaded}-{molecule}"
        pair = frozenset({(loaded, 1), (molecule, 1)})
        assert (pair, frozenset({(hybrid, 1)})) in keys        # hybridisation
        assert (frozenset({(hybrid, 1)}), pair) in keys        # dissociation
        assert reaction.rate is None          # no step time published for fig. 3B
    assert not any(r.catalysts for r in net.reactions)
    # FokI itself only ever assembles with software, never with an input
    for reaction in net.reactions:
        assert "FokI" not in reaction.reactants or set(reaction.reactants) == {"FokI", *[
            s for s in reaction.reactants if s.startswith("T")]}


def test_species_carry_the_published_sequences(net):
    for species in net.species:
        if species.id == "FokI":
            assert "GGATG" in species.structure
            continue
        strands = species.structure.split(" + ")[-1].split("/")
        assert len(strands) == 2 and all(set(s) <= set("ACGT") for s in strands)
        if species.id[:2] in ("S0", "S1"):
            sense, antisense = strands
            assert antisense == m.revcomp(sense[4:-1])
            state, symbols = species.id[:2], species.id[3:-1]
            assert (sense, antisense) == m.input_molecule(state, symbols)


def test_nucleotides_are_conserved(net):
    (law,) = net.extras["conservation"]
    assert law["name"] == "nucleotides"
    weight = dict(zip([s.id for s in net.species], law["vector"]))
    assert weight["FokI"] == 0
    for reaction in net.reactions:
        left = sum(weight[s] * n for s, n in reaction.reactants.items())
        right = sum(weight[s] * n for s, n in reaction.products.items())
        assert left == right, reaction.to_text()[:120]


def test_published_rate_constants():
    fast = chemart.generate_network("dna-automaton", conditions="fast")
    rates = {r.rate["k"] for r in fast.reactions if r.rate}
    assert rates == {0.05}                                  # 20 s per transition
    assert fast.initial_state["S0_abbbbaaaabbat"] == 0.01    # 10 nM input
    parallel = chemart.generate_network("dna-automaton", conditions="parallel")
    assert {r.rate["k"] for r in parallel.reactions if r.rate} == {round(1 / 45, 6)}
    assert all(r.rate["law"] == "mass-action" for r in parallel.reactions if r.rate)


def test_parallel_throughput_of_the_abstract():
    # abstract: 3e12 automata per ul performing 6.6e10 transitions per s per ul
    net = chemart.generate_network("dna-automaton", conditions="parallel")
    published = net.extras["published"]
    micromolar = net.extras["conditions"]["input_uM"]
    automata_per_ul = micromolar * 1e-6 * 1e-6 * AVOGADRO
    assert automata_per_ul == pytest.approx(published["automata_per_ul"], rel=0.02)
    per_second = automata_per_ul / net.extras["conditions"]["seconds_per_transition"]
    assert per_second == pytest.approx(published["transitions_per_second_per_ul"], rel=0.01)


def test_shared_intermediates_are_one_species():
    # I1 = abb and I2 = abba reach the same configurations until abba's last a
    net = chemart.generate_network("dna-automaton", program="A1", inputs=["I1", "I2"])
    ids = [s.id for s in net.species]
    assert len(ids) == len(set(ids))
    assert "S1_bbt" in ids and "S1_bbat" in ids
    assert net.initial_state["S0_abbt"] == 1.0 and net.initial_state["S0_abbat"] == 1.0


def test_seed_does_not_change_the_network():
    # the chemistry is deterministic: the wet experiment has no random choice
    a = chemart.generate_network("dna-automaton", seed=1)
    b = chemart.generate_network("dna-automaton", seed=99)
    assert a.to_dict()["species"] == b.to_dict()["species"]
    assert a.to_dict()["reactions"] == b.to_dict()["reactions"]


def test_invalid_parameters():
    with pytest.raises(ValueError, match="T9"):
        chemart.generate_network("dna-automaton", transitions=["T9"])
    with pytest.raises(ValueError, match="nondeterministic"):
        chemart.generate_network("dna-automaton", transitions=["T1", "T2"])
    with pytest.raises(ValueError, match="published input"):
        chemart.generate_network("dna-automaton", inputs=["abc"])
    with pytest.raises(ValueError, match="at least one"):
        chemart.generate_network("dna-automaton", inputs=[])


@pytest.mark.slow
def test_every_published_program_and_input_builds_a_network():
    for program in ("A1", "A2", "A3"):
        net = chemart.generate_network("dna-automaton", program=program,
                                       inputs=list(m.PUBLISHED_INPUTS))
        runs = net.extras["analysis"]["runs"]
        states = [int(run["final_state"][1]) for run in runs]
        assert states == FIG_3B[program]
        # 60 transitions over the eight inputs, but configurations shared by two
        # inputs are one species and therefore one cleavage reaction
        assert sum(len(run["transitions"]) for run in runs) == 60
        cleaved = {step["molecule"] for run in runs for step in run["transitions"]}
        counts = Counter(len(r.products) for r in net.reactions)
        assert counts[3] == len(cleaved) < 60
        weight = dict(zip([s.id for s in net.species],
                          net.extras["conservation"][0]["vector"]))
        for reaction in net.reactions:
            assert (sum(weight[s] * n for s, n in reaction.reactants.items())
                    == sum(weight[s] * n for s, n in reaction.products.items()))
