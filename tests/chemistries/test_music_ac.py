"""Music composition AC: book 16.4, Miura & Tominaga, GWAL-7 (2006) [592].

Every published item of the paper that can be checked: the notation and
recombination examples of section 2, the three design totals and the printed
rules of section 4, the algorithm of section 3, and the phrases of figures 1
and 3.
"""

import itertools

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.music_ac import (
    AVOID_PITCHES, CADENCES, CHORD_TONES, FUNCTION_CHORDS, NOTES, PITCH,
    SUBDOMINANT_TO_TONIC, Rule, describe, elements, initial_pool, is_phrase,
    match, object_text, parse_object, parse_pattern, read_bars, rule_set,
)

ID = "music-ac"

# Figure 1: the example object and the music it represents (C, G, C, G).
FIGURE_1 = ("0#StartE4D4C4D4E4C4G4G4D4E4F4E4F4G4F4E4G4C5G4A4B4C5B4C5G4F4G4F4G4F4E4D4Stop/"
            "1#TCCCCCCCDGGGGGGGTCCCCCCCDGGGGGGG/")
# Figure 3: the three generated phrases, by their chords.
FIGURE_3 = {1: ["Am", "Am", "C", "G"], 2: ["Em", "Dm", "G", "Em"], 3: ["C", "Bm", "Am", "F"]}

# The rules printed in section 4.3, for the default 8 notes per bar / 4 bars.
PRINTED_RULES = {
    "(1) note after D4":
        "0#<*1>D4/ + 0#Degree/0#D4<2>/ -> 0#<*1>D4<2>/ + 0#Degree/0#D4/",
    "(2) E4 -> T Em":
        "0#E4<1..7><8><9*>/ + 0#Chord/0#TEmEmEmEmEmEmEm/ "
        "-> 0#E4<1..7>/0#TEmEmEmEmEmEmEm/ + 0#<8><9*>/ + 0#Chord/",
    "(6) avoid F4 of C":
        "0#<*1>F4<2*>/0#<*3>C<4*>/ + 0#Avoid/0#Dummy/ "
        "-> 0#<*1>Dummy<2*>/0#<*3>C<4*>/ + 0#Avoid/0#F4/",
    "(7) Dummy -> a random note":
        "0#<*1>Dummy<2*>/0#<*3><4><5*>/ + 0#Avoid/0#<6>/ "
        "-> 0#<*1><6><2*>/0#<*3><4><5*>/ + 0#Avoid/0#Dummy/",
    "(8) T-D2":
        "0#<*1><2..9>/0#<*10>T<11..17>/ + 0#<18><19><20*>/0#D2<21*>/ "
        "-> 0#<*1><2..9><18><19><20*>/0#<*10>T<11..17>D2<21*>/",
    "(9) 4-bar phrase":
        "0#<1..32><33><34*>/0#T<35..65><66><67*>/ + 0#StartStop/ "
        "-> 0#Start<1..32>Stop/1#T<35..65>/ + 0#<33><34*>/0#<66><67*>/",
    "(10) detach leading D bar":
        "0#<1..8><9><10*>/0#D<11..17><18><19*>/ "
        "-> 0#<1..8>/0#D<11..17>/ + 0#<9><10*>/0#<18><19*>/",
}


def texts(outcome):
    return [tuple(object_text(o) for o in products) for products in outcome]


def test_notation_and_recombination_of_section_2():
    """Paper section 2: the matching examples and the worked recombination."""
    assert parse_object("0#CDE/") == ((0, ("C", "D", "E")),)
    assert parse_object("0#CDEF/3#EmAm/") == ((0, ("C", "D", "E", "F")), (3, ("Em", "Am")))
    assert object_text(parse_object("0#CBA/1#Am/")) == "0#CBA/1#Am/"
    # "the pattern 0#C<1>E/ matches the object 0#CDE/"
    assert match(parse_pattern("0#C<1>E/"), parse_object("0#CDE/")) == [{"1": ("D",)}]
    # "the pattern 0#<*1>EF/1#EmAm/ matches the object 0#CDEF/3#EmAm/"
    assert match(parse_pattern("0#<*1>EF/1#EmAm/"), parse_object("0#CDEF/3#EmAm/"))
    # displacements are meaningful: the same lines at another offset do not match
    assert not match(parse_pattern("0#<*1>EF/1#EmAm/"), parse_object("0#CDEF/2#EmAm/"))
    # "If this rule is applied to the objects 0#CDEF/3#EmAm/ and 0#GAB/, the object
    #  0#CDEFGAB/3#EmAm/ is produced and the reactants disappear."
    rule = Rule("0#<*1>EF/1#EmAm/ + 0#GA<2*>/ -> 0#<*1>EFGA<2*>/1#EmAm/")
    assert texts(rule.apply("0#CDEF/3#EmAm/", "0#GAB/")) == [("0#CDEFGAB/3#EmAm/",)]


def test_rules_reject_ill_formed_wildcards():
    """Section 2: a wildcard appears once in lhs and once and only once in rhs."""
    with pytest.raises(ValueError, match="twice on the left-hand side"):
        Rule("0#<1><1>/ -> 0#<1><1>/")
    with pytest.raises(ValueError, match="once on each side"):
        Rule("0#<1><2>/ -> 0#<1>/")
    with pytest.raises(ValueError, match="may only open a line"):
        Rule("0#<1><*2><3>/ -> 0#<1><*2><3>/")


def test_design_totals_of_section_4():
    """'24 kinds of elements, 43 kinds of initial objects (1826 objects in total),
    and 65 recombination rules'."""
    listed = [*NOTES, *CHORD_TONES, *FUNCTION_CHORDS, "Degree", "Chord", "Avoid", "Start", "Stop"]
    assert len(listed) == 24                       # the elements the paper lists (section 4.1)
    assert elements() == listed + ["Dummy"]        # its own initial objects also use Dummy
    assert len(elements()) == 25
    pool = initial_pool()
    assert len(pool) == 43
    assert sum(pool.values()) == 1826
    assert len(rule_set()) == 65


def test_rule_groups_add_up_to_65():
    """8 note rules + 38 chord rules + 8 avoid rules + 1 Dummy rule + 6 cadences + 1 + 3."""
    groups = {}
    for rule in rule_set():
        groups.setdefault(rule.label.split(")")[0] + ")", []).append(rule)
    sizes = {key: len(value) for key, value in groups.items()}
    assert sizes == {"(1)": 8, "(2)": 38, "(6)": 8, "(7)": 1, "(8)": 6, "(9)": 1, "(10)": 3}
    assert sum(sizes.values()) == 65
    # the avoid rules of G and Bm need one each for C4 and C5 ("C means both C4 and C5")
    assert sum(len(p) for p in AVOID_PITCHES.values()) == 8
    # dropping the avoid-note step leaves 56 rules; adding S-T makes 66
    assert len(rule_set(avoid_notes=False)) == 56
    assert len(rule_set(cadences="with-s-t")) == 66


def test_initial_multiset_of_section_4_2():
    """The published counts: 100, 20, 10, 10, 100 and 2 per kind."""
    pool = initial_pool()
    conjunct = {k: v for k, v in pool.items() if k.startswith("0#Degree/")}
    assert len(conjunct) == 14 and set(conjunct.values()) == {100}
    assert pool["0#Degree/0#C4D4/"] == 100        # "the motion from C4 to D4"
    chords = {k: v for k, v in pool.items() if k.startswith("0#Chord/")}
    assert len(chords) == 11 and set(chords.values()) == {20}
    assert pool["0#Chord/0#TCCCCCCC/"] == 20      # "the tonic chord C"
    assert pool["0#Chord/0#SAmAmAmAmAmAmAm/"] == 20   # Am is tonic and subdominant
    assert pool["0#Chord/0#TAmAmAmAmAmAmAm/"] == 20
    assert pool["0#Avoid/0#Dummy/"] == 10
    assert {pool[f"0#Avoid/0#{n}/"] for n in NOTES} == {10}
    assert pool["0#StartStop/"] == 100
    assert {pool[f"0#{n}/"] for n in NOTES} == {2}


@pytest.mark.parametrize("label", sorted(PRINTED_RULES))
def test_printed_rules_are_reproduced(label):
    """The rules printed in sections 4.3, character for character."""
    rule = next(r for r in rule_set() if r.label == label)
    assert rule.text == PRINTED_RULES[label]


def test_rule_1_grows_a_conjunct_melody():
    """Section 3.1 / rule (1): the next note is just above or below the previous one."""
    rule = next(r for r in rule_set() if r.label == "(1) note after D4")
    assert texts(rule.apply("0#C4D4/", "0#Degree/0#D4E4/")) == [("0#C4D4E4/", "0#Degree/0#D4/")]
    assert texts(rule.apply("0#C4D4/", "0#Degree/0#D4C4/")) == [("0#C4D4C4/", "0#Degree/0#D4/")]
    # the conjunct-motion objects are exactly the scale steps, so no rule can jump
    pool = initial_pool()
    for text in (t for t in pool if t.startswith("0#Degree/")):
        first, second = parse_object(text)[1][1]
        assert abs(NOTES.index(first) - NOTES.index(second)) == 1
    # "if E4 is chosen as the first note, the next note is either D4 or F4"
    rule_e4 = next(r for r in rule_set() if r.label == "(1) note after E4")
    grown = set()
    for text in (t for t in pool if t.startswith("0#Degree/")):
        for products in rule_e4.apply("0#E4/", text):
            grown.add(object_text(products[0]))
    assert grown == {"0#E4D4/", "0#E4F4/"}


def test_chord_is_chosen_by_the_first_note_of_the_bar():
    """Section 3.2: 'if the first note of the sequence is E4, a chord is chosen for this
    sequence from among C, Em and Am'."""
    rules = [r for r in rule_set() if r.label.startswith("(2)")]
    for_e4 = [r for r in rules if r.label.startswith("(2) E4 ->")]
    assert {r.label.split()[-1] for r in for_e4} == {"C", "Em", "Am"}
    # five rules, because Em and Am each carry two functions (the paper prints four of them)
    assert len(for_e4) == 5
    assert sorted(r.label for r in for_e4) == [
        "(2) E4 -> D Em", "(2) E4 -> S Am", "(2) E4 -> T Am", "(2) E4 -> T C", "(2) E4 -> T Em"]
    # every chord rule gives a bar whose chord contains the bar's first note
    for rule in rules:
        note, _, function, chord = rule.label.replace("(2) ", "").split()
        assert PITCH[note] in CHORD_TONES[chord]
        assert chord in FUNCTION_CHORDS[function]
    # rule (2) cuts the first eight notes into a bar and leaves the rest
    rule = next(r for r in rules if r.label == "(2) E4 -> T Em")
    melody = "0#" + "".join(["E4", "F4", "G4", "F4", "E4", "D4", "E4", "F4", "G4", "A4"]) + "/"
    assert texts(rule.apply(melody, "0#Chord/0#TEmEmEmEmEmEmEm/")) == [(
        "0#E4F4G4F4E4D4E4F4/0#TEmEmEmEmEmEmEm/", "0#G4A4/", "0#Chord/")]


def test_avoid_note_is_replaced_via_dummy():
    """Section 3.3 / rules (6) and (7): an avoid note becomes Dummy, then a random note."""
    bar = "0#C4D4E4F4G4A4B4C5/0#TCCCCCCC/"            # F4 is the avoid pitch of C
    six = next(r for r in rule_set() if r.label == "(6) avoid F4 of C")
    assert texts(six.apply(bar, "0#Avoid/0#Dummy/")) == [(
        "0#C4D4E4DummyG4A4B4C5/0#TCCCCCCC/", "0#Avoid/0#F4/")]
    seven = next(r for r in rule_set() if r.label == "(7) Dummy -> a random note")
    out = texts(seven.apply("0#C4D4E4DummyG4A4B4C5/0#TCCCCCCC/", "0#Avoid/0#E4/"))
    assert out == [("0#C4D4E4E4G4A4B4C5/0#TCCCCCCC/", "0#Avoid/0#Dummy/")]
    # an avoid pitch is never a tone of its chord, so a bar's first note is never replaced
    for chord, pitches in AVOID_PITCHES.items():
        for pitch in pitches:
            assert PITCH[pitch] not in CHORD_TONES[chord]
    # the chord F has no avoid pitch
    assert "F" not in AVOID_PITCHES


def test_rule_10_detaches_and_reuses_a_leading_non_tonic_bar():
    """Section 4.3 / figure 2 (10)-(12): a phrase cannot start with D, so the bar is cut off."""
    first = ["G4", "F4", "E4", "D4", "E4", "F4", "G4", "A4"]
    second = ["C4", "D4", "E4", "F4", "E4", "D4", "C4", "D4"]
    sequence = "0#" + "".join(first + second) + "/0#" + "DGGGGGGG" + "TCCCCCCC" + "/"
    rule = next(r for r in rule_set() if r.label == "(10) detach leading D bar")
    assert texts(rule.apply(sequence)) == [(
        "0#" + "".join(first) + "/0#DGGGGGGG/",
        "0#" + "".join(second) + "/0#TCCCCCCC/")]


def test_figure_1_phrase_is_a_final_product():
    """Figure 1 and section 5: the four-bar phrase C, G, C, G, terminated by Start and
    Stop, 'is no longer recombined by any rule'."""
    obj = parse_object(FIGURE_1)
    assert object_text(obj) == FIGURE_1
    assert is_phrase(obj)
    bars = read_bars(obj)
    assert [b["chord"] for b in bars] == ["C", "G", "C", "G"]
    assert [b["function"] for b in bars] == ["T", "D", "T", "D"]
    assert [len(b["notes"]) for b in bars] == [8, 8, 8, 8]
    assert bars[0]["notes"][0] == "E4" and bars[-1]["notes"][-1] == "D4"
    # every bar's chord contains its first note (section 3.2)
    assert all(PITCH[b["notes"][0]] in CHORD_TONES[b["chord"]] for b in bars)
    # it carries no avoid note, and no rule of the description matches it
    assert not [n for b in bars for n in b["notes"] if n in AVOID_PITCHES.get(b["chord"], ())]
    rules = rule_set()
    assert not [r.label for r in rules if any(match(p, obj) for p in r.lhs)]
    pool = initial_pool()
    for rule in rules:
        if rule.arity == 1:
            assert not rule.apply(obj)
        else:
            assert not [t for t in pool if rule.apply(obj, parse_object(t))]
    assert describe(obj).startswith("phrase: C(T) E4 D4 C4 D4 E4 C4 G4 G4 | G(D)")


def derivable(functions, pairs):
    return all(pair in pairs for pair in zip(functions, functions[1:]))


def test_figure_3_phrases_against_the_cadence_rules():
    """Figure 3 against section 3.4: 3(2) and 3(3) follow the six published cadence
    rules; 3(1) (Am, Am, C, G) needs S-T, which only the cadence T-S-T implies."""
    published, with_st = set(CADENCES), set(CADENCES) | {SUBDOMINANT_TO_TONIC}

    def assignments(chords, pairs):
        options = [[f for f, cs in FUNCTION_CHORDS.items() if c in cs] for c in chords]
        return [a for a in itertools.product(*options)
                if a[0] == "T" and derivable(list(a), pairs)]      # a phrase starts on a tonic

    assert assignments(FIGURE_3[2], published) == [("T", "D2", "D", "T")]   # a typical cadence
    assert ("T", "D", "T", "S") in assignments(FIGURE_3[3], published)
    assert assignments(FIGURE_3[1], published) == []
    assert assignments(FIGURE_3[1], with_st) == [("T", "S", "T", "D")]
    # figure 1 as well
    assert derivable(["T", "D", "T", "D"], published)
    # the typical cadences of section 3.4 are derivable once S-T is present
    for cadence in (["T", "D", "T"], ["T", "S", "T"], ["T", "D2", "D", "T"]):
        assert derivable(cadence, with_st)


def assert_phrase_is_well_formed(phrase, p_bars=4, p_notes=8, pairs=set(CADENCES)):
    assert len(phrase["cadence"]) == p_bars
    assert len(phrase["notes"]) == p_bars * p_notes
    assert phrase["cadence"][0] == "T"                       # section 3.5: starts on a tonic bar
    assert derivable(phrase["cadence"], pairs)               # section 3.4
    bars = read_bars(parse_object(phrase["object"]), p_notes)
    for bar, chord, function in zip(bars, phrase["chords"], phrase["cadence"]):
        assert chord in FUNCTION_CHORDS[function]
        assert PITCH[bar["notes"][0]] in CHORD_TONES[chord]   # section 3.2
        assert len(bar["notes"]) == p_notes


def test_default_run_composes_a_phrase():
    """A run of the system composes a phrase in the style of sections 3 and 6."""
    net = generate_network(ID, seed=1)
    assert net.status == "observed"
    analysis = net.extras["analysis"]
    assert (analysis["elements"], analysis["initial_object_kinds"],
            analysis["initial_objects"], analysis["rules"]) == (25, 43, 1826, 65)
    assert analysis["phrase_count"] >= 1
    for phrase in analysis["phrases"]:
        assert_phrase_is_well_formed(phrase)
    assert all(r.count >= 1 for r in net.reactions)
    assert sum(net.initial_state.values()) == 1826
    # the reactions that fired are applications of the description's rules
    labels = {r.label for r in rule_set()}
    assert set(net.extras["reaction_rules"]) <= labels
    assert len(net.extras["rules"]) == 65


def test_a_frame_per_collision_until_the_first_phrase():
    """chemart.evolve: a frame per collision; the run stops once a phrase is finished."""
    traj = evolve(ID, seed=1)
    assert traj.clock == "collisions"
    assert traj.times() == [float(i) for i in range(len(traj.frames))]
    assert traj.frames[0].fired == [] and traj.frames[0].observables == {"phrases": 0}
    assert all(sum(n for _, _, n in f.fired) <= 1 for f in traj.frames)
    finished = {}
    for f in traj.frames:
        count = sum(finished.setdefault(s, is_phrase(parse_object(s))) for s in f.state)
        assert f.observables["phrases"] == count
    phrases = traj.series("phrases")
    assert phrases[-1] == 1 and not any(phrases[:-1])
    assert traj.frames[-1].state == {s: float(n) for s, n in traj.network.extras["final_state"].items()}


def test_rules_conserve_elements():
    """Section 2: a recombination rule conserves elements 'just like a chemical
    reaction does', so m^T (P - R) = 0 for the count of every element."""
    net = generate_network(ID, seed=2)
    ids, R, P = net.matrices()
    net_stoich = (P - R).toarray()
    laws = {law["name"]: law["vector"] for law in net.extras["conservation"]}
    assert len(laws) >= 20
    for name, vector in laws.items():
        assert len(vector) == len(ids)
        assert not any(sum(v * s for v, s in zip(vector, column)) for column in net_stoich.T), name
    # the melody elements are among the conserved quantities
    assert "element:C4" in laws and "element:Start" in laws


def test_avoid_notes_and_cadences_change_the_description():
    """The two filters of the book's 'stricter rules' are switchable."""
    without = generate_network(ID, seed=3, avoid_notes=False)
    assert len(without.extras["rules"]) == 56
    assert not [r for r in without.extras["reaction_rules"] if r.startswith(("(6)", "(7)"))]
    assert without.extras["analysis"]["phrase_count"] >= 1
    for phrase in without.extras["analysis"]["phrases"]:
        assert_phrase_is_well_formed(phrase)
        # nothing perturbs the melody now, so every interval inside a bar is a
        # scale step: the jumps of section 5 come only from replaced avoid notes
        for bar in read_bars(parse_object(phrase["object"])):
            assert {abs(NOTES.index(a) - NOTES.index(b))
                    for a, b in zip(bar["notes"], bar["notes"][1:])} == {1}
    with_st = generate_network(ID, seed=3, cadences="with-s-t")
    assert len(with_st.extras["rules"]) == 66
    for phrase in with_st.extras["analysis"]["phrases"]:
        assert_phrase_is_well_formed(phrase, pairs=set(CADENCES) | {SUBDOMINANT_TO_TONIC})


def test_shorter_bars_and_phrases():
    """The bar size and the phrase length size rules (2), (8), (9) and (10)."""
    net = generate_network(ID, seed=0, notes_per_bar=4, bars_per_phrase=2)
    assert net.extras["analysis"]["initial_objects"] == 1826
    for phrase in net.extras["analysis"]["phrases"]:
        assert_phrase_is_well_formed(phrase, p_bars=2, p_notes=4)
    rule = next(r for r in net.extras["rules"] if r.startswith("(9)"))
    assert rule.startswith("(9) 2-bar phrase: 0#<1..8><9><10*>/0#T<11..17><18><19*>/")


@pytest.mark.parametrize("given, message", [
    # the polyphonic system of [854] is not reconstructed, so it is not a choice
    (dict(mode="polyphonic"), r"mode='polyphonic' is invalid"),
    (dict(cadences="all"), r"cadences='all' is invalid"),
    (dict(notes_per_bar=1), "notes_per_bar"),
])
def test_rejects_unsupported_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network(ID, **given)


@pytest.mark.slow
def test_longer_run_composes_several_phrases():
    """Section 5: 'a typical run generates three to five phrases'; several phrases grow
    in parallel in the same reactor, and each one observes the style."""
    net = generate_network(ID, seed=5, steps=3000, phrases=0)
    analysis = net.extras["analysis"]
    assert analysis["phrase_count"] >= 3
    for phrase in analysis["phrases"]:
        assert_phrase_is_well_formed(phrase)
    # "each phrase shown in the Figures has a jump in the melody. This is because
    #  avoid notes are replaced with random notes" (section 5)
    jumps = 0
    for phrase in analysis["phrases"]:
        steps = [abs(NOTES.index(a) - NOTES.index(b))
                 for a, b in zip(phrase["notes"], phrase["notes"][1:])]
        assert 1 in steps                                  # conjunct motion is still visible
        jumps += sum(1 for s in steps if s > 1)
    assert jumps
    assert any(r.startswith(("(6)", "(7)")) for r in net.extras["reaction_rules"])
    assert any(r.startswith("(10)") for r in net.extras["reaction_rules"])
