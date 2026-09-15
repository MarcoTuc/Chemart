"""The Chemical Abstract Machine reproduces Berry & Boudol's derivations.

Berry & Boudol, "The chemical abstract machine", POPL 1990 (book [113], journal
version TCS 96, 1992): the CCS- execution and non-determinism examples
(sec. 2.2), communication through a restriction membrane with an airlock and a
heavy ion (sec. 2.3), the laws (sec. 3.1), the TCCS cham rules (sec. 4.3) and
the external sum example (sec. 4.2). Book section 9.3.
"""

from collections import Counter

import pytest

from chemart import generate_network
from chemart.chemistries.cham import (
    CCS_MINUS_CLEANUP, CCS_MINUS_RULES, Machine, machine_rules, membranes, parse_molecule,
    parse_rules, parse_solution, solution_text, text, transitions,
)

ID = "cham"
TCCS = machine_rules("tccs")


def sol(src: str) -> str:
    return solution_text(parse_solution(src))


def step(before: str, after: str, cls: str, rules=TCCS) -> None:
    """`after` is a one-step transition of `before` (full relation) of class `cls`."""
    outs = {}
    for state, lab in transitions(before, rules):
        outs.setdefault(state, set()).add(lab["class"])
    assert sol(after) in outs, f"{after} is not one step from {before}: {sorted(outs)}"
    assert cls in outs[sol(after)]


def derivation(states, rules=TCCS) -> None:
    for (before, _), (after, cls) in zip(states, states[1:]):
        step(before, after, cls, rules)


def fires(net, before: str, after: str) -> bool:
    """Some reaction of the network turns solution `before` into `after`."""
    start, goal = Counter(parse_solution(before)), Counter(parse_solution(after))
    start = Counter({text(m): c for m, c in start.items()})
    goal = Counter({text(m): c for m, c in goal.items()})
    for r in net.reactions:
        left = Counter(r.reactants)
        if all(start[s] >= n for s, n in left.items()) and start - left + Counter(r.products) == goal:
            return True
    return False


# sec. 2.3: {|a.0 | (~a.p | q)\b|}; the paper's p, q are the process constants P, Q
HEAVY_ION = [
    (r"{| a.0|(~a.P|Q)\b |}", None),
    (r"{| a.0, (~a.P|Q)\b |}", "heating"),                  # parallel
    (r"{| a.0, {|~a.P|Q|}\b |}", "heating"),                # restriction membrane
    (r"{| a.0, {|~a.P, Q|}\b |}", "heating"),               # parallel inside: paper line 2
    (r"{| a.0, {|~a.P <| {|Q|}|}\b |}", "heating"),         # airlock
    (r"{| a.0, {|~a.(P <| {|Q|})|}\b |}", "heating"),       # heavy ion
    (r"{| a.0, (~a.(P <| {|Q|}))\b |}", "cooling"),         # restriction membrane
    (r"{| a.0, ~a.((P <| {|Q|})\b) |}", "heating"),         # restriction ion
    (r"{| 0, (P <| {|Q|})\b |}", "reaction"),               # reaction
    (r"{| (P <| {|Q|})\b |}", "heating"),                   # inaction cleanup (not listed in the paper)
    (r"{| {|P <| {|Q|}|}\b |}", "heating"),                 # restriction membrane (printed \c: erratum)
    (r"{| {|P, Q|}\b |}", "cooling"),                       # airlock
]

# sec. 4.2: {|a.~b.0 | ((~a.0 | b.0) [] q)|}
EXTERNAL_SUM = [
    (r"{| a.~b.0|((~a.0|b.0)[]Q) |}", None),
    (r"{| a.~b.0, (~a.0|b.0)[]Q |}", "heating"),            # parallel
    (r"{| a.~b.0, <{|~a.0|b.0|}, {|Q|}> |}", "heating"),    # []-expansion
    (r"{| a.~b.0, <{|~a.0, b.0|}, {|Q|}> |}", "heating"),   # parallel inside: paper line 2
    (r"{| a.~b.0, <{|~a.0 <| {|b.0|}|}, {|Q|}> |}", "heating"),     # airlock
    (r"{| a.~b.0, <{|~a.(0 <| {|b.0|})|}, {|Q|}> |}", "heating"),   # heavy ion: line 3
    (r"{| a.~b.0, ~a.l:<0 <| {|b.0|}, {|Q|}> |}", "heating"),       # left []-ion: line 4
    (r"{| ~b.0, l:<0 <| {|b.0|}, {|Q|}> |}", "reaction"),           # line 5
    (r"{| ~b.0, 0 <| {|b.0|} |}", "cooling"),                        # left projection: line 6
    (r"{| ~b.0, 0, b.0 |}", "cooling"),                              # airlock: line 7
]

# sec. 2.2: {|a.b.0 | ~a.0 | ~b.0|}
CCS_MINUS = [
    ("{| a.b.0|~a.0|~b.0 |}", None),
    ("{| a.b.0|~a.0, ~b.0 |}", "heating"),
    ("{| a.b.0, ~a.0, ~b.0 |}", "heating"),
    ("{| b.0, 0, ~b.0 |}", "reaction"),
    ("{| 0, 0, 0 |}", "reaction"),
    ("{| 0, 0 |}", "heating"),
    ("{| 0 |}", "heating"),
    ("{||}", "heating"),
]


# --- syntax -------------------------------------------------------------------------------
def test_notation_round_trips_and_precedence():
    assert parse_molecule(r"~a.(P<|{|Q|})\b") == parse_molecule(r"~a.((P <| {|Q|})\b)")   # postfix binds tightest
    assert parse_molecule("a.b.0|~a.0|~b.0") == parse_molecule("(a.b.0|~a.0)|~b.0")
    assert parse_molecule("~a.P <| {|Q|}") == parse_molecule("(~a.P) <| {|Q|}")
    assert sol("{| Q, ~a.P |}") == sol("{|~a.P,Q|}")                                         # a multiset
    for program in ("heavy-ion-communication", "ccs-minus-execution", "nondeterministic-choice",
                    "external-sum"):
        net = generate_network(ID, program=program)
        for s in net.species:
            assert " " not in s.id and text(parse_molecule(s.id)) == s.id
            assert text(parse_molecule(s.structure)) == s.id
    with pytest.raises(ValueError, match="agents"):
        parse_molecule("a.0 | {|P|}")                     # parallel composes agents only


# --- the published derivations ---------------------------------------------------------
def test_ccs_minus_execution_sec_2_2():
    derivation(CCS_MINUS, machine_rules("ccs-minus"))
    net = generate_network(ID, program="ccs-minus-execution", machine="ccs-minus")
    assert net.status == "complete"
    for (before, _), (after, _) in zip(CCS_MINUS, CCS_MINUS[1:]):
        assert fires(net, before, after)


def test_nondeterministic_choice_sec_2_2():
    reactions = {s for s, lab in transitions("{| a.0, ~a.b.0, ~a.c.0 |}", machine_rules("ccs-minus"))
                 if lab["class"] == "reaction"}
    assert reactions == {sol("{| 0, b.0, ~a.c.0 |}"), sol("{| 0, ~a.b.0, c.0 |}")}
    step("{| 0, b.0, ~a.c.0 |}", "{| b.0, ~a.c.0 |}", "heating")     # after cleanup
    step("{| 0, ~a.b.0, c.0 |}", "{| ~a.b.0, c.0 |}", "heating")


def test_communication_through_a_membrane_sec_2_3():
    derivation(HEAVY_ION, machine_rules("ccs-minus"))
    net = generate_network(ID)                                         # the default program
    assert net.status == "complete" and len(net.species) == 18 and len(net.reactions) == 30
    for (before, _), (after, _) in zip(HEAVY_ION, HEAVY_ION[1:]):
        assert fires(net, before, after), (before, after)
    labels = dict(zip((r.to_text() for r in net.reactions), net.extras["reaction_rules"]))
    airlock = labels[r"{|P<|{|Q|}|}\b -> {|P,Q|}\b"]
    assert airlock == {"class": "cooling", "rule": "airlock", "reversible": True, "depth": 1}
    assert net.extras["compartments"][r"{|P,Q|}\b"] == [
        {"depth": 1, "context": "restriction \\b", "molecules": ["P", "Q"]}]
    assert membranes(parse_molecule(r"{|~a.P <| {|Q|}|}\b")) == [
        {"depth": 1, "context": "restriction \\b", "molecules": ["~a.P<|{|Q|}"]},
        {"depth": 2, "context": "airlock", "molecules": ["Q"]}]
    assert net.initial_state == {r"a.0|(~a.P|Q)\b": 1.0}


def test_without_airlocks_the_membrane_blocks_communication():
    # sec. 2.3: the restriction ion alone does not work when p is compound
    net = generate_network(ID, airlock=False)
    assert net.status == "complete"
    assert not any(lab["class"] == "reaction" for lab in net.extras["reaction_rules"])
    assert "0" not in {s.id for s in net.species}


def test_external_sum_sec_4_2():
    derivation(EXTERNAL_SUM)
    net = generate_network(ID, program="external-sum")
    assert net.status == "complete"
    for (before, _), (after, _) in zip(EXTERNAL_SUM, EXTERNAL_SUM[1:]):
        assert fires(net, before, after), (before, after)


# --- the laws ------------------------------------------------------------------------------
def test_heating_then_cooling_returns_the_original_solution():
    for states in (CCS_MINUS[:3], HEAVY_ION, EXTERNAL_SUM):
        for state, _ in states:
            for after, lab in transitions(state, TCCS):
                if lab["reversible"]:
                    back = {s for s, l2 in transitions(after, TCCS) if l2["reversible"]}
                    assert sol(state) in back, (state, after, lab)
    # heat {|a.b.0|~a.0|~b.0|} completely, then cool it back
    step("{| a.b.0, ~a.0, ~b.0 |}", "{| a.b.0|~a.0, ~b.0 |}", "cooling")
    step("{| a.b.0|~a.0, ~b.0 |}", "{| a.b.0|~a.0|~b.0 |}", "cooling")
    net = generate_network(ID)
    pairs = {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in net.reactions}
    for r, lab in zip(net.reactions, net.extras["reaction_rules"]):
        if lab["reversible"]:
            assert (frozenset(r.products.items()), frozenset(r.reactants.items())) in pairs


def test_airlock_law_and_membrane_law():
    # {|m|} + S <-> {|m <| S|}, for any sub-multiset S (with the chemical law)
    outs = {s for s, lab in transitions("{| a.0, b.0, c.0 |}", []) if lab["rule"] == "airlock"}
    assert sol("{| a.0 <| {|b.0, c.0|} |}") in outs and sol("{| a.0 <| {|b.0|}, c.0 |}") in outs
    assert sol("{| a.0 <| {||}, b.0, c.0 |}") in outs
    # a subsolution evolves in any context, here inside a pair inside a prefix
    step("{| c.<{|a.0, ~a.0|}, {|P|}> |}", "{| c.<{|0, 0|}, {|P|}> |}", "reaction")
    # rules do not apply inside molecules that are not solutions (only the airlock law applies)
    assert not transitions("{| c.(a.0|~a.0) |}", TCCS, airlock=False)
    assert [s for s, _ in transitions("{| c.(a.0|~a.0) |}", TCCS)] == [sol("{| c.(a.0|~a.0) <| {||} |}")]


# --- the TCCS rules (sec. 4.3) -------------------------------------------------------------
def test_tccs_rule_schemata():
    step("{| (a.P)[b/a] |}", "{| b.(P[b/a]) |}", "heating")            # relabelling ion
    step("{| ~b.(P[b/a]) |}", "{| (~a.P)[b/a] |}", "cooling")          # its inverse: preimage of ~b
    step(r"{| (c.P)\a |}", r"{| c.(P\a) |}", "heating")                 # restriction ion
    assert sol(r"{| ~a.(P\a) |}") not in {s for s, _ in transitions(r"{| (~a.P)\a |}", TCCS)}
    assert {s for s, lab in transitions("{| P (+) Q |}", TCCS) if lab["class"] == "reaction"} == {
        sol("{| P |}"), sol("{| Q |}")}                                 # internal sum
    step("{| fix_X(X=a.X) |}", "{| a.fix_X(X=a.X) |}", "heating")      # fixpoint
    step("{| fix_X(X=a.Y;Y=b.X) |}", "{| a.fix_Y(X=a.Y;Y=b.X) |}", "heating")
    step("{| <{|P|}, {|b.Q|}> |}", "{| b.r:<{|P|}, Q> |}", "heating")  # right []-ion
    step("{| r:<{|P|}, Q> |}", "{| Q |}", "cooling")                   # right projection
    step(r"{| {||}\a, {||}[b/a] |}", r"{| {||}[b/a] |}", "heating")    # cleanup
    step(r"{| {||}[b/a] |}", "{||}", "heating")


# --- parameters ------------------------------------------------------------------------------
def test_custom_rules_and_machines():
    builtin = generate_network(ID, machine="ccs-minus")
    custom = generate_network(ID, machine="custom", rules=CCS_MINUS_RULES + CCS_MINUS_CLEANUP)
    assert custom.species == builtin.species and custom.reactions == builtin.reactions
    tccs = generate_network(ID)                                     # tccs contains the CCS- rules
    assert tccs.species == builtin.species and tccs.reactions == builtin.reactions
    assert len(parse_rules(CCS_MINUS_RULES)) == 5 and len(TCCS) == 18

    own = generate_network(ID, program="custom", solution="a.0, ~a.0, ~a.0", machine="custom",
                           rules="reaction: ?a.?m, ~?a.?n -> ?m, ?n")
    assert [r.to_text() for r in own.reactions] == ["a.0 + ~a.0 -> 2 0"]
    assert own.initial_state == {"a.0": 1.0, "~a.0": 2.0}

    looping = generate_network(ID, program="custom", solution=r"fix_X(X=(a.0|X)\b)", max_species=40)
    assert looping.status == "truncated" and len(looping.species) == 40

    with pytest.raises(ValueError, match="multiset matching"):
        generate_network(ID, machine="custom", rules="bad: {|?m, ?n|} -> ?m")
    with pytest.raises(ValueError, match="only on the right"):
        generate_network(ID, machine="custom", rules="bad: ?m -> ?m, ?n")
    with pytest.raises(ValueError, match="program=custom"):
        generate_network(ID, program="custom")
    with pytest.raises(ValueError, match="machine=custom"):
        generate_network(ID, rules="reaction: ?a.?m, ~?a.?n -> ?m, ?n")
    with pytest.raises(ValueError, match="needs rules"):
        generate_network(ID, machine="custom")


def test_closure_policy_is_a_subset_of_the_full_relation():
    full, closure = Machine(TCCS, full=True), Machine(TCCS, full=False)
    for states in (HEAVY_ION, EXTERNAL_SUM):
        for state, _ in states:
            items = parse_solution(state)
            assert {s for s, _ in closure.solution_steps(items, top=False)} <= {
                s for s, _ in full.solution_steps(items, top=False)}
