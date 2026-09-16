"""SOAS reproduces the published tape-roller case studies.

Frei, Di Marzo Serugendo & Serbanuta, "Ambient intelligence in self-organising assembly systems
using the chemical reaction model", JAmI 1(3):163-184, 2010 (book [296]): rule types 1-4
(sec. 5.1-5.4) and the simulation trace of the four-task case study (sec. 6.5, figs. 22-26).
Frei, Serbanuta & Di Marzo Serugendo, "Self-organising assembly systems formally specified in
Maude", JAmI 5:491-510 (online 2012): the six-task case study (sec. 3), the modules (sec. 5.1),
the layout rules (sec. 5.3) and the execution result (sec. 5.5, appendices A and B). Book 20.1.

Counts marked "soas-maude" come from running the authors' own specification (soas-maude 1.0,
sha1 90d77bf...) under Maude 3.5.1, which is the published form of the rules.
"""

import itertools
from collections import Counter

import pytest

from chemart import generate_network
from chemart.chemistries.soas import (
    GAP1, GAP2, GENERIC_HUMAN, PATTERNS, TAPE_ROLLER_MRAS, TAPE_ROLLER_PARTS, Skill, decompose,
    parse_skill, skill_match,
)

ID = "soas"

# 2012 appendix B, first solution: the coalitions and the order they sit in.
APPENDIX_B = {
    "t(1)": GENERIC_HUMAN, "t(2)": "f1+g1+r1", "t(3)": "f2+g2+r2", "t(4)": "f3+g3+r3",
    "t(5)": "a1+a4+a7+f4+g4", "t(6)": GENERIC_HUMAN,
}
# 2012 sec. 5.5: the 25th solution moves r3 to task 2, r1 to task 3 and r2 to task 4.
SOLUTION_25 = {**APPENDIX_B, "t(2)": "f1+g1+r3", "t(3)": "f2+g2+r1", "t(4)": "f3+g3+r2"}


def analysis(net):
    return net.extras["analysis"]


def lines(net):
    return [line["assignment"] for line in analysis(net)["assembly_lines"]]


@pytest.fixture(scope="module")
def gap2():
    return generate_network(ID)


@pytest.fixture(scope="module")
def gap1():
    return generate_network(ID, system="tape-roller-gap1")


# --- the published runs -----------------------------------------------------------------
def test_forty_eight_assembly_lines_for_the_four_part_tape_roller(gap2):
    # 2012 sec. 5.5 prints both 47 and 48; the specification itself gives 48 (soas-maude).
    assert analysis(gap2)["n_assembly_lines"] == 48
    assert len(lines(gap2)) == 48 == len({tuple(sorted(x.items())) for x in lines(gap2)})
    assert APPENDIX_B in lines(gap2)
    assert SOLUTION_25 in lines(gap2)

    # every line is one of the 24 permutations of the four workers over the four pick&place
    # tasks, times the two ways of splitting the body-case feeders f1 and f3
    workers = ("r1", "r2", "r3", "a1+a4+a7")
    expected = set()
    for order in itertools.permutations(workers):
        for body1, body3 in (("f1", "f3"), ("f3", "f1")):
            expected.add(tuple(
                "+".join(sorted(w.split("+") + [feeder, gripper]))
                for w, feeder, gripper in zip(
                    order, (body1, "f2", body3, "f4"), ("g1", "g2", "g3", "g4"))))
    assert {tuple(x[f"t({i})"] for i in (2, 3, 4, 5)) for x in lines(gap2)} == expected


def test_coalition_formation_matches_the_maude_specification(gap2):
    # soas-maude: 165 coalition molecules, 30 of them complete
    assert analysis(gap2)["coalitions"] == {"t(1)": 1, "t(2)": 44, "t(3)": 44, "t(4)": 44,
                                            "t(5)": 31, "t(6)": 1}
    complete = analysis(gap2)["complete_coalitions"]
    assert {task: len(v) for task, v in complete.items()} == {
        "t(1)": 1, "t(2)": 8, "t(3)": 8, "t(4)": 8, "t(5)": 4, "t(6)": 1}
    assert set(complete["t(5)"]) == {"f4+g4+r1", "f4+g4+r2", "f4+g4+r3", "a1+a4+a7+f4+g4"}
    assert complete["t(1)"] == complete["t(6)"] == [GENERIC_HUMAN]
    # a coalition molecule does not depend on the order in which its members joined
    assert analysis(gap2)["coalition_conflicts"] == 0
    assert analysis(gap2)["unassignable_tasks"] == []


def test_layout_of_a_line_matches_appendix_b(gap2):
    line = next(x for x in analysis(gap2)["assembly_lines"] if x["assignment"] == APPENDIX_B)

    def station(task, kind, subtype, position, angle, coalition=None):
        out = {"task": task, "kind": kind, "subtype": subtype, "position": position,
               "angle": [angle, 0, 0]}
        if coalition is not None:
            out["coalition"] = coalition
        return out

    # 2012 appendix B: the loading operator of t(1) is not placed; the line runs east, turns
    # north at the end of the shop floor and comes back west (t(6) from the soas-maude run)
    assert line["layout"] == [
        station("t(2)", "conveyor", "linear", [0, 0, 0], 0),
        station("t(2)", "coalition", None, [200, 50, 0], 0, "f1+g1+r1"),
        station("t(3)", "conveyor", "linear", [300, 0, 0], 0),
        station("t(3)", "coalition", None, [500, 50, 0], 0, "f2+g2+r2"),
        station("t(4)", "conveyor", "linear", [600, 0, 0], 0),
        station("t(4)", "coalition", None, [800, 50, 0], 0, "f3+g3+r3"),
        station("t(17/4)", "conveyor", "corner", [900, 0, 0], 0),
        station("t(9/2)", "conveyor", "linear", [900, 50, 0], -90),
        station("t(19/4)", "conveyor", "corner", [900, 350, 0], -90),
        station("t(5)", "conveyor", "linear", [600, 350, 0], 180),
        station("t(5)", "coalition", None, [600, 300, 0], 180, "a1+a4+a7+f4+g4"),
        station("t(6)", "conveyor", "linear", [300, 350, 0], 180),
        station("t(6)", "coalition", None, [300, 300, 0], 180, GENERIC_HUMAN),
    ]


def test_2010_trace_of_the_two_part_tape_roller(gap1):
    # soas-maude on the four-task order: 46 coalitions, 12 lines
    assert sum(analysis(gap1)["coalitions"].values()) == 46
    assert analysis(gap1)["n_assembly_lines"] == 12
    # 2010 sec. 6.5: "t1 by human1, t2 by coalition r3 - g1, t3 by r1 - g2, t4 by human2"
    trace = {"t(1)": GENERIC_HUMAN, "t(2)": "g1+r3", "t(3)": "g2+r1", "t(4)": GENERIC_HUMAN}
    assert trace in lines(gap1)

    # figs. 22-26: g1 starts the coalition for t(2) and r3 joins it; r1 starts the one for
    # t(3) and g2 joins it; both are complete, with the open interfaces of figs. 25 and 26
    structures = {s.id: s.structure for s in gap1.species}
    assert "open-interfaces=circular+" in structures["gap1->t(2)[g1]"]
    for done in ("gap1->t(2)[g1,r3]", "gap1->t(3)[g2,r1]"):
        assert "required-ops=none" in structures[done]
        assert "open-interfaces=diamond-, straight-, triangular-" in structures[done]
    assert "open-interfaces=circular-, diamond-, straight-, triangular-" in structures["gap1->t(3)[r1]"]
    def fired(reactants, products):
        want = (dict(Counter(reactants)), dict(Counter(products)))
        return any((r.reactants, r.products) == want for r in gap1.reactions)

    assert fired(["g1", "gap1"], ["g1", "gap1", "gap1->t(2)[g1]"])
    assert fired(["r3", "gap1->t(2)[g1]"], ["r3", "gap1->t(2)[g1]", "gap1->t(2)[g1,r3]"])
    assert fired(["r1", "gap1"], ["r1", "gap1", "gap1->t(3)[r1]"])
    assert fired(["g2", "gap1->t(3)[r1]"], ["g2", "gap1->t(3)[r1]", "gap1->t(3)[g2,r1]"])


# --- lines that must not appear ---------------------------------------------------------
def test_invalid_lines_are_not_produced(gap2):
    # p1 (body case) is gripped 60 mm apart, p2 (tape roll) 40 mm, p3 needs a vacuum gripper
    # and p4 a screw driver; each feeder feeds one part type (2012 sec. 5.1, parts of the GAP)
    grippers = {"t(2)": {"g1"}, "t(3)": {"g1", "g2"}, "t(4)": {"g3"}, "t(5)": {"g4"}}
    feeders = {"t(2)": {"f1", "f3"}, "t(3)": {"f2"}, "t(4)": {"f1", "f3"}, "t(5)": {"f4"}}
    complete = analysis(gap2)["complete_coalitions"]
    for task, allowed in grippers.items():
        used = {m for members in complete[task] for m in members.split("+")}
        assert used & {"g1", "g2", "g3", "g4"} <= allowed, task
        assert used & {"f1", "f2", "f3", "f4"} <= feeders[task], task
    # the 55 mm gripper can never take the 60 mm body case, and no feeder of the wrong part
    assert not any("g2" in c for c in complete["t(2)"] + complete["t(4)"])
    assert not any("f2" in c for c in complete["t(2)"] + complete["t(4)"])
    # positioning devices offer no skill any task needs, so they never join a coalition
    assert not any(f"pd{i}" in c for task in complete.values() for c in task for i in (1, 2, 3, 4))

    for line in lines(gap2):
        assert set(line) == {f"t({i})" for i in range(1, 7)}
        assert line["t(1)"] == line["t(6)"] == GENERIC_HUMAN
        used = [m for task, members in line.items() if members != GENERIC_HUMAN
                for m in members.split("+")]
        assert len(used) == len(set(used)), f"a module works at two stations: {line}"
        for task, members in line.items():
            assert members in complete[task], f"{members} does not complete {task}"


def test_missing_modules_leave_tasks_unassigned(gap1):
    # 2010 sec. 7 (Wermelinger's first property): the cham also terminates when the required
    # skills are not provided in sufficient number, and one or several tasks remain open
    without_gripper = generate_network(ID, unavailable=["g4"])
    assert without_gripper.status == "complete"
    assert analysis(without_gripper)["unassignable_tasks"] == ["t(5)"]
    assert analysis(without_gripper)["n_assembly_lines"] == 0

    # four pick&place tasks need four workers; with r1 gone only three are left, so every task
    # still has coalitions but no line covers the order
    without_robot = generate_network(ID, unavailable=["r1"])
    assert analysis(without_robot)["unassignable_tasks"] == []
    assert analysis(without_robot)["n_assembly_lines"] == 0
    # the two-task order survives the same failure (resilience, 2010 sec. 8.1)
    assert analysis(generate_network(ID, system="tape-roller-gap1", unavailable=["r1"]))[
        "n_assembly_lines"] == 6
    assert analysis(gap1)["n_assembly_lines"] == 12


def test_composition_patterns_gate_the_coalitions(gap1):
    # rule type 2: without "feed & move" a feeder never joins a robot, so no coalition can
    # cover a task that also requires feeding (the gap2 tasks start at the feeder)
    without = "\n".join(line for line in PATTERNS.splitlines() if line != "feed & move")
    crippled = generate_network(ID, patterns=without)
    assert analysis(crippled)["n_assembly_lines"] == 0
    assert analysis(crippled)["unassignable_tasks"] == ["t(2)", "t(3)", "t(4)", "t(5)"]
    # the gap1 tasks name their feeder instead, so they need no feed skill and are unaffected
    same = generate_network(ID, system="tape-roller-gap1", patterns=without)
    assert analysis(same)["n_assembly_lines"] == analysis(gap1)["n_assembly_lines"] == 12


# --- notation and matching ---------------------------------------------------------------
def test_skill_matching_and_composite_skills(gap2):
    wide = parse_skill("grip(gripper-type=2finger,range=70.0)")
    narrow = parse_skill("grip(gripper-type=2finger,range=55.0)")
    needed = parse_skill("grip(gripper-type=2finger,range=60.0)")
    assert skill_match(wide, needed) and not skill_match(narrow, needed)
    assert not skill_match(parse_skill("grip(gripper-type=vacuum,range=70.0)"), needed)
    assert skill_match(parse_skill("move(subtype=linear,direction=vertical,range=125.0)"),
                       parse_skill("move(direction=vertical)"))
    # rule type 3: a gripper and an axis together can pick and place
    assert decompose((parse_skill("pick&place"),)) == (
        Skill("move", ()), Skill("move", (("direction", "vertical"),)))

    # the parts complete the order and the grippers (2012 extension of rule 4): the body case
    # p1 is gripped 60 mm apart, the tape roll p2 40 mm, and g1 opens 70 mm
    structures = {s.id: s.structure for s in gap2.species}
    assert "grip(gripper-type=2finger,range=60.0)" in structures["gap2"]
    assert "grip(gripper-type=2finger,range=40.0)" in structures["gap2"]
    assert "feed(subtype=feeds(body-case))" in structures["gap2"]
    assert "grip(gripper-type=2finger,range=70.0)" in structures["g1"]
    assert "open-close" not in structures["g1"]
    assert "grip(gripper-type=vacuum)" in structures["g3"]


def test_custom_system_and_parameter_errors(gap2):
    custom = generate_network(ID, system="custom", mras=TAPE_ROLLER_MRAS, gap=GAP2,
                              parts=TAPE_ROLLER_PARTS)
    assert [s.id for s in custom.species] == [s.id for s in gap2.species]
    assert [r.to_text() for r in custom.reactions] == [r.to_text() for r in gap2.reactions]

    with pytest.raises(ValueError, match="only used with system=custom"):
        generate_network(ID, mras=TAPE_ROLLER_MRAS)
    with pytest.raises(ValueError, match="needs mras and gap"):
        generate_network(ID, system="custom")
    with pytest.raises(ValueError, match="unknown MRAs"):
        generate_network(ID, unavailable=["r9"])
    with pytest.raises(ValueError, match="4 fields"):
        generate_network(ID, system="custom", mras="r1 robot | move", gap=GAP2)
    with pytest.raises(ValueError, match="items are key=value"):
        generate_network(ID, system="custom", mras="r1 robot | move(speed=2) | | circular-",
                         gap=GAP2)
    with pytest.raises(ValueError, match="5 fields"):
        generate_network(ID, system="custom", mras=TAPE_ROLLER_MRAS, gap="gap0 1\nt(1) | other")
    with pytest.raises(ValueError, match="floor"):
        generate_network(ID, floor=[1000])


def test_network_shape(gap2):
    assert gap2.status == "complete"
    assert (len(gap2.species), len(gap2.reactions)) == (582, 1326)
    assert gap2.provides == ["catalysts", "initial-state", "stoichiometry", "topology"]
    assert all(r.rate is None for r in gap2.reactions)          # the rules have no kinetics
    assert all(" " not in s.id and s.structure for s in gap2.species)
    # the initial solution: the order and one of each module
    assert gap2.initial_state["gap2"] == 1.0
    assert set(gap2.initial_state) - {"gap2"} == {
        "r1", "r2", "r3", "a1", "a4", "a7", "f1", "f2", "f3", "f4", "g1", "g2", "g3", "g4",
        "pd1", "pd2", "pd3", "pd4", GENERIC_HUMAN}
    # modules and coalitions are catalysts of initiate and join; assignment consumes the line
    kinds = [label["rule"] for label in gap2.extras["reaction_rules"]]
    assert len(kinds) == len(gap2.reactions)
    assert set(kinds) == {"initiate", "join", "assign-generic", "assign"}
    for reaction, rule in zip(gap2.reactions, gap2.extras["reaction_rules"]):
        if rule["rule"] in ("initiate", "join"):
            assert len(reaction.catalysts) == 2 and len(reaction.products) == 3
        else:
            assert not reaction.catalysts and len(reaction.products) == 1
    assert analysis(gap2)["dead_end_lines"] == 0

    small = generate_network(ID, max_species=300)
    assert small.status == "truncated" and len(small.species) == 300
