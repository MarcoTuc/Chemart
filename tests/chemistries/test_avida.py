"""Avida: book 10.7.1; Ofria, Bryson & Wilke 2009 (book [629]); Lenski et al., Nature 423:139 (2003).

The machine is a port of the Avida 2.14 source (devosoft/avida). Reference values
marked "upstream" were produced by the compiled original (analyze mode: LOAD_SEQUENCE,
RECALCULATE, DETAIL; and population runs with support/config defaults), see the
catalog decisions.
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import avida as A

ID = "avida"
ANCESTOR = A.DEFAULT_ANCESTOR
NOT_GENOME = "wzcagc" + "yopcuyk" + "c" * 78 + "zvfcaxgab"       # IO push pop-C nand IO: outputs ~input
NO_MUTATIONS = dict(copy_mut_prob=0.0, divide_ins_prob=0.0, divide_del_prob=0.0)


def cpu(genome, cycles, **world):
    """Run a lone organism for a number of CPU cycles on the test inputs."""
    w = A.World(np.random.default_rng(0), width=3, height=3, record=False, age_limit=0, **world)
    org = A.Organism(genome, 0)
    org.cell = 0
    for _ in range(cycles):
        w.process(org)
    return org, w


# --- the ancestor ----------------------------------------------------------------
def test_default_ancestor_self_replicates_exactly():
    # default-heads.org; upstream tests/heads_default_100u: gestation 389, merit 97 (copied 100, executed 97)
    assert len(ANCESTOR) == 100 and ANCESTOR.count("c") == 88
    r = A.test_cpu(ANCESTOR)
    assert r["divided"] and r["offspring"] == ANCESTOR
    assert (r["gestation"], r["executed_size"], r["copied_size"], r["merit"]) == (389, 97, 100, 97.0)
    assert r["tasks"] == []


# instruction semantics of the heads CPU --------------------------------------------
def test_nop_modification_example_of_ofria_2009_figure_3():
    # pop; pop nop-A; inc; inc nop-A; inc nop-A; swap nop-C; add nop-A -> BX = 1, CX = 2, AX = 3
    org, _ = cpu("ppammamarcsa" + "q" * 10, 7)
    assert org.reg == [3, 1, 2]
    assert org.heads[A.HEAD_IP] == 12


def test_template_jump_example_of_ofria_2009_figure_2():
    # h-search nop-A nop-B puts the flow head after the complement nop-B nop-C; mov-head jumps there
    genome = "zabgpqqbcm" + "q" * 10
    org, w = cpu(genome, 2)
    assert org.heads[A.HEAD_FLOW] == 9 and org.heads[A.HEAD_IP] == 9
    assert org.reg[A.BX] == 6 and org.reg[A.CX] == 2           # distance to the label, label size
    w.process(org)                                               # executes inc at 9: pop at 4 was skipped
    assert org.reg[A.BX] == 7


def test_register_stack_and_head_instructions():
    org, _ = cpu("dmq" + "q" * 10, 1)                           # if-n-equ with BX == CX skips the next line
    assert org.heads[A.HEAD_IP] == 2
    org, _ = cpu("mdmq" + "q" * 10, 3)                          # BX = 1 != CX: the next inc runs
    assert org.reg[A.BX] == 2
    org, _ = cpu("uq" + "q" * 10, 1)                            # nand of BX = CX = 0
    assert org.reg[A.BX] == -1
    org, _ = cpu("ukq" + "q" * 10, 2)                           # shift-r is arithmetic
    assert org.reg[A.BX] == -1
    org, _ = cpu("m" + "l" * 31 + "qqqq", 32)                   # 32-bit registers wrap
    assert org.reg[A.BX] == -(2 ** 31)
    org, _ = cpu("mopp" + "q" * 10, 3)                          # push 1, pop 1 ...
    assert org.reg[A.BX] == 1
    org, _ = cpu("mopp" + "q" * 10, 4)                          # ... then a pop of the empty stack gives 0
    assert org.reg[A.BX] == 0
    org, _ = cpu("moqp" + "q" * 10, 4)                          # swap-stk: the other stack is empty
    assert org.reg[A.BX] == 0
    org, _ = cpu("moqpqp" + "q" * 10, 6)                        # and the first one still holds the value
    assert org.reg[A.BX] == 1
    # inc nop-C x3 (CX = 3); set-flow (flow = CX); jmp-head nop-C (write += CX); get-head nop-A (CX = IP)
    org, _ = cpu("mcmcmcjhciaq" + "q" * 10, 6)
    assert org.heads[A.HEAD_FLOW] == 3 and org.heads[A.HEAD_WRITE] == 3 and org.reg[A.CX] == 10
    org, _ = cpu("ww" + "q" * 20, 2)                             # h-alloc once: 22 -> 66, AX = 22
    assert len(org.mem) == 66 and org.reg[A.AX] == 22 and org.mal_active


def test_divide_needs_half_the_offspring_copied_and_half_the_parent_executed():
    # Ofria et al. 2009, 3.1.2: h-divide fails if less than half of the parent was executed
    # or less than half of the offspring's memory was copied into
    w = A.World(np.random.default_rng(0), width=3, height=3, **NO_MUTATIONS)
    org = w.inject(ANCESTOR, 4)
    w.process(org)                                               # h-alloc: 100 -> 300 lines
    assert len(org.mem) == 300
    org.heads[A.HEAD_READ], org.heads[A.HEAD_WRITE] = 100, 200   # offspring = lines 100..199
    assert not w.divide(org)
    org.copied[100:150] = bytes([1]) * 50
    assert not w.divide(org)                                     # only line 0 was executed
    org.executed[0:50] = bytes([1]) * 50
    assert w.divide(org) and w.stats["births"] == 1
    # without h-copy the ancestor never divides
    assert not A.test_cpu(ANCESTOR.replace("zvf", "zcf"))["divided"]


@pytest.mark.parametrize("genome, gestation, merit, executed, tasks", [
    (NOT_GENOME, 387, 194.0, 97, ["NOT"]),
    ("wzcagcyopcuyyopcuykcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab", 386, 194.0, 97, ["NOT"]),
    ("wzcagcccccccccccccccccccccccccccccceccccccccccccfcccccccccccccjccccccccccccpccccccccccccccczvfcaxgab", 375, 87.0, 87, []),
    ("wzcagczccccmicccccccccccccccceccccvqeacvbczcccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab", 364, 83.0, 83, []),
    ("wzcagccccccccccccykslybjybyrujyuylbcueciyckyccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab", 382, 192.0, 96, ["NAND"]),
    ("wzcagcccccccccccccccccccccccccccccccconqhayumnsskyreumyruayooyuytcccccccccccccccccccccccccczvfcaxgab", 384, 190.0, 95, ["NAND"]),
    ("wzcagcnmycdeuejtysnqyybooylmhmemocccccccccccccccccccccccccccccccccccccccccccccccccccccccccczvfcaxgab", 383, 188.0, 94, ["NOT"]),
    ("wzcagcccccccccccycccccccccccccccccccccccccccccccccchsoojqyctkbmjjmyjcuqyaeudyyccccccccccccczvfcaxgab", 380, 376.0, 94, ["ORN"]),
])
def test_test_cpu_matches_upstream_avida(genome, gestation, merit, executed, tasks):
    r = A.test_cpu(genome)
    assert (r["gestation"], r["merit"], r["executed_size"], r["tasks"]) == (gestation, merit, executed, tasks)


# --- merit and the environment ------------------------------------------------------
def test_logic_tasks_and_rewards_of_lenski_2003_table_1():
    a, b = A.TEST_INPUTS[1], A.TEST_INPUTS[0]              # most recent input first
    outputs = {
        "NOT": ~a, "NAND": ~(a & b), "AND": a & b, "ORN": a | ~b, "OR": a | b,
        "ANDN": a & ~b, "NOR": ~(a | b), "XOR": a ^ b, "EQU": ~(a ^ b),
    }
    table_1 = {"NOT": 2, "NAND": 2, "AND": 4, "ORN": 4, "OR": 8, "ANDN": 8, "NOR": 16, "XOR": 16, "EQU": 32}
    w = A.World(np.random.default_rng(0), width=3, height=3, record=False)
    for task, out in outputs.items():
        assert A.logic_id([a, b], out) in A.TASKS[task]
        org = A.Organism(ANCESTOR, 0)
        org.input_buf = [a, b]
        w.do_output(org, out)
        assert org.cur_bonus == table_1[task]
        w.do_output(org, out)                                # rewarded once per gestation (max_count=1)
        assert org.cur_bonus == table_1[task]
    assert A.logic_id([a, b], a) not in {i for ids in A.TASKS.values() for i in ids}   # ECHO is not rewarded
    # bits 0 and 1 of the input are both 0 but the output bits differ: no truth table
    assert A.logic_id([0x0F000000], 1) == -1


def test_not_genome_gets_the_published_merit_bonus():
    r = A.test_cpu(NOT_GENOME)
    assert r["tasks"] == ["NOT"] and r["merit"] == 2 * 97 and r["offspring"] == NOT_GENOME
    assert A.test_cpu(NOT_GENOME, rewards={})["merit"] == 97
    assert A.test_cpu(NOT_GENOME, rewards={"NOT": 3.0})["merit"] == 8 * 97


def test_merit_is_set_at_divide_and_inherited():
    # Ofria et al. 2009, 3.2.4: the bonus of a gestation sets the merit of parent and offspring at divide
    w = A.World(np.random.default_rng(1), width=5, height=5, **NO_MUTATIONS)
    parent = w.inject(NOT_GENOME, 12)
    assert parent.merit == 100.0                                # injected: merit = genome length
    while w.stats["births"] == 0:
        w.process(parent)
    child = next(o for o in w.cells if o is not None and o is not parent)
    assert parent.merit == child.merit == 194.0
    assert parent.cur_bonus == 1.0 and parent.last_tasks == {"NOT": 1}


def test_cpu_time_is_proportional_to_merit():
    # Ofria et al. 2009, 3.2.1: twice the merit, twice the instructions executed
    w = A.World(np.random.default_rng(2), width=4, height=4, record=False)
    w.set_merit(3, 97.0)
    w.set_merit(10, 194.0)
    counts = Counter(w.next_cell() for _ in range(30000))
    assert set(counts) == {3, 10}
    assert counts[10] / counts[3] == pytest.approx(2.0, rel=0.05)


# --- the observed network -------------------------------------------------------------
def classify(r):
    n_in, n_out = sum(r.reactants.values()), sum(r.products.values())
    return {(1, 0): "age-death", (1, 1): "parent-replaced", (1, 2): "birth", (2, 2): "overwrite"}[(n_in, n_out)]


def test_default_run_contains_replication_and_overwrite_events():
    net = generate_network(ID, seed=3)
    assert net.status == "observed"
    struct = {s.id: s.structure for s in net.species}
    assert all(A.genotype_id(g) == sid for sid, g in struct.items())
    kinds = Counter()
    for r in net.reactions:
        kind = classify(r)
        kinds[kind] += r.count
        if kind in ("birth", "overwrite"):
            assert r.catalysts                                  # the parent survives its own birth event
    a = net.extras["analysis"]
    assert kinds["birth"] > 0 and kinds["overwrite"] > 0
    assert kinds["overwrite"] == a["overwrites"] and kinds["age-death"] == a["age_deaths"]
    assert kinds["birth"] + kinds["overwrite"] + kinds["parent-replaced"] == a["births"]
    assert sum(a["per_update"]["births"]) == a["births"]
    ancestor = A.genotype_id(ANCESTOR)
    assert net.initial_state == {ancestor: 1.0} and struct[ancestor] == ANCESTOR
    grid = net.extras["space"]["final_grid"]
    assert len(grid) == 144 and sum(net.extras["final_state"].values()) == sum(g is not None for g in grid)
    # after the lattice fills, births overwrite: space is the limiting resource
    assert a["per_update"]["organisms"][-1] > 100


def test_without_mutations_the_network_is_pure_replication():
    net = generate_network(ID, seed=1, updates=120, **NO_MUTATIONS)
    anc = A.genotype_id(ANCESTOR)
    assert [s.id for s in net.species] == [anc]
    allowed = {"birth": ({anc: 1}, {anc: 2}), "overwrite": ({anc: 2}, {anc: 2}),
               "parent-replaced": ({anc: 1}, {anc: 1}), "age-death": ({anc: 1}, {})}
    for r in net.reactions:
        assert (r.reactants, r.products) == allowed[classify(r)]


def test_task_completions_are_recorded_in_extras():
    net = generate_network(ID, seed=4, ancestor=NOT_GENOME, updates=60, world_x=5, world_y=5, **NO_MUTATIONS)
    gid = A.genotype_id(NOT_GENOME)
    assert net.extras["task_events"][gid]["NOT"] > 0
    assert all(s.id == gid for s in net.species)
    a = net.extras["analysis"]
    assert a["tasks_performed"]["NOT"] == net.extras["task_events"][gid]["NOT"]
    # every organism has divided (or inherited the merit of a parent that did): 97 x 2 for the NOT bonus
    assert a["per_update"]["ave_merit"][-1] == 194.0
    quiet = generate_network(ID, seed=4, ancestor=NOT_GENOME, updates=60, world_x=5, world_y=5, rewards={},
                             **NO_MUTATIONS)
    assert quiet.extras["task_events"] == {}
    assert quiet.extras["analysis"]["per_update"]["ave_merit"][-1] == 97.0


def test_reproducible_with_a_seed():
    kw = dict(updates=60, world_x=6, world_y=6)
    assert generate_network(ID, seed=11, **kw).to_dict() == generate_network(ID, seed=11, **kw).to_dict()
    assert generate_network(ID, seed=11, **kw).to_dict() != generate_network(ID, seed=12, **kw).to_dict()


def test_mass_action_birth_method_places_offspring_anywhere():
    w = A.World(np.random.default_rng(5), width=10, height=10, birth_method="mass-action", **NO_MUTATIONS)
    w.inject(ANCESTOR, 0)
    w.run(40)
    cells = [c for c, o in enumerate(w.cells) if o is not None]
    near = set(w.neighbors[0]) | {0}
    assert len(cells) > 1 and any(c not in near for c in cells)


@pytest.mark.slow
def test_population_growth_matches_upstream_avida():
    # upstream Avida, support/config defaults on the 60x60 torus: 64.5 organisms at update 100
    # (6 seeds, 59-73; the shipped heads_default_100u test gives 53) and 386 +- 29 at update 200 (18 seeds)
    at100, at200 = [], []
    for seed in range(6):
        w = A.World(np.random.default_rng(seed), width=60, height=60, record=False)
        w.inject(ANCESTOR, 0)
        w.run(100)
        at100.append(w.num_orgs)
        w.run(100)
        at200.append(w.num_orgs)
    assert 50 < np.mean(at100) < 80
    assert 345 < np.mean(at200) < 430


def test_parameters_are_checked():
    with pytest.raises(ValueError, match="letters"):
        generate_network(ID, ancestor="WZC")
    with pytest.raises(ValueError, match="length"):
        generate_network(ID, ancestor="wzcag")
    with pytest.raises(ValueError, match="unknown tasks"):
        generate_network(ID, rewards={"FOO": 1.0})
    with pytest.raises(ValueError, match="number"):
        generate_network(ID, rewards={"NOT": "x"})
