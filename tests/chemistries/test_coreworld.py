"""Coreworld: Rasmussen, Feldberg, Hindsholm & Knudsen, LA-UR-89-2618 = Physica D 42:111 (1990).

Book 10.6.3, ref [697]. Reproduces the VENUS simulator of section 2: Table 1
(instructions and addressing modes), eqs. 1-3 (computational resources),
Table 2 (desert and jungle), fig. 3 (a self-propagating MOV $0,$1 and a copy
loop), the appendix program MICE (8 instructions, a cycle of 18 core
updates), the pointer fixed points of section 3, and the observations that
jungles are more active than deserts and that engineered programs are too
brittle to survive the noise.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import coreworld as C

ID = "coreworld"
BLANK = "DAT #0, #0"


def machine(lines, M=3584, pointers=(0,), fill=BLANK, key=True, seed=0, **kw):
    core = [C.parse(fill, M)] * M
    for k, line in enumerate(lines):
        core[k] = C.parse(line, M)
    kw.setdefault("resources", False)
    v = C.Venus(core, C._buffered(np.random.default_rng(seed)), **kw)
    if key:
        v.key = lambda w: w[0]
    v.queue = list(pointers)
    return v


def listing(v, n):
    return [C.text(w, v.M) for w in v.core[:n]]


def named_events(v):
    return {(tuple(C.OPS[k] for k in lhs), tuple(C.OPS[k] for k in rhs)): n
            for (lhs, rhs), (n,) in v.events.items()}


# --- Table 1: words ------------------------------------------------------------
def test_table_1_words_and_appendix_listing():
    M = 3584
    assert C.OPS == ("DAT", "JMP", "JMZ", "JMN", "DJN", "ADD", "SUB", "MOV", "CMP", "SPL")
    assert C.MODES == "#$@<"
    mice = [C.text(C.parse(s, M), M) for s in C.MICE]
    assert mice == ["DAT #0, #7", "MOV #7, $-1", "MOV @-2, <5", "DJN $-1, $-3", "SPL #0, @3",
                    "ADD #417, $2", "JMZ $-5, $-6", "DAT #0, #714"]
    assert C.parse("JMP $-2", M) == C.parse("JMP $-2, $0", M)
    assert C.word_id(C.parse("MOV $0, $1", M), M) == "MOV_$0_$1"
    with pytest.raises(ValueError, match="opcode"):
        C.parse("NOP $1", M)


# --- fig. 3 ---------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["venus-i", "venus-ii"])
def test_fig_3a_self_propagating_mov(mode):
    v = machine(["MOV $0, $1"], mode=mode)
    for t in range(3):
        v.update()
    assert listing(v, 5) == ["MOV $0, $1"] * 4 + [BLANK]
    assert v.queue == [3]
    # the network of this run: the executing MOV copies itself over a DAT
    assert named_events(v) == {(("DAT", "MOV"), ("MOV", "MOV")): 3}


@pytest.mark.parametrize("mode", ["venus-i", "venus-ii"])
def test_fig_3b_copy_loop_multiplies_an_instruction(mode):
    loop = ["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2", "DAT #2"]
    v = machine(loop, pointers=[1], mode=mode)
    frames = []
    for _ in range(4):
        frames.append((v.queue[0], listing(v, 8)))
        v.update()
    b = BLANK
    assert frames == [
        (1, ["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2, $0", "DAT #0, #2", b, b, b]),
        (2, ["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2, $0", "DAT #0, #2", b, "MOV #4, $10", b]),
        (3, ["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2, $0", "DAT #0, #3", b, "MOV #4, $10", b]),
        (1, ["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2, $0", "DAT #0, #3", b, "MOV #4, $10", b]),
    ]
    # "very powerful in multiplying any instruction": one copy per 3-update cycle
    for _ in range(3 * 100 - 4):
        v.update()
    copies = [i for i, w in enumerate(v.core) if w == C.parse("MOV #4, $10", v.M)]
    assert len(copies) == 1 + 100 and copies[1:] == list(range(6, 106))
    events = named_events(v)
    assert events[(("DAT", "MOV", "MOV"), ("MOV", "MOV", "MOV"))] == 100


# --- appendix: MICE ------------------------------------------------------------
@pytest.mark.parametrize("mode", ["venus-i", "venus-ii"])
def test_mice_replicates_with_a_cycle_of_18_core_updates(mode):
    M = 3584
    v = machine(C.MICE, pointers=[1], mode=mode)
    for t in range(1, 19):
        v.update()
        if t < 16:
            assert len(v.queue) == 1                    # the SPL is the 16th instruction of the cycle
    assert v.queue[0] == 1                              # back on its first MOV after 18 updates
    parent = [C.parse(s, M) for s in C.MICE]
    child = v.core[713:721]
    assert child[1:7] == parent[1:7]                    # the 6 executable words copied 713 cells away
    assert child[0] == C.parse("DAT #7", M)             # the copy's counter, set by its own MOV #7
    assert v.queue[1] == 716                            # the SPL started a pointer on the copy's MOV (714)
    assert v.core[7] == C.parse("DAT #1124", M)         # 714 - 7 copies + 417 for the next offspring
    # the replication keeps going: pointers and copies multiply
    for _ in range(5 * 18):
        v.update()
    assert len(v.queue) >= 8
    body = parent[1:7]
    intact = sum(all(v.core[(i + k) % M] == body[k] for k in range(6)) for i in range(M))
    assert intact >= 6


# --- section 3: pointer fixed points --------------------------------------------
@pytest.mark.parametrize("word", ["JMP #17", "JMZ #5, $0", "JMZ #5, #0", "JMN $0, #3", "DJN #9, $1"])
def test_fixed_points_trap_pointers(word):
    # "JMP(#X), JMZ(#X, 0), JMN($0, Y), and DJN(#X, Y) ... are all characterized by pointing to themselves"
    v = machine([word, "DAT #3000"])
    assert v.self_loop(0)
    for _ in range(50):
        v.update()
    assert v.queue == [0]
    assert not machine(["JMP #0, $0", "MOV $0, $1"], pointers=[1]).self_loop(1)


def test_dat_kills_spl_splits_and_the_queue_is_bounded():
    v = machine([], M=64, fill="SPL $0, $0", queue_len=10)
    for _ in range(12):
        v.update()
    assert len(v.queue) == 10                           # at most L pointers survive
    assert v.stats["splits_refused"] > 0
    w = machine(["DAT #0"], queue_len=10)
    w.update()
    assert w.queue == [] and w.stats["deaths"] == 1


def test_cmp_skips_if_unequal_and_spl_targets_b():
    # Table 1: "compare A and B and skip next statement if unequal"; "split execution between B and next"
    v = machine(["CMP #3, $1", "DAT #3"])
    assert v.execute(0, {}) == ((1,), None)
    v = machine(["CMP #4, $1", "DAT #3"])
    assert v.execute(0, {}) == ((2,), None)
    v = machine(["SPL $7, $5"])
    assert v.execute(0, {}) == ((1, 5), None)


# --- section 2: parallel update, operation radius, noise -----------------------
def test_venus_i_instructions_see_the_same_core_and_later_pointer_wins():
    prog = ["ADD #1, $2", "MOV $1, $2", "DAT #5"]
    par, seq = machine(prog, pointers=[0, 1]), machine(prog, pointers=[0, 1], mode="venus-ii")
    par.update()
    seq.update()
    assert C.text(par.core[3], 3584) == "DAT #0, #5"     # copied the core as it was at the start
    assert C.text(seq.core[3], 3584) == "DAT #0, #6"     # copied after the ADD
    for order, b in (([0, 1], 7), ([1, 0], 5)):
        v = machine(["MOV #5, $2", "MOV #7, $1"], pointers=order)
        v.update()
        assert v.core[2][4] == b
        assert v.stats["overwritten_writes"] == 1


def test_operation_radius_folds_relative_addresses():
    v = machine(["MOV $0, $150"], M=1024, operation_radius=100)
    v.update()
    assert v.core[-51] == v.core[0] and v.core[150] == C.parse(BLANK, 1024)   # (150 + 100) mod 201 - 100 = -51
    assert v.fold(100) == 100 and v.fold(-100) == -100 and v.fold(101) == -100
    v = machine(["JMP $300"], M=1024, operation_radius=100)
    v.update()
    assert v.queue == [(300 + 100) % 201 - 100]


def test_mov_mutation_writes_uniform_random_instructions():
    # "each word has an equal probability to change into any of the ten instructions; operands at random"
    v = machine(["MOV $0, $1"], mutation_rate=1.0)
    written = []
    for _ in range(20000):
        w = {}
        assert v.execute(0, w) == ((1,), 0)
        written.append(w[1])
    assert v.stats["mutations"] == 20000
    freq = np.bincount([w[0] for w in written], minlength=10) / len(written)
    assert np.all(np.abs(freq - 0.1) < 0.01)
    modes = np.bincount([w[1] for w in written] + [w[3] for w in written], minlength=4) / (2 * len(written))
    assert np.all(np.abs(modes - 0.25) < 0.01)
    # without mutation the MOV copies exactly
    clean = machine(["MOV $0, $1"], mutation_rate=0.0)
    w = {}
    clean.execute(0, w)
    assert w == {1: clean.core[0]}


def test_pointer_injection_rate():
    v = machine(["JMP #0"] * 3584, pointers=[], M=3584, injection_rate=0.05, queue_len=10000)
    for _ in range(4000):
        v.update()
    assert 160 < v.stats["injections"] < 240            # 0.05 per core update
    full = machine(["JMP #0"], pointers=[0], queue_len=1, injection_rate=1.0)
    full.update()
    assert full.stats["injections_refused"] == 1 and full.queue == [0]


# --- eqs. 1-3: computational resources -----------------------------------------
def test_execution_removes_one_exec_proportionally_eqs_2_3():
    v = machine(["JMP $0"], M=64, resources=True, resource_radius=2, resource_max=0.5, resource_influx=0.0)
    v.r = [0.5] * 64
    v.r[1], v.r[62] = 0.1, 0.3
    before = [v.r[i % 64] for i in range(-2, 3)]
    S = sum(before)                                     # 1.9 >= 1: the JMP executes
    v.update()
    after = [v.r[i % 64] for i in range(-2, 3)]
    assert sum(before) - sum(after) == pytest.approx(1.0)
    assert after == pytest.approx([x * (S - 1) / S for x in before])
    assert v.executed[C.JMP] == 1
    # the neighbourhood now holds 0.9 < 1 exec: the pointer waits in place
    v.update()
    assert v.executed[C.JMP] == 1 and v.stats["waits"] == 1 and v.queue == [0]


def test_refill_is_capped_at_r_max_eq_1():
    v = machine(["DAT #0"], M=32, pointers=[], resources=True, resource_max=0.5, resource_influx=0.3)
    v.r = [0.0] * 32
    v.low = set(range(32))
    v.update()
    assert v.r == pytest.approx([0.3] * 32)
    v.update()
    assert v.r == [0.5] * 32


@pytest.mark.parametrize("env, rate", [("desert", 11 * 0.1), ("jungle", 3.0)])
def test_resource_influx_limits_execution_rate_table_2(env, rate):
    # Table 2: desert R_res = 5, dr = 0.1; jungle R_res = 3, dr = 0.5; r_max = 0.5.
    R, dr = (5, 0.1) if env == "desert" else (3, 0.5)
    v = machine(["JMP $0"], M=256, pointers=[0] * 8, resources=True, resource_radius=R,
                resource_influx=dr, resource_max=0.5)
    for _ in range(100):
        v.update()
    start = v.executed[C.JMP]
    for _ in range(1000):
        v.update()
    # eight trapped pointers share one neighbourhood of 2 R_res + 1 addresses. Desert: it gains
    # 11 * 0.1 = 1.1 execs per update and never reaches r_max, so 1.1 executions per update on average.
    # Jungle: it holds at most 7 * 0.5 = 3.5 execs; 3 executions leave 0.5, the capped refill restores 3.5.
    assert (v.executed[C.JMP] - start) / 1000 == pytest.approx(rate, rel=0.02)


# --- the observed network ----------------------------------------------------------
def test_default_network_is_an_observed_conserving_run():
    net = generate_network(ID, seed=3)
    assert net.status == "observed" and net.reactions
    assert [s.id for s in net.species] == list(C.OPS)
    assert sum(net.initial_state.values()) == net.params["core_size"]
    for r in net.reactions:
        assert r.count > 0
        assert sum(r.reactants.values()) == sum(r.products.values())     # cells are conserved
        assert r.reactants != r.products
    assert any(r.catalysts for r in net.reactions)
    a = net.extras["analysis"]
    assert sum(a["final_composition"].values()) == net.params["core_size"]
    assert a["events"]["writes"] == a["events"]["elastic_writes"] + sum(r.count for r in net.reactions)
    assert a["final_pointers"] <= net.params["queue_len"]
    assert sum(n * int(k) for runs in a["final_run_lengths"].values() for k, n in runs.items()) == 1024


def test_word_species_carry_the_program_text():
    net = generate_network(ID, seed=2, species="word", updates=200, core_size=256, queue_len=32)
    ids = {s.id for s in net.species}
    for s in net.species:
        assert s.id == C.word_id(C.parse(s.structure, 256), 256)
    for r in net.reactions:
        assert set(r.reactants) | set(r.products) <= ids
        assert sum(r.reactants.values()) == sum(r.products.values())


def test_seed_is_reproducible_and_parameters_checked():
    a = generate_network(ID, seed=5, updates=300)
    b = generate_network(ID, seed=5, updates=300)
    assert a.to_dict() == b.to_dict()
    assert a.extras["coreworld"]["seed_program"] == "mice"
    with pytest.raises(ValueError, match="core_size"):
        C.generate(SimpleNamespace(core_size=4, seed_program="mice"), np.random.default_rng(0))


# --- observations of section 3 --------------------------------------------------------
def test_jungle_is_more_active_than_desert():
    # Table 2 parameters on a small core: a jungle sustains far more executions and splits
    common = dict(core_size=1024, queue_len=64, updates=1500, operation_radius=100, seed_program="mice", seed=0)
    jungle = generate_network(ID, resource_radius=3, resource_influx=0.5, **common)
    desert = generate_network(ID, resource_radius=5, resource_influx=0.1, **common)

    def activity(net):
        ex = net.extras["analysis"]["series"]["executions"]
        return np.mean(ex[len(ex) // 2:])

    assert activity(jungle) > 1.5 * activity(desert)
    assert jungle.extras["analysis"]["events"]["splits"] > desert.extras["analysis"]["events"]["splits"]


def test_engineered_mice_is_too_brittle_for_the_noisy_core():
    # "any human engineered organism ... were too brittle to survive in the noisy VENUS universe"
    M = 3584
    body = [C.parse(s, M) for s in C.MICE[1:7]]

    def intact_after(pmut, ppoint, updates=600):
        rand = C._buffered(np.random.default_rng(0))
        core = C.random_core(rand, M)
        for k, line in enumerate(C.MICE):
            core[100 + k] = C.parse(line, M)
        v = C.Venus(core, rand, queue_len=220, resource_radius=3, resource_influx=0.5, resource_max=0.5,
                    mutation_rate=pmut, injection_rate=ppoint, record=False)
        v.queue = [101]
        for _ in range(updates):
            v.update()
        return sum(all(v.core[(i + k) % M] == body[k] for k in range(6)) for i in range(M))

    assert intact_after(0.0, 0.0) >= 10        # without noise the copies persist and multiply
    assert intact_after(0.05, 0.05) == 0       # Table 2 noise destroys every copy
