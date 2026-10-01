"""Core War / Redcode (book 10.6.2): ICWS'94 as implemented by pMARS 0.9.2.

No warrior here is taken from pMARS: Imp, Dwarf and the Imp avalanche are
Dewdney's and the book's, the conformance probe and Quarry are written for
Chemart, and the three-warrior battles use randomly generated warriors.

Reference numbers marked "pMARS" come from pMARS 0.9.2 compiled with
-DEXT94 -DSERVER and instrumented to print, at the end of every round, the
steps executed, the survivors, the task counts and the whole core
(sha1 of the lines "opcode A_mode A_value B_mode B_value" in pMARS's own
encoding).
"""

import hashlib
from collections import Counter

import pytest

from chemart import generate_network
from chemart.chemistries import corewar as cw

IMP, DWARF, AVALANCHE = cw.WARRIORS["imp"], cw.WARRIORS["dwarf"], cw.WARRIORS["imp-avalanche"]

_POP = {cw.MOV: 0, cw.ADD: 1, cw.SUB: 2, cw.MUL: 3, cw.DIV: 4, cw.MOD: 5, cw.JMZ: 6, cw.JMN: 7, cw.DJN: 8,
        cw.CMP: 9, cw.SLT: 10, cw.SPL: 11, cw.DAT: 12, cw.JMP: 13, cw.SEQ: 14, cw.SNE: 15, cw.NOP: 16,
        cw.LDP: 17, cw.STP: 18}
_PMODE = {cw.IMM: 0, cw.DIR: 1, cw.BIND: 2, cw.BDEC: 3, cw.BINC: 4, cw.AIND: 0x82, cw.ADEC: 0x83, cw.AINC: 0x84}


def pmars_core_sha1(core):
    lines = (f"{_POP[c[0]] * 8 + c[1]} {_PMODE[c[2]]} {c[3]} {_PMODE[c[4]]} {c[5]}" for c in core)
    return hashlib.sha1("\n".join(lines).encode()).hexdigest()


def mars(sources, M, cycles, procs, max_length=100):
    return cw.Mars([cw.assemble(s, M, max_length) for s in sources], M, cycles, procs)


def cell(text, M=800):
    return cw.assemble(text, M).code[0]


# --- published behaviour of Imp and Dwarf ----------------------------------------------
def test_imp_moves_one_cell_per_cycle_leaving_a_trail():
    # Dewdney 1984 / book 10.6.2: MOV 0 1 copies itself to the next address and moves through memory
    m = mars([IMP], 800, 100, 800)
    res = m.play([0])
    imp = cell("MOV 0, 1")
    assert res["survivors"] == [0] and res["steps"] == 100
    assert m.core[:101] == [imp] * 101 and m.core[101:] == [cw.EMPTY] * 699


def test_dwarf_bombs_every_fourth_cell_and_spares_itself():
    # Karonen's guide: DAT bombs every 4 instructions; in a core divisible by 4 it loops around unharmed
    m = mars([DWARF], 800, 600, 800)                    # 200 loops of ADD, MOV, JMP
    res = m.play([0])
    assert res["survivors"] == [0]
    code = cw.assemble(DWARF, 800).code
    assert m.core[:4] == code                           # the 200th ADD brought the pointer back to #0
    bombs = {i: c for i, c in enumerate(m.core) if i >= 4 and c != cw.EMPTY}
    assert bombs == {3 + 4 * k: cell(f"DAT #0, #{4 * k}") for k in range(1, 200)}


@pytest.mark.parametrize("position, steps, survivors", [
    (300, 16000, [0, 1]),      # pMARS: tie after 8000 cycles
    (597, 1193, [1]),          # pMARS: Dwarf bombs Imp
    (598, 1187, [1]),
    (599, 1181, [1]),
])
def test_imp_vs_dwarf_matches_pmars(position, steps, survivors):
    res = mars([IMP, DWARF], 800, 8000, 800).play([0, position])
    assert res["steps"] == steps and res["survivors"] == survivors


def test_imp_subverts_dwarf_into_a_second_imp():
    # Dewdney 1984: "Dwarf will be subverted and become a second Imp ... the battle is a draw"
    net = generate_network("corewar", seed=1, positions=[300], rounds=1)
    imp, dwarf = (w["id"] for w in net.extras["warriors"])
    a = net.extras["analysis"]
    assert a["rounds"][0]["survivors"] == [imp, dwarf] and a["score"] == {imp: 1, dwarf: 1}
    subverted = [r for r in net.reactions if r.reactants.get(dwarf) and r.reactants.get("MOV.I$0,$1")
                 and r.products.get("MOV.I$0,$1", 0) == 2 and r.products.get(dwarf) == 1]
    assert subverted and a["events"][dwarf]["foreign_executions"] > 1000
    replication = {"reactants": {imp: 1, "DAT.F$0,$0": 1, "MOV.I$0,$1": 1}, "products": {imp: 1, "MOV.I$0,$1": 2}}
    assert any(r.reactants == replication["reactants"] and r.products == replication["products"]
               for r in net.reactions)


@pytest.mark.slow
def test_imp_vs_dwarf_all_positions_matches_pmars():
    # pMARS, KOTH rules scaled to core 800: over the 601 load addresses 100..700, Dwarf wins 104, 497 ties
    ws = [cw.assemble(IMP, 800), cw.assemble(DWARF, 800)]
    outcomes = Counter(tuple(cw.Mars(ws, 800, 8000, 800).play([0, f])["survivors"]) for f in range(100, 701))
    assert outcomes == {(0, 1): 497, (1,): 104}


def test_imp_avalanche_fills_the_task_queue():
    # book 10.6.2: SPL 2 / JMP -1 / MOV 0 1 "creates an avalanche of self-replicating MOV 0 1 programs"
    net = generate_network("corewar", warriors=["imp-avalanche"], max_processes=16, max_cycles=2000, rounds=1)
    wid = net.extras["warriors"][0]["id"]
    rnd = net.extras["analysis"]["rounds"][0]
    assert rnd["survivors"] == [wid] and rnd["tasks"][wid] == 16
    series = net.extras["analysis"]["series_round_1"]["tasks"][wid]
    assert max(series) == 16 and series == sorted(series)
    assert any(r.products.get(wid) == 2 for r in net.reactions)
    assert net.extras["final_state"]["MOV.I$0,$1"] > 16


# --- the ICWS'94 standard ------------------------------------------------------------------
def test_assembly_file_gives_the_standard_load_file():
    # ICWS'94 draft sec. 2.7 (assembly file) and sec. 3.5 (its load file)
    source = """;redcode
;name          Dwarf
;author        A. K. Dewdney
;assert        CORESIZE % 4 == 0
        ORG     start
step    EQU      4
target  DAT.F   #0,     #0
start   ADD.AB  #step,   target
        MOV.AB  #0,     @target
        JMP.A    start
        END
"""
    w = cw.assemble(source, 8000)
    assert (w.name, w.author) == ("Dwarf", "A. K. Dewdney")
    load = cw.assemble("ORG 1\nDAT.F #0, #0\nADD.AB #4, $-1\nMOV.AB #0, @-2\nJMP.A $-2, #0\n", 8000)
    assert w.start == load.start == 1
    # the draft's load file writes JMP.A $-2, #0; a missing operand of JMP assembles to $0 in pMARS
    assert w.code[:3] == load.code[:3] and w.code[3] == cell("JMP.A $-2, $0", 8000)


@pytest.mark.parametrize("text, expected", [
    # ICWS'94 app. A.2.1.2 default modifiers (pMARS asm.c)
    ("MOV 0, 1", "MOV.I $0, $1"), ("MOV #0, 1", "MOV.AB #0, $1"), ("MOV 0, #1", "MOV.B $0, #1"),
    ("ADD #4, 3", "ADD.AB #4, $3"), ("SUB 1, #2", "SUB.B $1, #2"), ("MUL 1, 2", "MUL.F $1, $2"),
    ("SLT 1, #2", "SLT.B $1, #2"), ("SLT #1, 2", "SLT.AB #1, $2"), ("DJN -1, <3", "DJN.B $-1, <3"),
    ("DAT -1", "DAT.F #0, $-1"), ("JMP -2", "JMP.B $-2, $0"), ("SPL 2", "SPL.B $2, $0"),
    ("seq *1, {2", "SEQ.I *1, {2"), ("nop }1, >2", "NOP.F }1, >2"),
    ("x EQU 3\n MOV x*2+(1-3)%2, <CORESIZE-x", "MOV.I $6, <-3"),
])
def test_default_modes_and_modifiers(text, expected):
    assert cw.instruction_text(cw.assemble(text, 8000).code[-1], 8000) == expected


def test_task_death_and_split_limit():
    # ICWS'94 sec. 5.5: DIV by zero removes the task; SPL queues only PC+1 when the queue is full
    m = mars(["DIV.AB #0, $1\nDAT.F #0, $5"], 800, 10, 800)     # A-value 0: B-target unchanged, the task dies
    res = m.play([0])
    assert res["survivors"] == [] and res["death_step"] == [1] and m.core[1] == cell("DAT.F #0, $5")
    m = mars(["DIV.F $1, $2\nDAT.F #0, #7\nDAT.F #3, #8"], 800, 10, 800)
    res = m.play([0])           # A-value (0, 7): A-number kept, B-number 8/7, then the task dies
    assert m.core[2] == cell("DAT.F #3, #1") and res["survivors"] == []
    m = mars(["DIV.X $1, $2\nDAT.F #0, #7\nDAT.F #3, #8"], 800, 10, 800)
    res = m.play([0])           # X crosses the pair: A-number 3/7 = 0, B-number kept, then the task dies
    assert m.core[2] == cell("DAT.F #0, #8") and res["survivors"] == []
    res = mars(["SPL.B $0, $0\nJMP.B $-1, $0"], 800, 50, 5).play([0])
    assert res["tasks"] == [5] and res["survivors"] == [0]      # the queue fills to the task limit


# --- a self-checking conformance probe ---------------------------------------------------
# Each check performs one ICWS'94 operation and compares the result against a cell holding
# the value the standard requires, jumping to a DAT when they differ. Surviving means every
# check passed, so the oracle is the standard rather than this MARS.
CONFORMANCE = """;redcode-94
;name Conformance
;author Chemart
;strategy Self-checking ICWS'94 probe: any failed check jumps to a DAT.
;assert CORESIZE >= 256
        ORG    start
start   MOV.X  xsrc, xdst
        SEQ.I  xdst, xwant
        JMP    fail
        MOV.I  isrc, idst
        SEQ.I  idst, iwant3
        JMP    fail
        SNE.I  nsrc, ndst
        JMP    fail
        MOV.AB #9, abdst
        SEQ.I  abdst, abwant
        JMP    fail
        MOV.AB #9, >iptr
        SEQ.I  iptr, iwant
        JMP    fail
        SEQ.I  itgt, iwant2
        JMP    fail
        MOV.AB #7, <dptr
        SEQ.I  dptr, dwant
        JMP    fail
        SEQ.I  dtgt, dwant2
        JMP    fail
        ADD.AB #1, wrap
        SEQ.I  wrap, wwant
        JMP    fail
        DJN.B  fail, dcount
        SEQ.I  dcount, dcwant
        JMP    fail
        SLT.AB #1, sltv
        JMP    fail
        SLT.AB #9, sltv
        JMP    done
        JMP    fail
done    JMP    done
fail    DAT.F  #0, #0
xsrc    DAT.F  #1, #2
xdst    DAT.F  #0, #0
xwant   DAT.F  #2, #1
isrc    SPL.B  #7, <3
idst    DAT.F  #0, #0
iwant3  SPL.B  #7, <3
nsrc    DAT.F  #0, #1
ndst    DAT.F  $0, #1
abdst   DAT.F  #5, #6
abwant  DAT.F  #5, #9
iptr    DAT.F  #0, $2
iwant   DAT.F  #0, $3
itgt    DAT.F  #0, #0
iwant2  DAT.F  #0, #9
dptr    DAT.F  #0, $3
dwant   DAT.F  #0, $2
dtgt    DAT.F  #0, #0
dwant2  DAT.F  #0, #7
wrap    DAT.F  #0, #799
wwant   DAT.F  #0, #0
dcount  DAT.F  #0, #1
dcwant  DAT.F  #0, #0
sltv    DAT.F  #0, #5
        END
"""

# (expected value, a wrong value) for every check the probe makes
CONFORMANCE_CHECKS = [
    ("xwant   DAT.F  #2, #1", "xwant   DAT.F  #1, #2", ".X crosses the operand pairs"),
    ("iwant3  SPL.B  #7, <3", "iwant3  DAT.F  #7, <3", ".I copies the opcode, not only the numbers"),
    ("ndst    DAT.F  $0, #1", "ndst    DAT.F  #0, #1", ".I comparison is sensitive to modes"),
    ("abwant  DAT.F  #5, #9", "abwant  DAT.F  #9, #9", ".AB writes the A-number into the B-number"),
    ("iwant   DAT.F  #0, $3", "iwant   DAT.F  #0, $2", "'>' increments after the operand is used"),
    ("iwant2  DAT.F  #0, #9", "iwant2  DAT.F  #0, #0", "'>' writes through the un-incremented pointer"),
    ("dwant   DAT.F  #0, $2", "dwant   DAT.F  #0, $3", "'<' decrements before the operand is used"),
    ("dwant2  DAT.F  #0, #7", "dwant2  DAT.F  #0, #0", "'<' writes through the decremented pointer"),
    ("wwant   DAT.F  #0, #0", "wwant   DAT.F  #0, #799", "arithmetic is modulo CORESIZE"),
    ("dcwant  DAT.F  #0, #0", "dcwant  DAT.F  #0, #1", "DJN decrements before it tests"),
    ("sltv    DAT.F  #0, #5", "sltv    DAT.F  #0, #0", "SLT skips when the A-value is lower"),
]


def test_conformance_probe_survives_only_on_a_standard_mars():
    m = mars([CONFORMANCE], 800, 5000, 100)
    assert m.play([0]) == {"steps": 5000, "survivors": [0], "tasks": [1], "death_step": [None]}


@pytest.mark.parametrize("expected, wrong, what", CONFORMANCE_CHECKS, ids=[c[2] for c in CONFORMANCE_CHECKS])
def test_every_conformance_check_is_live(expected, wrong, what):
    """Corrupting one expected value must kill the probe, or that check passes vacuously."""
    assert expected in CONFORMANCE
    res = mars([CONFORMANCE.replace(expected, wrong)], 800, 5000, 100).play([0])
    assert res["survivors"] == [], f"check passes vacuously: {what}"


# --- a stone whose stride is coprime to the core ----------------------------------------
QUARRY = """;redcode-94
;name Quarry
;author Chemart
;strategy A stone dropping DAT bombs in strides of 7. Unlike Dwarf's stride of 4,
;strategy 7 does not divide the core, so the walk comes around onto its own pointer.
;assert CORESIZE % 7 != 0
step    EQU    7
        ORG    start
bomb    DAT.F  #0, #0
start   MOV.I  bomb, @ptr
        ADD.AB #step, ptr
        JMP    start, $0
ptr     DAT.F  #0, #100
        END
"""
QUARRY_PTR = 4          # index of the `ptr` cell
QUARRY_LEN = 5


def test_quarry_drops_bombs_exactly_where_the_stride_predicts():
    # 100 loops of MOV/ADD/JMP, each dropping at ptr + ptr.B before ptr.B advances by 7
    m = mars([QUARRY], 800, 300, 800)
    assert m.play([0])["survivors"] == [0]
    empty_dat = cell("DAT.F #0, #0")
    bombs = sorted(i for i, c in enumerate(m.core) if c == empty_dat and i >= QUARRY_LEN)
    assert bombs == sorted((QUARRY_PTR + 100 + 7 * k) % 800 for k in range(100))
    assert m.core[:QUARRY_PTR] == cw.assemble(QUARRY, 800).code[:QUARRY_PTR]   # its code is untouched
    assert m.core[QUARRY_PTR] == cell("DAT.F #0, #0")                          # 100 + 700 = 0 mod 800


def test_quarry_bombs_itself_where_dwarf_never_does():
    # gcd(7, 800) = 1, so drop 100 lands on the pointer at address 4 and the walk collapses;
    # Dwarf's stride divides the core, so it loops around its own code unharmed forever.
    assert [k for k in range(800) if (QUARRY_PTR + 100 + 7 * k) % 800 < QUARRY_LEN][0] == 100
    res = mars([QUARRY], 800, 4000, 800).play([0])
    assert res["survivors"] == [] and res["death_step"][0] > 3 * 100
    assert mars([DWARF], 800, 4000, 800).play([0])["survivors"] == [0]


@pytest.mark.slow
def test_quarry_vs_dwarf_is_a_contest_over_every_load_position():
    ws = [cw.assemble(QUARRY, 800), cw.assemble(DWARF, 800)]
    outcomes = Counter(tuple(cw.Mars(ws, 800, 8000, 800).play([0, f])["survivors"]) for f in range(100, 701))
    assert outcomes == {(0,): 291, (1,): 293, (0, 1): 17}


RANDOM_BATTLES = [
    (["ORG 0\nSNE.BA #-76, {4\nJMN.A *1, >5\nMOV.I {4, {4\nLDP.A >1, <-1\nDJN.I #0, >-1\nMOV.I @1, *262\n"
      "SUB.AB *-2, {2\nSTP.I <-320, #-1\nNOP.BA {-6, <-2\nJMP.B $-9, $0\n",
      "ORG 0\nJMN.B @0, $3\nLDP.AB #3, #-183\nJMZ.B >6, }4\nJMN.X #-3, *4\nDJN.X *-2, *0\nJMP.A *-5, }-5\n"
      "DJN.A #-4, #-1\nSLT.B >-5, $1\nSPL.B *-4, {87\nJMP.B $-9, $0\n",
      "ORG 0\nDJN.B *5, $4\nSLT.B *5, @7\nJMN.AB }5, $1\nSPL.I *2, >5\nMOV.I *-3, *-4\nSPL.I >1, @-295\n"
      "SPL.F #-111, $0\nMUL.I *1, #-6\nMOV.AB *-1, #-5\nJMP.B $-9, $0\n"], 10,
     [([0, 129, 70], 4500, [1, 1, 16], "b7ca3a1548ea8f01e59203542dc02e20b2222564"),
      ([0, 235, 95], 4500, [1, 1, 16], "0fdc9bfe9b93f2b3d557b6cc4e67e041f5b27c04"),
      ([0, 226, 58], 4500, [1, 1, 16], "955093f97b9d96251945bba50135c64eccf5cc7f"),
      ([0, 230, 191], 4500, [1, 1, 16], "e1b35777b2f1ddd7f709af3fbbd7b4279f8c3b5d")]),
    (["ORG 0\nSTP.AB >5, }4\nJMN.F <5, <271\nJMP.F #-288, $5\nMOD.BA <0, @4\nJMZ.B *-95, }-3\nSPL.BA {2, {0\n"
      "JMP.F #-3, @-6\nSEQ.X >276, *-1\nJMP.B $-8, $0\n",
      "ORG 0\nSPL.X >1, >6\nSPL.F #1, <-121\nDJN.BA #1, >0\nMOV.I >-2, }-1\nSEQ.AB $1, *-2\nMOV.X #2, *1\n"
      "SPL.A }-1, }-351\nSPL.BA @-224, >-4\nJMP.B $-8, $0\n",
      "ORG 0\nDJN.A }5, *0\nSLT.A #6, >0\nJMP.F >5, {371\nADD.BA $-2, $4\nSEQ.X >-1, >-3\nDJN.B $1, }-1\n"
      "LDP.AB {-3, <0\nMOV.A $-5, *-2\nJMP.B $-8, $0\n"], 9,
     [([0, 129, 70], 4054, [1, 7, 0], "d7aeb5057333082b6539dfffce88d03b0e9ad038"),
      ([0, 235, 95], 4500, [1, 8, 1], "841814b1f68648db96e68bfa0726ce9ef2abba86"),
      ([0, 226, 58], 4110, [1, 7, 0], "95c36899e5b6930527df57513674d5e5feec8efd"),
      ([0, 230, 191], 4500, [1, 8, 1], "a44ad8a0605032f22d1618e62b0d3a3b961a3324")]),
    (["ORG 0\nSTP.B *4, $4\nSPL.A #0, $1\nADD.BA {-2, *-2\nSPL.B }0, *-3\nSUB.BA }0, #0\nSPL.A @1, }-5\n"
      "SUB.I *-1, *-5\nJMP.B $-7, $0\n",
      "ORG 0\nSPL.F {-346, $3\nMOV.A @-320, #3\nMOV.B $3, <-2\nMUL.F #-1, {2\nCMP.B >-4, *-1\nSPL.BA $-5, @1\n"
      "MUL.B @-5, #150\nJMP.B $-7, $0\n",
      "ORG 0\nNOP.B $3, {0\nSTP.BA #112, #5\nSNE.B <7, #1\nSLT.AB <139, >1\nSTP.X <-214, *1\nMOV.AB @-1, @-4\n"
      "MOV.BA <-5, <-4\nJMP.B $-7, $0\n"], 8,
     [([0, 129, 70], 3526, [16, 16, 0], "c50452134d0ddd25eb0edaa0729b9164589050a5"),
      ([0, 235, 95], 4500, [16, 7, 1], "1e3f89614ddd35cbc01319fbee8bc43b5165d3f5"),
      ([0, 226, 58], 4500, [16, 1, 1], "8e4e1b0797b3c1b45d5bd848943c1902f687ff22"),
      ([0, 230, 191], 3323, [16, 3, 0], "14dd2af39d944295e3cb644e5c461cd4e6af0d2b")]),
]


@pytest.mark.parametrize("sources, length, rounds", RANDOM_BATTLES)
def test_random_three_warrior_battles_match_pmars(sources, length, rounds):
    # every opcode, modifier and mode, 3 warriors, 4 rounds with rotating starter and persistent P-space;
    # pMARS -s 256 -c 1500 -p 16 -l <length+1> -d 20 -r 4
    m = mars(sources, 256, 1500, 16, max_length=length + 1)
    for r, (pos, steps, tasks, sha) in enumerate(rounds):
        res = m.play(pos, starter=r % 3)
        assert (res["steps"], res["tasks"], pmars_core_sha1(m.core)) == (steps, tasks, sha)


# --- the network -----------------------------------------------------------------------------
def test_observed_network_balances_and_conserves_cells():
    net = generate_network("corewar", seed=3)
    assert net.status == "observed" and all(r.count >= 1 and r.rate is None for r in net.reactions)
    rounds = len(net.extras["analysis"]["rounds"])
    state = Counter({s: rounds * int(v) for s, v in net.initial_state.items()})
    for r in net.reactions:
        state.update({s: r.count * n for s, n in r.products.items()})
        state.subtract({s: r.count * n for s, n in r.reactants.items()})
    assert +state == Counter(net.extras["final_state"]) and not -state
    (law,) = net.extras["conservation"]
    for r in net.reactions:
        assert sum(law["vector"][s] * n for s, n in r.reactants.items()) == \
               sum(law["vector"][s] * n for s, n in r.products.items())
    assert sum(law["vector"][s] * v for s, v in net.initial_state.items()) == 800
    ids = {w["id"] for w in net.extras["warriors"]}
    for s in net.species:
        if s.id not in ids:
            assert cw.instruction_id(cw.assemble(s.structure, 800).code[0], 800) == s.id


def test_multi_round_scores():
    net = generate_network("corewar", seed=5, rounds=6)
    a = net.extras["analysis"]
    imp, dwarf = (w["id"] for w in net.extras["warriors"])
    assert [r["starter"] for r in a["rounds"]] == [imp, dwarf] * 3
    assert a["wins"][imp] == 0                      # Imp never kills Dwarf (pMARS sweep above)
    assert a["score"][dwarf] == 3 * a["wins"][dwarf] + a["ties"][dwarf]
    assert a["wins"][dwarf] + a["ties"][dwarf] == 6


def test_bad_parameters():
    with pytest.raises(ValueError, match="built-in"):
        generate_network("corewar", warriors=["gemini"])
    with pytest.raises(ValueError, match="FOR"):
        generate_network("corewar", warriors=["x FOR 3\nMOV 0, 1\nROF\n"])
    with pytest.raises(ValueError, match="positions"):
        generate_network("corewar", positions=[5])
    with pytest.raises(ValueError, match="min_distance"):
        generate_network("corewar", min_distance=50)
    with pytest.raises(ValueError, match="max_length"):
        generate_network("corewar", warriors=["imp", "MOV 0, 1\n" * 5], max_length=4, min_distance=100)
    with pytest.raises(ValueError, match="unknown label"):
        generate_network("corewar", warriors=["imp", "JMP nowhere\n"])
