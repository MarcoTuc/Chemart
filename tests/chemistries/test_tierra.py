"""Tierra: Ray (1991), "An approach to the synthesis of life" (book 10.6.4, ref [701]).

Published results reproduced: the instruction set and the ancestor listing
(appendices B and C), the ancestor's 839 and 813 instructions per replication
(appendix C; 827 and 809 in the Tierra 6.02 genebank), its self-examination
(appendix C register comments), template addressing (book 10.6.4 and table
10.3), write protection, the parasite 0045aaa (instruction 42 mutated; no
self-replication, replication within the search limit of a host), the
0046aaa/0064aaa mutualists, and the soup filling to 80% in about 400000
instructions with about 375 creatures.
"""

from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import tierra as T

ID = "tierra"

#: Machine-code column of appendix C (hexadecimal opcodes of 0080aaa).
APPENDIX_C_HEX = """
01 01 01 01 04 02 03 03 18 1c 00 00 00 00 07 19 1d 00 00 00 01 08 06 01 01 00 01 1e 16 00
00 01 01 1f 14 00 00 01 00 05 01 01 00 00 0c 0d 0e 01 00 01 00 1a 0a 05 14 00 01 00 00 08
09 14 00 01 00 01 05 01 00 01 01 12 11 10 17 01 01 01 00 05
"""
#: gb0/0045aaa.tie (Tierra 6.02 names).
PARASITE_602 = """
nop1 nop1 nop1 nop1 zero not0 shl shl movDC adrb nop0 nop0 nop0 nop0 subAAC movBA adrf nop0
nop0 nop0 nop1 incA subCAB nop1 nop1 nop0 nop1 mal call nop0 nop0 nop1 nop1 divide jmpo nop0
nop0 nop1 nop0 ifz nop1 nop1 nop1 nop0 pushA
"""


def params(**kw):
    base = dict(soup_size=20000, instructions=20000, inoculum={"0080aaa": 1}, gap=0, search_limit=400,
                reaper_threshold=0.8, slicer="size-dependent", slice_power=1.0, slice_size=25,
                cosmic_ray_interval=0, copy_error_interval=0, flaw_interval=0, semantics="ray1991")
    base.update(kw)
    return SimpleNamespace(**base)


def run(seed=0, samples=0, **kw):
    p = params(**kw)
    soup, placed = T.build(p, np.random.default_rng(seed))
    history = soup.run(p.instructions, samples=samples)
    return soup, history


def events(soup):
    return [(lhs, rhs, n) for lhs, rhs, n in soup.events.values()]


# --- the instruction set and the published genomes ---------------------------------
def test_instruction_set_and_ancestor_listing():
    # appendix B: 32 instructions; a few anchors of the opcode table
    assert len(T.MNEMONICS) == len(set(T.MNEMONICS)) == 32
    assert (T.OPCODES["nop_0"], T.OPCODES["jmp"], T.OPCODES["mov_iab"], T.OPCODES["divide"]) == (0x00, 0x14, 0x1a, 0x1f)
    assert all(T.OPCODES[a] == T.OPCODES[b] for a, b in zip(T.MNEMONICS, T.MNEMONICS_602))
    # appendix C: the 80-instruction ancestor, 48 of which are no-operations
    assert T.ANCESTOR == bytes(int(h, 16) for h in APPENDIX_C_HEX.split())
    assert len(T.ANCESTOR) == 80 and sum(b < 2 for b in T.ANCESTOR) == 48
    # the parasite: low bit of instruction 42 flipped, 45 instructions (gb0/0045aaa.tie)
    assert T.PARASITE == T.assemble(PARASITE_602)
    assert T.PARASITE[:45] == bytes(b ^ (i == 42) for i, b in enumerate(T.ANCESTOR[:45]))
    assert [len(T.GENOMES[g]) for g in ("0079aab", "0046aaa", "0064aaa")] == [79, 46, 64]


def test_ancestor_self_examination_registers():
    # appendix C comments after instruction 22: ax = end of mother, bx = start, cx = size, dx = template size
    p = params(inoculum={"0080aaa": 1}, gap=0)
    soup, _ = T.build(p, np.random.default_rng(0))
    soup.release(0, 80)                           # move the ancestor to address 1000
    cell = next(iter(soup.cells.values()))
    soup.soup[0:80] = bytes(80)
    soup.soup[1000:1080] = T.ANCESTOR
    soup.claim(1000, 80, cell.id)
    cell.mp = cell.ip = 1000
    while cell.ip != 1023:
        soup.run_slice(cell, 1)
    assert (cell.ax, cell.bx, cell.cx, cell.dx) == (1080, 1000, 80, 4)
    assert cell.errors == 0


@pytest.mark.parametrize("semantics, counts", [("ray1991", [839, 813]), ("tierra602", [827, 809])])
def test_ancestor_replicates_exactly_in_published_cycles(semantics, counts):
    # Ray 1991 appendix C: 839 instructions for the first replication, 813 for each later one;
    # gb0/0080aaa.tie (Tierra 6.02): 827 and 809
    soup, _ = run(inoculum={"0080aaa": 1}, instructions=5000, semantics=semantics)
    assert soup.metabolism["0080aaa"] == counts
    assert list(soup.genomes) == ["0080aaa"] and soup.stats["errors"] == 0
    daughters = [c for c in soup.cells.values() if c.mp != 0]
    assert daughters and all(bytes(soup.soup[c.mp:c.mp + c.ms]) == T.ANCESTOR for c in daughters)
    assert events(soup)[0] == ({"0080aaa": 1}, {"0080aaa": 2}, soup.stats["births"])


# --- template addressing, protection, the reaper ------------------------------------
def test_template_addressing():
    soup = T.Soup(np.random.default_rng(0), 400, search_limit=100)
    # book 10.6.4: JMP NOP0 NOP1 NOP0 targets the nearest NOP1 NOP0 NOP1
    soup.soup[200:204] = T.assemble("jmp nop_0 nop_1 nop_0")
    soup.soup[230:233] = T.assemble("nop_1 nop_0 nop_1")         # forward, 25 steps from 205
    soup.soup[169:172] = T.assemble("nop_1 nop_0 nop_1")         # backward, 28 steps from 197
    assert soup.search(200, 3, "o") == 233                        # IP goes just after the pattern
    assert soup.search(200, 3, "b") == 172 and soup.search(200, 3, "f") == 233
    soup.soup[169:172] = bytes(3)
    soup.soup[172:175] = T.assemble("nop_1 nop_0 nop_1")         # backward, 25 steps: a tie
    assert soup.search(200, 3, "o") == 233                        # forward wins ties
    soup.soup[230:233] = bytes(3)
    assert soup.search(200, 3, "o") == 175
    assert soup.search(200, 3, "o", flaw=1) == 176                # a flaw finds it one place off
    # not complementary, or beyond the search limit: not found
    assert soup.search(200, 3, "f") == -1
    soup.soup[172:175] = bytes(3)
    soup.soup[320:323] = T.assemble("nop_1 nop_0 nop_1")
    assert soup.search(200, 3, "o") == -1
    soup.search_limit = 200
    assert soup.search(200, 3, "o") == 323
    # a template may match inside a longer run of nops
    soup.soup[320:323] = bytes(3)
    soup.soup[395:400] = T.assemble("nop_1 nop_1 nop_0 nop_1 nop_0")
    assert soup.search(200, 3, "f") == 399          # the match at 396..398
    # the soup is circular: a pattern may straddle its end
    soup.soup[395:400] = bytes(5)
    soup.soup[398:400], soup.soup[0] = T.assemble("nop_1 nop_0"), T.OPCODES["nop_1"]
    assert soup.search(200, 3, "f") == 1


@pytest.mark.parametrize("semantics, next_ip", [("ray1991", 201), ("tierra602", 204)])
def test_failed_jump_is_an_error_and_is_skipped(semantics, next_ip):
    p = params(soup_size=1000, inoculum={"jmp nop_0 nop_1 nop_0 " + "if_cz " * 20: 1}, semantics=semantics)
    soup, _ = T.build(p, np.random.default_rng(0))
    cell = next(iter(soup.cells.values()))
    soup.soup[0:24], soup.soup[200:224] = bytes(24), soup.soup[0:24]
    soup.release(0, 24)
    soup.claim(200, 24, cell.id)
    cell.mp = cell.ip = 200
    soup.run_slice(cell, 1)
    assert cell.errors == 1 and cell.ip == next_ip


def test_write_protection():
    # creatures may write in their own blocks and in free memory, never in another creature
    p = params(soup_size=1000, inoculum={"mov_iab " * 12: 1, "if_cz " * 12: 1})
    soup, _ = T.build(p, np.random.default_rng(0))
    a, b = soup.cells.values()
    a.bx = 0
    for ax, allowed in [(15, False), (500, True), (7, True)]:
        a.ip, a.ax, errors = 0, ax, a.errors
        before = soup.soup[ax]
        soup.run_slice(a, 1)
        assert (a.errors == errors) is allowed
        assert soup.soup[ax] == (T.OPCODES["mov_iab"] if allowed else before)


def test_reaper_queue_moves_on_errors_and_success():
    p = params(soup_size=2000, inoculum={"0080aaa": 3})
    soup, _ = T.build(p, np.random.default_rng(0))
    first, second, third = soup.reap
    third.errors = 1
    soup.up_reaper(third)                       # at least as many errors as the one above: up
    assert soup.reap == [first, third, second]
    soup.up_reaper(third)
    assert soup.reap == [third, first, second]
    soup.down_reaper(third)                     # more errors than the one below: stays
    assert soup.reap == [third, first, second]
    soup.down_reaper(first)                     # success with no more errors: down
    assert soup.reap == [third, second, first]


# --- ecology ---------------------------------------------------------------------------
def test_parasite_cannot_replicate_alone():
    # 0045aaa "is not able to self-replicate in isolated culture": its call finds no copy procedure
    soup, _ = run(inoculum={"0045aaa": 1}, instructions=20000)
    assert soup.stats["births"] == 0 and soup.stats["errors"] > 1000
    assert events(soup) == []


def test_parasite_replicates_with_the_copy_procedure_of_a_host():
    soup, _ = run(inoculum={"0080aaa": 1, "0045aaa": 1}, instructions=20000)
    rx = {(tuple(sorted(l.items())), tuple(sorted(r.items()))): n for l, r, n in events(soup)}
    parasitism = ((("0045aaa", 1), ("0080aaa", 1)), (("0045aaa", 2), ("0080aaa", 1)))
    assert rx[parasitism] == soup.stats["parasitic_births"] > 0
    assert rx[((("0080aaa", 1),), (("0080aaa", 2),))] > 0
    final = Counter(c.genotype for c in soup.cells.values())
    assert final["0045aaa"] > 1 and final["0080aaa"] > 1


def test_parasite_needs_a_host_within_the_search_limit():
    # "if it is within the search limit (generally set at 200-400 instructions) of the copy
    # procedure of a creature of genotype 0080aaa, it will match templates"
    far = dict(inoculum={"0080aaa": 1, "0045aaa": 1}, gap=1000, instructions=4000)
    near, _ = run(search_limit=2000, **far)
    out, _ = run(search_limit=400, **far)
    assert near.stats["parasitic_births"] >= 1
    assert out.stats["parasitic_births"] == 0 and out.stats["births"] >= 1


def test_dissected_ancestor_halves_replicate_only_together():
    # Ray 1991 (ECOLOGY, fig. 3): sizes 46 and 64, "neither could replicate when cultured alone,
    # but when cultured together, they both replicated"
    for g in ("0046aaa", "0064aaa"):
        alone, _ = run(inoculum={g: 1}, instructions=20000)
        assert alone.stats["births"] == 0
    both, _ = run(inoculum={"0046aaa": 1, "0064aaa": 1}, instructions=50000)
    rx = events(both)
    born = Counter()
    for lhs, rhs, n in rx:
        for g in rhs:
            born[g] += (rhs[g] - lhs.get(g, 0)) * n
    assert born["0046aaa"] > 5 and born["0064aaa"] > 5
    assert ({"0046aaa": 1, "0064aaa": 1}, {"0046aaa": 2, "0064aaa": 1}) in [(l, r) for l, r, _ in rx]
    assert ({"0064aaa": 1, "0046aaa": 1}, {"0064aaa": 2, "0046aaa": 1}) in [(l, r) for l, r, _ in rx]


def test_soup_fills_to_the_reaper_threshold_as_published():
    # "the system executes about 400,000 instructions in filling up the soup with about 375
    # individuals of size 80 ... the memory remains roughly 80% filled" (60000 soup)
    soup, history = run(seed=5, soup_size=60000, instructions=700000, samples=70,
                        cosmic_ray_interval=10000, copy_error_interval=1000, flaw_interval=26016)
    full = next(h for h in history if h["memory_fill"] > 0.79)
    assert 300000 <= full["instructions"] <= 500000
    assert 250 <= full["creatures"] <= 450
    assert all(h["memory_fill"] > 0.75 for h in history if h["instructions"] > full["instructions"])
    assert soup.stats["deaths"] > 0


# --- the observed network ----------------------------------------------------------------
def test_default_network_contains_ancestor_events_and_balances():
    net = generate_network(ID, seed=3)
    assert net.status == "observed" and net.initial_state == {"0080aaa": 1.0}
    text = [r.to_text() for r in net.reactions]
    assert any(t.startswith("0080aaa -> 2 0080aaa  (x") for t in text)
    assert any(t.startswith("0080aaa -> ∅") for t in text)
    genomes = {s.id: s.structure for s in net.species}
    assert genomes["0080aaa"] == T.disassemble(T.ANCESTOR)
    for sid, structure in genomes.items():
        assert len(structure.split()) == int(sid[:4])
    # every event is a replication, a replication using a partner, a death or a mutation,
    # and initial state + events = final population
    state = Counter({k: int(v) for k, v in net.initial_state.items()})
    for r in net.reactions:
        lhs, rhs = sum(r.reactants.values()), sum(r.products.values())
        assert (lhs, rhs) in {(1, 2), (2, 3), (1, 0), (1, 1)}
        if lhs == 2:
            assert all(r.products.get(s, 0) >= n for s, n in r.reactants.items())
        for s, n in r.reactants.items():
            state[s] -= n * r.count
        for s, n in r.products.items():
            state[s] += n * r.count
    a = net.extras["analysis"]
    assert {k: v for k, v in state.items() if v} == a["final_population"]
    assert a["instructions_executed"] == net.params["instructions"]
    assert a["births"] - a["deaths"] == sum(a["final_population"].values()) - 1
    # cosmic rays: about one per 10000 instructions
    assert 0.5 < a["cosmic_rays"] / (net.params["instructions"] / 10000) < 1.5
    assert a["replication_instructions"]["0080aaa"][0] in (838, 839, 840)


def test_same_seed_same_run_and_parameter_checks():
    kw = dict(soup_size=5000, instructions=60000)
    assert generate_network(ID, seed=4, **kw).to_dict() == generate_network(ID, seed=4, **kw).to_dict()
    assert generate_network(ID, seed=4, **kw).to_dict() != generate_network(ID, seed=5, **kw).to_dict()
    with pytest.raises(ValueError, match="unknown Tierra instruction"):
        generate_network(ID, inoculum={"0099zzz": 1})
    with pytest.raises(ValueError, match="positive integers"):
        generate_network(ID, inoculum={"0080aaa": 0})
    with pytest.raises(ValueError, match="cannot hold"):
        generate_network(ID, soup_size=300, inoculum={"0080aaa": 4})
    with pytest.raises(ValueError, match="minimum cell size"):
        generate_network(ID, inoculum={"nop_1 nop_0 divide": 1})
    net = generate_network(ID, seed=1, soup_size=2000, instructions=3000, cosmic_ray_interval=0,
                           copy_error_interval=0, flaw_interval=0,
                           inoculum={"nop1 nop1 nop1 nop1 " + T.disassemble(T.ANCESTOR[4:]): 1})
    assert list(net.initial_state) == ["0080aaa"]   # a listing of a known genome keeps its name
