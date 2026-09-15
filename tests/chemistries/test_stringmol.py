"""Stringmol reproduces the spec v0.2 (YCS-2010-458) and ALife XII (2010) results.

The machine is a port of the authors' C++ source (franticspider/stringmol); it was
cross-checked step by step against the compiled source (see the catalog decisions).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import stringmol as sm

SEED = sm.SEED_REPLICASE
SP9 = sm.ALIFE12_SPECIES_9
SP29 = SP9.replace("$BLUBO^", "$BLUBP^")      # ALife XII fig. 7: single point mutation of species 9
SP30 = ("OBEQBXUUUDYGRHBBOSEOLHHHRLUEUOBLROORE$BLUBP^B>C$=?>$$BLUBO%}OYHO"
        "OBEQBXUUUDYGRHBBOSEOLHHHRLUEUOBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB")
SP31 = "BBOSEOLHHHRLUEUOBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"


def machine(seed=0, **kw):
    return sm.Machine(sm._buffered(np.random.default_rng(seed)), **kw)


def react(a, b, traceback="matrix"):
    out = machine(traceback=traceback).react(a.encode(), b.encode(), "possible", 20000)
    return [s.decode() for s in out]


# --- symbols and alignment -----------------------------------------------------
def test_complements_are_spec_table_2():
    # spec table 2 / book fig. 11.5: template codes 13 letters on, function codes self-complementary
    assert sm.complement("ABC$DEF>GHIJ^KLM=NOP?QRS}TUV%WXYZ") == "NOP$QRS>TUVW^XYZ=ABC?DEF}GHI%JKLM"


def test_substitution_matrix_is_the_alxii_table():
    # spec eqs. 7-9 over the source symbol loop (alignment.cpp default_table, config/ALXII.mtx)
    row_a = sm.SUBSTITUTION[ord("A")]
    assert [round(row_a[ord(c)], 3) for c in "ABC$DEF%"] == [1.0, -0.121, -0.243, -0.364, -0.485, -0.607, -0.728, -0.849]
    assert round(row_a[ord("M")], 3) == round(row_a[ord("N")], 3) == -1.941               # 16 steps either way
    assert round(row_a[ord("R")], 3) == -1.213 and round(row_a[ord("Z")], 3) == -0.121   # cyclic loop
    assert all(sm.SUBSTITUTION[ord(c)][ord(c)] == 0.5 for c in sm.FUNCTION_CODES)
    assert round(sm.INDEL, 3) == -1.333


def report_table():
    """Spec v0.2 table 3: the report's symbol order, mismatch -d/8, match 1 (0.5 for function codes), indel -4/3."""
    key = "ABC$DEF>GHIJ^KLM=NOP?QRS}TUV%WXYZ"
    sub = [None] * 256
    for i, a in enumerate(key):
        row = [0.0] * 256
        for j, b in enumerate(key):
            d = abs(i - j)
            row[ord(b)] = sm.f32((1.0 if a.isalpha() else 0.5) if i == j else -min(d, 33 - d) / 8)
        sub[ord(a)] = row
    return sub, sm.f32(-4 / 3)


def test_spec_appendix_b_alignments():
    table = report_table()
    assert round(table[0][ord("A")][ord("B")], 2) == -0.12 and round(table[0][ord("$")][ord("$")], 2) == 0.5
    cs, s = b"AAABHG>$CDGBAAA", b"ABABABHG$CDBAAC"
    # table 6: maximum 7.72, alignment --AAABHG>$CDGBAAa / abABABHG-$CD-BAAc
    score, s1, e1, s2, e2 = sm.smith_waterman(cs, s, "matrix", table)
    assert score == pytest.approx(7.72, abs=0.015)
    assert (s1, e1, s2, e2) == (0, 14, 2, 14)
    # table 7: the v0.1 "highest score" trace gives ababA--BHG$C-D-BAAc instead
    assert sm.smith_waterman(cs, s, "highest-neighbour", table)[1:] == (0, 14, 4, 14)


# --- reactions -----------------------------------------------------------------
@pytest.mark.parametrize("traceback", ["matrix", "highest-neighbour"])
@pytest.mark.parametrize("replicase", [sm.SEED_REPLICASE, sm.CONFIG_REPLICASE, SP9])
def test_seed_replicases_copy_themselves(replicase, traceback):
    # spec 11.2.1 / app. B.1, ALife XII: two copies bind, the active one builds and cleaves a third
    assert react(replicase, replicase, traceback) == [replicase] * 3
    assert sm.bind_probability(replicase, replicase, traceback) > 0


@pytest.mark.parametrize("traceback", ["matrix", "highest-neighbour"])
def test_alife12_origin_of_species_31(traceback):
    # ALife XII fig. 7 (spec fig. 18): 29 active on 9 pastes a copy of 9 over its own last symbol
    assert len(SP29) == 65 and SP30 == SP29[:64] + SP9
    for a, b in [(SP29, SP9), (SP9, SP29)]:
        out = react(a, b, traceback)
        assert Counter(out) in (Counter([SP30, SP9]), Counter([SP9, SP29, SP29]))
    assert react(SP29, SP9, traceback) == [SP30, SP9]
    # 30 + 9: the bind site shifts and the first 14 symbols of 9 are not copied
    assert SP31 == SP9[14:]
    for a, b in [(SP30, SP9), (SP9, SP30)]:
        assert Counter(react(a, b, traceback)) == Counter([SP30, SP9, SP31])


def test_spec_b1_mutant_makes_a_double_length_molecule():
    # spec app. B.1: with the matrix traceback, BLUBO -> BLUBP in the seed replicase creates a double-length molecule
    mutant = SEED.replace("$BLUBO^", "$BLUBP^")
    out = react(SEED, mutant)
    assert Counter(out) == Counter([mutant[:-1] + SEED, SEED])
    assert len(mutant[:-1] + SEED) == 2 * len(SEED) - 1


def test_role_goes_to_the_later_bind_site():
    # spec 7.4.1: the executing string is the one whose bind site starts furthest from its start
    m = machine()
    a, b = sm.Molecule(SP30.encode(), m.size), sm.Molecule(SP9.encode(), m.size)
    sw = m.get_sw(a, b)
    act, pas = m.set_exec(a, b, sw)
    assert act is a and sw[1] > sw[3]
    assert act.i == [sw[3], sw[1]] and act.it == act.rt == act.wt == act.ft == 1


def test_copy_mutations_move_along_the_symbol_loop():
    # spec 10.6.1: a substitution writes the next or previous symbol in the loop of the substitution matrix
    net = generate_network("stringmol", seed=3, method="soup", molecules={SEED: 30}, steps=300,
                           substitution_rate=0.01, indel_rate=0.0)
    loop = sm.KEY
    seed_id = sm.species_id(SEED)
    structure = {s.id: s.structure for s in net.species}
    mutants = 0
    for r in net.reactions:
        if r.reactants != {seed_id: 2}:
            continue
        for sid in Counter(r.products).elements():
            copy = structure[sid]
            assert len(copy) == len(SEED)
            for x, y in zip(copy, SEED):
                if x != y:
                    mutants += 1
                    assert abs(loop.index(x) - loop.index(y)) in (1, len(loop) - 1)
    assert mutants > 0


# --- container, soup, closure --------------------------------------------------
def test_container_observed_network_balances_exactly():
    net = generate_network("stringmol", seed=5, molecules={SEED: 40, SP9: 10}, steps=1500, energy_per_step=40,
                           cell_radius=40.0, decay=0.001, substitution_rate=0.002, indel_rate=0.0005)
    assert net.status == "observed" and net.outflow == pytest.approx(0.001)
    ex = net.extras
    pop = Counter(net.initial_state)
    for r, active in zip(net.reactions, ex["active_counts"]):
        assert r.count >= 1 and sum(active.values()) == r.count
        assert set(active) <= set(r.reactants)
        pop.update({s: r.count * n for s, n in r.products.items()})
        pop.subtract({s: r.count * n for s, n in r.reactants.items()})
    for rec in ex["aborted"]:
        for s in rec["released"]:
            pop[s] += rec["count"]
        for s in rec["reactants"]:
            pop[s] -= rec["count"]
    for rec in ex["in_progress"]:
        pop.update(rec["released"])
        pop.subtract(rec["reactants"])
    pop.subtract(ex["decayed"])
    assert +pop == Counter(ex["final_state"]) and not -pop
    assert len(net.species) > 2, "mutation creates new species"
    ids = {s.id for s in net.species}
    assert all(sm.species_id(s.structure) == s.id for s in net.species) and len(ids) == len(net.species)


def test_energy_and_decay_maintain_the_population():
    # ALife XII: the balance between energy influx and decay maintains the population
    common = dict(molecules={SEED: 60}, steps=4000, cell_radius=30.0, decay=1 / 650, substitution_rate=0.0,
                  indel_rate=0.0)
    fed = generate_network("stringmol", seed=2, energy_per_step=25, **common)
    pops = fed.extras["analysis"]["population"]
    assert not fed.extras["extinct"] and 25 <= min(pops[len(pops) // 2:]) and max(pops) <= 200
    assert sum(r.count for r in fed.reactions) > 100
    starved = generate_network("stringmol", seed=2, energy_per_step=0, **common)
    assert not starved.reactions
    assert starved.extras["analysis"]["population"][-1] < 60 * 0.2


def test_closure_of_the_seed_replicase():
    net = generate_network("stringmol", method="closure", molecules={SEED: 1})
    seed_id = sm.species_id(SEED)
    assert net.status == "complete" and [s.id for s in net.species] == [seed_id]
    assert [(r.reactants, r.products) for r in net.reactions] == [({seed_id: 2}, {seed_id: 3})]


def test_closure_contains_the_alife12_cascade():
    net = generate_network("stringmol", method="closure", molecules={SP9: 1, SP29: 1}, max_species=8)
    ids = {s: sm.species_id(s) for s in (SP9, SP29, SP30, SP31)}
    found = {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in net.reactions}
    assert (frozenset({ids[SP29]: 1, ids[SP9]: 1}.items()), frozenset({ids[SP30]: 1, ids[SP9]: 1}.items())) in found
    assert (frozenset({ids[SP30]: 1, ids[SP9]: 1}.items()),
            frozenset({ids[SP30]: 1, ids[SP9]: 1, ids[SP31]: 1}.items())) in found


def test_soup_keeps_the_population_size():
    net = generate_network("stringmol", seed=4, method="soup", molecules={SEED: 20, SP31: 5}, steps=200)
    assert net.status == "observed" and net.outflow == "constant-total"
    assert sum(net.extras["final_state"].values()) == 25


def test_bad_parameters():
    with pytest.raises(ValueError, match="33 Stringmol symbols"):
        generate_network("stringmol", molecules={"ABC*": 3})
    with pytest.raises(ValueError, match="positive integer"):
        generate_network("stringmol", molecules={SEED: 0})
    with pytest.raises(ValueError, match="agent_radius"):
        generate_network("stringmol", agent_radius=20.0, cell_radius=10.0)
    with pytest.raises(ValueError, match="max_length"):
        generate_network("stringmol", molecules={SEED: 2}, max_length=10)
