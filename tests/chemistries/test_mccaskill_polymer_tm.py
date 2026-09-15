"""McCaskill (1988), 'Polymer chemistry on tape', MPI internal report: doublet codes and published strings."""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.mccaskill_polymer_tm import (
    OPPOSITE, PARASITE, REPLICATOR, SAME, can_collide, complement, decode_rules, process, recognizons,
    species_id,
)

REP, PAR = REPLICATOR, PARASITE


# --- the doublet codes (report section 3, items 4-5) ---------------------------
def test_rule_doublet_code():
    # the replicator's three rules, at its overlapping initiators 3, 4, 5
    rules = dict(decode_rules(REP))
    assert sorted(rules) == [3, 4, 5]
    assert rules[4] == {"read": (0, 1), "write": OPPOSITE, "state": (0, 1), "next": OPPOSITE,
                        "read_motion": 1, "write_motion": 1}
    assert rules[5] == {"read": (0, 1), "write": 1, "state": (0,), "next": OPPOSITE,
                        "read_motion": 0, "write_motion": 0}
    # the initiator 111 supplies READ (11, both) and the first bit of WRITE
    (_, r), = decode_rules("11" + "10" + "01" + "00" + "11" + "00")
    assert r == {"read": (0, 1), "write": 1, "state": (0,), "next": SAME, "read_motion": 1, "write_motion": -1}
    (_, r), = decode_rules("11" + "11" + "10" + "11" + "00" + "01")
    assert r == {"read": (0, 1), "write": OPPOSITE, "state": (1,), "next": OPPOSITE,
                 "read_motion": -1, "write_motion": 0}


def test_recognizon_doublet_code():
    assert recognizons("11" + "10" + "11" + "00" + "01" + "1") == ["#1##0"]   # 11 #, 10 1, 00/11 #, 01 0
    assert recognizons("0000000000") == [], "not all sequences have recognizons"
    assert recognizons("111" + "01" * 20, R=16) == ["#" + "1" * 15]
    assert can_collide("11101", "11110") and not can_collide("1110101", "11101")


# --- published strings (report section 4) --------------------------------------
def test_published_replicator_replicates_through_its_complement():
    assert process(REP, REP) == [complement(REP)]
    assert process(REP, complement(REP)) == [REP]


def test_published_parasite_is_replicated_but_does_not_replicate():
    assert decode_rules(PAR) == [], "the parasite does not encode the replication rule"
    assert process(REP, PAR) == [complement(PAR)]
    assert process(REP, complement(PAR)) == [PAR]
    assert process(PAR, REP) == [] and process(PAR, PAR) == []


def test_parasite_recognizon_collides_with_the_replicator():
    assert "#" in recognizons(REP) and "#" in recognizons(PAR)
    assert can_collide(REP, PAR) and can_collide(PAR, REP) and can_collide(REP, REP)


def test_halting_and_cutoff():
    writer = "11" + "10" + "00" + "00" + "10" + "10"    # read both, write 1, both heads right
    assert process(writer, "010") == ["111"], "halts on the blank after the tape"
    looping = "11" + "10" + "00" + "00" + "01" + "10"   # read head stays: never halts
    assert process(looping, "0", max_steps=50) == [], "uncompleted processes have no effect"


# --- networks ------------------------------------------------------------------
def reactions_of(net):
    return {(frozenset(Counter(r.reactants).items()), frozenset(Counter(r.products).items())) for r in net.reactions}


def key(lhs, rhs):
    return (frozenset(Counter(map(species_id, lhs)).items()), frozenset(Counter(map(species_id, rhs)).items()))


def test_default_closure_contains_replication_and_parasitism():
    net = generate_network("mccaskill-polymer-tm")
    assert net.status == "complete"
    ids = {s.id for s in net.species}
    assert {species_id(x) for x in (REP, PAR, complement(REP), complement(PAR))} <= ids
    found = reactions_of(net)
    assert key((REP, REP), (REP, REP, complement(REP))) in found
    assert key((REP, PAR), (REP, PAR, complement(PAR))) in found
    assert all(r.catalysts == r.reactants for r in net.reactions), "processor and tape both survive"
    assert all(s.structure == s.id[1:] for s in net.species)


def test_closure_without_recognition_closes_the_plus_minus_cycle():
    patterns = reactions_of(generate_network("mccaskill-polymer-tm", strings=[REP]))
    assert key((REP, complement(REP)), (REP, complement(REP), REP)) not in patterns   # open doubt, see decisions
    free = reactions_of(generate_network("mccaskill-polymer-tm", strings=[REP], recognition="none"))
    assert key((REP, REP), (REP, REP, complement(REP))) in free
    assert key((REP, complement(REP)), (REP, complement(REP), REP)) in free


def test_soup_products_are_processing_results():
    net = generate_network("mccaskill-polymer-tm", seed=3, method="soup", steps=3000)
    assert net.status == "observed" and net.outflow == "constant-total"
    assert sum(net.initial_state.values()) == 200
    assert net.initial_state[species_id(REP)] == 20 and net.initial_state[species_id(PAR)] == 20
    assert net.reactions and net.extras["analysis"]["recognition_collisions"] > 0
    for r in net.reactions:
        a, b = [s[1:] for s in Counter(r.reactants).elements()]
        new = sorted(s[1:] for s in (Counter(r.products) - Counter(r.reactants)).elements())
        assert new in (sorted(process(a, b)), sorted(process(b, a)))
        assert r.count >= 1


def test_soup_replication_fires_but_the_cycle_needs_recognition_none():
    # report section 4 inoculum (10% replicator, 10% parasite in random 19-mers). Under the
    # reconstructed recognizons the replicator never meets its complement, so it declines
    # (the report's survival is not reproduced; see decisions). Without recognition the
    # plus/minus reaction can fire (it does for seed 1; the replicator still dies out later).
    net = generate_network("mccaskill-polymer-tm", seed=0, method="soup", steps=20000)
    found = {k: r.count for k, r in zip(reactions_of_list(net), net.reactions)}
    assert found.get(key((REP, REP), (REP, REP, complement(REP))), 0) >= 1
    assert key((REP, complement(REP)), (REP, complement(REP), REP)) not in found
    assert net.extras["final_state"].get(species_id(REP), 0) < net.initial_state[species_id(REP)]

    free = generate_network("mccaskill-polymer-tm", seed=1, method="soup", steps=5000, recognition="none")
    fired = set(reactions_of_list(free))
    assert key((REP, REP), (REP, REP, complement(REP))) in fired
    assert key((REP, complement(REP)), (REP, complement(REP), REP)) in fired


def reactions_of_list(net):
    return [(frozenset(Counter(r.reactants).items()), frozenset(Counter(r.products).items())) for r in net.reactions]


def test_bad_parameters():
    with pytest.raises(ValueError, match="binary strings"):
        generate_network("mccaskill-polymer-tm", strings=["012"])
    with pytest.raises(ValueError, match="error_rate"):
        generate_network("mccaskill-polymer-tm", error_rate=0.01)
    with pytest.raises(ValueError, match="population"):
        generate_network("mccaskill-polymer-tm", method="soup", population=10, inoculum_fraction=0.9)
