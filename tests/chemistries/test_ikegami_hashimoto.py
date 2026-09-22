"""Machine-tape chemistry, reconstructed from Ikegami & Hashimoto.

[AL] Active mutation in self-reproducing networks of machines and tapes,
     Artificial Life 2(3):305-318, 1995 (author's manuscript, May 1996).
[ECAL] Coevolution of machines and tapes, ECAL 95, LNAI 929:234-245.
"""

from collections import Counter

import numpy as np
import pytest
from chemart.simulate import rhs

from chemart import generate_network
from chemart.chemistries import ikegami_hashimoto as ih


# --- the machine ---------------------------------------------------------------
def test_fig1_rewriting_and_translation():
    # [AL] fig. 1: M ebbd reads T 1b = [0011011] from state 1, writes [1110110], which folds back into M ebbd
    head, tprime, mprime = 0xB, [1, 1, 1, 0], [1, 0, 1, 1]
    assert ih.fields(0xEBBD) == (tprime, mprime, head, 0xD)     # table rows TM T'M': 0011 0110 1011 1101
    assert ih.bind(0xEBBD, 0b0011011) == (3, 6)                  # head 1011 at site 3, tail 1101 at site 2
    trace = [(f"{t:07b}", s) for t, s in ih.rewrite_trace(0xEBBD, 0b0011011, 1)]
    assert trace == [("0011011", 1), ("0010011", 1), ("0010111", 0), ("0010111", 1),
                     ("0010110", 1), ("1010110", 0), ("1110110", 1)]
    assert ih.translate(0b1110110) == 0xEBBD
    assert ih.react(0xEBBD, 0b0011011, 1) == (0xEBBD, 0b1110110, 3, 6, 5)


def test_description_tapes_of_the_published_machines():
    # [AL] figs. 2 and 7, [ECAL] fig. 5: M1002 with T1, M3006 with T5, M1222 with T3
    assert ih.translate(0x01) == 0x1002
    assert ih.translate(0x05) == 0x3006
    assert ih.translate(0x03) == 0x1222
    assert ih.necklace(0x41) == 0x03, "[ECAL] fig. 1 calls M1222's tape T41, a rotation of T3"


@pytest.mark.parametrize("machine, tape", [(0x1002, 0x01), (0x2004, 0x01), (0xDFFB, 0x3F), (0x9DD3, 0x1D),
                                           (0xBDD7, 0x37)])   # printed 'Mbdd1', see catalog decisions
def test_published_self_replicating_pairs(machine, tape):
    # [AL] sec. 4: the five self-replicating loops; the tape is named by its smallest rotation
    rotations = [ih.rotate(tape, r) for r in range(7)]
    assert any(ih.react(machine, t, 1)[:2] == (machine, t) for t in rotations if ih.bind(machine, t))


def test_bdd1_cannot_read_t37():
    assert all(ih.bind(0xBDD1, ih.rotate(0x37, r)) is None for r in range(7))


def test_m1002_is_produced_by_the_largest_variety_of_machines():
    # [AL] sec. 5: "the most large variety of machines generate the machine M1002"
    tables = np.arange(256)
    tp = np.array([(tables >> (7 - i)) & 1 for i in range(4)])
    mp = np.array([(tables >> (3 - i)) & 1 for i in range(4)])
    weights = 1 << (15 - np.arange(16))
    order = [i % 7 for i in ih._TRANSLATION]   # machine bit j <- tape site
    chunks = []
    for t in range(128):
        sites = np.array([ih.tape_bit(t, s) for s in range(7)])
        for head in range(16):
            for tail in range(16):
                found = ih.bind((head << 4) | tail, t)
                if found is None:
                    continue
                h, length = found
                for state in (0, 1):
                    tape = np.tile(sites, (256, 1))
                    m = np.full(256, state)
                    for i in range(length):
                        s = (h + i) % 7
                        k = 2 * tape[:, s] + m
                        tape[:, s], m = tp[k, tables], mp[k, tables]
                    produced = tape[:, order] @ weights
                    chunks.append(produced * 65536 + ((tables << 8) | (head << 4) | tail))
    pairs = np.unique(np.concatenate(chunks))
    variety = np.bincount(pairs // 65536, minlength=65536)
    assert variety[0x1002] == variety.max()
    assert ih.react(0x0000, 0x01, 0) and variety[0x1002] > variety[0x0000]


# --- networks ------------------------------------------------------------------
def test_fig7_parasitic_chain():
    # [AL] fig. 7a-c: the loop M1002/T1 is exploited by M3006 (made by M1002 from T5), then by M1222 (from T3)
    net = generate_network("ikegami-hashimoto", method="closure", machines=["1002", "3006", "1222"],
                           tapes=["01", "05", "03"])
    assert net.status == "complete"
    got = {(tuple(sorted(r.reactants.items())), tuple(sorted(r.products.items()))) for r in net.reactions}
    assert got == {
        ((("M1002", 1), ("T01", 1)), (("M1002", 2), ("T01", 2))),
        ((("M1002", 1), ("T05", 1)), (("M1002", 1), ("M3006", 1), ("T05", 2))),
        ((("M3006", 1), ("T03", 1)), (("M1222", 1), ("M3006", 1), ("T03", 2))),
    }
    assert all(r.rate == {"law": "mass-action", "k": 0.6, "frame_length": r.rate["frame_length"]}
               for r in net.reactions)


def test_closure_products_are_translations():
    net = generate_network("ikegami-hashimoto", seed=5, method="closure")
    assert net.status == "complete"
    tapes = [s for s in net.species if s.id.startswith("T")]
    assert len(tapes) <= 128
    for r in net.reactions:
        new = Counter(r.products) - Counter(r.reactants)
        (m2,) = [s for s in new if s.startswith("M")] or [s for s in r.reactants if s.startswith("M")]
        (t2,) = [s for s in new if s.startswith("T")] or [s for s in r.reactants if s.startswith("T")]
        assert ih.translate(int(t2[1:], 16)) == int(m2[1:], 16)


def test_rate_equations_are_eqs_4_to_6():
    # continuous-time eqs. 4-6 with c = d, no noise or integer parts:
    # df_i = c (sum over reading (k, j, state) -> i of f_k f_j / 2  -  f_i W),  W = sum over reading pairs of f_k f_j
    net = generate_network("ikegami-hashimoto", seed=2, method="closure", c=0.6)
    ids, f = rhs(net)
    rng = np.random.default_rng(0)
    x = rng.uniform(0.1, 1.0, len(ids))
    is_m = np.array([s.startswith("M") for s in ids])
    x[is_m] /= x[is_m].sum()
    x[~is_m] /= x[~is_m].sum()
    conc = dict(zip(ids, x))
    machines = [s for s in ids if s.startswith("M")]
    tapes = [s for s in ids if s.startswith("T")]
    production, W = Counter(), 0.0
    for m in machines:
        for t in tapes:
            if ih.bind(int(m[1:], 16), int(t[1:], 16)) is None:
                continue
            w = conc[m] * conc[t]
            W += w
            for state in (0, 1):
                m2, t2, *_ = ih.react(int(m[1:], 16), int(t[1:], 16), state)
                production[ih.machine_id(m2)] += w / 2
                production[ih.tape_id(t2)] += w / 2
    ours = dict(zip(ids, f(0.0, x)))
    for s in ids:
        assert ours[s] == pytest.approx(0.6 * (production[s] - conc[s] * W), abs=1e-12), s


# --- population dynamics -------------------------------------------------------
def test_minimal_loop_is_a_fixed_point_without_noise():
    # [AL] sec. 4, attractor 1: a minimal self-replicating loop with zero active mutation
    net = generate_network("ikegami-hashimoto", machines=["1002"], tapes=["01"], noise=0.0, generations=30)
    assert net.status == "observed" and len(net.reactions) == 1
    (r,) = net.reactions
    assert r.to_text() == "M1002 + T01 -> 2 M1002 + 2 T01  [mass-action k=0.6 frame_length=4]  (x30)"
    assert net.extras["final_state"] == {"M1002": 1000, "T01": 1000}
    a = net.extras["analysis"]
    assert set(a["active_mutation"]) == {0.0} and set(a["reading_length"]) == {4.0}


def test_machine_without_description_tape_is_washed_out():
    # [AL] sec. 3: "a machine without description tape is unstable and smoothly removed"
    net = generate_network("ikegami-hashimoto", machines=["1002", "3006"], tapes=["01"], noise=0.0, generations=12)
    assert net.initial_state == {"M1002": 500, "M3006": 500, "T01": 1000}
    final = net.extras["final_state"]
    assert set(final) == {"M1002", "T01"} and final["M1002"] >= 995   # integer parts stop just below N


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_noise_brings_the_published_parasites(seed):
    # [AL] fig. 2a / [ECAL] fig. 1a (noise 0.04): M3006 and M1222 with their tapes T5 and T3 invade the loop,
    # with bursts of active mutation
    net = generate_network("ikegami-hashimoto", seed=seed, machines=["1002"], tapes=["01"], noise=0.04,
                           generations=150)
    names = net.extras["paper_names"]
    made = {}   # parasite machine -> circular names of the tapes produced with it
    for r in net.reactions:
        new = Counter(r.products) - Counter(r.reactants)
        for m in ("M3006", "M1222"):
            if m in new:
                (tape,) = [s for s in new if s.startswith("T")]   # the written tape is always one new copy
                made.setdefault(m, set()).add(names[tape])
    assert made == {"M3006": {"T5"}, "M1222": {"T3"}}
    assert max(net.extras["analysis"]["active_mutation"]) > 0
    for i in net.extras["noise_induced"]:
        assert net.reactions[i].rate is None
    assert max(net.extras["analysis"]["distinct_machines"]) > 5


def test_turning_noise_off_stops_innovation():
    net = generate_network("ikegami-hashimoto", seed=0, machines=["1002"], tapes=["01"], noise=0.04,
                           generations=120, noise_off=80)
    a = net.extras["analysis"]
    after = a["distinct_machines"][90:]
    assert all(x >= y for x, y in zip(after, after[1:])), "no new machines without noise"


def test_observed_network_bookkeeping():
    net = generate_network("ikegami-hashimoto", seed=3)
    assert net.status == "observed" and net.outflow == "constant-total"
    assert sum(n for s, n in net.initial_state.items() if s.startswith("M")) == 1000
    assert sum(n for s, n in net.initial_state.items() if s.startswith("T")) == 1000
    for i, r in enumerate(net.reactions):
        assert r.count >= 1
        assert (r.rate is None) == (i in net.extras["noise_induced"])
        if r.rate:
            assert r.rate["k"] in (0.3, 0.6) and 1 <= r.rate["frame_length"] <= 6
    for s in net.species:
        width = 16 if s.id.startswith("M") else 7
        assert len(s.structure) == width and int(s.structure, 2) == int(s.id[1:], 16)
    final = net.extras["final_state"]
    assert sum(n for s, n in final.items() if s.startswith("M")) <= 1000


def test_bad_parameters():
    with pytest.raises(ValueError, match="16-bit"):
        generate_network("ikegami-hashimoto", machines=["12345"])
    with pytest.raises(ValueError, match="7-bit"):
        generate_network("ikegami-hashimoto", tapes=["80"])
    with pytest.raises(ValueError, match="noise_off"):
        generate_network("ikegami-hashimoto", generations=10, noise_off=20)
