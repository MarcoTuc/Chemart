"""Laing's molecular machines: book 10.5.1 / fig. 10.5, Freitas & Merkle (2004) 4.8, Sipper (1998).

Laing's own papers were not accessible (see the catalog decisions). The tests check the
instruction semantics stated by the sources, and universal computation ([485]) through
the published busy-beaver results (Rado 1962; Lin & Rado 1965).
"""

from collections import Counter

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.laing_molecular_machines import parse_program, run, turing_program

INCREMENT = "CT1.W1.H.TT1.W0.R.CT1.W1"


def ones(r):
    return sum(piece.count("1") for piece in r.pieces)


# --- the instructions -------------------------------------------------------------
def test_book_fig_10_5_instruction_set_parses():
    # w(1), w(0), TT, CT, L, R, NOP (book fig. 10.5); H and detach (Freitas & Merkle 4.8)
    assert parse_program("W1.W0.TT1.CT1.L.R.NOP.H.D") == ["W1", "W0", "TT1", "CT1", "L", "R", "NOP", "H", "D"]
    with pytest.raises(ValueError, match="TT"):
        parse_program("CT1.W1")
    with pytest.raises(ValueError, match="instructions"):
        parse_program("W2")


def test_write_and_slide_recruit_a_zero_at_the_ends():
    # W(0)/W(1) set the contacted constituent; sliding off the tape attaches a new 0 (KSRM 4.8)
    assert run(["W1", "R", "R"], "0").pieces == ["100"]
    assert run(["L", "W1"], "0").pieces == ["10"]
    assert run(["W0"], "1").pieces == ["0"]


def test_conditional_transfer_reads_the_contacted_constituent():
    program = parse_program("CT1.W1.H.TT1.W0")      # on 1 jump to TT1 and write 0, else write 1
    assert run(program, "1").pieces == ["0"]
    assert run(program, "0").pieces == ["1"]
    assert run(program, "0").outcome == "halt" and run(program, "1").outcome == "end"


def test_detach_severs_from_predecessors():
    r = run(parse_program("R.R.D.W0"), "1111")
    assert r.pieces == ["11", "01"]


def test_non_halting_run():
    r = run(parse_program("TT1.W1.CT1"), "0", max_steps=50)
    assert r.outcome == "steps" and not r.halted


# --- universal computation: Turing machines as Laing machines ---------------------
BB2 = {("A", 0): (1, "R", "B"), ("A", 1): (1, "L", "B"),
       ("B", 0): (1, "L", "A"), ("B", 1): (1, "R", "H")}
BB3_ONES = {("A", 0): (1, "R", "B"), ("A", 1): (1, "R", "H"), ("B", 0): (0, "R", "C"),
            ("B", 1): (1, "R", "B"), ("C", 0): (1, "L", "C"), ("C", 1): (1, "L", "A")}
BB3_STEPS = {("A", 0): (1, "R", "B"), ("A", 1): (1, "R", "H"), ("B", 0): (1, "L", "B"),
             ("B", 1): (0, "R", "C"), ("C", 0): (1, "L", "C"), ("C", 1): (1, "L", "A")}


@pytest.mark.parametrize("table, n_ones, n_moves", [
    (BB2, 4, 6),          # Sigma(2) = 4, S(2) = 6 (Rado 1962)
    (BB3_ONES, 6, 14),    # Sigma(3) = 6 (Lin & Rado 1965)
    (BB3_STEPS, 5, 21),   # S(3) = 21 (Lin & Rado 1965)
])
def test_busy_beavers_run_as_laing_machines(table, n_ones, n_moves):
    r = run(turing_program(table), "0", max_steps=10_000)
    assert r.halted
    assert ones(r) == n_ones
    assert r.counts["L"] + r.counts["R"] == n_moves      # one slide per Turing step


def test_increment_machine_counts_in_binary():
    program = parse_program(INCREMENT)
    tape = "0"
    for n in range(1, 40):
        (tape,) = run(program, tape).pieces
        assert int(tape[::-1], 2) == n


# --- networks ---------------------------------------------------------------------
def test_default_closure_is_a_4_bit_counter():
    net = generate_network("laing-molecular-machines")
    tapes = [s.id for s in net.species if s.id.startswith("t:")]
    assert len(tapes) == 16 and {int(t[2:][::-1], 2) for t in tapes} == set(range(16))
    assert len(net.reactions) == 15 and net.status == "truncated"        # 1111 would become 00001
    m = "m:" + INCREMENT
    for r in net.reactions:
        assert r.catalysts == {m: 1} and r.rate is None
    assert "m:" + INCREMENT in {s.id for s in net.species if s.structure.startswith("active")}


def test_binding_all_lists_one_reaction_per_site():
    net = generate_network("laing-molecular-machines", machines=["W1"], tapes=["000"], max_length=3,
                           binding="all")
    from_seed = [r for r in net.reactions if r.reactants == {"m:W1": 1, "t:000": 1}]
    assert sorted(next(iter(Counter(r.products) - Counter({"m:W1": 1}))) for r in from_seed) == \
        ["t:001", "t:010", "t:100"]
    assert net.status == "complete"


def test_detach_gives_several_products():
    net = generate_network("laing-molecular-machines", machines=["R.D"], tapes=["10"], max_length=2)
    # R then D cuts 10 into 1 and 0; on 0 and 1 the slide first recruits a 0 at the right end
    assert sorted(r.to_text() for r in net.reactions) == [
        "m:R.D + t:0 -> m:R.D + 2 t:0",
        "m:R.D + t:1 -> m:R.D + t:0 + t:1",
        "m:R.D + t:10 -> m:R.D + t:0 + t:1",
    ]
    assert net.status == "complete"


def test_soup_observes_counter_reactions():
    net = evolve("laing-molecular-machines", seed=2, copies=20, steps=400, max_length=8).network
    assert net.status == "observed" and net.outflow == "constant-total"
    assert net.initial_state == {"m:" + INCREMENT: 20, "t:0": 20}
    assert sum(net.extras["final_state"].values()) == 40
    assert all(r.count >= 1 and r.catalysts == {"m:" + INCREMENT: 1} for r in net.reactions)


def test_bad_parameters():
    with pytest.raises(ValueError, match="0 and 1"):
        generate_network("laing-molecular-machines", tapes=["012"])
    with pytest.raises(ValueError, match="CT<label>"):
        generate_network("laing-molecular-machines", machines=["CT7.H"])
    with pytest.raises(ValueError, match="binding"):
        evolve("laing-molecular-machines", binding="all")
    with pytest.raises(ValueError, match="max_length"):
        generate_network("laing-molecular-machines", tapes=["00000"])
