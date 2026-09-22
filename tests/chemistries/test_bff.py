"""BFF reproduces Agüera y Arcas et al. (2024), "Computational Life".

Source: arXiv:2406.19108v2, section 2 (the BFF language, the primordial soup,
high-order entropy, the Fig. 4 replicator, the seeded runs of Fig. 7 and the 2D
grid of section 2.2), checked against the authors' cubff code.
"""

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries.bff import (
    FIG4_REPLICATOR, TAPE, encode, high_order_entropy, run, show,
)

R = encode(FIG4_REPLICATOR)


def execute(tape: bytes, steps: int = 8192) -> bytearray:
    out = bytearray(tape)
    run(out, steps)
    return out


# --- the language ------------------------------------------------------------------
def test_heads_copy_and_arithmetic():
    # '+' on head0 (at 0), then '.' copies tape[head0] to tape[head1]; '}' moves head1
    tape = execute(encode("+}}}}}.") + bytes(1))
    assert tape[0] == ord("+") + 1 and tape[5] == ord("+") + 1


def test_arithmetic_wraps_and_loops_skip_on_zero():
    tape = execute(bytes([0]) + encode("-"))   # head0 at the zero byte: no-op, then 0 - 1
    assert tape[0] == 255
    # head0 is on the zero byte, so '[' jumps past the matching ']' ('+' is never run)
    tape = execute(bytes(1) + encode("[+]"))
    assert tape == bytearray(bytes(1) + encode("[+]"))


def test_unmatched_bracket_ends_the_program():
    # ']' with a non-zero byte under head0 and no '[' before it: execution stops at once
    tape = execute(encode("]+"))
    assert tape == bytearray(encode("]+"))


def test_heads_wrap_around_the_tape():
    # '{' moves head1 from 0 to the last byte; '.' copies tape[0] there
    tape = execute(encode("{.") + bytes(6))
    assert tape[-1] == ord("{")


def test_text_round_trip():
    tape = bytes(range(256))[:TAPE]
    assert encode(show(tape)) == tape
    assert show(encode("[0]")) == "[0]"


# --- the Fig. 4 replicator -------------------------------------------------------------
@pytest.mark.parametrize("row, suffix", [
    (4, "0000"),        # rows 1-4: nothing written yet
    (5, "000["),        # row 5: the first byte lands at the end of tape B
    (9, "00[["),
    (13, "0{[["),
])
def test_fig4_trace_end_of_tape_b(row, suffix):
    # Fig. 4 row k shows the tape after k - 1 steps; the write head moves left from the end
    tape = execute(R + bytes(TAPE), steps=row - 1)
    assert show(tape[-4:]) == suffix


def test_fig4_trace_first_byte_of_tape_b():
    # rows 255-256 show B starting '0[{.>]-]', row 257 '[[{.>]-]': the copy is complete
    assert show(execute(R + bytes(TAPE), steps=254)[TAPE:TAPE + 8]) == "0[{.>]-]"
    assert show(execute(R + bytes(TAPE), steps=256)[TAPE:TAPE + 8]) == "[[{.>]-]"


def test_fig4_replicator_copies_itself_exactly():
    # a palindrome copied in reverse is an exact copy; A is left intact
    tape = execute(R + bytes(TAPE))
    assert tape[:TAPE] == R and tape[TAPE:] == R


def test_replicator_turns_any_food_into_itself():
    # eq. 5: S + F -> split(exec(SF)) = 2 S, whatever F is
    rng = np.random.default_rng(0)
    for _ in range(20):
        food = rng.integers(0, 256, size=TAPE, dtype=np.uint8).tobytes()
        tape = execute(R + food)
        assert tape[:TAPE] == R and tape[TAPE:] == R


# --- high-order entropy -------------------------------------------------------------
def test_high_order_entropy_is_zero_for_noise_and_high_for_copies():
    rng = np.random.default_rng(0)
    noise = rng.integers(0, 256, size=64 * 1024, dtype=np.uint8).tobytes()
    assert abs(high_order_entropy(noise)) < 0.1
    copies = rng.integers(0, 256, size=64, dtype=np.uint8).tobytes() * 1024
    assert high_order_entropy(copies) > 5.0


# --- the soup ------------------------------------------------------------------------
def test_soup_records_reactions_and_mutations():
    net = generate_network("bff", seed=0)
    kinds = net.extras["reaction_kinds"]
    assert len(kinds) == len(net.reactions) and set(kinds) <= {"execution", "mutation"}
    for r, kind in zip(net.reactions, kinds):
        assert sum(r.reactants.values()) == sum(r.products.values()) == (2 if kind == "execution" else 1)
    assert sum(net.extras["final_state"].values()) == 128


def test_a_frame_per_epoch():
    traj = evolve("bff", seed=0, epochs=5)
    assert traj.clock == "epochs" and traj.times() == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    first = traj.frames[0]
    assert first.fired == [] and first.observables["ops_per_run"] == 0.0
    assert all(sum(f.state.values()) == 128 for f in traj.frames)
    assert traj.frames[-1].state == {s: float(n) for s, n in traj.network.extras["final_state"].items()}
    assert set(first.observables) == {"high_order_entropy", "top_tape_count", "ops_per_run", "zero_bytes"}


def test_no_mutation_no_mutation_reactions():
    net = generate_network("bff", seed=0, mutation_rate=0.0, epochs=4)
    assert "mutation" not in net.extras["reaction_kinds"]


def test_seeded_replicator_shows_as_catalyst():
    # S + F -> S + S: the replicator is on both sides, so the network carries catalysts
    net = generate_network("bff", seed=0, tapes=16, epochs=6, replicators=4, mutation_rate=0.0)
    rep = show(R)
    assert any(r.reactants.get(rep) == 1 and r.products.get(rep) == 2 for r in net.reactions)
    assert "catalysts" in net.provides


def test_grid_pairs_only_neighbours():
    # section 2.2: two programs interact only if |x0 - x1| <= 2 and |y0 - y1| <= 2
    from types import SimpleNamespace
    from chemart.chemistries.bff import _neighbours, _pairs
    width = 8
    pairs = _pairs(SimpleNamespace(tapes=64), np.random.default_rng(0), _neighbours(width, 8))
    assert pairs
    for i, j in pairs:
        assert i != j and abs(i % width - j % width) <= 2 and abs(i // width - j // width) <= 2
    assert len({k for pair in pairs for k in pair}) == 2 * len(pairs)
    net = generate_network("bff", seed=0, tapes=64, space="grid", width=8, epochs=2)
    space = net.extras["space"]
    assert (space["width"], space["height"], space["radius"]) == (8, 8, 2)
    assert len(space["final_grid"]) == 64


def test_grid_needs_a_multiple_of_the_width():
    with pytest.raises(ValueError, match="multiple of width"):
        generate_network("bff", tapes=10, space="grid", width=3)


@pytest.mark.slow
def test_seeded_runs_reach_a_state_transition_and_random_short_runs_do_not():
    # Fig. 7: 'seeded' runs (one Fig. 4 replicator, 128 epochs) reach a state transition
    # 22% of the time; 'short' runs from random programs almost never (3 of 1000).
    # A transition is high-order entropy >= 1 (Fig. 6), read at the start and after
    # 100 epochs (every=100 keeps those two frames). Here 64 programs.
    def transitions(**kw):
        out = 0
        for seed in range(12):
            traj = evolve("bff", seed=seed, tapes=64, epochs=100, every=100, **kw)
            out += max(traj.series("high_order_entropy")) >= 1.0
        return out

    seeded = transitions(replicators=1)
    assert 2 <= seeded <= 10
    assert transitions() == 0
