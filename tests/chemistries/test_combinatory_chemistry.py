"""Combinatory Chemistry reproduces Kruszewski & Mikolov (ALIFE 2020).

Source: "Combinatory Chemistry: Towards a Simple Model of Emergent Evolution",
ALIFE 2020, 411-419: eqs. 1-4, Algorithms 1-2 and the metabolic cycles of
Figs. 1-3 and of the (SSK) structure described in the text.
"""

from collections import Counter

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries.combinatory_chemistry import (
    Pool, Reactor, atom_counts, cleave, parse, reduce_at, redexes, show,
)


def reduce_once(text):
    """The single head redex of `text`, reduced: (result, by-products)."""
    seq = parse(text, variables=True)
    (path, _, _), = [r for r in redexes(seq) if r[0] == ()]
    new, by = reduce_at(seq, path)
    return show(new), sorted(show(b) for b in by)


# --- the reduce reactions (eqs. 1-3) ---------------------------------------------
def test_reduce_reactions_release_what_logic_destroys():
    assert reduce_once("Ifu") == ("fu", ["I"])                 # α(If)β -> αfβ + I
    assert reduce_once("Kfgu") == ("fu", ["K", "g"])           # α(Kfg)β -> αfβ + g + K
    assert reduce_once("Sfgxu") == ("fx(gx)u", ["S"])          # + x consumed as reactant


def test_s_redex_names_its_reactant():
    (_, kind, reactant), = redexes(parse("SII(SII)"))
    assert (kind, reactant) == ("S", "SII")


def test_i_sii_i_sii_has_two_redexes():
    # paper: "I(SII)(I(SII)) has two redexes ... the outermost or the innermost I"
    found = redexes(parse("I(SII)(I(SII))"))
    assert [(path, kind) for path, kind, _ in found] == [((), "I"), ((2,), "I")]


def test_sii_is_not_reducible():
    # paper: "(SII) is not reducible because S requires three arguments"
    assert redexes(parse("SII")) == []


def test_parentheses_are_left_associative():
    assert show(parse("((SI)I)")) == "SII"
    assert show(parse("S(I)((K)K)")) == "SI(KK)"


def test_cleave_is_the_top_split():
    x, y = cleave(parse("SII(SII)"))
    assert (show(x), show(y)) == ("SII", "SII")
    assert cleave(parse("S")) is None


def test_s_redex_needs_its_reactant_in_the_multiset():
    # "SII(SII) is reducible if and only if ... (SII) is also present in the set"
    alone = Reactor(Pool(["SII(SII)"]), F=1, max_reductions=100)
    assert alone.reductions("SII(SII)") == []
    fed = Reactor(Pool(["SII(SII)", "SII"]), F=1, max_reductions=100)
    assert len(fed.reductions("SII(SII)")) == 1


def test_reactant_assemblage_builds_small_reactants_from_free_atoms():
    # Algorithm 2: with F = 3, SII is built on the spot from S, I, I
    pool = Pool(["SII(SII)", "S", "I", "I"])
    reactor = Reactor(pool, F=3, max_reductions=100)
    assert reactor.reductions("SII(SII)")
    assert not Reactor(Pool(["SII(SII)", "S", "I", "I"]), F=2, max_reductions=100).reductions("SII(SII)")
    events = reactor.react("SII(SII)", np.random.default_rng(0))
    assert [e[0] for e in events] == ["assemblage", "S"]
    assert pool.count == Counter({"I(SII)(I(SII))": 1, "S": 1})


# --- metabolic cycles (Figs. 1-3) ---------------------------------------------------
def pathway(start, food, target, depth):
    """Is there a chain of reduce reactions from `start`, feeding on `food`, to `target`?

    `target` is (expression, food consumed, other released molecules, as a Counter).
    Reactants other than the food are not available, as in the paper's cycles.
    """
    goal_expr, goal_eaten, goal_released = target
    seen = set()

    def search(expr, eaten, released, left):
        key = (expr, eaten, tuple(sorted(released.items())))
        if key in seen:
            return False
        seen.add(key)
        if expr == goal_expr and eaten == goal_eaten and released == goal_released:
            return True
        if not left or eaten > goal_eaten:
            return False
        seq = parse(expr)
        for path, kind, reactant in redexes(seq):
            if kind == "S" and reactant != food:
                continue
            new, by = reduce_at(seq, path)
            if search(show(new), eaten + (kind == "S"), released + Counter(show(b) for b in by), left - 1):
                return True
        return False

    return search(start, 0, Counter(), depth)


def phi(a):
    """The atoms of a, as released molecules."""
    return Counter(atom_counts(a).elements())


def test_fig2_sii_sii_is_simple_autopoietic():
    # (AA) + A => (AA) + phi(A), A = SII: r1 takes A and releases S, two I-reactions follow
    a = "SII"
    assert pathway(a + f"({a})", a, (a + f"({a})", 1, phi(a)), depth=3)


def test_fig3_tail_recursive_growth():
    # (AA) + 2A => A(AA) + phi(A), A = S(SI)I
    a = "S(SI)I"
    aa = f"{a}({a})"
    assert pathway(aa, a, (f"{a}({aa})", 2, phi(a)), depth=8)


def test_fig1_self_reproduction():
    # (AA) + 3A => 2(AA) + phi(A), A = SI(S(SK)I): the second AA is released by a K-reaction
    a = "SI(S(SK)I)"
    aa = f"{a}({a})"
    assert pathway(aa, a, (aa, 3, phi(a) + Counter({aa: 1})), depth=10)


def test_ssk_cycle_releases_a_catalyst():
    # (AAA) + 2A => (AAA) + A + phi(A), A = SSK
    a = "SSK"
    aaa = f"{a}({a})({a})"
    assert pathway(aaa, a, (aaa, 2, phi(a) + Counter({a: 1})), depth=10)


# --- the reactor (Algorithm 1) ---------------------------------------------------------
def test_every_reaction_conserves_atoms():
    net = generate_network("combinatory-chemistry", seed=3)
    for r in net.reactions:
        assert atoms(r.reactants) == atoms(r.products), r.to_text()
    assert atoms(net.extras["final_state"]) == atoms(net.initial_state)


def atoms(multiset):
    out = Counter()
    for s, n in multiset.items():
        for a, k in atom_counts(s).items():
            out[a] += k * int(n)
    return out


def test_starts_from_atoms_and_diversity_explodes():
    traj = evolve("combinatory-chemistry", seed=0)
    assert set(traj.network.initial_state) == {"S", "K", "I"}
    diversity = [len(f.state) for f in traj.frames]
    assert diversity[0] == 3 and max(diversity) > 50


def test_a_frame_every_record_every_iterations():
    traj = evolve("combinatory-chemistry", seed=0, iterations=2500, record_every=1000)
    assert traj.clock == "iterations" and traj.times() == [0.0, 1000.0, 2000.0, 2500.0]
    first = traj.frames[0]
    assert first.fired == [] and first.state == {"I": 334.0, "K": 333.0, "S": 333.0}
    assert first.observables["free_atoms"] == {"S": 333, "K": 333, "I": 334}
    assert traj.frames[-1].state == {s: float(n) for s, n in traj.network.extras["final_state"].items()}


def test_reaction_kinds_are_recorded():
    net = generate_network("combinatory-chemistry", seed=0)
    kinds = net.extras["reaction_kinds"]
    assert len(kinds) == len(net.reactions)
    assert {"I", "K", "S", "cleave", "condense"} <= set(kinds)


@pytest.mark.slow
def test_paper_scale_diversity_and_sii_autopoiesis():
    # Fig. 4a: 10,000 atoms, diversity explodes to a few hundred within ~200k reactions;
    # Fig. 4c: reductions are 10-35% of reactions;
    # Fig. 5a: SII is consumed above binary reactants because SII(SII) structures emerge.
    traj = evolve("combinatory-chemistry", seed=0, n_I=3334, n_K=3333, n_S=3333,
                  iterations=2_000_000, record_every=500_000)
    assert 200 <= len(traj.frames[1].state) <= 500   # after the first 500,000 iterations
    assert all(0.08 <= r <= 0.35 for r in traj.series("reductions")[1:])
    late = traj.series("top_reactants")[-1]
    binary = max(late.get(x, 0) for x in ("KI", "KK", "II", "SI", "SK", "IK", "IS", "KS", "SS"))
    assert late.get("SII", 0) > binary
    assert traj.network.extras["final_state"].get("SII(SII)", 0) >= 1
