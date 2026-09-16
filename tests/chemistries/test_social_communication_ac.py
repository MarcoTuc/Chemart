"""Social communication AC (book 20.4; Dittrich, Kron & Banzhaf, JASSS 6(1):3, 2003).

Reproduces the paper's worked example (Sec. 2.13-2.27) number by number, and the
published order measures of figs. 3, 5, 7, 8, 9 and 11.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.social_communication_ac import Agent, certainty, memorize

EXAMPLE = dict(N=2, alpha=0.5, gamma=2.0, c_f=0.02, ee_memory="ego")   # Sec. 2.13
R_LEARN, R_FORGET = 0.1, 0.01


def analysis(**kwargs):
    return generate_network("social-communication-ac", **kwargs).extras["analysis"]


def mean_analysis(seeds, **kwargs):
    runs = [analysis(seed=s, **kwargs) for s in seeds]
    return {k: float(np.mean([r[k] for r in runs])) for k in ("O_P", "O_AV", "different_activities")}


def selection(agent, received):
    """The published decision table of one agent: EE, EC, f, w_AV, w_AP."""
    p = SimpleNamespace(**EXAMPLE)
    ee = (agent.alter if p.ee_memory == "alter" else agent.ego)[received]
    ec = np.array([certainty(row) for row in agent.alter])
    f = (1 - p.alpha) * ee + p.alpha * ec + p.c_f / p.N                 # eq. 2
    w_av, w_ap = agent.behaviour(p)
    return ee, ec, f, w_av[received], w_ap[received]


# --- the worked example, Sec. 2.13-2.27 --------------------------------------
def test_worked_example_reproduces_the_published_decision_tables():
    """The three decision tables printed for the first three steps of a dyadic run."""
    A, B = Agent(2, 0), Agent(2, 0)                      # both signs initialised to activity 1

    ee, ec, f, w_av, w_ap = selection(A, 0)              # Sec. 2.17-2.20: A reacts to activity 1
    assert ee == pytest.approx([0.5, 0.5])
    assert ec == pytest.approx([0.0, 0.0])
    assert f == pytest.approx([0.26, 0.26])
    assert w_av == pytest.approx([0.5, 0.5])
    assert w_ap == pytest.approx([0.5, 0.5])

    # Sec. 2.21: A selects activity 1; the pair (1, 1) goes to A's ego- and B's alter-memory.
    memorize(A.ego, 0, 0, R_LEARN, R_FORGET)
    memorize(B.alter, 0, 0, R_LEARN, R_FORGET)

    ee, ec, f, w_av, w_ap = selection(B, 0)              # Sec. 2.22-2.24: B reacts to activity 1
    assert ee == pytest.approx([0.5, 0.5])
    assert ec == pytest.approx([0.005863, 0.0], abs=5e-7)
    assert f == pytest.approx([0.262931, 0.26], abs=5e-7)
    assert w_av == pytest.approx([0.502803, 0.497197], abs=5e-7)
    assert w_ap == pytest.approx([0.505605, 0.494395], abs=5e-7)      # eq. 12

    memorize(B.ego, 0, 0, R_LEARN, R_FORGET)
    memorize(A.alter, 0, 0, R_LEARN, R_FORGET)

    ee, ec, f, w_av, w_ap = selection(A, 0)              # Sec. 2.26-2.27: A reacts again
    assert ee == pytest.approx([0.545045, 0.454955], abs=5e-7)
    assert ec == pytest.approx([0.005863, 0.0], abs=5e-7)
    assert f == pytest.approx([0.285454, 0.237477], abs=5e-7)
    assert w_av == pytest.approx([0.545872, 0.454128], abs=5e-7)
    assert w_ap == pytest.approx([0.590979, 0.409021], abs=5e-7)      # "about 59% / 41%"


def test_memorize_reproduces_the_published_memory_matrices():
    """Sec. 2.21: memorising (1, 1) turns the row (0.5, 0.5) into (0.545045, 0.454955)."""
    m = np.full((2, 2), 0.5)
    memorize(m, 0, 0, R_LEARN, R_FORGET)
    assert m[0] == pytest.approx([0.545045, 0.454955], abs=5e-7)
    assert m[1] == pytest.approx([0.5, 0.5]), "untouched rows stay uniform"
    assert m.sum(axis=1) == pytest.approx([1.0, 1.0]), "eq. 6: every row is a distribution"


def test_certainty_is_the_entropy_measure_of_eq_1():
    assert certainty([0.25] * 4) == 0.0                  # Sec. 2.5: no information
    assert certainty([1.0, 0.0, 0.0, 0.0]) == 1.0        # maximal certainty
    assert certainty([0.545045, 0.454955]) == pytest.approx(0.005863, abs=5e-7)   # Sec. 2.23


# --- the network ---------------------------------------------------------------
def test_default_network_is_the_observed_communication_events():
    net = generate_network("social-communication-ac", seed=1)
    assert net.status == "observed"
    assert [s.id for s in net.species] == [f"a{i}" for i in range(1, 11)]
    assert sum(r.count for r in net.reactions) == 1000, "every step is one communication event"
    assert sum(net.initial_state.values()) == 10, "one displayed sign per agent"
    for r in net.reactions:
        assert r.catalysts, "the message reacted to survives the event"
        assert sum(r.reactants.values()) == sum(r.products.values()) == 2
        assert r.rate is None
    assert len(net.extras["agents"]) == 10
    assert np.allclose([np.sum(a["ego_memory"], axis=1) for a in net.extras["agents"]], 1.0)


def test_same_seed_same_run():
    kwargs = dict(n_agents=6, steps=300, seed=11)
    assert generate_network("social-communication-ac", **kwargs).to_dict() == \
        generate_network("social-communication-ac", **kwargs).to_dict()


@pytest.mark.parametrize("given, message", [
    ({"interaction": "dyadic", "n_agents": 5}, "dyadic"),
    ({"n_agents": 4, "n_observers": 3}, "n_observers"),
])
def test_invalid_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        generate_network("social-communication-ac", **given)


# --- the published results -------------------------------------------------------
def test_expectation_certainty_destroys_systems_level_order():
    """Sec. 6.16 / fig. 9: with observers, order appears only for alpha = 0 (EE only)."""
    seeds = (1, 2, 3)
    assert mean_analysis(seeds, alpha=0.0)["O_P"] > 0.95
    for alpha in (0.5, 1.0):                             # published O_P at M = 10: ~0.26, ~0.30
        assert mean_analysis(seeds, alpha=alpha)["O_P"] < 0.6


@pytest.mark.slow
@pytest.mark.parametrize("M", [5, 10, 20, 30])
def test_order_is_scalable_with_alter_memory_and_observers(M):
    """Fig. 9: O_P is maximal and about three activities are used, for every M."""
    got = mean_analysis((1, 2), n_agents=M, steps=2000)
    assert got["O_P"] > 0.95                             # published: O_P = 1 from M = 2 to 30
    assert got["O_AV"] > 0.8                             # published: about 0.95, flat in M
    assert got["different_activities"] < 6               # published: about 3, flat in M


@pytest.mark.slow
def test_order_does_not_scale_without_observers():
    """Fig. 8: alter-memory alone gives maximal O_P only in small populations."""
    got = {M: mean_analysis((1, 2), n_agents=M, n_observers=0, steps=2000)["O_P"]
           for M in (5, 10, 30)}
    assert got[5] > 0.9                                  # published: O_P = 1 up to M = 6
    assert got[5] > got[10] > got[30]                    # published: 1.0, 0.93, 0.21
    assert got[30] < 0.5, "order has disappeared in the large population"


@pytest.mark.slow
@pytest.mark.parametrize("M, expected", [(10, 0.23), (20, 0.11)])
def test_ego_memory_gives_certain_agents_but_no_systems_level_order(M, expected):
    """Fig. 7 and Sec. 7.10: every agent is sure what to do, but all do something else."""
    got = mean_analysis((1, 2), n_agents=M, n_observers=0, ee_memory="ego", steps=4000)
    assert got["O_P"] == pytest.approx(expected, abs=0.08)
    assert got["O_AV"] > 0.7                             # published: 0.90 (M = 10), 0.83 (M = 20)
    assert got["different_activities"] > 9               # published: about 9.8 of 10


@pytest.mark.slow
def test_an_activity_system_emerges_and_is_shared():
    """Fig. 11: a small closed, self-maintaining set of activities, e.g. {1, 5, 6} of 10."""
    got = analysis(seed=1, n_observers=0, r_learn=0.2, steps=4000, edge_threshold=0.1)
    systems = got["activity_systems"]
    assert len(systems) == 1, "one common activity pattern, shared by all agents"
    assert 1 <= len(systems[0]) <= 6                     # published: 3 activities of 10
    assert set(systems[0]) <= set(got["activities_used"])
    assert got["O_P"] > 0.9, "transitions within the system are practically deterministic"


@pytest.mark.slow
def test_dyadic_order_matches_figure_5():
    """Sec. 5.14: at N = 20, alpha = 0, gamma = 1 agents use about 15 activities at O_AV ~ 0.58."""
    dyad = dict(N=20, n_agents=2, interaction="dyadic", ee_memory="ego", n_observers=0,
                r_learn=0.2, alpha=0.0)
    proportional = mean_analysis((1, 2, 3, 4, 5), gamma=1.0, **dyad)
    assert proportional["different_activities"] == pytest.approx(15, abs=3)
    assert proportional["O_AV"] == pytest.approx(0.58, abs=0.13)
    quadratic = mean_analysis((1, 2, 3), gamma=2.0, **dyad)
    assert quadratic["different_activities"] < proportional["different_activities"] / 2
    assert quadratic["O_AV"] > 0.9, "less randomness, much more order"


@pytest.mark.slow
def test_dyadic_order_survives_many_possible_messages():
    """Fig. 3: at gamma = 1, alpha = 0.5 the number of activities used saturates near 30."""
    got = mean_analysis((1, 2), N=64, n_agents=2, interaction="dyadic", ee_memory="ego",
                        n_observers=0, gamma=1.0, r_learn=0.2, alpha=0.5)
    assert got["different_activities"] == pytest.approx(27, abs=6)
    assert got["different_activities"] < 64 / 2, "much fewer than random choice would give"
