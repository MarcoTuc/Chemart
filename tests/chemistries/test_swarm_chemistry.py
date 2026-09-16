"""Swarm Chemistry: Sayama's update equations, his published recipes and their behaviour.

The reference step below is a literal transcription of
SwarmPopulationSimulator.simulateSwarmBehavior + SwarmIndividual.accelerate
(SwarmChemistry-1.3.0-src, https://bingweb.binghamton.edu/~sayama/SwarmChemistry/).
"""

import math

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import swarm_chemistry as S

ID = "swarm-chemistry"


def draws(seed, n, dim):
    """The three random arrays S.step consumes, in the order it consumes them."""
    rng = np.random.default_rng(seed)
    stray = (rng.random((n, dim)) - 0.5) * 1.0        # [-0.5, 0.5)
    whim = (rng.random((n, dim)) - 0.5) * 10.0        # [-5, 5)
    test = rng.random(n)
    return stray, whim, test


def reference_step(pos, vel, par, stray, whim, test):
    """One step of the Java simulator, particle by particle."""
    out = np.empty_like(vel)
    for i in range(len(pos)):
        R, v_normal, v_max, c1, c2, c3, c4, c5 = par[i]
        near = [j for j in range(len(pos))
                if j != i and float(np.sum((pos[j] - pos[i]) ** 2)) < R * R]
        if not near:
            acc = stray[i]
        else:
            acc = (c1 * (pos[near].mean(0) - pos[i])
                   + c2 * (vel[near].mean(0) - vel[i]))
            for j in near:
                d = float(np.sum((pos[i] - pos[j]) ** 2)) or 0.001
                acc = acc + c3 * (pos[i] - pos[j]) / d
            if test[i] < c4:
                acc = acc + whim[i]
        v = vel[i] + acc
        speed = float(np.linalg.norm(v))
        if speed > v_max:
            v = v * (v_max / speed)
        speed = float(np.linalg.norm(v)) or 0.001
        v = v + v * ((v_normal - speed) / speed * c5)
        speed = float(np.linalg.norm(v))
        if speed > v_max:
            v = v * (v_max / speed)
        out[i] = v
    return pos + out, out


# --- the update law --------------------------------------------------------
def test_step_reproduces_the_java_simulator_term_by_term():
    seed, n, dim = 11, 24, 2
    rng = np.random.default_rng(3)
    par = np.column_stack([
        rng.uniform(0, 300, n), rng.uniform(0, 20, n), rng.uniform(0, 40, n),
        rng.uniform(0, 1, n), rng.uniform(0, 1, n), rng.uniform(0, 100, n),
        rng.uniform(0, 0.5, n), rng.uniform(0, 1, n),
    ])
    pos = rng.uniform(0, 200, (n, dim))
    vel = rng.uniform(-5, 5, (n, dim))

    ours_pos, ours_vel = S.step(pos, vel, par, np.random.default_rng(seed))
    want_pos, want_vel = reference_step(pos, vel, par, *draws(seed, n, dim))
    assert ours_vel == pytest.approx(want_vel, abs=1e-12)
    assert ours_pos == pytest.approx(want_pos, abs=1e-12)


def test_three_terms_of_the_acceleration_by_hand():
    # two particles, 10 px apart: cohesion c1 (x_c - x), alignment c2 (v_c - v),
    # separation c3 (x_i - x_j)/|x_i - x_j|^2, then pace keeping toward V_normal.
    par = np.array([[50.0, 3.0, 100.0, 0.5, 0.25, 20.0, 0.0, 0.5]] * 2)
    pos = np.array([[0.0, 0.0], [10.0, 0.0]])
    vel = np.array([[1.0, 0.0], [0.0, 2.0]])
    acc = np.array([0.5 * 10 + 0.25 * (0.0 - 1.0) - 20 * 10 / 100.0,   # cohesion, align, separate
                    0.25 * (2.0 - 0.0)])                               # = (2.75, 0.5)
    v = vel[0] + acc
    speed = math.hypot(*v)
    expected = v + v * ((3.0 - speed) / speed * 0.5)

    new_pos, new_vel = S.step(pos, vel, par, np.random.default_rng(0))
    assert new_vel[0] == pytest.approx(expected, abs=1e-12)
    assert new_pos[0] == pytest.approx(pos[0] + expected, abs=1e-12)


def test_pace_keeping_equals_erskines_published_form():
    # Erskine & Herrmann (ECAL 2013) write v <- c5 (v_n/|v|) v + (1 - c5) v,
    # which is Sayama's v <- v + c5 (v_n - |v|) v/|v|.
    rng = np.random.default_rng(0)
    for _ in range(20):
        v = rng.normal(size=2)
        v_n, c5 = rng.uniform(0.5, 20), rng.uniform(0, 1)
        speed = np.linalg.norm(v)
        assert (v + c5 * (v_n - speed) * v / speed) == pytest.approx(
            c5 * (v_n / speed) * v + (1 - c5) * v, abs=1e-12)


def test_velocity_never_exceeds_v_max_and_lonely_particles_only_stray():
    par = np.array([[0.0, 4.0, 7.0, 1.0, 1.0, 100.0, 0.5, 1.0]] * 5)   # R = 0: no neighbours
    rng = np.random.default_rng(2)
    pos, vel = rng.uniform(0, 100, (5, 2)), rng.uniform(-30, 30, (5, 2))
    for _ in range(20):
        pos, vel = S.step(pos, vel, par, rng)
        assert np.all(np.linalg.norm(vel, axis=1) <= 7.0 + 1e-12)
    assert np.linalg.norm(vel, axis=1) == pytest.approx(4.0, abs=0.6)


def test_a_single_recipe_relaxes_to_its_normal_speed():
    # pace keeping (book 11.4.2, "approximate normal velocity of own type")
    net = generate_network(ID, seed=4, recipe="custom", particles=30, steps=200,
                           custom_recipe=[[30, 0.0, 12.0, 40.0, 0.5, 0.5, 10.0, 0.0, 0.7]])
    analysis = net.extras["analysis"]
    assert analysis["normal_speed"] == [12.0]
    assert analysis["speed"][0] == pytest.approx(12.0, rel=0.02)
    assert analysis["speed_ratio"] == pytest.approx(1.0, abs=0.02)


# --- the published recipes -------------------------------------------------
def test_recipes_are_within_the_ranges_of_book_table_11_6():
    assert S.BOUNDS == {"R": (0.0, 300.0), "V_normal": (0.0, 20.0), "V_max": (0.0, 40.0),
                        "c1": (0.0, 1.0), "c2": (0.0, 1.0), "c3": (0.0, 100.0),
                        "c4": (0.0, 0.5), "c5": (0.0, 1.0)}
    for name, rows in S.RECIPES.items():
        for row in rows:
            assert row[0] >= 1, name
            for field, value in zip(S.FIELDS, row[1:]):
                low, high = S.BOUNDS[field]
                assert low <= value <= high, f"{name}: {field}={value}"


def test_published_recipe_values():
    # Sayama's sample recipes, verbatim from his Swarm Chemistry page.
    assert S.RECIPES["blobs"] == [(300, 20.8, 1.95, 20.75, 0.95, 0.99, 9.31, 0.05, 0.68)]
    assert S.RECIPES["linear-oscillator"] == [
        (133, 214.41, 17.93, 35.14, 0.64, 0.13, 0.29, 0.08, 0.97),
        (24, 253.6, 7.19, 15.51, 0.82, 0.33, 32.65, 0.34, 0.56),
    ]
    assert S.RECIPES["rotary"][1] == (51, 299.13, 0.79, 38.71, 0.25, 0.18, 86.49, 0.38, 0.43)
    # Erskine & Herrmann (ECAL 2013) table 1, yellow then red, mixed 300:50.
    assert S.RECIPES["cell-division"] == [
        (300, 20.5, 1.94, 20.7, 1.0, 1.0, 18.6, 0.05, 1.0),
        (50, 300.0, 15.58, 37.08, 1.0, 0.05, 9.11, 0.47, 0.61),
    ]


def test_species_carry_the_recipe_tuple_and_the_counts_are_rescaled():
    net = generate_network(ID, seed=1, recipe="swinger", particles=60)
    assert [s.id for s in net.species] == [f"swinger-{i}" for i in (1, 2, 3, 4)]
    assert net.species[0].structure.endswith("(150.39, 15.89, 23.54, 0.74, 0.45, 62.65, 0.33, 0.13)")
    counts = [net.initial_state[s.id] for s in net.species]
    assert sum(counts) == 60 and min(counts) >= 1
    assert counts[1] == max(counts)          # 48 : 152 : 14 : 31 in the published recipe


def test_custom_recipes_are_checked_against_the_published_ranges():
    with pytest.raises(ValueError, match="custom_recipe"):
        generate_network(ID, recipe="custom", steps=1)
    with pytest.raises(ValueError, match="table 11.6"):
        generate_network(ID, recipe="custom", steps=1,
                         custom_recipe=[[5, 400.0, 1.0, 2.0, 0.5, 0.5, 1.0, 0.1, 0.5]])
    with pytest.raises(ValueError, match="9 numbers"):
        generate_network(ID, recipe="custom", steps=1, custom_recipe=[[5, 20.0, 1.0]])
    with pytest.raises(ValueError, match="particles"):
        generate_network(ID, recipe="swinger", particles=3, steps=1)


# --- what the network exports ---------------------------------------------
def test_the_default_network_has_no_reactions_only_an_interaction_law():
    net = generate_network(ID, seed=1)
    assert net.reactions == [] and net.status == "complete"
    assert "no reactions; see extras['interaction_law']" in net.summary()
    law = net.extras["interaction_law"]
    assert len(law["steps"]) == 7
    assert "c3_i sum_{j in N_i} (x_i - x_j) / |x_i - x_j|^2" in law["steps"][2]
    assert law["units"]["c3"] == "pixels^2/step^2"
    assert set(law["parameters"]) == {s.id for s in net.species}
    assert law["parameters"]["turbulent-runner-1"]["R"] == 177.1

    space = net.extras["space"]
    assert space["dimensions"] == 2 and len(space["particles"]) == 120
    assert space["boundary"].startswith("unbounded")
    assert sum(net.initial_state.values()) == 120


def test_same_seed_same_swarm_and_different_seeds_differ():
    a = generate_network(ID, seed=3, steps=50)
    b = generate_network(ID, seed=3, steps=50)
    c = generate_network(ID, seed=4, steps=50)
    assert a.to_dict() == b.to_dict()
    assert a.extras["space"]["particles"] != c.extras["space"]["particles"]


# --- the evolutionary variant ---------------------------------------------
def test_recipe_transmission_exports_the_observed_adoptions():
    net = generate_network(ID, seed=1, recipe_transmission=True, steps=200)
    assert net.status == "observed" and net.reactions
    ids = {s.id for s in net.species}
    for r in net.reactions:
        (winner, two), = [(s, n) for s, n in r.products.items()]
        assert two == 2 and set(r.reactants) <= ids and r.reactants[winner] == 1
        assert r.count >= 1 and r.rate is None
    assert sum(net.extras["analysis"]["final_counts"]) == 120
    assert "A + B -> 2 B" in net.extras["interaction_law"]["transmission"]


def test_the_default_box_keeps_sayamas_particle_density():
    # his simulator scatters up to 300 particles in a 300 px box
    domain = generate_network(ID, seed=1, steps=0).extras["space"]["initial_domain"]
    assert domain[0][1] == pytest.approx(300.0 * (120 / 300) ** 0.5, abs=1e-3)
    fixed = generate_network(ID, seed=1, steps=0, space=250.0).extras["space"]["initial_domain"]
    assert fixed == [[0.0, 250.0], [0.0, 250.0]]


# --- the published behaviour ----------------------------------------------
def test_turbulent_runner_segregates_the_two_types():
    for seed in (1, 2):
        a = generate_network(ID, seed=seed, recipe="turbulent-runner").extras["analysis"]
        assert a["segregation_index"] > 0.9      # each particle only sees its own kind
        assert a["cluster_count"] == 2
        assert a["polarization"] > 0.9           # both groups run in a straight line


def test_rotary_turns_as_a_collective():
    for seed in (1, 2):
        a = generate_network(ID, seed=seed, recipe="rotary").extras["analysis"]
        assert 0.15 < abs(a["rotation"]) < 0.5
        assert a["polarization"] > 0.7
    straight = generate_network(ID, seed=1, recipe="turbulent-runner").extras["analysis"]
    assert abs(straight["rotation"]) < 0.05      # a drifting swarm does not turn


def test_blobs_aggregates_and_keeps_its_normal_speed():
    a = generate_network(ID, seed=2, recipe="blobs").extras["analysis"]
    assert 3 <= a["cluster_count"] <= 20                     # blobs: neither a gas nor one ball
    assert a["largest_cluster_fraction"] >= 0.15
    assert a["speed_ratio"] == pytest.approx(1.0, abs=0.15)   # pace keeping


@pytest.mark.slow
def test_erskines_two_species_separate_into_a_core_and_a_divided_shell():
    # ECAL 2013: in 2D "the red particles travel to the inside of a yellow circle
    # of particles causing an inside out division to occur".
    for seed in (0, 1):
        a = generate_network(ID, seed=seed, recipe="cell-division",
                             particles=110, steps=800).extras["analysis"]
        yellow, red = a["mean_radius"]
        assert red < 0.5 * yellow
        assert a["clusters_per_species"][1] == 1             # the reds form one blob
        assert a["clusters_per_species"][0] >= 2             # the yellows have split
        assert a["segregation_index"] > 0.5


@pytest.mark.slow
def test_a_three_dimensional_swarm_behaves_like_the_two_dimensional_one():
    net = generate_network(ID, seed=1, recipe="jelly-fish", dimensions=3, steps=250)
    space = net.extras["space"]
    assert space["dimensions"] == 3 and len(space["particles"][0]["position"]) == 3
    assert net.extras["analysis"]["segregation_index"] > 0.8


def test_the_faster_particle_wins_a_collision():
    pos = np.array([[0.0, 0.0], [1.0, 0.0]])
    vel = np.array([[1.0, 0.0], [5.0, 0.0]])          # particle 1 is faster
    types = np.array([0, 1])
    rows = np.array([[50.0, 3.0, 40.0, 0.5, 0.5, 10.0, 0.0, 0.5]] * 2)
    new_types, events = S._transmit(pos, vel, types, rows, "faster")
    assert list(new_types) == [1, 1] and events == {(0, 1): 1}
    # with "slower" the winner is the other one
    new_types, events = S._transmit(pos, vel, types, rows, "slower")
    assert list(new_types) == [0, 0] and events == {(1, 0): 1}
