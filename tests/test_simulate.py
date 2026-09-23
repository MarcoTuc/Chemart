"""chemart.simulate: rate equations, Gillespie's SSA, and assigning rates and states."""

import json

import numpy as np
import pytest

import chemart
from chemart import simulate
from chemart.helpers.explicit import network
from chemart.network import CONSTANT_TOTAL
from chemart.simulate import NotSimulable
from chemart.trajectory import Trajectory


def birth_death(k_in=5.0, k_out=0.5, x=0.0):
    return network([(" -> X", k_in), ("X -> ", k_out)], initial_state={"X": x} if x else None)


# --- deterministic ------------------------------------------------------------
def test_brusselator_limit_cycle():
    # b = 3 > 1 + a^2 = 2: the fixed point is unstable and X settles on a limit cycle
    net = chemart.generate_network("brusselator", seed=1)
    traj = simulate.ode(net, t_end=60, points=600)
    ids, t, X = traj.array(["X"])
    late = X[t > 30, 0]
    assert late.max() - late.min() > 1.0
    assert traj.method == "ode" and traj.frames[0].t == 0.0 and traj.frames[-1].t == 60


def test_birth_death_reaches_k_in_over_k_out():
    traj = simulate.ode(birth_death(), t_end=40, x0=0.0)
    assert traj.frames[-1].state["X"] == pytest.approx(10.0, rel=1e-4)


def test_jacobian_matches_finite_differences():
    net = chemart.generate_network("brusselator", seed=1)
    ids, f = simulate.rhs(net)
    J = simulate.jacobian(net)
    x = np.array([1.3, 0.7, 2.0, 1.1, 0.4, 0.9])[:len(ids)]
    eps = 1e-6
    fd = np.array([(f(0, x + eps * e) - f(0, x - eps * e)) / (2 * eps) for e in np.eye(len(ids))]).T
    assert np.allclose(J(0, x), fd, atol=1e-6)


# --- stochastic ---------------------------------------------------------------
def test_ssa_mean_matches_the_ode_steady_state():
    # birth-death: the stationary distribution is Poisson with mean k_in V / k_out
    traj = simulate.ssa(birth_death(), t_end=400, volume=20, seed=0, points=801, x0=0.0)
    ids, t, X = traj.array(["X"])
    assert X[t > 50, 0].mean() == pytest.approx(10.0, rel=0.05)
    assert traj.method == "ssa" and traj.settings["events"] > 1000


def test_michaelis_menten_ssa_agrees_with_ode_at_large_volume():
    net = chemart.generate_network("michaelis-menten", seed=1)
    ode = simulate.ode(net, t_end=5, points=11)
    ssa = simulate.ssa(net, t_end=5, points=11, volume=2000, seed=1)
    for s in ("S", "P"):
        a = ode.frames[-1].state.get(s, 0.0)
        b = ssa.frames[-1].state.get(s, 0.0)
        assert b == pytest.approx(a, rel=0.05, abs=0.02), s


def test_ssa_is_reproducible_and_counts_fired_reactions():
    a = simulate.ssa(birth_death(), t_end=10, seed=3, x0=0.0)
    b = simulate.ssa(birth_death(), t_end=10, seed=3, x0=0.0)
    assert a == b
    fired = sum(n for f in a.frames for *_, n in f.fired)
    assert fired == a.settings["events"]


def test_ssa_constant_total_dilution_keeps_the_population():
    net = network([("X -> 2 X", 1.0), ("Y -> 2 Y", 0.5)], initial_state={"X": 50.0, "Y": 50.0},
                  outflow=CONSTANT_TOTAL)
    traj = simulate.ssa(net, t_end=15, seed=0, points=31)
    assert all(sum(f.state.values()) == 100 for f in traj.frames)
    assert traj.frames[-1].state.get("X", 0) > 80           # the faster replicator wins


def test_ssa_holds_buffered_species_and_applies_flows():
    net = network([("A + X -> 2 X", 0.01)], initial_state={"A": 100.0, "X": 5.0},
                  outflow={"X": 0.2}, extras={"buffered": ["A"]})
    traj = simulate.ssa(net, t_end=5, seed=0)
    assert all(f.state["A"] == 100.0 for f in traj.frames)
    # growth 0.01 * 100 = 1 against washout 0.2: X grows
    assert traj.frames[-1].state["X"] > 50


def test_max_events_stops_and_says_so():
    traj = simulate.ssa(birth_death(), t_end=1000, seed=0, max_events=50, x0=0.0)
    assert traj.settings["events"] == 50 and "max_events" in traj.settings["stopped"]


# --- assigning rates and states ----------------------------------------------------
def test_a_solution_that_escapes_stops_and_says_so():
    # x' = x^2 runs away at t = 1/x0; the solver would otherwise shrink its step forever
    net = network([("2 X -> 3 X", {"law": "mass-action", "k": 1.0})], initial_state={"X": 1.0})
    traj = simulate.ode(net, 10, points=11)
    assert "blew up" in traj.settings["stopped"]
    assert traj.frames[-1].t < 10


def test_unrated_networks_need_rates():
    net = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
    with pytest.raises(NotSimulable, match="no rate constants"):
        simulate.ode(net, t_end=1, x0=1.0)
    traj = simulate.ode(net, t_end=1, x0=1.0, rates={"dist": "lognormal", "mean": 0, "sigma": 1}, seed=0)
    assert traj.network.provides.count("rate-constants") == 1


def test_assign_is_deterministic_per_seed_and_takes_tables():
    net = network([("A -> B", None), ("B -> C", 2.0)])
    draw = {"dist": "uniform", "low": 1, "high": 2}
    a = simulate.assign(net, rates=draw, rng=np.random.default_rng(1))
    b = simulate.assign(net, rates=draw, rng=np.random.default_rng(1))
    assert a.reactions[0].rate == b.reactions[0].rate and 1 <= a.reactions[0].rate["k"] <= 2
    filled = simulate.assign(net, rates=7.0, fill_only=True)
    assert [r.rate["k"] for r in filled.reactions] == [7.0, 2.0]
    table = simulate.assign(net, rates={"A -> B": 3.0, "*": 0.5}, x0={"A": 1.0})
    assert [r.rate["k"] for r in table.reactions] == [3.0, 0.5]
    assert table.initial_state == {"A": 1.0}
    assert net.reactions[0].rate is None                   # the input is left alone


def test_rates_and_states_from_files(tmp_path):
    net = network([("A -> B", None), ("B -> A", None)])
    (tmp_path / "rates.csv").write_text("reaction,k\nA -> B,2.5\n1,0.5\n")
    (tmp_path / "x0.json").write_text(json.dumps({"A": 3.0}))
    out = simulate.assign(net, rates=tmp_path / "rates.csv", x0=str(tmp_path / "x0.json"))
    assert [r.rate["k"] for r in out.reactions] == [2.5, 0.5]
    assert out.initial_state == {"A": 3.0}
    with pytest.raises(ValueError, match="json or .csv"):
        simulate.assign(net, rates=tmp_path / "rates.txt")


def test_bad_specs_are_explained():
    net = network([("A -> B", None)])
    with pytest.raises(ValueError, match="unknown distribution"):
        simulate.assign(net, rates={"dist": "cauchy"})
    with pytest.raises(ValueError, match="needs"):
        simulate.assign(net, rates={"dist": "uniform", "low": 1})
    with pytest.raises(NotSimulable, match="no initial state"):
        simulate.ode(network([("A -> B", 1.0)]), t_end=1)


def test_arrhenius_needs_a_temperature():
    net = network([("A -> B", {"law": "arrhenius", "A": 1.0, "Ea": 1.0})], initial_state={"A": 1.0})
    with pytest.raises(NotSimulable, match="temperature"):
        simulate.ode(net, t_end=1)
    traj = simulate.ode(net, t_end=1, temperature=1.0, gas_constant=1.0)
    assert traj.frames[-1].state["A"] == pytest.approx(np.exp(-np.exp(-1.0)), rel=1e-5)


# --- the Trajectory record ---------------------------------------------------------------
def test_trajectory_round_trips_and_windows():
    traj = simulate.ssa(birth_death(), t_end=10, seed=0, points=11, x0=0.0)
    again = Trajectory.from_dict(json.loads(json.dumps(traj.to_dict())))
    assert again == traj
    whole = traj.window(-1, width=None)
    assert sum(r.count for r in whole.reactions) == traj.settings["events"]
    assert len(traj.turnover()) == len(traj.frames)
    with pytest.raises(KeyError, match="no observable"):
        traj.series("entropy")


# --- the command line ---------------------------------------------------------------
def test_cli_simulate(capsys):
    from chemart.cli import main

    assert main(["simulate", "brusselator", "--t-end", "5", "--points", "3", "--format", "csv"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0] == "t,A,B,X,Y,D,E" and len(lines) == 4
    assert main(["simulate", "brusselator", "--method", "ssa", "--t-end", "1", "--seed", "1"]) == 0
    assert "brusselator: 6 species" in capsys.readouterr().out
    assert main(["simulate", "kauffman-autocatalytic-sets", "--t-end", "1"]) == 2
    assert "no initial state" in capsys.readouterr().err
