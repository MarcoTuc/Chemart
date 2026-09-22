"""Embedded particles in cellular automata: the computational-mechanics reading.

Published results reproduced here:

- the particle catalogs of the two density-classification CAs - regular
  domains, the wall each particle is, and its velocity - exactly as printed in
  Crutchfield, Mitchell & Das, *The Evolutionary Design of Collective
  Computation in Cellular Automata* (arXiv:adap-org/9809001), Tables 3 and 4,
  together with their interaction tables;
- the GKL look-up table of Crutchfield & Mitchell (PNAS 1995, Table 1) expands
  bit for bit to the GKL rule as *defined* by Gacs, Kurdyumov and Levin
  (Mitchell, Hraber & Crutchfield 1993, Sec. 3);
- those published velocities, measured by running the automaton, and alpha's
  decay into gamma + mu;
- the density-classification performance P_149 of Table 1 (0.775 for phi_par^a,
  0.816 for GKL);
- the condensation time, about 12 steps at N = 149 (Hordijk, Crutchfield &
  Mitchell, PPSN V 1998, Sec. 4);
- every interaction of the published table fires in a real run.
"""

from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries import ca_embedded_particles as C

ID = "ca-embedded-particles"

#: Crutchfield, Mitchell & Das (1998), Table 3: phi_par^a's particle catalog.
TABLE_3 = {
    "domains": {"L0": "0", "L1": "1", "L2": "01"},
    "particles": {
        "alpha": ("L0", "L1", None),      # unstable, no velocity is listed
        "beta": ("L1", "L0", 0.0),
        "gamma": ("L0", "L2", -1.0),
        "delta": ("L2", "L0", -3.0),
        "eta": ("L1", "L2", 3.0),
        "mu": ("L2", "L1", 1.0),
    },
    "interactions": [
        "alpha -> gamma + mu",            # decay
        "beta + gamma -> eta",            # react
        "mu + beta -> delta",
        "eta + delta -> beta",
        "eta + mu -> ∅",                  # annihilate, leaving Lambda1
        "gamma + delta -> ∅",             # annihilate, leaving Lambda0
    ],
}

#: Table 4: phi_par^b's catalog, with the striped domain (011)+.
TABLE_4 = {
    "domains": {"L0": "0", "L1": "1", "L2": "011"},
    "particles": {
        "alpha": ("L1", "L0", 0.0),
        "beta": ("L0", "L1", 1.0),
        "gamma": ("L1", "L2", 0.0),
        "delta": ("L2", "L1", -3.0),
        "eta": ("L0", "L2", 3.0),
        "mu": ("L2", "L0", 1.5),
    },
    "interactions": TABLE_3["interactions"],
}

PUBLISHED = {"phi-par-a": TABLE_3, "phi-par-b": TABLE_4}

#: Table 1: the performance of each rule on the rho_c = 1/2 task at N = 149.
P149 = {"phi-par-a": 0.775, "phi-par-b": 0.766, "gkl": 0.816}


def run(**kw):
    """The observed network of a run of the automaton."""
    return evolve(ID, **kw).network


def multiset(stoich: dict) -> tuple:
    """A reaction side as a sorted tuple of species, e.g. ('gamma', 'delta')."""
    return tuple(sorted(Counter(stoich).elements()))


# --- the published catalogs ----------------------------------------------------------
@pytest.mark.parametrize("rule", sorted(PUBLISHED))
def test_published_particle_catalog_is_reproduced_exactly(rule):
    table, spec = PUBLISHED[rule], C.RULES[rule]
    assert spec["domains"] == table["domains"]
    assert spec["particles"] == table["particles"]
    assert spec["performance"]["149"] == P149[rule]


@pytest.mark.parametrize("rule", sorted(PUBLISHED))
def test_published_interaction_table_is_reproduced_exactly(rule):
    net = generate_network(ID, rule=rule)
    assert net.status == "complete"
    assert [r.to_text() for r in net.reactions] == PUBLISHED[rule]["interactions"]
    assert [s.id for s in net.species] == list(C.ORDER)
    # the species structure is the wall the particle is, and its published velocity
    structures = {s.id: s.structure for s in net.species}
    for name, (left, right, velocity) in PUBLISHED[rule]["particles"].items():
        shown = "unstable" if velocity is None else velocity
        assert structures[name] == f"Λ{left[1:]}Λ{right[1:]} (v={shown})"


@pytest.mark.parametrize("rule", sorted(PUBLISHED))
def test_every_published_interaction_composes_the_walls_of_its_reactants(rule):
    """A collision of Li|Lj with Lj|Lk must give Li|Lk - the table is consistent."""
    spec = C.RULES[rule]
    wall = {name: (left, right) for name, (left, right, _v) in spec["particles"].items()}
    for kind, lhs, rhs, leaves in spec["interactions"]:
        left = wall[lhs[0]][0]
        right = wall[lhs[-1]][1]
        for first, second in zip(lhs, lhs[1:]):        # they must share a domain
            assert wall[first][1] == wall[second][0], (kind, lhs)
        if rhs:
            assert wall[rhs[0]][0] == left and wall[rhs[-1]][1] == right
            for first, second in zip(rhs, rhs[1:]):
                assert wall[first][1] == wall[second][0], (kind, rhs)
        else:
            # an annihilation leaves one domain behind, and both sides agree
            assert left == right == leaves, (kind, lhs, leaves)


# --- the published look-up tables ----------------------------------------------------
def test_gkl_hex_expands_to_the_published_majority_rule():
    # PNAS 1995 Table 1 gives GKL as 005f005f 005f005f 005fff5f 005fff5f; the rule
    # itself is: s = 0 -> majority(s, s[-1], s[-3]); s = 1 -> majority(s, s[+1], s[+3])
    assert np.array_equal(C.lookup(C.RULES["gkl"]["hex"]), C.gkl_lookup())


@pytest.mark.parametrize("rule", sorted(C.RULES))
def test_every_lookup_table_is_128_bits(rule):
    table = C.lookup(C.RULES[rule]["hex"])
    assert table.shape == (128,) and set(table.tolist()) <= {0, 1}
    with pytest.raises(ValueError, match="32 hex digits"):
        C.lookup("00ff")


def test_gkl_and_the_majority_rule_differ_and_gkl_is_symmetric():
    table = C.lookup(C.RULES["gkl"]["hex"])
    # the task is symmetric under exchanging 0s and 1s, and so is GKL:
    # the output for a neighbourhood is the complement of the reversed one
    for index in range(128):
        mirror = 127 - int(f"{index:07b}"[::-1], 2)
        assert table[index] == 1 - table[mirror]


# --- the velocities, measured by running the automaton --------------------------------
def test_probe_measures_the_published_velocities():
    """Build a lattice out of the two domains a particle separates and follow it."""
    for name, (_left, _right, published) in TABLE_3["particles"].items():
        got = C.probe("phi-par-a", name)
        if published is None:
            continue
        assert got["steps_followed"] == 20, f"{name} did not survive"
        assert abs(got["velocity"] - published) <= 0.05, (name, got, published)


def test_alpha_is_unstable_and_decays_into_gamma_and_mu():
    # Table 3 lists alpha with no velocity; the text: "alpha is unstable and
    # immediately decays into two particles gamma and mu"
    got = C.probe("phi-par-a", "alpha")
    assert got["decays_into"] == ["gamma", "mu"]
    assert got["steps_followed"] < 20


def test_gkl_gets_no_catalog_because_its_walls_are_not_those_particles():
    # GKL has the same three domains but its walls move differently, so no
    # published name or velocity is borrowed and no interaction table is claimed
    spec = C.RULES["gkl"]
    assert spec["particles"] == {} and spec["catalog"] is None
    with pytest.raises(ValueError, match="no particle interaction table"):
        generate_network(ID, rule="gkl")
    net = run(rule="gkl", seed=2)
    assert all(s.id.startswith("w") for s in net.species)
    assert all("measured, not published" in s.structure for s in net.species)


# --- the domain filter ----------------------------------------------------------------
def test_domain_filter_finds_the_published_domains_and_the_wall_between_them():
    spec = C.RULES["phi-par-a"]
    lattice = C.two_domain_lattice(spec, "L0", "L2", 80)
    labels = C.domain_labels(lattice, spec["domains"], 3)
    assert set(labels.tolist()) <= {-1, 0, 2}, "only Lambda0 and Lambda2 are present"
    found = C.walls(labels, list(spec["domains"]))
    pairs = sorted((left, right) for _c, left, right, _w in found)
    assert pairs == [("L0", "L2"), ("L2", "L0")], "one wall each way round the ring"
    assert [C.name_wall(l, r, C._by_wall("phi-par-a")) for l, r in pairs] == ["gamma", "delta"]


def test_a_pure_domain_has_no_walls_at_all():
    spec = C.RULES["phi-par-a"]
    for name in spec["domains"]:
        lattice = C.two_domain_lattice(spec, name, name, 60)
        labels = C.domain_labels(lattice, spec["domains"], 3)
        assert (labels >= 0).all()
        assert C.walls(labels, list(spec["domains"])) == []


# --- the default run ------------------------------------------------------------------
def test_default_run_is_the_density_rule_classifying_a_low_density_lattice():
    net = run(seed=1)
    analysis = net.extras["analysis"]
    assert net.status == "observed"
    assert analysis["rule"] == "phi-par-a" and analysis["task"] == "density-classification"
    assert round(analysis["initial_density"], 3) == 0.483     # 72 of 149 cells
    assert analysis["classification"] == "correct"
    assert net.extras["space"]["final"] in {"0" * 149, "1" * 149}, "a fixed point"
    # the condensation time is a handful of steps (about 12 on average at N = 149)
    assert 0 < analysis["condensation_time"] <= 30
    # every catalogued particle is a species, and the reactions carry firing counts
    assert set(C.ORDER) <= {s.id for s in net.species}
    assert net.reactions and all(r.count >= 1 for r in net.reactions)
    assert net.initial_state and all(v > 0 for v in net.initial_state.values())
    space = net.extras["space"]
    assert (space["shape"], space["radius"], space["states"]) == ([149], 3, 2)
    assert len(space["initial"]) == 149 and set(space["initial"]) <= set("01")
    assert set("".join(space["filtered"])) <= {".", "#"}


def test_a_frame_per_iteration():
    traj = evolve(ID, seed=1)
    net = traj.network
    assert traj.clock == "iterations" and [f.t for f in traj.frames] == [float(t) for t in range(299)]
    assert not traj.frames[0].fired
    density = traj.series("density")
    assert density[0] == pytest.approx(72 / 149) and density[-1] in (0.0, 1.0)
    # the particles present at the condensation time are the network's initial state
    t_c = net.extras["analysis"]["condensation_time"]
    assert traj.frames[t_c].state == net.initial_state
    # the lattice ends at a fixed point with no particle left
    assert traj.frames[-1].state == {}
    fired = Counter()
    for f in traj.frames:
        for lhs, rhs, n in f.fired:
            fired[(tuple(sorted(lhs)), tuple(sorted(rhs)))] += n
    assert fired == Counter({(multiset(r.reactants), multiset(r.products)): r.count for r in net.reactions})
    # nothing is recorded before condensation
    assert not any(f.fired for f in traj.frames[:t_c + 1])


def test_observed_reactions_are_collisions_of_catalogued_particles():
    net = run(seed=3, lattice=299, steps=598)
    published = {tuple(sorted(lhs)): tuple(sorted(rhs))
                 for _k, lhs, rhs, _l in C.RULES["phi-par-a"]["interactions"]}
    two = [r for r in net.reactions if sum(r.reactants.values()) == 2
           and set(r.reactants) <= set(C.ORDER)]
    assert two, "collisions between two catalogued particles were observed"
    for r in two:
        key = multiset(r.reactants)
        assert key in published, f"{r.to_text()} is not in the published table"
        assert multiset(r.products) == published[key], r.to_text()


def test_the_published_table_does_not_run_the_automaton():
    net = generate_network(ID, seed=1)
    assert "final" not in net.extras["space"] and "initial" not in net.extras["space"], "nothing was iterated"
    assert net.to_dict()["reactions"] == generate_network(ID, seed=2).to_dict()["reactions"]
    assert net.extras["space"]["lookup_hex"] == C.RULES["phi-par-a"]["hex"]
    assert net.extras["analysis"]["catalog"]["source"].startswith("Crutchfield")


# --- parameters -----------------------------------------------------------------------
def test_same_seed_same_run_and_a_different_seed_moves_the_particles():
    a, b = run(seed=5), run(seed=5)
    assert a.to_dict() == b.to_dict()
    c = run(seed=6)
    assert a.extras["space"]["initial"] != c.extras["space"]["initial"]


def test_bad_parameters():
    with pytest.raises(ValueError, match="lattice"):
        run(lattice=11)
    with pytest.raises(ValueError, match="rule"):
        run(rule="rule-110")
    # phi_par^b's published look-up table does not reproduce its published
    # behaviour, so only its catalog is offered
    with pytest.raises(ValueError, match="generate_network"):
        run(rule="phi-par-b")
    # the run's parameters belong to the evolve face, and the old switch is gone
    with pytest.raises(ValueError, match="belongs to the evolve face"):
        generate_network(ID, lattice=75)
    with pytest.raises(ValueError, match="'reactions' is gone"):
        generate_network(ID, reactions="published")


def test_the_filter_window_widens_the_walls():
    narrow = run(seed=2, filter_window=3)
    wide = run(seed=2, filter_window=6)
    assert narrow.extras["space"]["initial"] == wide.extras["space"]["initial"]
    assert narrow.extras["analysis"]["particles_seen"] != \
        wide.extras["analysis"]["particles_seen"]


# --- published performance (slow) -----------------------------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("rule, published", [("phi-par-a", 0.775), ("gkl", 0.816)])
def test_density_classification_performance_is_in_the_published_range(rule, published):
    # Table 1: P_149 over 10^4 unbiased ICs; here 250 ICs, so the sampling
    # standard deviation is about 0.026
    got = C.performance(rule, n=149, trials=250, rng=np.random.default_rng(7))
    assert abs(got - published) < 0.08, f"{rule}: {got} against a published {published}"


@pytest.mark.slow
def test_performance_at_the_published_lattice_sizes():
    # Table 1: phi_par^a scores 0.775 at N = 149, 0.740 at N = 599 and 0.728 at
    # N = 999. The gaps between those are smaller than the sampling deviation of
    # a run this size (about 0.04), so each size is checked against its own value
    # rather than against the others.
    for n, trials, published in ((149, 150, 0.775), (599, 120, 0.740)):
        got = C.performance("phi-par-a", n=n, trials=trials, rng=np.random.default_rng(n))
        assert abs(got - published) < 0.10, f"N={n}: {got} against a published {published}"


@pytest.mark.slow
def test_a_real_run_fires_every_interaction_of_the_published_table():
    """The whole point: the published table is what the automaton actually does."""
    fired: Counter = Counter()
    condensation = []
    for seed in range(6):
        net = run(seed=seed, lattice=599, steps=1198,
                  density=0.48 if seed % 2 else 0.52)
        analysis = net.extras["analysis"]
        condensation.append(analysis["condensation_time"])
        for r in net.reactions:
            fired[(multiset(r.reactants), multiset(r.products))] += r.count
    for _kind, lhs, rhs, _leaves in C.RULES["phi-par-a"]["interactions"]:
        key = (tuple(sorted(lhs)), tuple(sorted(rhs)))
        assert fired[key] > 0, f"{' + '.join(lhs)} -> {' + '.join(rhs) or '∅'} never fired"
    # the published collisions are the bulk of what happens between particles
    between = {k: v for k, v in fired.items() if len(k[0]) == 2}
    known = {(tuple(sorted(lhs)), tuple(sorted(rhs)))
             for _k, lhs, rhs, _l in C.RULES["phi-par-a"]["interactions"]}
    assert sum(v for k, v in between.items() if k in known) > 0.8 * sum(between.values())
    # Hordijk, Crutchfield & Mitchell (1998): t_c averages about 12 at N = 149
    assert max(condensation) < 60
