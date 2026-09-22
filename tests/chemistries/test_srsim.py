"""SRSim: what the 2010 paper and its supplement publish about spatial rule-based kinetics.

Gruenert, Ibrahim, Lenser, Lohel, Hinze & Dittrich, "Rule-based spatial modeling with
diffusing, geometrically constrained molecules", BMC Bioinformatics 11:307 (2010):
figure 3 (the geometry of the binding sites decides which complexes a polymerisation
rule can build), figures 5-6 with the input files and data of its additional file 3
(scaffold proteins co-localise their ligands in the spatial simulation but not in the
well-mixed one: 107 against 81 phosphorylated A in SRSim, 61 against 79 in BioNetGen),
and additional file 2 ("Calculation of Kinetic Parameters from Macroscopic Values"):
the reactive volume (2 pi / 3)(3 d0^2 t_dist + t_dist^3)(1 - cos t_ang)^2, the
equilibrium constant K = k1/k-1 it implies, and the refractory molecules
[A*] = [B*] = t_ref k_-1 [C].
"""

import math

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries.srsim import (SCAFFOLD_PREFACTORS, SCAFFOLD_RATES, canonical_text,
                                       parse_types, reactive_volume)

ID = "srsim"

#: a custom elementary species: two binding sites and one modifiable component
TYPES = [{"name": "M", "radius": 0.5, "mass": 1.0, "diffusion": 1.0,
          "sites": [{"name": "a", "dist": 1.0, "theta": 90.0, "phi": 0.0},
                    {"name": "b", "dist": 1.0, "theta": 90.0, "phi": 90.0},
                    {"name": "c", "dist": 1.0, "theta": 0.0, "phi": 0.0, "states": ["u", "p"]}]}]


@pytest.fixture(scope="module")
def net():
    return generate_network(ID, seed=1)


def dimer_counts(network):
    """Mean free A, free B and AB complexes over the second half of a dimerization run."""
    mean = network.extras["analysis"]["mean_counts"]
    free_a = sum(v for s, v in mean.items() if s.startswith("A(") and "," not in s)
    free_b = sum(v for s, v in mean.items() if s.startswith("B(") and "," not in s)
    return free_a, free_b, sum(v for s, v in mean.items() if "," in s)


def equilibrium_constant(network):
    """K = [C] / ([A][B]) over the molecules that are not waiting out a refractory time."""
    free_a, free_b, complexes = dimer_counts(network)
    waiting = network.extras["analysis"]["mean_refractory"] / 2.0      # per species
    volume = network.extras["kinetics"]["reactor_volume"]
    return complexes * volume / ((free_a - waiting) * (free_b - waiting))


def predicted_constant(network):
    """k1 / k-1 with the macroscopic k1 = k2_mic V_react of additional file 2."""
    rules = network.extras["kinetics"]["rules"]
    return rules["bind"]["k_macro"] / rules["bind_rev"]["k_macro"]


def well_mixed_dimers(seed, **params):
    return generate_network(ID, seed=seed, model="dimerization", well_mixed=True, n_molecules=60,
                            box=12.0, steps=3000, dt=0.02, k_on=0.5, k_off=0.05,
                            refractory_time=0.0, **params)


# ---------------------------------------------------------------------------
# Additional file 2: reactive volumes and the macroscopic rate constants
# ---------------------------------------------------------------------------
def test_reactive_volume_is_the_supplements_spherical_cone():
    d0, t_dist, t_ang = 2.0, 0.5, 60.0
    shell = (8.0 * math.pi / 3.0) * (3.0 * d0 ** 2 * t_dist + t_dist ** 3)
    cone = ((2.0 * math.pi / 3.0) * (3.0 * d0 ** 2 * t_dist + t_dist ** 3)
            * (1.0 - math.cos(math.radians(t_ang))) ** 2)
    assert reactive_volume(d0, t_dist, t_ang) == pytest.approx(shell)
    assert reactive_volume(d0, t_dist, t_ang, "sampled") == pytest.approx(cone)
    # at an angular tolerance of 180 degrees the orientation factors are 1 and both agree
    assert reactive_volume(d0, t_dist, 180.0, "sampled") == pytest.approx(shell)

    quiet = generate_network(ID, seed=1, model="dimerization", steps=1, well_mixed=True)
    bind = quiet.extras["kinetics"]["rules"]["bind"]
    assert bind["ideal_bond_length"] == pytest.approx(2.0)        # d_ij + d_kl = 1 + 1
    assert bind["reactive_volume"] == pytest.approx(shell)
    assert bind["k_macro"] == pytest.approx(0.5 * shell)          # k2_mac = k2_mic V_react


def test_well_mixed_equilibrium_matches_the_published_conversion():
    """Well stirred, the reactor reaches K = [C]/([A][B]) = k1/k-1 with k1 = k2_mic V_react."""
    nets = [well_mixed_dimers(seed) for seed in (1, 2, 3)]
    ratios = [equilibrium_constant(n) / predicted_constant(n) for n in nets]
    assert sum(ratios) / len(ratios) == pytest.approx(1.0, rel=0.2)
    assert all(0.6 < r < 1.6 for r in ratios)


def test_the_angular_tolerance_enters_as_one_minus_cos_squared():
    """With both orientations sampled, V_react carries (1 - cos t_ang)^2 (additional file 2)."""
    narrow = well_mixed_dimers(1, angle_tolerance=60.0, orientation="sampled")
    wide = well_mixed_dimers(1, angle_tolerance=180.0, orientation="sampled")
    factor = ((1.0 - math.cos(math.radians(60.0))) / 2.0) ** 2
    assert predicted_constant(narrow) == pytest.approx(predicted_constant(wide) * factor)
    assert equilibrium_constant(narrow) / predicted_constant(narrow) == pytest.approx(1.0, rel=0.35)


# ---------------------------------------------------------------------------
# The spatial simulation in the fast-diffusion limit
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def spatial_dimers():
    return [generate_network(ID, seed=seed, model="dimerization", n_molecules=50, box=12.0,
                             steps=2000, dt=0.02, k_on=0.5, k_off=0.05, refractory_time=2.0,
                             integrator="brownian") for seed in (1, 2)]


@pytest.mark.slow
def test_spatial_kinetics_reproduce_well_stirred_mass_action(spatial_dimers):
    """Diffusing particles reach the mass-action equilibrium constant of the supplement."""
    ratios = [equilibrium_constant(n) / predicted_constant(n) for n in spatial_dimers]
    assert sum(ratios) / len(ratios) == pytest.approx(1.0, rel=0.3)


@pytest.mark.slow
def test_refractory_molecules_follow_the_supplement(spatial_dimers):
    """[A*] = [B*] = t_ref k_-1 [C]: what the refractory time takes out of the reactor."""
    for network in spatial_dimers:
        _a, _b, complexes = dimer_counts(network)
        predicted = 2.0 * network.params["refractory_time"] * network.params["k_off"] * complexes
        assert network.extras["analysis"]["mean_refractory"] == pytest.approx(predicted, rel=0.25)


def test_free_diffusion_matches_six_d_t():
    """The paper's rule of thumb rests on <r^2> = 6 D t; with no rules that is all there is."""
    for integrator in ("brownian", "langevin"):
        free = generate_network(ID, seed=4, model="dimerization", n_molecules=40, k_on=0.0,
                                k_off=0.0, k_repulsion=0.0, integrator=integrator,
                                steps=2000, dt=0.02, box=40.0)
        walk = free.extras["analysis"]
        assert walk["mean_squared_displacement"] == pytest.approx(walk["expected_msd_6Dt"], rel=0.2)
        assert not free.reactions


# ---------------------------------------------------------------------------
# Figure 3: the geometry decides which complexes can form
# ---------------------------------------------------------------------------
def closed_sizes(network):
    """Sizes of the complexes with no free site left, i.e. the closed rings."""
    return sorted(s.structure.count("M(") for s in network.species if "[.]" not in s.structure)


@pytest.mark.slow
def test_geometry_decides_which_complexes_can_form():
    """Fig. 3: one rule, three geometries - triangles, no small ring at all, and rods."""
    triangles = generate_network(ID, seed=1, bond_angle=60.0, angle_tolerance=15.0)
    assert 3 in closed_sizes(triangles)

    # A closed ring of n monomers turns by 360 degrees in total (Fenchel), so its interior
    # angles cannot exceed 180 - 360/n: with 120 +/- 15 degrees no ring below
    # 360 / (180 - 120 + 15) = 4.8 monomers can ever close.
    for seed in (1, 2):
        wide = generate_network(ID, seed=seed, bond_angle=120.0, angle_tolerance=15.0)
        assert all(n >= 5 for n in closed_sizes(wide))
        assert wide.extras["analysis"]["largest_complex"] > 3        # chains still grow

    rods = generate_network(ID, seed=1, bond_angle=180.0, angle_tolerance=15.0)
    assert closed_sizes(rods) == [] and rods.extras["analysis"]["closed_complexes"] == 0


# ---------------------------------------------------------------------------
# Figures 5-6: scaffold proteins
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_scaffold_raises_phosphorylation_only_in_the_spatial_simulation():
    """Fig. 6: with the scaffold, phosphorylated A rises in SRSim (107 against 81) but not in
    BioNetGen (61 against 79), because the effect is co-localisation, not a different rate."""
    def run(**params):
        network = generate_network(ID, seed=1, model="scaffold", n_molecules=100, n_scaffolds=10,
                                   box=69.3, steps=1500, dt=0.5, diffusion=2.0,
                                   angle_tolerance=180.0, k_angle=0.0, rate_scale=100.0, **params)
        return network.extras["analysis"]["mean_state_counts"].get("A.r~P", 0.0), network

    spatial_with, network = run(scaffold_binding=1.0)
    spatial_without, _ = run(scaffold_binding=0.0)
    mixed_with, _ = run(scaffold_binding=1.0, well_mixed=True)
    mixed_without, _ = run(scaffold_binding=0.0, well_mixed=True)

    assert spatial_with > 1.3 * spatial_without
    assert mixed_with < 1.3 * mixed_without
    # the scaffold does bind its four ligands in both cases (S plus four A is a complex of five)
    assert network.extras["analysis"]["size_histogram"].get("5", 0) > 0
    assert SCAFFOLD_RATES == {"d1_on": 1.0e-4, "d1_off": 0.1, "d2_on": 0.05, "d2_off": 1.0e-3,
                              "d3_on": 1.0e-5, "d3_off": 0.9}
    phosphorylate = network.extras["kinetics"]["rules"]["phosphorylate"]
    assert phosphorylate["k_micro"] == pytest.approx(
        SCAFFOLD_RATES["d1_on"] * SCAFFOLD_PREFACTORS["modify2"] * 100.0)
    # the unmodified partner is a catalyst: it stands on both sides of the reaction
    assert any(r.catalysts for r in network.reactions)


# ---------------------------------------------------------------------------
# Complexes are species up to isomorphism
# ---------------------------------------------------------------------------
def test_canonical_ids_are_invariant_under_isomorphism():
    monomer = parse_types([{"name": "M", "radius": 0.5, "sites": [
        {"name": "a", "dist": 1.0, "theta": 90.0, "phi": 0.0},
        {"name": "b", "dist": 1.0, "theta": 90.0, "phi": 90.0}]}])[0]

    def chain(ids, close):
        kind = {i: monomer for i in ids}
        states = {i: {0: 0, 1: 0} for i in ids}
        links = {}
        for k in range(len(ids) if close else len(ids) - 1):
            i, j = ids[k], ids[(k + 1) % len(ids)]
            links[(i, 0)], links[(j, 1)] = (j, 1), (i, 0)
        return canonical_text(list(ids), kind, states, links)

    # the same 4-ring, written with other molecule numbers and started elsewhere
    assert chain([0, 1, 2, 3], True) == chain([7, 3, 9, 1], True) == chain([2, 0, 3, 1], True)
    assert chain([0, 1, 2, 3], False) != chain([0, 1, 2, 3], True)
    assert chain([0, 1, 2, 3], False) == chain([5, 6, 7, 8], False)

    # components with the same name are interchangeable, as in BNGL: it does not matter
    # which of the scaffold's four t sites carries which ligand
    scaffold, ligand = parse_types([
        {"name": "S", "radius": 8.0, "sites": [{"name": "t", "dist": 8.0, "theta": 0.0, "phi": phi}
                                               for phi in (0.0, 72.0, 144.0, 216.0)]},
        {"name": "A", "radius": 1.0, "sites": [{"name": "s", "dist": 1.0, "theta": 100.0,
                                                "phi": 0.0}]}])

    def loaded(slots):
        kind = {0: scaffold, 1: ligand, 2: ligand}
        states = {0: {i: 0 for i in range(4)}, 1: {0: 0}, 2: {0: 0}}
        links = {}
        for em, slot in zip((1, 2), slots):
            links[(0, slot)], links[(em, 0)] = (em, 0), (0, slot)
        return canonical_text([0, 1, 2], kind, states, links)

    assert loaded((0, 1)) == loaded((0, 2)) == loaded((2, 3))


def test_elementary_molecules_are_conserved(net):
    ids, reactants, products = net.matrices()
    stoichiometry = (products - reactants).toarray()
    assert net.extras["conservation"]
    for law in net.extras["conservation"]:
        vector = np.array([law["vector"][s] for s in ids])
        assert np.all(vector @ stoichiometry == 0), law["name"]
    assert sum(net.initial_state.values()) == net.params["n_molecules"]
    assert net.extras["analysis"]["size_histogram"]


def test_the_rule_set_is_kept_unflattened(net):
    assert net.status == "observed" and net.reactions
    assert all(r.count and r.count > 0 for r in net.reactions)
    assert net.extras["rules"] == ["'grow' M(a) + M(b) <-> M(a!1).M(b!1) @ 0.5, 0.01"]
    assert "begin reaction rules" in net.extras["bngl"]
    assert net.extras["molecule_types"][0]["sites"][1]["phi"] == 90.0
    assert net.extras["space"]["ideal_angles"]["M"][0][1] == pytest.approx(90.0)
    assert net.extras["space"]["tolerances"] == {"distance": 0.5, "angle_degrees": 30.0,
                                                 "orientation": "bound-only"}
    assert len(net.extras["space"]["particles"]) == net.params["n_molecules"]
    assert all(s.structure and " " not in s.id for s in net.species)


def test_a_custom_model_runs_its_own_rules():
    custom = generate_network(ID, seed=2, model="custom", molecule_types=TYPES,
                              rules=["'bind' M(a) + M(b) -> M(a!1).M(b!1) @ 0.5",
                                     "'phos' M(c~u) -> M(c~p) @ 0.5"],
                              init={"M": 30}, steps=500)
    assert custom.extras["analysis"]["events_by_rule"]["phos"] > 0
    assert custom.extras["analysis"]["events_by_rule"]["bind"] > 0
    assert any("{p}" in s.structure for s in custom.species)


def test_a_frame_per_step_at_the_simulated_time():
    traj = evolve(ID, seed=1, steps=100)
    assert traj.clock == "time"
    assert traj.times() == pytest.approx([0.02 * k for k in range(101)])
    net = traj.network
    assert traj.frames[0].state == net.initial_state and traj.frames[0].fired == []
    final = net.extras["analysis"]["final_counts"]
    assert traj.frames[-1].state == {s: float(n) for s, n in final.items()}
    # every frame holds the 50 monomers, bound or not
    (law,) = net.extras["conservation"]
    for frame in traj.frames:
        assert sum(law["vector"][s] * n for s, n in frame.state.items()) == 50


def test_the_seed_reproduces_the_run():
    first = generate_network(ID, seed=11, steps=400)
    assert first.to_dict() == generate_network(ID, seed=11, steps=400).to_dict()
    assert first.to_dict() != generate_network(ID, seed=12, steps=400).to_dict()


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError, match="only used with model='custom'"):
        generate_network(ID, rules=["'x' M(a) + M(b) -> M(a!1).M(b!1) @ 1"])
    with pytest.raises(ValueError, match="exactly one bond"):
        generate_network(ID, model="custom", molecule_types=TYPES, init={"M": 2},
                         rules=["'x' M(a,c~u) + M(b) -> M(a!1,c~p).M(b!1) @ 1"])
    with pytest.raises(ValueError, match="no site"):
        generate_network(ID, model="custom", molecule_types=TYPES, init={"M": 2},
                         rules=["'x' M(z) -> M(z) @ 1"])
    with pytest.raises(ValueError, match="one or two elementary molecules"):
        generate_network(ID, model="custom", molecule_types=TYPES, init={"M": 3},
                         rules=["'x' M(a) + M(b) + M(a) -> M(a!1).M(b!1).M(a) @ 1"])
    with pytest.raises(ValueError, match="init"):
        generate_network(ID, model="custom", molecule_types=TYPES,
                         rules=["'x' M(a) + M(b) -> M(a!1).M(b!1) @ 1"])
    with pytest.raises(ValueError, match="box"):
        generate_network(ID, box=1.5)
