"""Ono & Ikegami autopoietic protocells: the book's reaction set, the published
constants of the 1D predecessor, and the membrane phenomena of book 6.3.2.

The numeric source is Ono & Ikegami, J. theor. Biol. 206:243-253 (2000),
appendices A and B. The spatial claims come from book section 6.3.2 and its
figure 6.14 (isotropic M gives clusters that die, anisotropic M gives membrane
filaments).
"""

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.ono_ikegami_protocell import JTB2000, anisotropy_field

SMALL = dict(width=12, height=12, steps=8, relaxation=4)
NOCHEM = dict(P_A=0.0, P_M=0.0, P_decay=0.0, X_supply=0.0)


def analysis(**kw):
    return evolve("ono-ikegami-protocell", **kw).network.extras["analysis"]


# --------------------------------------------------------------------------
# The reaction set (book 6.3.2)
# --------------------------------------------------------------------------
def test_reaction_set_is_the_book_s_four_reactions():
    net = generate_network("ono-ikegami-protocell", seed=1)
    assert net.status == "complete"
    assert [s.id for s in net.species] == ["A", "M_a", "X", "Y", "W"]
    written = {r.to_text() for r in net.reactions}
    assert written == {
        "A + X -> 2 A",          # autocatalytic replication
        "A + X -> A + M_a",      # membrane production, catalysed by A
        "A -> Y", "M_a -> Y", "X -> Y",   # everything but water decays
        "Y -> X",                # recycling by the external energy source
    }
    # A is a catalyst of membrane production, and water never decays.
    (membrane,) = [r for r in net.reactions if "M_a" in r.products]
    assert membrane.catalysts == {"A": 1}
    assert all("W" not in r.reactants for r in net.reactions)


def test_isotropic_membrane_is_a_distinct_species():
    net = generate_network("ono-ikegami-protocell", seed=1, membrane="isotropic")
    assert [s.id for s in net.species] == ["A", "M_i", "X", "Y", "W"]


def test_no_rate_constant_is_invented():
    """No accessible source publishes rates for this 2D scheme, so every rate is None."""
    net = generate_network("ono-ikegami-protocell", seed=1)
    assert all(r.rate is None for r in net.reactions)
    assert "rate-constants" not in net.provides


# --------------------------------------------------------------------------
# Published numbers: Ono & Ikegami (2000), appendices A and B
# --------------------------------------------------------------------------
def test_published_constants_of_the_one_dimensional_predecessor():
    """Appendix A: P_E = 0.5e-6, P_A' = 2e-6, P_R = P_W = 100e-6,
    P_D = 7e-3, P_R0 = 5e-3, P_R1 = 0.5 P_R0."""
    assert JTB2000["P_E"] == pytest.approx(0.5e-6)
    assert JTB2000["P_A_spontaneous"] == pytest.approx(2e-6)
    assert JTB2000["P_R_recycling"] == pytest.approx(100e-6)
    assert JTB2000["P_W_decay"] == pytest.approx(100e-6)
    assert JTB2000["P_D_diffusion"] == pytest.approx(7e-3)
    assert JTB2000["P_R0_repulsion_same_site"] == pytest.approx(5e-3)
    assert JTB2000["P_R1_repulsion_neighbour_site"] == pytest.approx(
        0.5 * JTB2000["P_R0_repulsion_same_site"])


def test_appendix_b_waste_fixed_point():
    """W_0 = 100/(1 + P_R/P_W) = 50 at the paper's P_R = P_W = 100e-6."""
    w0 = 100.0 / (1.0 + JTB2000["P_R_recycling"] / JTB2000["P_W_decay"])
    assert w0 == pytest.approx(50.0)


def test_phase_diagram_points_of_figures_4_7_and_10():
    """The (P_A, P_M) values labelling each regime of the paper's figures."""
    points = JTB2000["phase_points"]
    assert points["fig4b_periodic_membranes"] == {"P_A": 9e-6, "P_M": 7e-6}
    assert points["fig7b_region_Ia_cell_dies"] == {"P_A": 5e-6, "P_M": 3e-6}
    assert points["fig7c_region_Ib_cell_stable"] == {"P_A": 7e-6, "P_M": 9e-6}
    assert points["fig7d_region_Ic_cell_reproduces"] == {"P_A": 11e-6, "P_M": 5e-6}
    # Recursive division needs the reduced autocatalyst mobility m = 5.
    assert points["fig10_recursive_division"]["m"] == 5


# --------------------------------------------------------------------------
# The repulsion field (book figure 6.13)
# --------------------------------------------------------------------------
def test_anisotropic_field_redistributes_but_does_not_add_repulsion():
    """M_i and M_a carry the same total repulsion: F has mean 1 for both."""
    iso = anisotropy_field(0.0, True)
    ani = anisotropy_field(2.4, False)
    assert iso.tolist() == [[1.0] * 6] * 6
    for o in range(6):
        assert ani[:, o].mean() == pytest.approx(1.0)
    # Strongest along the +-o axis, weakest (attractive past a = 2) on the flanks.
    assert ani[0, 0] == pytest.approx(3.4)
    assert ani[1, 0] == pytest.approx(-0.2)
    # Nematic: the two ends are equivalent, which is what leaves the mutual
    # energy of an exchanged pair unchanged.
    for k in range(6):
        for o in range(6):
            assert ani[(k + 3) % 6, o] == pytest.approx(ani[k, o])


# --------------------------------------------------------------------------
# Spatial phenomena
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_isotropic_membranes_form_droplets_and_anisotropic_ones_form_filaments():
    """Book 6.3.2 and figure 6.14: with only isotropic M particles the membrane
    collapses into compact clusters, while anisotropic M builds thin filaments.

    Measured by the mean number of membrane neighbours per membrane particle
    (about 5 inside a droplet, about 2 along a one-particle-thick filament) and
    by the fraction of membrane contacts aligned along both particles' axes
    (1/9 for random orientations).
    """
    for seed in (1, 2):
        common = dict(seed=seed, repulsion=5.0, relaxation=25, steps=60,
                      M_fraction=0.12, A_fraction=0.0, X_fraction=0.20, **NOCHEM)
        ani = analysis(membrane="anisotropic", **common)
        iso = analysis(membrane="isotropic", **common)

        assert iso["mean_membrane_coordination"] > ani["mean_membrane_coordination"] + 1.0
        assert iso["mean_membrane_coordination"] > 2.5, "isotropic M should ball up"
        assert ani["mean_membrane_coordination"] < 2.5, "anisotropic M should stay thin"
        # The isotropic droplet is one big cluster; the filaments are not.
        assert iso["largest_membrane_cluster"] > ani["largest_membrane_cluster"]
        # Only M_a organises its orientations; M_i has none to organise.
        assert ani["membrane_alignment"] > 0.30
        assert iso["membrane_alignment"] < 0.20


@pytest.mark.slow
def test_metabolism_collapses_without_recycling_of_waste():
    """The external energy source is what keeps the cell alive: with Y -> X
    switched off every particle ends up as waste and the metabolism dies."""
    common = dict(seed=1, repulsion=5.0, relaxation=20, steps=100,
                  initial="cell", cell_radius=5, width=24, height=24)
    fed_run = evolve("ono-ikegami-protocell", X_supply=0.30, **common)
    starved_run = evolve("ono-ikegami-protocell", X_supply=0.0, **common)
    fed, starved = (run.network.extras["analysis"] for run in (fed_run, starved_run))

    # the A count after each sweep
    a_fed = [f.state.get("A", 0) for f in fed_run.frames[1:]]
    a_starved = [f.state.get("A", 0) for f in starved_run.frames[1:]]
    # Fed: food is recycled, so the autocatalyst grows and holds.
    assert fed["counts"]["X"] > 15
    assert a_fed[-1] > 2 * a_fed[0]
    # Starved: food drains away and the autocatalyst falls back from its peak.
    assert starved["counts"]["X"] <= 10
    assert starved["counts"]["Y"] > 3 * fed["counts"]["Y"]
    assert a_starved[-1] < max(a_starved), "A should decay once the food is gone"
    assert a_starved[-1] < a_fed[-1]


@pytest.mark.slow
def test_spatial_run_reports_the_lattice_energies_and_measurements():
    traj = evolve("ono-ikegami-protocell", seed=1)
    net = traj.network
    assert net.status == "observed" and traj.clock == "sweeps"
    assert all(r.count is not None and r.count > 0 for r in net.reactions)
    assert {"space", "energies", "analysis"} <= set(net.extras)

    space = net.extras["space"]
    assert space["lattice"] == "hexagonal" and space["dimensions"] == 2
    assert len(space["neighbour_directions"]) == 6
    grid = space["final"]
    assert len(grid) == space["shape"][0]
    assert all(len(row) == space["shape"][1] for row in grid)
    assert set("".join(grid)) <= set("AMXYW")
    # Every cell holds exactly one particle, so the lattice is always full.
    assert sum(net.extras["analysis"]["counts"].values()) == space["shape"][0] * space["shape"][1]

    energies = net.extras["energies"]
    assert energies["hydrophilic"] == ["A", "W"] and energies["neutral"] == ["X", "Y"]
    assert energies["hydrophobic"] == ["M_a"]


def test_a_frame_per_sweep_counts_the_lattice():
    traj = evolve("ono-ikegami-protocell", seed=3, **SMALL)
    net = traj.network
    assert [f.t for f in traj.frames] == [float(t) for t in range(SMALL["steps"] + 1)]
    assert traj.frames[0].state == {s: n for s, n in net.initial_state.items() if n}
    assert all(sum(f.state.values()) == 144 for f in traj.frames)      # one particle per cell
    counts = net.extras["analysis"]["counts"]                           # keyed by lattice code
    assert traj.frames[-1].state == {("M_a" if c == "M" else c): float(n) for c, n in counts.items() if n}
    # the counts change exactly as the fired reactions say
    ids, _, _ = net.matrices()
    for before, after in zip(traj.frames, traj.frames[1:]):
        change = {s: after.state.get(s, 0) - before.state.get(s, 0) for s in ids}
        net_change = dict.fromkeys(ids, 0)
        for lhs, rhs, n in after.fired:
            for s in lhs:
                net_change[s] -= n
            for s in rhs:
                net_change[s] += n
        assert change == net_change


def test_same_seed_gives_the_same_lattice():
    kw = dict(seed=5, **SMALL)
    assert evolve("ono-ikegami-protocell", **kw).to_dict() == evolve("ono-ikegami-protocell", **kw).to_dict()
    assert analysis(**kw) == analysis(**kw)
    assert analysis(seed=6, **SMALL) != analysis(**kw)


# --------------------------------------------------------------------------
# Parameter guards
# --------------------------------------------------------------------------
def test_rejects_a_cell_too_big_for_its_lattice():
    with pytest.raises(ValueError, match="does not fit"):
        evolve("ono-ikegami-protocell", seed=1, initial="cell",
               cell_radius=20, width=24, height=24)


def test_rejects_an_overfull_initial_composition():
    with pytest.raises(ValueError, match="exceeds 1"):
        evolve("ono-ikegami-protocell", seed=1, A_fraction=0.6,
               X_fraction=0.6, **SMALL)


def test_the_reaction_set_takes_no_lattice_parameters():
    with pytest.raises(ValueError, match="belongs to the evolve face"):
        generate_network("ono-ikegami-protocell", steps=10)
    with pytest.raises(ValueError, match="is gone"):
        generate_network("ono-ikegami-protocell", mode="spatial")
