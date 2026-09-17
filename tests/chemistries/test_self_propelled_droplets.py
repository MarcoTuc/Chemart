"""Self-propelled oil droplets: the published interfacial chemistry and the
measured behaviour of the droplets.

Numeric sources:
  - Hanczyc & Ikegami, Artificial Life 16(3):233-243 (2010): the chemical
    system, the Marangoni mechanism, the tensiometry of figure 2 and the
    materials and methods of the appendix.
  - Horibe, Hanczyc & Ikegami, Entropy 13(3):709-719 (2011): the behavioural
    modes by droplet size, the size series and the two-droplet attraction.
  - Banzhaf & Yamamoto (2015) section 19.2.5: chemotaxis toward the highest pH.
Everything asserted here is stated in one of those; nothing is fitted.
"""

import math

import pytest

from chemart import generate_network
from chemart.chemistries.self_propelled_droplets import (
    ALIFE2010, HORIBE, TENSION, droplet_geometry, fuel_budget,
)

ID = "self-propelled-droplets"


def net(**kw):
    return generate_network(ID, seed=1, **kw)


# --------------------------------------------------------------------------
# The published reaction set
# --------------------------------------------------------------------------
def test_hydrolysis_stoichiometry_of_the_oleic_anhydride_system():
    """One anhydride plus water gives TWO oleic acids, which deprotonate to
    oleate and protons (Hanczyc & Ikegami 2010, section 3.1)."""
    n = net()
    written = {r.to_text() for r in n.reactions}
    assert written == {
        "oleic_anhydride + H2O -> 2 oleic_acid",
        "oleic_acid -> oleate + H",
    }
    (hydrolysis,) = [r for r in n.reactions if "oleic_anhydride" in r.reactants]
    assert hydrolysis.reactants == {"oleic_anhydride": 1, "H2O": 1}
    assert hydrolysis.products == {"oleic_acid": 2}       # the 2 is the whole point
    # The proton is a product, not a catalyst: it is what lowers the local pH.
    (deprotonation,) = [r for r in n.reactions if "H" in r.products]
    assert deprotonation.products == {"oleate": 1, "H": 1}
    assert not deprotonation.catalysts


def test_overall_reaction_makes_two_oleates_and_two_protons():
    """The sum of the two steps is the sentence in the paper: the precursor
    'is hydrolyzed at the oil-water interface to produce more oleate and
    protons'."""
    n = net()
    ids, R, P = n.matrices()
    S = (P - R).toarray()
    index = {s: i for i, s in enumerate(ids)}
    # Fire the hydrolysis once and the deprotonation twice.
    extent = [0, 0]
    for j, r in enumerate(n.reactions):
        extent[j] = 1 if "oleic_anhydride" in r.reactants else 2
    net_change = S @ extent
    assert net_change[index["oleic_anhydride"]] == -1
    assert net_change[index["H2O"]] == -1
    assert net_change[index["oleic_acid"]] == 0      # fully consumed again
    assert net_change[index["oleate"]] == 2
    assert net_change[index["H"]] == 2


def test_oleoyl_groups_are_conserved():
    """One anhydride carries two oleoyl groups, one acid and one oleate carry
    one each, so m = (2, 1, 1) is a conservation law: S^T m = 0."""
    n = net()
    (law,) = n.extras["conservation"]
    assert law["name"] == "oleoyl groups"
    assert law["vector"] == {"oleic_anhydride": 2, "oleic_acid": 1, "oleate": 1}

    ids, R, P = n.matrices()
    S = (P - R).toarray()
    m = [law["vector"].get(s, 0) for s in ids]
    assert list(S.T @ m) == [0, 0]
    assert "mass-conservation" in n.provides


def test_no_rate_constant_is_invented():
    """No accessible source publishes a rate for the interfacial hydrolysis."""
    for system in ("oleic-anhydride", "fuel-surfactant", "maze-chemotaxis"):
        n = net(system=system)
        assert all(r.rate is None for r in n.reactions), system
        assert "rate-constants" not in n.provides
        assert "rate-law" not in n.provides


# --------------------------------------------------------------------------
# Published numbers: tensiometry and composition
# --------------------------------------------------------------------------
def test_published_tensiometry_of_figure_2():
    """Figure 2 caption: a bare nitrobenzene droplet at pH 11 with no oleate has
    an interfacial tension of 27 mN/m. Section 3.1: the tension rises as the pH
    falls, reaching a maximum at pH 9, and the local pH falls as low as 7."""
    assert TENSION["bare_nitrobenzene_mN_per_m"] == pytest.approx(27.0)
    assert TENSION["maximum_at_pH"] == pytest.approx(9.0)
    assert TENSION["bulk_pH"] == pytest.approx(11.0)
    assert TENSION["local_pH_minimum"] == pytest.approx(7.0)
    assert TENSION["oleate_mM"] == pytest.approx(10.0)
    # The droplets run above the tension maximum, which is why it takes a local
    # acidification - not extra surfactant - to break the symmetry.
    assert TENSION["bulk_pH"] > TENSION["maximum_at_pH"]
    assert TENSION["local_pH_minimum"] < TENSION["maximum_at_pH"]


def test_default_composition_is_the_published_one():
    """Appendix A.1/A.2: 0.5 M oleic anhydride in nitrobenzene, added to 10 mM
    oleate micelles at pH 11."""
    n = net()
    assert n.params["precursor_M"] == pytest.approx(0.5)
    assert n.params["surfactant_mM"] == pytest.approx(10.0)
    assert n.params["pH"] == pytest.approx(11.0)
    assert n.initial_state["oleic_anhydride"] == pytest.approx(500.0)   # mM
    assert n.initial_state["oleate"] == pytest.approx(10.0)
    # pH 11 is 1e-11 M = 1e-8 mM of protons.
    assert n.initial_state["H"] == pytest.approx(1e-8)
    assert n.extras["analysis"]["above_tension_maximum"] is True


def test_the_pH_of_the_tension_maximum_is_the_threshold():
    """Below the pH 9 maximum the droplet is no longer on the low-tension side."""
    assert net(pH=11.0).extras["analysis"]["above_tension_maximum"] is True
    assert net(pH=8.0).extras["analysis"]["above_tension_maximum"] is False
    assert net(pH=8.0).initial_state["H"] == pytest.approx(1e-5)       # mM


# --------------------------------------------------------------------------
# Droplet geometry and the interfacial fuel budget
# --------------------------------------------------------------------------
def test_droplet_geometry_is_a_sphere_of_the_published_volume():
    g = droplet_geometry(20.0)
    assert (4.0 / 3.0) * math.pi * g["radius_mm"] ** 3 == pytest.approx(20.0)
    assert g["surface_area_mm2"] == pytest.approx(4 * math.pi * g["radius_mm"] ** 2)
    assert g["surface_to_volume_per_mm"] == pytest.approx(3.0 / g["radius_mm"])
    # Hanczyc & Ikegami's convective droplets are about 0.1 mm across; the 20 uL
    # droplets of the behavioural experiments are centimetre-scale by comparison.
    assert ALIFE2010["convective_droplet_diameter_mm"] == pytest.approx(0.1)
    assert g["diameter_mm"] > 3.0


def test_surface_to_volume_falls_with_size():
    """The reaction is interfacial, which is why Horibe et al. say the reaction
    rate is 'largely determined by the size of droplet'."""
    ratios = [droplet_geometry(v)["surface_to_volume_per_mm"] for v in HORIBE["volumes_uL"]]
    assert ratios == sorted(ratios, reverse=True)
    areas = [droplet_geometry(v)["surface_area_mm2"] for v in HORIBE["volumes_uL"]]
    assert areas == sorted(areas)          # bigger droplets: more interface...
    assert ratios[0] > ratios[-1]          # ...but less interface per unit fuel


def test_one_anhydride_yields_two_surfactants():
    b = fuel_budget(20.0, 0.5, 2)
    assert b["precursor_mol"] == pytest.approx(0.5 * 20e-6)
    assert b["surfactant_mol"] == pytest.approx(2 * b["precursor_mol"])
    assert net().extras["analysis"]["fuel_budget"]["surfactant_per_precursor"] == 2


# --------------------------------------------------------------------------
# Published behaviour: Horibe et al. (2011)
# --------------------------------------------------------------------------
def test_published_size_series_and_modes():
    """Figure 1 and section 2.4: 1 uL droplets are circular and fluctuating,
    20 uL add a directional mode, and only the 50 uL droplets vibrate."""
    assert HORIBE["volumes_uL"] == [1.0, 3.0, 10.0, 20.0, 30.0, 50.0]
    assert set(HORIBE["modes"]) == {"directional", "circular", "fluctuating", "vibrating"}
    assert HORIBE["modes_by_volume_uL"]["1"] == ["circular", "fluctuating"]
    assert HORIBE["modes_by_volume_uL"]["20"] == ["directional", "circular", "fluctuating"]
    assert HORIBE["modes_by_volume_uL"]["50"] == ["vibrating", "circular", "fluctuating"]
    # "Both small (1 uL) and middle size (20 uL) droplets did not show vibrating mode."
    for volume, modes in HORIBE["modes_by_volume_uL"].items():
        assert ("vibrating" in modes) == (float(volume) in HORIBE["vibrating_only_at_uL"])


def test_modes_are_reported_for_the_droplet_that_was_asked_for():
    assert net(droplet_volume_uL=20.0).extras["analysis"]["modes_at_this_volume"] == [
        "directional", "circular", "fluctuating"]
    assert net(droplet_volume_uL=50.0).extras["analysis"]["modes_at_this_volume"] == [
        "vibrating", "circular", "fluctuating"]
    # 3, 10 and 30 uL were run but their mode lists were not printed per size.
    assert net(droplet_volume_uL=3.0).extras["analysis"]["modes_at_this_volume"] is None


def test_collective_attraction_at_3_and_20_but_not_50_microlitres():
    """Section 2.5: two 20 uL droplets in one dish stay closer than two in
    separate dishes; 3 uL behaves the same way; 50 uL shows no attraction."""
    series = {row["volume_uL"]: row for row in net().extras["analysis"]["size_series"]}
    assert series[3.0]["collective_attraction"] is True
    assert series[20.0]["collective_attraction"] is True
    assert series[50.0]["collective_attraction"] is False
    assert series[10.0]["collective_attraction"] is None      # not reported
    assert HORIBE["observation_minutes"] == pytest.approx(60.0)
    assert HORIBE["window_minutes"] == pytest.approx(20.0)


def test_qualitative_velocity_findings_are_recorded_not_invented():
    """No source gives a droplet speed, so velocity survives only as the two
    published relative statements."""
    findings = " ".join(HORIBE["findings"])
    assert "velocity falls and turning angle rises as the droplet ages" in findings
    assert "negatively correlated" in findings
    analysis = net().extras["analysis"]
    assert analysis["motion_is_simulated"] is False
    assert "no accessible source gives the force law or a droplet speed" in \
        analysis["why_no_trajectory"]


# --------------------------------------------------------------------------
# The propulsion law, recorded as data
# --------------------------------------------------------------------------
def test_the_marangoni_mechanism_is_recorded_as_an_ordered_causal_chain():
    law = net().extras["interaction_law"]
    assert law["configuration"].startswith("onboard fuel")
    steps = law["steps"]
    assert len(steps) == 6
    assert "hydrolyses to" in steps[0]
    assert "local pH is not even" in steps[1]
    assert "tension rises as local pH falls" in steps[2]
    assert "Marangoni force tangential to the surface" in steps[3]
    assert "sustains the imbalance" in steps[4]
    assert "low-pH trail" in steps[5]
    assert law["sensor"] == "the oil-water interface itself"
    assert law["motor"] == "the convective flow structure inside the droplet"


def test_chemotaxis_runs_up_the_pH_gradient_for_the_hanczyc_droplet():
    """Book 19.2.5: the droplet follows a pH gradient, 'moving towards the
    highest pH'."""
    chemotaxis = net().extras["interaction_law"]["chemotaxis"]
    assert chemotaxis["stimulus"].startswith("pH gradient")
    assert chemotaxis["direction"] == "up the gradient, toward the highest pH"


def test_the_maze_solver_runs_the_other_way_down_the_pH_gradient():
    """Lagzi et al. (2010): these droplets are chemotactic toward LOW pH. The
    two published systems genuinely disagree in sign."""
    chemotaxis = net(system="maze-chemotaxis").extras["interaction_law"]["chemotaxis"]
    assert chemotaxis["direction"] == "down the gradient, toward the low-pH region"
    assert chemotaxis["direction"] != \
        net().extras["interaction_law"]["chemotaxis"]["direction"]


def test_no_direction_is_asserted_where_none_was_published():
    """The accessible abstract of Cejkova et al. reports induction time and
    migration velocity, but not the sign of the salt response."""
    for system in ("decanol-salt", "fuel-surfactant"):
        chemotaxis = net(system=system).extras["interaction_law"]["chemotaxis"]
        assert chemotaxis["direction"] is None, system


# --------------------------------------------------------------------------
# The other published droplet chemistries
# --------------------------------------------------------------------------
def test_the_fuel_surfactant_droplet_keeps_its_catalyst():
    """Toyota et al. (2009): 5 mol % of an amphiphilic catalyst hydrolyses an
    exogenous fuel, so the droplet itself is not consumed."""
    n = net(system="fuel-surfactant")
    (reaction,) = n.reactions
    assert reaction.catalysts == {"catalyst": 1}
    assert "catalysts" in n.provides
    assert reaction.reactants["precursor"] == 1
    assert reaction.products["octylaniline"] == 1
    assert "not consumed" in n.extras["interaction_law"]["configuration"]
    assert "5 mol %" in n.extras["interaction_law"]["steps"][0]
    # The catalyst is conserved because it is never used up.
    names = {law["name"] for law in n.extras["conservation"]}
    assert names == {"octylaniline moiety", "catalyst"}


def test_the_maze_droplet_is_an_acid_base_pair():
    n = net(system="maze-chemotaxis")
    written = {r.to_text() for r in n.reactions}
    assert written == {
        "hexyldecanoic_acid + OH -> hexyldecanoate + H2O",
        "hexyldecanoate + H -> hexyldecanoic_acid",
    }
    (law,) = n.extras["conservation"]
    assert law["vector"] == {"hexyldecanoic_acid": 1, "hexyldecanoate": 1}


def test_the_decanol_droplet_has_no_reactive_chemistry_at_all():
    """Hanczyc (2014) classifies the salt-gradient droplets as type 1: 'no
    reactive chemistry'. The reaction list is therefore empty."""
    n = net(system="decanol-salt")
    assert n.reactions == [] and n.status == "complete"
    assert "no reactions; see extras['interaction_law']" in n.summary()
    assert "stoichiometry" not in n.provides
    assert {s.id for s in n.species} == {"decanol", "decanoate", "Na", "Cl", "H2O"}
    assert n.extras["analysis"]["reactive"] is False
    assert "no reactive chemistry" in n.extras["interaction_law"]["configuration"]
    # Nothing reacts, so there is no moiety to conserve either.
    assert "conservation" not in n.extras
    assert "mass-conservation" not in n.provides


# --------------------------------------------------------------------------
# What the network exports
# --------------------------------------------------------------------------
def test_space_holds_the_published_assay_geometry_not_a_trajectory():
    space = net().extras["space"]
    assert space["units"] == "mm"
    assert "not simulated" in space["positions"]
    vessels = {v["diameter_mm"] for v in space["vessels"]}
    # The 27 mm base of the 35 mm dish is the "diameter of the liquid recipient"
    # in the book's figure 19.11; the mode experiments used a 100 mm dish.
    assert {15.0, 35.0, 60.0, 100.0} <= vessels
    base = [v for v in space["vessels"] if v.get("base_diameter_mm")]
    assert all(v["base_diameter_mm"] == 27.0 for v in base)


def test_phases_say_which_species_is_in_which_liquid():
    phases = net().extras["phases"]
    assert phases["oil"] == ["nitrobenzene", "oleic_anhydride"]
    assert set(phases["aqueous"]) == {"H2O", "oleate", "H"}
    assert phases["interface"] == ["oleic_acid"]        # where the reaction happens
    assert net().extras["buffered"] == ["H2O"]


def test_same_seed_same_network_and_parameters_change_it():
    assert net().to_dict() == net().to_dict()
    assert net(droplet_volume_uL=50.0).to_dict() != net().to_dict()
    assert net(system="decanol-salt").to_dict() != net().to_dict()


def test_rejects_bad_parameters():
    with pytest.raises(ValueError, match="system"):
        generate_network(ID, seed=1, system="water-in-oil")
    with pytest.raises(ValueError, match="droplet_volume_uL"):
        generate_network(ID, seed=1, droplet_volume_uL=0.0)
    with pytest.raises(ValueError, match="pH"):
        generate_network(ID, seed=1, pH=15.0)
