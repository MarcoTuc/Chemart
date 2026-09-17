"""Self-propelled oil droplets (book 19.2.5). Catalog id: self-propelled-droplets.

A droplet is not a reaction: it is a body that moves. What *is* a chemistry
here is the small interfacial reaction set that breaks the droplet's symmetry,
so this module returns the in-silico network of the published experiment - the
hydrolysis that makes surfactant at the oil-water interface - and records the
motion as data.

The canonical system (``system="oleic-anhydride"``, Hanczyc, Toyota, Ikegami,
Packard & Sugawara 2007; Hanczyc & Ikegami 2010) is a nitrobenzene droplet
loaded with oleic anhydride, placed in 10 mM aqueous oleate micelles at pH 11.
At the interface the anhydride hydrolyses to two oleic acids, which deprotonate
to oleate and protons::

    oleic_anhydride + H2O -> 2 oleic_acid
    oleic_acid -> oleate + H

The protons lower the local pH; interfacial tension *rises* as pH falls
(maximum at pH 9), so an uneven reaction leaves an uneven tension, which is a
Marangoni force tangential to the surface, which drives a pair of convective
vortices, which carries fresh precursor to one pole and expels product at the
other - and the droplet swims. Nothing in this module integrates that motion:
**no accessible source publishes the force law**, so the mechanism, the
tension relation and the observed behaviour are recorded as data in
``extras["interaction_law"]``, ``extras["analysis"]`` and ``extras["space"]``,
and every rate is ``None``.

Three further published droplet chemistries are available through ``system``:
Toyota et al. (2009) ``fuel-surfactant`` (onboard catalyst, exogenous fuel),
Lagzi et al. (2010) ``maze-chemotaxis`` (acid/base, solves a maze), and
Cejkova et al. (2014) ``decanol-salt`` (no reactive chemistry at all, hence an
empty reaction list). See the catalog's `decisions` for what each is based on.
"""

from __future__ import annotations

import math

from chemart.helpers.explicit import network
from chemart.network import Species

SYSTEMS = ("oleic-anhydride", "fuel-surfactant", "maze-chemotaxis", "decanol-salt")

#: Interfacial tension facts of the oleic-anhydride system, from the pendant-drop
#: tensiometry of Hanczyc & Ikegami (2010), figure 2 and its caption. The curve
#: itself was published only as a picture, so only the stated values are kept.
TENSION = {
    "source": "Hanczyc & Ikegami, Artificial Life 16(3):233-243 (2010), figure 2 "
              "and section 3.1; pendant-drop tensiometry, each point in triplicate",
    "bare_nitrobenzene_mN_per_m": 27.0,   # nitrobenzene at pH 11 with no oleate
    "maximum_at_pH": 9.0,                 # tension rises to a maximum at pH 9
    "bulk_pH": 11.0,                      # the pH the droplets are added to
    "local_pH_minimum": 7.0,              # measured with pH-sensitive dyes
    "oleate_mM": 10.0,                    # the aqueous phase the droplets run in
    "relation": "interfacial tension increases as pH decreases from the bulk pH 11, "
                "reaching a maximum at pH 9; at pH 11 the tension is already low, so "
                "producing more oleate has little further effect and it is the local "
                "drop in pH, not the extra surfactant, that moves the droplet",
}

#: The behaviour that Horibe, Hanczyc & Ikegami (Entropy 13(3):709-719, 2011)
#: measured on this system, droplet size by droplet size.
HORIBE = {
    "source": "Horibe, Hanczyc & Ikegami, Entropy 13(3):709-719 (2011), "
              "doi:10.3390/e13030709 (open access)",
    "oil": "fresh oleic anhydride mixed with nitrobenzene at 1:1 v/v",
    "aqueous": "10 mM oleate micelles at pH 11",
    "volumes_uL": [1.0, 3.0, 10.0, 20.0, 30.0, 50.0],
    "modes": ["directional", "circular", "fluctuating", "vibrating"],
    # Figure 1: the modes each size was seen in, over the whole experiment.
    "modes_by_volume_uL": {
        "1": ["circular", "fluctuating"],
        "20": ["directional", "circular", "fluctuating"],
        "50": ["vibrating", "circular", "fluctuating"],
    },
    "vibrating_only_at_uL": [50.0],
    "collective_attraction_uL": [3.0, 20.0],
    "no_collective_attraction_uL": [50.0],
    "observation_minutes": 60.0,
    "window_minutes": 20.0,
    "sampling_interval_s": 1.0,
    "binning_interval_s": 20.0,
    "mean_coefficient_of_restitution": 1.0,
    "findings": [
        "velocity falls and turning angle rises as the droplet ages",
        "turning angle is negatively correlated with velocity: droplets change "
        "direction while stopped",
        "the reaction happens at the surface, so the reaction rate is set by droplet "
        "size, and it weakens as the oil volume is consumed",
        "spontaneous (non-collision) events give a coefficient of restitution with "
        "mean around 1, and values above 1 when the reaction is re-activated",
        "two droplets of 3 or 20 uL stay closer together than two droplets in "
        "separate dishes, and the attraction weakens over time; 50 uL droplets, "
        "whose convective flow is unstable, show no attraction",
    ],
}

#: Sizes and timings reported by Hanczyc & Ikegami (2010).
ALIFE2010 = {
    "source": "Hanczyc & Ikegami, Artificial Life 16(3):233-243 (2010)",
    "oil": "oleic anhydride in nitrobenzene at 0.5 M (appendix A.2)",
    "aqueous": "oleate micelles at pH 11, NaOH-adjusted (appendix A.1)",
    "diameter_range": "a few hundred microns to a few centimetres",
    "convective_droplet_diameter_mm": 0.1,
    "no_shape_fluctuation_below_um": 100.0,
    "shape_fluctuation_up_to_cm": 0.5,
    "shape_series_uL": [1.0, 5.0, 10.0, 30.0],
    "frame_interval_s": 8.0,
    "onset": "the droplet starts moving within seconds of the oil meeting the water",
    "without_precursor": "droplets still move chemotactically in a pH gradient, but the "
                         "motion is not sustained and stops within a few seconds once "
                         "the tension imbalance is equilibrated",
    "mechanism_shift": "above a few hundred microns self-movement shifts from "
                       "convection-driven to shape-driven, and a horseshoe shape "
                       "best supports straight directional motion",
}


# ---------------------------------------------------------------------------
# Geometry and fuel budget: published parameters, plain arithmetic
# ---------------------------------------------------------------------------
def droplet_geometry(volume_uL: float) -> dict:
    """Sphere of the given volume. 1 uL = 1 mm^3, so radii come out in mm."""
    radius = (3.0 * volume_uL / (4.0 * math.pi)) ** (1.0 / 3.0)
    return {
        "volume_uL": float(volume_uL),
        "radius_mm": radius,
        "diameter_mm": 2.0 * radius,
        "surface_area_mm2": 4.0 * math.pi * radius * radius,
        # 3/r: the reaction is interfacial, so this is what sets the rate per
        # unit of fuel, and it is why small droplets burn out first.
        "surface_to_volume_per_mm": 3.0 / radius,
    }


def fuel_budget(volume_uL: float, precursor_M: float, surfactant_per_precursor: int) -> dict:
    """Moles of precursor in the droplet and of surfactant it can make."""
    moles = precursor_M * volume_uL * 1e-6          # mol/L * L
    return {
        "precursor_mol": moles,
        "precursor_umol": moles * 1e6,
        "surfactant_per_precursor": int(surfactant_per_precursor),
        "surfactant_mol": moles * surfactant_per_precursor,
        "surfactant_umol": moles * surfactant_per_precursor * 1e6,
    }


# ---------------------------------------------------------------------------
# The four published droplet chemistries
# ---------------------------------------------------------------------------
def _oleic_anhydride(p):
    species = [
        Species("oleic_anhydride", "oleic anhydride, the precursor (fuel) dissolved "
                                   "in the nitrobenzene oil phase"),
        Species("H2O", "bulk aqueous phase"),
        Species("oleic_acid", "oleic acid, the protonated surfactant"),
        Species("oleate", "oleate, the surfactant anion that lowers interfacial tension"),
        Species("H", "proton released at the interface; it is the local drop in pH "
                     "that raises the tension and drives the Marangoni flow"),
        Species("nitrobenzene", "the oil solvent; it carries the precursor and does not react"),
    ]
    reactions = [
        # Hydrolysis at the oil-water interface. No accessible source publishes a
        # rate constant for it, so the rate is None (see the catalog's decisions).
        ("oleic_anhydride + H2O -> 2 oleic_acid", None),
        ("oleic_acid -> oleate + H", None),
    ]
    conservation = [{
        "name": "oleoyl groups",
        "vector": {"oleic_anhydride": 2, "oleic_acid": 1, "oleate": 1},
    }]
    initial = {
        "oleic_anhydride": p.precursor_M * 1000.0,
        "oleate": float(p.surfactant_mM),
        "H": 1000.0 * 10.0 ** (-p.pH),
    }
    return species, reactions, conservation, initial


def _fuel_surfactant(p):
    species = [
        Species("precursor", "amphiphilic precursor of 4-octylaniline, dispersed in the "
                             "aqueous phase as the exogenous 'fuel'"),
        Species("H2O", "bulk aqueous phase"),
        Species("catalyst", "amphiphilic catalyst, 5 mol % of the droplet; it is not "
                            "consumed, which is why the motion outlasts the droplet"),
        Species("octylaniline", "4-octylaniline, the oil the droplet is made of"),
        Species("byproduct", "the other hydrolysis product; the accessible abstract does "
                             "not name it (see the catalog's decisions)"),
    ]
    reactions = [
        ("precursor + H2O + catalyst -> octylaniline + byproduct + catalyst", None),
    ]
    conservation = [{
        "name": "octylaniline moiety",
        "vector": {"precursor": 1, "octylaniline": 1},
    }, {
        "name": "catalyst",
        "vector": {"catalyst": 1},
    }]
    initial = {
        "octylaniline": p.precursor_M * 1000.0,
        "precursor": float(p.surfactant_mM),
    }
    return species, reactions, conservation, initial


def _maze_chemotaxis(p):
    species = [
        Species("hexyldecanoic_acid", "2-hexyldecanoic acid, the protonated surfactant "
                                      "in the dichloromethane droplet"),
        Species("hexyldecanoate", "2-hexyldecanoate, the surface-active anion the "
                                  "droplet emits into the aqueous phase"),
        Species("OH", "hydroxide of the alkaline aqueous phase"),
        Species("H", "proton; the acid placed at one exit of the maze is its source"),
        Species("H2O", "bulk aqueous phase"),
        Species("dichloromethane", "the oil the droplet is made of; it does not react"),
    ]
    reactions = [
        ("hexyldecanoic_acid + OH -> hexyldecanoate + H2O", None),
        ("hexyldecanoate + H -> hexyldecanoic_acid", None),
    ]
    conservation = [{
        "name": "hexyldecanoyl groups",
        "vector": {"hexyldecanoic_acid": 1, "hexyldecanoate": 1},
    }]
    initial = {
        "hexyldecanoic_acid": p.precursor_M * 1000.0,
        "hexyldecanoate": float(p.surfactant_mM),
        "H": 1000.0 * 10.0 ** (-p.pH),
    }
    return species, reactions, conservation, initial


def _decanol_salt(p):
    species = [
        Species("decanol", "1-decanol, the droplet"),
        Species("decanoate", "sodium decanoate, the aqueous surfactant"),
        Species("Na", "sodium ion of the NaCl gradient"),
        Species("Cl", "chloride ion of the NaCl gradient"),
        Species("H2O", "bulk aqueous phase"),
    ]
    # Type 1 in Hanczyc's own classification: no reactive chemistry. The droplet
    # moves because an externally imposed gradient makes its tension uneven, so
    # there is genuinely nothing to put in the reaction list.
    return species, [], [], {"decanoate": float(p.surfactant_mM)}


_BUILD = {
    "oleic-anhydride": _oleic_anhydride,
    "fuel-surfactant": _fuel_surfactant,
    "maze-chemotaxis": _maze_chemotaxis,
    "decanol-salt": _decanol_salt,
}

#: Per-system provenance, propulsion mechanism and chemotaxis, all as recorded
#: data. `direction` is None wherever no source I could read states the sign.
LAWS = {
    "oleic-anhydride": {
        "name": "Marangoni flow driven by interfacial hydrolysis (onboard fuel)",
        "configuration": "onboard fuel: the droplet carries its own precursor and moves "
                         "for minutes to hours, until the fuel runs out or waste "
                         "accumulates",
        "steps": [
            "oleic anhydride meets water at the oil-water interface and hydrolyses to "
            "two oleic acids, which deprotonate to oleate and protons",
            "the reaction is not even over the surface, so the local pH is not even",
            "interfacial tension rises as local pH falls (maximum at pH 9), so the "
            "tension is not even either",
            "the tension gradient is a Marangoni force tangential to the surface, which "
            "drives Marangoni flow; once the initial symmetry breaks by fluctuation a "
            "pair of convective vortices organises inside the droplet",
            "the convective flow brings fresh precursor to one pole and releases product "
            "at the other, which sustains the imbalance instead of equilibrating it",
            "the droplet moves along the axis of the convection, leaving a low-pH trail "
            "of expelled surfactant behind it",
        ],
        "sensor": "the oil-water interface itself",
        "motor": "the convective flow structure inside the droplet",
        "chemotaxis": {
            "stimulus": "pH gradient, self-generated or externally imposed",
            "direction": "up the gradient, toward the highest pH",
            "source": "Banzhaf & Yamamoto (2015) section 19.2.5; Hanczyc & Ikegami (2010) "
                      "section 3.1",
        },
        "tension": dict(TENSION),
        "refs": ["[365]", "[364]", "[363]"],
    },
    "fuel-surfactant": {
        "name": "Marangoni flow driven by catalysed hydrolysis of an exogenous fuel",
        "configuration": "onboard catalyst: the fuel is supplied by the environment and "
                         "the droplet is not consumed, so in principle the motion can be "
                         "sustained indefinitely",
        "steps": [
            "the droplet of 4-octylaniline carries 5 mol % of an amphiphilic catalyst",
            "the catalyst hydrolyses the amphiphilic precursor dispersed in the aqueous "
            "phase into 4-octylaniline",
            "tiny oil droplets form on the surface, are conveyed to the posterior surface "
            "and are released into the aqueous solution",
            "the interfacial tension around the droplet is left asymmetric, which is the "
            "authors' hypothesis for the unidirectional motion",
        ],
        "chemotaxis": {"stimulus": None, "direction": None,
                       "source": "not reported in the accessible abstract"},
        "refs": ["[858]"],
    },
    "maze-chemotaxis": {
        "name": "Marangoni flow driven by acid/base chemistry at the droplet surface",
        "configuration": "no onboard reaction fuel: the droplet emits surface-active "
                         "chemicals and follows an externally imposed pH gradient",
        "steps": [
            "the dichloromethane droplet carries 2-hexyldecanoic acid and sits in an "
            "alkaline aqueous phase, where the acid deprotonates to the surface-active "
            "carboxylate",
            "acid placed at one exit of the maze re-protonates the carboxylate nearby, "
            "which raises the interfacial tension on that side",
            "the resulting tension gradient propels the droplet, and it finds the "
            "shortest path through the maze",
        ],
        "chemotaxis": {
            "stimulus": "pH gradient",
            "direction": "down the gradient, toward the low-pH region",
            "source": "Lagzi, Soh, Wesson, Browne & Grzybowski (2010), abstract",
        },
        "refs": ["[483]"],
    },
    "decanol-salt": {
        "name": "externally imposed tension gradient, with no reactive chemistry",
        "configuration": "no reactive chemistry: the droplet moves as long as the "
                         "external gradient lasts, and the reaction list is empty",
        "steps": [
            "decanol droplets sit in an aqueous solution of sodium decanoate",
            "a NaCl concentration gradient makes the interfacial tension uneven around "
            "the droplet, with no chemical transformation involved",
            "the droplet migrates in the gradient after an induction time; induction time "
            "and migration velocity both depend on the decanoate and NaCl concentrations",
        ],
        "chemotaxis": {
            "stimulus": "NaCl concentration gradient",
            "direction": None,
            "source": "Cejkova, Novak, Stepanek & Hanczyc (2014); the accessible abstract "
                      "reports induction time and migration velocity but not the sign of "
                      "the response, so no direction is asserted here",
        },
        "demonstrated": [
            "migration in a linear chemotactic assay and in a topologically complex "
            "environment",
            "reversing the direction of movement repeatedly",
            "carrying and releasing a chemically reactive cargo",
            "selecting the stronger of two concentration gradients",
            "starting chemotaxis from an external temperature stimulus",
        ],
        "refs": ["[885]"],
    },
}

#: The published assay vessels, as reported in the methods sections.
VESSELS = {
    "oleic-anhydride": [
        {"what": "glass slide with a concave depression", "diameter_mm": 15.0,
         "aqueous_uL": 100.0, "source": "Hanczyc & Ikegami (2010), appendix A.3"},
        {"what": "glass dish with a quartz base", "diameter_mm": 35.0,
         "base_diameter_mm": 27.0, "aqueous_uL": 800.0,
         "source": "Hanczyc & Ikegami (2010), appendix A.4; the 27 mm base is the "
                   "'diameter of the liquid recipient' in the book's figure 19.11"},
        {"what": "glass dish, behavioural modes", "diameter_mm": 100.0,
         "aqueous_uL": 20000.0, "source": "Horibe et al. (2011), experimental"},
        {"what": "glass dish with a triangular glass wall, restitution",
         "diameter_mm": 60.0, "wall_length_mm": 48.0, "aqueous_uL": 5000.0,
         "source": "Horibe et al. (2011), experimental"},
        {"what": "glass dish with a glass base, two-droplet collective behaviour",
         "diameter_mm": 35.0, "base_diameter_mm": 27.0, "aqueous_uL": 2000.0,
         "source": "Horibe et al. (2011), experimental"},
    ],
    "fuel-surfactant": [],
    "maze-chemotaxis": [],
    "decanol-salt": [],
}

#: How many surfactant molecules one precursor molecule yields.
SURFACTANT_PER_PRECURSOR = {
    "oleic-anhydride": 2,       # one anhydride -> two oleic acids -> two oleates
    "fuel-surfactant": 1,
    "maze-chemotaxis": 1,
    "decanol-salt": 0,
}

PHASES = {
    "oleic-anhydride": {
        "oil": ["nitrobenzene", "oleic_anhydride"],
        "aqueous": ["H2O", "oleate", "H"],
        "interface": ["oleic_acid"],
    },
    "fuel-surfactant": {
        "oil": ["octylaniline"],
        "aqueous": ["H2O", "precursor", "byproduct"],
        "interface": ["catalyst"],
    },
    "maze-chemotaxis": {
        "oil": ["dichloromethane", "hexyldecanoic_acid"],
        "aqueous": ["H2O", "OH", "H", "hexyldecanoate"],
        "interface": [],
    },
    "decanol-salt": {
        "oil": ["decanol"],
        "aqueous": ["H2O", "decanoate", "Na", "Cl"],
        "interface": [],
    },
}


# ---------------------------------------------------------------------------
def generate(p, rng):
    if p.system not in SYSTEMS:
        raise ValueError(f"system must be one of {list(SYSTEMS)}, got {p.system!r}")

    species, reactions, conservation, initial = _BUILD[p.system](p)
    geometry = droplet_geometry(p.droplet_volume_uL)
    per = SURFACTANT_PER_PRECURSOR[p.system]

    analysis = {
        "system": p.system,
        "droplet": geometry,
        "concentration_units": "mM",
        "pH": float(p.pH),
        "surfactant_mM": float(p.surfactant_mM),
        "reactive": bool(reactions),
        "motion_is_simulated": False,
        "why_no_trajectory": (
            "the propulsion mechanism is published qualitatively (a Marangoni force "
            "set by the interfacial-tension gradient) but no accessible source gives "
            "the force law or a droplet speed, so Chemart records the mechanism and "
            "the measured behaviour instead of integrating a trajectory"
        ),
    }
    if per:
        analysis["fuel_budget"] = fuel_budget(p.droplet_volume_uL, p.precursor_M, per)

    if p.system == "oleic-anhydride":
        analysis["tension"] = dict(TENSION)
        # Published threshold: tension peaks at pH 9, so a bulk phase above it sits
        # on the low-tension side and only a local acidification moves the droplet.
        analysis["above_tension_maximum"] = bool(p.pH > TENSION["maximum_at_pH"])
        analysis["overall_reaction"] = (
            "oleic_anhydride + H2O -> 2 oleate + 2 H (hydrolysis followed by "
            "deprotonation of both oleic acids)"
        )
        analysis["measured_behaviour"] = dict(HORIBE)
        analysis["sizes_and_timings"] = dict(ALIFE2010)
        analysis["modes_at_this_volume"] = HORIBE["modes_by_volume_uL"].get(
            f"{p.droplet_volume_uL:g}"
        )
        analysis["size_series"] = [
            {**droplet_geometry(v),
             "modes": HORIBE["modes_by_volume_uL"].get(f"{v:g}"),
             "collective_attraction": (
                 True if v in HORIBE["collective_attraction_uL"]
                 else False if v in HORIBE["no_collective_attraction_uL"] else None),
             **fuel_budget(v, p.precursor_M, per)}
            for v in HORIBE["volumes_uL"]
        ]

    extras = {
        "interaction_law": dict(LAWS[p.system]),
        "analysis": analysis,
        "space": {
            "kind": "continuous",
            "dimensions": 2,
            "units": "mm",
            "droplet": geometry,
            "vessels": [dict(v) for v in VESSELS[p.system]],
            "positions": (
                "not simulated: this entry returns the published interfacial chemistry "
                "and the published assay geometry, not a droplet trajectory"
            ),
        },
        "phases": dict(PHASES[p.system]),
        "buffered": ["H2O"] if any(s.id == "H2O" for s in species) else [],
        "units": {
            "initial_state": "mM; oil-phase and aqueous-phase species are listed "
                             "together, see extras['phases'] for which is which",
            "droplet_volume": "uL (1 uL = 1 mm^3)",
        },
    }
    if conservation:
        extras["conservation"] = conservation
    return network(reactions, species=species, initial_state=initial,
                   status="complete", extras=extras)
