"""Swarm Chemistry (Sayama, 2009-2012; book 11.4.2). Catalog id: swarm-chemistry.

A species here is a kinetic recipe -- the eight-parameter tuple
(R, V_normal, V_max, c1..c5) of book table 11.6 -- and particles are never
created or destroyed, so the network carries NO reactions: the whole
mechanism is the Boid-like force law in ``extras["interaction_law"]``, the
simulated configuration is in ``extras["space"]`` and the order parameters
measured from the run are in ``extras["analysis"]``.

The only transformation Sayama defines is at the level of recipes: in the
evolutionary simulator a particle that loses a collision adopts its
neighbour's recipe. With ``recipe_transmission`` that variant is run and the
observed adoption events become catalytic reactions ``A + B -> 2 B``.

Update law and constants transcribed from Sayama's own source release
(SwarmChemistry-1.3.0-src, SwarmPopulationSimulator.simulateSwarmBehavior and
SwarmIndividual.accelerate) and from the evolutionary simulator (2015).
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from chemart.helpers.params import apportion
from chemart.network import Network, Reaction, Species

#: Recipe fields, in the order Sayama's recipe text uses them.
FIELDS = ("R", "V_normal", "V_max", "c1", "c2", "c3", "c4", "c5")

#: Ranges of book table 11.6 (= the maxima in SwarmParameters.java).
BOUNDS = {
    "R": (0.0, 300.0), "V_normal": (0.0, 20.0), "V_max": (0.0, 40.0),
    "c1": (0.0, 1.0), "c2": (0.0, 1.0), "c3": (0.0, 100.0),
    "c4": (0.0, 0.5), "c5": (0.0, 1.0),
}

#: Constants of the released simulator.
INITIAL_SPEED = 5.0        # Recipe.createPopulation: velocity components in [-5, 5)
WHIM = 5.0                 # random steering: acceleration components in [-5, 5)
STRAY = 0.5                # lone particle: acceleration components in [-0.5, 0.5)
MIN_D2 = 0.001             # "if (d == 0) d = 0.001" in separation and pace keeping
TRANSMISSION_RADIUS = 10.0  # collision range of recipe transmission (2011 variant)
METRIC_RADIUS = 30.0       # radius used for the segregation and cluster measures
COMPETITIONS = ("faster", "slower", "behind", "majority-relative")

#: Sayama's published sample recipes (bingweb.binghamton.edu/~sayama/SwarmChemistry/),
#: each a list of (count, R, V_normal, V_max, c1, c2, c3, c4, c5).
RECIPES: dict[str, list[tuple[float, ...]]] = {
    "blobs": [(300, 20.8, 1.95, 20.75, 0.95, 0.99, 9.31, 0.05, 0.68)],
    "linear-oscillator": [
        (133, 214.41, 17.93, 35.14, 0.64, 0.13, 0.29, 0.08, 0.97),
        (24, 253.6, 7.19, 15.51, 0.82, 0.33, 32.65, 0.34, 0.56),
    ],
    "turbulent-runner": [
        (131, 177.1, 9.71, 30.06, 0.8, 0.43, 19.65, 0.45, 0.91),
        (169, 277.3, 14.67, 37.71, 0.68, 0.23, 77.01, 0.02, 0.31),
    ],
    "playing-catch": [
        (76, 84.06, 0.09, 9.89, 0.33, 0.32, 15.66, 0.22, 0.68),
        (100, 158.86, 18.4, 24.98, 0.3, 0.3, 1.72, 0.06, 0.37),
    ],
    "recombining-blobs": [
        (132, 45.91, 10.82, 21.11, 0.86, 0.13, 42.48, 0.32, 0.74),
        (84, 113.26, 3.41, 25.71, 0.4, 0.39, 49.53, 0.13, 0.24),
    ],
    "wedding-ring": [
        (24, 220.51, 13.88, 3.47, 0.46, 0.38, 6.23, 0.19, 0.68),
        (13, 64.07, 1.4, 19.7, 0.88, 0.27, 0.36, 0.47, 0.72),
        (35, 117.53, 7.31, 21.72, 0.3, 0.5, 98.69, 0.03, 0.29),
    ],
    "rotary": [
        (29, 122.13, 19.19, 17.98, 0.65, 0.44, 19.88, 0.46, 0.2),
        (51, 299.13, 0.79, 38.71, 0.25, 0.18, 86.49, 0.38, 0.43),
        (10, 252.92, 19.99, 10.21, 0.23, 0.17, 1.22, 0.28, 0.92),
    ],
    "swinger": [
        (48, 150.39, 15.89, 23.54, 0.74, 0.45, 62.65, 0.33, 0.13),
        (152, 217.14, 12.13, 12.42, 0.59, 0.98, 14.06, 0.04, 0.65),
        (14, 248.54, 5.85, 22.26, 0.43, 0.11, 17.14, 0.06, 0.68),
        (31, 141.53, 2.91, 4.86, 0.92, 0.03, 21.87, 0.28, 0.2),
    ],
    "jelly-fish": [
        (134, 262.65, 12.01, 25.87, 0.97, 1.0, 56.35, 0.26, 0.61),
        (67, 288.17, 6.19, 23.37, 0.95, 1.0, 1.31, 0.1, 0.9),
        (68, 150.5, 12.97, 15.87, 0.46, 0.39, 57.95, 0.17, 0.48),
    ],
    "no-wait-this-way": [
        (60, 262.68, 2.82, 38.32, 0.21, 0.01, 54.93, 0.11, 0.19),
        (40, 78.58, 5.7, 33.23, 0.89, 0.18, 45.44, 0.04, 0.05),
        (40, 257.27, 14.96, 35.66, 0.2, 0.8, 47.81, 0.13, 0.13),
    ],
    "pulsating-eye": [
        (102, 293.86, 17.06, 38.3, 0.81, 0.05, 0.83, 0.2, 0.9),
        (124, 226.18, 19.27, 24.57, 0.95, 0.84, 13.09, 0.07, 0.8),
        (74, 49.98, 8.44, 4.39, 0.92, 0.14, 96.92, 0.13, 0.51),
    ],
    "chaos-cells": [
        (144, 109.03, 6.71, 12.7, 0.47, 0.6, 61.43, 0.02, 0.21),
        (89, 117.15, 16.33, 31.88, 0.39, 0.13, 12.96, 0.48, 0.8),
        (67, 76.3, 8.59, 26.57, 0.7, 0.64, 28.39, 0.3, 0.35),
    ],
    "fast-walker-slow-follower": [
        (67, 216.35, 11.75, 7.7, 0.83, 0.97, 97.31, 0.02, 0.38),
        (29, 254.64, 7.28, 7.0, 0.95, 0.11, 22.41, 0.43, 0.31),
        (13, 105.4, 3.55, 5.24, 0.34, 0.18, 23.53, 0.39, 0.24),
    ],
    "aggressive-predator": [
        (18, 211.92, 12.59, 19.37, 0.09, 0.21, 57.92, 0.0, 0.95),
        (41, 257.27, 14.96, 35.66, 0.2, 0.8, 47.81, 0.13, 0.13),
        (35, 262.68, 2.82, 38.32, 0.21, 0.01, 54.93, 0.11, 0.19),
        (31, 78.58, 5.7, 33.23, 0.89, 0.18, 45.44, 0.04, 0.05),
        (7, 194.21, 12.88, 21.68, 0.97, 0.19, 99.21, 0.5, 0.13),
    ],
    "cell-with-two-nuclei": [
        (41, 249.84, 4.85, 28.73, 0.34, 0.45, 14.44, 0.09, 0.82),
        (26, 277.87, 15.02, 35.48, 0.68, 0.05, 82.96, 0.46, 0.9),
        (30, 277.87, 15.02, 24.44, 0.68, 0.05, 82.96, 0.43, 0.9),
        (28, 110.8, 16.12, 38.6, 0.18, 0.34, 14.3, 0.01, 0.01),
        (48, 83.79, 13.29, 7.54, 0.08, 0.79, 1.07, 0.15, 0.45),
        (74, 269.64, 6.62, 34.69, 0.36, 0.5, 30.2, 0.03, 0.23),
    ],
    "multicellularity": [
        (99, 19.8, 15.73, 2.61, 0.85, 0.64, 10.51, 0.17, 0.06),
        (48, 300.0, 14.63, 0.0, 0.48, 0.81, 90.27, 0.25, 0.78),
        (37, 275.18, 16.9, 7.05, 0.48, 0.81, 90.27, 0.17, 0.85),
        (8, 159.59, 2.09, 24.19, 0.96, 0.59, 76.03, 0.01, 0.07),
        (42, 73.07, 1.82, 2.36, 0.27, 0.61, 40.55, 0.22, 0.86),
    ],
    "insurmountable-wall": [
        (42, 52.57, 9.91, 20.42, 0.32, 0.76, 1.8, 0.01, 0.64),
        (25, 84.87, 8.82, 24.98, 0.91, 0.44, 40.97, 0.18, 0.6),
        (45, 220.42, 4.65, 7.53, 0.96, 0.35, 46.18, 0.25, 1.0),
        (49, 279.64, 10.29, 35.95, 0.37, 0.49, 38.09, 0.32, 0.89),
    ],
    # Erskine & Herrmann (ECAL 2013) table 1, with their 300:50 population mix:
    # two species that each form a single blob alone but divide like a cell together.
    "cell-division": [
        (300, 20.5, 1.94, 20.7, 1.0, 1.0, 18.6, 0.05, 1.0),
        (50, 300.0, 15.58, 37.08, 1.0, 0.05, 9.11, 0.47, 0.61),
    ],
    "uzushio": [
        (49, 184.07, 12.64, 38.49, 0.67, 0.31, 63.78, 0.11, 0.21),
        (67, 177.92, 15.66, 11.05, 0.85, 0.41, 51.77, 0.09, 0.86),
        (58, 177.92, 15.66, 11.05, 0.85, 0.41, 51.77, 0.03, 0.86),
        (10, 298.48, 0.65, 19.27, 0.84, 0.11, 83.18, 0.24, 0.5),
        (37, 81.6, 18.34, 27.3, 0.94, 0.17, 3.41, 0.2, 0.9),
        (1, 274.4, 3.33, 7.42, 0.9, 0.19, 7.37, 0.43, 0.17),
        (12, 203.64, 3.33, 7.42, 0.9, 0.19, 7.37, 0.43, 0.17),
    ],
}

LAW = [
    "neighbourhood: N_i = {j != i : |x_j - x_i| < R_i}",
    "if N_i is empty: a_i = u, u componentwise uniform in [-0.5, 0.5) (stray randomly)",
    "else: a_i = c1_i (<x>_N - x_i) + c2_i (<v>_N - v_i) "
    "+ c3_i sum_{j in N_i} (x_i - x_j) / |x_i - x_j|^2, "
    "with |x_i - x_j|^2 replaced by 0.001 when it is 0",
    "with probability c4_i: a_i += w, w componentwise uniform in [-5, 5) (steer randomly)",
    "v_i <- v_i + a_i, then rescaled to |v_i| = V_max_i if it exceeds it",
    "pace keeping: v_i <- v_i + c5_i (V_normal_i - |v_i|) v_i / |v_i|, "
    "then rescaled to |v_i| = V_max_i if it exceeds it (|v_i| = 0.001 when it is 0)",
    "x_i <- x_i + v_i; all particles are updated synchronously, one step at a time",
]


# ---------------------------------------------------------------------------
# the kinetic law
# ---------------------------------------------------------------------------
def step(pos, vel, par, rng):
    """One synchronous time step; `par` is (N, 8) in the order of FIELDS."""
    n_particles, dim = pos.shape
    R, v_normal, v_max, c1, c2, c3, c4, c5 = (par[:, i] for i in range(8))

    diff = pos[:, None, :] - pos[None, :, :]                 # x_i - x_j
    d2 = (diff * diff).sum(-1)
    near = d2 < (R * R)[:, None]
    np.fill_diagonal(near, False)
    n = near.sum(1)
    has = n > 0
    weight = near / np.maximum(n, 1)[:, None]

    centre = weight @ pos
    mean_v = weight @ vel
    inv = np.where(near, 1.0 / np.where(d2 == 0.0, MIN_D2, d2), 0.0)
    separation = (inv[:, :, None] * diff).sum(1)

    acc = (c1[:, None] * (centre - pos) + c2[:, None] * (mean_v - vel)
           + c3[:, None] * separation)
    stray = (rng.random((n_particles, dim)) - 0.5) * (2 * STRAY)
    whim = (rng.random((n_particles, dim)) - 0.5) * (2 * WHIM)
    steer = (rng.random(n_particles) < c4) & has
    acc = np.where(steer[:, None], acc + whim, acc)
    acc = np.where(has[:, None], acc, stray)

    new_v = _limit(vel + acc, v_max)
    speed = np.linalg.norm(new_v, axis=1)
    speed = np.where(speed == 0.0, MIN_D2, speed)
    new_v = _limit(new_v + new_v * ((v_normal - speed) / speed * c5)[:, None], v_max)
    return pos + new_v, new_v


def _limit(vel, v_max):
    """SwarmIndividual.accelerate: rescale to V_max only when it is exceeded."""
    speed = np.linalg.norm(vel, axis=1)
    over = speed > v_max
    factor = np.where(over, v_max / np.where(speed == 0.0, 1.0, speed), 1.0)
    return vel * factor[:, None]


def _losing(pos, vel, types, par, defender, attacker, mode):
    """Sayama's competition functions: does `defender` lose to `attacker`?"""
    dv, av = vel[defender], vel[attacker]
    if mode == "faster":
        return ~((dv * dv).sum(1) > (av * av).sum(1))
    if mode == "slower":
        return ~((dv * dv).sum(1) < (av * av).sum(1))
    if mode == "behind":
        threshold = np.cos(0.75 * np.pi)
        a = pos[attacker] - pos[defender]
        dot = (a * dv).sum(1)
        return ~(dot > threshold * np.linalg.norm(a, axis=1) * np.linalg.norm(dv, axis=1))
    if mode == "majority-relative":
        radius = np.maximum(METRIC_RADIUS, par[:, 0])
        d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
        near = d2 < (radius * radius)[:, None]
        np.fill_diagonal(near, False)
        same = near & (types[:, None] == types[None, :])
        total = near.sum(1)
        share = np.where(total > 0, same.sum(1) / np.maximum(total, 1), 0.0)
        return share[defender] < share[attacker]
    raise ValueError(f"unknown competition function {mode!r}; use one of {list(COMPETITIONS)}")


def _transmit(pos, vel, types, rows, mode):
    """One round of recipe transmission by collision (evolutionary variant)."""
    d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    other = types[:, None] != types[None, :]
    candidates = np.where((d2 < TRANSMISSION_RADIUS ** 2) & other, d2, np.inf)
    nearest = candidates.argmin(1)
    touching = np.isfinite(candidates.min(1))
    index = np.arange(len(types))
    losing = _losing(pos, vel, types, rows[types], index, nearest, mode)
    changed = touching & losing
    events = Counter(
        (int(types[i]), int(types[nearest[i]])) for i in np.flatnonzero(changed)
    )
    types = types.copy()
    types[changed] = types[nearest[changed]]
    return types, events


# ---------------------------------------------------------------------------
# a run and its order parameters
# ---------------------------------------------------------------------------
def run(rows, counts, steps, dim, extent, rng, transmission=None):
    """Simulate a recipe; return the final state, the metrics and the events."""
    rows = np.asarray(rows, dtype=float)
    types = np.repeat(np.arange(len(counts)), counts)
    n_particles = len(types)
    pos = rng.random((n_particles, dim)) * extent
    vel = (rng.random((n_particles, dim)) - 0.5) * (2 * INITIAL_SPEED)

    window = max(1, steps // 10)
    samples: list[dict] = []
    events: Counter = Counter()
    for t in range(steps):
        if transmission:
            types, fired = _transmit(pos, vel, types, rows, transmission)
            events.update(fired)
        pos, vel = step(pos, vel, rows[types], rng)
        if t >= steps - window:
            samples.append(_metrics(pos, vel, types, rows))
    if not samples:                                    # steps == 0
        samples.append(_metrics(pos, vel, types, rows))

    analysis = {k: float(np.mean([s[k] for s in samples])) for k in samples[0]}
    analysis.update(_structure(pos, types, len(rows)))
    analysis["speed"] = _speeds(vel, types, len(rows))
    analysis["steps"] = int(steps)
    analysis["averaging_window"] = int(len(samples))
    return pos, vel, types, analysis, events


def _metrics(pos, vel, types, rows):
    speed = np.linalg.norm(vel, axis=1)
    moving = speed > 0
    direction = np.zeros_like(vel)
    direction[moving] = vel[moving] / speed[moving, None]
    polarization = float(np.linalg.norm(direction.mean(0)))

    radius = pos - pos.mean(0)
    length = np.linalg.norm(radius, axis=1)
    off = length > 0
    unit = np.zeros_like(radius)
    unit[off] = radius[off] / length[off, None]
    if pos.shape[1] == 2:
        spin = unit[:, 0] * direction[:, 1] - unit[:, 1] * direction[:, 0]
        rotation = float(spin.mean())
    else:
        rotation = float(np.linalg.norm(np.cross(unit, direction).mean(0)))

    return {
        "polarization": polarization,
        "rotation": rotation,
        "mean_speed": float(speed.mean()),
        "radius_of_gyration": float(np.sqrt((length ** 2).mean())),
        "speed_ratio": float(np.mean(speed / np.maximum(rows[types][:, 1], 1e-12))),
        **_mixing(pos, types),
        **_entropy(pos),
    }


def _entropy(pos):
    """Erskine & Herrmann's spatial entropy: patches of 0.1 of the swarm's extent.

    H = -sum_k P(k) log P(k) over the 10^d patches of the swarm's bounding box;
    D_KL against an evenly dispersed swarm is then d log 10 - H, and grows when
    the swarm breaks into separate clumps.
    """
    dim = pos.shape[1]
    low, high = pos.min(0), pos.max(0)
    span = np.where(high > low, high - low, 1.0)
    cell = np.clip(((pos - low) / span * 10).astype(int), 0, 9)
    counts = np.bincount(np.ravel_multi_index(cell.T, (10,) * dim), minlength=10 ** dim)
    share = counts[counts > 0] / len(pos)
    h = float(-(share * np.log(share)).sum())
    return {"spatial_entropy": h, "kl_divergence": float(dim * np.log(10.0) - h)}


def _mixing(pos, types):
    d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    near = d2 < METRIC_RADIUS ** 2
    total = near.sum(1)
    same = (near & (types[:, None] == types[None, :])).sum(1)
    has = total > 0
    observed = float((same[has] / total[has]).mean()) if has.any() else 1.0
    shares = np.array([np.mean(types == t) for t in np.unique(types)])
    expected = float((shares ** 2).sum())
    index = (observed - expected) / (1.0 - expected) if expected < 1.0 else 0.0
    return {
        "same_type_neighbour_fraction": observed,
        "mixed_expectation": expected,
        "segregation_index": float(index),
    }


def _structure(pos, types, n_types):
    from scipy.sparse import csr_array
    from scipy.sparse.csgraph import connected_components

    d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
    near = d2 < METRIC_RADIUS ** 2
    np.fill_diagonal(near, False)
    count, labels = connected_components(csr_array(near), directed=False)
    sizes = np.bincount(labels)
    radius = np.linalg.norm(pos - pos.mean(0), axis=1)
    out = {
        "cluster_count": int(count),
        "largest_cluster_fraction": float(sizes.max() / len(types)),
        "mean_radius": [float(radius[types == t].mean()) if (types == t).any() else 0.0
                        for t in range(n_types)],
    }
    for t in range(n_types):
        members = types == t
        if members.sum() > 1:
            sub = np.bincount(labels[members])
            out.setdefault("clusters_per_species", []).append(int((sub > 0).sum()))
        else:
            out.setdefault("clusters_per_species", []).append(int(members.sum()))
    return out


def _speeds(vel, types, n_types):
    speed = np.linalg.norm(vel, axis=1)
    return [float(speed[types == t].mean()) if (types == t).any() else 0.0
            for t in range(n_types)]


# ---------------------------------------------------------------------------
# the catalogued generator
# ---------------------------------------------------------------------------
def _rows(p):
    if p.recipe == "custom":
        if not p.custom_recipe:
            raise ValueError("recipe='custom' needs custom_recipe: a list of "
                             "[count, R, V_normal, V_max, c1, c2, c3, c4, c5] rows")
        rows = []
        for i, row in enumerate(p.custom_recipe):
            if (not isinstance(row, (list, tuple)) or len(row) != 9
                    or any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in row)):
                raise ValueError(
                    f"custom_recipe[{i}] must be 9 numbers "
                    f"[count, R, V_normal, V_max, c1, c2, c3, c4, c5], got {row!r}")
            if row[0] < 1:
                raise ValueError(f"custom_recipe[{i}]: count must be at least 1, got {row[0]!r}")
            for name, value in zip(FIELDS, row[1:]):
                low, high = BOUNDS[name]
                if not low <= value <= high:
                    raise ValueError(
                        f"custom_recipe[{i}]: {name}={value!r} is outside the published "
                        f"range [{low}, {high}] (book table 11.6)")
            rows.append([float(v) for v in row])
        return rows
    return [[float(v) for v in row] for row in RECIPES[p.recipe]]


def _counts(rows, total):
    """Rescale the recipe's counts to `total`, keeping at least one of each type."""
    if total < len(rows):
        raise ValueError(f"particles={total} is fewer than the {len(rows)} ingredients "
                         f"of the recipe; use particles >= {len(rows)}")
    weights = [row[0] - 1 for row in rows]
    if sum(weights) <= 0:
        weights = [1.0] * len(rows)
    return [n + 1 for n in apportion(total - len(rows), weights)]


def _extent(p):
    """`space` = 0 keeps the density of Sayama's 300 particles in a 300 px box.

    The recipes are tuned to that density: a particle whose perception radius is
    20 px finds no neighbours at all in a sparser box, and the swarm never
    aggregates.
    """
    if p.space > 0:
        return float(p.space)
    return 300.0 * (p.particles / 300.0) ** (1.0 / p.dimensions)


def generate(p, rng):
    rows = _rows(p)
    counts = _counts(rows, p.particles)
    extent = _extent(p)
    slug = p.recipe
    ids = [f"{slug}-{i + 1}" for i in range(len(rows))]
    species = [
        Species(sid, structure="{} * ({})".format(
            n, ", ".join(_short(v) for v in row[1:])))
        for sid, n, row in zip(ids, counts, rows)
    ]

    mode = p.competition if p.recipe_transmission else None
    pos, vel, types, analysis, events = run(
        [row[1:] for row in rows], counts, p.steps, p.dimensions, extent, rng,
        transmission=mode,
    )

    reactions = [
        Reaction({ids[a]: 1, ids[b]: 1}, {ids[b]: 2}, count=int(n))
        for (a, b), n in sorted(events.items())
    ]
    final = Counter(int(t) for t in types)

    analysis["normal_speed"] = [row[2] for row in rows]
    analysis["final_counts"] = [final.get(i, 0) for i in range(len(rows))]
    analysis["definitions"] = {
        "polarization": "|mean of v_i/|v_i||, 1 = one common heading",
        "rotation": "mean of (x_i - <x>)^ x v_i^ about the centre of mass "
                    "(signed in 2D, magnitude of the mean vector in 3D)",
        "segregation_index": f"(f - f0) / (1 - f0), f = mean fraction of same-species "
                             f"particles within {METRIC_RADIUS} px, f0 = sum of squared "
                             "species shares (0 = well mixed, 1 = fully segregated)",
        "cluster_count": f"connected components of the {METRIC_RADIUS} px proximity graph",
        "mean_radius": "per species, the mean distance from the swarm's centre of mass at the "
                       "final step; a core-shell structure shows as two different radii",
        "spatial_entropy": "H = -sum_k P(k) log P(k) over the 10^d patches of the swarm's "
                           "bounding box (Erskine & Herrmann 2013)",
        "kl_divergence": "d log 10 - H, the divergence from an evenly dispersed swarm; "
                         "it grows when the swarm breaks into separate clumps",
        "speed_ratio": "mean |v_i| / V_normal_i (1 = pace keeping satisfied)",
        "averaging": "scalar order parameters are averaged over the last 10% of the steps; "
                     "the structural ones are measured at the final step",
    }

    return Network(
        species=species,
        reactions=reactions,
        status="observed" if mode else "complete",
        initial_state={sid: float(n) for sid, n in zip(ids, counts)},
        extras={
            "space": {
                "kind": "continuous",
                "dimensions": int(p.dimensions),
                "units": "pixels",
                "initial_domain": [[0.0, round(extent, 4)]] * p.dimensions,
                "boundary": "unbounded: no walls and no wrapping",
                "steps": int(p.steps),
                "particles": [
                    {"species": ids[int(t)],
                     "position": [round(float(v), 4) for v in x],
                     "velocity": [round(float(v), 4) for v in u]}
                    for t, x, u in zip(types, pos, vel)
                ],
            },
            "interaction_law": {
                "name": "Sayama swarm chemistry (Boid-like kinetic interaction)",
                "state": "each particle has a position x_i and a velocity v_i; its "
                         "species fixes the eight kinetic parameters",
                "parameters": {
                    sid: dict(zip(FIELDS, [float(v) for v in row[1:]]))
                    for sid, row in zip(ids, rows)
                },
                "units": {"R": "pixels", "V_normal": "pixels/step", "V_max": "pixels/step",
                          "c1": "1/step^2", "c2": "1/step", "c3": "pixels^2/step^2",
                          "c4": "probability per step", "c5": "dimensionless"},
                "bounds": {k: list(v) for k, v in BOUNDS.items()},
                "steps": list(LAW),
                "initial_condition": f"positions uniform in [0, {round(extent, 4)})^{p.dimensions}, "
                                     f"velocity components uniform in [-{INITIAL_SPEED}, {INITIAL_SPEED})",
                "transmission": _transmission_law(mode),
                "source": "H. Sayama, Swarm Chemistry, Artificial Life 15(1):105-114 (2009); "
                          "equations as released in SwarmChemistry-1.3.0-src",
            },
            "analysis": analysis,
        },
    )


def _transmission_law(mode):
    if not mode:
        return "off: particles never change their recipe, so the network has no reactions"
    return (
        f"on ({mode}): a particle that has a particle of another species within "
        f"{TRANSMISSION_RADIUS} px adopts the species of the nearest such particle "
        "when it loses the competition; the observed adoptions are the reactions "
        "A + B -> 2 B with their firing counts"
    )


def _short(value):
    """Sayama's recipe text rounds every parameter to two decimals."""
    return f"{round(float(value), 2):g}"
