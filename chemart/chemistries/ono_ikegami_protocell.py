"""Ono & Ikegami autopoietic protocells (book 6.3.2). Catalog id: ono-ikegami-protocell.

Five species live on a hexagonal torus, one particle per cell: A (autocatalyst,
hydrophilic), M (membrane, hydrophobic; isotropic ``M_i`` or anisotropic ``M_a``
carrying one of the six lattice orientations), X (food, neutral), Y (waste,
neutral) and W (water, hydrophilic). The reaction set is the book's:

    A + X -> 2 A            autocatalytic replication
    A + X -> A + M          membrane production, catalysed by A
    A -> Y, M -> Y, X -> Y  every particle but water decays to waste
    Y -> X                  recycling by an external energy source

Two faces. ``generate`` returns exactly that network. ``evolve`` runs the
lattice and yields a frame per sweep (the count of each species and the
reactions fired in the sweep); it returns the reactions that actually fired
with their counts, the final lattice in ``extras["space"]``, the interaction
energies in ``extras["energies"]`` and the membrane/protocell measurements of
the final lattice in ``extras["analysis"]``.

Spatial dynamics (reconstructed; see the catalog's `decisions`): hydrophilic and
hydrophobic particles repel each other, neutral particles couple weakly to both,
and an anisotropic M_a spreads that repulsion unevenly over the six directions.
Particles move by Kawasaki exchange with a neighbour and M_a rotates, both
accepted with the Metropolis rule min(1, exp(-dE/T)). Each step runs
`relaxation` exchange passes before one chemistry pass, which is how the
paper's separation between mobility (7e-3) and reaction (1e-4) rates enters.
"""

from __future__ import annotations

from collections import deque

import numpy as np

from chemart.network import Network, Reaction, Species
from chemart.trajectory import Frame

#: Species order used inside the lattice arrays.
A, MEM, X, Y, W = range(5)
CODE = "AMXYW"

#: The six axial directions of the hexagonal lattice, 60 degrees apart.
DIRS = ((0, 1), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1))

#: Published constants of the one-dimensional predecessor, Ono & Ikegami
#: (2000), appendix A. They belong to that model's enzyme-mediated scheme
#: (A + A -> A + E, E + R -> E + A, E + R -> E + M), not to the book's
#: two-dimensional scheme implemented here, so they are reported as a
#: reference rather than used as this network's rate constants.
JTB2000 = {
    "source": "Ono & Ikegami, J. theor. Biol. 206:243-253 (2000), appendix A",
    "P_E": 0.5e-6,
    "P_A_spontaneous": 2e-6,
    "P_R_recycling": 100e-6,
    "P_W_decay": 100e-6,
    "P_D_diffusion": 7e-3,
    "P_R0_repulsion_same_site": 5e-3,
    "P_R1_repulsion_neighbour_site": 2.5e-3,
    "phase_points": {
        "fig4b_periodic_membranes": {"P_A": 9e-6, "P_M": 7e-6},
        "fig7b_region_Ia_cell_dies": {"P_A": 5e-6, "P_M": 3e-6},
        "fig7c_region_Ib_cell_stable": {"P_A": 7e-6, "P_M": 9e-6},
        "fig7d_region_Ic_cell_reproduces": {"P_A": 11e-6, "P_M": 5e-6},
        "fig10_recursive_division": {"P_A": 7e-6, "P_M": 5e-6, "m": 5},
    },
}


# ---------------------------------------------------------------------------
# Lattice helpers
# ---------------------------------------------------------------------------
def neighbour(arr: np.ndarray, k: int) -> np.ndarray:
    """Value held by the direction-k neighbour of every cell."""
    dq, dr = DIRS[k]
    return np.roll(arr, (-dq, -dr), axis=(0, 1))


def anisotropy_field(anisotropy: float, isotropic: bool) -> np.ndarray:
    """F[k, o]: how strongly a membrane particle of orientation o repels in direction k.

    F = 1 + a cos(2 * 60deg * (k - o)), which is nematic (180 degree symmetric),
    so the particle has two equivalent ends. Repulsion is strongest along the
    +-o axis (F = 1 + a) and weakest on the four flanking directions
    (F = 1 - a/2); the mean over the six directions is 1, so M_i and M_a carry
    the same total repulsion and differ only in how it is distributed.

    Past a = 2 the flank term turns negative, which is what makes a membrane
    rather than a droplet: the particle becomes effectively amphiphilic, so it
    wants membrane neighbours along its +-o axis and water on the four flanks.
    That is a one-particle-thick sheet with water on both sides, and it beats
    the compact droplet that a purely repulsive M_i collapses into.
    F[k + 3, o] == F[k, o], so exchanging two particles leaves their mutual
    energy unchanged.
    """
    k = np.arange(6)[:, None]
    o = np.arange(6)[None, :]
    if isotropic:
        return np.ones((6, 6))
    return 1.0 + anisotropy * np.cos(2.0 * np.pi / 3.0 * (k - o))


def _fields(sp, ori, F, nc):
    """Neighbour lookups plus the two potentials the energy is built from.

    hb[k] is the hydrophilic load seen in direction k (neutral particles count
    nc as much); psi is the membrane repulsion felt at a cell by whatever
    hydrophilic particle sits there.
    """
    nb_sp = [neighbour(sp, k) for k in range(6)]
    nb_ori = [neighbour(ori, k) for k in range(6)]
    hb, psi = [], np.zeros(sp.shape)
    for k in range(6):
        ns = nb_sp[k]
        hb.append(((ns == A) | (ns == W)) + nc * ((ns == X) | (ns == Y)))
        psi += (ns == MEM) * F[(k + 3) % 6][nb_ori[k]]
    return nb_sp, nb_ori, hb, psi


def _energy_at(sp_self, ori_self, hb, psi, F, eps, nc) -> np.ndarray:
    """Interaction energy of every cell, given the particle it is assumed to hold."""
    phi = np.zeros(sp_self.shape)
    for k in range(6):
        phi += hb[k] * F[k][ori_self]
    phob = sp_self == MEM
    phil = (sp_self == A) | (sp_self == W)
    neut = (sp_self == X) | (sp_self == Y)
    return eps * (phob * phi + (phil + nc * neut) * psi)


def _components(mask: np.ndarray) -> list[dict]:
    """Connected components of `mask` on the hexagonal torus.

    Each component reports its cells, their unwrapped coordinates, and whether
    it winds around the torus; `wraps` is how the unbounded "outside" is told
    apart from an enclosed compartment.
    """
    H, Wd = mask.shape
    seen = np.full((H, Wd), -1, dtype=np.int64)
    out: list[dict] = []
    for q0 in range(H):
        for r0 in range(Wd):
            if not mask[q0, r0] or seen[q0, r0] >= 0:
                continue
            label = len(out)
            cells: list[tuple[int, int]] = []
            wraps = False
            unwrapped = {(q0, r0): (q0, r0)}
            seen[q0, r0] = label
            queue = deque([(q0, r0)])
            while queue:
                q, r = queue.popleft()
                cells.append((q, r))
                uq, ur = unwrapped[(q, r)]
                for dq, dr in DIRS:
                    nq, nr = (q + dq) % H, (r + dr) % Wd
                    if not mask[nq, nr]:
                        continue
                    target = (uq + dq, ur + dr)
                    if seen[nq, nr] < 0:
                        seen[nq, nr] = label
                        unwrapped[(nq, nr)] = target
                        queue.append((nq, nr))
                    elif unwrapped.get((nq, nr)) != target:
                        wraps = True
            out.append({
                "cells": cells,
                "unwrapped": [unwrapped[c] for c in cells],
                "wraps": wraps,
            })
    return out


def _gyration(unwrapped: list[tuple[int, int]]) -> float:
    """Radius of gyration of a component, in lattice spacings.

    Divided by sqrt(size) this separates shapes: it stays near 0.37 for a
    compact blob at any size and grows like 0.29 sqrt(size) for a filament.
    """
    q = np.array([c[0] for c in unwrapped], dtype=float)
    r = np.array([c[1] for c in unwrapped], dtype=float)
    x, y = q + 0.5 * r, 0.8660254 * r
    return float(np.sqrt(((x - x.mean()) ** 2 + (y - y.mean()) ** 2).mean()))


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
def _hex_distance(dq: np.ndarray, dr: np.ndarray) -> np.ndarray:
    return (np.abs(dq) + np.abs(dq + dr) + np.abs(dr)) // 2


def _initial(p, rng) -> tuple[np.ndarray, np.ndarray]:
    H, Wd = p.height, p.width
    ori = rng.integers(0, 6, size=(H, Wd)).astype(np.int8)
    seeded = p.X_fraction + p.A_fraction + p.M_fraction
    if seeded > 1.0:
        raise ValueError(
            f"X_fraction + A_fraction + M_fraction = {seeded:.3f} exceeds 1; "
            "the remainder of the lattice is water W"
        )
    if p.initial == "random":
        u = rng.random((H, Wd))
        sp = np.full((H, Wd), W, dtype=np.int8)
        sp[u < p.X_fraction] = X
        sp[(u >= p.X_fraction) & (u < p.X_fraction + p.A_fraction)] = A
        sp[(u >= p.X_fraction + p.A_fraction) & (u < seeded)] = MEM
        return sp, ori

    # "cell": a ring of membrane enclosing autocatalyst and food, the paper's
    # section 3.2 experiment lifted to two dimensions.
    if 2 * p.cell_radius + 2 >= min(H, Wd):
        raise ValueError(
            f"cell_radius = {p.cell_radius} does not fit a {H}x{Wd} lattice; "
            f"use cell_radius <= {min(H, Wd) // 2 - 2} or a larger lattice"
        )
    q = np.broadcast_to(np.arange(H)[:, None] - H // 2, (H, Wd))
    r = np.broadcast_to(np.arange(Wd)[None, :] - Wd // 2, (H, Wd))
    dist = _hex_distance(q, r)
    sp = np.full((H, Wd), W, dtype=np.int8)
    sp[(dist > p.cell_radius) & (rng.random((H, Wd)) < p.X_fraction)] = X
    inside = dist < p.cell_radius
    u = rng.random((H, Wd))
    sp[inside] = np.where(u[inside] < 0.5, A, X)
    sp[dist == p.cell_radius] = MEM
    return sp, ori


# ---------------------------------------------------------------------------
# The lattice run
# ---------------------------------------------------------------------------
def _motion_pass(sp, ori, k, phase, F, p, rng):
    """One Metropolis exchange pass along direction k, on one parity sublattice."""
    eps, nc, T = p.repulsion, p.neutral_coupling, p.temperature
    nb_sp, nb_ori, hb, psi = _fields(sp, ori, F, nc)
    partner_sp, partner_ori = nb_sp[k], nb_ori[k]

    e_self = _energy_at(sp, ori, hb, psi, F, eps, nc)
    e_self_new = _energy_at(partner_sp, partner_ori, hb, psi, F, eps, nc)
    hb_n = [neighbour(h, k) for h in hb]
    psi_n = neighbour(psi, k)
    e_partner = neighbour(e_self, k)
    e_partner_new = _energy_at(sp, ori, hb_n, psi_n, F, eps, nc)
    delta = (e_self_new + e_partner_new) - (e_self + e_partner)

    # Pairs must not overlap: take every other row (or column) along the axis
    # the direction moves in, dropping the wrap-around pair on an odd side.
    axis = 0 if DIRS[k][0] else 1
    n = sp.shape[axis]
    coord = np.arange(n)
    coord = coord[:, None] if axis == 0 else coord[None, :]
    mask = (coord % 2 == phase)
    if n % 2:
        mask = mask & (coord != n - 1)
    mask = np.broadcast_to(mask, sp.shape) & (sp != partner_sp)
    if p.mobility_ratio > 1.0:
        mask = mask & ~(((sp == A) | (partner_sp == A))
                        & (rng.random(sp.shape) >= 1.0 / p.mobility_ratio))
    accept = mask & ((delta <= 0) | (rng.random(sp.shape) < np.exp(-delta / T)))

    dq, dr = DIRS[k]
    back = np.roll(accept, (dq, dr), axis=(0, 1))
    sp_new, ori_new = sp.copy(), ori.copy()
    sp_new[accept] = partner_sp[accept]
    ori_new[accept] = partner_ori[accept]
    sp_new[back] = np.roll(sp, (dq, dr), axis=(0, 1))[back]
    ori_new[back] = np.roll(ori, (dq, dr), axis=(0, 1))[back]
    return sp_new, ori_new


def _sweeps(p, rng):
    """Run the lattice: yield (sp, ori, fired) after each sweep, fired being the
    reactions of that sweep by name. The first yield is the initial lattice."""
    sp, ori = _initial(p, rng)
    F = anisotropy_field(p.anisotropy, p.membrane == "isotropic")
    eps, nc, T = p.repulsion, p.neutral_coupling, p.temperature
    yield sp, ori, {}

    for _ in range(p.steps):
        for _ in range(p.relaxation):
            sp, ori = _motion_pass(sp, ori, int(rng.integers(6)),
                                   int(rng.integers(2)), F, p, rng)

        if p.membrane == "anisotropic":
            _, _, hb, psi = _fields(sp, ori, F, nc)
            proposal = rng.integers(0, 6, size=sp.shape).astype(np.int8)
            d = (_energy_at(sp, proposal, hb, psi, F, eps, nc)
                 - _energy_at(sp, ori, hb, psi, F, eps, nc))
            take = (sp == MEM) & ((d <= 0) | (rng.random(sp.shape) < np.exp(-d / T)))
            ori = np.where(take, proposal, ori).astype(np.int8)

        # --- chemistry -------------------------------------------------------
        n_A = sum((neighbour(sp, k) == A) for k in range(6)).astype(float)
        u = rng.random(sp.shape)
        p_rep = p.P_A * n_A
        p_mem = p_rep + p.P_M * n_A
        p_dec = p_mem + p.P_decay

        is_X, is_A, is_M, is_Y = sp == X, sp == A, sp == MEM, sp == Y
        to_A = is_X & (u < p_rep)
        to_M = is_X & (u >= p_rep) & (u < p_mem)
        x_to_Y = is_X & (u >= p_mem) & (u < p_dec)
        a_to_Y = is_A & (u < p.P_decay)
        m_to_Y = is_M & (u < p.P_decay)
        y_to_X = is_Y & (u < p.X_supply)

        fired = {
            "replicate": int(to_A.sum()),
            "membrane": int(to_M.sum()),
            "decay_X": int(x_to_Y.sum()),
            "decay_A": int(a_to_Y.sum()),
            "decay_M": int(m_to_Y.sum()),
            "recycle": int(y_to_X.sum()),
        }

        sp = np.where(to_A, np.int8(A), sp)
        sp = np.where(to_M, np.int8(MEM), sp)
        sp = np.where(x_to_Y | a_to_Y | m_to_Y, np.int8(Y), sp)
        sp = np.where(y_to_X, np.int8(X), sp).astype(np.int8)
        if to_M.any():
            fresh = rng.integers(0, 6, size=sp.shape).astype(np.int8)
            ori = np.where(to_M, fresh, ori).astype(np.int8)

        yield sp, ori, fired


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------
def _analyse(sp, ori, membrane_id) -> dict:
    counts = {CODE[s]: int((sp == s).sum()) for s in range(5)}
    clusters = sorted(_components(sp == MEM), key=lambda c: len(c["cells"]), reverse=True)
    sizes = [len(c["cells"]) for c in clusters]
    biggest = clusters[0] if clusters else None
    gyration = _gyration(biggest["unwrapped"]) if biggest and not biggest["wraps"] else 0.0

    protocells, outside_cells = [], 0
    for comp in _components(sp != MEM):
        if comp["wraps"]:
            outside_cells += len(comp["cells"])
            continue
        content = {CODE[s]: 0 for s in range(5)}
        for q, r in comp["cells"]:
            content[CODE[sp[q, r]]] += 1
        protocells.append({"size": len(comp["cells"]), "contents": content})
    protocells.sort(key=lambda c: c["size"], reverse=True)

    # Shape of the membrane phase. Coordination is the mean number of membrane
    # neighbours per membrane particle: ~5 inside a droplet, ~2 along a
    # one-particle-thick filament. Alignment is the fraction of membrane
    # contacts that run along the +-o axis of both particles, which is how a
    # filament is held together: 1 for a well-formed membrane, 1/9 for random
    # orientations.
    is_m = sp == MEM
    n_m = int(is_m.sum())
    coordination = contacts = in_plane = 0
    o_self = ori.astype(int)
    for k in range(6):
        both = is_m & neighbour(is_m, k)
        coordination += int(both.sum())
        if k < 3:
            o_nb = neighbour(ori, k).astype(int)
            ok = ((o_self - k) % 3 == 0) & ((o_nb - k) % 3 == 0)
            contacts += int(both.sum())
            in_plane += int((both & ok).sum())

    return {
        "counts": counts,
        "membrane_species": membrane_id,
        "membrane_clusters": len(clusters),
        "largest_membrane_cluster": sizes[0] if sizes else 0,
        "membrane_cluster_sizes": sizes[:20],
        "largest_cluster_gyration": round(gyration, 4),
        # ~0.37 for a compact blob at any size, growing like 0.29 sqrt(size)
        # for a filament: the shape difference between M_i and M_a.
        "largest_cluster_elongation": (
            round(gyration / np.sqrt(sizes[0]), 4) if sizes and sizes[0] > 1 else 0.0),
        "protocells": len(protocells),
        "protocell_details": protocells[:20],
        "enclosed_cells": int(sum(c["size"] for c in protocells)),
        "enclosed_A": int(sum(c["contents"]["A"] for c in protocells)),
        "outside_cells": outside_cells,
        "mean_membrane_coordination": round(coordination / n_m, 4) if n_m else 0.0,
        "membrane_alignment": round(in_plane / contacts, 4) if contacts else 0.0,
    }


# ---------------------------------------------------------------------------
def _species(membrane_id: str) -> list[Species]:
    return [
        Species("A", structure="hydrophilic autocatalyst"),
        Species(membrane_id, structure=(
            "hydrophobic membrane particle, anisotropic: one of six lattice orientations"
            if membrane_id == "M_a" else "hydrophobic membrane particle, isotropic")),
        Species("X", structure="neutral food"),
        Species("Y", structure="neutral waste"),
        Species("W", structure="hydrophilic water"),
    ]


def _scheme(membrane_id: str) -> list[tuple[str, dict, dict]]:
    return [
        ("replicate", {"A": 1, "X": 1}, {"A": 2}),
        ("membrane", {"A": 1, "X": 1}, {"A": 1, membrane_id: 1}),
        ("decay_A", {"A": 1}, {"Y": 1}),
        ("decay_M", {membrane_id: 1}, {"Y": 1}),
        ("decay_X", {"X": 1}, {"Y": 1}),
        ("recycle", {"Y": 1}, {"X": 1}),
    ]


def _membrane_id(p) -> str:
    return "M_a" if p.membrane == "anisotropic" else "M_i"


def generate(p, rng):
    """The book's reaction set, as a complete network."""
    membrane_id = _membrane_id(p)
    return Network(
        species=_species(membrane_id),
        # No rate: no accessible source publishes a rate constant for this
        # two-dimensional scheme (see the catalog's `decisions`).
        reactions=[Reaction(lhs, rhs) for _, lhs, rhs in _scheme(membrane_id)],
        status="complete",
        extras={"reference_rates": JTB2000},
    )


def evolve(p, rng):
    """Run the lattice: a frame per sweep."""
    membrane_id = _membrane_id(p)
    scheme = _scheme(membrane_id)
    ids = [membrane_id if c == "M" else c for c in CODE]

    def census(sp):
        return {ids[k]: float(n) for k, n in enumerate(np.bincount(sp.ravel(), minlength=5).tolist()) if n}

    totals = {name: 0 for name, _, _ in scheme}
    initial_counts = None
    for t, (sp, ori, fired) in enumerate(_sweeps(p, rng)):
        if initial_counts is None:
            initial_counts = {ids[k]: int((sp == k).sum()) for k in range(5)}
        for name, n in fired.items():
            totals[name] += n
        yield Frame(t=float(t), state=census(sp), fired=[
            [sorted(_elements(lhs)), sorted(_elements(rhs)), fired[name]]
            for name, lhs, rhs in scheme if fired.get(name)])

    reactions = [
        Reaction(lhs, rhs, count=totals[name]) for name, lhs, rhs in scheme if totals[name]
    ]
    grid = ["".join(CODE[s] for s in row) for row in sp.tolist()]
    orientations = [
        "".join(str(o) if s == MEM else "." for s, o in zip(rs, ro))
        for rs, ro in zip(sp.tolist(), ori.tolist())
    ]
    energies = {
        "law": (
            "E = sum over neighbouring pairs; a hydrophilic/hydrophobic pair costs "
            "repulsion * F[k, o], a neutral/hydrophobic pair neutral_coupling * "
            "repulsion * F[k, o], every other pair 0"
        ),
        "repulsion": float(p.repulsion),
        "neutral_coupling": float(p.neutral_coupling),
        "temperature": float(p.temperature),
        "hydrophilic": ["A", "W"],
        "hydrophobic": [membrane_id],
        "neutral": ["X", "Y"],
        "anisotropy": 0.0 if membrane_id == "M_i" else float(p.anisotropy),
        "anisotropy_field": (
            "F[k, o] = 1 (isotropic M_i)" if membrane_id == "M_i" else
            "F[k, o] = 1 + anisotropy * cos(2 * 60deg * (k - o)); mean 1 over the "
            "six directions, so M_a redistributes the same total repulsion onto "
            "the two faces normal to its orientation o"
        ),
        "update": (
            "Kawasaki exchange with a neighbour and (for M_a) rotation, both "
            "accepted with min(1, exp(-dE/temperature))"
        ),
        "reference_rates": JTB2000,
    }
    return Network(
        species=_species(membrane_id),
        reactions=reactions,
        status="observed",
        initial_state={s: float(n) for s, n in initial_counts.items()},
        extras={
            "space": {
                "dimensions": 2,
                "lattice": "hexagonal",
                "coordinates": "axial (q, r) on a torus",
                "shape": [int(p.height), int(p.width)],
                "boundary": "periodic",
                "occupancy": "one particle per cell",
                "neighbour_directions": [list(d) for d in DIRS],
                "codes": {c: (membrane_id if c == "M" else c) for c in CODE},
                "final": grid,
                "orientations": orientations,
            },
            "energies": energies,
            "analysis": _analyse(sp, ori, membrane_id),
        },
    )


def _elements(side: dict) -> list[str]:
    return [s for s, n in side.items() for _ in range(n)]
