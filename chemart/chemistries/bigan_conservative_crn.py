"""Random conservative chemical reaction networks (Bigan, Steyaert & Douady 2013). Catalog id: bigan-conservative-crn.

Topology (paper sec. 2.1): candidate elementary reactions (association
Ai + Aj <-> Ak, which is also the dissociation, and transformation Ai <-> Ak)
are tried in random order and kept while a strictly positive mass vector m with
m^T S = 0 still exists, until the network has maximum size.

Kinetics (sec. 2.2): random formation energies G_i/RT ~ U(0, G_max); the
downhill direction gets k = k_avg 10^U(-s/2, s/2) and the uphill rate follows
from detailed balance, with standard concentration c0 = 1.
"""

import math
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from chemart.helpers.explicit import network, term


def _positive_mass(S: np.ndarray) -> np.ndarray | None:
    """A mass vector m >= 1 with m^T S = 0, or None if the network is not conservative."""
    n = S.shape[0]
    result = linprog(np.zeros(n), A_eq=S.T, b_eq=np.zeros(S.shape[1]), bounds=[(1, None)] * n, method="highs")
    return result.x if result.status == 0 else None


def _candidates(n: int):
    """(reactants, products) of every association/transformation, one direction each."""
    out = []
    for i in range(n):
        for k in range(i + 1, n):
            out.append(((i,), (k,)))
    for i in range(n):
        for j in range(i, n):
            for k in range(n):
                out.append(((i, j), (k,)))
    return out


def _column(n, lhs, rhs):
    c = np.zeros(n)
    for i in lhs:
        c[i] -= 1
    for k in rhs:
        c[k] += 1
    return c


def _integer_masses(m: np.ndarray, S: np.ndarray) -> list:
    scaled = [Fraction(v / m.min()).limit_denominator(1000) for v in m]
    lcm = math.lcm(*(f.denominator for f in scaled))
    ints = [int(f * lcm) for f in scaled]
    return ints if not np.any(np.array(ints) @ S) else [float(v) for v in m]


def generate(p, rng):
    n = p.N
    if not -1 <= p.nutrient < n:
        raise ValueError(f"nutrient must be a species index in 0..{n - 1}, or -1 for no nutrient flux; got {p.nutrient}")

    candidates = _candidates(n)
    columns, kept, mass, unique = [], [], np.ones(n), False
    for index in rng.permutation(len(candidates)):
        lhs, rhs = candidates[index]
        c = _column(n, lhs, rhs)
        if not c.any():
            continue
        if unique:
            if abs(c @ mass) < 1e-9:
                columns.append(c)
                kept.append((lhs, rhs))
        else:
            trial = np.column_stack(columns + [c])
            m = _positive_mass(trial)
            if m is None:
                continue
            columns.append(c)
            kept.append((lhs, rhs))
            mass = m
            unique = n - np.linalg.matrix_rank(trial) == 1
        if p.max_reactions and len(kept) >= p.max_reactions:
            break

    A = [f"A{i}" for i in range(n)]
    G = rng.uniform(0.0, p.G_max, n)

    def side(indices):
        return term(2, A[indices[0]]) if len(indices) == 2 and indices[0] == indices[1] else " + ".join(A[i] for i in indices)

    def rate(k):
        if p.kinetics == "mass-action":
            return float(k)
        return {"law": "saturating", "k": float(k), "K": float(p.K_avg * 10 ** rng.uniform(-p.p / 2, p.p / 2))}

    reactions, pairs = [], []
    for lhs, rhs in kept:
        dG = sum(G[k] for k in rhs) - sum(G[i] for i in lhs)      # in units of RT
        if dG > 0:
            lhs, rhs, dG = rhs, lhs, -dG                             # forward = downhill
        k_avg = p.k_avg_bi if len(lhs) == 2 else p.k_avg_mono
        k_forward = k_avg * 10 ** rng.uniform(-p.s / 2, p.s / 2)
        # K = k_forward / k_reverse = exp(|dG|/RT) c0^(n_products - n_reactants), c0 = 1
        k_reverse = k_forward * math.exp(-abs(dG))
        pairs.append([len(reactions), len(reactions) + 1])
        reactions += [(f"{side(lhs)} -> {side(rhs)}", rate(k_forward)), (f"{side(rhs)} -> {side(lhs)}", rate(k_reverse))]

    net = network(
        reactions,
        species=A,
        initial_state={a: p.initial_concentration for a in A},
        inflow={A[p.nutrient]: p.nutrient_flux} if p.nutrient >= 0 else None,
    )
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    masses = _integer_masses(mass, S) if kept else [1] * n
    net.extras = {
        "energies": {"formation_free_energy_RT": {a: float(g) for a, g in zip(A, G)}},
        "conservation": [{"name": "mass", "vector": dict(zip(A, masses)), "unique": bool(unique)}],
        "reaction_pairs": pairs,
    }
    return net
