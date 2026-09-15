"""Ecolab (Standish): an evolving generalised Lotka-Volterra ecosystem. Catalog id: ecolab.

Population dynamics (Standish 2004, eq. 1):

    n_i' = r_i n_i - n_i sum_j beta_ij n_j

integrated as Ecolab does, one explicit timestep at a time with populations
kept integral by probabilistic rounding. Periodically each species produces
round(sp_sep r_i mu_i dt n_i) mutant species whose r, beta_ii, mu and
interaction row/column are varied from the parent by the operators of the
Ecolab technical report and `ecolab_model.cc`; extinct species are removed
(`condense`). After `cycles` generate/mutate/condense cycles the surviving
ecosystem is returned as its mass-action Lotka-Volterra network, the same
term-by-term mapping as `lotka_volterra.py`.

Internally beta is kept in Ecolab's source sign (n' = r n + n (b n), b_ii < 0);
extras and the docstrings use the paper sign beta = -b (beta_ii > 0).
"""

import numpy as np

from chemart.helpers.explicit import network

INT_MAX = 2**31 - 1
_MAX_ROUND = float(INT_MAX - 1)
_THRESHOLD = 10   # Ecolab's hardwired population for starting a lifetime


def _round(rng, x):
    """Ecolab's ROUND: floor(x) plus one with probability frac(x), clipped to [0, INT_MAX - 1]."""
    x = np.minimum(np.asarray(x, dtype=float), _MAX_ROUND)
    base = np.floor(x)
    out = base + (rng.random(x.shape) < x - base)
    return np.where(x <= 0, 0, out).astype(np.int64)


def _row_or_col(tmp, rng, lo, hi, gdist, gen_bias):
    """do_row_or_col: add or delete floor(1/|r|) - 1 links, then spread the nonzero values."""
    size = len(tmp)
    if size == 0:
        return tmp
    r = 2.0 * rng.random() - 1.0 + gen_bias
    if r > 0:
        r /= 1.0 + gen_bias
    elif r < 0:
        r /= 1.0 - gen_bias
    count = size if r == 0 else min(int(1.0 / abs(r)) - 1, size)
    if count > 0:
        pos = ((size - 1) * rng.random(count) + 0.5).astype(np.int64)
        if r > 0:
            # sequential semantics: only the first draw at an empty position sets a value
            values = (hi - lo) * rng.random(count) + lo
            upos, first = np.unique(pos, return_index=True)
            empty = tmp[upos] == 0.0
            tmp[upos[empty]] = values[first[empty]]
        else:
            tmp[pos] = 0.0
    nonzero = tmp != 0.0
    tmp[nonzero] += (hi - lo) * gdist * rng.standard_normal(int(nonzero.sum()))
    return tmp


def _bound_pairs(row, col):
    """Enforce b_ij + b_ji <= 0 (source sign) by removing the excess from the nonzero entries."""
    for j in np.flatnonzero(row + col > 0):
        if row[j] != 0.0 and col[j] != 0.0:
            row[j] = (row[j] - col[j]) / 2.0
            col[j] = -row[j]
        else:
            row[j] = col[j] = 0.0


class _Ecosystem:
    def __init__(self, p, rng):
        self.p, self.rng = p, rng
        n = p.nsp
        self.label = np.arange(n)
        self.next_label = n
        self.density = np.full(n, p.density, dtype=np.int64)
        self.r = rng.uniform(p.repro_min, p.repro_max, n)
        self.mu = np.full(n, p.mut_max)
        self.create = np.zeros(n, dtype=np.int64)
        B = np.zeros((n, n))
        for i in range(n):
            cols = rng.choice(n, size=p.conn, replace=False)
            cols = cols[cols != i]
            B[i, cols] = rng.uniform(p.odiag_min, p.odiag_max, len(cols))
        for i in range(n):
            row, col = B[i, :i], B[:i, i]
            _bound_pairs(row, col)
        np.fill_diagonal(B, -p.beta_diag)
        self.B = B
        self.tstep = 0
        self.last_mut = 0
        self.phylogeny = {
            int(s): {"parent": None, "born": 0, "r": float(ri), "died": None}
            for s, ri in zip(self.label, self.r)
        }
        self.history = {"timestep": [], "diversity": [], "speciations": [], "extinctions": []}
        self.lifetimes = []

    def generate(self, steps):
        n = self.density
        for _ in range(steps):
            n = _round(self.rng, n + n * (self.r + self.B @ n))
        self.density = n
        self.tstep += steps

    def mutate(self):
        p, rng = self.p, self.rng
        scale = p.sp_sep * self.r * self.mu * (self.tstep - self.last_mut)
        self.last_mut = self.tstep
        births = np.minimum(_round(rng, scale * self.density), self.density)
        parents = np.repeat(np.arange(len(births)), births)
        k, m0 = len(parents), len(self.density)
        if k == 0:
            return 0
        self.density = self.density - births

        gdist = rng.exponential(1.0, k) * self.mu[parents]
        r_new = self.r[parents] + (p.repro_max - p.repro_min) * gdist * rng.standard_normal(k)
        mu_new = np.minimum(self.mu[parents] * np.exp(gdist * rng.standard_normal(k)), p.mut_max)
        diag = np.diag(self.B)[parents] * np.exp(gdist * rng.standard_normal(k))   # stays negative
        diag = np.minimum(diag, -np.abs(r_new) / (0.1 * INT_MAX))

        B = np.zeros((m0 + k, m0 + k))
        B[:m0, :m0] = self.B
        B[range(m0, m0 + k), range(m0, m0 + k)] = diag
        for i, parent in enumerate(parents):
            q = m0 + i
            row, col = B[parent, :q].copy(), B[:q, parent].copy()
            row[parent] = col[parent] = B[parent, parent]
            _row_or_col(row, rng, p.odiag_min, p.odiag_max, gdist[i], p.gen_bias)
            _row_or_col(col, rng, p.odiag_min, p.odiag_max, gdist[i], p.gen_bias)
            _bound_pairs(row, col)
            B[q, :q], B[:q, q] = row, col
        self.B = B

        labels = np.arange(self.next_label, self.next_label + k)
        self.next_label += k
        for s, parent, ri in zip(labels, parents, r_new):
            self.phylogeny[int(s)] = {"parent": int(self.label[parent]), "born": self.tstep,
                                      "r": float(ri), "died": None}
        self.label = np.concatenate([self.label, labels])
        self.r = np.concatenate([self.r, r_new])
        self.mu = np.concatenate([self.mu, mu_new])
        self.create = np.concatenate([self.create, np.zeros(k, dtype=np.int64)])
        self.density = np.concatenate([self.density, np.ones(k, dtype=np.int64)])
        return k

    def record_lifetimes(self):
        started = (self.create == 0) & (self.density > _THRESHOLD)
        ended = (self.create > 0) & (self.density == 0)
        self.lifetimes += [int(t) for t in self.tstep - self.create[ended]]
        self.create[started] = self.tstep
        self.create[ended] = 0

    def condense(self):
        alive = self.density != 0
        for s in self.label[~alive]:
            self.phylogeny[int(s)]["died"] = self.tstep
        self.label, self.density = self.label[alive], self.density[alive]
        self.r, self.mu, self.create = self.r[alive], self.mu[alive], self.create[alive]
        self.B = self.B[np.ix_(alive, alive)]
        return int((~alive).sum())

    def cycle(self):
        self.generate(self.p.steps_per_cycle)
        births = self.mutate()
        self.record_lifetimes()
        deaths = self.condense()
        h = self.history
        h["timestep"].append(self.tstep)
        h["diversity"].append(len(self.label))
        h["speciations"].append(births)
        h["extinctions"].append(deaths)


def _check(p):
    if p.repro_min > p.repro_max:
        raise ValueError(f"repro_min ({p.repro_min}) must not exceed repro_max ({p.repro_max})")
    if p.odiag_min > p.odiag_max:
        raise ValueError(f"odiag_min ({p.odiag_min}) must not exceed odiag_max ({p.odiag_max})")
    if p.conn >= p.nsp:
        raise ValueError(f"conn ({p.conn}) must be smaller than nsp ({p.nsp})")
    if p.beta_diag <= 0:
        raise ValueError(f"beta_diag must be > 0 for boundedness, got {p.beta_diag}")


def _reactions(X, r, B):
    """Mass-action form of n_i' = n_i (r_i + sum_j b_ij n_j), b = -beta (as lotka_volterra.py)."""
    out = []
    for i, x in enumerate(X):
        if r[i] > 0:
            out.append((f"{x} -> 2 {x}", float(r[i])))
        elif r[i] < 0:
            out.append((f"{x} -> ", float(-r[i])))
    for i, x in enumerate(X):
        out.append((f"2 {x} -> {x}", float(-B[i, i])))

    def one_way(i, j, b):
        if b > 0:
            return [(f"{X[i]} + {X[j]} -> 2 {X[i]} + {X[j]}", float(b))]
        if b < 0:
            return [(f"{X[i]} + {X[j]} -> {X[j]}", float(-b))]
        return []

    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            bij, bji = B[i, j], B[j, i]
            if bij > 0 and bji == -bij:
                out.append((f"{X[j]} + {X[i]} -> 2 {X[i]}", float(bij)))
            elif bji > 0 and bij == -bji:
                out.append((f"{X[i]} + {X[j]} -> 2 {X[j]}", float(bji)))
            else:
                out += one_way(i, j, bij) + one_way(j, i, bji)
    return out


def generate(p, rng):
    _check(p)
    eco = _Ecosystem(p, rng)
    for _ in range(p.cycles):
        eco.cycle()

    X = [f"S{s}" for s in eco.label]
    beta = -eco.B
    m = len(X)
    offdiag = [[X[i], X[j], float(beta[i, j])]
               for i, j in zip(*np.nonzero(beta)) if i != j]
    if m:
        sign, _ = np.linalg.slogdet(beta)
        max_eig = float(np.max(np.linalg.eigvals(eco.B).real))
        diag = np.diag(beta)
        # May's connectivity as in ModelData::connectivity: sum b_ij^2 / (b_ii b_jj) over n(n-1)
        ratio = beta**2 / np.outer(diag, diag)
        np.fill_diagonal(ratio, 0.0)
        connectivity = float(ratio.sum() / (m * (m - 1))) if m > 1 else 0.0
    else:
        sign, max_eig, diag, connectivity = 0.0, None, np.array([]), 0.0
    phylogeny = {f"S{s}": {**v, "parent": None if v["parent"] is None else f"S{v['parent']}"}
                 for s, v in eco.phylogeny.items()}
    return network(
        _reactions(X, eco.r, eco.B),
        species=X,
        initial_state={x: float(n) for x, n in zip(X, eco.density)},
        extras={
            "ecosystem": {
                "equation": "n_i' = r_i n_i - n_i sum_j beta_ij n_j (one Ecolab timestep = one time unit)",
                "species": X,
                "r": [float(v) for v in eco.r],
                "beta_diag": [float(v) for v in diag],
                "beta_offdiag": offdiag,
                "mu": [float(v) for v in eco.mu],
            },
            "analysis": {
                "timestep": eco.tstep,
                "diversity": m,
                "speciations": int(sum(eco.history["speciations"])),
                "extinctions": int(sum(eco.history["extinctions"])),
                "history": eco.history,
                "lifetimes": eco.lifetimes,
                "phylogeny": phylogeny,
                "connectivity": connectivity,
                "max_eigenvalue_b": max_eig,
                "det_beta_positive": bool(sign > 0),
            },
        },
    )
