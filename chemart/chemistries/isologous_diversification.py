"""Isologous diversification (Kaneko & Yomo 1997). Catalog id: isologous-diversification.

Each cell holds k + 1 chemicals: a source (nutrient) x(0) and x(1)..x(k).
The intracellular network is a set of catalysed paths m -> l aided by j
(eq. 1), the source feeds every chemical (S(l) = 1), and the chemicals with
a path to the "division factor" decay at gamma (P(l) = 1); the integral of
that decay is the division condition (eq. 8).

The reactions of the returned Network are exactly those terms::

    X0 + Xl -> 2 Xl                   e0 x(0) x(l)                  (source)
    Xm + Xj -> Xl + Xj                e1 x(j) x(m) / (1 + x(m)/x_M) (eq. 1)
    Xl -> DF                          gamma x(l)                    (eq. 1)

Cells are compartments (``extras["compartments"]``: division at DF > R with
an almost equal split, death at sum_l x(l) < S), and the coupling through the
well-stirred medium -- active transport p (sum_l x(l)) X(m) (eq. 2), diffusion
D (X(m) - x(m)) (eq. 3) and the medium's own flow (eqs. 6-7) -- is in
``extras["interaction_law"]``: it mixes the cell and the medium with the
volume factor 1/V, so it is not a reaction of one cell's network.

``Model`` integrates the whole cell society (eqs. 1-9) and the four published
stages measured from that run are reported in ``extras["analysis"]``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from chemart.helpers.explicit import network, term

#: Kaneko & Yomo (1997), caption of Fig. 6: the parameters of the main
#: simulation (network of Fig. 5a, k = 8, 3 autocatalytic paths per chemical).
PAPER = {
    "transport": 10.0,           # p
    "e0": 1.0, "e1": 1.0,        # reaction coefficients
    "nutrient": 40.0,            # X0, the external source concentration
    "diffusion": 0.02,           # D
    "gamma": 0.2,                # decay to the division factor
    "x_M": 10.0,                 # Michaelis-Menten constant
    "division_threshold": 2000.0,  # R
    "death_threshold": 0.05,     # S
    "flow": 0.005,               # f
    "washout": 0.005,            # D_out
    "medium_volume": 1000.0,     # V
    "split_noise": 0.001,        # epsilon, uniform over [-1e-3, 1e-3]
}

#: Caption of Fig. 18: the second published parameter set (network of Fig. 5b,
#: 2 paths per chemical), used for the simulations with cell death.
PAPER_FIG18 = {**PAPER, "nutrient": 10.0, "division_threshold": 100.0,
               "death_threshold": 0.01}

ENZYMES = ("michaelis-menten", "quadratic")

#: Screening a candidate network the way the paper does: it integrates a single
#: cell in its medium with neither division nor death (its Fig. 4) and keeps the
#: network when the concentrations still oscillate after the transient --
#: "Only for medium number of reaction paths, non-trivial oscillations of
#: chemicals appear as in Fig.1. We use such network for our simulation"
#: (Kaneko & Yomo 1998, section 2.4).
SCREEN_SETTLE = 60.0      # transient discarded before measuring
SCREEN_WINDOW = 20.0      # window the oscillation is measured over
SCREEN_SAMPLES = 200
SCREEN_AMPLITUDE = 0.05   # (max - min)/mean of the most strongly varying chemical
SCREEN_DIVISION_SHARE = 1e-3   # the activated chemicals must feed the division factor

#: The four stages of the isologous diversification scenario (1997, section 5);
#: stage 5 (successive differentiation) is stage 4 repeated.
STAGES = {
    0: "single cell: the attractor of the intracellular dynamics",
    1: "synchronous oscillation of identical cells (cell number 1, 2, 4, 8, ...)",
    2: "clustering by the phase of the oscillations; the temporal averages are still equal",
    3: "fixed differentiation: the temporal averages differ, giving distinct cell types",
    4: "determination: the differentiated type is inherited by the daughter cells",
}


# ---------------------------------------------------------------------------
# the reaction network
# ---------------------------------------------------------------------------
def draw_paths(k: int, connections: int, autocatalytic: bool, rng) -> list[list[int]]:
    """`connections` outgoing paths per chemical, as (m, l, j) with m -> l aided by j.

    With `autocatalytic` the catalyst of a path is its own product (j = l),
    which is the choice of the paper (appendix 2); otherwise it is drawn at
    random among all k chemicals.
    """
    others = np.arange(1, k + 1)
    paths = []
    for m in range(1, k + 1):
        targets = rng.choice(others[others != m], size=connections, replace=False)
        for l in sorted(int(t) for t in targets):
            j = l if autocatalytic else int(rng.integers(1, k + 1))
            paths.append([m, l, j])
    return paths


def check_paths(value, k: int) -> list[list[int]]:
    out = []
    for path in value:
        if (not isinstance(path, (list, tuple)) or len(path) != 3
                or any(isinstance(v, bool) or not isinstance(v, int) for v in path)):
            raise ValueError(f"paths: every entry must be a triple of integers "
                             f"[m, l, j] (substrate, product, catalyst), got {path!r}")
        m, l, j = (int(v) for v in path)
        if not all(1 <= v <= k for v in (m, l, j)):
            raise ValueError(f"paths: {path!r} uses a chemical outside 1..{k} "
                             f"(the source chemical 0 is fed in separately)")
        if m == l:
            raise ValueError(f"paths: {path!r} is a path from a chemical to itself")
        out.append([m, l, j])
    if not out:
        raise ValueError("paths must list at least one reaction path [m, l, j]")
    return out


def reactions_of(k, paths, source, division, enzyme, e0, e1, gamma, x_M):
    """The terms of eq. (1), one reaction each."""
    def name(i):
        return f"X{i}"

    out = []
    for l in source:
        out.append((f"X0 + {name(l)} -> 2 {name(l)}", {"law": "mass-action", "k": float(e0)}))
    order = 2 if enzyme == "quadratic" else 1
    for m, l, j in paths:
        rate = {"law": "mass-action", "k": float(e1)}
        if enzyme == "michaelis-menten":
            # e1 x(j) x(m) / (1 + x(m)/x_M): the Michaelis-Menten factor is not
            # one of Chemart's rate laws, so it is recorded as extra keys.
            rate.update(saturation="michaelis-menten",
                        saturated_species=name(m), x_M=float(x_M))
        left = f"{name(m)} + {term(order, name(j))}"
        right = (term(order + 1, name(l)) if j == l
                 else f"{name(l)} + {term(order, name(j))}")
        out.append((f"{left} -> {right}", rate))
    for l in division:
        out.append((f"{name(l)} -> DF", {"law": "mass-action", "k": float(gamma)}))
    return out


# ---------------------------------------------------------------------------
# the cell society (eqs. 1-9)
# ---------------------------------------------------------------------------
@dataclass
class Model:
    """Eqs. (1)-(9) for N cells in one well-stirred medium.

    The state of a cell is ``[x(0..k), DF, A(0..k)]``: the concentrations, the
    division factor DF = integral of sum_l gamma P(l) x(l) since the cell's
    birth (eq. 8), and A, the integral of the concentrations since that birth,
    which gives the temporal averages the paper classifies cells by.
    """

    k: int
    paths: list[list[int]]
    source: list[int]
    division: list[int]
    enzyme: str = "michaelis-menten"
    e0: float = PAPER["e0"]
    e1: float = PAPER["e1"]
    gamma: float = PAPER["gamma"]
    x_M: float = PAPER["x_M"]
    transport: float = PAPER["transport"]
    diffusion: float = PAPER["diffusion"]
    medium_volume: float = PAPER["medium_volume"]
    flow: float = PAPER["flow"]
    washout: float = PAPER["washout"]
    nutrient: float = PAPER["nutrient"]
    division_threshold: float = PAPER["division_threshold"]
    death_threshold: float = PAPER["death_threshold"]
    split_noise: float = PAPER["split_noise"]
    rtol: float = 1e-5
    atol: float = 1e-7

    def __post_init__(self):
        arr = np.array(self.paths, dtype=int).reshape(-1, 3)
        self._m, self._l, self._j = arr[:, 0], arr[:, 1], arr[:, 2]
        self._source = np.array(self.source, dtype=int)
        self._division = np.array(self.division, dtype=int)
        self.width = 2 * self.k + 3          # x(0..k), DF, A(0..k)

    # -- layout --------------------------------------------------------
    def unpack(self, y, n_cells):
        cells = y[: n_cells * self.width].reshape(n_cells, self.width)
        return cells, y[n_cells * self.width:]

    def pack(self, cells, medium):
        return np.concatenate([cells.ravel(), medium])

    def initial(self, rng):
        cells = np.zeros((1, self.width))
        cells[0, 1: self.k + 1] = rng.random(self.k)
        medium = np.zeros(self.k + 1)
        medium[0] = self.nutrient
        return cells, medium

    # -- the right-hand side -------------------------------------------
    def derivative(self, cells, medium):
        """d/dt of (cells, medium): eqs. (1)-(3) with (4)-(7)."""
        k = self.k
        x = np.maximum(cells[:, : k + 1], 0.0)
        X = np.maximum(medium, 0.0)
        dx = np.zeros_like(x)

        # eq. (1), reaction paths m -> l aided by j
        xm, xj = x[:, self._m], x[:, self._j]
        if self.enzyme == "quadratic":
            flux = self.e1 * xm * xj ** 2
        else:
            flux = self.e1 * xj * xm / (1.0 + xm / self.x_M)
        np.add.at(dx.T, self._l, flux.T)
        np.subtract.at(dx.T, self._m, flux.T)

        # eq. (1), paths from the source chemical, and eq. (4)
        if len(self._source):
            fed = self.e0 * x[:, :1] * x[:, self._source]
            np.add.at(dx.T, self._source, fed.T)
            dx[:, 0] -= fed.sum(1)

        # eq. (1), paths to the division factor
        d_factor = np.zeros(len(x))
        if len(self._division):
            decay = self.gamma * x[:, self._division]
            np.subtract.at(dx.T, self._division, decay.T)
            d_factor = decay.sum(1)

        # eqs. (2) and (3): active transport and diffusion through the membrane
        activity = x[:, 1:].sum(1)
        transport = self.transport * activity[:, None] * X[None, :]
        diffusion = self.diffusion * (X[None, :] - x)
        exchange = transport + diffusion
        dx += exchange

        # eqs. (6) and (7): the medium
        dX = -exchange.sum(0) / self.medium_volume
        dX[0] += self.flow * (self.nutrient - X[0])
        dX[1:] -= self.washout * X[1:]

        d_cells = np.empty_like(cells)
        d_cells[:, : k + 1] = dx
        d_cells[:, k + 1] = d_factor
        d_cells[:, k + 2:] = x
        return d_cells, dX

    def rhs(self, n_cells):
        def f(t, y):
            cells, medium = self.unpack(y, n_cells)
            d_cells, d_medium = self.derivative(cells, medium)
            return self.pack(d_cells, d_medium)
        return f

    # -- the run -------------------------------------------------------
    def run(self, max_cells, t_max, rng):
        """Integrate until `max_cells` cells exist or t reaches `t_max`.

        Cells that reach the threshold together divide together: while the cells
        are synchronised (stage 1) that is the rule rather than the exception.
        """
        cells, medium = self.initial(rng)
        start = cells[0, : self.k + 1].copy()
        k, t = self.k, 0.0
        info = [{"id": 0, "parent": -1, "born": 0.0, "mother": None, "mother_life": 0.0}]
        next_id, divisions, deaths = 1, 0, 0
        history: list[tuple[np.ndarray, np.ndarray]] = []   # (mother, daughter) averages

        while t < t_max and len(cells):
            n = len(cells)
            events = [self._division_event(n)] if n < max_cells else []
            if self.death_threshold > 0:
                events.append(self._death_event(n))
            if not events:
                break
            sol = solve_ivp(self.rhs(n), (t, t_max), self.pack(cells, medium),
                            method="LSODA", events=events,
                            rtol=self.rtol, atol=self.atol)
            if not sol.success:
                break
            t = float(sol.t[-1])
            cells, medium = self.unpack(sol.y[:, -1], n)
            cells, medium = cells.copy(), medium.copy()

            ready = (np.flatnonzero(cells[:, k + 1] >= self.division_threshold * (1 - 1e-9))
                     if n < max_cells else np.empty(0, dtype=int))
            starved = np.flatnonzero(
                np.maximum(cells[:, 1: k + 1], 0.0).sum(1)
                <= self.death_threshold * (1 + 1e-9)
            ) if self.death_threshold > 0 else np.empty(0, dtype=int)
            if not len(ready) and not len(starved):
                break

            for i in ready:                                        # division (eq. 8)
                if len(cells) >= max_cells:
                    break
                life = t - info[i]["born"]
                average = self._average(cells[i], life)
                if info[i]["mother"] is not None:
                    # the paper's return map (Fig. 12b): a full life of the
                    # mother against a full life of the daughter
                    history.append((info[i]["mother"], average))
                noise = float(rng.uniform(-self.split_noise, self.split_noise))
                other = cells[i].copy()
                cells[i, : k + 1] *= 0.5 + noise
                other[: k + 1] *= 0.5 - noise
                cells[i, k + 1:] = 0.0
                other[k + 1:] = 0.0
                cells = np.vstack([cells, other])
                child = {"parent": info[i]["id"], "born": t,
                         "mother": average, "mother_life": life}
                info[i] = {"id": next_id, **child}
                info.append({"id": next_id + 1, **child})
                next_id += 2
                divisions += 1

            for i in sorted(starved, reverse=True):                # death (eq. 9)
                medium = medium + np.maximum(cells[i, : k + 1], 0.0) / self.medium_volume
                cells = np.delete(cells, i, axis=0)
                info.pop(i)
                deaths += 1

        age = [t - c["born"] for c in info]
        averages = np.array(
            [self._average(row, a) for row, a in zip(cells, age)]
        ).reshape(-1, k + 1)
        # A cell that has just been born has no meaningful average yet: it counts
        # only once it has lived at least half as long as its mother did.
        mature = [a >= 0.5 * c["mother_life"] for a, c in zip(age, info)]
        for row, c, ripe in zip(averages, info, mature):
            if c["mother"] is not None and ripe:
                history.append((c["mother"], row))
        return {
            "t": t, "cells": cells, "medium": medium, "averages": averages,
            "initial": start, "lineage": info, "age": age, "mature": mature,
            "return_map": history, "divisions": divisions, "deaths": deaths,
        }

    def _average(self, row, age):
        if age <= 1e-9:
            return np.maximum(row[: self.k + 1], 0.0)
        return row[self.k + 2:] / age

    def _division_event(self, n_cells):
        k = self.k

        def division(t, y):
            cells, _ = self.unpack(y, n_cells)
            return float(self.division_threshold - cells[:, k + 1].max())
        division.terminal, division.direction = True, -1
        return division

    def _death_event(self, n_cells):
        k = self.k

        def death(t, y):
            cells, _ = self.unpack(y, n_cells)
            return float(np.maximum(cells[:, 1: k + 1], 0.0).sum(1).min()
                         - self.death_threshold)
        death.terminal, death.direction = True, -1
        return death


# ---------------------------------------------------------------------------
# cell types: the clustering of the temporal averages (1998, eqs. 5-6)
# ---------------------------------------------------------------------------
def composition(vectors: np.ndarray) -> np.ndarray:
    """Normalised composition of the k chemicals, dropping the source x(0).

    A cell's type is its chemical composition, not the amount it holds: while
    the population grows every cell's concentrations fall together (the medium
    feeds N cells from one flow), so the composition is what can be compared
    across divisions. Kaneko & Yomo (1998) build the normalisation into the
    model itself, as sum_l x(l) = 1 (their eqs. 1-2).
    """
    x = np.maximum(np.atleast_2d(vectors)[:, 1:], 0.0)
    total = x.sum(1, keepdims=True)
    return x / np.where(total > 0.0, total, 1.0)


def classify(compositions: np.ndarray, tolerance: float) -> list[int]:
    """Group cells whose averaged compositions are within `tolerance`.

    The paper's criterion (Kaneko & Yomo 1998, eqs. 5-6) is the Euclidean
    distance between the time-averaged concentration vectors: it is much
    smaller within a type than between types. Chemart makes that a
    single-linkage clustering of the compositions at `tolerance`.
    """
    n = len(compositions)
    if n == 0:
        return []
    label = list(range(n))

    def find(i):
        while label[i] != i:
            label[i] = label[label[i]]
            i = label[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if float(np.linalg.norm(compositions[i] - compositions[j])) < tolerance:
                label[find(i)] = find(j)
    roots = {}
    return [roots.setdefault(find(i), len(roots)) for i in range(n)]


def spread(compositions: np.ndarray) -> float:
    """Largest pairwise distance between compositions."""
    if len(compositions) < 2:
        return 0.0
    d = np.linalg.norm(compositions[:, None, :] - compositions[None, :, :], axis=-1)
    return float(d.max())


def analyse(run: dict, k: int, tolerance: float, network: dict | None = None) -> dict:
    """The stage reached, the cell types and the recursivity of the run."""
    cells, averages = run["cells"], run["averages"]
    # Only cells that have lived long enough to have an average are classified:
    # right after a division a cell holds half of its mother's concentrations.
    keep = [i for i, ripe in enumerate(run["mature"]) if ripe] or list(range(len(cells)))
    grown = averages[keep] if len(averages) else averages
    snapshots = composition(cells[keep]) if len(cells) else np.empty((0, k))
    mean_composition = composition(grown) if len(grown) else np.empty((0, k))
    average_spread = spread(mean_composition)
    types = classify(mean_composition, tolerance)
    n_types = len(set(types))

    pairs = run["return_map"]
    recursivity = inherited = None
    if pairs:
        distance = [float(np.linalg.norm(composition(a) - composition(b)))
                    for a, b in pairs]
        recursivity = float(np.mean(distance))
        inherited = float(np.mean([d < tolerance for d in distance]))

    if len(cells) <= 1:
        stage = 0
    elif spread(snapshots) < tolerance and average_spread < tolerance:
        stage = 1
    elif average_spread < tolerance:
        stage = 2
    elif inherited is not None and inherited >= 0.5:
        stage = 4
    else:
        stage = 3

    signatures = []
    for t in range(n_types):
        members = [i for i, label in enumerate(types) if label == t]
        signatures.append({
            "type": t,
            "cells": len(members),
            "activity": float(np.maximum(grown[members][:, 1:], 0.0).sum(1).mean()),
            "composition": [round(float(v), 6) for v in mean_composition[members].mean(0)],
        })
    signatures.sort(key=lambda s: -s["activity"])
    for rank, s in enumerate(signatures):
        s["type"] = rank

    out = {
        "t_end": round(float(run["t"]), 6),
        "cells": int(len(cells)),
        "cells_classified": int(len(grown)),
        "divisions": int(run["divisions"]),
        "deaths": int(run["deaths"]),
        "n_types": int(n_types),
        "types": signatures,
        "stage": int(stage),
        "stage_name": STAGES[stage],
        "snapshot_spread": round(spread(snapshots), 6),
        "average_spread": round(average_spread, 6),
        "recursivity": None if recursivity is None else round(recursivity, 6),
        "inherited_fraction": None if inherited is None else round(inherited, 6),
        "divisions_measured": len(pairs),
        "type_tolerance": float(tolerance),
    }
    out.update(network or {})
    out["definitions"] = {
        "composition": "a cell's averaged concentrations of X1..Xk, divided by their "
                       "sum: while the population grows every cell's concentrations "
                       "fall together, so the composition is what is comparable",
        "stage": "the stage of the isologous diversification scenario "
                 "(Kaneko & Yomo 1997, section 5) that the run reached; "
                 + "; ".join(f"{i}: {text}" for i, text in STAGES.items()),
        "snapshot_spread": "largest pairwise distance between the instantaneous "
                           "compositions of the cells",
        "average_spread": "the same for the compositions averaged since each cell's "
                          "last division (the paper's Fig. 8)",
        "n_types": "number of groups of the averaged compositions, single-linkage at "
                   "type_tolerance (1998, eqs. 5-6)",
        "activity": "sum_l x(l) of a type's cells, the paper's measure of how strong "
                    "a cell is; strong cells divide faster",
        "recursivity": "mean distance between a mother's averaged composition over her "
                       "whole life and her daughter's; 0 is the diagonal of the "
                       "paper's return map (Fig. 12b)",
        "inherited_fraction": "fraction of those mother-daughter pairs that stay "
                              "within type_tolerance, i.e. that keep the type",
        "oscillatory": "whether the chosen network's single-cell dynamics oscillates, "
                       "which is what the paper selects its networks for",
    }
    return out


# ---------------------------------------------------------------------------
# the catalogued generator
# ---------------------------------------------------------------------------
def single_cell(model: Model, rng, settle=SCREEN_SETTLE, window=SCREEN_WINDOW,
                samples=SCREEN_SAMPLES):
    """One cell in its medium, with neither division nor death (the paper's Fig. 4).

    Returns the concentrations x(0..k) sampled over `window` after the
    transient `settle`, or None if the integration fails.
    """
    cells, medium = model.initial(rng)
    sol = solve_ivp(model.rhs(1), (0.0, settle + window), model.pack(cells, medium),
                    method="LSODA", t_eval=np.linspace(settle, settle + window, samples),
                    rtol=model.rtol, atol=model.atol)
    if not sol.success or sol.y.shape[1] < 2:
        return None
    return np.maximum(sol.y[: model.k + 1, :].T, 0.0)


def oscillation(trace) -> tuple[float, float]:
    """(relative amplitude of the most strongly varying chemical, division share)."""
    if trace is None:
        return 0.0, 0.0
    x = trace[:, 1:]
    mean = x.mean(0)
    amplitude = np.where(mean > 0, (x.max(0) - x.min(0)) / np.where(mean > 0, mean, 1.0), 0.0)
    return float(amplitude.max(initial=0.0)), float(mean.sum())


def screen(params: dict, rng, attempts: int):
    """Draw networks until one oscillates as a single cell, as the paper selects them.

    The paper keeps only networks whose intracellular dynamics oscillates (1997,
    section 3 and appendix 2; 1998, section 2.4), because the whole scenario
    rests on coupled oscillators. Chemart additionally requires that the
    chemicals that stay active feed the division factor, since a cell that
    never reaches the division condition cannot show the scenario either.
    """
    k, best = params["n_chemicals"], None
    for attempt in range(1, attempts + 1):
        paths = draw_paths(k, params["connections"], params["autocatalytic"], rng)
        division = draw_division_paths(k, params["n_division_paths"], rng)
        if attempts == 1:
            return paths, division, {
                "network_attempts": 1, "screened": False, "oscillatory": None,
                "single_cell_amplitude": None, "division_factor_share": None}
        model = model_of(params, paths, division)
        trace = single_cell(model, rng)
        amplitude, total = oscillation(trace)
        share = 0.0
        if trace is not None and total > 0 and division:
            share = float(trace[:, division].mean(0).sum() / total)
        feeds = share >= SCREEN_DIVISION_SHARE
        report = {"network_attempts": attempt, "screened": True,
                  "accepted": bool(amplitude >= SCREEN_AMPLITUDE and feeds),
                  "oscillatory": bool(amplitude >= SCREEN_AMPLITUDE),
                  "single_cell_amplitude": round(amplitude, 6),
                  "division_factor_share": round(share, 6)}
        if report["accepted"]:
            return paths, division, report
        # keep the best candidate so far: one that at least feeds the division
        # factor, and among those the one that varies most
        key = (feeds, amplitude)
        if best is None or key > best[0]:
            best = (key, paths, division, report)
    _, paths, division, report = best
    return paths, division, {**report, "network_attempts": attempts, "accepted": False}


def draw_division_paths(k: int, count: int, rng) -> list[int]:
    """The chemicals with a path to the division factor, P(l) = 1."""
    if not count:
        return []
    return sorted(int(v) for v in rng.choice(np.arange(1, k + 1), size=count, replace=False))


def model_of(params: dict, paths: list[list[int]], division: list[int]) -> Model:
    """Build the ODE model from resolved catalog parameters and a drawn network."""
    k = params["n_chemicals"]
    return Model(
        k=k, paths=paths, source=list(range(1, k + 1)), division=division,
        enzyme=params["enzyme"],
        **{name: float(params[name]) for name in (
            "e0", "e1", "gamma", "x_M", "transport", "diffusion", "medium_volume",
            "flow", "washout", "nutrient", "division_threshold", "death_threshold",
            "split_noise")},
    )


def model_from_network(net) -> Model:
    """Rebuild the exact ODE model of a generated network."""
    return model_of(net.params, net.extras["reaction_network"]["paths"],
                    net.extras["reaction_network"]["division_factor_paths"])


def generate(p, rng):
    k = p.n_chemicals
    if not p.paths and p.connections > k - 1:
        raise ValueError(f"connections = {p.connections} needs at least "
                         f"{p.connections + 1} chemicals, but n_chemicals = {k}")
    if p.n_division_paths > k:
        raise ValueError(f"n_division_paths = {p.n_division_paths} is more than "
                         f"n_chemicals = {k}")
    params = vars(p)
    if p.paths:
        paths = check_paths(p.paths, k)
        division = draw_division_paths(k, p.n_division_paths, rng)
        report = {"network_attempts": 0, "screened": False, "oscillatory": None,
                  "single_cell_amplitude": None, "division_factor_share": None}
    else:
        paths, division, report = screen(params, rng, p.network_attempts)
    source = list(range(1, k + 1))                       # S(l) = 1 for every chemical

    model = model_of(params, paths, division)
    run = model.run(p.max_cells, p.t_max, rng)
    analysis = analyse(run, k, p.type_tolerance, report)

    reactions = reactions_of(k, paths, source, division, p.enzyme,
                             p.e0, p.e1, p.gamma, p.x_M)
    species = ["X0"] + [f"X{i}" for i in range(1, k + 1)] + ["DF"]
    initial = {f"X{i}": round(float(v), 9) for i, v in enumerate(run["initial"])}

    return network(
        reactions,
        species=species,
        initial_state=initial,
        extras={
            "reaction_network": {
                "paths": paths,
                "source_paths": source,
                "division_factor_paths": division,
                "notation": "paths are [m, l, j]: Con(m, l, j) = 1, the reaction from "
                            "chemical m to l catalysed by chemical j (eq. 1); "
                            "source_paths are the l with S(l) = 1, "
                            "division_factor_paths the l with P(l) = 1",
            },
            "compartments": {
                "cell": {
                    "chemicals": [f"X{i}" for i in range(k + 1)],
                    "volume": "constant, except for the short span of a division in "
                              "which it doubles, so the concentrations are halved",
                    "division": f"when DF > R = {float(p.division_threshold)}, i.e. when "
                                "the integral of sum_l gamma P(l) x(l) since the cell's "
                                "birth passes R (eq. 8); the cell splits into two, with "
                                f"(1/2 +- eps) of every concentration, eps uniform over "
                                f"[-{float(p.split_noise)}, {float(p.split_noise)}]",
                    "death": f"when sum_l x(l) < S = {float(p.death_threshold)} (eq. 9); "
                             "the cell's chemicals are released into the medium, "
                             "each divided by the medium volume V",
                    "cells_simulated": int(analysis["cells"]),
                }
            },
            "interaction_law": {
                "name": "cell-cell interaction through the well-stirred medium "
                        "(Kaneko & Yomo 1997, eqs. 2-3 and 6-7)",
                "state": "x(m) is the concentration of chemical m in a cell, X(m) its "
                         "concentration in the medium, V the volume of the medium in "
                         "units of one cell, N the number of cells",
                "active_transport": "Transp_i(m) = p (sum_{l=1..k} x_i(l)) X(m), "
                                    f"p = {float(p.transport)} (eq. 2): a cell with more "
                                    "chemicals takes up more",
                "diffusion": f"Diff_i(m) = D (X(m) - x_i(m)), D = {float(p.diffusion)} (eq. 3)",
                "cell": "dx_i(m)/dt = (network terms) + Transp_i(m) + Diff_i(m) (eqs. 4-5)",
                "medium": "dX(m)/dt = -(1/V) sum_i {Transp_i(m) + Diff_i(m)} - D_out X(m) "
                          "for m > 0 (eq. 6), and dX(0)/dt = f (X0 - X(0)) "
                          "- (1/V) sum_i {Transp_i(0) + Diff_i(0)} for the source (eq. 7)",
                "constants": {"p": float(p.transport), "D": float(p.diffusion),
                              "V": float(p.medium_volume), "f": float(p.flow),
                              "D_out": float(p.washout), "X0": float(p.nutrient)},
                "medium_state": [round(float(v), 9) for v in run["medium"]],
                "reason": "the exchange couples a cell to the medium with the volume "
                          "factor 1/V, so it is not a reaction of one cell's network",
            },
            "analysis": analysis,
        },
    )
