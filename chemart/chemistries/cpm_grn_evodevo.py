"""Cellular Potts + GRN evo-devo (book 18.6.1). Catalog id: cpm-grn-evodevo.

Hogeweg's evo-devo model: a Boolean gene-regulation network, read through
lock-and-key surface receptors, sets each cell's adhesion, and a Cellular
Potts Model (Glazier & Graner 1993) turns those adhesions into movement,
growth, division and death. Evolution selects the network for cell
differentiation alone; morphogenesis appears as a side-effect.

**The Cellular Potts Model.** A lattice site carries the identity sigma of
the cell occupying it (0 = medium), so one cell is many sites. The energy is

    H = sum over neighbouring site pairs with different sigma of J(type, type)
        + lambda * sum over cells of (v - V)^2

(Savill & Hogeweg 1997 eq. 1; Hogeweg 2000 eq. 1), with v the cell's volume
in sites and V its target volume. A copy attempt takes a random interior site
i and a random Moore neighbour j and copies sigma(j) into i with

    P = 1                                if dH < -0.1
    P = exp(-(dH + 0.1) / T)             otherwise

(Savill & Hogeweg 1997 eq. 2). One time step is one copy attempt per site.

**Adhesion.** Of the GRN's nodes, ``adhesion_nodes`` set the surface
receptors: the first half are "locks", the second half "keys", matched
complementarily (Hogeweg 2000 eqs. 2a-2b):

    J_ij = (sum_k 2^k [key k of i != lock k of j]
            + sum_k 2^k [key k of j != lock k of i]) / 2
    J_im = sum_l 2^l (node l of i),  l a window straddling the locks and keys

**Development.** One zygote undergoes ``divisions`` pre-scheduled cleavages,
each followed by ``steps_per_stage`` time steps (``final_steps`` after the
last). Both daughters inherit the mother's GRN state; after the first and the
second cleavage one maternal node is flipped in one daughter for one step.
A cell whose volume reaches 0 dies; when stretching takes its volume above
V + ``growth_threshold`` its target grows by 1, and at twice the reference
target it divides across the middle of its longest axis.

**Evolution.** A population of networks is scored by the amount of cell
differentiation - the summed Hamming distance between the expression patterns
of the distinct settled cell types, minimised over the last ``fitness_window``
steps - and bred by tournament selection with point mutations. Growth and
division are disabled during evolution, as in the paper, and enabled for the
final development of the best network, which is what the network reports.

The network is the observed events of that development, each with its count:

    Ta -> Tb      a cell of type a differentiates into type b
    Ta -> 2 Ta    a cell divides (both daughters inherit the mother's state)
    Ta ->         a cell is squeezed to zero volume and dies

Species are cell types, ``T<k>`` in order of first appearance, with the
expression pattern as their structure. ``mode="cell-sorting"`` runs the
Glazier & Graner differential-adhesion experiment instead: a mixed aggregate
of two fixed cell types sorts out, and the network carries the CPM itself in
``extras["interaction_law"]``.
"""

from __future__ import annotations

import math

import numpy as np

from chemart.network import Network, Reaction, Species

MEDIUM = 0

#: Savill & Hogeweg (1997) eq. 2: a copy whose dH is below -0.1 is certain.
#: The offset keeps cells still when nothing pushes them.
DELTA_H_OFFSET = 0.1

#: Moore neighbourhood; the second half alone visits every unordered pair once.
OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
HALF_OFFSETS = ((0, 1), (1, -1), (1, 0), (1, 1))

#: The 16 two-input Boolean functions, as 4-bit truth tables.
FUNCTIONS = 16


# ---------------------------------------------------------------------------
# the Cellular Potts Model (Glazier & Graner 1993; Savill & Hogeweg 1997)

class CPM:
    """A Potts lattice of cells, its energy and its Metropolis copy dynamics.

    ``bond(sigma_a, sigma_b)`` gives the surface energy J between two lattice
    states (0 is the medium); it is called only for sigma_a != sigma_b and the
    result is cached until :meth:`invalidate` is called.
    """

    def __init__(self, width: int, height: int, *, temperature: float,
                 inelasticity: float, bond, offset: float = DELTA_H_OFFSET):
        if width < 3 or height < 3:
            raise ValueError(f"the lattice needs at least 3 sites a side, got {width}x{height}")
        self.width, self.height = width, height
        self.size = width * height
        self.temperature = float(temperature)
        self.inelasticity = float(inelasticity)
        self.offset = float(offset)
        self.bond = bond
        self.sigma = [MEDIUM] * self.size
        self.volume: dict[int, int] = {}
        self.target: dict[int, int] = {}
        self.occupied: set[int] = set()
        self._cache: dict[tuple[int, int], float] = {}
        self.nbrs = [
            tuple((r + dr) * width + (c + dc)
                  for dr, dc in OFFSETS
                  if 0 <= r + dr < height and 0 <= c + dc < width)
            for r in range(height) for c in range(width)
        ]
        self.half_nbrs = [
            tuple((r + dr) * width + (c + dc)
                  for dr, dc in HALF_OFFSETS
                  if 0 <= r + dr < height and 0 <= c + dc < width)
            for r in range(height) for c in range(width)
        ]
        #: copies only ever target interior sites, so the border stays medium
        self.sites = [r * width + c
                      for r in range(1, height - 1) for c in range(1, width - 1)]

    # -- energies ----------------------------------------------------------
    def invalidate(self) -> None:
        """Forget the cached J values (their cell states have changed)."""
        self._cache.clear()

    def J(self, a: int, b: int) -> float:
        if a == b:
            return 0.0
        value = self._cache.get((a, b))
        if value is None:
            value = float(self.bond(a, b))
            self._cache[(a, b)] = self._cache[(b, a)] = value
        return value

    def hamiltonian(self) -> float:
        """H = sum of J over unordered neighbour pairs + lambda sum (v - V)^2."""
        total = 0.0
        sigma, J = self.sigma, self.J
        for x in range(self.size):
            sx = sigma[x]
            for y in self.half_nbrs[x]:
                sy = sigma[y]
                if sy != sx:
                    total += J(sx, sy)
        for cell, v in self.volume.items():
            d = v - self.target[cell]
            total += self.inelasticity * d * d
        return total

    def delta_h(self, site: int, new: int) -> float:
        """The change in H if ``site`` were taken over by cell ``new``."""
        old = self.sigma[site]
        if old == new:
            return 0.0
        dh = 0.0
        sigma, J = self.sigma, self.J
        for y in self.nbrs[site]:
            sy = sigma[y]
            if sy != new:
                dh += J(new, sy)
            if sy != old:
                dh -= J(old, sy)
        lam = self.inelasticity
        if new != MEDIUM:                      # (v+1-V)^2 - (v-V)^2
            dh += lam * (2 * (self.volume[new] - self.target[new]) + 1)
        if old != MEDIUM:                      # (v-1-V)^2 - (v-V)^2
            dh += lam * (1 - 2 * (self.volume[old] - self.target[old]))
        return dh

    def acceptance(self, dh: float) -> float:
        if dh < -self.offset:
            return 1.0
        return math.exp(-(dh + self.offset) / self.temperature)

    # -- dynamics ----------------------------------------------------------
    def sweep(self, rng) -> int:
        """One time step: one copy attempt per interior site. Returns the copies."""
        sites, nbrs, sigma = self.sites, self.nbrs, self.sigma
        volume, occupied = self.volume, self.occupied
        n = len(sites)
        picks = rng.integers(0, n, n).tolist()
        dirs = rng.integers(0, len(OFFSETS), n).tolist()
        draws = rng.random(n).tolist()
        temperature, offset = self.temperature, self.offset
        copies = 0
        for a, d, u in zip(picks, dirs, draws):
            x = sites[a]
            new = sigma[nbrs[x][d]]
            old = sigma[x]
            if new == old:
                continue
            dh = self.delta_h(x, new)
            if dh >= -offset and u >= math.exp(-(dh + offset) / temperature):
                continue
            sigma[x] = new
            if old != MEDIUM:
                volume[old] -= 1
            else:
                occupied.add(x)
            if new != MEDIUM:
                volume[new] += 1
            else:
                occupied.discard(x)
            copies += 1
        return copies

    # -- cells -------------------------------------------------------------
    def place(self, cell: int, sites, target: int) -> None:
        for x in sites:
            self.sigma[x] = cell
            self.occupied.add(x)
        self.volume[cell] = len(sites)
        self.target[cell] = target

    def pixels(self, cell: int) -> list[int]:
        return [x for x in self.occupied if self.sigma[x] == cell]

    def remove(self, cell: int) -> None:
        for x in self.pixels(cell):
            self.sigma[x] = MEDIUM
            self.occupied.discard(x)
        self.volume.pop(cell, None)
        self.target.pop(cell, None)

    def divide(self, cell: int, daughter: int) -> bool:
        """Split ``cell`` across the middle of and perpendicular to its longest axis."""
        px = self.pixels(cell)
        if len(px) < 2:
            return False
        rows = np.array([x // self.width for x in px], dtype=float)
        cols = np.array([x % self.width for x in px], dtype=float)
        rows -= rows.mean()
        cols -= cols.mean()
        cov = np.cov(np.vstack([rows, cols]))
        values, vectors = np.linalg.eigh(cov)
        axis = vectors[:, int(np.argmax(values))]          # the longest axis
        projection = rows * axis[0] + cols * axis[1]
        order = np.argsort(projection, kind="stable")
        moved = [px[int(k)] for k in order[len(px) // 2:]]
        for x in moved:
            self.sigma[x] = daughter
        self.volume[cell] -= len(moved)
        self.volume[daughter] = len(moved)
        self.target[daughter] = self.target[cell]
        return True

    # -- observables -------------------------------------------------------
    def contacts(self) -> dict[int, set[int]]:
        """Which cells touch which (the medium excluded)."""
        adj: dict[int, set[int]] = {c: set() for c in self.volume}
        sigma = self.sigma
        for x in self.occupied:
            sx = sigma[x]
            for y in self.nbrs[x]:
                sy = sigma[y]
                if sy != sx and sy != MEDIUM:
                    adj[sx].add(sy)
        return adj

    def contact_counts(self, type_of) -> tuple[int, int]:
        """(cell-cell neighbour pairs, those joining two different cell types)."""
        sigma = self.sigma
        total = hetero = 0
        for x in self.occupied:
            sx = sigma[x]
            for y in self.half_nbrs[x]:
                sy = sigma[y]
                if sy != MEDIUM and sy != sx:
                    total += 1
                    if type_of[sx] != type_of[sy]:
                        hetero += 1
        return total, hetero

    def surface(self) -> int:
        """Neighbour pairs between a cell and the medium."""
        sigma = self.sigma
        out = 0
        for x in self.occupied:
            for y in self.nbrs[x]:
                if sigma[y] == MEDIUM:
                    out += 1
        return out

    def centroid(self, cell: int) -> tuple[float, float]:
        px = self.pixels(cell)
        return (sum(x // self.width for x in px) / len(px),
                sum(x % self.width for x in px) / len(px))


# ---------------------------------------------------------------------------
# the gene regulation network (Hogeweg 2000, section 2.2)

class Genome:
    """A Boolean network: two inputs and one of the 16 two-input functions per node.

    An input >= 0 is a node of the cell itself; an input -1 - k reads entry k
    of the environment vector, which is the OR over the neighbouring cells of
    their node k. Only the first ``signals`` entries carry a signal; the rest
    are a constant 0, which is how a variable functional connectivity and a
    genotypic redundancy arise.
    """

    __slots__ = ("inputs", "functions", "nodes", "signals")

    def __init__(self, inputs, functions, nodes: int, signals: int):
        self.inputs = inputs
        self.functions = functions
        self.nodes = nodes
        self.signals = signals

    @staticmethod
    def _draw(nodes: int, value: int) -> int:
        """0..nodes-1 -> the environment entry -1-k; nodes..2*nodes-1 -> own node."""
        return value - nodes if value >= nodes else -1 - value

    @classmethod
    def random(cls, nodes: int, signals: int, rng) -> Genome:
        raw = rng.integers(0, 2 * nodes, size=(nodes, 2)).tolist()
        inputs = [[cls._draw(nodes, int(v)) for v in row] for row in raw]
        functions = [int(f) for f in rng.integers(0, FUNCTIONS, size=nodes).tolist()]
        return cls(inputs, functions, nodes, signals)

    def mutate(self, rng) -> Genome:
        """One point mutation: a new input for a node, or a new Boolean function."""
        inputs = [list(row) for row in self.inputs]
        functions = list(self.functions)
        i = int(rng.integers(self.nodes))
        if rng.random() < 0.5:
            slot = int(rng.integers(2))
            inputs[i][slot] = self._draw(self.nodes, int(rng.integers(2 * self.nodes)))
        else:
            functions[i] = int(rng.integers(FUNCTIONS))
        return Genome(inputs, functions, self.nodes, self.signals)

    def step(self, state: tuple[int, ...], env: tuple[int, ...]) -> tuple[int, ...]:
        """One synchronous update of every node."""
        out = []
        for function, (a, b) in zip(self.functions, self.inputs):
            va = state[a] if a >= 0 else env[-1 - a]
            vb = state[b] if b >= 0 else env[-1 - b]
            out.append((function >> (va * 2 + vb)) & 1)
        return tuple(out)

    def to_json(self) -> dict:
        return {"nodes": self.nodes, "signal_nodes": self.signals,
                "inputs": [list(row) for row in self.inputs],
                "functions": list(self.functions),
                "encoding": ("input >= 0 is the cell's own node; input -1-k is entry k of the "
                             "environment vector (the OR over neighbouring cells of their node k), "
                             "and is a constant 0 for k >= signal_nodes. A function is the 4-bit "
                             "truth table read as bit (2 a + b)")}


class Adhesion:
    """Lock-and-key surface receptors (Hogeweg 2000 eqs. 2a-2b).

    The first ``adhesion // 2`` adhesion nodes are locks, the rest keys; key k
    of one cell matches lock k of the other when they are complementary, and
    the matches are summed with binary place value.
    """

    def __init__(self, adhesion: int):
        self.half = adhesion // 2
        self.locks = tuple(range(self.half))
        self.keys = tuple(range(self.half, 2 * self.half))
        start = self.half - self.half // 2          # straddles locks and keys
        self.medium_nodes = tuple(range(start, start + self.half))
        self.maximum = float((1 << self.half) - 1)
        self._between: dict[tuple, float] = {}
        self._medium: dict[tuple, float] = {}

    def _match(self, i: tuple[int, ...], j: tuple[int, ...]) -> int:
        """sum_k 2^k [key k of i complements lock k of j]."""
        return sum(1 << k for k in range(self.half) if i[self.keys[k]] != j[self.locks[k]])

    def between(self, a: tuple[int, ...], b: tuple[int, ...]) -> float:
        value = self._between.get((a, b))
        if value is None:
            value = (self._match(a, b) + self._match(b, a)) / 2.0
            self._between[(a, b)] = self._between[(b, a)] = value
        return value

    def medium(self, a: tuple[int, ...]) -> float:
        value = self._medium.get(a)
        if value is None:
            value = float(sum(1 << k for k, node in enumerate(self.medium_nodes) if a[node]))
            self._medium[a] = value
        return value


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x != y for x, y in zip(a, b))


def differentiation(patterns) -> int:
    """The fitness: summed Hamming distance between all distinct cell types."""
    distinct = sorted(set(patterns))
    return sum(hamming(distinct[i], distinct[j])
               for i in range(len(distinct)) for j in range(i + 1, len(distinct)))


# ---------------------------------------------------------------------------
# development

class Development:
    """One critter: a lattice, its cells, their GRN states and what happened."""

    def __init__(self, genome: Genome, p, rng, *, growth: bool):
        self.genome = genome
        self.rng = rng
        self.growth = growth
        self.reference = int(p.target_volume)
        self.growth_threshold = int(p.growth_threshold)
        self.divisions = int(p.divisions)
        self.steps_per_stage = int(p.steps_per_stage)
        self.final_steps = int(p.final_steps)
        self.max_cycle = int(p.max_cycle)
        self.fitness_window = int(p.fitness_window)
        self.adhesion = Adhesion(int(p.adhesion_nodes))
        self.maternal = (genome.nodes - 5, genome.nodes - 4)
        self.state: dict[int, tuple[int, ...]] = {}
        self.history: dict[int, list[tuple[int, ...]]] = {}
        self.events: dict[tuple[tuple[int, ...], tuple], int] = {}
        self.per_step: list[dict] = []
        self.samples: list[int] = []
        self.next_id = 1
        self.step_count = 0
        self.cpm = CPM(int(p.grid), int(p.grid), temperature=p.temperature,
                       inelasticity=p.inelasticity, bond=self._bond)
        self._seed_zygote()

    # -- setup -------------------------------------------------------------
    def _bond(self, a: int, b: int) -> float:
        if a == MEDIUM:
            return self.adhesion.medium(self.state[b])
        if b == MEDIUM:
            return self.adhesion.medium(self.state[a])
        return self.adhesion.between(self.state[a], self.state[b])

    def _seed_zygote(self) -> None:
        side = max(1, round(math.sqrt(self.reference)))
        r0 = (self.cpm.height - side) // 2
        c0 = (self.cpm.width - side) // 2
        block = [(r0 + dr) * self.cpm.width + (c0 + dc)
                 for dr in range(side) for dc in range(side)]
        zygote = self.next_id
        self.next_id += 1
        self.cpm.place(zygote, block, self.reference)
        # "the zygote is initiated in state 0 (no genes expressed)"
        self.state[zygote] = (0,) * self.genome.nodes
        self.history[zygote] = [self.state[zygote]]

    # -- one time step -----------------------------------------------------
    def _environment(self) -> dict[int, tuple[int, ...]]:
        nodes, signals = self.genome.nodes, self.genome.signals
        blank = (0,) * nodes
        out = {}
        for cell, neighbours in self.cpm.contacts().items():
            if not neighbours or signals == 0:
                out[cell] = blank
                continue
            env = [0] * nodes
            for other in neighbours:
                other_state = self.state[other]
                for k in range(signals):
                    env[k] |= other_state[k]
            out[cell] = tuple(env)
        return out

    def _record(self, before: tuple[int, ...], after: tuple) -> None:
        key = (before, after)
        self.events[key] = self.events.get(key, 0) + 1

    def _update_genes(self) -> None:
        env = self._environment()
        blank = (0,) * self.genome.nodes
        updated = {}
        for cell, state in self.state.items():
            new = self.genome.step(state, env.get(cell, blank))
            updated[cell] = new
            if new != state:
                self._record(state, (new,))
            window = self.history[cell]
            window.append(new)
            if len(window) > self.max_cycle + 1:
                del window[0]
        self.state = updated
        self.cpm.invalidate()

    def _settled(self, cell: int) -> bool:
        """True when the cell's expression pattern repeats with period <= max_cycle."""
        window = self.history[cell]
        return window[-1] in window[:-1]

    def _reap(self) -> None:
        for cell in [c for c, v in self.cpm.volume.items() if v == 0]:
            self._record(self.state[cell], ())
            self.cpm.remove(cell)
            del self.state[cell], self.history[cell]

    def _grow(self) -> None:
        for cell in list(self.cpm.volume):
            target = self.cpm.target[cell]
            if self.cpm.volume[cell] > target + self.growth_threshold:
                target += 1
                self.cpm.target[cell] = target
            if target >= 2 * self.reference:
                self.cpm.target[cell] = max(1, target // 2)
                self._cleave(cell, self.cpm.target[cell])

    def _cleave(self, cell: int, target: int) -> None:
        daughter = self.next_id
        if not self.cpm.divide(cell, daughter):
            return
        self.next_id += 1
        self.cpm.target[cell] = self.cpm.target[daughter] = target
        self.state[daughter] = self.state[cell]
        self.history[daughter] = list(self.history[cell])
        self._record(self.state[cell], (self.state[cell], self.state[cell]))

    def _sample(self) -> None:
        cells = sorted(self.cpm.volume)
        patterns = [self.state[c] for c in cells if self._settled(c)]
        self.samples.append(differentiation(patterns))
        total, hetero = self.cpm.contact_counts(self.state)
        self.per_step.append({
            "step": self.step_count,
            "cells": len(cells),
            "cell_types": len(set(self.state[c] for c in cells)),
            "settled_types": len(set(patterns)),
            "differentiation": self.samples[-1],
            "heterotypic_contact_fraction": (hetero / total) if total else 0.0,
            "surface": self.cpm.surface(),
            "volume": int(sum(self.cpm.volume.values())),
        })

    def step(self) -> None:
        self.step_count += 1
        self._update_genes()
        self.cpm.sweep(self.rng)
        self._reap()
        if self.growth:
            self._grow()
        self._sample()

    # -- the whole development --------------------------------------------
    def run(self) -> None:
        for stage in range(self.divisions + 1):
            steps = self.final_steps if stage == self.divisions else self.steps_per_stage
            for _ in range(steps):
                self.step()
                if not self.cpm.volume:                    # every cell died
                    return
            if stage < self.divisions:
                self._cleavage(stage)

    def _cleavage(self, stage: int) -> None:
        """A pre-scheduled cleavage: every cell divides at once."""
        mothers = sorted(self.cpm.volume)
        born = []
        for cell in mothers:
            before = self.next_id
            self._cleave(cell, self.reference)
            born.append(self.next_id - 1 if self.next_id > before else None)
        # "after the first and the second cell division unequal cell division can
        # occur ... by flipping the state of a pre-defined node during 1 time-step
        # in one of the daughters" - the maternal factors
        if stage < len(self.maternal):
            index = min(stage, len(born) - 1)
            daughter = born[index]
            if daughter is not None:
                node = self.maternal[stage]
                state = list(self.state[daughter])
                state[node] ^= 1
                before = self.state[daughter]
                self.state[daughter] = tuple(state)
                self._record(before, (self.state[daughter],))
        self.cpm.invalidate()

    @property
    def fitness(self) -> int:
        """The minimum differentiation seen over the last fitness_window steps."""
        if not self.samples:
            return 0
        return min(self.samples[-self.fitness_window:])


# ---------------------------------------------------------------------------
# the generator

def _check_sorting_energies(value) -> dict[str, float]:
    keys = ("aa", "ab", "bb", "am", "bm")
    if not isinstance(value, dict):
        raise ValueError(f"sorting_energies must be an object with the keys {list(keys)}, got {value!r}")
    missing = [k for k in keys if k not in value]
    if missing:
        raise ValueError(f"sorting_energies is missing {missing}; it needs {list(keys)}")
    out = {}
    for k in keys:
        v = value[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise ValueError(f"sorting_energies[{k!r}] must be a number, got {v!r}")
        if v < 0:
            raise ValueError(f"sorting_energies[{k!r}] must be >= 0 (J is a bond energy), got {v!r}")
        out[k] = float(v)
    return out


def _validate(p) -> None:
    if p.adhesion_nodes % 2:
        raise ValueError(
            f"adhesion_nodes must be even (half locks, half keys), got {p.adhesion_nodes}")
    if p.adhesion_nodes > p.grn_nodes - 5:
        raise ValueError(
            f"adhesion_nodes must leave room for the two maternal nodes: "
            f"adhesion_nodes {p.adhesion_nodes} > grn_nodes - 5 = {p.grn_nodes - 5}")
    if p.signal_nodes > p.adhesion_nodes:
        raise ValueError(
            f"signal_nodes must be at most adhesion_nodes (the signalling nodes are "
            f"adhesion nodes), got {p.signal_nodes} > {p.adhesion_nodes}")
    if p.grn_nodes < 6:
        raise ValueError(f"grn_nodes must be at least 6, got {p.grn_nodes}")


def generate(p, rng):
    _validate(p)
    if p.mode == "cell-sorting":
        return _cell_sorting(p, rng)
    return _evo_devo(p, rng)


# --- the Glazier & Graner cell-sorting experiment ---------------------------

def _cell_sorting(p, rng):
    energies = _check_sorting_energies(p.sorting_energies)
    n = int(p.sorting_cells)
    types = ["A"] * n + ["B"] * n
    rng.shuffle(types)
    type_of = {MEDIUM: None}

    def bond(a: int, b: int) -> float:
        ta, tb = type_of[a], type_of[b]
        if ta is None:
            return energies[f"{tb.lower()}m"]
        if tb is None:
            return energies[f"{ta.lower()}m"]
        return energies["ab"] if ta != tb else energies[f"{ta.lower() * 2}"]

    cpm = CPM(int(p.grid), int(p.grid), temperature=p.temperature,
              inelasticity=p.inelasticity, bond=bond)
    volume = int(p.sorting_volume)
    side = max(1, round(math.sqrt(volume)))
    across = max(1, round(math.sqrt(len(types))))
    r0 = (cpm.height - across * side) // 2
    c0 = (cpm.width - across * side) // 2
    if r0 < 1 or c0 < 1:
        raise ValueError(
            f"a {p.grid}x{p.grid} lattice is too small for {2 * n} cells of "
            f"{volume} sites; raise grid, or lower sorting_cells or sorting_volume")
    for k, kind in enumerate(types):
        br, bc = divmod(k, across)
        block = [(r0 + br * side + dr) * cpm.width + (c0 + bc * side + dc)
                 for dr in range(side) for dc in range(side)]
        cpm.place(k + 1, block, volume)
        type_of[k + 1] = kind

    trace = []
    deaths: dict[str, int] = {}
    for step in range(int(p.sorting_steps) + 1):
        total, hetero = cpm.contact_counts(type_of)
        trace.append({"step": step, "cells": len(cpm.volume),
                      "contacts": total, "heterotypic": hetero,
                      "heterotypic_contact_fraction": (hetero / total) if total else 0.0,
                      "surface": cpm.surface(),
                      "energy": cpm.hamiltonian()})
        if step == int(p.sorting_steps):
            break
        cpm.sweep(rng)
        for cell in [c for c, v in cpm.volume.items() if v == 0]:
            deaths[type_of[cell]] = deaths.get(type_of[cell], 0) + 1
            cpm.remove(cell)

    reactions = [Reaction({kind: 1}, {}, count=n) for kind, n in sorted(deaths.items())]
    alive = {"A": 0, "B": 0}
    for cell in cpm.volume:
        alive[type_of[cell]] += 1
    inner = _inner_type(cpm, type_of)
    return Network(
        species=[Species("A", structure="cell type A"), Species("B", structure="cell type B")],
        reactions=reactions,
        status="observed",
        initial_state={"A": float(n), "B": float(n)},
        extras={
            "space": _space(cpm, type_of),
            "compartments": _compartments(cpm, type_of),
            "energies": {
                "units": "dimensionless surface bond energy J; H is minimised",
                "hamiltonian": "H = sum over neighbouring site pairs with different sigma of "
                               "J(type, type) + lambda sum over cells of (v - V)^2",
                "acceptance": "P = 1 if dH < -0.1, else exp(-(dH + 0.1) / T)",
                "temperature": float(p.temperature),
                "inelasticity": float(p.inelasticity),
                "delta_h_offset": DELTA_H_OFFSET,
                "J": {"A-A": energies["aa"], "A-B": energies["ab"], "B-B": energies["bb"],
                      "A-medium": energies["am"], "B-medium": energies["bm"]},
                "engulfment": "Hogeweg (2000) sec. 3.1 after Glazier & Graner: type a engulfs "
                              "type b when J_ab < J_mb and J_am < J_bm",
            },
            "interaction_law": (
                "Cellular Potts Model: a random interior site takes the state of a random Moore "
                "neighbour with probability 1 if dH < -0.1 and exp(-(dH + 0.1) / T) otherwise; "
                "one time step is one copy attempt per site. Cell types are fixed here, so the "
                "only reaction is the death of a cell squeezed to zero volume."
            ),
            "analysis": {
                "experiment": "Glazier & Graner (1993) differential-adhesion cell sorting: a mixed "
                              "aggregate of two cell types sorts out",
                "steps": int(p.sorting_steps),
                "cells_alive": alive,
                "engulfed_type": inner,
                "heterotypic_contact_fraction": {
                    "initial": trace[0]["heterotypic_contact_fraction"],
                    "final": trace[-1]["heterotypic_contact_fraction"],
                },
                "energy": {"initial": trace[0]["energy"], "final": trace[-1]["energy"]},
                "per_step": trace,
            },
        },
    )


def _inner_type(cpm: CPM, type_of) -> str | None:
    """The type whose cells sit further from the medium: the engulfed one."""
    exposure: dict[str, list[float]] = {}
    for cell in cpm.volume:
        px = cpm.pixels(cell)
        touching = sum(1 for x in px for y in cpm.nbrs[x] if cpm.sigma[y] == MEDIUM)
        exposure.setdefault(type_of[cell], []).append(touching / len(px))
    if len(exposure) < 2:
        return None
    return min(exposure, key=lambda k: sum(exposure[k]) / len(exposure[k]))


# --- the Hogeweg evo-devo model ---------------------------------------------

def _evo_devo(p, rng):
    population = [Genome.random(int(p.grn_nodes), int(p.signal_nodes), rng)
                  for _ in range(int(p.population))]
    generations: list[dict] = []
    best_genome, best_fitness = population[0], -1

    for generation in range(int(p.generations) + 1):
        scored = []
        for genome in population:
            dev = Development(genome, p, rng, growth=False)
            dev.run()
            scored.append((dev.fitness, genome))
            if dev.fitness > best_fitness:
                best_fitness, best_genome = dev.fitness, genome
        fits = [f for f, _ in scored]
        generations.append({
            "generation": generation,
            "best_fitness": max(fits),
            "mean_fitness": sum(fits) / len(fits),
            "differentiating": sum(1 for f in fits if f > 0),
        })
        if generation == int(p.generations):
            break
        population = _breed(scored, p, rng)

    # "During evolution cell growth/division was disabled ... the development of
    # selected critters was studied with cell growth/division enabled."
    final = Development(best_genome, p, rng, growth=True)
    final.run()
    return _network(final, best_fitness, generations, p)


def _breed(scored, p, rng) -> list[Genome]:
    """Tournament selection plus point mutation; sterile winners are redrawn."""
    out = []
    size = min(int(p.tournament), len(scored))
    for _ in range(len(scored)):
        picks = rng.integers(0, len(scored), size).tolist()
        fitness, genome = max((scored[i] for i in picks), key=lambda s: s[0])
        if fitness == 0:
            # "as long as none of the networks shows differentiation, new networks
            # are generated randomly"
            out.append(Genome.random(int(p.grn_nodes), int(p.signal_nodes), rng))
        elif rng.random() < float(p.mutation_rate):
            out.append(genome.mutate(rng))
        else:
            out.append(genome)
    return out


def _network(dev: Development, best_fitness: int, generations, p) -> Network:
    names: dict[tuple[int, ...], str] = {}

    def name(pattern: tuple[int, ...]) -> str:
        if pattern not in names:
            names[pattern] = f"T{len(names)}"
        return names[pattern]

    zygote = (0,) * dev.genome.nodes
    name(zygote)                                   # T0 is always the zygote
    reactions = []
    for (before, after), count in dev.events.items():
        products: dict[str, int] = {}
        for pattern in after:
            key = name(pattern)
            products[key] = products.get(key, 0) + 1
        reactions.append(Reaction({name(before): 1}, products, count=count))
    for cell in dev.cpm.volume:                    # types alive at the end
        name(dev.state[cell])

    species = [Species(sid, structure="".join(str(b) for b in pattern))
               for pattern, sid in names.items()]
    type_of = {c: dev.state[c] for c in dev.cpm.volume}
    total, hetero = dev.cpm.contact_counts(type_of)
    final = {}
    for cell in dev.cpm.volume:
        key = name(dev.state[cell])
        final[key] = final.get(key, 0) + 1

    adhesion = dev.adhesion
    alive = sorted(set(dev.state[c] for c in dev.cpm.volume))
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={"T0": 1.0},
        extras={
            "space": _space(dev.cpm, {c: name(dev.state[c]) for c in dev.cpm.volume}),
            "compartments": _compartments(dev.cpm, {c: name(dev.state[c]) for c in dev.cpm.volume}),
            "energies": {
                "units": "dimensionless surface bond energy J, 0 to "
                         f"{adhesion.maximum:g}; H is minimised",
                "hamiltonian": "H = sum over neighbouring site pairs with different sigma of "
                               "J(type, type) + lambda sum over cells of (v - V)^2",
                "acceptance": "P = 1 if dH < -0.1, else exp(-(dH + 0.1) / T)",
                "temperature": float(p.temperature),
                "inelasticity": float(p.inelasticity),
                "delta_h_offset": DELTA_H_OFFSET,
                "receptors": {
                    "locks": list(adhesion.locks), "keys": list(adhesion.keys),
                    "medium_nodes": list(adhesion.medium_nodes),
                    "law": "J_ij = (sum_k 2^k [key k of i != lock k of j] + sum_k 2^k "
                           "[key k of j != lock k of i]) / 2; J_im = sum_l 2^l (node l of i)",
                },
                "J": {f"{name(a)}-{name(b)}": adhesion.between(a, b)
                      for i, a in enumerate(alive) for b in alive[i:]},
                "J_medium": {name(a): adhesion.medium(a) for a in alive},
                "engulfment": "Hogeweg (2000) sec. 3.1 after Glazier & Graner: type a engulfs "
                              "type b when J_ab < J_mb and J_am < J_bm",
            },
            "grn": dev.genome.to_json(),
            "maternal_nodes": list(dev.maternal),
            "events": ("Ta -> Tb: a cell differentiates; Ta -> 2 Ta: a cell divides and both "
                       "daughters inherit the mother's expression pattern; Ta -> (nothing): a "
                       "cell squeezed to zero volume dies"),
            "final_state": {k: float(v) for k, v in sorted(final.items())},
            "analysis": {
                "fitness": {
                    "law": "the summed Hamming distance between the expression patterns of the "
                           "distinct settled cell types, minimised over the last fitness_window "
                           "steps of development",
                    "best_evolved": int(best_fitness),
                    "final_development": int(dev.fitness),
                    "per_generation": generations,
                },
                "cells_alive": len(dev.cpm.volume),
                "cell_types": len(alive),
                "types_seen": len(names),
                "divisions": sum(c for (_, after), c in dev.events.items() if len(after) == 2),
                "deaths": sum(c for (_, after), c in dev.events.items() if not after),
                "differentiations": sum(c for (_, after), c in dev.events.items() if len(after) == 1),
                "heterotypic_contact_fraction": (hetero / total) if total else 0.0,
                "surface": dev.cpm.surface(),
                "volume": int(sum(dev.cpm.volume.values())),
                "steps": dev.step_count,
                "per_step": dev.per_step,
            },
        },
    )


# --- shared extras ----------------------------------------------------------

def _space(cpm: CPM, label) -> dict:
    return {
        "dimensions": 2,
        "shape": [cpm.width, cpm.height],
        "lattice": "Cellular Potts: a site carries the id of the cell occupying it, 0 = medium",
        "neighbourhood": "Moore (8 sites); the outermost ring of sites stays medium",
        "site_index": "row * width + column",
        "sigma": list(cpm.sigma),
        "cell_type": {str(c): label.get(c) for c in sorted(cpm.volume)},
    }


def _compartments(cpm: CPM, label) -> dict:
    cells = {}
    for cell in sorted(cpm.volume):
        row, col = cpm.centroid(cell)
        cells[str(cell)] = {
            "type": label.get(cell),
            "volume": int(cpm.volume[cell]),
            "target_volume": int(cpm.target[cell]),
            "centroid": [row, col],
        }
    return {
        "cells": cells,
        "medium": "sigma 0; it has no volume constraint",
        "constraint": "lambda (v - V)^2 per cell; a cell reduced to v = 0 dies",
    }
