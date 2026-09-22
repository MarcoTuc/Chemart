"""Autopoiesis of Varela, Maturana & Uribe (1974), in McMullin's SCL reconstruction.

Catalog id: autopoiesis-vmu. Book 6.1.5 and figure 6.3.

A 2D toroidal lattice holds one particle per site: hole, substrate (O),
catalyst (*) or link ([O]). A link may additionally hold one absorbed
substrate and may carry up to two bonds to neighbouring links, so links form
chains; a closed chain is a membrane, and a membrane with a catalyst inside is
a cell. The three reactions of the book are

    production       * + 2 O -> * + [O]      (one site is left empty: a hole)
    bonding          links bond to neighbouring links (at most 2 bonds each)
    disintegration   [O] -> O + O            (spontaneous)

plus the absorption and emission of substrate by a link, which is what makes a
membrane permeable to substrate while bonded links, being immobile, make it
impermeable to links and catalysts.

The rule that the 1974 paper omitted, and that McMullin & Varela (1997)
identified as essential, is **chain-based bond inhibition**: a free link may
not bond while a chain link (a doubly bonded link, i.e. a membrane
constituent) is in its neighbourhood. It keeps the links inside a cell free
and mobile, hence available to repair a rupture of the membrane; without it
they bond to each other and to the membrane and the cell cannot heal
(`bond_inhibition=False` reproduces that failure).

Two faces. ``generate`` returns the reactions the chemistry defines, as a
complete network. ``evolve`` runs the lattice and yields a frame per time step
(the species counts, the reactions fired in the step, and the membrane
measurements closed_chains, membranes and enclosed_catalysts); it returns the
reactions that fired, with counts (status observed), the final grid in
``extras["space"]`` and the whole-run membrane measurements in
``extras["analysis"]``.

Species: ``S`` substrate, ``C`` catalyst, ``L0``/``L1``/``L2`` a link with 0, 1
or 2 bonds, and ``L0S``/``L1S``/``L2S`` the same link holding an absorbed
substrate. Holes are empty space, not a species.
"""

from __future__ import annotations

import math
from collections import Counter

from chemart.network import Network, Reaction, Species
from chemart.trajectory import Frame

HOLE, SUBSTRATE, CATALYST, LINK = 0, 1, 2, 3

#: extras["space"]["grid"] codes.
CODE = {HOLE: ".", SUBSTRATE: "o", CATALYST: "*"}
LINK_CODE = {(0, 0): "L", (1, 0): "l", (2, 0): "=", (0, 1): "M", (1, 1): "m", (2, 1): "#"}

#: Mobility factors of Von Kamp (2002), table 1; bonded links are immobile.
DEFAULT_MOBILITY = {"substrate": 0.5, "catalyst": 0.1, "link": 0.1, "hole": 0.5}
MOBILITY_KEYS = tuple(DEFAULT_MOBILITY)

#: A twelve-link closed chain around one site, the cell of Von Kamp (2002),
#: figure 2; consecutive links are 8-adjacent and it encloses the 3x3 interior.
RING = ((2, 0), (2, 1), (1, 2), (0, 2), (-1, 2), (-2, 1),
        (-2, 0), (-2, -1), (-1, -2), (0, -2), (1, -2), (2, -1))

#: Motion directions, in the order of World.von_neumann.
VON_NEUMANN_DIRS = ((1, 0), (0, 1), (-1, 0), (0, -1))

SPECIES_STRUCTURE = {
    "S": "substrate O",
    "C": "catalyst *",
    "L0": "link [O], free (no bonds)",
    "L1": "link [O] with one bond (chain end)",
    "L2": "link [O] with two bonds (chain link, membrane constituent)",
    "L0S": "free link holding an absorbed substrate",
    "L1S": "singly bonded link holding an absorbed substrate",
    "L2S": "chain link holding an absorbed substrate",
}


def link_species(bonds: int, loaded: bool) -> str:
    return f"L{bonds}S" if loaded else f"L{bonds}"


class _Uniform:
    """Buffered uniform [0, 1) draws from a numpy Generator."""

    def __init__(self, rng, block: int = 8192):
        self.rng, self.block = rng, block
        self.buf, self.i = rng.random(block).tolist(), 0

    def __call__(self) -> float:
        if self.i == self.block:
            self.buf, self.i = self.rng.random(self.block).tolist(), 0
        u = self.buf[self.i]
        self.i += 1
        return u


class _Ints:
    """Buffered integer draws in [0, n) for one fixed n."""

    def __init__(self, rng, n: int, block: int = 8192):
        self.rng, self.n, self.block = rng, n, block
        self.buf, self.i = rng.integers(n, size=block).tolist(), 0

    def __call__(self) -> int:
        if self.i == self.block:
            self.buf, self.i = self.rng.integers(self.n, size=self.block).tolist(), 0
        v = self.buf[self.i]
        self.i += 1
        return v


class World:
    """The SCL lattice: particles, bonds, and one time step of the algorithm.

    Sites are indexed ``i = y * width + x`` on a torus. Interactions (catalysis,
    bonding, absorption, emission, bond inhibition) use the 8-cell Moore
    neighbourhood; motion swaps with one of the 4 von Neumann neighbours, so an
    8-connected closed chain of links is impermeable to moving particles.
    """

    def __init__(self, rng, width: int = 30, height: int = 30, *,
                 production_probability: float = 1.0,
                 disintegration_probability: float = 0.01,
                 bond_probability: float = 1.0,
                 absorption_probability: float = 0.5,
                 emission_probability: float = 0.5,
                 bond_inhibition: bool = True,
                 chain_inhibit_count: int = 1,
                 mobility: dict | None = None,
                 record: bool = True):
        self.rng = rng
        self.width, self.height = width, height
        self.size = width * height
        self.p_production = production_probability
        self.p_disintegration = disintegration_probability
        self.p_bond = bond_probability
        self.p_absorption = absorption_probability
        self.p_emission = emission_probability
        self.bond_inhibition = bond_inhibition
        self.chain_inhibit_count = chain_inhibit_count
        self.mobility = dict(DEFAULT_MOBILITY if mobility is None else mobility)

        self.cell = [SUBSTRATE] * self.size
        self.absorbed = bytearray(self.size)
        self.bonds: dict[int, set[int]] = {}
        self.moore, self.von_neumann = _neighbourhoods(width, height)
        self.uniform = _Uniform(rng)
        self.pick8 = _Ints(rng, 8)
        self.pick4 = _Ints(rng, 4)
        self.pick2 = _Ints(rng, 2)

        self.time = 0
        self.record = record
        self.fired: dict[tuple, list] = {}
        self.window: dict[tuple, list] = {}      # fired since the last flush()
        self.events = Counter()
        self.enclosure: list[bool] = []
        self.crossings = Counter()
        self._origin_side: dict[int, bool] = {}     # link site -> interior at absorption
        self._interior: frozenset[int] = frozenset()

    # ------------------------------------------------------------------
    # geometry and set-up
    def index(self, x: int, y: int) -> int:
        return (y % self.height) * self.width + (x % self.width)

    def put(self, i: int, kind: int) -> None:
        self.cell[i] = kind
        if kind != LINK:
            self.bonds.pop(i, None)
            self.absorbed[i] = 0
        else:
            self.bonds.setdefault(i, set())

    def bond(self, a: int, b: int) -> None:
        self.bonds[a].add(b)
        self.bonds[b].add(a)

    def n_bonds(self, i: int) -> int:
        return len(self.bonds.get(i, ()))

    def catalysts(self) -> list[int]:
        return [i for i, k in enumerate(self.cell) if k == CATALYST]

    def seed_catalysts(self, n: int) -> list[int]:
        """One catalyst at the centre (book figure 6.3, t = 0), the rest at random."""
        placed: list[int] = []
        if n > 0:
            placed.append(self.index(self.width // 2, self.height // 2))
        pick = _Ints(self.rng, self.size)
        while len(placed) < n:
            i = pick()
            if i not in placed:
                placed.append(i)
        for i in placed:
            self.put(i, CATALYST)
        return placed

    def seed_ring(self, centre: int) -> list[int]:
        """The twelve-link membrane of Von Kamp (2002), figure 2, around `centre`."""
        cx, cy = centre % self.width, centre // self.width
        sites = [self.index(cx + dx, cy + dy) for dx, dy in RING]
        for i in sites:
            self.put(i, LINK)
        for a, b in zip(sites, sites[1:] + sites[:1]):
            self.bond(a, b)
        return sites

    def scatter(self, kind: int, n: int, sites: list[int]) -> None:
        if n <= 0:
            return
        pick = _Ints(self.rng, len(sites)) if sites else None
        placed = 0
        while placed < n and pick is not None:
            i = sites[pick()]
            if self.cell[i] == SUBSTRATE:
                self.put(i, kind)
                placed += 1

    # ------------------------------------------------------------------
    # the observed network
    def _fire(self, reactants: Counter, products: Counter) -> None:
        if not self.record or reactants == products:
            return
        key = (frozenset(reactants.items()), frozenset(products.items()))
        for book in (self.fired, self.window):
            entry = book.setdefault(key, [dict(reactants), dict(products), 0])
            entry[2] += 1

    def flush(self) -> list[list]:
        """The reactions fired since the last flush, as a frame's `fired`; then forget them."""
        out = [[sorted(Counter(lhs).elements()), sorted(Counter(rhs).elements()), n]
               for lhs, rhs, n in self.window.values()]
        self.window = {}
        return out

    def census(self) -> dict[str, float]:
        """How many of each species the lattice holds; holes are not a species."""
        counts = Counter(self.species_at(i) for i in range(self.size) if self.cell[i] != HOLE)
        return {s: float(n) for s, n in sorted(counts.items())}

    def species_at(self, i: int) -> str:
        kind = self.cell[i]
        if kind == SUBSTRATE:
            return "S"
        if kind == CATALYST:
            return "C"
        return link_species(self.n_bonds(i), bool(self.absorbed[i]))

    # ------------------------------------------------------------------
    # one time step
    def step(self) -> None:
        order = self.rng.permutation(self.size).tolist()
        done = bytearray(self.size)
        for i in order:
            if done[i]:
                continue
            done[i] = 1
            j = self.move(i)
            if j is not None:
                done[j] = 1
                i = j
            kind = self.cell[i]
            if kind == CATALYST:
                self.produce(i)
            elif kind == LINK:
                if self.uniform() < self.p_disintegration:
                    self.disintegrate(i)
                    continue
                self.try_bond(i)
                self.transport(i)
        self.time += 1

    def mobility_of(self, i: int) -> float:
        kind = self.cell[i]
        if kind == LINK:
            return 0.0 if self.bonds.get(i) else self.mobility["link"]
        return self.mobility[("hole", "substrate", "catalyst")[kind]]

    def move(self, i: int) -> int | None:
        """Swap with a random von Neumann neighbour; return the new site, if moved."""
        j = self.von_neumann[i][self.pick4()]
        if self.cell[i] == self.cell[j] and self.absorbed[i] == self.absorbed[j]:
            return None
        m = self.mobility_of(i) * self.mobility_of(j)
        if m <= 0.0 or self.uniform() >= math.sqrt(m):
            return None
        self.cell[i], self.cell[j] = self.cell[j], self.cell[i]
        self.absorbed[i], self.absorbed[j] = self.absorbed[j], self.absorbed[i]
        for site in (i, j):
            if self.cell[site] == LINK:
                self.bonds.setdefault(site, set())
            else:
                self.bonds.pop(site, None)
        self.events["motion"] += 1
        if self._interior and (i in self._interior) != (j in self._interior):
            for site in (i, j):
                kind = self.cell[site]
                if kind == CATALYST:
                    self.crossings["catalyst"] += 1
                elif kind == LINK:
                    self.crossings["link"] += 1
        return j

    def produce(self, i: int) -> None:
        """* + 2 O -> * + [O]: two neighbouring substrates become a link and a hole."""
        if self.uniform() >= self.p_production:
            return
        found = []
        start = self.pick8()
        for k in range(8):
            j = self.moore[i][(start + k) % 8]
            if self.cell[j] == SUBSTRATE:
                found.append(j)
                if len(found) == 2:
                    break
        if len(found) < 2:
            return
        link, hole = found if self.pick2() else found[::-1]
        self.put(link, LINK)
        self.put(hole, HOLE)
        self.events["production"] += 1
        self._fire(Counter({"C": 1, "S": 2}), Counter({"C": 1, "L0": 1}))

    def disintegrate(self, i: int) -> None:
        """[O] -> O + O: the link loses its bonds and becomes substrate."""
        partners = sorted(self.bonds.get(i, ()))
        reactants = Counter({self.species_at(i): 1})
        for b in partners:
            reactants[self.species_at(b)] += 1
        loaded = bool(self.absorbed[i])
        for b in partners:
            self.bonds[b].discard(i)
        products = Counter()
        for b in partners:
            products[self.species_at(b)] += 1
        self.put(i, SUBSTRATE)
        released = 1
        for _ in range(1 + int(loaded)):
            hole = self.free_neighbour(i)
            if hole is None:
                break
            self.put(hole, SUBSTRATE)
            released += 1
        products["S"] += released
        self.events["disintegration"] += 1
        self.events["substrates_released"] += released
        self._fire(reactants, products)

    def free_neighbour(self, i: int) -> int | None:
        start = self.pick8()
        for k in range(8):
            j = self.moore[i][(start + k) % 8]
            if self.cell[j] == HOLE:
                return j
        return None

    def inhibited(self, a: int, b: int) -> bool:
        """Chain-based bond inhibition: a free link may not bond next to a chain link."""
        if not self.bond_inhibition:
            return False
        for x in (a, b):
            if self.bonds.get(x):
                continue                      # only free links are inhibited
            chain = sum(1 for n in self.moore[x] if self.n_bonds(n) == 2)
            if chain >= self.chain_inhibit_count:
                return True
        return False

    def try_bond(self, i: int) -> None:
        """Bond to one randomly chosen neighbouring link with a free bond site."""
        if len(self.bonds[i]) >= 2:
            return
        j = self.moore[i][self.pick8()]
        if self.cell[j] != LINK or len(self.bonds[j]) >= 2 or j in self.bonds[i]:
            return
        if self.inhibited(i, j) or self.uniform() >= self.p_bond:
            return
        # a list, not a dict: two links of the same species are two reactants
        reactants = Counter([self.species_at(i), self.species_at(j)])
        self.bond(i, j)
        products = Counter([self.species_at(i), self.species_at(j)])
        self.events["bond"] += 1
        self._fire(reactants, products)

    def transport(self, i: int) -> None:
        """Absorb a neighbouring substrate, or emit the absorbed one into a hole."""
        if not self.absorbed[i]:
            if self.uniform() >= self.p_absorption:
                return
            j = self.moore[i][self.pick8()]
            if self.cell[j] != SUBSTRATE:
                return
            reactants = Counter({self.species_at(i): 1, "S": 1})
            self.put(j, HOLE)
            self.absorbed[i] = 1
            self.events["absorption"] += 1
            self._origin_side[i] = j in self._interior
            self._fire(reactants, Counter({self.species_at(i): 1}))
        else:
            if self.uniform() >= self.p_emission:
                return
            j = self.moore[i][self.pick8()]
            if self.cell[j] != HOLE:
                return
            reactants = Counter({self.species_at(i): 1})
            self.absorbed[i] = 0
            self.put(j, SUBSTRATE)
            self.events["emission"] += 1
            if self._interior and self._origin_side.get(i) is not None:
                if self._origin_side.pop(i) != (j in self._interior):
                    self.crossings["substrate"] += 1
            self._fire(reactants, Counter({self.species_at(i): 1, "S": 1}))

    # ------------------------------------------------------------------
    # membrane analysis
    def closed_chains(self) -> list[list[int]]:
        """Every closed chain of links (a cycle in which each link has two bonds)."""
        seen: set[int] = set()
        cycles = []
        for i, partners in self.bonds.items():
            if i in seen or len(partners) != 2:
                continue
            chain, previous, current = [i], i, min(partners)
            seen.add(i)
            while current != i and len(self.bonds.get(current, ())) == 2:
                chain.append(current)
                seen.add(current)
                nxt = [n for n in self.bonds[current] if n != previous]
                previous, current = current, nxt[0]
            if current == i and len(chain) > 2:
                cycles.append(chain)
            else:
                seen.update(chain)
        return cycles

    def interior_of(self, start: int, walls: set[int]) -> frozenset[int] | None:
        """The 4-connected region of `start` with `walls` blocked, or None.

        None means `start` is not enclosed. The region is followed on the
        unwrapped plane: one that is not enclosed wraps around the torus and so
        reaches the same site at two different offsets.
        """
        if start in walls:
            return None
        offset = {start: (0, 0)}
        stack = [(start, 0, 0)]
        while stack:
            i, ox, oy = stack.pop()
            for (dx, dy), j in zip(VON_NEUMANN_DIRS, self.von_neumann[i]):
                if j in walls:
                    continue
                here = (ox + dx, oy + dy)
                if j in offset:
                    if offset[j] != here:
                        return None
                    continue
                offset[j] = here
                stack.append((j, here[0], here[1]))
        return frozenset(offset)

    def enclosed(self) -> tuple[list[list[int]], dict[int, frozenset[int]]]:
        """Closed chains, and the interior of each catalyst that one encloses."""
        cycles = self.closed_chains()
        walls = {i for chain in cycles for i in chain}
        interiors: dict[int, frozenset[int]] = {}
        if walls:
            for c in self.catalysts():
                region = self.interior_of(c, walls)
                if region is not None:
                    interiors[c] = region
        return cycles, interiors

    def observe(self) -> dict[str, int]:
        """Measure the membranes now, as a frame's observables, and refresh the
        interior used to count crossings."""
        cycles, interiors = self.enclosed()
        self.enclosure.append(bool(interiors))
        self._interior = next(iter(interiors.values()), frozenset())
        return {
            "closed_chains": len(cycles),
            # Von Kamp (2002), 2.2: a closed chain of six links or more is a membrane
            "membranes": sum(1 for c in cycles if len(c) >= 6),
            "enclosed_catalysts": len(interiors),
        }

    def run(self, steps: int) -> None:
        self.observe()
        for _ in range(steps):
            self.step()
            self.observe()


def _neighbourhoods(width: int, height: int):
    moore, von_neumann = [], []
    for y in range(height):
        for x in range(width):
            moore.append(tuple(
                ((y + dy) % height) * width + (x + dx) % width
                for dx, dy in ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))
            ))
            von_neumann.append(tuple(
                ((y + dy) % height) * width + (x + dx) % width
                for dx, dy in ((1, 0), (0, 1), (-1, 0), (0, -1))
            ))
    return moore, von_neumann


# ----------------------------------------------------------------------------
# the defined network (the generate face)
def defined_reactions() -> list[Reaction]:
    out = [Reaction({"C": 1, "S": 2}, {"C": 1, "L0": 1})]
    for a, b in ((0, 0), (0, 1), (1, 1)):
        left = Counter({f"L{a}": 1})
        left[f"L{b}"] += 1
        right = Counter({f"L{a + 1}": 1})
        right[f"L{b + 1}"] += 1
        out.append(Reaction(dict(left), dict(right)))
    for k in range(3):
        out.append(Reaction({f"L{k}": 1}, {"S": 2}))
        out.append(Reaction({f"L{k}S": 1}, {"S": 3}))
        out.append(Reaction({f"L{k}": 1, "S": 1}, {f"L{k}S": 1}))
        out.append(Reaction({f"L{k}S": 1}, {f"L{k}": 1, "S": 1}))
    return out


# ----------------------------------------------------------------------------
def _mobility(value) -> dict[str, float]:
    if not isinstance(value, dict):
        raise ValueError(f"mobility must be an object with the keys {list(MOBILITY_KEYS)}, got {value!r}")
    unknown = sorted(set(value) - set(MOBILITY_KEYS))
    if unknown:
        raise ValueError(f"mobility: unknown particle type(s) {unknown}; valid keys: {list(MOBILITY_KEYS)}")
    out = dict(DEFAULT_MOBILITY)
    for key, v in value.items():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not 0.0 <= v <= 1.0:
            raise ValueError(f"mobility[{key!r}] must be a number in [0, 1], got {v!r}")
        out[key] = float(v)
    return out


def _species() -> list[Species]:
    return [Species(i, structure=s) for i, s in SPECIES_STRUCTURE.items()]


def generate(p, rng):
    """The reactions the chemistry defines (book 6.1.5), as a complete network."""
    return Network(species=_species(), reactions=defined_reactions(), status="complete",
                   extras={"rules": _RULES})


def evolve(p, rng):
    """Run the SCL lattice (McMullin's reconstruction): a frame per time step."""
    mobility = _mobility(p.mobility)
    if p.n_catalysts > p.width * p.height:
        raise ValueError(f"n_catalysts = {p.n_catalysts} exceeds the {p.width * p.height} sites of the lattice")
    world = World(
        rng, p.width, p.height,
        production_probability=p.production_probability,
        disintegration_probability=p.disintegration_probability,
        bond_probability=p.bond_probability,
        absorption_probability=p.absorption_probability,
        emission_probability=p.emission_probability,
        bond_inhibition=p.bond_inhibition,
        chain_inhibit_count=p.chain_inhibit_count,
        mobility=mobility,
    )
    holes = int(round((1.0 - p.substrate_density) * world.size))
    world.scatter(HOLE, holes, list(range(world.size)))
    centres = world.seed_catalysts(p.n_catalysts)
    if p.initial == "cell" and centres:
        world.seed_ring(centres[0])
    initial = world.census()

    yield Frame(t=0.0, state=initial, observables=world.observe())
    for _ in range(p.steps):
        world.step()
        yield Frame(t=float(world.time), state=world.census(), fired=world.flush(),
                    observables=world.observe())
    return _network(world, initial, p)


_RULES = (
    "Lattice rules (McMullin's SCL reconstruction). Each time step every particle, "
    "in random order: (1) swaps with one random von Neumann neighbour with probability "
    "sqrt(m_i m_j) of the two mobility factors, bonded links being immobile; (2) a "
    "catalyst turns two substrates of its Moore neighbourhood into a link and a hole; "
    "(3) a link disintegrates, losing its bonds and becoming substrate; else it bonds to "
    "one random neighbouring link that has a free bond site, unless chain-based bond "
    "inhibition applies (a free link with a chain link in its Moore neighbourhood may not "
    "bond), and then absorbs a neighbouring substrate or emits the absorbed one into a hole."
)


def _network(world: World, initial: dict[str, float], p) -> Network:
    reactions = [
        Reaction(reactants, products, count=count)
        for reactants, products, count in world.fired.values()
    ]
    grid = []
    for y in range(world.height):
        row = ""
        for x in range(world.width):
            i = y * world.width + x
            kind = world.cell[i]
            row += (LINK_CODE[(world.n_bonds(i), int(world.absorbed[i]))]
                    if kind == LINK else CODE[kind])
        grid.append(row)
    cycles, interiors = world.enclosed()
    enclosure = world.enclosure
    first = enclosure.index(True) if True in enclosure else None
    ruptures = sum(1 for a, b in zip(enclosure, enclosure[1:]) if a and not b)
    repairs = sum(
        1 for t in range(1, len(enclosure))
        if enclosure[t] and not enclosure[t - 1] and first is not None and t > first
    )
    return Network(
        species=_species(),
        reactions=reactions,
        status="observed",
        initial_state=initial,
        extras={
            "space": {
                "dimensions": 2,
                "shape": [world.width, world.height],
                "topology": "torus",
                "neighbourhood": {"interaction": "moore-8", "motion": "von-neumann-4"},
                "grid": grid,
                "legend": {".": "hole", "o": "substrate", "*": "catalyst", "L": "free link",
                           "l": "link with one bond", "=": "chain link (two bonds)",
                           "M": "free link + absorbed substrate",
                           "m": "one-bond link + absorbed substrate",
                           "#": "chain link + absorbed substrate"},
                "bonds": sorted([a, b] for a, partners in world.bonds.items() for b in partners if a < b),
                "mobility": world.mobility,
                "site_index": "i = y * width + x",
            },
            "analysis": {
                "steps": p.steps,
                "events": {k: int(v) for k, v in sorted(world.events.items())},
                "membrane": {
                    "closed_chains": len(cycles),
                    "chain_lengths": sorted(len(c) for c in cycles),
                    # Von Kamp (2002), 2.2: a closed chain of fewer than six
                    # links is a cluster, one of six or more a membrane.
                    "membranes": sum(1 for c in cycles if len(c) >= 6),
                    "clusters": sum(1 for c in cycles if len(c) < 6),
                    "enclosed_catalysts": len(interiors),
                    "interior_sizes": sorted(len(r) for r in interiors.values()),
                    "first_enclosure_step": first,
                    "steps_enclosed": sum(enclosure),
                    "ruptures": ruptures,
                    "repairs": repairs,
                },
                "permeability": {
                    "substrate_crossings": int(world.crossings["substrate"]),
                    "link_crossings": int(world.crossings["link"]),
                    "catalyst_crossings": int(world.crossings["catalyst"]),
                },
            },
            "rules": _RULES,
        },
    )
