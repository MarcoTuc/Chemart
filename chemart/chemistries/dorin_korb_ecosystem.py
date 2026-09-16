"""Dorin & Korb's virtual ecosystem chemistry (book 8.2.3, ref [240]).

Virtual square atoms from {A, B, C, O} plus four catalyst types float on a 2D
grid, bond covalently edge to edge, and every bond made or broken moves energy
around. Reconstructed from Dorin & Korb, *Building Virtual Ecosystems from
Artificial Chemistry* (Monash TR 2007/212, the full version of the ECAL 2007
paper): table 1 is the bond table, the appendix gives the electron shells that
fix the valences, and sections 2.1-2.2 give the reactor.

The three reactions the paper writes out are consequences of the bond table,
not extra rules:

    AO + BO --(chlorophyll & sunlight)--> AB + 2 O    photosynthesis (3.1.1)
    O + AB  --(enzyme)--> A + BO + energy             respiration    (3.1.2)
    C + C   --(energy)--> C2                          biosynthesis   (3.1.3)

Energy convention (table 1 caption): a bond's energy is released when the bond
is made and must be supplied to break it, so the negative entries (A-B, C-C)
are the bonds that *store* energy. Released energy is pooled in the
spatially-connected cluster of atoms around the reaction site, is spendable
only in the same time step, and is otherwise lost.

The returned network is the set of reactions that actually fired in one run,
with firing counts.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter

from chemart.network import Network, Reaction, Species

#: Appendix "Known atoms": atom : shell 1 : shell 2 : shell 3.
SHELLS = {
    "A": (1,), "B": (2, 3), "O": (2, 1), "C": (2, 4, 4), "K": (2, 2),
    "EAB": (2, 2), "ECC": (2, 2), "EO": (2, 2),
}
#: Appendix bonding rule 3: electron shells from inner to outer hold 2, 4, 8.
CAPACITY = (2, 4, 8)

#: The four catalyst types of section 2.1 ("four types of catalyst are required").
CATALYSTS = {
    "K": "chlorophyll: makes A-B from sunlight, breaks A-O and B-O",
    "EAB": "sugar-breaking enzyme: breaks A-B and releases its stored energy",
    "ECC": "organic decomposer: breaks C-C, makes C-O",
    "EO": "inorganic decomposer: breaks A-O and B-O",
}
ATOM_ORDER = ("A", "B", "C", "O", "K", "EAB", "ECC", "EO")
BUILDING_BLOCKS = ("A", "B", "C", "O")

#: Table 1, verbatim: bond -> (make, break, catalysed make, catalysed break,
#: bond energy as (sign, magnitude)). "-" is written as None.
TABLE = {
    "A-B": ("low", "low", {"K": "high"}, {"EAB": "high"}, (-1, "high")),
    "C-C": ("moderate", "low", None, {"ECC": "high"}, (-1, "low")),
    "A-O": ("high", "low", None, {"K": "high", "EO": "high"}, (+1, "low")),
    "B-O": ("high", "low", None, {"K": "high", "EO": "high"}, (+1, "low")),
    "C-O": ("low", "moderate", {"ECC": "high"}, None, (+1, "low")),
}
MAKE, BREAK, MAKE_CAT, BREAK_CAT, ENERGY = range(5)

DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def valence(atom: str) -> int:
    """Appendix rule 2: a bond needs a free outer-shell electron *and* a free slot."""
    shells = SHELLS[atom]
    electrons = shells[-1]
    return min(electrons, CAPACITY[len(shells) - 1] - electrons)


VALENCE = {a: valence(a) for a in SHELLS}


def bond_key(a: str, b: str) -> str:
    return "-".join(sorted((a, b)))


# ---------------------------------------------------------------------------
# Canonical molecule identity: the bond graph, up to isomorphism.
# ---------------------------------------------------------------------------
def _refine(colour: dict[int, int], adj: dict[int, set[int]]) -> dict[int, int]:
    """Weisfeiler-Lehman colour refinement to a stable partition."""
    while True:
        sig = {v: (colour[v], tuple(sorted(colour[u] for u in adj[v]))) for v in adj}
        rank = {s: i for i, s in enumerate(sorted(set(sig.values())))}
        new = {v: rank[sig[v]] for v in adj}
        if len(set(new.values())) == len(set(colour.values())):
            return new
        colour = new


def _code(colour: dict[int, int], adj: dict[int, set[int]], types: list[str]) -> str:
    order = sorted(colour, key=lambda v: colour[v])
    bits = "".join(
        "1" if order[j] in adj[order[i]] else "0"
        for i in range(len(order)) for j in range(i + 1, len(order))
    )
    return "".join(types[v] for v in order) + (f"|{bits}" if bits else "")


def _canon(colour: dict[int, int], adj: dict[int, set[int]], types: list[str]) -> str:
    """Canonical certificate by refinement plus individualisation of one cell."""
    colour = _refine(colour, adj)
    cells: dict[int, list[int]] = {}
    for v, c in colour.items():
        cells.setdefault(c, []).append(v)
    target = next((sorted(vs) for _, vs in sorted(cells.items()) if len(vs) > 1), None)
    if target is None:
        return _code(colour, adj, types)
    best = None
    for v in target:
        split = {u: colour[u] * 2 + (0 if u == v else 1) for u in colour}
        code = _canon(split, adj, types)
        if best is None or code < best:
            best = code
    return best


def certificate(types: list[str], adj: dict[int, set[int]]) -> str:
    """Isomorphism-invariant code of a molecule: atom types plus its bond graph."""
    colour = {v: ATOM_ORDER.index(types[v]) for v in adj}
    return _canon(colour, adj, types)


def formula(types) -> str:
    count = Counter(types)
    return "".join(f"{a}{count[a] if count[a] > 1 else ''}" for a in ATOM_ORDER if count[a])


class _Registry:
    """Molecule certificate -> species id, structure and atom/bond content."""

    def __init__(self) -> None:
        self.id: dict[str, str] = {}
        self.atoms: dict[str, dict[str, int]] = {}
        self.bonds: dict[str, dict[str, int]] = {}
        self.structure: dict[str, str] = {}
        self._taken: dict[str, str] = {}

    def species(self, types: list[str], adj: dict[int, set[int]]) -> str:
        cert = certificate(types, adj)
        if cert in self.id:
            return self.id[cert]
        name = formula(types)
        if len(types) > 2:
            name += "_" + hashlib.blake2s(cert.encode(), digest_size=2).hexdigest()
        if self._taken.setdefault(name, cert) != cert:
            raise ValueError(
                f"species id {name!r} would name two different molecules "
                f"({self._taken[name]!r} and {cert!r})"
            )
        bonds = Counter(
            bond_key(types[v], types[u]) for v in adj for u in adj[v] if u > v
        )
        self.id[cert] = name
        self.atoms[name] = dict(Counter(types))
        self.bonds[name] = dict(bonds)
        self.structure[name] = cert
        return name


# ---------------------------------------------------------------------------
# The reactor
# ---------------------------------------------------------------------------
class _World:
    def __init__(self, p, rng):
        self.p, self.rng = p, rng
        self.width, self.height = p.width, p.height
        self.type: list[str] = []
        self.pos: list[tuple[int, int]] = []
        self.bond: list[set[int]] = []
        self.static: list[set[int]] = []       # bonds that never break (anchors)
        self.occ: dict[tuple[int, int], int] = {}
        self.registry = _Registry()
        self.prob = {"low": p.p_low, "moderate": p.p_moderate, "high": p.p_high}
        self.energy = {"low": p.energy_low, "high": p.energy_high}
        self.fired: dict[tuple, list] = {}
        self.ledger = Counter()
        self.history: dict[str, list[int]] = {
            "sugar_bonds": [], "biomass_bonds": [], "inorganic_bonds": [],
            "free_atoms": [], "molecules": [],
        }
        self.events = Counter()                # (action, bond, catalyst) -> count

    # -- construction -------------------------------------------------------
    def add(self, atom: str, cell: tuple[int, int]) -> int:
        if cell in self.occ:
            raise ValueError(f"cell {cell} is already occupied")
        i = len(self.type)
        self.type.append(atom)
        self.pos.append(cell)
        self.bond.append(set())
        self.static.append(set())
        self.occ[cell] = i
        return i

    def anchor(self, i: int, j: int) -> None:
        """A bond that is part of an organism's body and never reacts."""
        self.bond[i].add(j)
        self.bond[j].add(i)
        self.static[i].add(j)
        self.static[j].add(i)

    def join(self, i: int, j: int) -> None:
        self.bond[i].add(j)
        self.bond[j].add(i)

    def split(self, i: int, j: int) -> None:
        self.bond[i].discard(j)
        self.bond[j].discard(i)

    # -- geometry -----------------------------------------------------------
    def neighbours(self, cell: tuple[int, int]):
        x, y = cell
        for dx, dy in DIRECTIONS:
            yield ((x + dx) % self.width, (y + dy) % self.height)

    def molecule(self, i: int) -> frozenset[int]:
        seen, stack = {i}, [i]
        while stack:
            v = stack.pop()
            for u in self.bond[v]:
                if u not in seen:
                    seen.add(u)
                    stack.append(u)
        return frozenset(seen)

    def molecules(self) -> list[frozenset[int]]:
        out, seen = [], set()
        for i in range(len(self.type)):
            if i not in seen:
                mol = self.molecule(i)
                seen |= mol
                out.append(mol)
        return out

    def clusters(self) -> dict[int, int]:
        """Atom -> id of its contact cluster (atoms touching directly or indirectly).

        Section 2.1: released energy reaches "any continuous atomic structure that
        contacts directly or indirectly (through intermediate neighbours) the
        reaction site", so the pool is spatial, not per molecule.
        """
        label: dict[int, int] = {}
        for start in range(len(self.type)):
            if start in label:
                continue
            cid = len(set(label.values()))
            stack = [start]
            label[start] = cid
            while stack:
                v = stack.pop()
                for cell in self.neighbours(self.pos[v]):
                    u = self.occ.get(cell)
                    if u is not None and u not in label:
                        label[u] = cid
                        stack.append(u)
        return label

    def species(self, atoms: frozenset[int]) -> str:
        order = sorted(atoms)
        index = {v: k for k, v in enumerate(order)}
        types = [self.type[v] for v in order]
        adj = {index[v]: {index[u] for u in self.bond[v]} for v in order}
        return self.registry.species(types, adj)

    # -- one time step ------------------------------------------------------
    def move(self) -> None:
        """Molecules drift one square; bonded atoms move identically, no collisions."""
        mols = self.molecules()
        for k in self.rng.permutation(len(mols)):
            if self.rng.random() >= self.p.move_probability:
                continue
            atoms = mols[k]
            dx, dy = DIRECTIONS[int(self.rng.integers(4))]
            target = {
                i: ((self.pos[i][0] + dx) % self.width, (self.pos[i][1] + dy) % self.height)
                for i in atoms
            }
            if any(self.occ.get(c, i) not in atoms for i, c in target.items()):
                continue
            for i in atoms:
                del self.occ[self.pos[i]]
            for i, c in target.items():
                self.pos[i] = c
                self.occ[c] = i

    def catalyst_of(self, u: int, v: int, key: str, action: str) -> str | None:
        """Section 2.1: a catalyst is present when an atom of the reaction neighbours it."""
        allowed = TABLE[key][MAKE_CAT if action == "make" else BREAK_CAT]
        if not allowed:
            return None
        for site in (u, v):
            for cell in self.neighbours(self.pos[site]):
                w = self.occ.get(cell)
                if w is not None and w not in (u, v) and self.type[w] in allowed:
                    return self.type[w]
        return None

    def candidates(self) -> list[tuple[int, int]]:
        pairs = []
        for cell, u in self.occ.items():
            for other in self.neighbours(cell):
                v = self.occ.get(other)
                if v is not None and u < v and bond_key(self.type[u], self.type[v]) in TABLE:
                    pairs.append((u, v))
        return pairs

    def react(self, t: int) -> None:
        label = self.clusters()
        size = Counter(label.values())
        sun = max(0.0, math.sin(2.0 * math.pi * t / self.p.light_period))
        light = {c: int(self.p.light_amplitude * sun * n) for c, n in size.items()}
        self.ledger["light_incident"] += sum(light.values())
        pool = dict.fromkeys(size, 0)

        pairs = self.candidates()
        draws = self.rng.random(len(pairs))
        for k in self.rng.permutation(len(pairs)):
            u, v = pairs[k]
            key = bond_key(self.type[u], self.type[v])
            bonded = v in self.bond[u]
            if bonded and v in self.static[u]:
                continue
            if not bonded and (len(self.bond[u]) >= VALENCE[self.type[u]]
                               or len(self.bond[v]) >= VALENCE[self.type[v]]):
                continue
            action = "break" if bonded else "make"
            catalyst = self.catalyst_of(u, v, key, action)
            level = TABLE[key][MAKE_CAT if action == "make" else BREAK_CAT]
            rate = level[catalyst] if catalyst else TABLE[key][MAKE if action == "make" else BREAK]
            if draws[k] >= self.prob[rate]:
                continue

            sign, magnitude = TABLE[key][ENERGY]
            released = sign * self.energy[magnitude] * (1 if action == "make" else -1)
            cid = label[u]
            if released < 0:
                need = -released
                from_light = min(light[cid], need) if catalyst == "K" else 0
                if from_light + pool[cid] < need:
                    continue                     # cannot proceed (section 2.2)
                light[cid] -= from_light
                pool[cid] -= need - from_light
                self.ledger["light_spent"] += from_light
                self.ledger["consumed"] += need
            else:
                pool[cid] += released
                self.ledger["released"] += released
            self.record(u, v, action, key, catalyst, released)

        self.ledger["dissipated"] += sum(pool.values())
        self.ledger["light_lost"] += sum(light.values())

    def record(self, u: int, v: int, action: str, key: str,
               catalyst: str | None, released: int) -> None:
        mol_u = self.molecule(u)
        before = [mol_u] if v in mol_u else [mol_u, self.molecule(v)]
        left = [self.species(m) for m in before]
        catalysing = None
        if catalyst is not None:
            cat_atoms = next(
                self.molecule(w)
                for site in (u, v) for cell in self.neighbours(self.pos[site])
                for w in [self.occ.get(cell)]
                if w is not None and w not in (u, v) and self.type[w] == catalyst
            )
            if not any(cat_atoms == m for m in before):
                catalysing = self.species(cat_atoms)

        (self.join if action == "make" else self.split)(u, v)
        mol_u = self.molecule(u)
        after = [mol_u] if v in mol_u else [mol_u, self.molecule(v)]
        right = [self.species(m) for m in after]
        if catalysing is not None:
            left.append(catalysing)
            right.append(catalysing)
        reactants, products = Counter(left), Counter(right)
        rkey = (frozenset(reactants.items()), frozenset(products.items()))
        entry = self.fired.setdefault(
            rkey, [dict(reactants), dict(products), 0,
                   {"action": action, "bond": key, "energy_released": released,
                    "catalysts": Counter()}],
        )
        entry[2] += 1
        if catalyst:
            entry[3]["catalysts"][catalyst] += 1
        self.events[(action, key, catalyst)] += 1

    def sample(self) -> None:
        bonds = Counter()
        for i, partners in enumerate(self.bond):
            for j in partners:
                if j > i:
                    bonds[bond_key(self.type[i], self.type[j])] += 1
        self.history["sugar_bonds"].append(bonds["A-B"])
        self.history["biomass_bonds"].append(bonds["C-C"])
        self.history["inorganic_bonds"].append(bonds["A-O"] + bonds["B-O"])
        self.history["free_atoms"].append(sum(1 for b in self.bond if not b))
        self.history["molecules"].append(len(self.molecules()))

    def run(self) -> None:
        self.initial = Counter(self.species(m) for m in self.molecules())
        self.sample()
        for t in range(self.p.steps):
            self.move()
            self.react(t)
            self.sample()
        self.final = Counter(self.species(m) for m in self.molecules())


# ---------------------------------------------------------------------------
# Setting up a world
# ---------------------------------------------------------------------------
def _counts(name: str, value, allowed) -> dict[str, int]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a dict of atom -> count, got {value!r}")
    out = {}
    for atom, n in value.items():
        if atom not in allowed:
            raise ValueError(f"{name}: unknown atom {atom!r}; choose from {sorted(allowed)}")
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise ValueError(f"{name}[{atom!r}] must be a non-negative integer, got {n!r}")
        out[atom] = n
    return out


def _photoautotroph(world: _World, x0: int, y0: int) -> None:
    """Figure 5: a box of C with chlorophyll on one inner wall and the sugar
    enzyme on the opposite one, enclosing a vacuole."""
    wall = {}
    for x in range(x0, x0 + 5):
        for y in range(y0, y0 + 5):
            if x in (x0, x0 + 4) or y in (y0, y0 + 4):
                wall[(x, y)] = world.add("C", ((x) % world.width, y % world.height))
    for (x, y), i in wall.items():
        for cell in ((x + 1, y), (x, y + 1)):
            if cell in wall:
                world.anchor(i, wall[cell])
    chl = world.add("K", ((x0 + 1) % world.width, (y0 + 2) % world.height))
    enz = world.add("EAB", ((x0 + 3) % world.width, (y0 + 2) % world.height))
    world.anchor(chl, wall[(x0, y0 + 2)])
    world.anchor(enz, wall[(x0 + 4, y0 + 2)])
    for atom, dx in (("A", 1), ("O", 2)):
        world.add(atom, ((x0 + dx) % world.width, (y0 + 1) % world.height))
    world.join(world.occ[((x0 + 1) % world.width, (y0 + 1) % world.height)],
               world.occ[((x0 + 2) % world.width, (y0 + 1) % world.height)])
    for atom, dx in (("B", 1), ("O", 2)):
        world.add(atom, ((x0 + dx) % world.width, (y0 + 3) % world.height))
    world.join(world.occ[((x0 + 1) % world.width, (y0 + 3) % world.height)],
               world.occ[((x0 + 2) % world.width, (y0 + 3) % world.height)])


def _decomposer(world: _World, x0: int, y0: int) -> None:
    """Figure 6: a C-C body carrying the C-C enzyme held off its own bonds by an O."""
    chain = [world.add("C", ((x0 + k) % world.width, y0 % world.height)) for k in range(4)]
    for i, j in zip(chain, chain[1:]):
        world.anchor(i, j)
    spacer = world.add("O", ((x0 + 1) % world.width, (y0 + 1) % world.height))
    enzyme = world.add("ECC", ((x0 + 1) % world.width, (y0 + 2) % world.height))
    world.anchor(chain[1], spacer)
    world.anchor(spacer, enzyme)


FOOTPRINTS = {"photoautotroph": (5, 5), "decomposer": (4, 3)}
BUILD = {"photoautotroph": _photoautotroph, "decomposer": _decomposer}


def _place(world: _World, kind: str, rng) -> None:
    w, h = FOOTPRINTS[kind]
    for _ in range(400):
        x0 = int(rng.integers(world.width))
        y0 = int(rng.integers(world.height))
        cells = [((x0 + dx) % world.width, (y0 + dy) % world.height)
                 for dx in range(-1, w + 1) for dy in range(-1, h + 1)]
        if all(c not in world.occ for c in cells):
            BUILD[kind](world, x0, y0)
            return
    raise ValueError(f"no free room for the {kind} structure; enlarge width/height")


def _check(p) -> None:
    if p.energy_high <= p.energy_low:
        raise ValueError(
            f"energy_high ({p.energy_high}) must exceed energy_low ({p.energy_low}): "
            "table 1 stores energy in the A-B sugar bond"
        )
    for name in ("p_low", "p_moderate", "p_high"):
        value = getattr(p, name)
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be a probability in [0, 1], got {value}")
    if not p.p_low <= p.p_moderate <= p.p_high:
        raise ValueError(
            f"table 1 is ordered: p_low ({p.p_low}) <= p_moderate ({p.p_moderate}) "
            f"<= p_high ({p.p_high})"
        )


def generate(p, rng):
    _check(p)
    atoms = _counts("atoms", p.atoms, BUILDING_BLOCKS)
    catalysts = _counts("catalysts", p.catalysts, CATALYSTS)
    world = _World(p, rng)

    kinds = ([] if p.structures == "none"
             else list(FOOTPRINTS) if p.structures == "all" else [p.structures])
    for kind in kinds:
        _place(world, kind, rng)

    free = [(x, y) for x in range(p.width) for y in range(p.height) if (x, y) not in world.occ]
    wanted = {**atoms, **catalysts}
    total = sum(wanted.values())
    if total > len(free):
        raise ValueError(
            f"{total} atoms do not fit in the {p.width}x{p.height} grid "
            f"({len(free)} cells are free); lower atoms/catalysts or enlarge the grid"
        )
    chosen = rng.choice(len(free), size=total, replace=False)
    order = [a for a in ATOM_ORDER for _ in range(wanted.get(a, 0))]
    for atom, k in zip(order, chosen):
        world.add(atom, free[int(k)])

    world.run()
    return _network(world, p)


def _network(world: _World, p) -> Network:
    reg = world.registry
    names = sorted(reg.atoms)
    reactions, events = [], []
    for reactants, products, count, detail in world.fired.values():
        reactions.append(Reaction(reactants, products, None, count))
        events.append({
            "action": detail["action"], "bond": detail["bond"],
            "energy_released": detail["energy_released"],
            "catalysts": dict(sorted(detail["catalysts"].items())),
            "count": count,
        })

    energy = {name: world.energy[level] * sign for name, (sign, level)
              in ((b, TABLE[b][ENERGY]) for b in TABLE)}
    ledger = {k: int(v) for k, v in sorted(world.ledger.items())}
    ledger["balanced"] = (ledger.get("released", 0) + ledger.get("light_spent", 0)
                          == ledger.get("consumed", 0) + ledger.get("dissipated", 0))
    by_event = {f"{a}:{b}:{c or 'none'}": n for (a, b, c), n in sorted(
        world.events.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2] or ""))}

    return Network(
        species=[Species(n, structure=reg.structure[n]) for n in names],
        reactions=reactions,
        status="observed",
        initial_state={n: float(c) for n, c in sorted(world.initial.items())},
        extras={
            "energies": {
                "units": "integer energy units; energy_low and energy_high scale table 1",
                "convention": "a bond's energy is released when the bond is made and must be "
                              "supplied to break it (table 1 caption), so negative bonds store energy",
                "bond_energy": energy,
                "light": {
                    "amplitude": p.light_amplitude, "period": p.light_period,
                    "law": "floor(amplitude * max(0, sin(2 pi t / period)) * atoms in the cluster)",
                    "usable_by": "chlorophyll-catalysed reactions only (section 2.2)",
                },
                "species_bonds": {n: reg.bonds[n] for n in names},
                "ledger": ledger,
            },
            "conservation": [
                {"name": f"atom {a}", "vector": {n: reg.atoms[n][a] for n in names if a in reg.atoms[n]}}
                for a in ATOM_ORDER if any(a in reg.atoms[n] for n in names)
            ],
            "space": {
                "dimensions": 2, "lattice": "torus", "width": p.width, "height": p.height,
                "neighbourhood": "von-neumann", "cell": "one square atom, bonds are shared edges",
                "motion": "each molecule moves one square in a random direction with probability "
                          f"{p.move_probability}; bonded atoms move identically",
                "collisions": "not permitted: a move is rejected if any target cell is occupied",
            },
            "chemistry": {
                "valence": VALENCE,
                "shells": {a: list(s) for a, s in SHELLS.items()},
                "table": {b: {"make": r[MAKE], "break": r[BREAK],
                              "make_catalyst": r[MAKE_CAT], "break_catalyst": r[BREAK_CAT],
                              "energy": ("-" if r[ENERGY][0] < 0 else "+") + " " + r[ENERGY][1]}
                          for b, r in TABLE.items()},
                "probabilities": {"low": p.p_low, "moderate": p.p_moderate, "high": p.p_high},
                "catalysts": CATALYSTS,
                "static_bonds": "bonds inside a seeded organism body are anchors and never react",
            },
            "events": events,
            "species_encoding": "id = atom formula for one- and two-atom molecules, else the "
                                "formula and a hash of the certificate; structure = certificate, "
                                "the canonical (isomorphism-invariant) bond graph as atom types "
                                "in canonical order followed by the upper adjacency bits",
            "final_state": {n: int(c) for n, c in sorted(world.final.items())},
            "analysis": {
                "steps": p.steps,
                "event_counts": by_event,
                "history": world.history,
                "trophic": {
                    "photosynthesis": "AO + BO -(K, sunlight)-> AB + 2 O (section 3.1.1)",
                    "respiration": "O + AB -(EAB)-> A + BO + energy (section 3.1.2)",
                    "biosynthesis": "C + C -(energy)-> C2 (section 3.1.3)",
                    "sugar_made": world.events[("make", "A-B", "K")],
                    "sugar_respired": world.events[("break", "A-B", "EAB")],
                    "biomass_built": sum(n for (a, b, _), n in world.events.items()
                                         if a == "make" and b == "C-C"),
                    "biomass_decomposed": world.events[("break", "C-C", "ECC")],
                    "inorganic_split": (world.events[("break", "A-O", "K")]
                                        + world.events[("break", "B-O", "K")]
                                        + world.events[("break", "A-O", "EO")]
                                        + world.events[("break", "B-O", "EO")]),
                },
            },
        },
    )
