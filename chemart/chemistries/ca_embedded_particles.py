"""Embedded particles in cellular automata (book 10.7.2). Catalog id: ca-embedded-particles.

The book's cleanest "a CA *is* a chemistry" reading: in a one-dimensional
binary CA the space-time behaviour organises into *regular domains*; the walls
between domains are spatially localised, temporally periodic structures -
*particles* - and a collision of two particles is literally a reaction
(book 10.7.2, figure 10.14, after Hordijk, Crutchfield & Mitchell).

So: **domains are the solvent, particles are the molecules, and the published
particle interaction table is the reaction network.**

    alpha -> gamma + mu          (an unstable wall decays)
    beta + gamma -> eta          (two signals react)
    eta + mu -> ∅                (annihilation, leaving domain Lambda1)

Everything here is run honestly. The CA is the published look-up table, applied
synchronously on a periodic lattice. Domains are detected by matching every
site against the phases of each domain's spatial pattern over a window of the
neighbourhood's width; sites that match nothing are walls; each maximal wall is
a particle, named by the pair of domains it separates; particles are followed
from step to step and every appearance/disappearance event is recorded as a
reaction with its firing count.

Two faces. ``generate`` returns the particle catalog's own interaction table,
exactly as printed (status ``complete``); nothing is run. ``evolve`` runs the
automaton and yields a frame per CA iteration: the particles present, the
interactions recorded in that iteration, and the density of 1s; it returns the
interactions that fired with their counts (status ``observed``).

Rules (`rule`), with the lookup tables embedded verbatim from their sources:

    phi-par-a   the density-classification CA of Crutchfield, Mitchell & Das,
                with the particle catalog of their Table 3
    phi-par-b   the second density CA, particle catalog of their Table 4
    gkl         the hand-designed Gacs-Kurdyumov-Levin rule; the same three
                domains, but its particle catalog is not published, so its
                walls are named by the same domain-pair convention

extras["analysis"] holds the catalog (domains, particles, velocities, the
published interaction table) and, for a run, the measured particle velocities,
the observed interaction counts and the classification outcome;
extras["space"] holds the rule and, for a run, the lattice and the filtered
space-time diagram.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from chemart.network import Network, Reaction, Species
from chemart.trajectory import Frame, ticks

#: CA radius of every rule here: neighbourhood = 2r + 1 = 7 cells.
RADIUS = 3

#: Greek particle names, in the order the sources list them.
ORDER = ("alpha", "beta", "gamma", "delta", "eta", "mu")


# ---------------------------------------------------------------------------
# The published rules and particle catalogs
# ---------------------------------------------------------------------------
#: Lookup tables as hexadecimal, verbatim from the sources (see the catalog
#: entry). Expanding each hex digit to binary, left to right, gives the 128
#: output bits in lexicographic order of neighbourhood, the leftmost bit being
#: the output for neighbourhood 0000000.
RULES: dict[str, dict] = {
    "phi-par-a": {
        "task": "density-classification",
        "runnable": True,
        "hex": "0504058605000F77037755877BFFB77F",
        "source": "Crutchfield, Mitchell & Das (1998), Table 1 (phi_par^a)",
        "performance": {"149": 0.775, "599": 0.740, "999": 0.728},
        "domains": {"L0": "0", "L1": "1", "L2": "01"},
        "catalog": "Crutchfield, Mitchell & Das (1998), Table 3",
        # name: (left domain, right domain, velocity or None for unstable)
        "particles": {
            "alpha": ("L0", "L1", None),
            "beta": ("L1", "L0", 0.0),
            "gamma": ("L0", "L2", -1.0),
            "delta": ("L2", "L0", -3.0),
            "eta": ("L1", "L2", 3.0),
            "mu": ("L2", "L1", 1.0),
        },
        # (type, reactants, products, domain left behind by an annihilation)
        "interactions": [
            ("decay", ["alpha"], ["gamma", "mu"], None),
            ("react", ["beta", "gamma"], ["eta"], None),
            ("react", ["mu", "beta"], ["delta"], None),
            ("react", ["eta", "delta"], ["beta"], None),
            ("annihilate", ["eta", "mu"], [], "L1"),
            ("annihilate", ["gamma", "delta"], [], "L0"),
        ],
    },
    "phi-par-b": {
        # The lookup table as printed does not reproduce this rule's published
        # performance (see the catalog entry's `decisions`), so only its
        # published particle catalog is offered, not a run of the automaton.
        "task": "density-classification",
        "runnable": False,
        "hex": "00240066A0A0224676EFEFFFFBFFAAFE",
        "source": "Crutchfield, Mitchell & Das (1998), Table 1 (phi_par^b)",
        "performance": {"149": 0.766, "599": 0.687, "999": 0.641},
        "domains": {"L0": "0", "L1": "1", "L2": "011"},
        "catalog": "Crutchfield, Mitchell & Das (1998), Table 4",
        "particles": {
            "alpha": ("L1", "L0", 0.0),
            "beta": ("L0", "L1", 1.0),
            "gamma": ("L1", "L2", 0.0),
            "delta": ("L2", "L1", -3.0),
            "eta": ("L0", "L2", 3.0),
            "mu": ("L2", "L0", 1.5),
        },
        "interactions": [
            ("decay", ["alpha"], ["gamma", "mu"], None),
            ("react", ["beta", "gamma"], ["eta"], None),
            ("react", ["mu", "beta"], ["delta"], None),
            ("react", ["eta", "delta"], ["beta"], None),
            ("annihilate", ["eta", "mu"], [], "L0"),
            ("annihilate", ["gamma", "delta"], [], "L1"),
        ],
    },
    "gkl": {
        "task": "density-classification",
        "runnable": True,
        "hex": "005F005F005F005F005FFF5F005FFF5F",
        "source": "Crutchfield & Mitchell (1995), Table 1 (GKL); expands exactly "
                  "to the majority definition of Gacs, Kurdyumov & Levin",
        "performance": {"149": 0.816, "599": 0.766, "999": 0.757},
        "domains": {"L0": "0", "L1": "1", "L2": "01"},
        # No particle catalog is published for GKL in these sources, and its
        # walls do not behave like the ones phi-par-a names with the same
        # domain pair (GKL's Lambda0Lambda2 wall moves +3 where phi-par-a's
        # moves -1), so no name and no velocity is borrowed: every wall is
        # reported by the domains it separates and its velocity is measured.
        "catalog": None,
        "particles": {},
        "interactions": [],
    },
}

#: The particles a wall can be, keyed by rule and by (left, right) domain.
def _by_wall(rule: str) -> dict[tuple[str, str], str]:
    return {(l, r): name for name, (l, r, _v) in RULES[rule]["particles"].items()}


# ---------------------------------------------------------------------------
# The automaton
# ---------------------------------------------------------------------------
def lookup(hex_code: str) -> np.ndarray:
    """The 128 output bits of an r = 3 rule, indexed by the neighbourhood value.

    Neighbourhood (s[i-3] ... s[i+3]) read as a binary number with s[i-3] the
    most significant bit, which is the "lexicographic order of neighbourhood,
    leftmost bit = neighbourhood 0000000" convention of the sources.
    """
    digits = "".join(hex_code.split())
    if len(digits) != 32:
        raise ValueError(f"an r = 3 lookup table needs 32 hex digits, got {len(digits)}")
    bits = "".join(bin(int(d, 16))[2:].zfill(4) for d in digits)
    return np.array([int(b) for b in bits], dtype=np.int8)


def gkl_lookup() -> np.ndarray:
    """The GKL rule from its definition (Gacs, Kurdyumov & Levin).

    ``s[i] = 0 -> majority(s[i], s[i-1], s[i-3]);
    s[i] = 1 -> majority(s[i], s[i+1], s[i+3])``
    (Mitchell, Hraber & Crutchfield 1993, Sec. 3).
    """
    table = np.zeros(128, dtype=np.int8)
    for index in range(128):
        cell = [(index >> (6 - k)) & 1 for k in range(7)]   # cell[k] = s[i-3+k]
        centre = cell[3]
        vote = (centre, cell[2], cell[0]) if centre == 0 else (centre, cell[4], cell[6])
        table[index] = int(sum(vote) >= 2)
    return table


def step(lattice: np.ndarray, table: np.ndarray) -> np.ndarray:
    """One synchronous update of a periodic lattice."""
    index = np.zeros(lattice.shape, dtype=np.int64)
    for k in range(-RADIUS, RADIUS + 1):
        index = index * 2 + np.roll(lattice, -k)
    return table[index]


def initial_lattice(n: int, density: float, rng) -> np.ndarray:
    """An IC with exactly round(density * n) ones, in random positions."""
    ones = int(round(density * n))
    lattice = np.zeros(n, dtype=np.int8)
    lattice[:ones] = 1
    return lattice[rng.permutation(n)]


def classify(lattice: np.ndarray, density: float) -> str:
    """'correct', 'incorrect' or 'undecided' for the density task."""
    if lattice.min() != lattice.max():
        return "undecided"
    want = 1 if density > 0.5 else 0
    return "correct" if int(lattice[0]) == want else "incorrect"


def performance(rule: str, n: int = 149, trials: int = 100, rng=None,
                steps: int | None = None) -> float:
    """P_{N,I}: the fraction of unbiased random ICs classified correctly.

    The task of the sources: rho_c = 1/2, T_max = 2N, no partial credit for
    configurations that have not reached an all-0s or all-1s fixed point.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    table = lookup(RULES[rule]["hex"])
    steps = 2 * n if steps is None else steps
    correct = 0
    for _ in range(trials):
        lattice = rng.integers(0, 2, n).astype(np.int8)
        want = 1 if lattice.mean() > 0.5 else 0
        for _t in range(steps):
            nxt = step(lattice, table)
            if np.array_equal(nxt, lattice):
                break
            lattice = nxt
        correct += int(lattice.min() == lattice.max() == want)
    return correct / trials


# ---------------------------------------------------------------------------
# The domain filter
# ---------------------------------------------------------------------------
def domain_labels(lattice: np.ndarray, domains: dict[str, str], window: int) -> np.ndarray:
    """Index of the domain each site participates in, or -1 for a wall.

    A site is in domain Lambda when the whole window of 2*window + 1 sites
    around it agrees with one spatial phase of Lambda's pattern - the same
    criterion as the published domain-recognising transducer, which leaves the
    sites it cannot assign as wall sites.
    """
    n = lattice.size
    names = list(domains)
    out = np.full(n, -1, dtype=np.int8)
    position = np.arange(n)
    for d, name in enumerate(names):
        pattern = np.array([int(c) for c in domains[name]], dtype=np.int8)
        period = pattern.size
        covered = np.zeros(n, dtype=bool)
        for phase in range(period):
            agrees = lattice == pattern[(position - phase) % period]
            # every site of the window agrees: circular all() over 2*window + 1
            runs = agrees.copy()
            for shift in range(-window, window + 1):
                if shift:
                    runs &= np.roll(agrees, -shift)
            covered |= runs
        out[(out < 0) & covered] = d
    return out


def walls(labels: np.ndarray, names: list[str]) -> list[tuple[float, str, str, int]]:
    """Every wall as (centre position, left domain, right domain, width).

    A wall is a maximal circular run of unassigned sites; a direct change of
    domain with no unassigned site between is a wall of width 0.
    """
    n = labels.size
    if (labels < 0).all():
        return []
    start = int(np.argmax(labels >= 0))          # a site that is in a domain
    order = (np.arange(n) + start) % n
    seq = labels[order]
    out = []
    i = 0
    while i < n:
        if seq[i] >= 0:
            j = (i + 1) % n
            if seq[j] >= 0 and seq[j] != seq[i]:
                out.append((float(order[i]) + 0.5, names[seq[i]], names[seq[j]], 0))
            i += 1
            continue
        j = i
        while j < n and seq[j] < 0:
            j += 1
        left = names[seq[i - 1]]
        right = names[seq[j % n]]
        centre = (float(order[i]) + (j - i - 1) / 2.0) % n
        out.append((centre, left, right, j - i))
        i = j
    return out


def name_wall(left: str, right: str, by_wall: dict[tuple[str, str], str]) -> str:
    """The catalog name of a wall, or a `w<left><right>` id for an unlisted one."""
    return by_wall.get((left, right), f"w{left}{right}")


def two_domain_lattice(spec: dict, left: str, right: str, n: int) -> np.ndarray:
    """A periodic lattice made of domain `left` on one half and `right` on the other.

    It carries exactly two walls: the `left|right` one at the middle and the
    `right|left` one at the wrap-around.
    """
    pattern = {name: [int(c) for c in word] for name, word in spec["domains"].items()}
    half = n // 2
    cells = [pattern[left][i % len(pattern[left])] for i in range(half)]
    cells += [pattern[right][i % len(pattern[right])] for i in range(n - half)]
    return np.array(cells, dtype=np.int8)


def probe(rule: str, particle: str, n: int = 199, steps: int = 20,
          window: int = RADIUS) -> dict:
    """Measure one particle by building the two domains it separates.

    Returns its measured velocity (displacement per step of the wall at the
    middle of the lattice) or, when the wall does not survive, the walls that
    replaced it - which is how an unstable particle's decay shows up.
    """
    spec = RULES[rule]
    left, right, _v = spec["particles"][particle]
    by_wall = _by_wall(rule)
    table = lookup(spec["hex"])
    lattice = two_domain_lattice(spec, left, right, n)
    names = list(spec["domains"])
    here = float(n // 2) - 0.5
    moved, lived, products = 0.0, 0, []
    for t in range(1, steps + 1):
        lattice = step(lattice, table)
        found = [(centre, name_wall(l, r, by_wall))
                 for centre, l, r, _w in walls(domain_labels(lattice, spec["domains"], window), names)]
        same = [(abs(_gap(here, centre, n)), centre) for centre, name in found
                if name == particle and abs(_gap(here, centre, n)) <= MAX_MOVE]
        if not same:
            products = sorted(name for centre, name in found
                              if abs(_gap(here, centre, n)) <= 2 * REACH)
            break
        distance, centre = min(same)
        moved += _gap(here, centre, n)
        here = centre
        lived = t
    return {"particle": particle, "wall": f"Λ{left[1:]}Λ{right[1:]}",
            "velocity": round(moved / lived, 3) if lived else None,
            "steps_followed": lived, "decays_into": products}


# ---------------------------------------------------------------------------
# Following the particles and recording their interactions
# ---------------------------------------------------------------------------
def _gap(a: float, b: float, n: int) -> float:
    """Signed circular displacement from a to b, in (-n/2, n/2]."""
    return (b - a + n / 2) % n - n / 2


#: A particle interaction is not instantaneous: for a few steps around a
#: collision the lattice cannot be decomposed into domains and particles
#: (the fourth simplifying assumption of the embedded-particle model,
#: Hordijk, Crutchfield & Mitchell 1998). What vanishes and what appears
#: within WIDEST_WALL cells and SETTLE steps is collected into one event: the
#: radius has to be the width of a wall, not of the neighbourhood, because the
#: centres of two colliding wedges are that far apart when they meet.
SETTLE = 4
REACH = 2 * RADIUS + 1

#: The fastest published particle moves 3 cells per step, so a particle
#: followed from one step to the next may not move further than this.
MAX_MOVE = 4

#: A wall wider than this is not a particle but a region the description does
#: not cover yet. Particle wedges are wider than the neighbourhood - phi-par-a's
#: beta wedge is a steady 8 cells - so the bound is two neighbourhoods.
WIDEST_WALL = 2 * REACH


class _Event:
    """One collision in progress: what vanished and what appeared near a place."""

    def __init__(self, pos: float, t: int):
        self.pos = pos
        self.reactants: Counter = Counter()
        self.products: Counter = Counter()
        self.first = t
        self.last = t


class Run:
    """One run of the CA, read as a chemistry of embedded particles."""

    def __init__(self, rule: str, window: int):
        self.rule = rule
        self.spec = RULES[rule]
        self.by_wall = _by_wall(rule)
        self.window = window
        self.domains = list(self.spec["domains"])
        self.events: Counter = Counter()
        self.recent: Counter = Counter()      # events recorded since the last flush()
        self.seen: Counter = Counter()
        self.condensation: int | None = None
        self.filtered: list[str] = []
        self.open: list[_Event] = []
        self.tracks: dict[int, dict] = {}     # one particle followed over time
        self.live: list[int] = []             # track id of each current particle

    # -- one configuration ------------------------------------------------
    def observe(self, lattice: np.ndarray, t: int) -> list[tuple[float, str]]:
        labels = domain_labels(lattice, self.spec["domains"], self.window)
        found = walls(labels, self.domains)
        particles = [(centre, name_wall(left, right, self.by_wall))
                     for centre, left, right, _w in found]
        # Condensation time: every site is in a regular domain or in a narrow
        # wall between two of them. A wall between two patches of the *same*
        # domain is a dislocation - the published catalogs leave those out
        # ("the four dislocations within Lambda2 are not listed", Crutchfield &
        # Mitchell 1995, Table 2) - so they count as condensed too.
        # A rule with no published catalog (GKL) has no names to demand, so for
        # it condensation is the width criterion alone.
        condensed = bool((labels >= 0).any()) and all(
            width <= WIDEST_WALL and (not self.by_wall
                                      or (left, right) in self.by_wall or left == right)
            for _c, left, right, width in found)
        if condensed and self.condensation is None:
            self.condensation = t
        self.seen.update(name for _c, name in particles)
        self.filtered.append("".join("." if v >= 0 else "#" for v in labels))
        return sorted(particles)

    # -- follow the particles from one step to the next --------------------
    def transition(self, before: list[tuple[float, str]],
                   after: list[tuple[float, str]], t: int, n: int) -> None:
        """Match same-named particles by proximity; the rest are interactions.

        The match uses no knowledge of the published velocities - it takes the
        nearest particle of the same type that is within MAX_MOVE cells - so
        the velocity it measures is an independent measurement.
        """
        if len(self.live) != len(before):
            self.live = [self.new_track(name, t - 1) for _pos, name in before]
        matched: set[int] = set()
        taken: dict[int, int] = {}                     # index in after -> in before
        candidates = []
        for i, (pos, name) in enumerate(before):
            for j, (pos2, name2) in enumerate(after):
                if name2 != name:
                    continue
                move = _gap(pos, pos2, n)
                if abs(move) <= MAX_MOVE:
                    candidates.append((abs(move), i, j, move))
        for _cost, i, j, move in sorted(candidates):
            if i in matched or j in taken:
                continue
            matched.add(i)
            taken[j] = i
            track = self.tracks[self.live[i]]
            track["moved"] += move
            track["until"] = t
        live = []
        for j, (_pos, name) in enumerate(after):
            live.append(self.live[taken[j]] if j in taken else self.new_track(name, t))
        self.absorb(t, n,
                    [before[i] for i in range(len(before)) if i not in matched],
                    [after[j] for j in range(len(after)) if j not in taken])
        self.live = live

    def new_track(self, name: str, t: int) -> int:
        tid = len(self.tracks)
        self.tracks[tid] = {"name": name, "since": t, "until": t, "moved": 0.0}
        return tid

    def velocities(self, least: int = 5) -> dict[str, dict]:
        """Measured velocity of each particle type: displacement over lifetime.

        Only tracks that survive `least` steps are used, and the matching that
        built them never looks at the published velocities.
        """
        by_name: dict[str, list[float]] = {}
        for track in self.tracks.values():
            span = track["until"] - track["since"]
            if span >= least:
                by_name.setdefault(track["name"], []).append(track["moved"] / span)
        return {name: {"velocity": round(float(np.mean(v)), 3), "tracks": len(v)}
                for name, v in sorted(by_name.items())}

    # -- interactions ------------------------------------------------------
    def find(self, pos: float, t: int, n: int) -> _Event:
        near = [(abs(_gap(e.pos, pos, n)), k) for k, e in enumerate(self.open)
                if abs(_gap(e.pos, pos, n)) <= WIDEST_WALL]
        if near:
            return self.open[min(near)[1]]
        event = _Event(pos, t)
        self.open.append(event)
        return event

    def absorb(self, t: int, n: int, gone: list[tuple[float, str]],
               born: list[tuple[float, str]]) -> None:
        """Add what vanished and what appeared to the event at that place.

        A particle that vanishes and comes back within the settling window was
        only hidden by the collision it is passing through, so it cancels out
        instead of being read as a reaction.
        """
        self.close(t - SETTLE)
        for pos, name in gone:
            event = self.find(pos, t, n)
            if event.products[name]:
                event.products[name] -= 1
                event.products += Counter()          # drop zero counts
            else:
                event.reactants[name] += 1
            event.last = t
        for pos, name in born:
            event = self.find(pos, t, n)
            if event.reactants[name]:
                event.reactants[name] -= 1
                event.reactants += Counter()
            else:
                event.products[name] += 1
            event.last = t

    def close(self, before: int) -> None:
        """Record every event that has settled, and keep the rest open."""
        still_open = []
        for event in self.open:
            if event.last > before:
                still_open.append(event)
                continue
            self.emit(event)
        self.open = still_open

    def emit(self, event: _Event) -> None:
        if self.condensation is None or event.first < self.condensation:
            return                                    # pre-condensation debris
        left = tuple(sorted(event.reactants.elements()))
        right = tuple(sorted(event.products.elements()))
        if left or right:
            self.events[(left, right)] += 1
            self.recent[(left, right)] += 1

    def flush(self) -> list[list]:
        """The events recorded since the last flush, as a frame's `fired`; then forget them."""
        out = [[list(lhs), list(rhs), n] for (lhs, rhs), n in self.recent.items()]
        self.recent = Counter()
        return out

    def finish(self) -> None:
        self.close(float("inf"))


# ---------------------------------------------------------------------------
# The chemistry
# ---------------------------------------------------------------------------
def _structure(name: str, spec) -> str:
    """A species' structure: the wall it is, and its published velocity if any."""
    if name in spec["particles"]:
        left, right, velocity = spec["particles"][name]
        return (f"Λ{left[1:]}Λ{right[1:]} "
                f"(v={'unstable' if velocity is None else velocity})")
    left, right = name[1:3], name[3:5]                   # wL0L2 -> L0, L2
    kind = "dislocation" if left == right else "wall, no published catalog"
    return f"Λ{left[1:]}Λ{right[1:]} ({kind}; velocity measured, not published)"


def _catalog(spec) -> dict:
    return {
        "domains": {f"Λ{k[1:]}": f"({v})+" for k, v in spec["domains"].items()},
        "particles": {
            name: {"wall": f"Λ{left[1:]}Λ{right[1:]}",
                   "velocity": velocity, "stable": velocity is not None}
            for name, (left, right, velocity) in spec["particles"].items()
        },
        "interactions": [
            {"type": kind, "reaction": " + ".join(lhs) + " -> " + (" + ".join(rhs) or "∅"),
             "leaves": f"Λ{left[1:]}" if left else None}
            for kind, lhs, rhs, left in spec["interactions"]
        ],
        "source": spec["catalog"],
    }


def _analysis(p, spec) -> dict:
    return {
        "rule": p.rule,
        "task": spec["task"],
        "lookup_hex": spec["hex"],
        "published_performance": spec["performance"],
        "catalog": _catalog(spec),
    }


def _space(p, spec) -> dict:
    return {
        "dimensions": 1, "boundary": "periodic", "states": 2,
        "radius": RADIUS, "neighbourhood": 2 * RADIUS + 1, "rule": p.rule,
        "lookup_hex": spec["hex"],
        "lookup_convention": "hex digits left to right give the 128 output bits in "
                             "lexicographic order of neighbourhood; the leftmost bit "
                             "is the output for neighbourhood 0000000",
    }


def generate(p, rng):
    """The rule's published particle interaction table, as a complete network."""
    spec = RULES[p.rule]
    if not spec["interactions"]:
        raise ValueError(
            f"no particle interaction table is published for rule {p.rule!r}; "
            "use chemart.evolve to record the interactions of a run"
        )
    species = [Species(name, structure=_structure(name, spec))
               for name in ORDER if name in spec["particles"]]
    reactions = [Reaction(dict(Counter(lhs)), dict(Counter(rhs)))
                 for _kind, lhs, rhs, _left in spec["interactions"]]
    return Network(species=species, reactions=reactions, status="complete",
                   extras={"space": _space(p, spec), "analysis": _analysis(p, spec)})


def evolve(p, rng):
    """Run the automaton: a frame per CA iteration."""
    spec = RULES[p.rule]
    table = lookup(spec["hex"])
    n = p.lattice
    if n < 4 * RADIUS:
        raise ValueError(f"lattice must hold a few neighbourhoods; got {n} < {4 * RADIUS}")
    if not spec["runnable"]:
        raise ValueError(
            f"the published look-up table of rule {p.rule!r} does not reproduce its "
            f"published behaviour (see the catalog entry's decisions), so only its "
            f"particle catalog is offered: use chemart.generate_network"
        )

    lattice = initial_lattice(n, p.density, rng)
    space = {**_space(p, spec), "shape": [n], "initial": "".join(str(v) for v in lattice)}

    def frame(t, particles):
        return Frame(t=float(t), state={name: float(k) for name, k in
                                        sorted(Counter(name for _c, name in particles).items())},
                     fired=run.flush(), observables={"density": float(lattice.mean())})

    run = Run(p.rule, p.filter_window)
    before = run.observe(lattice, 0)
    at_condensation = before if run.condensation == 0 else None
    yield frame(0, before)
    for t in ticks(p.steps):
        lattice = step(lattice, table)
        after = run.observe(lattice, t)
        run.transition(before, after, t, n)
        before = after
        if run.condensation == t and at_condensation is None:
            at_condensation = after
        if t == p.steps:
            run.finish()                       # the events still open end with the run
        yield frame(t, after)

    start = at_condensation if at_condensation is not None else before
    names = sorted({name for _kind, lhs, rhs, _l in spec["interactions"] for name in lhs + rhs}
                   | set(run.seen) | set(spec["particles"]),
                   key=lambda s: (ORDER.index(s) if s in ORDER else len(ORDER), s))
    species = [Species(name, structure=_structure(name, spec)) for name in names]
    reactions = [Reaction(dict(Counter(lhs)), dict(Counter(rhs)), count=count)
                 for (lhs, rhs), count in sorted(run.events.items(),
                                                 key=lambda kv: (-kv[1], kv[0]))]
    analysis = _analysis(p, spec)
    analysis.update({
        "steps": p.steps,
        "initial_density": float(np.mean([int(c) for c in space["initial"]])),
        "condensation_time": run.condensation,
        "classification": classify(lattice, p.density),
        "particles_seen": dict(run.seen),
        "measured_velocities": run.velocities(),
        "interactions_observed": [
            {"reaction": " + ".join(lhs) + " -> " + (" + ".join(rhs) or "∅"), "count": count}
            for (lhs, rhs), count in sorted(run.events.items(), key=lambda kv: (-kv[1], kv[0]))
        ],
    })
    space["final"] = "".join(str(v) for v in lattice)
    space["filtered"] = run.filtered[:: max(1, len(run.filtered) // 200)]
    space["filtered_key"] = ". = site in a regular domain, # = site in a wall"
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={name: float(count)
                       for name, count in sorted(Counter(n for _c, n in start).items())},
        extras={"space": space, "analysis": analysis},
    )
