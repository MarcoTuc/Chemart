"""Squirm3: Hutton's 2D artificial chemistry of typed, stateful atoms. Catalog id: squirm3.

Book 11.4.1. An atom has a fixed type in {a..f} and a mutable integer state.
Reactions are local rules on two (or three) colliding atoms that change their
states and make or break bonds; a bonded group of atoms is a molecule and moves
as a unit. Rules are *specific* (fixed types) or *general* (type variables
x, y, z: the same letter means the same type, different letters are
unconstrained). Written in the papers' own notation, where juxtaposition means
bonded and " + " means not bonded:

    e8 + e0 -> e4e3     association     (R1 of table 11.5)
    x4y1 -> x2y5        transformation  (R2)
    f4f3 -> f8 + f8     dissociation    (R6)

Species are molecules, identified by a canonical string over their bonded atoms
("e8-a1-b1-f1") that is invariant under the molecule's symmetries. Atoms are
never created or destroyed, so the count of each atom type is an exact
conservation law of every reaction (extras["conservation"]).

Rule sets (parameter `rules`):
  replicator  the eight rules of Hutton (2002) table 1 = book table 11.5, which
              replicate any string of a1..d1 flanked by e8 and f1.
  membrane    the 41 rules of Hutton (2007) for self-reproducing cells, which
              add membrane growth, gene read-out into enzymes and mutation.
  custom      pair rules given as text in `rule_text`.

Enzymatic reactions (2007 R40, book eq. 11.15): a catalyst z(i) whose state i
encodes a whole two-atom reaction,

    i = 2(2(T(T(S(S(S g + h) + j) + k) + x) + y) + b1) + b2 + S

with S the number of states, T the number of types, g,h the reactant states,
j,k the product states, x,y the reactant types and b1,b2 the bond before and
after. i is read off a gene of base atoms in base |bases| (2007 R37).

Faces:
  generate  chemart.expand.expand over molecules: every reaction reachable from
            the seed molecule and one food atom of each type, without a
            spatial constraint (the closure).
  evolve    a 2D lattice (2002 CA physics) or continuous space (2007), with
            optional periodic flooding; yields a frame about every
            steps // 200 steps and returns the reactions that fired, with
            counts.
"""

from __future__ import annotations

import math
import re
from collections import Counter, deque
from dataclasses import dataclass

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame, ticks

#: Atom types of the papers; |T| enters the enzyme encoding.
ATOM_TYPES = "abcdef"
#: Type variables usable in a rule ("the same letter means the same type").
TYPE_VARS = "xyz"

_ATOM = re.compile(r"([a-fxyz])(\d+)")

# --- Hutton (2002) table 1 = book table 11.5 ------------------------------------
#: The self-replicator: any string of a1..d1 with e8 at one end and f1 at the
#: other replicates in a soup of atoms in state 0.
REPLICATOR_RULES = """
R1: e8 + e0 -> e4e3
R2: x4y1 -> x2y5
R3: x5 + x0 -> x7x6
R4: x3 + y6 -> x2y3
R5: x7y3 -> x4y3
R6: f4f3 -> f8 + f8
R7: x2y8 -> x9y1
R8: x9y9 -> x8 + y8
"""

# --- Hutton (2007) table 1, R1-R34 (the pair rules of cell reproduction) --------
MEMBRANE_PAIR_RULES = """
R1: e1a37 -> e5a10
R2: a10 + e6 -> a37e3
R3: e6e3 -> e2e3
R4: x2y1 -> x7y4
R5: x4 + y3 -> x5y7
R6: x5 + x0 -> x6x6
R7: x6 + y7 -> x3y4
R8: x6y4 -> x1 + y2
R9: x7y1 -> x2y2
R10: f2a37 -> f9a11
R11: a11 + f3 -> a11f9
R12: x2y8 -> x9y1
R13: x9y9 -> x8 + y8
R14: a11a36 -> a11a12
R15: f1 + a12 -> f13a37
R16: x13y1 -> x14y15
R17: a11 + x15 -> a11x16
R18: x14y16 -> x27y16
R19: x27a11 -> x17 + a11
R20: x17y16 -> x17y13
R21: x13e8 -> x14e15
R22: e13a37 -> e18a19
R23: e13a19 -> e18a20
R24: a20 + a11 -> a21a22
R25: e18a22 -> e32 + a23
R26: e18a21 -> e32 + a24
R27: a24 + a37 -> a26a27
R28: a27a23 -> a37 + a28
R29: a26a36 -> a29a30
R30: a29a36 -> a31a30
R31: a30 + a28 -> a25a33
R32: a31a25 -> a32 + a36
R33: a32a30 -> a34a36
R34: a34a33 -> a37 + a37
"""

#: Membrane states of Hutton (2007): a36 is plain membrane, a37 an anchor point.
MEMBRANE_STATES = (36, 37)
#: Number of states of the 2007 system; also the offset of the enzyme states.
MEMBRANE_N_STATES = 38


# --- rules ---------------------------------------------------------------------
@dataclass(frozen=True)
class Rule:
    """One reaction rule over 2 or 3 atoms.

    types/states describe the reactants in slot order (A, B[, C]); a type is a
    concrete letter or a variable in TYPE_VARS. bonds/new_bonds are (AB,) for a
    pair and (AB, BC, AC) for a triple. kind "plain" uses new_states as given;
    "readout" is 2007 R37 (the enzyme state is computed from the base it
    leaves) and "enzyme" is 2007 R40 (the catalyst's state encodes the whole
    reaction).
    """

    name: str
    types: tuple[str, ...]
    states: tuple[int, ...]
    bonds: tuple[bool, ...]
    new_states: tuple[int, ...]
    new_bonds: tuple[bool, ...]
    kind: str = "plain"
    cases: int = 1              # 1 in `cases` chance of firing (C++ TestProb)

    @property
    def n(self) -> int:
        return len(self.types)

    def to_text(self) -> str:
        return f"{self.name}: {_rule_text(self)}" if self.name else _rule_text(self)


def _atom_text(t: str, s: int) -> str:
    return f"{t}{s}"


def _rule_text(r: Rule) -> str:
    def side(states, bonds):
        parts = [_atom_text(t, s) for t, s in zip(r.types, states)]
        if r.n == 2:
            return parts[0] + ("" if bonds[0] else " + ") + parts[1]
        ab, bc, ac = bonds
        joined = parts[0] + ("-" if ab else "+") + parts[1] + ("-" if bc else "+") + parts[2]
        return joined + ("/ac" if ac else "")

    text = f"{side(r.states, r.bonds)} -> {side(r.new_states, r.new_bonds)}"
    if r.kind != "plain":
        text += f"  [{r.kind}]"
    if r.cases != 1:
        text += f"  [1 in {r.cases}]"
    return text


def _parse_side(side: str, text: str) -> tuple[list[tuple[str, int]], bool]:
    side = side.strip()
    if "+" in side:
        left, right = side.split("+", 1)
        atoms = _ATOM.findall(left.strip()) + _ATOM.findall(right.strip())
        bonded = False
    else:
        atoms = _ATOM.findall(side.replace("-", ""))
        bonded = True
    if len(atoms) != 2:
        raise ValueError(
            f"rule {text!r}: each side needs exactly two atoms written as a type "
            f"letter (a-f, or a variable x/y/z) and a state, e.g. 'e8 + e0' or 'x4y1'"
        )
    return [(t, int(s)) for t, s in atoms], bonded


def parse_rules(text: str) -> list[Rule]:
    """Parse pair rules in the notation of the papers, one per line.

    "R1: e8 + e0 -> e4e3"; juxtaposition (or '-') means bonded, ' + ' means not
    bonded; '#' starts a comment. Types never change, so both sides must name
    the same two types.
    """
    out: list[Rule] = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        name = ""
        if ":" in line.split("->")[0]:
            name, line = (part.strip() for part in line.split(":", 1))
        if line.count("->") != 1:
            raise ValueError(f"rule {raw.strip()!r} needs exactly one '->'")
        lhs, rhs = line.split("->")
        left, lb = _parse_side(lhs, raw.strip())
        right, rb = _parse_side(rhs, raw.strip())
        if [t for t, _ in left] != [t for t, _ in right]:
            raise ValueError(
                f"rule {raw.strip()!r}: an atom's type never changes, but the two "
                f"sides name different types"
            )
        rule = Rule(
            name=name,
            types=tuple(t for t, _ in left),
            states=tuple(s for _, s in left),
            bonds=(lb,),
            new_states=tuple(s for _, s in right),
            new_bonds=(rb,),
        )
        if rule.states == rule.new_states and rule.bonds == rule.new_bonds:
            raise ValueError(f"rule {raw.strip()!r} changes nothing")
        out.append(rule)
    return out


def _triple(name, types, states, bonds, new_states, new_bonds, kind="plain", cases=1) -> Rule:
    return Rule(name, tuple(types), tuple(states), tuple(bonds),
                tuple(new_states), tuple(new_bonds), kind, cases)


def membrane_rules(n_states: int = MEMBRANE_N_STATES, n_types: int = 6,
                   mutation_cases: int = 1000000) -> list[Rule]:
    """The 41 rules of Hutton (2007): R1-R34 pairs plus R35-R41."""
    rules = parse_rules(MEMBRANE_PAIR_RULES)
    m36 = MEMBRANE_STATES[0]
    # R35: a membrane spontaneously acquires an atom from, or loses one to, the
    # soup (three atoms; runs in both directions at equal rates).
    for i in MEMBRANE_STATES:
        for j in MEMBRANE_STATES:
            rules.append(_triple(
                "R35", ("a", "a", "a"), (i, j, 0), (True, False, False),
                (i, j, m36), (False, True, True), cases=100))
            rules.append(_triple(
                "R35", ("a", "a", "a"), (i, m36, j), (True, True, False),
                (i, 0, j), (False, False, True), cases=100))
    # R36: attach the seed of a new enzyme to the top of the gene string.
    rules.append(_triple("R36", ("x", "y", "z"), (32, 17, 0), (True, False, False),
                         (1, 35, n_states), (True, False, True)))
    # R37: the seed works down the gene, reading each base in base-|bases|.
    rules.append(_triple("R37", ("x", "y", "z"), (35, 17, -1), (True, False, True),
                         (1, 35, -1), (True, False, True), kind="readout"))
    # R38: release the finished enzyme as a bonded pair, so it stays in the cell.
    rules.append(_triple("R38", ("f", "x", "y"), (35, -1, 0), (True, False, False),
                         (32, -1, -1), (False, True, False), kind="release"))
    # R39: stop when every gene has been read out.
    rules += parse_rules("R39: f32a37 -> f1a37")
    # R40: the enzyme z(i) drives the reaction its state i encodes.
    rules.append(_triple("R40", ("z", "x", "y"), (-1, -1, -1), (False, False, False),
                         (-1, -1, -1), (False, False, False), kind="enzyme"))
    # R41: a chain of atoms in state 1 gains or loses an atom (mutation), with
    # the low probability the paper leaves to the user.
    rules.append(_triple("R41", ("x", "y", "z"), (1, 1, 0), (True, False, False),
                         (1, 1, 1), (False, True, True), cases=mutation_cases))
    rules.append(_triple("R41", ("x", "y", "z"), (1, 1, 1), (True, True, False),
                         (1, 0, 1), (False, False, True), cases=mutation_cases))
    return rules


def rule_set(name: str, text: str, n_states: int, n_types: int,
             mutation_cases: int = 1000000) -> list[Rule]:
    if name == "replicator":
        return parse_rules(REPLICATOR_RULES)
    if name == "membrane":
        return membrane_rules(n_states, n_types, mutation_cases)
    if not text.strip():
        raise ValueError("rules='custom' needs a non-empty rule_text")
    return parse_rules(text)


def explicit_rules(rules, types: str = ATOM_TYPES) -> list[tuple]:
    """Expand type variables: every rule written out without variables.

    Returns (types, states, bonds, new_states, new_bonds) tuples, one per
    distinct reaction. Different variable letters are unconstrained, so x and y
    may take the same type (Hutton 2002, which counts 188 of these for table 1).
    """
    out: dict[tuple, tuple] = {}
    for r in rules:
        variables = sorted({t for t in r.types if t in TYPE_VARS})
        for combo in _assignments(variables, types):
            concrete = tuple(combo.get(t, t) for t in r.types)
            key = (concrete, r.states, r.bonds)
            out.setdefault(key, (concrete, r.states, r.bonds, r.new_states, r.new_bonds))
    return list(out.values())


def _assignments(variables, types):
    if not variables:
        yield {}
        return
    head, rest = variables[0], variables[1:]
    for t in types:
        for tail in _assignments(rest, types):
            yield {head: t, **tail}


# --- the enzyme encoding (2007 R40 / book eq. 11.15) ---------------------------
def encode_enzyme(g: int, h: int, j: int, k: int, x: int, y: int, b1: int, b2: int,
                  *, n_states: int, n_types: int = 6) -> int:
    """The state i of a catalyst z(i) that drives x(g) [b1] y(h) -> x(j) [b2] y(k).

    i = 2(2(T(T(S(S(S g + h) + j) + k) + x) + y) + b1) + b2 + S, with types
    encoded a = 0 .. f = 5 and bonds 0 unbonded, 1 bonded.
    """
    for value, limit, what in ((g, n_states, "g"), (h, n_states, "h"),
                               (j, n_states, "j"), (k, n_states, "k"),
                               (x, n_types, "x"), (y, n_types, "y")):
        if not 0 <= value < limit:
            raise ValueError(f"{what} = {value} is out of range [0, {limit})")
    if b1 not in (0, 1) or b2 not in (0, 1):
        raise ValueError(f"bonds must be 0 or 1, got {b1!r} and {b2!r}")
    v = (((g * n_states + h) * n_states + j) * n_states + k) * n_types + x
    return ((v * n_types + y) * 2 + b1) * 2 + b2 + n_states


def decode_enzyme(i: int, *, n_states: int, n_types: int = 6) -> dict:
    """Invert `encode_enzyme`: the reaction that the catalyst state i encodes."""
    if i < n_states:
        raise ValueError(f"enzyme states start at {n_states}, got {i}")
    rest, b2 = divmod(i - n_states, 2)
    rest, b1 = divmod(rest, 2)
    rest, y = divmod(rest, n_types)
    rest, x = divmod(rest, n_types)
    rest, k = divmod(rest, n_states)
    rest, j = divmod(rest, n_states)
    g, h = divmod(rest, n_states)
    if g >= n_states:
        raise ValueError(f"state {i} is not a valid enzyme for S = {n_states}")
    return {"g": g, "h": h, "j": j, "k": k, "x": x, "y": y, "b1": b1, "b2": b2}


def gene_to_enzyme(bases: str, *, n_states: int, alphabet: str = "abcd") -> int:
    """Read a gene of base atoms into an enzyme state (2007 R37).

    The seed starts in state S and each base x moves it to
    j = |alphabet| (i - S) + val(x) + S, i.e. the base-|alphabet| value of the
    string offset by S.
    """
    i = n_states
    for base in bases:
        if base not in alphabet:
            raise ValueError(f"gene base {base!r} is not one of {alphabet!r}")
        i = len(alphabet) * (i - n_states) + alphabet.index(base) + n_states
    return i


# --- molecules ------------------------------------------------------------------
@dataclass(frozen=True)
class Mol:
    """A molecule: atom labels in canonical order plus its bonds."""

    labels: tuple[tuple[str, int], ...]
    edges: tuple[tuple[int, int], ...]

    @property
    def id(self) -> str:
        text = "-".join(_atom_text(t, s) for t, s in self.labels)
        path = tuple((i, i + 1) for i in range(len(self.labels) - 1))
        if self.edges == path:
            return text
        ring = path + ((0, len(self.labels) - 1),) if len(self.labels) > 2 else ()
        if ring and self.edges == tuple(sorted(ring)):
            return text + "/ring"
        return text + "/" + ",".join(f"{u}.{v}" for u, v in self.edges)

    @property
    def atom_counts(self) -> Counter:
        return Counter(t for t, _ in self.labels)


def canonical(labels: list[tuple[str, int]], edges) -> Mol:
    """Canonical form of a labelled molecule, invariant under its symmetries."""
    n = len(labels)
    adj: list[list[int]] = [[] for _ in range(n)]
    clean = set()
    for u, v in edges:
        if u == v:
            raise ValueError("an atom cannot be bonded to itself")
        clean.add((min(u, v), max(u, v)))
    for u, v in clean:
        adj[u].append(v)
        adj[v].append(u)
    if n == 1:
        return Mol((labels[0],), ())
    degrees = [len(a) for a in adj]
    order = None
    if max(degrees) <= 2:
        if len(clean) == n - 1:                         # a simple chain
            order = min(_walk(adj, e, labels) for e in range(n) if degrees[e] == 1)
        elif len(clean) == n:                           # a ring
            order = min(_ring(adj, labels, s, d) for s in range(n) for d in (0, 1))
    if order is None:
        order = _canonical_order(labels, adj)
    index = {v: i for i, v in enumerate(order)}
    new_edges = sorted(
        (min(index[u], index[v]), max(index[u], index[v])) for u, v in clean
    )
    return Mol(tuple(labels[v] for v in order), tuple(new_edges))


def _walk(adj, start, labels):
    """Traverse a chain from an end; returns the order keyed by its labels."""
    order, previous, current = [start], None, start
    while True:
        nxt = [u for u in adj[current] if u != previous]
        if not nxt:
            break
        previous, current = current, nxt[0]
        order.append(current)
    return _Keyed(order, labels)


def _ring(adj, labels, start, direction):
    order, previous, current = [start], None, start
    for _ in range(len(labels) - 1):
        nxt = [u for u in adj[current] if u != previous]
        previous, current = current, nxt[direction % len(nxt)]
        order.append(current)
    return _Keyed(order, labels)


class _Keyed(list):
    """A vertex order compared by the label sequence it produces."""

    def __init__(self, order, labels):
        super().__init__(order)
        self._key = tuple(labels[v] for v in order)

    def __lt__(self, other):
        return self._key < other._key


def _refine(colours, adj):
    while True:
        signature = [(colours[v], tuple(sorted(colours[u] for u in adj[v])))
                     for v in range(len(colours))]
        ranks = {c: i for i, c in enumerate(sorted(set(signature)))}
        new = [ranks[c] for c in signature]
        if new == colours:
            return colours
        colours = new


def _canonical_order(labels, adj):
    """Canonical vertex order by individualization-refinement (any graph)."""
    n = len(labels)
    ranks = {c: i for i, c in enumerate(sorted(set(labels)))}
    best: list = [None, None]

    def search(colours):
        colours = _refine(colours, adj)
        cells: dict[int, list[int]] = {}
        for v, c in enumerate(colours):
            cells.setdefault(c, []).append(v)
        target = next((cells[c] for c in sorted(cells) if len(cells[c]) > 1), None)
        if target is None:
            order = sorted(range(n), key=lambda v: colours[v])
            index = {v: i for i, v in enumerate(order)}
            certificate = (
                tuple(labels[v] for v in order),
                tuple(sorted((min(index[u], index[v]), max(index[u], index[v]))
                             for u in range(n) for v in adj[u] if u < v)),
            )
            if best[0] is None or certificate < best[0]:
                best[0], best[1] = certificate, order
            return
        for v in target:
            branch = list(colours)
            branch[v] = -1
            search(branch)

    search([ranks[c] for c in labels])
    return best[1]


# --- the world -------------------------------------------------------------------
class _Rand:
    """Buffered draws from a numpy Generator (one call per few thousand draws)."""

    def __init__(self, rng, size: int = 1 << 15):
        self._rng, self._size, self._buf = rng, size, []

    def random(self) -> float:
        if not self._buf:
            self._buf = self._rng.random(self._size).tolist()
        return self._buf.pop()

    def below(self, n: int) -> int:
        return int(self.random() * n)

    def choice(self, items):
        return items[self.below(len(items))]

    def shuffled(self, items: list) -> list:
        out = list(items)
        for i in range(len(out) - 1, 0, -1):
            j = self.below(i + 1)
            out[i], out[j] = out[j], out[i]
        return out


MOORE = [(-1, -1), (0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0)]
VON_NEUMANN = [(0, -1), (1, 0), (0, 1), (-1, 0)]


class World:
    """Atoms with types, states, bonds and positions, and the rules acting on them."""

    def __init__(self, rules, rand: _Rand, *, n_states: int, bases: str,
                 lattice: bool, width: float, height: float,
                 reaction_neighbourhood: str = "von-neumann", bond_range: int = 1,
                 reaction_radius: float = 2.5, radius: float = 6.0):
        self.rules = list(rules)
        self.rand = rand
        self.n_states = n_states
        self.bases = bases
        self.lattice = lattice
        self.width, self.height = width, height
        self.bond_range = bond_range
        self.reaction_radius = reaction_radius
        self.radius = radius
        self.moves = MOORE
        self.react_offsets = VON_NEUMANN if reaction_neighbourhood == "von-neumann" else MOORE
        # atoms
        self.type: list[str] = []
        self.state: list[int] = []
        self.x: list[float] = []
        self.y: list[float] = []
        self.vx: list[float] = []
        self.vy: list[float] = []
        self.bonds: list[set[int]] = []
        self.occupied: dict[tuple[int, int], int] = {}
        # rules indexed by the state of the first reactant (-1 = any)
        self.by_state: dict[int, list[Rule]] = {}
        for rule in self.rules:
            self.by_state.setdefault(rule.states[0], []).append(rule)
        self.wildcard = self.by_state.get(-1, [])
        self.fired: dict[tuple, list] = {}
        self.tally = Tally()            # the same firings, flushed into frames
        self.rule_counts: Counter = Counter()
        self.seen: dict[str, Mol] = {}
        self.floods = 0
        self.dissolved = 0

    # -- construction --------------------------------------------------------
    def add(self, t: str, s: int, x: float, y: float) -> int:
        i = len(self.type)
        self.type.append(t)
        self.state.append(s)
        self.x.append(x)
        self.y.append(y)
        self.bonds.append(set())
        if self.lattice:
            self.occupied[(int(x), int(y))] = i
            self.vx.append(0.0)
            self.vy.append(0.0)
        else:
            angle = self.rand.random() * 2 * math.pi
            speed = self.radius * 0.4
            self.vx.append(speed * math.cos(angle))
            self.vy.append(speed * math.sin(angle))
        return i

    def bond(self, i: int, j: int) -> None:
        self.bonds[i].add(j)
        self.bonds[j].add(i)

    # -- molecules -----------------------------------------------------------
    def component(self, start: int) -> list[int]:
        seen, queue = {start}, deque([start])
        while queue:
            v = queue.popleft()
            for u in self.bonds[v]:
                if u not in seen:
                    seen.add(u)
                    queue.append(u)
        return sorted(seen)

    def molecule(self, atoms: list[int]) -> Mol:
        index = {a: i for i, a in enumerate(atoms)}
        labels = [(self.type[a], self.state[a]) for a in atoms]
        edges = [(index[a], index[b]) for a in atoms for b in self.bonds[a] if a < b]
        mol = canonical(labels, edges)
        self.seen.setdefault(mol.id, mol)
        return mol

    def molecules_of(self, atoms) -> list[Mol]:
        out, done = [], set()
        for a in atoms:
            if a in done:
                continue
            comp = self.component(a)
            done.update(comp)
            out.append(self.molecule(comp))
        return out

    def population(self) -> list[Mol]:
        out, done = [], set()
        for a in range(len(self.type)):
            if a not in done:
                comp = self.component(a)
                done.update(comp)
                out.append(self.molecule(comp))
        return out

    # -- matching ------------------------------------------------------------
    def near(self, i: int) -> list[int]:
        if self.lattice:
            xi, yi = int(self.x[i]), int(self.y[i])
            out = []
            for dx, dy in self.react_offsets:
                j = self.occupied.get((xi + dx, yi + dy))
                if j is not None:
                    out.append(j)
            return out
        limit = (self.reaction_radius * self.radius) ** 2
        xi, yi = self.x[i], self.y[i]
        return [j for j in range(len(self.type))
                if j != i and (self.x[j] - xi) ** 2 + (self.y[j] - yi) ** 2 < limit]

    def _close(self, i: int, j: int) -> bool:
        if self.lattice:
            dx, dy = int(self.x[i]) - int(self.x[j]), int(self.y[i]) - int(self.y[j])
            return (dx, dy) in self.react_offsets
        limit = (self.reaction_radius * self.radius) ** 2
        return (self.x[i] - self.x[j]) ** 2 + (self.y[i] - self.y[j]) ** 2 < limit

    def _types_ok(self, rule: Rule, atoms: tuple[int, ...]) -> bool:
        bound: dict[str, str] = {}
        for pattern, a in zip(rule.types, atoms):
            actual = self.type[a]
            if pattern in TYPE_VARS:
                if bound.setdefault(pattern, actual) != actual:
                    return False
            elif pattern != actual:
                return False
        return True

    def _bonds_ok(self, rule: Rule, atoms: tuple[int, ...]) -> bool:
        pairs = ((0, 1),) if rule.n == 2 else ((0, 1), (1, 2), (0, 2))
        for want, (u, v) in zip(rule.bonds, pairs):
            if (atoms[v] in self.bonds[atoms[u]]) != want:
                return False
        return True

    def matches(self, i: int) -> list[tuple[Rule, tuple[int, ...], tuple]]:
        """Every (rule, atoms, outcome) that applies with atom i in the first slot."""
        candidates = self.by_state.get(self.state[i], ())
        out = []
        neighbours = None
        for rule in (*candidates, *self.wildcard):
            if rule.types[0] not in TYPE_VARS and rule.types[0] != self.type[i]:
                continue
            if rule.kind == "plain" and rule.states[0] != self.state[i]:
                continue
            if neighbours is None:
                neighbours = self.near(i)
            for j in neighbours:
                if j == i:
                    continue
                if rule.kind == "plain" and rule.states[1] != self.state[j]:
                    continue
                if rule.n == 2:
                    atoms = (i, j)
                    if not self._types_ok(rule, atoms) or not self._bonds_ok(rule, atoms):
                        continue
                    outcome = self._outcome(rule, atoms)
                    if outcome is not None:
                        out.append((rule, atoms, outcome))
                    continue
                for k in neighbours:
                    if k == i or k == j or not self._close(j, k):
                        continue
                    atoms = (i, j, k)
                    if not self._types_ok(rule, atoms):
                        continue
                    if rule.kind == "plain" and rule.states[2] != self.state[k]:
                        continue
                    if not self._bonds_ok(rule, atoms):
                        continue
                    outcome = self._outcome(rule, atoms)
                    if outcome is not None:
                        out.append((rule, atoms, outcome))
        return out

    def _outcome(self, rule: Rule, atoms):
        """(new_states, new_bonds) for this match, or None if it does not apply."""
        if rule.kind == "plain":
            return rule.new_states, rule.new_bonds
        if rule.kind == "readout":
            # R37: x35 y17 z(i) -> x1 y35 z(j), i >= S, x a base type
            a, b, c = atoms
            if self.state[a] != 35 or self.state[b] != 17 or self.state[c] < self.n_states:
                return None
            if self.type[a] not in self.bases:
                return None
            i = self.state[c]
            j = (len(self.bases) * (i - self.n_states)
                 + self.bases.index(self.type[a]) + self.n_states)
            return (1, 35, j), (True, False, True)
        if rule.kind == "release":
            # R38: f35 x(i) y0 -> f32 x(i) y(i), the enzyme leaves as a pair
            a, b, c = atoms
            if self.state[a] != 35 or self.state[b] < self.n_states or self.state[c] != 0:
                return None
            return (32, self.state[b], self.state[b]), (False, True, False)
        if rule.kind == "enzyme":
            # R40: z(i) catalyses x(g) [b1] y(h) -> x(j) [b2] y(k)
            z, a, b = atoms
            i = self.state[z]
            if i < self.n_states:
                return None
            try:
                spec = decode_enzyme(i, n_states=self.n_states, n_types=len(ATOM_TYPES))
            except ValueError:
                return None
            if self.state[a] != spec["g"] or self.state[b] != spec["h"]:
                return None
            if ATOM_TYPES.index(self.type[a]) != spec["x"]:
                return None
            if ATOM_TYPES.index(self.type[b]) != spec["y"]:
                return None
            if (b in self.bonds[a]) != bool(spec["b1"]):
                return None
            if spec["g"] == spec["j"] and spec["h"] == spec["k"] and spec["b1"] == spec["b2"]:
                return None
            return (i, spec["j"], spec["k"]), (False, bool(spec["b2"]), False)
        raise ValueError(f"unknown rule kind {rule.kind!r}")

    # -- firing --------------------------------------------------------------
    def apply(self, rule: Rule, atoms, outcome) -> None:
        before = self.molecules_of(atoms)
        new_states, new_bonds = outcome
        for a, s in zip(atoms, new_states):
            self.state[a] = s
        pairs = ((0, 1),) if rule.n == 2 else ((0, 1), (1, 2), (0, 2))
        for want, (u, v) in zip(new_bonds, pairs):
            i, j = atoms[u], atoms[v]
            if want and j not in self.bonds[i]:
                self.bond(i, j)
            elif not want and j in self.bonds[i]:
                self.bonds[i].discard(j)
                self.bonds[j].discard(i)
        after = self.molecules_of(atoms)
        self.record(before, after, rule)

    def record(self, before, after, rule) -> None:
        left = Counter(m.id for m in before)
        right = Counter(m.id for m in after)
        if left == right:
            return
        key = (frozenset(left.items()), frozenset(right.items()))
        entry = self.fired.setdefault(key, [left, right, 0, Counter()])
        entry[2] += 1
        self.tally.add(tuple(left.elements()), tuple(right.elements()))
        entry[3][rule.name or _rule_text(rule)] += 1
        self.rule_counts[rule.name or _rule_text(rule)] += 1

    # -- dynamics ------------------------------------------------------------
    def react_phase(self) -> None:
        for i in self.rand.shuffled(range(len(self.type))):
            if not self.by_state.get(self.state[i]) and not self.wildcard:
                continue
            found = self.matches(i)
            if not found:
                continue
            rule, atoms, outcome = self.rand.choice(found)
            if rule.cases > 1 and self.rand.random() * rule.cases >= 1.0:
                continue
            self.apply(rule, atoms, outcome)

    def move_phase(self) -> None:
        if self.lattice:
            self._move_lattice()
        else:
            self._move_continuous()

    def _move_lattice(self) -> None:
        rand, occupied, span = self.rand, self.occupied, self.bond_range
        w, h = int(self.width), int(self.height)
        for i in rand.shuffled(range(len(self.type))):
            dx, dy = self.moves[rand.below(8)]
            nx, ny = int(self.x[i]) + dx, int(self.y[i]) + dy
            if not (0 <= nx < w and 0 <= ny < h) or (nx, ny) in occupied:
                continue
            ok = True
            for j in self.bonds[i]:
                if max(abs(nx - int(self.x[j])), abs(ny - int(self.y[j]))) > span:
                    ok = False
                    break
            if not ok:
                continue
            del occupied[(int(self.x[i]), int(self.y[i]))]
            self.x[i], self.y[i] = float(nx), float(ny)
            occupied[(nx, ny)] = i

    def _move_continuous(self) -> None:
        """The 2007 continuous physics: bond springs, volume exclusion, walls."""
        r = self.radius
        phys2 = (2.0 * r) ** 2
        top = r * 0.4
        n = len(self.type)
        for i in range(n):
            fx = fy = 0.0
            free_food = not self.bonds[i] and self.state[i] == 0
            if not free_food:
                for j in range(n):
                    if j == i:
                        continue
                    dx, dy = self.x[i] - self.x[j], self.y[i] - self.y[j]
                    d2 = dx * dx + dy * dy
                    if d2 >= phys2 or d2 == 0.0:
                        continue
                    if not self.bonds[j] and self.state[j] == 0:
                        continue          # unbonded food passes through (semi-permeable)
                    push = 0.4 * (2.0 * r / math.sqrt(d2) - 1.0)
                    fx += dx * push
                    fy += dy * push
            for j in self.bonds[i]:
                dx, dy = self.x[j] - self.x[i], self.y[j] - self.y[i]
                d2 = dx * dx + dy * dy
                if d2 > phys2:
                    pull = 0.4 * (math.sqrt(d2) / (2.0 * r) - 1.0)
                    fx += dx * pull
                    fy += dy * pull
            if self.x[i] < r:
                fx += max(0.0, r / max(self.x[i], 0.1) - 1.0) * self.x[i]
            elif self.x[i] > self.width - r:
                fx += (self.x[i] - self.width) * max(0.0, r / max(self.width - self.x[i], 0.1) - 1.0)
            if self.y[i] < r:
                fy += max(0.0, r / max(self.y[i], 0.1) - 1.0) * self.y[i]
            elif self.y[i] > self.height - r:
                fy += (self.y[i] - self.height) * max(0.0, r / max(self.height - self.y[i], 0.1) - 1.0)
            vx, vy = self.vx[i] + fx, self.vy[i] + fy
            speed = math.hypot(vx, vy)
            if speed > top:
                vx, vy = vx * top / speed, vy * top / speed
            self.vx[i], self.vy[i] = vx, vy
        for i in range(n):
            nx, ny = self.x[i] + self.vx[i], self.y[i] + self.vy[i]
            if not (0.0 <= nx <= self.width and 0.0 <= ny <= self.height):
                self.vx[i] = self.vy[i] = 0.0
                nx = min(max(nx, 0.0), self.width)
                ny = min(max(ny, 0.0), self.height)
            self.x[i], self.y[i] = nx, ny

    def cosmic_rays(self, probability: float) -> None:
        """2002 experiment 3: perturb an atom's state, leaving type and bonds."""
        for i in range(len(self.type)):
            if self.rand.random() < probability:
                self.state[i] = self.rand.below(self.n_states)

    def flood(self, sectors: int) -> int:
        """Dissolve every atom in one sector: break its bonds, reset it to state 0."""
        quarter = self.floods % sectors
        if sectors == 2:
            x0, y0 = (self.width / 2.0) * quarter, 0.0
            w, h = self.width / 2.0, self.height
        else:
            x0 = (self.width / 2.0) * (quarter % 2)
            y0 = (self.height / 2.0) * (quarter // 2)
            w, h = self.width / 2.0, self.height / 2.0
        hit = [i for i in range(len(self.type))
               if x0 <= self.x[i] <= x0 + w and y0 <= self.y[i] <= y0 + h]
        dissolved = [m.id for m in self.molecules_of(hit) if len(m.labels) > 1]
        for i in hit:
            for j in list(self.bonds[i]):
                self.bonds[i].discard(j)
                self.bonds[j].discard(i)
            self.state[i] = 0
        self.floods += 1
        return len(dissolved)

    def run(self, steps: int, *, flood_period: int, flood_sectors: int,
            cosmic_ray: float, samples: int = 200):
        """Run `steps` time steps, yielding the index of the steps after which
        to sample: the first one and every steps // samples from there.
        `steps=0` runs the world on until the caller stops reading."""
        every = max(1, steps // samples) if steps else 1
        for step in (t - 1 for t in ticks(steps)):
            self.react_phase()
            self.move_phase()
            if cosmic_ray > 0.0:
                self.cosmic_rays(cosmic_ray)
            if flood_period and step > 0 and step % flood_period == 0:
                self.dissolved += self.flood(flood_sectors)
            if step % every == 0:
                yield step


# --- molecule-level reactions for the closure -----------------------------------
def _mol_world(rules, rand, mols, n_states, bases) -> World:
    world = World(rules, rand, n_states=n_states, bases=bases, lattice=False,
                  width=1e9, height=1e9, reaction_radius=1e9)
    for mol in mols:
        base = len(world.type)
        for t, s in mol.labels:
            world.add(t, s, 0.0, 0.0)
        for u, v in mol.edges:
            world.bond(base + u, base + v)
    return world


def _react_molecules(rules, rand, n_states, bases, mols, require_all: bool):
    """Every distinct outcome of letting these molecules react once."""
    outcomes: dict[tuple, tuple] = {}
    owner: list[int] = []
    for k, mol in enumerate(mols):
        owner += [k] * len(mol.labels)
    world = _mol_world(rules, rand, mols, n_states, bases)
    states, bonds = list(world.state), [set(b) for b in world.bonds]
    found = [match for i in range(len(world.type)) for match in world.matches(i)]
    for rule, atoms, outcome in found:
        if require_all and len({owner[a] for a in atoms}) < len(mols):
            continue
        world.state, world.bonds = list(states), [set(b) for b in bonds]
        new_states, new_bonds = outcome
        for a, s in zip(atoms, new_states):
            world.state[a] = s
        pairs = ((0, 1),) if rule.n == 2 else ((0, 1), (1, 2), (0, 2))
        for want, (u, v) in zip(new_bonds, pairs):
            a, b = atoms[u], atoms[v]
            if want:
                world.bond(a, b)
            else:
                world.bonds[a].discard(b)
                world.bonds[b].discard(a)
        products = tuple(sorted(m.id for m in world.population()))
        if products == tuple(sorted(m.id for m in mols)):
            continue
        outcomes.setdefault(products, products)
    world.state, world.bonds = states, bonds
    return list(outcomes.values())


# --- generate ---------------------------------------------------------------------
def parse_molecule(text: str) -> Mol:
    """Parse a chain like 'e8-a1-b1-f1' (atoms bonded in the order given)."""
    tokens = [t for t in text.replace(" ", "").split("-") if t]
    if not tokens:
        raise ValueError("seed_molecule must name at least one atom, e.g. 'e8-a1-b1-f1'")
    labels = []
    for token in tokens:
        match = re.fullmatch(r"([a-f])(\d+)", token)
        if not match:
            raise ValueError(
                f"seed_molecule: {token!r} is not an atom (a type letter a-f and a state, e.g. 'a1')"
            )
        labels.append((match.group(1), int(match.group(2))))
    edges = [(i, i + 1) for i in range(len(labels) - 1)]
    return canonical(labels, edges)


def _parse_seed(text: str) -> tuple[bool, Mol]:
    """"cell:<gene>" asks for the 2007 cell; anything else is a bare molecule."""
    body = text.strip()
    if body.lower().startswith("cell:"):
        return True, parse_molecule(body.split(":", 1)[1])
    return False, parse_molecule(body)


def _chain_atoms(mol: Mol):
    """Layout of a bare molecule: a straight chain of lattice points."""
    points = [(k, 0) for k in range(len(mol.labels))]
    return points, list(mol.labels), [(u, v) for u, v in mol.edges]


def _cell_atoms(gene: Mol):
    """The 2007 cell (fig. 2): a membrane loop with the gene anchored at two a37."""
    n = len(gene.labels)
    right = n + 1
    ring = ([(x, -2) for x in range(right + 1)]
            + [(right, y) for y in range(-1, 3)]
            + [(x, 2) for x in range(right - 1, -1, -1)]
            + [(0, y) for y in range(1, -2, -1)])
    labels = [("a", 37) if pt in ((0, 0), (right, 0)) else ("a", 36) for pt in ring]
    bonds = [(i, (i + 1) % len(ring)) for i in range(len(ring))]
    at = {pt: i for i, pt in enumerate(ring)}
    points = list(ring)
    base = len(points)
    for k, label in enumerate(gene.labels):
        points.append((k + 1, 0))
        labels.append(label)
        if k:
            bonds.append((base + k - 1, base + k))
    bonds.append((at[(0, 0)], base))
    bonds.append((at[(right, 0)], base + n - 1))
    return points, labels, bonds


def _setup(p, mutation_cases: int = 1000000):
    """Validate the chemistry parameters: (types, bases, rules, cell?, seed molecule)."""
    if p.n_types < 2:
        raise ValueError(f"n_types must be at least 2, got {p.n_types}")
    types = ATOM_TYPES[:p.n_types]
    bases = types[:-2] if p.n_types > 2 else types
    rules = rule_set(p.rules, p.rule_text, p.n_states, p.n_types, mutation_cases)
    # Ordinary states run 0..n_states-1; enzyme states start at n_states, so only
    # the states a rule matches on bound n_states (R36 *produces* the seed n_states).
    needed = max(s for r in rules for s in r.states if s >= 0)
    if p.n_states <= needed:
        raise ValueError(
            f"rules={p.rules!r} reacts with atoms in state {needed}, so n_states "
            f"must be greater than {needed}, got {p.n_states}"
        )
    cell, seed = _parse_seed(p.seed_molecule)
    for t, s in seed.labels:
        if t not in types:
            raise ValueError(f"seed_molecule uses type {t!r}, but n_types = {p.n_types}")
        if s >= p.n_states:
            raise ValueError(f"seed_molecule uses state {s}, but n_states = {p.n_states}")
    return types, bases, rules, cell, seed


def generate(p, rng):
    """The closure: every molecule-level reaction reachable from the seed
    molecule and one food atom of each type, cut off by max_species."""
    # Rule probabilities (R35, R41) do not matter here: every outcome is a reaction.
    types, bases, rules, cell, seed = _setup(p)
    if cell:
        raise ValueError(
            "the closure starts from a molecule, not a cell; drop the 'cell:' prefix, "
            "or run the cell with chemart.evolve"
        )
    return _closure(p, rules, seed, types, bases, _Rand(rng))


def evolve(p, rng):
    """The 2D world: a frame at the start, after the first step, every
    steps // 200 steps from there, and at the end."""
    types, bases, rules, cell, seed = _setup(p, p.mutation_cases)
    return (yield from _spatial(p, rules, seed, types, bases, _Rand(rng), cell))


def _conservation(species: list[Mol], types: str) -> list[dict]:
    out = []
    for t in types:
        vector = {m.id: m.atom_counts.get(t, 0) for m in species}
        if any(vector.values()):
            out.append({"name": f"atoms_{t}", "vector": vector})
    out.append({"name": "atoms", "vector": {m.id: len(m.labels) for m in species}})
    return out


def _space(p) -> dict:
    if p.space == "continuous":
        return {"model": "continuous-2d", "width": float(p.width), "height": float(p.height),
                "atom_radius": 6.0, "reaction_radius": float(p.reaction_radius) * 6.0,
                "semi_permeable": "unbonded atoms in state 0 feel no volume exclusion"}
    return {"model": "lattice-2d", "width": int(p.width), "height": int(p.height),
            "movement_neighbourhood": "moore",
            "reaction_neighbourhood": "von-neumann" if p.space == "lattice-vn" else "moore",
            "bond_range": 1 if p.space == "lattice-vn" else 2,
            "occupancy": "at most one atom per lattice point"}


def _spatial(p, rules, seed, types, bases, rand, cell=False):
    lattice = p.space != "continuous"
    world = World(
        rules, rand, n_states=p.n_states, bases=bases, lattice=lattice,
        width=p.width, height=p.height,
        reaction_neighbourhood="von-neumann" if p.space == "lattice-vn" else "moore",
        bond_range=1 if p.space == "lattice-vn" else 2,
        reaction_radius=p.reaction_radius,
    )
    _place(world, _cell_atoms(seed) if cell else _chain_atoms(seed), types, p, rand, lattice)

    def frame(t: int) -> Frame:
        # population() also records every molecule present as a species
        mols = world.population()
        return Frame(t=float(t), state={i: float(n) for i, n in Counter(m.id for m in mols).items()},
                     fired=world.tally.flush(),
                     observables={"molecules": sum(1 for m in mols if len(m.labels) > 1)})

    first = frame(0)
    initial = dict(first.state)
    yield first
    done = 0
    for step in world.run(p.steps, flood_period=p.flood_period,
                          flood_sectors=2 if p.flood_sectors == "halves" else 4,
                          cosmic_ray=p.cosmic_ray):
        done = step + 1
        yield frame(done)
    if done < p.steps:
        yield frame(p.steps)
    final = Counter(m.id for m in world.population())

    reactions = [
        Reaction(dict(left), dict(right), count=count)
        for left, right, count, _ in world.fired.values()
    ]
    species = [world.seen[i] for i in sorted(world.seen)]
    return Network(
        species=[Species(m.id, structure=m.id) for m in species],
        reactions=reactions,
        status="observed",
        initial_state={i: n for i, n in sorted(initial.items())},
        extras={
            "space": _space(p),
            "conservation": _conservation(species, types),
            "rules": [r.to_text() for r in rules],
            "rule_counts": dict(world.rule_counts.most_common()),
            "reaction_rules": [
                {"reactants": dict(left), "products": dict(right), "count": count,
                 "rules": dict(by_rule.most_common())}
                for left, right, count, by_rule in world.fired.values()
            ],
            "final_state": {i: n for i, n in sorted(final.items())},
            "floods": world.floods,
            "dissolved_by_flood": world.dissolved,
        },
    )


def _place(world, layout, types, p, rand, lattice):
    """Put the seed structure in the middle of the world, then scatter the food."""
    points, labels, bonds = layout
    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    x0, y0 = min(xs), min(ys)
    span_x, span_y = max(xs) - x0 + 1, max(ys) - y0 + 1
    if lattice:
        w, h = int(p.width), int(p.height)
        if span_x > w or span_y > h:
            raise ValueError(
                f"the seed needs a {span_x}x{span_y} world, but width x height is {w}x{h}"
            )
        if len(points) + p.food > w * h:
            raise ValueError(
                f"the lattice holds {w * h} atoms but the seed plus food needs "
                f"{len(points) + p.food}; raise width/height or lower food"
            )
        dx, dy = (w - span_x) // 2 - x0, (h - span_y) // 2 - y0
        index = [world.add(t, s, x + dx, y + dy)
                 for (x, y), (t, s) in zip(points, labels)]
    else:
        step = 2.0 * world.radius
        dx = p.width / 2.0 - (x0 + span_x / 2.0) * step
        dy = p.height / 2.0 - (y0 + span_y / 2.0) * step
        index = [world.add(t, s, x * step + dx, y * step + dy)
                 for (x, y), (t, s) in zip(points, labels)]
    for u, v in bonds:
        world.bond(index[u], index[v])
    if lattice:
        placed, attempts = 0, 0
        while placed < p.food and attempts < 200 * (p.food + 1):
            attempts += 1
            x, y = rand.below(int(p.width)), rand.below(int(p.height))
            if (x, y) in world.occupied:
                continue
            world.add(types[rand.below(len(types))], 0, x, y)
            placed += 1
    else:
        for _ in range(p.food):
            world.add(types[rand.below(len(types))], 0,
                      rand.random() * p.width, rand.random() * p.height)


def _closure(p, rules, seed, types, bases, rand):
    food = [canonical([(t, 0)], []) for t in types]
    seeds = [seed] + food
    cache = {m.id: m for m in seeds}

    def molecule(i: str) -> Mol:
        mol = cache.get(i)
        if mol is None:
            mol = cache[i] = _mol_from_id(i)
        return mol

    def react(*ids_in):
        return _react_molecules(
            rules, rand, p.n_states, bases, [molecule(i) for i in ids_in],
            require_all=len(ids_in) > 1,
        ) or None

    ids, found, status = expand(
        react, [m.id for m in seeds], arity=[1, 2],
        max_species=p.max_species, ordered=False, alternatives=True,
    )
    species = [molecule(i) for i in ids]
    kept = set(ids)                 # a budget below the seed count cuts seeds too
    return Network(
        species=[Species(m.id, structure=m.id) for m in species],
        reactions=[Reaction.of(left, right) for left, right in found],
        status=status,
        initial_state={m.id: 1.0 for m in seeds if m.id in kept},
        extras={
            "conservation": _conservation(species, types),
            "rules": [r.to_text() for r in rules],
            "seed": [m.id for m in seeds],
            "note": "closure: molecules may meet in any way, with no spatial constraint",
        },
    )


def _mol_from_id(text: str) -> Mol:
    """Rebuild a molecule from its canonical id."""
    body, _, edge_text = text.partition("/")
    labels = []
    for token in body.split("-"):
        match = re.fullmatch(r"([a-f])(\d+)", token)
        labels.append((match.group(1), int(match.group(2))))
    if not edge_text:
        edges = [(i, i + 1) for i in range(len(labels) - 1)]
    elif edge_text == "ring":
        edges = [(i, (i + 1) % len(labels)) for i in range(len(labels))]
    else:
        edges = [tuple(int(v) for v in pair.split(".")) for pair in edge_text.split(",")]
    return canonical(labels, edges)
