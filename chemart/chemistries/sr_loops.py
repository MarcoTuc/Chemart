"""Self-replicating loops in cellular automata (book 10.4, 10.7.2; ecology 8.2.3).

Catalog id: sr-loops.

The cellular automaton is implemented honestly: a square lattice of cells with
`n_states` states, the von Neumann neighbourhood and a published transition
table, applied synchronously. Five tables are offered as named choices:

    langton    Langton (1984), 8 states, 219 transitions, the 86-cell loop
    byl        Byl (1989), 6 states, the 12-cell loop
    reggia-1   Chou & Reggia (Reggia et al. 1993), 8 states, the 6-cell loop
    reggia-2   the same paper's 5-cell loop (the book's "size 5")
    sdsr       Sayama (1998, 1999), 9 states: Langton's loop plus the
               dissolving state 8 (structural dissolution)
    evoloop    Sayama (1999), 9 states: the SDSR loop made evolvable

The tables are the @TABLE sections of Golly's rule files, which come from
Sayama's own `loops.java`; they are reproduced verbatim at the end of this
module and each one is checked against the published replication period in
the tests.

Two readings of "what is a molecule" are exported (book 10.7.2), as `mode`:

- **macro** (the default): an emergent loop is a molecule. A loop is a
  connected group of at least `min_loop_cells` non-quiescent cells; loops are
  followed from step to step by overlap, and two loops belong to the same
  species when they ever show the same configuration (up to rotation). The
  observed reactions are the events of the run:

      L -> 2 L        replication (the daughter is the same species)
      L -> L + L'     replication with variation (evoloop)
      L -> ∅          dissolution / death

- **micro**: a cell state is a molecule and the quiescent state is the absence
  of one; every state transition that fires is a reaction with the four
  neighbours as catalysts,

      s2 + s1 + s7 -> s1 + s1 + s7      (a cell in state 2 becomes a 1)

  which is the book's reading of movement and reaction as state transitions.

The grid is in extras["space"] and the population time series in
extras["analysis"].
"""

from __future__ import annotations

import itertools
import re
from collections import Counter

import numpy as np
from scipy import ndimage

from chemart.network import Network, Reaction, Species

#: 8-connectivity: a loop is one connected blob of non-quiescent cells.
_STRUCT8 = np.ones((3, 3), dtype=int)

#: Replication period of each ancestor (steps until a complete copy of it stands
#: on the lattice). The published tables quote 15 for "the Chou-Reggia loop",
#: which is the 5-cell reggia-2; the 6-cell reggia-1 takes 13.
PERIODS = {"langton": 151, "byl": 25, "reggia-1": 13, "reggia-2": 15,
           "sdsr": 151, "evoloop": 363}

#: Meaning of the states in the Langton family (Sayama's legend).
STATE_MEANING = {0: "background (quiescent: no molecule)", 1: "core", 2: "sheath",
                 3: "signal", 4: "signal (left turn)", 5: "signal", 6: "signal",
                 7: "signal (straight growth)", 8: "dissolver (SDSR, evoloop)"}


# --------------------------------------------------------------------------
# Golly rule tables and RLE patterns
# --------------------------------------------------------------------------
def parse_table(text: str):
    """(n_states, symmetries, variables, transitions) from a Golly @TABLE body."""
    n_states = symmetries = neighbourhood = None
    variables: dict[str, list[int]] = {}
    transitions: list[list[str]] = []
    for raw in text.splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line.startswith("n_states:"):
            n_states = int(line.split(":")[1])
        elif line.startswith("neighborhood:"):
            neighbourhood = line.split(":")[1].strip()
        elif line.startswith("symmetries:"):
            symmetries = line.split(":")[1].strip()
        elif line.startswith("var "):
            name, values = line[4:].split("=", 1)
            variables[name.strip()] = [int(v) for v in values.strip(" {}").split(",")]
        else:
            toks = [t.strip() for t in line.split(",")] if "," in line else list(line)
            if len(toks) != 6:
                raise ValueError(f"transition {line!r} must have 6 entries (C,N,E,S,W,C')")
            transitions.append(toks)
    if neighbourhood != "vonNeumann" or symmetries != "rotate4":
        raise ValueError(f"expected a vonNeumann/rotate4 table, got {neighbourhood}/{symmetries}")
    return n_states, symmetries, variables, transitions


def dense_table(n_states: int, variables: dict[str, list[int]],
                transitions: list[list[str]]) -> np.ndarray:
    """Expand a rule table into a lookup array over (C, N, E, S, W).

    Variables are bound per transition, every transition stands for its four
    rotations, and the first matching transition wins (Golly's semantics; this
    expansion reproduces the compiled @TREE of every rule file exactly). -1
    marks a neighbourhood no transition covers: the cell then keeps its state.
    """
    table = np.full(n_states ** 5, -1, dtype=np.int8)
    for toks in transitions:
        names = [t for t in toks if not t.isdigit()]
        unknown = [t for t in names if t not in variables]
        if unknown:
            raise ValueError(f"transition {','.join(toks)} uses undeclared variable {unknown[0]!r}")
        unique = list(dict.fromkeys(names))
        for combination in itertools.product(*[variables[u] for u in unique]):
            bound = dict(zip(unique, combination))
            c, n, e, s, w, out = [int(t) if t.isdigit() else bound[t] for t in toks]
            for _ in range(4):                       # rotate4
                index = ((((c * n_states + n) * n_states + e) * n_states + s) * n_states + w)
                if table[index] < 0:
                    table[index] = out
                n, e, s, w = e, s, w, n
    return table


_CACHE: dict[str, tuple[int, np.ndarray, int]] = {}


def rule(name: str) -> tuple[int, np.ndarray, int]:
    """(n_states, lookup table, number of transitions) of a named rule, cached."""
    if name not in _CACHE:
        n_states, _, variables, transitions = parse_table(RULE_TABLES[name])
        _CACHE[name] = (n_states, dense_table(n_states, variables, transitions), len(transitions))
    return _CACHE[name]


def parse_rle(text: str) -> np.ndarray:
    """A Golly RLE pattern as an array of states ('.'/'b' = 0, 'A' = 1, ...)."""
    header, data = None, ""
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if header is None and line.startswith("x"):
            values = [int(v) for v in re.findall(r"=\s*(\d+)", line)]
            header = (values[1], values[0])
            continue
        data += line
    if header is None:
        raise ValueError("RLE pattern needs an 'x = ..., y = ...' header line")
    grid = np.zeros(header, dtype=np.int8)
    x = y = 0
    count = ""
    for ch in data:
        if ch.isdigit():
            count += ch
            continue
        n, count = (int(count) if count else 1), ""
        if ch == "!":
            break
        if ch == "$":
            y, x = y + n, 0
        elif ch in ".b":
            x += n
        elif "A" <= ch <= "X":
            if y >= grid.shape[0] or x + n > grid.shape[1]:
                raise ValueError("RLE pattern does not fit its own header")
            grid[y, x:x + n] = ord(ch) - ord("A") + 1
            x += n
        else:
            raise ValueError(f"unexpected character {ch!r} in the RLE pattern")
    return grid


# --------------------------------------------------------------------------
# the automaton
# --------------------------------------------------------------------------
def neighbours(grid: np.ndarray, periodic: bool):
    """(N, E, S, W) of every cell; outside a non-periodic lattice is quiescent."""
    north, south = np.roll(grid, 1, 0), np.roll(grid, -1, 0)
    west, east = np.roll(grid, 1, 1), np.roll(grid, -1, 1)
    if not periodic:
        north[0], south[-1], west[:, 0], east[:, -1] = 0, 0, 0, 0
    return north, east, south, west


def step(grid: np.ndarray, table: np.ndarray, n_states: int, periodic: bool = True) -> np.ndarray:
    """One synchronous update of the whole lattice."""
    north, east, south, west = neighbours(grid, periodic)
    index = grid.astype(np.int32)
    for side in (north, east, south, west):
        index = index * n_states + side
    out = table[index]
    return np.where(out < 0, grid, out).astype(np.int8)


def evolve(grid: np.ndarray, table: np.ndarray, n_states: int, steps: int,
           periodic: bool = True):
    """Yield the lattice after each of `steps` updates."""
    for _ in range(steps):
        grid = step(grid, table, n_states, periodic)
        yield grid


def rotations(pattern: np.ndarray) -> list[np.ndarray]:
    """The four rotations of a pattern, without duplicates."""
    out, seen = [], set()
    for k in range(4):
        r = np.rot90(pattern, k)
        key = (r.shape, r.tobytes())
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def occurrences(grid: np.ndarray, pattern: np.ndarray) -> list[tuple[int, int]]:
    """Top-left positions where `pattern` occurs in `grid` exactly."""
    h, w = pattern.shape
    if h > grid.shape[0] or w > grid.shape[1]:
        return []
    windows = np.lib.stride_tricks.sliding_window_view(grid, (h, w))
    hits = np.all(windows == pattern, axis=(2, 3))
    return [(int(y), int(x)) for y, x in zip(*np.nonzero(hits))]


def copies(grid: np.ndarray, pattern: np.ndarray, home: tuple[int, int] | None = None):
    """Every place where the ancestor (in any rotation) sits, other than `home`."""
    out = []
    for r in rotations(pattern):
        for pos in occurrences(grid, r):
            if home is None or (pos, r.shape) != (home, pattern.shape):
                out.append((pos, r.shape))
    return out


# --------------------------------------------------------------------------
# loops as molecules (the macro reading)
# --------------------------------------------------------------------------
def canonical(grid: np.ndarray, labels: np.ndarray, index: int) -> tuple:
    """A loop's configuration, made independent of its position and rotation."""
    ys, xs = np.nonzero(labels == index)
    box = np.where(labels[ys.min():ys.max() + 1, xs.min():xs.max() + 1] == index,
                   grid[ys.min():ys.max() + 1, xs.min():xs.max() + 1], 0).astype(np.uint8)
    return min((r.shape, r.tobytes()) for r in rotations(box))


def label_name(index: int) -> str:
    """Base-26 label: 0 -> aaa, 1 -> aab (as in the Tierra genebank)."""
    out = ""
    for _ in range(3):
        index, r = divmod(index, 26)
        out = chr(97 + r) + out
    return out if index == 0 else f"{out}{index}"


class _Species:
    """One loop species, registered by the configuration a loop was born with."""

    def __init__(self, sid: int, step_seen: int, cells: int, shape: tuple, pattern: bytes):
        self.id = sid
        self.parent = sid                 # union-find
        self.step = step_seen
        self.cells = cells
        self.shape = shape
        self.pattern = pattern
        self.name = ""


class Colony:
    """Follows the loops of a run and records the events between them."""

    def __init__(self, min_cells: int):
        self.min_cells = min_cells
        self.species: dict[int, _Species] = {}
        self.by_pattern: dict[tuple, int] = {}
        self.loops: dict[int, dict] = {}
        self.owner: np.ndarray | None = None
        self.ancestor_key: tuple | None = None    # the seed loop's configuration
        self.initial_loops: set[int] = set()
        self.ancestor_seen: dict[int, int] = {}   # loop -> step it showed that configuration
        self.next_loop = 1
        self.events: dict[tuple, list] = {}
        self.births: list[dict] = []
        self.history: list[dict] = []
        self.deaths = 0

    # -- species union-find --------------------------------------------
    def root(self, sid: int) -> int:
        while self.species[sid].parent != sid:
            self.species[sid].parent = self.species[self.species[sid].parent].parent
            sid = self.species[sid].parent
        return sid

    def merge(self, a: int, b: int) -> int:
        a, b = self.root(a), self.root(b)
        if a == b:
            return a
        old, new = (a, b) if self.species[a].step <= self.species[b].step else (b, a)
        self.species[new].parent = old
        return old

    def new_species(self, step_seen: int, cells: int, key: tuple) -> int:
        sid = len(self.species) + 1
        self.species[sid] = _Species(sid, step_seen, cells, key[0], key[1])
        self.by_pattern[key] = sid
        return sid

    # -- the run --------------------------------------------------------
    def observe(self, grid: np.ndarray, t: int, record: bool = True) -> None:
        labels, count = ndimage.label(grid > 0, structure=_STRUCT8)
        sizes = np.bincount(labels.ravel(), minlength=count + 1)
        blobs = [i for i in range(1, count + 1) if sizes[i] >= self.min_cells]
        if self.owner is None:
            self.owner = np.zeros(grid.shape, dtype=np.int32)

        claims: dict[int, list[tuple[int, int]]] = {}
        orphans: list[int] = []
        for i in blobs:
            overlap = Counter(self.owner[labels == i].tolist())
            overlap.pop(0, None)
            if overlap:
                best = max(overlap, key=lambda k: (overlap[k], -k))
                claims.setdefault(best, []).append((overlap[best], i))
            else:
                orphans.append(i)

        alive: dict[int, int] = {}            # loop id -> blob
        newborn: list[tuple[int, int | None]] = []
        for loop_id, claimed in claims.items():
            claimed.sort(reverse=True)
            alive[loop_id] = claimed[0][1]
            newborn += [(blob, loop_id) for _, blob in claimed[1:]]
        newborn += [(blob, None) for blob in orphans]

        for loop_id in [i for i in self.loops if i not in alive]:
            self.record([self.loops[loop_id]["species"]], [])
            self.deaths += 1
            del self.loops[loop_id]

        centres = {i: ndimage.center_of_mass(labels == blob) for i, blob in alive.items()}
        for blob, parent in newborn:
            if parent is None and centres:
                where = ndimage.center_of_mass(labels == blob)
                parent = min(centres, key=lambda i: (centres[i][0] - where[0]) ** 2
                             + (centres[i][1] - where[1]) ** 2)
            loop_id = self.next_loop
            self.next_loop += 1
            self.loops[loop_id] = {"species": None, "born": t, "parent": parent, "matched": None}
            alive[loop_id] = blob

        owner = np.zeros(grid.shape, dtype=np.int32)
        for loop_id, blob in alive.items():
            mask = labels == blob
            owner[mask] = loop_id
            loop = self.loops[loop_id]
            key = canonical(grid, labels, blob)
            if t and key == self.ancestor_key and loop_id not in self.ancestor_seen:
                self.ancestor_seen[loop_id] = t
            known = self.by_pattern.get(key)
            if loop["species"] is None:
                # a loop is registered by the configuration it is born with
                if known is None:
                    loop["species"] = self.new_species(t, int(sizes[blob]), key)
                else:
                    loop["species"] = self.root(known)
                    loop["matched"] = t
            elif known is not None and self.root(known) != self.root(loop["species"]):
                if loop["matched"] is None:
                    loop["matched"] = t
                loop["species"] = self.merge(known, loop["species"])
        self.owner = owner

        for blob, parent in newborn:
            if not record:                       # the ancestors are not an event
                continue
            loop_id = next(i for i, b in alive.items() if b == blob)
            child = self.loops[loop_id]
            if parent is None:
                self.record([], [child["species"]])
            else:
                mother = self.loops[parent]["species"]
                self.record([mother], [mother, child["species"]])
                self.births.append({"step": t, "loop": loop_id, "parent_loop": parent,
                                    "mother": mother, "daughter": child["species"]})

        self.history.append({
            "step": t, "loops": len(alive),
            "cells": sorted(int(sizes[b]) for b in alive.values()),
            "species": Counter(self.root(self.loops[i]["species"]) for i in alive),
        })

    def record(self, reactants: list[int], products: list[int]) -> None:
        key = (tuple(sorted(reactants)), tuple(sorted(products)))
        entry = self.events.get(key)
        if entry is None:
            self.events[key] = entry = [list(reactants), list(products), 0]
        entry[2] += 1

    # -- the network -----------------------------------------------------
    def name_species(self) -> dict[int, str]:
        roots = sorted({self.root(s) for s in self.species}, key=lambda s: (self.species[s].step, s))
        for i, sid in enumerate(roots):
            self.species[sid].name = f"L{self.species[sid].cells:03d}{label_name(i)}"
        return {s: self.species[self.root(s)].name for s in self.species}


def _structure(sp: _Species) -> str:
    rows = np.frombuffer(sp.pattern, dtype=np.uint8).reshape(sp.shape)
    return "/".join("".join(str(v) for v in row) for row in rows)


def _grid_rows(grid: np.ndarray) -> list[str]:
    return ["".join(str(v) for v in row) for row in grid]


# --------------------------------------------------------------------------
# the chemistry
# --------------------------------------------------------------------------
def ancestor(p) -> np.ndarray:
    """The seed pattern: the published ancestor of the rule, or the user's RLE."""
    text = p.ancestor_pattern.strip() or ANCESTORS[p.rule]
    return parse_rle(text)


def start(p, rng) -> np.ndarray:
    """The initial lattice, with `ancestors` copies of the seed pattern."""
    pattern = ancestor(p)
    h, w = pattern.shape
    if h > p.grid or w > p.grid:
        raise ValueError(f"grid {p.grid} is smaller than the {h}x{w} ancestor of rule {p.rule!r}")
    grid = np.zeros((p.grid, p.grid), dtype=np.int8)
    if p.ancestors == 1:
        y, x = (p.grid - h) // 2, (p.grid - w) // 2
        grid[y:y + h, x:x + w] = pattern
        return grid
    placed = []
    for _ in range(p.ancestors):
        for _try in range(200):
            y = int(rng.integers(0, p.grid - h + 1))
            x = int(rng.integers(0, p.grid - w + 1))
            if all(abs(y - b) > h or abs(x - a) > w for b, a in placed):
                placed.append((y, x))
                grid[y:y + h, x:x + w] = pattern
                break
        else:
            raise ValueError(f"cannot place {p.ancestors} ancestors of size {h}x{w} "
                             f"on a {p.grid}x{p.grid} lattice; use a larger grid")
    return grid


def generate(p, rng):
    if p.rule not in RULE_TABLES:
        raise ValueError(f"unknown rule {p.rule!r}; known: {', '.join(RULE_TABLES)}")
    n_states, table, n_transitions = rule(p.rule)
    grid = start(p, rng)
    periodic = p.boundary == "periodic"
    space = {
        "dimensions": 2, "shape": [p.grid, p.grid], "neighbourhood": "von-neumann",
        "boundary": p.boundary, "states": n_states, "rule": p.rule,
        "transitions": n_transitions,
        "state_meaning": {str(s): STATE_MEANING.get(s, "signal") for s in range(n_states)},
        "ancestor": _grid_rows(ancestor(p)),
    }
    if p.mode == "micro":
        return _micro(p, grid, table, n_states, periodic, space)
    return _macro(p, grid, table, n_states, periodic, space)


def _macro(p, grid, table, n_states, periodic, space):
    colony = Colony(p.min_loop_cells)
    colony.observe(grid, 0, record=False)
    colony.initial_loops = set(colony.loops)
    colony.ancestor_key = next(iter(colony.by_pattern), None)
    initial = Counter(colony.root(colony.loops[i]["species"]) for i in colony.loops)
    if not initial:
        raise ValueError(f"the ancestor of rule {p.rule!r} has fewer than min_loop_cells="
                         f"{p.min_loop_cells} cells; lower min_loop_cells")
    for t in range(1, p.steps + 1):
        grid = step(grid, table, n_states, periodic)
        if t % p.track_every == 0 or t == p.steps:
            colony.observe(grid, t)

    names = colony.name_species()
    used = sorted({names[s] for s in colony.species}, key=lambda n: (int(n[1:4]), n))
    species = []
    for name in used:
        sid = next(s for s in colony.species if names[s] == name and colony.root(s) == s)
        species.append(Species(name, structure=_structure(colony.species[sid])))
    merged: dict[tuple, list] = {}          # species merge after the fact: aggregate by name
    for lhs, rhs, count in colony.events.values():
        left, right = Counter(names[s] for s in lhs), Counter(names[s] for s in rhs)
        key = (tuple(sorted(left.items())), tuple(sorted(right.items())))
        entry = merged.setdefault(key, [dict(left), dict(right), 0])
        entry[2] += count
    reactions = [Reaction(left, right, count=count) for left, right, count in merged.values()]
    final = colony.history[-1]
    analysis = {
        "steps": p.steps,
        "tracked_every": p.track_every,
        "births": len(colony.births),
        "deaths": colony.deaths,
        "species_seen": len(used),
        "loops_final": final["loops"],
        "replications": [
            {"step": b["step"], "identified": colony.loops.get(b["loop"], {}).get("matched"),
             "copy": colony.ancestor_seen.get(b["loop"]),
             "mother": names[b["mother"]], "daughter": names[b["daughter"]]}
            for b in colony.births[:200]
        ],
        # the published replication period: when a new loop first shows the exact
        # configuration of the ancestor, and when the ancestor itself shows it again
        "ancestor_copy_steps": sorted(t for i, t in colony.ancestor_seen.items()
                                      if i not in colony.initial_loops)[:20],
        "ancestor_recurrence": sorted(t for i, t in colony.ancestor_seen.items()
                                      if i in colony.initial_loops)[:20],
        "population": [
            {"step": h["step"], "loops": h["loops"], "cells": h["cells"],
             "by_species": {names[s]: c for s, c in h["species"].items()}}
            for h in colony.history[:: max(1, len(colony.history) // 50)]
        ],
        "final_population": {names[s]: c for s, c in final["species"].items()},
    }
    space["final"] = _grid_rows(grid)
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={names[s]: float(c) for s, c in initial.items()},
        extras={"space": space, "analysis": analysis,
                "events": "L -> 2 L replication; L -> L + L' replication with variation; "
                          "L -> ∅ dissolution or death"},
    )


def _micro(p, grid, table, n_states, periodic, space):
    """Every state transition that fires, as a reaction with the neighbours as catalysts."""
    initial_cells = Counter(grid[grid > 0].tolist())
    fired: Counter = Counter()
    states = Counter()
    history = []
    for t in range(1, p.steps + 1):
        north, east, south, west = neighbours(grid, periodic)
        after = step(grid, table, n_states, periodic)
        changed = after != grid
        if changed.any():
            hood = np.sort(np.stack([north[changed], east[changed], south[changed],
                                     west[changed]]).astype(np.int64), axis=0)
            key = grid[changed].astype(np.int64)
            for row in hood:
                key = key * n_states + row
            key = key * n_states + after[changed]
            values, counts = np.unique(key, return_counts=True)
            fired.update(dict(zip(values.tolist(), counts.tolist())))
        grid = after
        if t % max(1, p.steps // 50) == 0 or t == p.steps:
            cells = Counter(grid[grid > 0].tolist())
            history.append({"step": t, "cells": {f"s{k}": int(v) for k, v in sorted(cells.items())}})
        states.update(grid[grid > 0].tolist())

    reactions = []
    seen = set()
    for key, count in sorted(fired.items()):
        rest, out = divmod(key, n_states)
        hood = []
        for _ in range(4):
            rest, value = divmod(rest, n_states)
            hood.append(value)
        before = rest
        lhs = [before] + hood
        rhs = [out] + hood
        seen.update(v for v in lhs + rhs if v)
        reactions.append(Reaction(
            dict(Counter(f"s{v}" for v in lhs if v)),
            dict(Counter(f"s{v}" for v in rhs if v)),
            count=int(count)))
    species = [Species(f"s{v}", structure=STATE_MEANING.get(v, "signal"))
               for v in sorted(seen)]
    space["final"] = _grid_rows(grid)
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={f"s{v}": float(c) for v, c in sorted(initial_cells.items()) if v},
        extras={"space": space,
                "analysis": {"steps": p.steps, "transitions_fired": len(reactions),
                             "firings": int(sum(fired.values())), "cells": history},
                "events": "a cell in state C with neighbours N,E,S,W becomes C': the "
                          "neighbours are catalysts, the quiescent state is the absence "
                          "of a molecule"},
    )


# --------------------------------------------------------------------------
# The published transition tables and ancestors follow (generated verbatim
# from Golly's rule and pattern files; see the catalog entry's `sources`).
# --------------------------------------------------------------------------
#: Ancestor patterns, in Golly RLE (Patterns/Loops/*.rle).
ANCESTORS: dict[str, str] = {
    "langton": """\
x = 15, y = 10, rule = Langtons-Loops
.8B$BAG.AD.ADB$B.6B.B$BGB4.BAB$BAB4.BAB$B.B4.BAB$BGB4.BAB$BA6BA5B$B.G
A.GA.G5AB$.13B!""",
    "byl": """\
x = 4, y = 4, rule = Byl-Loop
.2B$BCAB$BCDB$.BE!""",
    "reggia-1": """\
x = 4, y = 2, rule = Chou-Reggia-1
2A$CD2A!""",
    "reggia-2": """\
x = 3, y = 2, rule = Chou-Reggia-2
2A$CDA!""",
    "sdsr": """\
x = 15, y = 10, rule = SDSR-Loop
.8B$BAG.AD.ADB$B.6B.B$BGB4.BAB$BAB4.BAB$B.B4.BAB$BGB4.BAB$BA6BA5B$B.G
A.GA.G5AB$.13B!""",
    "evoloop": """\
x = 17, y = 17, rule = Evoloop
.15B$BG.AG.AG.AG.4AB$BA13BAB$B.B11.BAB$BGB11.BAB$BAB11.BAB$B.B11.BAB$
BGB11.BAB$BAB11.BAB$B.B11.BAB$BGB11.BGB$BAB11.B.B$B.B11.BAB$BGB11.BGB
$BA13B.B$B.GA.GA.GA.DA.DAB$.14BE!""",
}

#: Transition tables, verbatim from the @TABLE section of Golly's rule files
#: (comments and blank lines removed); see `sources` in the catalog entry.
RULE_TABLES: dict[str, str] = {
    "langton": """\
n_states:8
neighborhood:vonNeumann
symmetries:rotate4
000000
000012
000020
000030
000050
000063
000071
000112
000122
000132
000212
000220
000230
000262
000272
000320
000525
000622
000722
001022
001120
002020
002030
002050
002125
002220
002322
005222
012321
012421
012525
012621
012721
012751
014221
014321
014421
014721
016251
017221
017255
017521
017621
017721
025271
100011
100061
100077
100111
100121
100211
100244
100277
100511
101011
101111
101244
101277
102026
102121
102211
102244
102263
102277
102327
102424
102626
102644
102677
102710
102727
105427
111121
111221
111244
111251
111261
111277
111522
112121
112221
112244
112251
112277
112321
112424
112621
112727
113221
122244
122277
122434
122547
123244
123277
124255
124267
125275
200012
200022
200042
200071
200122
200152
200212
200222
200232
200242
200250
200262
200272
200326
200423
200517
200522
200575
200722
201022
201122
201222
201422
201722
202022
202032
202052
202073
202122
202152
202212
202222
202272
202321
202422
202452
202520
202552
202622
202722
203122
203216
203226
203422
204222
205122
205212
205222
205521
205725
206222
206722
207122
207222
207422
207722
211222
211261
212222
212242
212262
212272
214222
215222
216222
217222
222272
222442
222462
222762
222772
300013
300022
300041
300076
300123
300421
300622
301021
301220
302511
401120
401220
401250
402120
402221
402326
402520
403221
500022
500215
500225
500232
500272
500520
502022
502122
502152
502220
502244
502722
512122
512220
512422
512722
600011
600021
602120
612125
612131
612225
700077
701120
701220
701250
702120
702221
702251
702321
702525
702720""",
    "byl": """\
n_states:6
neighborhood:vonNeumann
symmetries:rotate4
000000
000010
000031
000420
000110
000120
000311
000330
000320
000200
000242
000233
004330
005120
001400
001120
001305
001310
001320
003020
003400
003420
003120
002040
002050
002010
002020
002420
002500
002520
002200
002250
002220
040000
050000
050500
010052
010022
010400
010100
020055
400233
405235
401233
442024
452020
415233
411233
412024
412533
432024
421433
422313
501302
502230
540022
542002
512024
530025
520025
520442
523242
522020
100000
100010
100033
100330
101233
103401
103244
111244
113244
112404
133244
121351
123444
123543
122414
122434
300100
300300
300211
300233
304233
301100
301211
303211
303233
302233
344233
343233
351202
353215
314233
311233
313251
313211
312211
335223
333211
320533
324433
325415
321433
321511
321321
323411
200000
200042
200032
200022
200442
200515
200112
200122
200342
200332
200242
200212
201502
203202
202302
244022
245022
242042
242022
254202
255042
252025
210042
214022
215022
212055
212022
230052
234022
235002
235022
232042
232022
220042
220020
220533
221552""",
    "reggia-1": """\
n_states:8
neighborhood:vonNeumann
symmetries:rotate4
000000
000200
000110
004420
006301
001447
003030
046000
060400
062102
051060
010420
010360
012152
032000
401033
410011
431206
640406
660630
616600
636000
000440
000260
000330
004512
005103
001540
003100
045002
060600
061110
051710
010506
010320
030400
076010
460313
410633
600000
640106
610006
616100
636100
000420
000231
000705
006030
002020
001210
040000
045600
060100
050006
020320
010140
016000
036000
400035
466633
410100
606100
660060
610106
613001
676600
000600
000100
004040
006440
001010
001110
040160
041040
060166
050306
024200
010154
016100
036606
400313
423106
411033
606313
660460
614004
630160
500000
000540
000152
004010
006110
001030
003000
044400
060066
062030
054040
023010
010340
012000
036300
405033
410061
431033
640006
660615
616000
630360
500044
500011
230000
100111
101044
146004
143321
166611
121121
117044
300111
340066
366611
374017
200100
100044
104414
101601
141403
160414
161101
112244
130414
304011
340611
314001
703010
203010
100011
104674
101301
141424
160111
161301
111044
133001
301471
340261
314676
703330
240400
100441
106101
103011
141141
164104
120414
113011
133011
301303
340211
313023
733630
210102
100414
105011
103103
141121
166644
124104
113601
300101
307141
344011
331001""",
    "reggia-2": """\
n_states:8
neighborhood:vonNeumann
symmetries:rotate4
000000
000440
000547
000100
000110
000330
004040
004445
004100
001040
001010
001740
003000
003010
003030
007040
007030
007104
007110
040000
040070
047100
050007
017000
070000
070077
071010
400101
400313
401033
407103
411033
431033
500033
503330
100045
100011
100414
101044
101301
103011
107131
140414
141044
111044
133011
177711
304011
305011
305141
307141
344011
345011
354010
314001
377711
700000
700337
707100
707141
707110
717000
730007
770070
770711""",
    "sdsr": """\
n_states:9
neighborhood:vonNeumann
symmetries:rotate4
var a={0,1,2}
var b={0,1}
var c={0,1}
var d={0,1}
var e={2,3}
var f={0,1,2,3,4,5,6,7}
var g={0,1,2,3,4,6,7}
var h={4,6,7}
var i={6,7}
var j={0,1,2,3,4,5,6,7,8}
var k={0,1,2,3,4,5,6,7,8}
var l={2,3,4,5,6,7}
var m={0,1,2,3,5}
var n={0,1,8}
var o={1,4,6,7}
var p={0,1,3,5}
var q={1,2,4,6,7}
var r={3,5}
var s={2,3,4,5,6}
var t={0,3,5}
var u={1,3,5}
var v={0,5}
var w={2,3,5}
var x={0,2}
var y={1,4,7}
var z={0,3,5}
var A={0,3,5}
var B={0,2,3,4,5}
var C={0,2,4,6,7}
var D={0,1,2,4,6,7}
var E={0,1,2,4,6,7}
var F={1,3,4,6,7}
var G={0,1,2,3,4,5}
var H={0,1,2,3,4}
var I={1,2,4}
var J={0,3,5,6}
var K={0,1,3,4,5,6}
var L={1,2,4,6}
var M={0,3}
var N={0,2,3,5}
var O={1,3,4,5,6,7}
var P={0,1,3,4}
var Q={1,2}
var R={2,4,6,7}
var S={0,1,2,3,4,5,6}
var T={0,1,3,4,5,6,7}
var U={1,5}
var V={0,1,3,4,5}
var W={1,3,5}
var X={1,3,5}
var Y={0,3,4,5,6}
var Z={0,6}
var aa={0,6}
var ab={1,3}
var ac={2,5}
var ad={2,3,5}
var ae={0,7}
var af={2,7}
var ag={1,4,6,7,8}
var ah={1,3,5}
var ai={1,2,3,4,5}
var aj={1,4,5,6,7}
var ak={1,2,3,5}
var al={1,2,3,5}
var am={4,7}
var an={0,5,6}
var ao={0,5,6}
var ap={0,5,6}
var aq={0,5,6}
var ar={0,3}
var as={0,3}
var at={3,7}
var au={2,3,5,8}
var av={0,1,2,3,4,5,6,7,8}
var aw={0,1,4,5,6,7}
var ax={0,1,4,5,6,7}
var ay={1,4,6,7}
var az={0,1,2,3,4,5,6,7}
var aA={2,8}
var aB={4,5,6,7}
var aC={4,5,6}
var aD={0,1,5}
var aE={0,1,5}
var aF={1,4,5,6,7}
var aG={0,1,2,4,5,6,7}
var aH={5,6,7}
var aI={0,2,3,4,5,6,7}
var aJ={0,4,6,7}
var aK={2,3,5,7}
var aL={1,2,3,4,6,7}
var aM={1,2,3,4,5,6,7}
var aN={3,4,6,7,8}
var aO={1,3,4,5,6,7}
var aP={2,4,6,7}
var aQ={4,5}
var aR={0,1,3,4,6,7}
var aS={4,6,7,8}
var aT={1,2,4,6,7}
var aU={4,6,8}
var aV={1,2,3,4,5,6,7}
var aW={1,2,3,4,5,6,7}
var aX={1,2,3,4,5,6,7}
var aY={0,3,7}
var aZ={0,1,2,3}
var ba={0,3,4,5,6,7}
var bb={5,8}
var bc={0,1,7}
var bd={0,1,3,7}
var be={1,2,7}
var bf={1,3,6,7}
var bg={0,1,2,3,7}
var bh={3,6}
var bi={5,6}
var bj={2,3,4,5,6,7}
var bk={2,3,4,5,6,7}
var bl={2,3,4,5,6,7}
var bm={6,8}
var bn={6,7,8}
var bo={3,4,6,7}
var bp={1,6}
var bq={0,1,4,6,7}
var br={0,1,2,3,4,5,6,7}
var bs={0,1,2,3,4,5,6,7}
0,0,0,a,1,2
0,0,0,0,6,3
0,b,c,d,7,1
0,0,0,1,e,2
0,f,g,1,h,1
0,0,0,2,i,2
b,j,k,l,8,8
0,m,f,h,1,1
0,0,0,5,2,5
0,0,0,i,2,2
b,f,n,8,l,8
0,m,1,f,o,1
0,0,1,0,2,2
0,p,1,q,r,1
0,m,1,s,2,1
b,c,l,d,8,8
0,t,2,1,u,1
0,v,2,1,2,5
0,0,2,w,1,1
0,0,2,3,2,2
0,0,r,1,q,1
0,0,r,2,1,1
0,0,5,2,2,2
x,1,q,3,f,1
0,1,7,2,5,5
0,2,5,2,7,1
y,t,z,A,B,8
1,C,D,E,7,7
F,A,t,f,r,8
1,G,H,I,4,4
1,J,K,L,6,6
1,M,p,N,w,8
O,A,t,3,f,8
1,m,P,4,Q,4
F,0,A,5,R,8
1,S,K,6,L,6
1,f,T,7,q,7
U,M,u,m,w,8
1,p,I,A,4,4
1,A,L,V,6,6
u,A,W,w,X,8
1,p,q,f,7,7
1,p,L,7,Y,7
1,Z,2,aa,2,6
1,A,w,U,W,8
ab,N,w,ac,ad,8
1,0,2,2,6,3
1,ae,af,3,2,7
1,B,2,6,W,6
1,0,2,6,4,4
ag,0,2,7,1,0
O,A,r,v,D,8
1,b,s,L,7,7
1,p,5,I,4,4
1,A,5,4,1,4
1,x,5,4,2,7
1,0,7,r,L,7
O,W,X,u,ah,8
1,1,Q,2,7,7
ab,m,Q,5,x,8
1,ai,2,6,4,6
aj,ak,5,ad,al,8
1,W,5,4,2,4
1,2,ad,2,6,6
1,2,am,2,5,5
1,2,4,2,6,7
ad,an,ao,ap,aq,8
2,M,ar,as,at,1
au,j,k,av,8,0
ad,aw,ax,o,ay,8
2,f,az,O,3,1
aA,x,0,2,5,0
2,aw,x,3,aj,1
2,a,0,3,2,6
2,0,0,4,2,3
ad,aw,ax,aB,aC,8
2,0,0,5,1,7
2,0,b,5,aC,8
2,x,0,5,7,5
ad,aD,aj,aE,aF,8
2,0,O,x,3,1
2,0,2,0,7,3
2,aG,l,2,3,1
af,0,2,3,2,1
2,0,3,2,aB,1
2,0,3,aH,2,1
2,0,5,5,2,1
2,1,1,2,6,1
3,0,0,Z,2,2
3,aI,az,f,r,8
3,x,0,0,4,1
3,0,0,0,7,6
3,D,aJ,aK,aL,8
3,v,1,0,2,1
3,az,aM,F,aL,8
aN,0,1,2,2,0
r,T,O,az,aO,8
3,0,R,0,aP,8
aQ,A,0,v,aR,8
h,av,j,k,8,1
aS,0,aw,q,aT,0
aS,0,ay,A,aT,0
aU,0,aT,q,O,0
aU,0,2,aw,aT,0
aS,0,2,aM,ay,0
4,0,e,2,2,1
4,0,2,3,2,6
aS,0,e,ay,aT,0
aS,0,3,2,ay,0
am,aM,aV,aW,aX,8
5,aY,0,0,2,2
5,aZ,T,F,az,8
5,az,ba,e,aC,8
bb,0,0,5,2,0
5,b,2,bc,2,2
5,0,e,bd,h,8
5,v,e,be,bf,8
5,0,2,1,5,2
bb,b,2,2,2,0
5,0,2,2,4,4
5,0,2,af,5,8
5,bg,2,bh,2,8
5,1,2,4,2,2
bi,l,bj,bk,bl,8
6,0,0,0,aJ,8
6,0,0,0,Q,1
i,0,0,5,1,8
bm,0,2,e,2,0
bn,0,3,2,2,0
6,O,aF,aM,bj,8
6,1,2,Q,2,5
6,1,2,1,3,1
6,1,e,bo,2,8
6,1,3,2,2,8
7,0,0,0,bp,8
7,0,ay,aT,r,0
7,0,2,bq,2,0
7,0,2,ay,r,0
7,0,2,2,ac,1
7,0,2,2,3,0
7,0,2,5,2,5
8,az,f,br,bs,0""",
    "evoloop": """\
n_states:9
neighborhood:vonNeumann
symmetries:rotate4
var a={0,2,5}
var b={0,1,2,3,4,5,6,7}
var c={0,1,2,3,4,5,6,7}
var d={1,4,6,7}
var e={1,4}
var f={0,1}
var g={0,1,2,3,4,5,6,7,8}
var h={0,1,2,3,4,5,6,7,8}
var i={2,3,4,5,6,7}
var j={0,2,3,5}
var k={0,2,3,4,5,6,7}
var l={4,6,7}
var m={2,5}
var n={0,1,8}
var o={0,3,5}
var p={2,3,5}
var q={3,5}
var r={2,4,6,7}
var s={0,1}
var t={0,1}
var u={0,1,2,3,4,7}
var v={5,6}
var w={1,6,7}
var x={0,3,5}
var y={0,3,5}
var z={1,3,5}
var A={0,1,3,5}
var B={0,1,3,5}
var C={1,3,5}
var D={1,3,5}
var E={0,1,2,3,4,5}
var F={0,1,2,4,5}
var G={1,2,4}
var H={0,1,3,4,5,6}
var I={0,1,2,3,4,5,6}
var J={1,2,4,6}
var K={0,1,2,4,5,6,7}
var L={0,1,2,3,4,5,7}
var M={1,2,4,6,7}
var N={0,3}
var O={0,1,2,3,5}
var P={0,1,2,3,4}
var Q={0,2,3,5,6}
var R={0,1,3,4,5,6,7}
var S={1,3}
var T={1,2,3,5}
var U={0,1,3,4,5}
var V={0,1,4,5,6}
var W={1,4,6}
var X={2,3,5}
var Y={0,5}
var Z={1,2}
var aa={0,1,5}
var ab={0,2}
var ac={2,3,4,5}
var ad={2,3,4,5,6}
var ae={0,3}
var af={0,3}
var ag={3,7}
var ah={2,8}
var ai={0,2,8}
var aj={0,8}
var ak={0,2,8}
var al={6,8}
var am={0,1,4,5,6,7}
var an={0,1,4,5,6,7}
var ao={1,4,6,7}
var ap={1,3,4,5,6,7}
var aq={2,3,5,8}
var ar={0,1,2,3,4,5,6,7,8}
var as={0,1,4,6,7}
var at={1,5,6}
var au={4,7}
var av={1,3,4,5,7}
var aw={0,4,5,7}
var ax={1,4,5,6,7}
var ay={0,4,5,7}
var az={1,4,5,6,7}
var aA={0,1,3,5,6,7}
var aB={1,2,3,4,5,7}
var aC={0,1,4,7}
var aD={2,4,6,7,8}
var aE={3,5,6}
var aF={2,3,5,6}
var aG={2,3,5,6}
var aH={2,3}
var aI={0,1,3,4,7}
var aJ={2,5,6,7}
var aK={1,2,3,4,6,7}
var aL={1,3,4,7}
var aM={0,2,3}
var aN={1,2,4,7}
var aO={3,6}
var aP={1,2,3,4,7}
var aQ={3,4,7}
var aR={1,2,3,4,5,6,7}
var aS={3,4,6,7,8}
var aT={0,3,6}
var aU={3,4,5,7}
var aV={2,3,6}
var aW={1,2,3,4,5,6,7}
var aX={0,8}
var aY={0,8}
var aZ={4,5,6,7}
var ba={4,6,7,8}
var bb={1,2,4,6,7}
var bc={4,5}
var bd={4,7,8}
var be={0,1,4,5,7}
var bf={4,5,7}
var bg={1,2,4,5,7}
var bh={1,2,4,5,7}
var bi={0,1,2,3,4,6,7}
var bj={2,5,6}
var bk={1,3,4,6,7}
var bl={0,4,6,7}
var bm={6,7,8}
var bn={1,3,4,5,6,7}
var bo={1,4,7}
var bp={2,3,5,6}
var bq={2,3,5,6}
var br={0,1,2,3,4,5,6,7}
var bs={0,1,2,3,4,5,6,7}
0,a,0,0,1,2
0,0,0,0,4,3
0,b,c,1,d,1
0,0,0,2,e,2
f,g,h,i,8,8
0,j,k,l,1,1
0,0,0,4,m,2
0,0,0,7,5,2
f,b,n,8,i,8
0,o,1,j,d,1
0,0,1,0,2,2
0,o,1,2,p,1
0,j,1,q,r,1
f,s,i,t,8,8
0,u,2,1,p,1
0,0,r,p,1,1
0,0,2,3,2,2
0,0,q,1,2,1
0,o,q,2,1,1
0,1,2,v,2,6
w,o,x,y,k,8
z,A,B,C,D,8
1,E,F,G,4,4
1,H,I,J,6,6
1,K,L,M,7,7
1,N,A,j,p,8
1,O,P,4,G,4
1,Q,I,6,J,6
1,b,R,7,M,7
1,N,S,f,T,8
1,O,e,o,4,4
1,U,J,o,6,6
1,V,M,I,7,7
1,I,W,7,3,7
1,o,2,N,4,4
S,j,p,O,5,8
C,j,p,2,X,8
1,0,2,3,2,4
1,f,2,5,2,7
1,f,2,5,4,3
1,f,2,7,3,5
1,Y,3,Z,4,4
C,aa,q,ab,D,8
1,o,5,4,Z,4
1,f,6,2,4,4
1,o,7,3,M,7
C,f,2,X,q,8
1,e,ac,6,2,6
C,f,A,5,ab,8
1,1,4,3,3,4
1,X,2,5,4,4
1,ad,2,7,3,7
1,2,4,3,3,3
1,2,6,2,7,6
X,0,0,0,0,8
2,N,ae,af,ag,1
ah,ai,aj,ak,al,0
X,am,an,d,ao,8
2,b,R,ap,3,1
aq,g,h,ar,8,0
2,as,ab,3,ap,1
2,0,0,3,2,4
2,0,0,4,2,3
X,am,an,v,at,8
2,ab,0,5,au,5
2,0,0,8,av,0
X,aw,ax,ay,az,8
2,aA,ap,ab,3,1
2,0,aB,0,8,0
2,aC,2,0,6,5
2,0,2,0,7,3
2,am,2,ax,3,1
2,aC,2,3,2,3
2,0,2,5,2,5
aD,0,2,6,Y,0
2,0,3,2,ax,1
2,ab,3,aE,2,1
2,2,aF,aG,3,1
2,2,3,4,5,1
3,0,0,0,aH,2
3,R,k,b,v,8
3,0,0,0,7,4
3,aI,R,aJ,aK,8
3,0,0,3,2,2
3,u,U,aL,ao,8
3,0,0,4,2,1
3,aM,b,aN,aO,8
3,0,1,0,2,1
3,f,aP,s,aQ,8
q,R,ax,aR,aK,8
aS,0,1,2,5,0
3,0,2,aT,2,8
3,aC,2,5,2,1
3,0,3,3,2,1
aU,ap,aR,aV,aW,8
4,aj,aX,aY,ai,1
aZ,o,x,y,ap,8
ba,0,am,M,bb,0
l,ar,g,h,8,1
4,o,x,2,q,8
bc,0,o,q,2,8
4,0,0,8,ap,1
ba,0,ao,o,M,0
4,0,ap,0,8,1
ba,0,M,bb,ap,0
bd,0,2,be,M,0
4,0,2,0,q,8
ba,0,2,C,ao,0
4,0,2,aH,2,1
au,0,2,6,2,6
ba,0,3,ao,M,0
ba,0,3,2,ao,0
4,0,3,2,2,1
bf,bg,bh,aR,aW,8
5,b,bi,aF,bj,8
5,0,0,2,3,2
5,aI,bk,bl,bi,8
5,0,1,2,1,8
5,0,2,0,m,2
5,0,2,1,5,2
5,0,2,au,5,8
5,0,3,1,2,0
6,0,2,aC,2,2
al,0,2,aF,2,0
bm,0,3,2,2,0
6,ap,bn,aR,aW,8
6,ap,2,bn,2,8
al,bo,2,2,2,0
6,aF,aG,bp,bq,8
7,0,2,2,2,1
7,0,2,3,2,0
8,b,c,br,bs,0""",
}
