"""Reflexive artificial chemistry: finite state machines composing finite state machines.

Catalog id: reflexive-ac. Book 9.8; Salzberg, "A graph-based reflexive artificial
chemistry", BioSystems 87:1-12 (2007) [737]. That paper was not accessible; the
chemistry is reconstructed from its abstract, from Salzberg, Sayama & Ikegami, "A
tangled hierarchy of graph-constructing graphs" (ALife IX, 2004), and from
Salzberg & Sayama, "Reflexive composition of elementary state machines"
(arXiv:2505.07186, 2025), which restates the formulation of the 2006-2007 papers.

A molecule is a Mealy machine drawn as a labelled directed graph: every node
(state) has at most one outgoing link per input symbol, labelled
``input/output``, and one node is the current state (the pointer). Only the
part reachable from the pointer is kept, so a molecule is a pointed, accessible,
deterministic graph over the alphabet {0, ..., k-1}. Partial machines (missing
links) are allowed, so input data such as a tape can be written as a graph too.

The reaction is reflexive composition (arXiv:2505.07186, eqs. 1-2). A sender
machine s and a receiver machine r read one input i together: the sender's
output is piped into the receiver, and the two states advance with their
roles swapped,

    o(t)      = phi_r(phi_s(i))
    q_r(t+1)  = delta_s(i)                (now in the sender position)
    q_s(t+1)  = delta_r(phi_s(i))         (now in the receiver position)

so neither machine is only a processor or only a carrier of information. The
composite is itself a machine whose nodes are (sender state, receiver state)
pairs; its reachable part from the pair of pointers is the product graph:

    s1 + s2 -> s1 + s2 + s1 o s2          (s1 sends, s2 receives)

Composing a molecule with itself (s + s) is the graph reading its own
structure. Repeated composition gives the n-machine lattice of eqs. 3-4.
Products with more than max_states nodes are elastic.

Species id: the canonical text of the graph, nodes numbered breadth-first from
the pointer (links followed in input order), so isomorphic pointed graphs are
one species. Node j is written as its links for inputs 0..k-1, each
``destination/output`` or ``-`` when missing; nodes are separated by ``;``.
E.g. the adding machine M45 at S0 is ``0/0,1/0;1/1,0/1``. Elementary two-state
two-symbol machines can also be given by their number (Table 1 of
arXiv:2505.07186), e.g. ``M45`` or ``M45@S1`` for the pointer at S1.
"""

from __future__ import annotations

import re
from collections import Counter, deque

from chemart.expand import expand
from chemart.helpers.params import apportion
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import soup

# A machine is a tuple of nodes; a node is a tuple over input symbols of
# (destination, output) pairs or None for a missing link. Canonical machines
# have their pointer at node 0 and nodes numbered breadth-first.


# --- canonical form -----------------------------------------------------------
def canonical(nodes, start: int = 0) -> tuple:
    """The accessible part of a machine from `start`, renumbered breadth-first."""
    order = {start: 0}
    queue = deque([start])
    out = []
    while queue:
        node = nodes[queue.popleft()]
        row = []
        for link in node:
            if link is None:
                row.append(None)
                continue
            dest, sym = link
            if dest not in order:
                order[dest] = len(order)
                queue.append(dest)
            row.append((order[dest], sym))
        out.append(tuple(row))
    return tuple(out)


def to_id(m: tuple) -> str:
    return ";".join(",".join("-" if t is None else f"{t[0]}/{t[1]}" for t in node) for node in m)


def symbols(m: tuple) -> int:
    return len(m[0])


def to_text(m: tuple) -> str:
    """Readable structure: q<j>[<input>:q<dest>/<output> ...], prefixed by the elementary name if any."""
    body = " ".join(
        f"q{j}[" + " ".join(f"{i}:q{t[0]}/{t[1]}" for i, t in enumerate(node) if t is not None) + "]"
        for j, node in enumerate(m)
    )
    n = number(m) if len(m) == 2 and symbols(m) == 2 else None
    return f"M{n}@S0 {body}" if n is not None else body


# --- elementary machines (arXiv:2505.07186, Table 1) ---------------------------
def elementary(n: int) -> tuple:
    """Two-state two-symbol machine number n: bits 7..0 are delta/phi for (S0,0), (S0,1), (S1,0), (S1,1)."""
    if not isinstance(n, int) or not 0 <= n <= 255:
        raise ValueError(f"elementary machine numbers run from 0 to 255, got {n!r}")
    bits = [(n >> (7 - k)) & 1 for k in range(8)]
    return tuple(tuple((bits[4 * q + 2 * i], bits[4 * q + 2 * i + 1]) for i in range(2)) for q in range(2))


def number(nodes) -> int | None:
    """Inverse of elementary() for a complete machine with states S0 = node 0 and S1 = node 1."""
    if len(nodes) != 2 or any(len(node) != 2 or None in node for node in nodes):
        return None
    n = 0
    for q in range(2):
        for i in range(2):
            d, o = nodes[q][i]
            n = (n << 2) | (d << 1) | o
    return n


# --- notation -------------------------------------------------------------------
_ELEMENTARY = re.compile(r"^M(\d{1,3})(?:@S?([01]))?$")
_LINK = re.compile(r"^(\d+)/(\d+)$")


def parse(text: str) -> tuple:
    """Parse ``M45``, ``M45@S1``, or a machine text ``0/0,1/0;1/1,0/1`` (optionally ``@j`` for the pointer)."""
    if not isinstance(text, str):
        raise ValueError(f"a machine must be a string, got {text!r}")
    text = text.strip()
    m = _ELEMENTARY.match(text)
    if m:
        return canonical(elementary(int(m.group(1))), int(m.group(2) or 0))
    body, _, at = text.partition("@")
    rows = [r.strip() for r in body.split(";")]
    nodes = []
    for j, row in enumerate(rows):
        links = []
        for item in row.split(","):
            item = item.strip()
            if item == "-":
                links.append(None)
                continue
            lm = _LINK.match(item)
            if not lm:
                raise ValueError(f"machine {text!r}: node {j} has the link {item!r}; "
                                 "write destination/output (e.g. 1/0) or - for a missing link")
            links.append((int(lm.group(1)), int(lm.group(2))))
        nodes.append(tuple(links))
    k = len(nodes[0])
    for j, node in enumerate(nodes):
        if len(node) != k:
            raise ValueError(f"machine {text!r}: node {j} has {len(node)} links, node 0 has {k}; "
                             "every node needs one entry per input symbol")
        for link in node:
            if link and (link[0] >= len(nodes) or link[1] >= k):
                raise ValueError(f"machine {text!r}: link {link[0]}/{link[1]} of node {j} is out of range "
                                 f"({len(nodes)} nodes, {k} symbols)")
    start = 0
    if at:
        if not at.strip().lstrip("S").isdigit() or int(at.strip().lstrip("S")) >= len(nodes):
            raise ValueError(f"machine {text!r}: pointer {at!r} is not a node index")
        start = int(at.strip().lstrip("S"))
    return canonical(tuple(nodes), start)


# --- the reaction -----------------------------------------------------------------
def compose(s: tuple, r: tuple, max_states: int | None = None, same: bool | None = None) -> tuple | None:
    """Reflexive composition s o r (s sends, r receives): the reachable composite machine.

    A composite node is (tag, x, y): x is the state in the sender position, y the
    state in the receiver position, and tag says which machine is sending (0: s,
    1: r). When s and r are the same species the two machines are identical and
    the tag is dropped (the uniform lattice of arXiv:2505.07186). Returns None
    when the product would exceed max_states.
    """
    if symbols(s) != symbols(r):
        raise ValueError("machines must share the alphabet to compose")
    if same is None:
        same = s == r
    pair = (s, r)
    start = (0, 0, 0)
    index = {start: 0}
    queue = deque([start])
    out = []
    while queue:
        tag, x, y = queue.popleft()
        sender, receiver = pair[tag], pair[1 - tag]
        row = []
        for link in sender[x]:
            if link is None or receiver[y][link[1]] is None:
                row.append(None)
                continue
            dx, msg = link
            dy, o = receiver[y][msg]
            nxt = (tag if same else 1 - tag, dy, dx)
            if nxt not in index:
                if max_states is not None and len(index) >= max_states:
                    return None
                index[nxt] = len(index)
                queue.append(nxt)
            row.append((index[nxt], o))
        out.append(tuple(row))
    return tuple(out)


class Chemistry:
    """Composition rule plus the size boundary, with outcomes cached by species id."""

    def __init__(self, p):
        self.max_states = p.max_states
        self.machines: dict[str, tuple] = {}
        self.cache: dict[tuple[str, str], str | None] = {}
        self.stats = Counter()

    def add(self, m: tuple) -> str:
        sid = to_id(m)
        self.machines.setdefault(sid, m)
        return sid

    def product(self, a: str, b: str) -> str | None:
        """Species id of a o b, or None when the product exceeds max_states (elastic)."""
        key = (a, b)
        if key not in self.cache:
            m = compose(self.machines[a], self.machines[b], self.max_states, same=a == b)
            if m is None:
                self.stats["too_large"] += 1
                self.cache[key] = None
            else:
                self.cache[key] = self.add(m)
        return self.cache[key]

    def react(self, a: str, b: str):
        c = self.product(a, b)
        return None if c is None else (a, b, c)

    def species(self, ids) -> list[Species]:
        return [Species(s, structure=to_text(self.machines[s])) for s in ids]


def _random_machine(rng, states: int, k: int) -> tuple:
    nodes = tuple(
        tuple((int(rng.integers(states)), int(rng.integers(k))) for _ in range(k)) for _ in range(states)
    )
    return canonical(nodes, 0)


def _seeds(p, chem: Chemistry, rng) -> list[str]:
    if p.machines:
        if not isinstance(p.machines, list):
            raise ValueError(f"machines must be a list of machine texts, got {p.machines!r}")
        ids = [chem.add(parse(t)) for t in p.machines]
        if len({symbols(chem.machines[s]) for s in ids}) > 1:
            raise ValueError("machines must all use the same number of input symbols (links per node)")
        return ids
    return [chem.add(_random_machine(rng, p.states, p.symbols)) for _ in range(p.M)]


def generate(p, rng):
    chem = Chemistry(p)
    drawn = _seeds(p, chem, rng)
    seeds = list(dict.fromkeys(drawn))
    extras = {
        "notation": "species id: nodes numbered breadth-first from the pointer (node 0), separated by ';'; "
                    "node j lists its links for inputs 0..k-1 as destination/output or '-'",
        "seed": seeds,
    }
    if p.method == "closure":
        found, pairs, status = expand(chem.react, seeds, arity=2, max_species=p.max_species, ordered=True)
        reactions = [Reaction.of(lhs, rhs) for lhs, rhs in pairs]
        extras["analysis"] = {
            "elastic": dict(chem.stats),
            "states": {s: len(chem.machines[s]) for s in found},
        }
        return Network(species=chem.species(found), reactions=reactions, status=status, extras=extras)
    return _soup(p, chem, drawn, seeds, rng, extras)


def _soup(p, chem, drawn, seeds, rng, extras):
    if p.machines:
        if p.M < len(seeds):
            raise ValueError(f"M={p.M} is smaller than the {len(seeds)} distinct machines given")
        pop = [s for s, n in zip(seeds, apportion(p.M, [1.0] * len(seeds))) for _ in range(n)]
    else:
        pop = list(drawn)
    start = Counter(pop)
    size = len(pop)
    fired: dict[tuple, list] = {}
    diversity = [len(start)]
    done = 0
    while done < p.collisions:
        n = min(size, p.collisions - done)
        chunk, pop = soup(chem.react, pop, n, rng, arity=2, dilution="constant")
        done += n
        for lhs, rhs, count in chunk:
            key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
            fired.setdefault(key, [lhs, rhs, 0])[2] += count
        diversity.append(len(set(pop)))
    reactions = [Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in fired.values()]
    ids = list(dict.fromkeys([*start, *(c for _, rhs, _ in fired.values() for c in rhs)]))
    extras["analysis"] = {
        "collisions_per_sample": size,
        "distinct_species": diversity,
        "elastic": dict(chem.stats),
        "states": {s: len(chem.machines[s]) for s in ids},
    }
    extras["final_state"] = dict(Counter(pop).most_common())
    return Network(species=chem.species(ids), reactions=reactions, status="observed",
                   initial_state={s: float(n) for s, n in start.items()}, outflow=CONSTANT_TOTAL,
                   extras=extras)
