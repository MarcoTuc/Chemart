"""RBN World (book 18.4.2, 10.7.3). Catalog id: rbn-world.

Faulconbridge, Stepney, Miller & Caves's sub-symbolic chemistry. An atom is a
bonding RBN (bRBN, b = 2 bonding sites): a random Boolean network from
chemart.chemistries.rbn with two of its inputs rewired to bonding sites. A
molecule is a tree of bRBNs whose bonds are chosen, kept and broken by an
emergent property of the bRBNs' attractors.

The chemistry has one face: `evolve` runs the well-stirred collision reactor
(chemart.soup.stir), a frame every `elements * copies` collisions, and returns
the reactions observed in the run (status "observed").
"""

from __future__ import annotations

from collections import Counter
from itertools import groupby, permutations, product
from math import factorial, prod
from string import ascii_uppercase

from chemart.chemistries.rbn import bits, random_rbn
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

SITES = 2                     # bonding sites per atom (b = 2 in every RBN World paper)
SETTLE_LIMIT = 1 << 16        # steps allowed to find an attractor before a bRBN counts as undefined
PERMUTATION_LIMIT = 720       # tie-breaking budget of the species canonical form


# ============================================================================
# RBN World: atoms, bRBN dynamics and emergent bonding properties
# ============================================================================
class Element:
    """A bRBN atom: inputs >= 0 are nodes, -1 - s is bonding site s; tables are row bitmasks."""

    def __init__(self, inputs, tables, site_target, initial):
        self.inputs, self.tables, self.site_target, self.initial = inputs, tables, site_target, initial


def random_element(N: int, K: int, bias: float, rng) -> Element:
    inputs, rows = random_rbn(N, [K] * N, bias, rng)
    tables = [sum(b << r for r, b in enumerate(t)) for t in rows]
    site_target = []
    for s, slot in enumerate(rng.choice(N * K, size=SITES, replace=False)):
        node, t = divmod(int(slot), K)
        inputs[node][t] = -1 - s             # redirect one input of a random node to bonding site s
        site_target.append(node)
    initial = sum(int(b) << i for i, b in enumerate(rng.integers(0, 2, size=N)))
    return Element(inputs, tables, site_target, initial)


def cycle_properties(cycle: list[list[int]]) -> dict[str, float]:
    """Bonding properties of an attractor (Faulconbridge et al. 2010, table 1; thesis eqs. 7.1.1-7.1.8).

    cycle[j][i] is the state (0/1) of Boolean node i at step j of the cycle.
    """
    c, n = len(cycle), len(cycle[0])
    ones = sum(map(sum, cycle))
    columns = list(zip(*cycle))
    return {
        "cycle-length": c,
        "flashing": sum(1 for col in columns if 0 < sum(col) < c),
        "flashes": sum(col[j] != col[j - 1] for col in columns for j in range(c)),
        "total": 2 * ones - n * c,
        "magnitude": max(ones, n * c - ones),
        "proportion": ones / (n * c),
    }


def _satisfied(rule: str, a, b) -> bool:
    if a is None or b is None:
        return False
    if rule == "cycle-length-equal":
        return a == b
    if rule == "proportion-sum-one":
        return abs(a + b - 1) <= 0.001
    return abs(a + b) <= 0.001                       # total-sum-zero


class Node:
    """A bRBN in a molecule's structure tree: a leaf atom or a composite of child bRBNs.

    state maps each atom id under the node to that atom's node states as seen at this level.
    """

    __slots__ = ("children", "atoms", "state")

    def __init__(self, children=(), state=None, atom=None):
        self.children = list(children)
        self.atoms = frozenset([atom]) if atom is not None else frozenset().union(*(c.atoms for c in self.children))
        self.state = dict(state or {})


def _parents(roots):
    parent, stack = {}, list(roots)
    while stack:
        node = stack.pop()
        for c in node.children:
            parent[id(c)] = node
            stack.append(c)
    return parent


def _preorder(node):
    yield node
    for c in node.children:
        yield from _preorder(c)


def link_trees(roots: list[Node], x: Node, y: Node) -> list[Node]:
    """Join the trees holding x and y after a bond forms between bRBNs x and y (thesis fig. 6.1.3).

    The parents of x and y are merged level by level ('zipped'); a side that runs
    out of parents carries the merged node up its own tree. Returns the new roots.
    """
    parent = _parents(roots)
    X, Y = parent.get(id(x)), parent.get(id(y))
    if X is None and Y is None:
        merged = Node([x, y], {**x.state, **y.state})
        return [r for r in roots if r is not x and r is not y] + [merged]
    cur_a, cur_b, merged = x, y, None
    while True:
        X, Y = parent.get(id(cur_a)), parent.get(id(cur_b))
        if X is None and Y is None:
            break
        kids = [c for c in (X.children if X else []) if c is not cur_a]
        kids += [c for c in (Y.children if Y else []) if c is not cur_b]
        kids += [cur_a, cur_b] if merged is None else [merged]
        merged = Node(kids, {**(X or cur_a).state, **(Y or cur_b).state})
        cur_a, cur_b = X or cur_a, Y or cur_b
    return [r for r in roots if r is not cur_a and r is not cur_b] + [merged]


def split_tree(node: Node, bonds: dict) -> list[Node]:
    """Rebuild a tree after bonds were removed: composites whose children are no longer
    connected split, and composites of a single bRBN disappear (thesis fig. 6.1.4)."""
    if not node.children:
        return [node]
    pieces = [q for c in node.children for q in split_tree(c, bonds)]
    owner = {a: k for k, q in enumerate(pieces) for a in q.atoms}
    root = list(range(len(pieces)))

    def find(k):
        while root[k] != k:
            root[k] = root[root[k]]
            k = root[k]
        return k

    for (a, _), (b, _) in bonds.items():
        if a in owner and b in owner:
            root[find(owner[a])] = find(owner[b])
    groups: dict[int, list[Node]] = {}
    for k, q in enumerate(pieces):
        groups.setdefault(find(k), []).append(q)
    if len(groups) == 1 and len(pieces) == len(node.children) and all(p is c for p, c in zip(pieces, node.children)):
        return [node]
    return [g[0] if len(g) == 1 else Node(g, {a: node.state[a] for q in g for a in q.atoms})
            for g in groups.values()]


class World:
    """Molecules under reaction: atom elements, bonds (both directions) and structure-tree roots."""

    def __init__(self, elements, N, rule):
        self.elements, self.N, self.rule = elements, N, rule
        self.elem: dict[int, int] = {}
        self.bonds: dict[tuple[int, int], tuple[int, int]] = {}
        self.roots: list[Node] = []

    # -- dynamics -----------------------------------------------------------
    def settle(self, node: Node):
        """Run the node's bRBN from its state to an attractor; keep the first repeated state; return the property."""
        N = self.N
        atoms = sorted(node.atoms)
        pos = {a: k for k, a in enumerate(atoms)}
        M = N * len(atoms)
        zero, one = M, M + 1
        wiring = []
        for a in atoms:
            el = self.elements[self.elem[a]]
            base = pos[a] * N
            for i in range(N):
                ins = []
                for ref in el.inputs[i]:
                    if ref >= 0:
                        ins.append(base + ref)
                        continue
                    partner = self.bonds.get((a, -1 - ref))
                    if partner is None:
                        ins.append(zero)                         # empty site reads 0
                    elif partner[0] in pos:                      # both ends inside: reciprocal input
                        b, t = partner
                        ins.append(pos[b] * N + self.elements[self.elem[b]].site_target[t])
                    else:
                        ins.append(one)                          # filled site reads 1
                wiring.append((ins, el.tables[i]))
        s = sum(node.state[a] << (pos[a] * N) for a in atoms)
        seen, seq, high = {}, [], 1 << one
        while s not in seen:
            if len(seq) >= SETTLE_LIMIT:
                return None
            seen[s] = len(seq)
            seq.append(s)
            x, s = s | high, 0
            for i, (ins, table) in enumerate(wiring):
                r = 0
                for g in ins:
                    r = (r << 1) | ((x >> g) & 1)
                if (table >> r) & 1:
                    s |= 1 << i
        mask = (1 << N) - 1
        for a in atoms:
            node.state[a] = (s >> (pos[a] * N)) & mask
        cycle = seq[seen[s]:]
        if self.rule == "cycle-length-equal":
            return len(cycle)
        ones = sum(v.bit_count() for v in cycle)
        return ones / (M * len(cycle)) if self.rule == "proportion-sum-one" else 2 * ones - M * len(cycle)

    # -- reaction (thesis fig. 6.1.1) ---------------------------------------
    def react(self, site_a, site_b):
        rule = self.rule
        parent = _parents(self.roots)
        leaf = {next(iter(n.atoms)): n for r in self.roots for n in _preorder(r) if not n.children}
        x, y = leaf[site_a[0]], leaf[site_b[0]]
        pa, pb = self.settle(x), self.settle(y)
        while not _satisfied(rule, pa, pb):                   # TestForInteraction: climb to larger bRBNs
            X, Y = parent.get(id(x)), parent.get(id(y))
            if X is None and Y is None:
                x = None
                break
            if X is not None and Y is not None:
                if len(X.atoms) >= len(Y.atoms):
                    x, pa = X, self.settle(X)
                if len(Y.atoms) >= len(X.atoms):
                    y, pb = Y, self.settle(Y)
            elif X is not None:
                x, pa = X, self.settle(X)
            else:
                y, pb = Y, self.settle(Y)
        if x is not None:
            self.bonds[site_a], self.bonds[site_b] = site_b, site_a      # fill both sites
            if _satisfied(rule, self.settle(x), self.settle(y)):
                self.roots = link_trees(self.roots, x, y)                 # MakeStableBond
            else:
                del self.bonds[site_a], self.bonds[site_b]                # empty them again
        self.decompose()

    def decompose(self):
        """Break every bond whose two bRBNs no longer satisfy the rule, until none breaks."""
        while True:
            broken = []
            for root in self.roots:
                for T in _preorder(root):
                    if not T.children:
                        continue
                    owner = {a: c for c in T.children for a in c.atoms}
                    for a in sorted(T.atoms):
                        for s in range(SITES):
                            partner = self.bonds.get((a, s))
                            if partner is None or (a, s) > partner or partner[0] not in owner:
                                continue
                            cx, cy = owner[a], owner[partner[0]]
                            if cx is not cy and not _satisfied(self.rule, self.settle(cx), self.settle(cy)):
                                broken.append((a, s))
            if not broken:
                return
            for site in broken:
                del self.bonds[self.bonds.pop(site)]
            self.roots = [q for r in self.roots for q in split_tree(r, self.bonds)]


# ============================================================================
# RBN World: species identity
# ============================================================================
class Registry:
    """Canonical molecules. A species is a structure tree with its bonds and the states at every level."""

    def __init__(self, elements, N, rule):
        self.elements, self.N, self.rule = elements, N, rule
        self.molecules: dict[str, tuple] = {}     # code -> (elem, bonds, root) with atom ids 0..m-1
        self.ids: dict[str, str] = {}
        self.names: Counter = Counter()

    def _canon(self, node, elem, bonds):
        N = self.N
        if not node.children:
            (a,) = node.atoms
            letter = ascii_uppercase[elem[a]]
            return f"{letter}:{bits(node.state[a], N)}", letter, [a]
        subs = sorted((self._canon(c, elem, bonds) for c in node.children), key=lambda s: s[0])
        groups = [list(g) for _, g in groupby(range(len(subs)), key=lambda k: subs[k][0])]
        if prod(factorial(len(g)) for g in groups) > PERMUTATION_LIMIT:
            orderings = [list(range(len(subs)))]
        else:
            orderings = [[k for part in combo for k in part] for combo in product(*(permutations(g) for g in groups))]
        best = None
        for order in orderings:
            rank = {a: (r, i) for r, k in enumerate(order) for i, a in enumerate(subs[k][2])}
            links = sorted(
                min((rank[a] + (s,), rank[b] + (t,)), (rank[b] + (t,), rank[a] + (s,)))
                for (a, s), (b, t) in bonds.items()
                if a in rank and b in rank and rank[a][0] != rank[b][0] and (a, s) < (b, t)
            )
            atoms = [a for k in order for a in subs[k][2]]
            key = (";".join(f"{u[0]}.{u[1]}.{u[2]}={v[0]}.{v[1]}.{v[2]}" for u, v in links),
                   "".join(bits(node.state[a], N) for a in atoms), atoms)
            best = key if best is None or key[:2] < best[:2] else best
        links, state, atoms = best
        code = f"({','.join(s[0] for s in subs)}|{links}|{state})"
        return code, f"({'-'.join(s[1] for s in subs)})", atoms

    def add(self, root, elem, bonds) -> str:
        code, name, atoms = self._canon(root, elem, bonds)
        if code not in self.molecules:
            new = {a: k for k, a in enumerate(atoms)}

            def copy(n):
                if not n.children:
                    (a,) = n.atoms
                    return Node(state={new[a]: n.state[a]}, atom=new[a])
                return Node([copy(c) for c in n.children], {new[a]: v for a, v in n.state.items()})

            self.molecules[code] = ({new[a]: elem[a] for a in atoms},
                                    {(new[a], s): (new[b], t) for (a, s), (b, t) in bonds.items() if a in new},
                                    copy(root))
            self.names[name] += 1
            self.ids[code] = f"{name}.{self.names[name]}"
        return code

    def load(self, world: World, code: str, offset: int):
        elem, bonds, root = self.molecules[code]

        def copy(n):
            if not n.children:
                (a,) = n.atoms
                return Node(state={a + offset: n.state[a]}, atom=a + offset)
            return Node([copy(c) for c in n.children], {a + offset: v for a, v in n.state.items()})

        world.elem.update({a + offset: e for a, e in elem.items()})
        world.bonds.update({(a + offset, s): (b + offset, t) for (a, s), (b, t) in bonds.items()})
        world.roots.append(copy(root))
        return sorted((a + offset, s) for a in elem for s in range(SITES) if (a, s) not in bonds)

    def products(self, world: World) -> tuple[str, ...]:
        out = []
        for root in world.roots:
            elem = {a: world.elem[a] for a in root.atoms}
            out.append(self.add(root, elem, {k: v for k, v in world.bonds.items() if k[0] in root.atoms}))
        return tuple(out)


# ============================================================================
# RBN World: the collision run
# ============================================================================
SPECIES_ENCODING = (
    "id = readable structure (children joined by '-', one letter per atom) + '.' + index of first "
    "observation; structure = (child codes|bonds as rank.atom.site=rank.atom.site|node states at "
    "this level), atoms as LETTER:states, node 0 first"
)


def evolve(p, rng):
    """Iterative mixing of random pairs (thesis fig. 8.2.1): a frame per elements * copies collisions."""
    if p.K > p.N:
        raise ValueError(f"K must satisfy 1 <= K <= N, got K = {p.K} with N = {p.N}")
    if p.N * p.K < SITES:
        raise ValueError(f"an rbn-world atom needs N * K >= {SITES} input slots for its bonding sites")
    N, rule = p.N, p.bonding_rule
    elements = [random_element(N, p.K, p.function_bias, rng) for _ in range(p.elements)]
    registry = Registry(elements, N, rule)

    atoms = []
    for e in range(p.elements):
        world = World(elements, N, rule)
        world.elem[0] = e
        leaf = Node(state={0: elements[e].initial}, atom=0)
        world.roots = [leaf]
        world.settle(leaf)                       # atoms enter the reactor on their attractor
        atoms.append(registry.add(leaf, {0: e}, {}))

    cache: dict[tuple, tuple] = {}

    def react(a, b):
        world = World(elements, N, rule)
        sites_a = registry.load(world, a, 0)
        sites_b = registry.load(world, b, len(registry.molecules[a][0]))
        i, j = int(rng.integers(len(sites_a))), int(rng.integers(len(sites_b)))
        key = (a, b, i, j)
        if key not in cache:
            world.react(sites_a[i], sites_b[j])
            cache[key] = registry.products(world)
        return cache[key]

    ids = registry.ids
    size = {}                                    # code -> atoms in the molecule

    def atoms_in(code):
        if code not in size:
            size[code] = len(registry.molecules[code][0])
        return size[code]

    population = [code for code in atoms for _ in range(p.copies)]
    tally = Tally()
    final = population
    for step, final, tally in stir(react, population, p.collisions, rng, arity=2, tally=tally):
        present = Counter(final)
        yield Frame(
            t=float(step),
            state={ids[c]: float(n) for c, n in present.items()},
            fired=[[[ids[c] for c in lhs], [ids[c] for c in rhs], n] for lhs, rhs, n in tally.flush()],
            observables={"largest_molecule_atoms": max(map(atoms_in, present), default=0)},
        )
    fired = tally.reactions()

    used = dict.fromkeys(atoms)
    for lhs, rhs, _ in fired:
        used.update(dict.fromkeys(lhs + rhs))
    species = [Species(ids[c], structure=c) for c in used]
    reactions = [Reaction.of([ids[c] for c in lhs], [ids[c] for c in rhs], count=int(n)) for lhs, rhs, n in fired]
    letters = ascii_uppercase[: p.elements]
    conservation = [{"name": f"atoms of element {letters[e]}",
                     "vector": {ids[c]: sum(1 for x in registry.molecules[c][0].values() if x == e) for c in used}}
                    for e in range(p.elements)]
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={ids[c]: float(p.copies) for c in atoms},
        extras={
            "conservation": conservation,
            "elements": [{
                "name": letters[e],
                "inputs": [[int(r) if r >= 0 else f"site{-1 - r}" for r in ins] for ins in el.inputs],
                "functions": ["".join(str((t >> r) & 1) for r in range(2 ** p.K)) for t in el.tables],
                "bonding_site_targets": el.site_target,
                "initial_state": bits(el.initial, N),
            } for e, el in enumerate(elements)],
            "species_encoding": SPECIES_ENCODING,
            "analysis": {
                "final_population": {ids[c]: n for c, n in Counter(final).items()},
                "largest_molecule_atoms": max(atoms_in(c) for c in used),
                "distinct_collision_outcomes": len(cache),
            },
        },
    )
