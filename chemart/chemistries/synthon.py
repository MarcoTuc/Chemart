"""Synthon: molecules as graphs of atoms and electrons, reactions as isomerisations.

Catalog id: synthon. Book 18.3.2 (one paragraph); Lenaerts & Bersini, "A synthon
approach to artificial chemistry", Artificial Life 15(1):89-103 (2009) (book
[505]); reconstructed from the authors' own description of the same framework in
Lenaerts & Bersini, "On the generation and analysis of complex reaction networks
in interstellar chemistry" (2005), and from the Synthon model of Koca (1988) and
the Dugundji-Ugi model it extends.

A *synthon* S(A) = <W, A, R, E> is a set of atoms A (plus virtual atoms W, not
used here), explicitly modelled electrons R and edges E: covalent bonds between
atoms, links between an atom and a free (unpaired) electron, and loops standing
for lone pairs. A reaction is an *isomerisation* S(A) => S'(A) of one ensemble of
molecules into another: both sides carry the same atoms and the same electrons,
only the edges change. A reaction class is a reaction graph G_T = <W, A, R, E,
psi> with psi(e) = +1 for a formed edge and -1 for a deleted one; a *reaction
object* is one instance of a class on concrete molecules.

Chemart's molecules follow that representation:

- an atom carries its element, its lone pairs and its unpaired electrons
  (radicals), so radicals and charges are explicit;
- bonds are covalent (an integer order, two electrons per order) or *ionic* -
  the extension the authors added for H2+ and H3+, which they represent as two
  disjoint parts of one molecule (H+ . H and H+ . H2);
- the free electron e- is a species of its own.

The formal charge of an atom follows the Dugundji-Ugi bookkeeping,
q = v - 2*lone_pairs - radicals - sum of covalent bond orders, with v the
valence electrons of the free atom; a molecule's charge is the sum over its
atoms. Every reaction conserves the atoms of each element, the electrons and
therefore the total charge (extras["conservation"]).

Species identity is the canonical code of the graph (the role of the CANGEN
algorithm in the papers): atoms are ordered canonically and written as
element + ":"*lone_pairs + "*"*radicals, covalent bonds by juxtaposition
("=" double, "#" triple), an ionic link as ".", branches in parentheses and ring
closures as "%n" - so H2 is "HH", the H atom "H*", the proton "H", H2+ "H.H*",
H3+ "H.HH", water "HO::H" and the free electron "e-". extras["formulas"] gives
the chemical formula of each species.

Reaction classes (`templates`):

- "interstellar": the eleven classes G1-G11 of table 1 of the 2005 paper, with
  the rate constants it takes from Duley & Williams, "Interstellar chemistry"
  (1984) p. 143;
- "fig3": the three classes drawn in its figure 3 (radiative association, charge
  transfer, dissociative recombination), for which no rates are published.

Methods (`method`):

- "closure" is the deterministic network generator (DNG): every reaction object
  reachable from the initial species under the size constraints N (atoms per
  species), es (lone pairs or radicals per species) and ep (per atom);
- "kinetic" is the MC-sampling generator (MCNG): Gillespie's SSA over the
  molecules present, with the network grown from the species that actually have
  a non-zero amount, returning the reactions that fired with their counts.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from itertools import permutations

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

METHODS = ("closure", "kinetic")
TEMPLATE_SETS = ("interstellar", "fig3")

#: The molecules of figure 4 of the 2005 paper, in the notation of its caption
#: (H=H+, *H=H, HH=H2, H.*H=H2+, H.HH=H3+, *O*=O, HO=OH+, *OH=OH, HOH=H2O,
#: *OO*=O2, *O(H)H=OH2+, HO(H)H=OH3+, *=e-), as Chemart formulas. This is the
#: "observational constraint" the paper used to keep that network drawable.
FIGURE_4 = ("H^+", "H", "H2", "H2^+", "H3^+", "O", "HO^+", "HO", "H2O", "O2",
            "H2O^+", "H3O^+", "e-")


@dataclass(frozen=True)
class Limits:
    """The growth constraints of the papers: N, es, ep and the allowed molecules."""

    max_atoms: int
    max_pairs_species: int
    max_pairs_atom: int
    max_charge: int = 1
    allowed: frozenset[str] | None = None

#: element -> (valence electrons, electrons that fill the shell)
ELEMENTS: dict[str, tuple[int, int]] = {
    "H": (1, 2), "He": (2, 2), "C": (4, 8), "N": (5, 8), "O": (6, 8),
    "F": (7, 8), "Ne": (8, 8), "S": (6, 8), "Cl": (7, 8), "Ar": (8, 8),
}

#: 2005 paper, table 1: class, chemical class name, reaction template, k, units.
#: (c.r. = cosmic ray, hv = photon; both carry no mass and are not species.)
INTERSTELLAR = (
    ("G1", "cosmic ray ionization", "A + c.r. -> A+ + e-", 1e-17, "s^-1"),
    ("G2", "cosmic ray ionization", "AB + c.r. -> AB+ + e-", 1e-17, "s^-1"),
    ("G3", "cosmic ray ionization", "AB + c.r. -> A + B+ + e-", 1e-19, "s^-1"),
    ("G4", "ion-molecule exchange", "A+ + BC -> AB+ + C", 1e-9, "cm^3 s^-1"),
    ("G5", "charge rearrangement", "A+ + B -> A + B+", 1e-9, "cm^3 s^-1"),
    ("G6", "dissociative recombination", "AB+ + e- -> A + B", 1e-6, "cm^3 s^-1"),
    ("G7", "dissociative recombination", "A+ + e- -> A", 1e-11, "cm^3 s^-1"),
    ("G8", "neutral reactions", "A + BC -> AB + C", 1e-11, "cm^3 s^-1"),
    ("G9", "photodissociation", "AB + hv -> A + B", 1e-11, "s^-1"),
    ("G10", "photodissociation", "H2 + hv -> 2 H", 1e-14, "s^-1"),
    ("G11", "grain surface reaction", "H + H :g -> H2 + g", 1e-17, "cm^3 s^-1"),
)

#: 2005 paper, figure 3: the three reaction graphs drawn there. No rates given.
FIG3 = (
    ("ra", "radiative association", "A + B -> AB + hv", None, None),
    ("ct", "charge transfer", "A+ + B -> A + B+", None, None),
    ("dr", "dissociative recombination", "AB+ + e- -> A + B", None, None),
)

BOND_TEXT = {1: "", 2: "=", 3: "#"}


# ============================================================================
# Molecules
# ============================================================================
@dataclass(frozen=True)
class Mol:
    """A species: atoms (element, lone pairs, radicals) joined by bonds.

    `bonds` are covalent (i, j, order) with i < j and order >= 1; `ionic` are
    (i, j) links that hold charged parts together without sharing electrons.
    `kind` is "molecule" or "electron" (the free electron e-, no atoms).
    """

    atoms: tuple[tuple[str, int, int], ...] = ()
    bonds: tuple[tuple[int, int, int], ...] = ()
    ionic: tuple[tuple[int, int], ...] = ()
    kind: str = "molecule"

    # -- bookkeeping ---------------------------------------------------------
    @property
    def charge(self) -> int:
        if self.kind == "electron":
            return -1
        orders = Counter()
        for i, j, o in self.bonds:
            orders[i] += o
            orders[j] += o
        return sum(ELEMENTS[el][0] - 2 * lp - rad - orders[i]
                   for i, (el, lp, rad) in enumerate(self.atoms))

    @property
    def electrons(self) -> int:
        if self.kind == "electron":
            return 1
        return (sum(2 * lp + rad for _, lp, rad in self.atoms)
                + 2 * sum(o for _, _, o in self.bonds))

    @property
    def atom_counts(self) -> Counter:
        return Counter(el for el, _, _ in self.atoms)

    @property
    def size(self) -> int:
        return len(self.atoms)

    @property
    def id(self) -> str:
        return "e-" if self.kind == "electron" else _text(self.atoms, self.bonds, self.ionic)

    @property
    def formula(self) -> str:
        if self.kind == "electron":
            return "e-"
        counts = self.atom_counts
        order = ([e for e in ("C", "H") if e in counts]
                 + sorted(e for e in counts if e not in ("C", "H")))
        text = "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)
        q = self.charge
        if q:
            text += "^" + (str(abs(q)) if abs(q) > 1 else "") + ("+" if q > 0 else "-")
        return text


ELECTRON = Mol(kind="electron")


def ground_radicals(element: str) -> int:
    """Unpaired electrons of the neutral free atom (H 1, O 2, N 3, C 4)."""
    v, shell = ELEMENTS[element]
    return min(v, shell - v)


def _part_charges(atoms, bonds) -> list[int]:
    """Charge of each covalently connected part of one molecule."""
    n = len(atoms)
    root = list(range(n))

    def find(i):
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    for i, j in bonds:
        root[find(i)] = find(j)
    orders: Counter = Counter()
    for (i, j), o in bonds.items():
        orders[i] += o
        orders[j] += o
    parts: dict[int, list[int]] = {}
    for i in range(n):
        parts.setdefault(find(i), []).append(i)
    return [sum(ELEMENTS[atoms[k][0]][0] - 2 * atoms[k][1] - atoms[k][2] - orders[k] for k in group)
            for group in parts.values()]


def _ionic_ok(atoms, bonds, ionic) -> bool:
    """An ionic link holds one cation and one neutral part together (H+ . H2, H+ . H)."""
    if len(ionic) != 1:
        return False
    charges = _part_charges(atoms, bonds)
    return len(charges) == 2 and max(charges) > 0 and min(charges) == 0


def atom(element: str) -> Mol:
    """The neutral free atom of an element: unpaired electrons up to the octet."""
    if element not in ELEMENTS:
        raise ValueError(f"unknown element {element!r}; known: {', '.join(ELEMENTS)}")
    v, shell = ELEMENTS[element]
    radicals = min(v, shell - v)
    return Mol(((element, (v - radicals) // 2, radicals),))


def molecule(atoms, bonds, ionic) -> Mol:
    """Canonical Mol from raw atoms/bonds/ionic links (any vertex order)."""
    order = _canonical_order(atoms, bonds, ionic)
    index = {v: i for i, v in enumerate(order)}
    return Mol(
        tuple(atoms[v] for v in order),
        tuple(sorted((min(index[i], index[j]), max(index[i], index[j]), o)
                     for (i, j), o in bonds.items())),
        tuple(sorted((min(index[i], index[j]), max(index[i], index[j]))
                     for i, j in ionic)),
    )


# -- canonical form (the role of CANGEN in the papers) ------------------------
def _neighbours(n: int, bonds: dict, ionic) -> list[list[tuple[int, int]]]:
    """adj[i] = [(neighbour, label)], label 0 for an ionic link, else bond order."""
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    for (i, j), o in bonds.items():
        adj[i].append((j, o))
        adj[j].append((i, o))
    for i, j in ionic:
        adj[i].append((j, 0))
        adj[j].append((i, 0))
    return adj


def _refine(atoms, adj) -> list[int]:
    """Morgan-style refinement: a rank per atom, equal ranks = same environment."""
    n = len(atoms)
    keys = [(atoms[i], tuple(sorted(lab for _, lab in adj[i]))) for i in range(n)]
    rank = _rank(keys)
    for _ in range(n):
        keys = [(rank[i], tuple(sorted((lab, rank[j]) for j, lab in adj[i]))) for i in range(n)]
        new = _rank(keys)
        if new == rank:
            break
        rank = new
    return rank


def _rank(keys) -> list[int]:
    order = {k: i for i, k in enumerate(sorted(set(keys)))}
    return [order[k] for k in keys]


def _score(atoms, bonds, ionic, order) -> tuple:
    index = {v: i for i, v in enumerate(order)}
    return (
        tuple(atoms[v] for v in order),
        tuple(sorted((min(index[i], index[j]), max(index[i], index[j]), o)
                     for (i, j), o in bonds.items())),
        tuple(sorted((min(index[i], index[j]), max(index[i], index[j]))
                     for i, j in ionic)),
    )


def _canonical_order(atoms, bonds, ionic) -> list[int]:
    """Atom order that minimises the graph's code; ties broken by search."""
    n = len(atoms)
    if n <= 1:
        return list(range(n))
    rank = _refine(atoms, _neighbours(n, bonds, ionic))
    cells: dict[int, list[int]] = {}
    for i, r in enumerate(rank):
        cells.setdefault(r, []).append(i)
    sizes = [len(c) for c in cells.values()]
    total = math.prod(math.factorial(s) for s in sizes)
    base = [cells[r] for r in sorted(cells)]
    if total > 5040:                      # symmetric and large: keep the refinement order
        return [v for cell in base for v in cell]
    best = None
    for combo in _orderings(base):
        score = _score(atoms, bonds, ionic, combo)
        if best is None or score < best[0]:
            best = (score, combo)
    return best[1]


def _orderings(cells):
    if not cells:
        yield []
        return
    head, rest = cells[0], cells[1:]
    for first in permutations(head):
        for tail in _orderings(rest):
            yield list(first) + tail


def _text(atoms, bonds, ionic) -> str:
    """Linear code of a canonically ordered molecule (see the module docstring)."""
    n = len(atoms)
    bond_map = {}
    if isinstance(bonds, dict):
        bonds = tuple((i, j, o) for (i, j), o in bonds.items())
    for i, j, o in bonds:
        bond_map[(min(i, j), max(i, j))] = o
    for i, j in ionic:
        bond_map[(min(i, j), max(i, j))] = 0
    adj: list[list[int]] = [[] for _ in range(n)]
    for (i, j) in bond_map:
        adj[i].append(j)
        adj[j].append(i)
    for a in adj:
        a.sort()

    def label(o):
        return "." if o == 0 else BOND_TEXT.get(o, f"~{o}~")

    def token(i):
        return f"{atoms[i][0]}{':' * atoms[i][1]}{'*' * atoms[i][2]}"

    # First pass: the same traversal, only to number the ring-closure edges.
    used: set[tuple[int, int]] = set()
    visited: set[int] = set()
    marks: dict[int, list[str]] = {}
    counter = 0

    def scan(i: int) -> None:
        nonlocal counter
        visited.add(i)
        for j in adj[i]:
            edge = (min(i, j), max(i, j))
            if edge in used:
                continue
            used.add(edge)
            if j in visited:
                counter += 1
                marks.setdefault(i, []).append(f"%{counter}")
                marks.setdefault(j, []).append(f"%{counter}")
            else:
                scan(j)

    scan(0)

    # Second pass: build the code, ring marks now known for every atom.
    used.clear()
    visited.clear()

    def walk(i: int) -> str:
        visited.add(i)
        parts = []
        for j in adj[i]:
            edge = (min(i, j), max(i, j))
            if edge in used:
                continue
            used.add(edge)
            if j in visited:
                continue                                    # ring closure, marked above
            parts.append(label(bond_map[edge]) + walk(j))
        body = "".join(f"({p})" for p in parts[:-1]) + (parts[-1] if parts else "")
        return token(i) + "".join(marks.get(i, ())) + body

    return walk(0)


# ============================================================================
# Ensembles and the reaction graphs G_T
# ============================================================================
class Ensemble:
    """A mutable synthon: the merged reactants, edited by a reaction graph."""

    def __init__(self, mols):
        self.atoms: list[tuple[str, int, int]] = []
        self.bonds: dict[tuple[int, int], int] = {}
        self.ionic: set[tuple[int, int]] = set()
        self.loose: list[Mol] = []
        self.parts: list[list[int]] = []
        for m in mols:
            if m.kind != "molecule":
                self.loose.append(m)
                self.parts.append([])
                continue
            off = len(self.atoms)
            self.atoms.extend(m.atoms)
            for i, j, o in m.bonds:
                self.bonds[(i + off, j + off)] = o
            for i, j in m.ionic:
                self.ionic.add((i + off, j + off))
            self.parts.append(list(range(off, off + len(m.atoms))))

    def copy(self) -> "Ensemble":
        clone = Ensemble(())
        clone.atoms = list(self.atoms)
        clone.bonds = dict(self.bonds)
        clone.ionic = set(self.ionic)
        clone.loose = list(self.loose)
        clone.parts = [list(p) for p in self.parts]
        return clone

    # -- edits -------------------------------------------------------------
    def set_atom(self, i: int, lp: int, rad: int) -> None:
        el, _, _ = self.atoms[i]
        self.atoms[i] = (el, lp, rad)

    def lone_pairs(self, i: int) -> int:
        return self.atoms[i][1]

    def radicals(self, i: int) -> int:
        return self.atoms[i][2]

    def break_edge(self, i: int, j: int, homolytic: bool = True) -> None:
        key = (min(i, j), max(i, j))
        if key in self.ionic:
            self.ionic.discard(key)
            return
        order = self.bonds[key] - 1
        if order:
            self.bonds[key] = order
        else:
            del self.bonds[key]
        if homolytic:
            for k in (i, j):
                el, lp, rad = self.atoms[k]
                self.atoms[k] = (el, lp, rad + 1)
        else:                                   # heterolytic: i keeps the pair
            el, lp, rad = self.atoms[i]
            self.atoms[i] = (el, lp + 1, rad)

    def bind(self, i: int, j: int) -> bool:
        """Covalent bond if both atoms carry a radical, else an ionic link."""
        key = (min(i, j), max(i, j))
        if self.radicals(i) and self.radicals(j) and i != j:
            for k in (i, j):
                el, lp, rad = self.atoms[k]
                self.atoms[k] = (el, lp, rad - 1)
            self.bonds[key] = self.bonds.get(key, 0) + 1
            return True
        if key in self.bonds or key in self.ionic:
            return False
        self.ionic.add(key)
        return True

    # -- products ------------------------------------------------------------
    def species(self, limits) -> tuple[Mol, ...] | None:
        """Split into molecules; None if a size constraint is violated."""
        n = len(self.atoms)
        root = list(range(n))

        def find(i):
            while root[i] != i:
                root[i] = root[root[i]]
                i = root[i]
            return i

        for i, j in list(self.bonds) + list(self.ionic):
            root[find(i)] = find(j)
        groups: dict[int, list[int]] = {}
        for i in range(n):
            groups.setdefault(find(i), []).append(i)

        out = list(self.loose)
        for comp in groups.values():
            if len(comp) > limits.max_atoms:
                return None
            index = {v: k for k, v in enumerate(comp)}
            atoms = [self.atoms[v] for v in comp]
            if sum(lp + rad for _, lp, rad in atoms) > limits.max_pairs_species:
                return None
            if any(lp + rad > limits.max_pairs_atom for _, lp, rad in atoms):
                return None
            bonds = {(index[i], index[j]): o for (i, j), o in self.bonds.items() if i in index}
            ionic = {(index[i], index[j]) for i, j in self.ionic if i in index}
            orders: Counter = Counter()
            for (a, b), o in bonds.items():
                orders[a] += o
                orders[b] += o
            for k, (el, lp, rad) in enumerate(atoms):
                if rad > ground_radicals(el):        # no more unpaired electrons than the free atom
                    return None
                if 2 * lp + rad + 2 * orders[k] > ELEMENTS[el][1]:      # duet / octet
                    return None
            if ionic and not _ionic_ok(atoms, bonds, ionic):
                return None
            mol = molecule(atoms, bonds, ionic)
            if abs(mol.charge) > limits.max_charge:
                return None
            if limits.allowed and mol.formula not in limits.allowed:
                return None
            out.append(mol)
        return tuple(out)


def _ionisations(ens: Ensemble):
    """X -> X+ + e-: one electron leaves an atom (from a radical or a lone pair)."""
    for i in range(len(ens.atoms)):
        if ens.radicals(i):
            out = ens.copy()
            out.set_atom(i, out.lone_pairs(i), out.radicals(i) - 1)
            out.loose.append(ELECTRON)
            yield out
        if ens.lone_pairs(i):
            out = ens.copy()
            out.set_atom(i, out.lone_pairs(i) - 1, out.radicals(i) + 1)
            out.loose.append(ELECTRON)
            yield out


def _cleavages(ens: Ensemble):
    """Every way to break one edge: a covalent bond homolytically, or an ionic link.

    The eleven published classes only ever split a molecule into neutral radical
    fragments (AB + hv -> A + B) or undo an ionic link, so heterolytic cleavage,
    which would manufacture ion pairs out of neutral bonds, is not offered.
    """
    for (i, j), _ in list(ens.bonds.items()):
        out = ens.copy()
        out.break_edge(i, j, True)
        yield out
    for i, j in list(ens.ionic):
        out = ens.copy()
        out.break_edge(i, j)
        yield out


def _electron_capture(ens: Ensemble, atoms):
    """An atom takes up a free electron: it pairs with a radical or becomes one."""
    for i in atoms:
        out = ens.copy()
        if out.radicals(i):
            out.set_atom(i, out.lone_pairs(i) + 1, out.radicals(i) - 1)
        else:
            out.set_atom(i, out.lone_pairs(i), out.radicals(i) + 1)
        yield out


def _electron_moves(ens: Ensemble, donor, acceptor):
    """One electron moves from an atom of the donor to an atom of the acceptor.

    Only an unpaired electron moves: breaking a lone pair to hand over one of its
    electrons is an ionisation, which is what classes G1-G3 are for.
    """
    for i in donor:
        givers = []
        if ens.radicals(i):
            givers.append((ens.lone_pairs(i), ens.radicals(i) - 1))
        for lp, rad in givers:
            for j in acceptor:
                out = ens.copy()
                out.set_atom(i, lp, rad)
                if out.radicals(j):
                    out.set_atom(j, out.lone_pairs(j) + 1, out.radicals(j) - 1)
                else:
                    out.set_atom(j, out.lone_pairs(j), out.radicals(j) + 1)
                yield out


def _transfers(ens: Ensemble, donor, acceptor):
    """A fragment leaves the donor molecule and binds to an atom of the acceptor."""
    for (i, j), _ in list(ens.bonds.items()):
        if i not in donor or j not in donor:
            continue
        for moving, staying in ((i, j), (j, i)):
            broken = ens.copy()
            broken.break_edge(moving, staying, True)
            for a in acceptor:
                out = broken.copy()
                if out.bind(a, moving):
                    yield out


# ============================================================================
# Reaction classes
# ============================================================================
def _outcomes(mols: tuple[Mol, ...], classes: set[str], limits):
    """Yield (products, class) for every reaction object of the given classes."""
    kinds = [m.kind for m in mols]
    charges = [m.charge for m in mols]
    ens = Ensemble(mols)
    parts = ens.parts

    def emit(edited: Ensemble, cls: str, want: int | None = None):
        products = edited.species(limits)
        if products is None:
            return None
        if want is not None and len(products) != want:
            return None
        return (products, cls)

    if len(mols) == 1 and kinds[0] == "molecule":
        mol = mols[0]
        # G1/G2: ionisation of an atom (G1) or a molecule (G2)
        cls = "G1" if mol.size == 1 else "G2"
        if cls in classes:
            for edited in _ionisations(ens):
                got = emit(edited, cls, want=2)
                if got:
                    yield got
        # G3: ionisation with dissociation
        if "G3" in classes:
            for ionised in _ionisations(ens):
                for edited in _cleavages(ionised):
                    got = emit(edited, "G3", want=3)
                    if got:
                        yield got
        # G9/G10: photodissociation (G10 is H2 specifically)
        cls = "G10" if mol.formula == "H2" else "G9"
        if cls in classes:
            for edited in _cleavages(ens):
                got = emit(edited, cls, want=2)
                if got:
                    yield got
        return

    if len(mols) != 2:
        return

    for x, y in ((0, 1), (1, 0)):
        if kinds[x] == "electron" and kinds[y] == "molecule" and charges[y] > 0:
            # G7: A+ + e- -> A;  G6/dr: AB+ + e- -> A + B. The electron is consumed.
            base = ens.copy()
            base.loose = [m for m in base.loose if m.kind != "electron"]
            if "G7" in classes:
                for edited in _electron_capture(base, parts[y]):
                    got = emit(edited, "G7", want=1)
                    if got:
                        yield got
            cls = "G6" if "G6" in classes else ("dr" if "dr" in classes else None)
            if cls:
                for captured in _electron_capture(base, parts[y]):
                    for edited in _cleavages(captured):
                        got = emit(edited, cls, want=2)
                        if got:
                            yield got
        if kinds[x] != "molecule" or kinds[y] != "molecule":
            continue
        # G5 / ct: charge rearrangement, an electron moves to the cation
        cls = "G5" if "G5" in classes else ("ct" if "ct" in classes else None)
        if cls and charges[x] > 0:
            for edited in _electron_moves(ens, parts[y], parts[x]):
                got = emit(edited, cls, want=2)
                if got:
                    yield got
        # G4 (ion + molecule) and G8 (two neutrals): a fragment changes molecule
        cls = None
        if "G4" in classes and charges[x] > 0 and charges[y] == 0:
            cls = "G4"
        elif "G8" in classes and charges[x] == 0 and charges[y] == 0:
            cls = "G8"
        if cls:
            for edited in _transfers(ens, parts[y], parts[x]):
                got = emit(edited, cls, want=2)
                if got:
                    yield got
        # ra: radiative association, two radicals bind
        if "ra" in classes:
            for i in parts[x]:
                for j in parts[y]:
                    if ens.radicals(i) and ens.radicals(j):
                        edited = ens.copy()
                        edited.bind(i, j)
                        got = emit(edited, "ra", want=1)
                        if got:
                            yield got
    # G11: two hydrogen atoms make H2 on a grain (the grain is implicit, see notes)
    if "G11" in classes and all(m.kind == "molecule" and m.formula == "H" for m in mols):
        edited = ens.copy()
        if edited.bind(parts[0][0], parts[1][0]):
            got = emit(edited, "G11", want=1)
            if got:
                yield got


# ============================================================================
# Generator
# ============================================================================
def templates(name: str) -> tuple[tuple, ...]:
    return INTERSTELLAR if name == "interstellar" else FIG3


def _seeds(p) -> list[Mol]:
    if not isinstance(p.initial, list) or not p.initial or not all(
            isinstance(e, str) for e in p.initial):
        raise ValueError(f"initial must be a non-empty list of element symbols, got {p.initial!r}")
    seeds = [atom(e) for e in p.initial]
    if len(set(m.id for m in seeds)) != len(seeds):
        raise ValueError(f"initial lists an element twice: {p.initial!r}")
    return seeds


def _allowed(p) -> frozenset[str] | None:
    """The observational constraint: the formulas a molecule may have, or None."""
    if not p.allowed:
        return None
    if not isinstance(p.allowed, list) or not all(isinstance(f, str) and f for f in p.allowed):
        raise ValueError(
            f"allowed must be a list of chemical formulas like {list(FIGURE_4[:3])}, "
            f"got {p.allowed!r}; an empty list removes the constraint")
    return frozenset(p.allowed)


def _amounts(p, seeds) -> list[float]:
    values = p.densities if p.densities else [1.0] * len(seeds)
    if (not isinstance(values, list) or len(values) != len(seeds)
            or not all(isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0
                       for v in values)):
        raise ValueError(
            f"densities must be a list of {len(seeds)} positive numbers, one per initial "
            f"species, got {p.densities!r}")
    return [float(v) for v in values]


def _reaction_table(p):
    table = {name: (cls, text, k, units) for name, cls, text, k, units in templates(p.templates)}
    return table, set(table)


def generate(p, rng):
    if p.method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {p.method!r}")
    if p.templates not in TEMPLATE_SETS:
        raise ValueError(f"templates must be one of {TEMPLATE_SETS}, got {p.templates!r}")
    seeds = _seeds(p)
    amounts = _amounts(p, seeds)
    table, classes = _reaction_table(p)
    limits = Limits(p.max_atoms, p.max_pairs_species, p.max_pairs_atom,
                    p.max_charge, _allowed(p))
    if any(m.size > p.max_atoms for m in seeds):
        raise ValueError(f"max_atoms={p.max_atoms} is smaller than an initial species")
    if limits.allowed is not None:
        outside = [m.formula for m in seeds if m.formula not in limits.allowed]
        if outside:
            raise ValueError(
                f"the initial species {outside} are not in allowed={sorted(limits.allowed)}; "
                "add them or pass allowed=[] to remove the constraint")

    found: dict[tuple, list] = {}

    def react(*mols):
        outcomes: dict[tuple, list] = {}
        left_key = tuple(sorted(m.id for m in mols))
        for products, cls in _outcomes(tuple(mols), classes, limits):
            key = tuple(sorted(m.id for m in products))
            if key == left_key:                     # an isomerisation that changes nothing
                continue
            entry = outcomes.setdefault(key, [products, []])
            if cls not in entry[1]:
                entry[1].append(cls)
        if not outcomes:
            return None
        left = tuple(sorted(m.id for m in mols))
        for key, (products, labels) in outcomes.items():
            rate = _rate(table, labels)
            found[(left, key)] = (rate, list(labels))
        return [entry[0] for entry in outcomes.values()]

    if p.method == "closure":
        return _closure(p, seeds, amounts, react, found, table)
    return _kinetic(p, rng, seeds, amounts, react, found, table)


def _rate(table, labels) -> dict | None:
    """Rate constants of the classes that give this reaction, summed (see notes)."""
    ks = [table[c][2] for c in labels if table[c][2] is not None]
    if not ks:
        return None
    return {"law": "mass-action", "k": float(sum(ks)),
            "classes": ",".join(labels), "units": table[labels[0]][3]}


def _network(p, species, reactions, status, amounts, seeds, table, extra) -> Network:
    by_id = {m.id: m for m in species}
    elements = sorted({el for m in species for el in m.atom_counts})
    conservation = [
        {"name": f"atoms of {el}", "vector": {i: int(m.atom_counts.get(el, 0)) for i, m in by_id.items()}}
        for el in elements
    ]
    conservation.append({"name": "electrons",
                         "vector": {i: int(m.electrons) for i, m in by_id.items()}})
    conservation.append({"name": "charge",
                         "vector": {i: int(m.charge) for i, m in by_id.items()}})
    initial = {m.id: float(a) for m, a in zip(seeds, amounts) if m.id in by_id}
    return Network(
        species=[Species(m.id, structure=m.id) for m in species],
        reactions=reactions,
        status=status,
        initial_state=initial or None,
        extras={
            "method": p.method,
            "templates": p.templates,
            "reaction_classes": [
                {"class": name, "name": cls, "template": text, "k": k, "units": units}
                for name, (cls, text, k, units) in table.items()
            ],
            "constraints": {"max_atoms_per_species": p.max_atoms,
                            "max_pairs_per_species": p.max_pairs_species,
                            "max_pairs_per_atom": p.max_pairs_atom},
            "formulas": {m.id: m.formula for m in species},
            "charges": {m.id: int(m.charge) for m in species},
            "electrons": {m.id: int(m.electrons) for m in species},
            "conservation": conservation,
            "species_encoding": "id = structure = canonical code of the molecular graph: "
                                "element + ':'*lone_pairs + '*'*radicals per atom, covalent bonds "
                                "by juxtaposition ('=' double, '#' triple), '.' for an ionic link, "
                                "branches in parentheses, ring closures '%n'",
            **extra,
        },
    )


def _closure(p, seeds, amounts, react, found, table) -> Network:
    species, pairs, status = expand(react, seeds, arity=[1, 2],
                                    max_species=p.max_species, ordered=False, alternatives=True)
    reactions, labels = [], []
    for lhs, rhs in pairs:
        rate, classes = found[(tuple(sorted(m.id for m in lhs)), tuple(sorted(m.id for m in rhs)))]
        reactions.append(Reaction(dict(Counter(m.id for m in lhs)),
                                  dict(Counter(m.id for m in rhs)), rate))
        labels.append(classes)
    return _network(p, species, reactions, status, amounts, seeds, table,
                    {"reaction_classes_used": labels,
                     "note": "deterministic network generator (DNG): every reaction object "
                             "reachable from the initial species under the size constraints"})


def _kinetic(p, rng, seeds, amounts, react, found, table) -> Network:
    """MC-sampling network generator: Gillespie's SSA over the molecules present."""
    if any(k is None for _, _, k, _ in table.values()):
        raise ValueError("method='kinetic' needs rate constants; templates='fig3' publishes none")

    known: dict[str, Mol] = {}
    counts: Counter = Counter()
    total = sum(amounts)
    for m, a in zip(seeds, amounts):
        known[m.id] = m
        counts[m.id] = max(1, int(round(p.molecules * a / total)))

    cache: dict[tuple, list] = {}
    fired: dict[tuple, list] = {}
    time = 0.0

    def channels(ids: tuple[str, ...]):
        if ids not in cache:
            out = []
            products = react(*[known[i] for i in ids])
            for group in products or ():
                for m in group:
                    known.setdefault(m.id, m)
                key = tuple(sorted(m.id for m in group))
                rate, classes = found[(tuple(sorted(ids)), key)]
                if rate is not None:
                    out.append((key, rate, classes))
            cache[ids] = out
        return cache[ids]

    for _ in range(p.steps):
        present = sorted(i for i, n in counts.items() if n > 0)
        active = []
        for a, i in enumerate(present):
            for channel in channels((i,)):
                active.append(((i,), channel, counts[i]))
            for j in present[a:]:
                multiplicity = (counts[i] * (counts[i] - 1) // 2 if i == j
                                else counts[i] * counts[j])
                if multiplicity <= 0:
                    continue
                for channel in channels(tuple(sorted((i, j)))):
                    active.append(((i, j), channel, multiplicity))
        weights = [channel[1]["k"] * mult for _, channel, mult in active]
        total_rate = sum(weights)
        if total_rate <= 0:
            break
        time += float(rng.exponential(1.0 / total_rate))
        pick = int(rng.choice(len(active), p=[w / total_rate for w in weights]))
        lhs, (rhs, rate, classes), _ = active[pick]
        for i in lhs:
            counts[i] -= 1
        for i in rhs:
            counts[i] += 1
        key = (tuple(sorted(lhs)), rhs)
        entry = fired.setdefault(key, [lhs, rhs, rate, classes, 0])
        entry[4] += 1

    species = [known[i] for i in sorted(known)]
    reactions = [
        Reaction(dict(Counter(lhs)), dict(Counter(rhs)), rate, count)
        for lhs, rhs, rate, _, count in fired.values()
    ]
    labels = [classes for _, _, _, classes, _ in fired.values()]
    return _network(p, species, reactions, "observed", amounts, seeds, table,
                    {"reaction_classes_used": labels,
                     "final_state": {i: int(n) for i, n in sorted(counts.items()) if n},
                     "elapsed_time": float(time),
                     "molecules": int(p.molecules),
                     "analysis": {"species_seen": len(species),
                                  "reactions_fired": sum(r.count for r in reactions)},
                     "note": "MC-sampling network generator (MCNG): Gillespie's SSA over the "
                             "molecules present; only species with a non-zero amount react, so "
                             "the recorded network is the kinetically relevant part of the closure"})
