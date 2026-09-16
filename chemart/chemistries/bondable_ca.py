"""Bondable cellular automata (BCA): Hatcher, Banzhaf & Yu (ECAL 2011). Book 10.7.3.

A *subsymbolic* artificial chemistry: reactivity is not declared, it is read off
an observable of an internal dynamical system. An atom is a one-dimensional
binary cellular automaton of width `w` with periodic (closed) boundaries, and
its *type* is the Wolfram rule number 0..255 that runs it. A molecule is a set
of bonded atoms.

Three things are read off the dynamics:

- **Mean polarity.** Each cell contributes +1 when it is 1 and -1 when it is 0,
  so a configuration has polarity `ones - zeros` in `-w .. +w` (the scale of
  book figure 10.16). A CA on `w` cells is a finite deterministic map, so its
  trajectory settles onto a cycle; the *settled mean polarity* is the mean of
  the polarity over that cycle and the *settling iteration* is the length of
  the transient before it ("the transition time of the dynamical system made up
  of that CA", caption of fig. 10.16). Ordering the 256 rules by this value is
  the "periodic table" of figure 10.16, which `periodic_table()` reproduces.
- **Bonding strength**, the "green bars" of book figure 10.15: the longest
  contiguous run of complementary cells between the two atoms, maximised over
  the `w` cyclic alignments of one ring against the other. It is symmetric.
- **Coupling.** "Once bonded, the neighborhood in each CA is influenced by the
  other CA, so as to potentially change the overall dynamics" (book 10.7.3):
  inside the bonded run a cell takes the neighbour that faces the partner from
  the partner's aligned cell instead of from its own ring. That can strengthen
  or weaken the bond, which is what makes molecules decay.

Reactions are the association and dissociation events these three produce, and
the network is the closure (`chemart.expand.expand`) of a seed set of atom
types under them.

The CA step is done bitwise on the whole ring at once (eight neighbourhood
patterns, independent of `w`), so settling a molecule is cheap.
"""

from __future__ import annotations

from collections import Counter

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

RULES = 256                 # Wolfram's indexing of 1-d binary nearest-neighbour CAs
SUM_TOLERANCE = 0.5         # tolerance of the "sum-zero" polarity criterion
COUPLINGS = ("replace", "xor", "none")
POLARITY_RULES = ("opposite-sign", "sum-zero", "any")


# ============================================================================
# Ring arithmetic
# ============================================================================
def _rotl(x: int, k: int, w: int) -> int:
    """Bit i of the result is bit i - k of x."""
    k %= w
    return x if k == 0 else ((x << k) | (x >> (w - k))) & ((1 << w) - 1)


def _rotr(x: int, k: int, w: int) -> int:
    """Bit i of the result is bit i + k of x."""
    return _rotl(x, -k, w)


def polarity(config: int, w: int) -> int:
    """ones - zeros of one configuration: the -w .. +w scale of figure 10.16."""
    return 2 * config.bit_count() - w


def initial_config(w: int) -> int:
    """The balanced start: floor(w / 2) contiguous ones on the ring, the rest zeros."""
    return (1 << (w // 2)) - 1


def bits(config: int, w: int) -> str:
    """Cell 0 first."""
    return "".join(str((config >> i) & 1) for i in range(w))


# ============================================================================
# The cellular automaton
# ============================================================================
def step(config: int, rule: int, w: int, left: tuple[int, int] = (0, 0),
         right: tuple[int, int] = (0, 0), coupling: str = "replace") -> int:
    """One synchronous CA step on a ring of `w` cells.

    `left` and `right` are ``(mask, partner_bits)`` pairs aligned to this ring's
    cell indices: on the cells selected by the mask, that side of the
    neighbourhood comes from the partner instead of from this ring (`replace`)
    or is exclusive-ored with it (`xor`).
    """
    mask = (1 << w) - 1
    a = _rotl(config, 1, w)          # a's bit i is cell i - 1
    b = _rotr(config, 1, w)          # b's bit i is cell i + 1
    if coupling != "none":
        for side, (m, partner) in (("l", left), ("r", right)):
            if not m:
                continue
            value = a if side == "l" else b
            value = (value & ~m & mask) | (partner & m) if coupling == "replace" else value ^ (partner & m)
            if side == "l":
                a = value
            else:
                b = value
    out = 0
    for pattern in range(8):
        if not (rule >> pattern) & 1:
            continue
        sel = a if pattern & 4 else ~a
        sel &= config if pattern & 2 else ~config
        sel &= b if pattern & 1 else ~b
        out |= sel
    return out & mask


def _trajectory(advance, start, cap: int):
    """Iterate `advance` from `start` until a state repeats or `cap` steps are used.

    Returns (states to average over, settling iteration, settled state, settled?).
    A trajectory that has not closed within `cap` steps is reported unsettled and
    averaged over its second half.
    """
    seen: dict = {}
    seq: list = []
    s = start
    while s not in seen and len(seq) < cap:
        seen[s] = len(seq)
        seq.append(s)
        s = advance(s)
    if s in seen:
        return seq[seen[s]:], seen[s], s, True
    return seq[len(seq) // 2:], len(seq), s, False


def settle(rule: int, w: int, config: int | None = None, cap: int = 2048) -> dict:
    """Run one free atom to its attractor: the observable of figure 10.16."""
    start = initial_config(w) if config is None else config
    states, transient, final, ok = _trajectory(lambda s: step(s, rule, w), start, cap)
    return {
        "mean_polarity": sum(polarity(x, w) for x in states) / len(states),
        "settling_iteration": transient,
        "cycle_length": len(states) if ok else 0,
        "config": final,
        "settled": ok,
    }


def periodic_table(w: int, cap: int = 2048) -> list[dict]:
    """The 256 rules ordered by settled mean polarity: book figure 10.16."""
    rows = []
    for rule in range(RULES):
        r = settle(rule, w, cap=cap)
        rows.append({
            "rule": rule,
            "mean_polarity": r["mean_polarity"],
            "settling_iteration": r["settling_iteration"],
            "cycle_length": r["cycle_length"],
            "settled_config": bits(r["config"], w),
        })
    return sorted(rows, key=lambda r: (r["mean_polarity"], r["rule"]))


# ============================================================================
# Bonding strength: the longest contiguous complementary run (figure 10.15)
# ============================================================================
def _longest_run(x: int, w: int) -> tuple[int, int]:
    """(length, mask) of the longest contiguous cyclic run of set bits in `x`."""
    full = (1 << w) - 1
    if x == full:
        return w, full
    if x == 0:
        return 0, 0
    zero = next(i for i in range(w) if not (x >> i) & 1)
    best_len = best_start = cur = cur_start = 0
    for k in range(w):
        i = (zero + 1 + k) % w
        if (x >> i) & 1:
            if cur == 0:
                cur_start = i
            cur += 1
            if cur > best_len:
                best_len, best_start = cur, cur_start
        else:
            cur = 0
    mask = 0
    for t in range(best_len):
        mask |= 1 << ((best_start + t) % w)
    return best_len, mask


def run_at(a: int, b: int, w: int, offset: int) -> tuple[int, int]:
    """(length, mask) of the longest complementary run when cell i of `a` faces cell i + offset of `b`."""
    return _longest_run(a ^ _rotr(b, offset, w), w)


def best_bond(a: int, b: int, w: int) -> tuple[int, int, int]:
    """(strength, offset, mask): the strongest alignment of two atoms.

    Strength is the longest contiguous run of complementary cells; ties go to
    the smallest offset, so the choice is deterministic.
    """
    best = (0, 0, 0)
    for offset in range(w):
        length, mask = run_at(a, b, w, offset)
        if length > best[0]:
            best = (length, offset, mask)
    return best


def bond_strength(a: int, b: int, w: int) -> int:
    """Bonding strength of two atom configurations; symmetric in a and b."""
    return best_bond(a, b, w)[0]


# ============================================================================
# Molecules: a chain of bonded atoms
# ============================================================================
class Molecule:
    """A chain of atoms: rules, current configurations and one offset per bond.

    The bond between atom k and atom k + 1 stores only its alignment `offset`;
    the bonded run itself is recomputed from the current configurations at every
    step, so a bond can grow or shrink while the coupled dynamics run.
    """

    __slots__ = ("rules", "configs", "offsets", "w")

    def __init__(self, rules, configs, offsets, w):
        self.rules, self.configs = tuple(rules), tuple(configs)
        self.offsets, self.w = tuple(offsets), w

    def __len__(self) -> int:
        return len(self.rules)

    def reversed(self) -> Molecule:
        w = self.w
        return Molecule(self.rules[::-1], self.configs[::-1],
                        tuple((-d) % w for d in reversed(self.offsets)), w)

    def code(self) -> str:
        w = self.w
        parts = [f"r{self.rules[0]}:{bits(self.configs[0], w)}"]
        for k, d in enumerate(self.offsets):
            parts += [f"+{d}", f"r{self.rules[k + 1]}:{bits(self.configs[k + 1], w)}"]
        return "~".join(parts)

    def canonical(self) -> tuple[str, str, Molecule]:
        """(code, readable name, this molecule in its canonical orientation)."""
        other = self.reversed()
        best = self if self.code() <= other.code() else other
        return best.code(), "-".join(f"r{r}" for r in best.rules), best

    def runs(self) -> list[tuple[int, int]]:
        """(length, mask) of every bond's current run, in the left atom's coordinates."""
        return [run_at(self.configs[k], self.configs[k + 1], self.w, d)
                for k, d in enumerate(self.offsets)]

    def polarity(self) -> float:
        """ones - zeros averaged over the atoms, on the same -w .. +w scale."""
        return sum(polarity(c, self.w) for c in self.configs) / len(self.configs)


def step_molecule(m: Molecule, coupling: str) -> tuple[int, ...]:
    """One synchronous step of every atom, with bonded cells coupled to their partner."""
    w, n = m.w, len(m)
    sides = [[(0, 0), (0, 0)] for _ in range(n)]      # per atom: [left, right]
    if coupling != "none":
        for k, d in enumerate(m.offsets):
            _, mask = run_at(m.configs[k], m.configs[k + 1], w, d)
            if not mask:
                continue
            # atom k's cell i faces atom k + 1's cell i + d, and vice versa.
            sides[k][1] = (mask, _rotr(m.configs[k + 1], d, w))
            sides[k + 1][0] = (_rotl(mask, d, w), _rotl(m.configs[k], d, w))
    return tuple(
        step(m.configs[k], m.rules[k], w, sides[k][0], sides[k][1], coupling)
        for k in range(n)
    )


def settle_molecule(m: Molecule, coupling: str, cap: int = 2048) -> tuple[Molecule, int, bool]:
    """Run the coupled molecule to its attractor: (settled molecule, settling iteration, settled?)."""
    def advance(state):
        return step_molecule(Molecule(m.rules, state, m.offsets, m.w), coupling)

    _, transient, final, ok = _trajectory(advance, m.configs, cap)
    return Molecule(m.rules, final, m.offsets, m.w), transient, ok


# ============================================================================
# The chemistry: association and dissociation
# ============================================================================
def _polarity_ok(criterion: str, p1: float, p2: float) -> bool:
    if criterion == "any":
        return True
    if criterion == "opposite-sign":
        return p1 * p2 < 0
    return abs(p1 + p2) <= SUM_TOLERANCE


class Chemistry:
    """The reaction rules, over canonical molecule codes."""

    def __init__(self, p):
        self.w, self.threshold = p.w, p.bond_threshold
        self.coupling, self.criterion = p.coupling, p.polarity_rule
        self.max_atoms, self.cap = p.max_atoms, p.settle_iterations
        self.molecules: dict[str, Molecule] = {}
        self.names: Counter = Counter()
        self.ids: dict[str, str] = {}
        self.unsettled: set[str] = set()

    # -- registry ----------------------------------------------------------
    def register(self, m: Molecule, settled: bool = True) -> str:
        code, name, canon = m.canonical()
        if code not in self.molecules:
            self.molecules[code] = canon
            self.names[name] += 1
            self.ids[code] = name if self.names[name] == 1 else f"{name}#{self.names[name]}"
        if not settled:
            self.unsettled.add(code)
        return code

    def atom(self, rule: int) -> str:
        r = settle(rule, self.w, cap=self.cap)
        return self.register(Molecule([rule], [r["config"]], [], self.w), r["settled"])

    def settled(self, m: Molecule) -> str:
        """Settle a molecule under coupling and register the result."""
        s, _, ok = settle_molecule(m, self.coupling, self.cap)
        return self.register(s, ok)

    # -- reactions ---------------------------------------------------------
    def associate(self, c1: str, c2: str):
        m1, m2 = self.molecules[c1], self.molecules[c2]
        if len(m1) + len(m2) > self.max_atoms:
            return None
        if not _polarity_ok(self.criterion, m1.polarity(), m2.polarity()):
            return None
        best = None
        for a in (m1, m1.reversed()):
            for b in (m2, m2.reversed()):
                strength, offset, _ = best_bond(a.configs[-1], b.configs[0], self.w)
                joined = Molecule(a.rules + b.rules, a.configs + b.configs,
                                  a.offsets + (offset,) + b.offsets, self.w)
                key = (-strength, joined.canonical()[0])
                if best is None or key < best[0]:
                    best = (key, strength, joined)
        _, strength, joined = best
        if strength < self.threshold:
            return None
        return (self.settled(joined),)

    def dissociate(self, code: str):
        """A molecule decays when the coupled dynamics leave one of its bonds too weak."""
        m = self.molecules[code]
        if len(m) < 2:
            return None
        weak = [k for k, (length, _) in enumerate(m.runs()) if length < self.threshold]
        if not weak:
            return None
        k = weak[0]
        left = Molecule(m.rules[:k + 1], m.configs[:k + 1], m.offsets[:k], m.w)
        right = Molecule(m.rules[k + 1:], m.configs[k + 1:], m.offsets[k + 1:], m.w)
        return (self.settled(left), self.settled(right))

    def react(self, *codes):
        return self.dissociate(*codes) if len(codes) == 1 else self.associate(*codes)


# ============================================================================
def _atom_types(value) -> list[int]:
    if (not isinstance(value, list) or not value
            or not all(isinstance(v, int) and not isinstance(v, bool) for v in value)):
        raise ValueError(f"atom_types must be a non-empty list of integers, got {value!r}")
    if any(not 0 <= v < RULES for v in value):
        raise ValueError(f"atom_types must be Wolfram rule numbers 0..{RULES - 1}, got {value!r}")
    return list(dict.fromkeys(value))


def generate(p, rng):
    if p.bond_threshold > p.w:
        raise ValueError(
            f"bond_threshold must be at most the CA width w, got {p.bond_threshold} with w = {p.w}"
        )
    types = _atom_types(p.atom_types)

    chem = Chemistry(p)
    seed = [chem.atom(rule) for rule in types]
    codes, found, status = expand(chem.react, seed, arity=(1, 2),
                                  max_species=p.max_species, ordered=True)

    ids = chem.ids
    species = [Species(ids[c], structure=c) for c in codes]
    reactions = [Reaction.of([ids[c] for c in lhs], [ids[c] for c in rhs]) for lhs, rhs in found]

    counts = {c: Counter(chem.molecules[c].rules) for c in codes}
    conservation = [
        {"name": f"atoms of type r{rule}",
         "vector": {ids[c]: int(counts[c][rule]) for c in codes}}
        for rule in types if any(counts[c][rule] for c in codes)
    ]
    return Network(
        species=species,
        reactions=reactions,
        status=status,
        extras={
            "conservation": conservation,
            "species_encoding": (
                "id = the chain of atom types (r<rule>) joined by '-', with '#k' when two "
                "molecules share a skeleton but not their cell states; structure = "
                "r<rule>:<cells>~+<offset>~r<rule>:<cells>, cell 0 first, the offset being the "
                "alignment of the next atom's ring against this one"
            ),
            "interaction_law": {
                "name": "bondable cellular automata (Hatcher, Banzhaf & Yu)",
                "atom": f"a 1-d binary CA of width {p.w} with periodic boundaries, run by its "
                        f"Wolfram rule; every atom starts from {bits(initial_config(p.w), p.w)}",
                "observable": "mean polarity = (ones - zeros) averaged over the attractor cycle, "
                              f"on the scale -{p.w} .. +{p.w} of book figure 10.16",
                "bond_strength": "the longest contiguous run of complementary cells between the "
                                 "two rings, maximised over the w cyclic alignments (the green "
                                 "bars of book figure 10.15); symmetric",
                "association": f"two molecules bond end to end when their mean polarities satisfy "
                               f"'{p.polarity_rule}' and the strongest alignment of the two facing "
                               f"atoms has strength >= {p.bond_threshold}",
                "coupling": {
                    "mode": p.coupling,
                    "law": "inside the bonded run a cell takes the neighbour that faces the "
                           "partner from the partner's aligned cell instead of from its own ring "
                           "(replace), or exclusive-ors the two (xor); every other cell is "
                           "unaffected, and the run is recomputed from the current cells at every "
                           "step, so a bond can grow or shrink",
                },
                "dissociation": f"a molecule whose coupled dynamics have left a bond with a run "
                                f"shorter than {p.bond_threshold} splits at that bond, and each "
                                f"fragment settles again",
                "settle_iterations": p.settle_iterations,
            },
            "analysis": {
                "periodic_table": periodic_table(p.w, p.settle_iterations),
                "periodic_table_scale": (
                    "settled mean polarity in -w .. +w, rules ordered by it: book figure 10.16 "
                    "(drawn there for w = 12)"
                ),
                "molecule_polarity": {ids[c]: chem.molecules[c].polarity() for c in codes},
                "bond_strengths": {
                    ids[c]: [int(length) for length, _ in chem.molecules[c].runs()]
                    for c in codes if len(chem.molecules[c]) > 1
                },
                "largest_molecule_atoms": max(len(chem.molecules[c]) for c in codes),
                "unsettled_molecules": sorted(ids[c] for c in chem.unsettled),
            },
        },
    )
