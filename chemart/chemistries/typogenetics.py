"""Typogenetics (Hofstadter 1979; Morris 1989; Varetto 1993). Catalog id: typogenetics.

Strands are strings over A, C, G, T. A strand is translated duplet by duplet
into typoenzymes (AA is the punctuation between enzymes, an odd last base is
ignored). Each amino acid carries a kink (s, l, r); the fold of the enzyme
fixes the base it binds to. An enzyme then runs its amino acids on a complex
of two aligned rows: the lower strand (read left to right) and the upper,
complementary strand (read right to left), with gaps. When the enzymes are
done the rows fall apart into single strands, the daughters.

    self  (book 10.5.2, Hofstadter's puzzle):  s  -> daughters(s)
    pair  (Chemart addition):          g + s  -> g + daughters of s under g's enzymes

The machine follows the unambiguous specification of Snare (1999, Monash
honours thesis, chapter 2 and its Typogenetics.py), except that an inserting
amino acid moves the enzyme onto the inserted base, which is what Hofstadter's
worked example (GEB p. 508; see the tests) needs. method "closure" returns
every reaction reachable from the seed strands (chemart.expand.expand);
method "soup" draws molecules at random (chemart.soup.soup) and returns the
reactions that fired.
"""

from __future__ import annotations

import re
from bisect import bisect_right
from collections import Counter
from itertools import product as cartesian

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import soup

BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")

# (amino acid, kink) for duplets AA, AC, ..., TT; None is the AA punctuation.
# Hofstadter: GEB p. 510, book table 10.1 (left), Snare table 2.1(a).
# Varetto (1993, p. 187) per Snare table 2.1(b): row G carries inc, ing, int, ina
# while the kinks stay with the columns. Book table 10.1 prints the Hofstadter
# table twice (erratum).
_H = [None, ("cut", "s"), ("del", "s"), ("swi", "r"),
      ("mvr", "s"), ("mvl", "s"), ("cop", "r"), ("off", "l"),
      ("ina", "s"), ("inc", "r"), ("ing", "r"), ("int", "l"),
      ("rpy", "r"), ("rpu", "l"), ("lpy", "l"), ("lpu", "l")]
_V = _H[:8] + [("inc", "s"), ("ing", "r"), ("int", "r"), ("ina", "l")] + _H[12:]
CODE_TABLES = {
    name: {a + b: table[4 * i + j] for i, a in enumerate(BASES) for j, b in enumerate(BASES)}
    for name, table in (("hofstadter", _H), ("varetto", _V))
}
_TURN = {"s": 0, "r": 1, "l": -1}
# Direction of the last segment relative to the first (east), in right turns:
# east A, south G, west T, north C (GEB p. 511; Snare table 2.2).
_PREFERENCE = "AGTC"

_STRAND = re.compile(r"^[ACGT]+$")


def translate(strand: str, code_table: str = "hofstadter") -> list[list[tuple[str, str]]]:
    """Enzymes coded by a strand: lists of (amino acid, kink); empty enzymes are dropped."""
    code = CODE_TABLES[code_table]
    enzymes: list[list[tuple[str, str]]] = [[]]
    for i in range(0, len(strand) - 1, 2):
        acid = code[strand[i:i + 2]]
        if acid is None:
            enzymes.append([])
        else:
            enzymes[-1].append(acid)
    return [e for e in enzymes if e]


def enzyme(names: str, code_table: str = "hofstadter") -> list[tuple[str, str]]:
    """An enzyme written by name, e.g. enzyme("rpu-inc-cop"), with the kinks of the code table."""
    kinks = {acid: kink for acid, kink in filter(None, CODE_TABLES[code_table].values())}
    return [(a, kinks[a]) for a in names.split("-")]


def binding_preference(enzyme, rule: str = "hofstadter") -> str:
    """The base an enzyme binds to.

    hofstadter: the kinks between the first and the last segment count, i.e.
        those of all amino acids but the first and the last (Snare 2.4).
    morris: every amino acid's kink counts, as if hanging links preceded the
        first and followed the last (Morris 1989, Varetto 1993; book fig. 10.6).
    """
    kinks = [k for _, k in enzyme]
    if rule == "hofstadter":
        kinks = kinks[1:-1]
    return _PREFERENCE[sum(_TURN[k] for k in kinks) % 4]


def enzyme_text(enzyme, rule: str = "hofstadter") -> str:
    return "-".join(a for a, _ in enzyme) + ":" + binding_preference(enzyme, rule)


class Complex:
    """Two aligned rows (row 0 lower, row 1 upper) with gaps (None) and an enzyme on them."""

    def __init__(self, strand: str):
        self.rows = [list(strand), [None] * len(strand)]
        self.row = 0          # the row the enzyme is bound to
        self.pos = 0
        self.copy = False
        self.done = False

    def clone(self) -> Complex:
        c = Complex.__new__(Complex)
        c.rows = [self.rows[0][:], self.rows[1][:]]
        c.row, c.pos, c.copy, c.done = self.row, self.pos, self.copy, self.done
        return c

    def key(self) -> tuple:
        return (tuple(self.rows[0]), tuple(self.rows[1]), self.row, self.pos)

    # -- binding ------------------------------------------------------------
    def sites(self, base: str) -> list[int]:
        """Units holding `base` in the stretch of the bound row around the last position,
        in that row's reading direction (Snare, Applicator.apply)."""
        row = self.rows[self.row]
        blanks = [i for i, b in enumerate(row) if b is None]
        hi = bisect_right(blanks, self.pos)
        lo = blanks[hi - 1] if hi else -1
        up = blanks[hi] if hi < len(blanks) else len(row)
        out = [i for i in range(lo + 1, up) if row[i] == base]
        return out if self.row == 0 else out[::-1]

    def bind(self, site: int) -> None:
        self.pos, self.copy, self.done = site, False, False

    # -- the machine --------------------------------------------------------
    def _fill(self) -> None:
        lower, upper = self.rows
        if lower[self.pos] is None:
            lower[self.pos] = COMPLEMENT[upper[self.pos]]
        else:
            upper[self.pos] = COMPLEMENT[lower[self.pos]]

    def _moveto(self, pos: int, through_gaps: bool = False) -> None:
        self.pos = pos
        n = len(self.rows[0])
        if not 0 <= pos < n:
            self.done = True
        elif self.rows[self.row][pos] is None and (
                not through_gaps or self.rows[1 - self.row][pos] is None):
            self.done = True

    def _step(self, right: bool, through_gaps: bool = False) -> None:
        delta = 1 if right == (self.row == 0) else -1
        self._moveto(self.pos + delta, through_gaps)

    def _move(self, right: bool) -> None:
        self._step(right)
        if self.copy and not self.done:
            self._fill()

    def _search(self, right: bool, bases: frozenset) -> None:
        # Searches cross a gap in one row but not in both; bases filled in by
        # copy mode on the way do not stop the search (Morris; Snare 2.5).
        while True:
            self._step(right, through_gaps=True)
            if self.done:
                return
            found = self.rows[self.row][self.pos] in bases
            if self.copy:
                self._fill()
            if found:
                return

    def _insert(self, base: str) -> None:
        at = self.pos + 1 if self.row == 0 else self.pos      # right of the unit
        self.rows[self.row].insert(at, base)
        self.rows[1 - self.row].insert(at, COMPLEMENT[base] if self.copy else None)
        self.pos = at                                         # onto the new base (GEB p. 508)

    def run(self, acid: str) -> None:
        if acid == "cut":
            at = self.pos + 1 if self.row == 0 else self.pos
            self.rows[0].insert(at, None)
            self.rows[1].insert(at, None)
            if self.row == 1:
                self.pos += 1
        elif acid == "del":
            self.rows[self.row][self.pos] = None
            self._move(right=True)
        elif acid == "swi":
            self.row = 1 - self.row
            if self.rows[self.row][self.pos] is None:
                self.done = True
        elif acid in ("mvr", "mvl"):
            self._move(right=acid == "mvr")
        elif acid == "cop":
            self.copy = True
            self._fill()
        elif acid == "off":
            self.copy = False
        elif acid in ("ina", "inc", "ing", "int"):
            self._insert(acid[2].upper())
        else:  # rpy, rpu, lpy, lpu
            self._search(acid[0] == "r", PYRIMIDINES if acid[1:] == "py" else PURINES)

    def apply(self, enzyme, site: int) -> None:
        self.bind(site)
        for acid, _ in enzyme:
            self.run(acid)
            if self.done:
                break

    def strands(self) -> list[str]:
        """The single strands the complex falls apart into (upper ones read right to left)."""
        lower = "".join(b or " " for b in self.rows[0]).split()
        upper = "".join(b or " " for b in self.rows[1]).split()
        return lower + [s[::-1] for s in upper]


def apply_enzyme(enzyme, strand: str, site: int) -> list[str]:
    """Run one enzyme bound at unit `site` of a fresh strand; return the daughters."""
    c = Complex(strand)
    c.apply(enzyme, site)
    return c.strands()


class TooManyBranches(Exception):
    pass


def outcomes(enzymes, strand: str, fold: str, tiebreak: str, rng=None,
             max_branches: int = 1000) -> list[tuple[str, ...]]:
    """Distinct daughter multisets (sorted tuples) of running `enzymes` in order on `strand`.

    Each enzyme binds in the stretch where the previous one stopped (Snare 2.6);
    an enzyme with no site is skipped. tiebreak picks among equal sites: all
    (every choice, one outcome each), random, leftmost or rightmost (in the
    bound strand's reading direction; rightmost is Varetto's rule).
    """
    states = [Complex(strand)]
    for enzyme in enzymes:
        base = binding_preference(enzyme, fold)
        nxt: dict[tuple, Complex] = {}
        for c in states:
            sites = c.sites(base)
            if not sites:
                nxt.setdefault(c.key(), c)
                continue
            if tiebreak == "all":
                chosen = sites
            elif tiebreak == "random":
                chosen = [sites[int(rng.integers(len(sites)))]]
            else:
                chosen = [sites[0] if tiebreak == "leftmost" else sites[-1]]
            for i, site in enumerate(chosen):
                d = c if i == len(chosen) - 1 else c.clone()
                d.apply(enzyme, site)
                nxt.setdefault(d.key(), d)
            if len(nxt) > max_branches:
                raise TooManyBranches
        states = list(nxt.values())
    return sorted({tuple(sorted(c.strands())) for c in states})


# ---------------------------------------------------------------------------
def _check_strands(strands) -> list[str]:
    if not isinstance(strands, list) or not all(isinstance(s, str) and _STRAND.match(s) for s in strands):
        raise ValueError(f"strands must be a list of non-empty strings over A, C, G, T, got {strands!r}")
    return strands


class _Chemistry:
    def __init__(self, p, rng):
        self.p, self.rng = p, rng
        self.truncated = False
        self.cache: dict[tuple, list[tuple]] = {}

    def enzymes(self, strand):
        return translate(strand, self.p.code_table)

    def results(self, lhs, fresh: bool = False) -> list[tuple]:
        """Right-hand sides (complete, catalysts included) of the reactions of `lhs`."""
        if not fresh and lhs in self.cache:
            return self.cache[lhs]
        gene, target = (lhs[0], lhs[0]) if len(lhs) == 1 else lhs
        enzymes = self.enzymes(gene)
        out = []
        if enzymes:
            try:
                found = outcomes(enzymes, target, self.p.fold, self.p.binding_tiebreak,
                                 self.rng, self.p.max_branches)
            except TooManyBranches:
                self.truncated = True
                found = []
            for daughters in found:
                if any(len(s) > self.p.max_length for s in daughters):
                    self.truncated = True
                    continue
                rhs = daughters if len(lhs) == 1 else (gene, *daughters)
                if Counter(rhs) != Counter(lhs):
                    out.append(tuple(rhs))
        if not fresh:
            self.cache[lhs] = out
        return out

    def species(self, ids) -> list[Species]:
        rule = self.p.fold
        return [Species(s, structure="; ".join(enzyme_text(e, rule) for e in self.enzymes(s)) or None)
                for s in ids]


def generate(p, rng) -> Network:
    strands = _check_strands(p.strands)
    if not strands:
        strands = ["".join(BASES[int(i)] for i in rng.integers(4, size=p.strand_length))
                   for _ in range(p.n_random)]
    if any(len(s) > p.max_length for s in strands):
        raise ValueError(f"seed strands must not be longer than max_length={p.max_length}")
    chem = _Chemistry(p, rng)
    arity = 1 if p.reaction == "self" else 2
    if p.method == "closure":
        return _closure(p, chem, strands, arity)
    if p.binding_tiebreak == "all":
        raise ValueError("method 'soup' needs one outcome per collision: use binding_tiebreak "
                         "random, leftmost or rightmost (or method 'closure' for 'all')")
    return _soup(p, chem, strands, arity, rng)


def _closure(p, chem, strands, arity) -> Network:
    def discover(*lhs):
        found = [s for rhs in chem.results(lhs) for s in rhs]
        return tuple(dict.fromkeys(found)) if found else None

    seed = list(dict.fromkeys(strands))
    ids, _, status = expand(discover, seed, arity=arity, max_species=p.max_species, ordered=True)
    known = set(ids)
    reactions = []
    combos = ((s,) for s in ids) if arity == 1 else cartesian(ids, repeat=2)
    for lhs in combos:
        for rhs in chem.results(tuple(lhs)):
            if set(rhs) <= known:
                reactions.append(Reaction.of(lhs, rhs))
            else:
                status = "truncated"
    if chem.truncated:
        status = "truncated"
    return Network(species=chem.species(ids), reactions=reactions, status=status,
                   extras={"seed": seed})


def _soup(p, chem, strands, arity, rng) -> Network:
    def react(*lhs):
        found = chem.results(lhs, fresh=p.binding_tiebreak == "random")
        return found[0] if found else None

    start = [s for s in strands for _ in range(p.copies)]
    if len(start) < arity:
        raise ValueError(f"the soup needs at least {arity} molecules, got {len(start)}")
    fired, final = soup(react, start, p.steps, rng, arity=arity, dilution="constant")
    ids = list(dict.fromkeys([*start, *(s for _, rhs, _ in fired for s in rhs)]))
    return Network(
        species=chem.species(ids),
        reactions=[Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in fired],
        status="observed",
        initial_state={s: c for s, c in Counter(start).items()},
        outflow=CONSTANT_TOTAL,
        extras={"final_state": dict(Counter(final).most_common())},
    )
