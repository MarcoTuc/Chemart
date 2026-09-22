"""Fraglets: Tschudin (2003), book 16.2.1 and 17.3.1. Catalog id: fraglets.

A fraglet is a string of symbols [s1 s2 ... sn] stored in a network node
(a well-stirred vessel). The head symbols are the instruction, the rest is
data. Transformations rewrite one fraglet; match/matchp join two fraglets
whose tags agree, and matchp survives the reaction (a catalyst). `send`
moves the tail of a fraglet to another node.

Two dialects, each ported from its source:

- "pycellchem": the book's reference interpreter (PyCellChemistry
  Fraglets.py / NetFraglets.py, L. Yamamoto 2013) and book table 16.1:
  match matchp dup exch pop nop nul split send fork, with ``node[send dest
  tail] -> dest[tail]``.
- "fraglets-2007": the upstream C interpreter (fraglets0.32 src/frag.c,
  C. Tschudin et al.) and its instruction set document (2007-09-24): the
  core set, ``[send seg dest tail] -> dest[tail]`` over segments, and the
  logic and arithmetic instructions (empty eq lt length sum sub mult div mod
  pow abs pop pop2 copy).

A program uses the upstream .fra syntax: ``f node[s1 s2 ...]mult`` puts
fraglets in a node, ``a node seg`` attaches a node to a segment, ``#``
starts a comment and ``e`` ends the input. Symbols are separated by spaces
or ':' as in the AINS paper.

Two faces:
- ``generate``: chemart.expand.expand of the distinct program fraglets
  (status complete or truncated).
- ``evolve``: the PyCellChemistry scheduler (Fraglets.propensity/react and
  Cell.gillespie): all pending transformations fire at once, then one
  match is chosen with propensity n_active * n_passive, over all nodes; a
  frame after every match (status observed).
"""

from __future__ import annotations

import re
from collections import Counter

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame

DIALECTS = ("pycellchem", "fraglets-2007")
MATCH_OPS = frozenset({"match", "matchp"})

#: PyCellChemistry Fraglets.instr, minus the two match rules.
PYCELLCHEM_TRANSFORMS = frozenset({"dup", "exch", "pop", "nop", "nul", "split", "send", "fork"})

#: frag.c f_execFromQueue transformations that are ported.
FRAGLETS_2007_TRANSFORMS = frozenset({
    "nop", "copy", "fork", "dup", "empty", "exch", "lt", "eq", "length",
    "sum", "sub", "mult", "div", "mod", "pow", "abs", "nul", "pop", "pop2",
    "send", "split",
})
#: frag.c instructions that depend on time, node creation, printing or the
#: host application; a program using them is rejected.
FRAGLETS_2007_UNSUPPORTED = frozenset({
    "delay", "wait", "newname", "newnode", "printsym", "broadcast", "anycast",
    "inject", "expel", "splitat", "_", "out",
})
ARITHMETIC = frozenset({"sum", "sub", "mult", "div", "mod", "pow"})

STAR = "*"
#: A cascade of transformations longer than this within one step stops the run.
MAX_TRANSFORMATIONS = 1_000_000

_NUMBER = re.compile(r"^-?\d+")
_FRAGLET = re.compile(r"^([^\[\]\s]*)\s*\[([^\[\]]*)\]\s*(\S+)?\s*(#.*)?$")

#: The active confirmed delivery protocol of Tschudin (2003) sec. IV.B and book figure 16.10.
CDP = """\
a a net
a b net
f a[matchp cdp send b split send a ack *]
f a[cdp data]
"""


# --- molecules ---------------------------------------------------------------
def species_id(mol) -> str:
    node, syms = mol
    return f"{node}[{','.join(syms)}]"


def structure(mol) -> str:
    node, syms = mol
    return f"{node}[{' '.join(syms)}]"


def transforms(dialect: str) -> frozenset:
    return PYCELLCHEM_TRANSFORMS if dialect == "pycellchem" else FRAGLETS_2007_TRANSFORMS


def is_transformation(syms, dialect: str) -> bool:
    return bool(syms) and syms[0] in transforms(dialect)


def is_match(syms) -> bool:
    return bool(syms) and syms[0] in MATCH_OPS


def is_passive(syms, dialect: str) -> bool:
    return bool(syms) and not is_match(syms) and not is_transformation(syms, dialect)


def keep(syms, dialect: str) -> bool:
    """Whether a new fraglet enters the store (Fraglets.inject / frag.c rx_fraglet)."""
    if not syms:
        return False
    if is_match(syms):
        return len(syms) > 1
    if dialect == "pycellchem" and is_transformation(syms, dialect):
        return len(syms) > 1
    return True


# --- transformations -----------------------------------------------------------
def _pycellchem(node, s, links):
    """PyCellChemistry Fraglets.r_* (one symbol per character there, whole symbols here)."""
    op, n = s[0], len(s)
    if op == "nul":
        return []
    if n < 2:
        return []
    if op == "nop":
        return [s[1:]]
    if op == "dup":
        return [s[1:2]] if n < 3 else [(s[1], s[2]) + s[2:]]
    if op == "exch":
        return [s[1:]] if n < 4 else [(s[1], s[3], s[2]) + s[4:]]
    if op == "pop":
        return [s[1:2]] if n < 4 else [(s[1],) + s[3:]]
    if op == "fork":
        return [s[1:]] if n < 3 else [(s[1],) + s[3:], s[2:]]
    if op == "split":
        return _split(s)
    if op == "send":
        if n < 3 or s[1] not in links.get(node, ()):
            return []
        return [(s[1], s[2:])]
    raise AssertionError(op)


def _split(s):
    tail = s[1:]
    if STAR not in tail:
        return [tail]
    i = tail.index(STAR)
    return [tail[:i], tail[i + 1:]]


def _isnum(sym: str) -> bool:
    return bool(re.fullmatch(r"-?\d+", sym))


def _fraglets2007(node, s, segments, nodes):
    """frag.c f_nop, f_copy, f_fork, f_dup, f_empty, f_exch, f_lt, f_eq, f_length,
    f_op2, f_abs, f_nul, f_pop, f_pop2, f_send, f_split. A dropped fraglet gives []."""
    op, n = s[0], len(s)
    if op == "nul":
        return []
    if op == "nop":
        return [s[1:]]
    if op == "split":
        return _split(s)
    if op == "copy":
        return [] if n < 2 else [s[1:], s[1:]]
    if op == "fork":
        return [] if n < 3 else [(s[1],) + s[3:], s[2:]]
    if op == "dup":
        return [] if n < 3 else [(s[1], s[2]) + s[2:]]
    if op == "empty":
        if n < 3:
            return []
        return [(s[1],)] if n == 3 else [(s[2],) + s[3:]]
    if op == "exch":
        return [] if n < 4 else [(s[1], s[3], s[2]) + s[4:]]
    if op in ("lt", "eq"):
        if n < 5:
            return []
        if op == "lt":
            if not (_isnum(s[3]) and _isnum(s[4])):
                return []
            yes = int(s[3]) < int(s[4])
        else:
            yes = s[3] == s[4]
        return [((s[1] if yes else s[2]),) + s[3:]]
    if op == "length":
        return [] if n < 2 else [(s[1], str(n - 2)) + s[2:]]
    if op in ARITHMETIC:
        if n < 4 or not (_isnum(s[2]) and _isnum(s[3])):
            return []
        a, b = int(s[2]), int(s[3])
        if op == "sum":
            v = a + b
        elif op == "sub":
            v = a - b
        elif op == "mult":
            v = a * b
        elif op in ("div", "mod"):
            if b == 0:
                return []
            q = abs(a) // abs(b) * (1 if (a < 0) == (b < 0) else -1)   # C: truncate toward 0
            v = q if op == "div" else a - b * q
        else:
            if b < 0 or (a == 0 and b == 0):
                return []
            v = a ** b
        return [(s[1], str(v)) + s[4:]]
    if op == "abs":
        if n < 3 or not _isnum(s[2]):
            return []
        return [(s[1], str(abs(int(s[2])))) + s[3:]]
    if op == "pop":
        return [] if n < 3 else [(s[1],) + s[3:]]
    if op == "pop2":
        if n < 3:
            return []
        if n < 4:
            return [(s[1],), (s[2],)]
        return [(s[1], s[3]), (s[2],) + s[4:]]
    if op == "send":
        if n < 3 or s[1] in ("stdout", "stderr"):
            return []
        if n < 4:
            return []
        dest, seg = s[2], s[1]
        if dest not in nodes or node not in segments.get(seg, ()) or dest not in segments[seg]:
            return []
        return [(dest, s[3:])]
    raise AssertionError(op)


class Chemistry:
    """The reaction rules of one dialect over one network topology."""

    def __init__(self, dialect: str, segments: dict[str, list[str]], nodes=()):
        if dialect not in DIALECTS:
            raise ValueError(f"dialect must be one of {DIALECTS}, got {dialect!r}")
        self.dialect = dialect
        self.segments = {k: list(v) for k, v in segments.items()}
        self.nodes = set(nodes) | {x for v in segments.values() for x in v}
        # PyCellChemistry add_cnx: a node reaches the other nodes on its segments.
        self.links = {
            x: {y for v in segments.values() if x in v for y in v if y != x}
            for x in self.nodes
        }

    def transform(self, mol):
        """Products of a transformation, or None if mol is not one."""
        node, s = mol
        if not is_transformation(s, self.dialect):
            return None
        if self.dialect == "pycellchem":
            out = _pycellchem(node, s, self.links)
        else:
            out = _fraglets2007(node, s, self.segments, self.nodes)
        mols = []
        for p in out:
            if s[0] == "send":
                dest, syms = p
            else:
                dest, syms = node, p
            if keep(syms, self.dialect):
                mols.append((dest, tuple(syms)))
        return mols

    def match(self, a, b):
        """Products of [match(p) t tail1] + [t tail2] in the same node, or None."""
        (na, sa), (nb, sb) = a, b
        if na != nb or not is_match(sa) or len(sa) < 2 or not is_passive(sb, self.dialect):
            return None
        if sa[1] != sb[0]:
            return None
        joined = sa[2:] + sb[1:]
        out = [a] if sa[0] == "matchp" else []
        if keep(joined, self.dialect):
            out.append((na, joined))
        return out

    def react(self, *mols):
        return self.transform(mols[0]) if len(mols) == 1 else self.match(*mols)


# --- programs ------------------------------------------------------------------
def parse_program(text: str, dialect: str):
    """Parse .fra text: returns (initial Counter of molecules, segments {seg: [nodes]}, nodes)."""
    if not isinstance(text, str):
        raise ValueError(f"program must be a string, got {text!r}")
    initial: Counter = Counter()
    segments: dict[str, list[str]] = {}
    nodes: list[str] = []
    ended = False
    for lno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if ended:
            raise ValueError(f"program line {lno}: nothing may follow the end marker 'e': {raw!r}")
        head, rest = line[0], line[1:].strip()
        if head == "e" and not rest:
            ended = True
        elif head in "fF":
            m = _FRAGLET.match(rest)
            if not m:
                raise ValueError(f"program line {lno}: expected 'f node[symbols]mult', got {raw!r}")
            node, body, mult = m.group(1), m.group(2), m.group(3)
            if mult is not None and not mult.isdigit():
                raise ValueError(f"program line {lno}: multiplicity must be a positive integer, got {mult!r}")
            count = 1 if mult is None else int(mult)
            if count < 1:
                raise ValueError(f"program line {lno}: multiplicity must be > 0")
            syms = tuple(_symbol(tok, dialect, lno) for tok in re.split(r"[\s:]+", body) if tok)
            if node not in nodes:
                nodes.append(node)
            if keep(syms, dialect):
                initial[(node, syms)] += count
        elif head in "aA":
            parts = rest.split("#")[0].split()
            if len(parts) != 2:
                raise ValueError(f"program line {lno}: expected 'a node segment', got {raw!r}")
            node, seg = parts
            if node not in nodes:
                nodes.append(node)
            if node not in segments.setdefault(seg, []):
                segments[seg].append(node)
        else:
            raise ValueError(f"program line {lno}: lines start with 'f', 'a', 'e' or '#', got {raw!r}")
    if not initial:
        raise ValueError("program contains no fraglet (lines 'f node[symbols]')")
    return initial, segments, nodes


def _symbol(tok: str, dialect: str, lno: int) -> str:
    if "," in tok:
        raise ValueError(f"program line {lno}: symbol {tok!r} may not contain ','")
    if dialect == "fraglets-2007":
        if tok in FRAGLETS_2007_UNSUPPORTED:
            raise ValueError(
                f"program line {lno}: instruction {tok!r} is not supported "
                f"(timing, node creation, printing and application hooks are not ported)")
        if tok[0].isdigit() or (tok[0] == "-" and len(tok) > 1 and tok[1].isdigit()):
            return str(int(_NUMBER.match(tok).group(0)))     # frag.c name_add: atol
    return tok


# --- runs ------------------------------------------------------------------------
def closure(chem: Chemistry, seeds, max_species: int):
    return expand(chem.react, seeds, arity=(1, 2), max_species=max_species, ordered=True)


def ssa(chem: Chemistry, initial: Counter, steps: int, rng, tally: Tally, info: dict):
    """PyCellChemistry scheduling, as a generator of recording points.

    Yields (matches so far, live population) once the program's transformations
    have settled and again after every match and the transformations it sets
    off. `tally` records the reactions that fire; `info` is filled with the
    run's summary (bimolecular_events, inert, stopped, elastic) when it ends.
    """
    pop = Counter(initial)
    elastic: set = set()
    info.update({"bimolecular_events": 0, "inert": False, "stopped": None})

    def settle():
        done = 0
        while True:
            pending = [m for m, c in pop.items() if c > 0 and m not in elastic
                       and is_transformation(m[1], chem.dialect)]
            if not pending:
                return True
            for m in pending:
                c = pop.pop(m)
                out = chem.transform(m)
                if Counter(out) == Counter([m]):
                    elastic.add(m)          # e.g. [dup dup dup]: the reference would loop forever
                    pop[m] += c
                    continue
                for p in out:
                    pop[p] += c
                tally.add((m,), tuple(out), c)
                done += c
                if done > MAX_TRANSFORMATIONS:
                    info["stopped"] = f"more than {MAX_TRANSFORMATIONS} transformations in one step"
                    return False

    def propensities():
        pairs, weights = [], []
        passive: dict = {}
        for m, c in pop.items():
            if c > 0 and is_passive(m[1], chem.dialect):
                passive.setdefault((m[0], m[1][0]), []).append(m)
        for m, c in pop.items():
            if c > 0 and is_match(m[1]):
                for b in passive.get((m[0], m[1][1]), ()):
                    pairs.append((m, b))
                    weights.append(c * pop[b])
        return pairs, weights

    ok = settle()
    yield 0, pop
    for _ in range(steps if ok else 0):
        pairs, weights = propensities()
        total = sum(weights)
        if total <= 0:
            break
        w = rng.random() * total
        for (a, b), wt in zip(pairs, weights):
            if w < wt:
                break
            w -= wt
        out = chem.match(a, b)
        pop[a] -= 1
        pop[b] -= 1
        for p in out:
            pop[p] += 1
        tally.add((a, b), tuple(out))
        info["bimolecular_events"] += 1
        ok = settle()
        yield info["bimolecular_events"], pop
        if not ok:
            break
    info["inert"] = info["stopped"] is None and not propensities()[0]
    info["elastic"] = sorted(species_id(m) for m in elastic)


# --- faces ---------------------------------------------------------------------
SCHEDULING = {
    "generate": "every reaction reachable from the program fraglets; no kinetics",
    "evolve": (
        "PyCellChemistry Fraglets.propensity/react and Cell.gillespie: before each step all "
        "transformations fire at once (instantaneous); then one match is chosen with "
        "propensity n_active * n_passive (mass action with k = 1) over all nodes; "
        "steps counts the matches"
    ),
}


def _state(pop) -> dict[str, float]:
    return {species_id(m): float(c) for m, c in pop.items() if c > 0}


def _network(mols, rows, status, initial, extras) -> Network:
    species = [Species(species_id(m), structure(m)) for m in mols]
    reactions = [
        Reaction(dict(Counter(species_id(m) for m in lhs)), dict(Counter(species_id(m) for m in rhs)),
                 None, count)
        for lhs, rhs, count in rows
    ]
    compartments: dict[str, list[str]] = {}
    for m in mols:
        compartments.setdefault(m[0], []).append(species_id(m))
    extras["compartments"] = compartments
    return Network(
        species=species,
        reactions=reactions,
        status=status,
        initial_state={species_id(m): float(c) for m, c in initial.items()},
        extras=extras,
    )


def generate(p, rng):
    """The closure of the distinct program fraglets, cut off by max_species."""
    initial, segments, nodes = parse_program(p.program, p.dialect)
    chem = Chemistry(p.dialect, segments, nodes)
    mols, rxns, status = closure(chem, list(initial), p.max_species)
    extras = {"segments": segments, "scheduling": SCHEDULING["generate"]}
    return _network(mols, [(lhs, rhs, None) for lhs, rhs in rxns], status, initial, extras)


def evolve(p, rng):
    """A PyCellChemistry-scheduled run of the program multiset: a frame after every match.

    The first frame is the program as written; if its transformations fire
    before the first match, a second frame at t = 0 holds the settled state.
    """
    initial, segments, nodes = parse_program(p.program, p.dialect)
    chem = Chemistry(p.dialect, segments, nodes)
    tally, info = Tally(), {}
    yield Frame(t=0.0, state=_state(initial))
    final = initial
    for matches, pop in ssa(chem, initial, p.steps, rng, tally, info):
        final = pop
        fired = [[[species_id(m) for m in lhs], [species_id(m) for m in rhs], n]
                 for lhs, rhs, n in tally.flush()]
        if matches or fired:
            yield Frame(t=float(matches), state=_state(pop), fired=fired)
    final = Counter({m: c for m, c in final.items() if c > 0})
    rows = tally.reactions()
    seen = dict.fromkeys(initial)
    for lhs, rhs, _ in rows:
        seen.update(dict.fromkeys(lhs))
        seen.update(dict.fromkeys(rhs))
    seen.update(dict.fromkeys(final))
    extras = {"segments": segments, "scheduling": SCHEDULING["evolve"],
              "final_state": {species_id(m): int(c) for m, c in final.items()}}
    extras.update(info)
    return _network(list(seen), rows, "observed", initial, extras)
