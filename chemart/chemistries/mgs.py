"""MGS: transformations of topological collections (Giavitto & Michel).

Catalog id: mgs. Book 9.6 (and 18.3.2); Giavitto & Michel, UMC 2002 (LNCS 2509)
and Cohen, ENTCS 86(2) 2003 for the pattern language, the collection types and
the published programs.

A topological collection is a set of values organised by a neighbourhood
relation. A transformation is a list of rules ``pattern => replacement``; a rule
selects a path of neighbouring elements matching its pattern and substitutes the
replacement for it. Rules are prioritised by their order.

Collection types (``collection``):

- ``set``, ``bag``: every element is a neighbour of every other (Leibnizian,
  no undefined positions). A rule is a multiset rewriting; species are values.
- ``seq``: ``x, y`` means y is the right neighbour of x.
- ``grid``: the GBF ``<north, east>`` (von Neumann neighbourhood); ``hexagon``:
  the GBF ``<east, north, northeast; east + north = northeast>``. Both are finite
  windows (optionally a torus) of the Newtonian space: positions exist a priori
  and may hold ``<undef>``. ``x, y`` means any neighbour; ``x |d> y`` means y is
  the d-neighbour of x, where d is a sum of generators (``north-east``).
- ``tube``: a bag of sequences (the restriction-enzyme program): a rule matches a
  whole strand and its replacement is ``strand :: strand :: ():TUBE``.

On positional collections (seq, grid, hexagon) species are position-qualified,
``value@position`` (``3@2``, ``true@1_4``, ``undef@0_0``), a replacement must
have as many elements as the pattern, and every reaction is a position-specific
rewriting that conserves one species per site.

The rule language (never evaluated as Python):

    rule     ::= [Label =] pattern => replacement
    pattern  ::= element (("," | "|dir>") element)*
    element  ::= var | literal | <undef> | X+ | X* | x+ as X | x* as X,
                 each optionally followed by "/ guard"
    literal  ::= integer | "string" | true | false
    guard    ::= expression over the pattern variables with
                 or and not = == != < <= > >= + - * mod div size( ) and parentheses
    replacement ::= expr, expr, ...   |  [expr, ...]  |  expr :: expr :: empty_seq  |  ()

Guards are conjoined and evaluated once every variable is bound. Sequence
variables (X+ one or more, X* zero or more elements) are only for tube strands.

The network is the closure of the initial collection under the transformation
(every local rewriting reachable species by species). ``extras["run"]`` holds
one run of the transformation iterated to its fixpoint (no rule changes the
collection) with the maximal-parallel strategy of MGS (each rule in priority
order selects a maximal set of non-intersecting matches) or the asynchronous one
(one match per step).
"""

from __future__ import annotations

import re
from collections import Counter
from itertools import permutations, product

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

COLLECTIONS = ("set", "bag", "seq", "grid", "hexagon", "tube")
POSITIONAL = ("seq", "grid", "hexagon")
STRATEGIES = ("maximal-parallel", "asynchronous")
MAX_STEPS = 10000

GENERATORS = {
    "grid": {"north": (-1, 0), "east": (0, 1), "south": (1, 0), "west": (0, -1)},
    "hexagon": {"east": (0, 1), "north": (-1, 0), "northeast": (-1, 1),
                "west": (0, -1), "south": (1, 0), "southwest": (1, -1)},
}
GBF = {"grid": "<north, east>", "hexagon": "<east, north, northeast; east + north = northeast>"}
_UNITS = {"grid": [(-1, 0), (0, 1), (1, 0), (0, -1)],
          "hexagon": [(-1, 0), (-1, 1), (0, 1), (1, 0), (1, -1), (0, -1)]}


def bead_rows(numbers: list[int], width: int) -> list[list[bool]]:
    """Cohen 2003 fig. 4c: each number is a row of beads (true) justified to the west."""
    return [[k < n for k in range(width)] for n in numbers]


def _eden_start(size: int) -> list[list]:
    rows = [[None] * size for _ in range(size)]
    rows[size // 2][size // 2] = True
    return rows


PROGRAMS = {
    # book 9.6 [317]
    "sort": {"collection": "seq", "rules": ["x, y / x > y => y, x"], "initial": [4, 2, 5, 1, 3]},
    # Cohen 2003 sec. 2.1
    "sieve": {"collection": "set", "rules": ["x, y/(y mod x = 0) => [x]", "x => [x]"],
              "initial": list(range(2, 31))},
    # book 9.6 (EcoRI) and Giavitto & Michel 2002 sec. 3.1 (Void, the input tube)
    "restriction-enzymes": {
        "collection": "tube",
        "rules": ['EcoRI = X+, "G","A","A","T","T","C", Y+ => (X,"G") :: ("A","A","T","T","C",Y) :: ():TUBE',
                  "Void = x+ as X => X :: ():TUBE"],
        "initial": ["CCCGAATTCAA", "TTGAATTCGGG"]},
    # Giavitto & Michel 2002 sec. 3.2
    "eden": {"collection": "grid", "rules": ["x, <undef> / x => x, true"], "initial": _eden_start(7)},
    # Cohen 2003 sec. 5.2, fig. 4
    "bead-sort": {"collection": "grid",
                  "rules": ["x/x=false |north> y/y=true => y::x::empty_seq", "x=>[x]"],
                  "initial": bead_rows([3, 2, 4, 2], 4)},
    # Giavitto & Michel 2002 sec. 3.2, fig. 2 (X = east, Y = north)
    "turn": {"collection": "grid",
             "rules": ["a |east> b |north-east> c |-east-north> d |east-north> e => a, e, b, c, d"],
             "initial": [[None, 2, None], [3, 0, 1], [None, 4, None]]},
}

# --- values ---------------------------------------------------------------------------
_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_RESERVED = {"true", "false", "undef", "and", "or", "not", "mod", "div", "size", "as"}


def _is_int(v) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def value_id(v) -> str:
    """The text of a value in species ids: 3, -1, true, false, undef, G."""
    if v is None:
        return "undef"
    if isinstance(v, bool):
        return "true" if v else "false"
    if _is_int(v):
        return str(v)
    if isinstance(v, str) and _NAME.match(v) and v not in _RESERVED:
        return v
    raise ValueError(f"collection values must be integers, booleans, null (<undef>) or names "
                     f"[A-Za-z_][A-Za-z0-9_]* (not a keyword), got {v!r}")


def _eq(a, b) -> bool:
    if isinstance(a, tuple) or isinstance(b, tuple):
        return (isinstance(a, tuple) and isinstance(b, tuple) and len(a) == len(b)
                and all(_eq(x, y) for x, y in zip(a, b)))
    return type(a) is type(b) and a == b


# --- expressions -------------------------------------------------------------------------
_TOKEN = re.compile(r"""\s*(?:(<undef>)|(\d+)|("[^"]*"|'[^']*')|(==|!=|<=|>=|&&|\|\||[-+*()=<>!])|([A-Za-z_]\w*))""")
_PREC = {"or": 1, "and": 2, "=": 4, "!=": 4, "<": 4, "<=": 4, ">": 4, ">=": 4,
         "+": 5, "-": 5, "*": 6, "mod": 6, "div": 6}
_ALIAS = {"==": "=", "&&": "and", "||": "or", "!": "not"}


def _tokens(text: str, where: str) -> list[tuple[str, object]]:
    out, i, text = [], 0, text.rstrip()
    while i < len(text):
        m = _TOKEN.match(text, i)
        if not m or m.end() == i:
            raise ValueError(f"{where}: cannot read {text[i:].strip()!r}")
        undef, num, string, op, name = m.groups()
        if undef:
            out.append(("lit", None))
        elif num:
            out.append(("lit", int(num)))
        elif string:
            out.append(("lit", string[1:-1]))
        elif op:
            out.append(("op", _ALIAS.get(op, op)))
        elif name in ("true", "false"):
            out.append(("lit", name == "true"))
        elif name in ("and", "or", "not", "mod", "div"):
            out.append(("op", name))
        else:
            out.append(("name", name))
        i = m.end()
    return out


class _Parser:
    def __init__(self, text: str, where: str):
        self.where, self.toks, self.i = where, _tokens(text, where), 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else (None, None)

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def fail(self, msg):
        raise ValueError(f"{self.where}: {msg}")

    def parse(self):
        if not self.toks:
            self.fail("empty expression")
        node = self.expr(0)
        if self.i != len(self.toks):
            self.fail(f"unexpected {self.peek()[1]!r}")
        return node

    def expr(self, rbp):
        left = self.nud()
        while True:
            kind, op = self.peek()
            if kind != "op" or op not in _PREC or _PREC[op] <= rbp:
                return left
            self.take()
            left = (op, left, self.expr(_PREC[op]))

    def nud(self):
        kind, v = self.take()
        if kind == "lit":
            return ("lit", v)
        if kind == "name":
            if v == "size":
                if self.take() != ("op", "("):
                    self.fail("size needs parentheses: size(X)")
                node = ("size", self.expr(0))
                if self.take() != ("op", ")"):
                    self.fail("missing ')' after size(")
                return node
            return ("var", v)
        if (kind, v) == ("op", "("):
            node = self.expr(0)
            if self.take() != ("op", ")"):
                self.fail("missing ')'")
            return node
        if (kind, v) == ("op", "-"):
            return ("neg", self.expr(7))
        if (kind, v) == ("op", "not"):
            return ("not", self.expr(3))
        self.fail(f"unexpected {v!r}" if kind else "expression ends too early")


def _names(node) -> set[str]:
    if node[0] == "var":
        return {node[1]}
    if node[0] == "lit":
        return set()
    return set().union(*(_names(n) for n in node[1:]))


def evaluate(node, env: dict, where: str = "expression"):
    """Evaluate a parsed guard or replacement expression (no Python eval)."""
    k = node[0]
    if k == "lit":
        return node[1]
    if k == "var":
        return env[node[1]]
    if k == "size":
        v = evaluate(node[1], env, where)
        if not isinstance(v, tuple):
            raise ValueError(f"{where}: size() applies to a sequence variable")
        return len(v)
    if k == "not":
        v = evaluate(node[1], env, where)
        if not isinstance(v, bool):
            raise ValueError(f"{where}: 'not' needs a boolean, got {v!r}")
        return not v
    if k == "neg":
        v = evaluate(node[1], env, where)
        if not _is_int(v):
            raise ValueError(f"{where}: '-' needs an integer, got {v!r}")
        return -v
    if k in ("and", "or"):
        a = evaluate(node[1], env, where)
        if not isinstance(a, bool):
            raise ValueError(f"{where}: '{k}' needs booleans, got {a!r}")
        if (k == "and" and not a) or (k == "or" and a):
            return a
        b = evaluate(node[2], env, where)
        if not isinstance(b, bool):
            raise ValueError(f"{where}: '{k}' needs booleans, got {b!r}")
        return b
    a, b = evaluate(node[1], env, where), evaluate(node[2], env, where)
    if k == "=":
        return _eq(a, b)
    if k == "!=":
        return not _eq(a, b)
    if k in ("<", "<=", ">", ">="):
        if not ((_is_int(a) and _is_int(b)) or (isinstance(a, str) and isinstance(b, str))):
            raise ValueError(f"{where}: cannot compare {a!r} {k} {b!r}")
        return {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b}[k]
    if not (_is_int(a) and _is_int(b)):
        raise ValueError(f"{where}: '{k}' needs integers, got {a!r} and {b!r}")
    if k in ("mod", "div") and b == 0:
        raise ValueError(f"{where}: {k} by zero")
    return {"+": a + b, "-": a - b, "*": a * b, "mod": a % b if k == "mod" else 0,
            "div": a // b if k == "div" else 0}[k]


# --- rules -------------------------------------------------------------------------------
_HOP = re.compile(r"\|\s*([-+]?\s*[A-Za-z_]\w*(?:\s*[-+]\s*[A-Za-z_]\w*)*)\s*>")
_LABEL = re.compile(r"^\s*([A-Za-z_]\w*)\s*=(?![=>])\s*(.*)$", re.S)
_SEQVAR = re.compile(r"^([A-Za-z_]\w*)\s*([+*])(?:\s+as\s+([A-Za-z_]\w*))?$")
_EMPTY = re.compile(r"^(\(\s*\)(\s*:\s*\w+)?|\[\s*\]|empty[_ ]seq)$")


def _scan(text: str, seps: tuple[str, ...], hops: bool = False):
    """Split at top-level separators (outside quotes and brackets); returns (parts, separators)."""
    parts, found, buf, depth, quote, i = [], [], [], 0, None, 0
    while i < len(text):
        ch = text[i]
        if quote:
            quote = None if ch == quote else quote
        elif ch in "\"'":
            quote = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif depth == 0:
            if hops and ch == "|":
                m = _HOP.match(text, i)
                if m:
                    parts.append("".join(buf))
                    found.append(re.sub(r"\s+", "", m.group(1)))
                    buf, i = [], m.end()
                    continue
            sep = next((s for s in seps if text.startswith(s, i)), None)
            if sep:
                parts.append("".join(buf))
                found.append(None)
                buf, i = [], i + len(sep)
                continue
        buf.append(ch)
        i += 1
    parts.append("".join(buf))
    return [p.strip() for p in parts], found


def _unwrap(s: str, open_: str, close: str) -> str:
    s = s.strip()
    if not (s.startswith(open_) and s.endswith(close)):
        return s
    depth = 0
    for i, ch in enumerate(s):
        depth += ch == open_
        depth -= ch == close
        if depth == 0 and i < len(s) - 1:
            return s
    return s[1:-1].strip()


class Element:
    def __init__(self, kind: str, name: str | None = None, value=None, minimum: int = 1):
        self.kind, self.name, self.value, self.minimum = kind, name, value, minimum


class Rule:
    """A parsed MGS rule for one collection type."""

    def __init__(self, text: str, collection: str, label: str = ""):
        self.text = " ".join(text.split())
        m = _LABEL.match(text)
        body = text
        if m:
            label, body = m.group(1), m.group(2)
        self.label = label
        where = f"rule {self.text!r}"
        if body.count("=>") != 1:
            raise ValueError(f"{where} needs exactly one '=>'")
        lhs, rhs = body.split("=>")
        parts, hops = _scan(lhs, (",",), hops=True)
        self.elements, self.guards, bound = [], [], {}
        for part in parts:
            head, guards = _scan(part, ("/",))
            head = head.strip() if isinstance(head, str) else head
            el = self._element(head[0], where)
            if el.name:
                if el.name in bound:
                    raise ValueError(f"{where}: variable {el.name} appears twice (patterns are linear)")
                bound[el.name] = el.kind
            self.elements.append(el)
            self.guards += [_Parser(g, where).parse() for g in head[1:]]
        self.hops = hops
        self.bound = bound
        for g in self.guards:
            self._check_names(g, where)
        self._check_collection(collection, where)
        self.collection = collection
        self.rhs = self._replacement(rhs, collection, where)

    @staticmethod
    def _element(head: str, where: str) -> Element:
        if not head:
            raise ValueError(f"{where}: empty pattern element")
        if head == "<undef>":
            return Element("undef")
        m = _SEQVAR.match(head)
        if m:
            return Element("seqvar", m.group(3) or m.group(1), minimum=1 if m.group(2) == "+" else 0)
        toks = _tokens(head, where)
        if len(toks) == 2 and toks[0] == ("op", "-") and toks[1][0] == "lit" and _is_int(toks[1][1]):
            return Element("lit", value=-toks[1][1])
        if len(toks) == 1 and toks[0][0] == "lit" and toks[0][1] is not None:
            value_id(toks[0][1])
            return Element("lit", value=toks[0][1])
        if len(toks) == 1 and toks[0][0] == "name" and toks[0][1] not in _RESERVED:
            return Element("var", name=toks[0][1])
        raise ValueError(f"{where}: cannot read pattern element {head!r}")

    def _check_names(self, node, where):
        unbound = _names(node) - set(self.bound)
        if unbound:
            raise ValueError(f"{where}: unbound variable(s) {sorted(unbound)}")

    def _check_collection(self, collection, where):
        kinds = {e.kind for e in self.elements}
        if collection != "tube" and "seqvar" in kinds:
            raise ValueError(f"{where}: sequence variables (X+, X*) are only supported in a tube")
        if collection in ("set", "bag", "tube") and "undef" in kinds:
            raise ValueError(f"{where}: <undef> only exists in seq, grid and hexagon collections")
        if any(self.hops) and collection not in GENERATORS:
            raise ValueError(f"{where}: directions |d> need a grid or hexagon collection")
        self.vectors = []
        for h in self.hops:
            if h is None:
                self.vectors.append(None)
                continue
            vec, gens = (0, 0), GENERATORS[collection]
            for sign, name in re.findall(r"([-+]?)([A-Za-z_]\w*)", h):
                if name not in gens:
                    raise ValueError(f"{where}: unknown direction {name!r} for {collection} "
                                     f"(use {', '.join(gens)})")
                s = -1 if sign == "-" else 1
                vec = (vec[0] + s * gens[name][0], vec[1] + s * gens[name][1])
            self.vectors.append(vec)

    def _replacement(self, rhs, collection, where):
        groups, _ = _scan(rhs.strip(), ("::",))
        groups = [g for g in groups if g and not _EMPTY.match(g)]
        if collection == "tube":
            strands = []
            for g in groups:
                items, _ = _scan(_unwrap(_unwrap(g, "(", ")"), "[", "]"), (",",))
                strand = []
                for item in items:
                    if not item:
                        raise ValueError(f"{where}: empty item in replacement strand {g!r}")
                    if item in self.bound and self.bound[item] == "seqvar":
                        strand.append(("splice", item))
                    else:
                        node = _Parser(item, where).parse()
                        self._check_names(node, where)
                        strand.append(("expr", node))
                strands.append(strand)
            return strands
        exprs = []
        for g in groups:
            items, _ = _scan(_unwrap(_unwrap(g, "[", "]"), "(", ")"), (",",))
            for item in items:
                if not item:
                    if len(items) == 1:
                        continue
                    raise ValueError(f"{where}: empty item in the replacement")
                node = _Parser(item, where).parse()
                self._check_names(node, where)
                exprs.append(node)
        if collection in POSITIONAL and len(exprs) != len(self.elements):
            raise ValueError(f"{where}: on a {collection} the replacement must have as many elements "
                             f"as the pattern ({len(self.elements)}), got {len(exprs)}")
        return exprs

    @property
    def arity(self) -> int:
        return len(self.elements)

    def match(self, values) -> dict | None:
        """Bindings if the pattern (elements and guards) matches these values, else None."""
        env = {}
        for el, v in zip(self.elements, values):
            if el.kind == "undef":
                if v is not None:
                    return None
            elif v is None:
                return None
            elif el.kind == "lit":
                if not _eq(v, el.value):
                    return None
            else:
                env[el.name] = v
        return env if self._guards(env) else None

    def _guards(self, env) -> bool:
        for g in self.guards:
            v = evaluate(g, env, f"rule {self.text!r}")
            if not isinstance(v, bool):
                raise ValueError(f"rule {self.text!r}: a guard must be boolean, got {v!r}")
            if not v:
                return False
        return True

    def match_strand(self, atoms: tuple) -> list[dict]:
        """Every binding under which the pattern matches the whole strand."""
        found, els = [], self.elements

        def rec(k, pos, env):
            if k == len(els):
                if pos == len(atoms) and self._guards(env):
                    found.append(dict(env))
                return
            el = els[k]
            if el.kind == "seqvar":
                for end in range(pos + el.minimum, len(atoms) + 1):
                    env[el.name] = tuple(atoms[pos:end])
                    rec(k + 1, end, env)
                env.pop(el.name, None)
            elif pos < len(atoms):
                if el.kind == "lit":
                    if _eq(atoms[pos], el.value):
                        rec(k + 1, pos + 1, env)
                else:
                    env[el.name] = atoms[pos]
                    rec(k + 1, pos + 1, env)
                    env.pop(el.name)

        rec(0, 0, {})
        return found

    def values(self, env) -> list:
        return [evaluate(e, env, f"rule {self.text!r}") for e in self.rhs]

    def strands(self, env) -> list[tuple]:
        out = []
        for strand in self.rhs:
            atoms = []
            for kind, x in strand:
                if kind == "splice":
                    atoms.extend(env[x])
                else:
                    v = evaluate(x, env, f"rule {self.text!r}")
                    atoms.extend(v if isinstance(v, tuple) else (v,))
            out.append(tuple(atoms))
        return out


# --- topology ----------------------------------------------------------------------------
class Topology:
    """Positions and neighbourhood of a seq, grid or hexagon collection."""

    def __init__(self, kind: str, shape: tuple[int, ...], torus: bool = False):
        self.kind, self.shape, self.torus = kind, shape, torus
        if kind == "seq":
            self.positions = [(i,) for i in range(shape[0])]
        else:
            self.positions = [(r, c) for r in range(shape[0]) for c in range(shape[1])]
        self.index = {p: i for i, p in enumerate(self.positions)}

    def label(self, p) -> str:
        return "_".join(str(x) for x in p)

    def shift(self, p, vec):
        if self.kind == "seq":
            q = p[0] + vec[0]
            return (q,) if 0 <= q < self.shape[0] else None
        r, c = p[0] + vec[0], p[1] + vec[1]
        if self.torus:
            return (r % self.shape[0], c % self.shape[1])
        return (r, c) if 0 <= r < self.shape[0] and 0 <= c < self.shape[1] else None

    def neighbours(self, p) -> list:
        out = []
        for vec in ([(1,)] if self.kind == "seq" else _UNITS[self.kind]):
            q = self.shift(p, vec)
            if q is not None and q != p and q not in out:
                out.append(q)
        return out

    def paths(self, vectors: list) -> list[tuple]:
        """Tuples of distinct positions linked by the pattern's hops (None: any neighbour)."""
        out = []

        def rec(path):
            if len(path) == len(vectors) + 1:
                out.append(tuple(path))
                return
            vec = vectors[len(path) - 1]
            for q in (self.neighbours(path[-1]) if vec is None else [self.shift(path[-1], vec)]):
                if q is not None and q not in path:
                    rec(path + [q])

        for p in self.positions:
            rec([p])
        return out

    def edges(self) -> list[list[str]]:
        seen = []
        for p in self.positions:
            for q in self.neighbours(p):
                pair = sorted([self.index[p], self.index[q]])
                if pair not in seen:
                    seen.append(pair)
        return [[self.label(self.positions[a]), self.label(self.positions[b])] for a, b in seen]

    def space(self) -> dict:
        out = {"collection": self.kind, "species": "value@position"}
        if self.kind == "seq":
            out.update(positions=self.shape[0], neighbourhood="x, y: y is the right neighbour of x")
        else:
            out.update(gbf=GBF[self.kind], shape=list(self.shape), torus=self.torus,
                       position="row_column, row 0 to the north",
                       directions={k: list(v) for k, v in GENERATORS[self.kind].items()},
                       neighbourhood="x, y: y is any neighbour of x; x |d> y: y is the d-neighbour of x")
        out["edges"] = self.edges()
        return out


# --- the machine: closure and runs -------------------------------------------------------
class Machine:
    """A transformation (list of rules) acting on one collection."""

    def __init__(self, collection: str, rules: list[str], initial: list, torus: bool = False):
        if collection not in COLLECTIONS:
            raise ValueError(f"collection must be one of {COLLECTIONS}, got {collection!r}")
        if not isinstance(rules, list) or not rules or not all(isinstance(r, str) for r in rules):
            raise ValueError("rules must be a non-empty list of rule strings 'pattern => replacement'")
        self.collection = collection
        self.rules = [Rule(t, collection, f"r{i + 1}") for i, t in enumerate(rules)]
        self.values: dict[str, object] = {}
        self.strands: dict[str, tuple] = {}
        if torus and collection not in GENERATORS:
            raise ValueError("torus only applies to grid and hexagon collections")
        if not isinstance(initial, list) or not initial:
            raise ValueError(f"initial must be a non-empty list, got {initial!r}")
        if collection in POSITIONAL:
            if collection == "seq":
                if any(isinstance(v, list) or v is None for v in initial):
                    raise ValueError("a seq initial collection is a flat list of values without null")
                self.topology, flat = Topology("seq", (len(initial),)), initial
            else:
                if not all(isinstance(r, list) and r for r in initial) or len({len(r) for r in initial}) != 1:
                    raise ValueError(f"a {collection} initial collection is a list of equal-length rows")
                self.topology = Topology(collection, (len(initial), len(initial[0])), torus)
                flat = [v for row in initial for v in row]
            self.start = tuple(self._vid(v) for v in flat)
            self.paths = [self.topology.paths(r.vectors) for r in self.rules]
            self.pathsets = [set(ps) for ps in self.paths]
        elif collection == "tube":
            strands = []
            for s in initial:
                if isinstance(s, str):
                    s = list(s)
                if not isinstance(s, list):
                    raise ValueError(f"a tube holds strands (strings or lists of values), got {s!r}")
                strands.append(self._sid(tuple(s)))
            self.start = tuple(strands)
        else:
            if any(isinstance(v, list) or v is None for v in initial):
                raise ValueError(f"a {collection} initial collection is a flat list of values without null")
            ids = [self._vid(v) for v in initial]
            self.start = tuple(dict.fromkeys(ids)) if collection == "set" else tuple(ids)

    # ids ---------------------------------------------------------------------------------
    def _vid(self, v) -> str:
        i = value_id(v)
        self.values.setdefault(i, v)
        return i

    def _sid(self, atoms: tuple) -> str:
        ids = [self._vid(a) for a in atoms]
        sid = "".join(ids) if ids and all(len(i) == 1 for i in ids) else "(" + ",".join(ids) + ")"
        self.strands.setdefault(sid, tuple(atoms))
        return sid

    def species_id(self, pos, vid: str) -> str:
        return f"{vid}@{self.topology.label(pos)}"

    def _shadowed(self, j: int, vals, path=None) -> bool:
        """A higher-priority rule matches the same sub-collection."""
        return any(r.arity == len(vals) and (path is None or path in self.pathsets[i])
                   and r.match(vals) is not None for i, r in enumerate(self.rules[:j]))

    def _outcome(self, j: int, ids: tuple):
        """Replacement ids of rule j on these element ids (None if it does not apply)."""
        r = self.rules[j]
        if self.collection == "tube":
            for i in range(j):
                if self.rules[i].match_strand(self.strands[ids[0]]):
                    return None
            outs = []
            for env in r.match_strand(self.strands[ids[0]]):
                out = tuple(self._sid(s) for s in r.strands(env))
                if out not in outs:
                    outs.append(out)
            return outs or None
        vals = [self.values[i] for i in ids]
        if self.collection == "set" and len(set(ids)) < len(ids):
            return None
        env = r.match(vals)
        if env is None or self._shadowed(j, vals):
            return None
        out = [self._vid(v) for v in r.values(env)]
        if self.collection in ("set", "bag") and any(v is None for v in r.values(env)):
            raise ValueError(f"rule {r.text!r}: <undef> cannot be put in a {self.collection}")
        return [tuple(dict.fromkeys(out)) if self.collection == "set" else tuple(out)]

    # closure -----------------------------------------------------------------------------
    def closure(self, max_species: int = 5000):
        """(species ids, reactions [(lhs ids, rhs ids, rule label)], status)."""
        if self.collection in POSITIONAL:
            return self._positional_closure(max_species)
        labels: dict = {}

        def react(*ids):
            outs = []
            for j, r in enumerate(self.rules):
                if r.arity != len(ids) and self.collection != "tube":
                    continue
                for out in self._outcome(j, ids) or []:
                    if Counter(out) != Counter(ids):
                        labels.setdefault(_key(ids, out), r.label)
                        outs.append(out)
            return outs or None

        arity = 1 if self.collection == "tube" else sorted({r.arity for r in self.rules})
        species, found, status = expand(react, self.start, arity=arity, max_species=max_species,
                                        ordered=True, alternatives=True)
        return species, [(lhs, rhs, labels[_key(lhs, rhs)]) for lhs, rhs in found], status

    def _positional_closure(self, max_species):
        topo = self.topology
        species, index, site = [], {}, {p: [] for p in topo.positions}

        def add(pos, vid):
            index[(pos, vid)] = len(species)
            site[pos].append(len(species))
            species.append((pos, vid))

        for pos, vid in zip(topo.positions, self.start):
            add(pos, vid)
        reactions, truncated, done = {}, False, 0
        while done < len(species):
            old, n = done, len(species)
            done = n
            for j, r in enumerate(self.rules):
                for path in self.paths[j]:
                    lists = [[i for i in site[p] if i < n] for p in path]
                    for first in range(len(path)):
                        choice = ([[i for i in lists[q] if i < old] for q in range(first)]
                                  + [[i for i in lists[first] if i >= old]] + lists[first + 1:])
                        for combo in product(*choice):
                            vids = [species[i][1] for i in combo]
                            vals = [self.values[v] for v in vids]
                            env = r.match(vals)
                            if env is None or self._shadowed(j, vals, path):
                                continue
                            out = [self._vid(v) for v in r.values(env)]
                            if out == vids:
                                continue
                            novel = [(p, o) for p, o in zip(path, out) if (p, o) not in index]
                            if len(species) + len(novel) > max_species:
                                truncated = True
                                continue
                            for p, o in novel:
                                add(p, o)
                            lhs = tuple(self.species_id(p, v) for p, v in zip(path, vids))
                            rhs = tuple(self.species_id(p, v) for p, v in zip(path, out))
                            reactions.setdefault(_key(lhs, rhs), (lhs, rhs, r.label))
        ids = [self.species_id(p, v) for p, v in species]
        return ids, list(reactions.values()), "truncated" if truncated else "complete"

    # runs --------------------------------------------------------------------------------
    def matches(self, state: tuple) -> list[list[tuple]]:
        """Per rule, the non-elastic matches (elements, replacement) in this state."""
        per_rule = [[] for _ in self.rules]
        if self.collection in POSITIONAL:
            idx = self.topology.index
            for j, r in enumerate(self.rules):
                for path in self.paths[j]:
                    vids = [state[idx[p]] for p in path]
                    vals = [self.values[v] for v in vids]
                    env = r.match(vals)
                    if env is None or self._shadowed(j, vals, path):
                        continue
                    out = [self._vid(v) for v in r.values(env)]
                    if out != vids:
                        per_rule[j].append((path, tuple(out)))
            return per_rule
        for j, r in enumerate(self.rules):
            k = 1 if self.collection == "tube" else r.arity
            for elems in permutations(range(len(state)), k):
                ids = tuple(state[e] for e in elems)
                for out in self._outcome(j, ids) or []:
                    if Counter(out) != Counter(ids):
                        per_rule[j].append((elems, out))
        return per_rule

    def step(self, state: tuple, strategy: str, rng):
        """One transformation step; None when no rule changes the collection."""
        per_rule = self.matches(state)
        chosen, used = [], set()
        if strategy == "asynchronous":
            ms = next((m for m in per_rule if m), None)
            if ms:
                chosen.append(ms[int(rng.integers(len(ms)))])
        else:
            for ms in per_rule:
                for k in (rng.permutation(len(ms)) if ms else []):
                    elems, out = ms[int(k)]
                    if used.isdisjoint(elems):
                        used.update(elems)
                        chosen.append((elems, out))
        if not chosen:
            return None
        if self.collection in POSITIONAL:
            new = list(state)
            for path, out in chosen:
                for p, o in zip(path, out):
                    new[self.topology.index[p]] = o
            return tuple(new)
        gone = {e for elems, _ in chosen for e in elems}
        new = [x for i, x in enumerate(state) if i not in gone]
        for _, out in chosen:
            new.extend(out)
        return tuple(dict.fromkeys(new)) if self.collection == "set" else tuple(new)

    def run(self, strategy: str, rng, max_steps: int = MAX_STEPS):
        """Iterate to the fixpoint: (final state, steps, reached fixpoint)."""
        state = self.start
        for steps in range(max_steps):
            new = self.step(state, strategy, rng)
            if new is None:
                return state, steps, True
            state = new
        return state, max_steps, False

    def to_json(self, state: tuple):
        if self.collection == "seq":
            return [self.values[v] for v in state]
        if self.collection in GENERATORS:
            cols = self.topology.shape[1]
            return [[self.values[v] for v in state[r * cols:(r + 1) * cols]]
                    for r in range(self.topology.shape[0])]
        if self.collection == "set":
            return [self.values[v] for v in state]
        return dict(Counter(state))


def _key(lhs, rhs):
    return frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items())


# --- generator ----------------------------------------------------------------------------
def program(p) -> tuple[str, list[str], list]:
    """(collection, rules, initial collection) selected by the parameters."""
    if p.program == "custom":
        if p.collection == "default":
            raise ValueError("program custom needs a collection (set, bag, seq, grid, hexagon or tube)")
        if not p.rules or not p.initial:
            raise ValueError("program custom needs rules (a list of 'pattern => replacement') and initial")
        return p.collection, list(p.rules), list(p.initial)
    if p.rules:
        raise ValueError("rules are only used with program='custom'")
    spec = PROGRAMS[p.program]
    collection = spec["collection"] if p.collection == "default" else p.collection
    return collection, list(spec["rules"]), list(p.initial) if p.initial else spec["initial"]


def generate(p, rng):
    collection, rules, initial = program(p)
    m = Machine(collection, rules, initial, p.torus)
    ids, found, status = m.closure(p.max_species)
    final, steps, fixpoint = m.run(p.strategy, rng)

    extras = {
        "program": p.program,
        "rules": [f"{r.label}: {r.text}" for r in m.rules],
        "reaction_rules": [label for _, _, label in found],
        "run": {"strategy": p.strategy, "steps": steps, "fixpoint": fixpoint, "final": m.to_json(final)},
    }
    if collection in POSITIONAL:
        topo = m.topology
        species = [Species(i) for i in ids]
        initial_state = {m.species_id(pos, v): 1.0 for pos, v in zip(topo.positions, m.start)}
        extras["space"] = topo.space()
        by_site: dict = {}
        for i in ids:
            by_site.setdefault(i.rsplit("@", 1)[1], {})[i] = 1
        extras["conservation"] = [{"name": f"site {s}", "vector": v} for s, v in by_site.items()]
    else:
        if collection == "tube":
            species = [Species(i, structure=",".join(value_id(a) for a in m.strands[i])) for i in ids]
            hood = "a bag of sequences: every strand neighbours every other; a rule matches a whole strand"
        else:
            species = [Species(i) for i in ids]
            hood = "complete: every element is a neighbour of every other"
        initial_state = {i: float(c) for i, c in Counter(m.start).items()}
        extras["space"] = {"collection": collection, "neighbourhood": hood}
    reactions = [Reaction.of(lhs, rhs) for lhs, rhs, _ in found]
    return Network(species=species, reactions=reactions, status=status,
                   initial_state=initial_state, extras=extras)
