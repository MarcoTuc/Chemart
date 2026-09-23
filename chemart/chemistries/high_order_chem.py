"""High-order chemistry: reaction rules are molecules (book appendix; PyCellChemistry HighOrderChem.py).

Catalog id: high-order-chem.

Two multisets float in the reactor: data molecules and rule molecules. One
iteration of the algorithm (book appendix, figure 3) is:

1. draw one rule molecule at random from the rule multiset (stop if empty);
2. count its binding sites, nsites (the comma-separated parameters of
   ``function(m1, m2, ...)``);
3. if the data multiset holds at least nsites molecules, draw nsites of them
   at random, without replacement, execute the rule on them and inject the
   returned products into the data multiset;
4. reinject the rule molecule.

The rule is therefore a catalyst of every reaction it executes:
``rule:divrule + n12 + n3 -> rule:divrule + n4 + n3``.

No string is ever executed. PyCellChemistry builds a Python command from the
rule string and runs it with ``exec``; here a rule is either

- a **named rule** from the built-in registry `RULES`, written as its name
  (``"divrule"``) or as the rule molecule of the reference code
  (``"divrule(m1, m2)"``, whose binding sites must match the rule's arity):

  ``divrule(m1, m2)``          NumberChemHO.py: the larger of two integers is
                               divided by the smaller when it is a multiple of it
  ``exchangeMachine(m)``       MolecularTSP.py E-machine: swap two cities
  ``cutMachine(m)``            MolecularTSP.py C-machine: move a segment
  ``invertMachine(m)``         MolecularTSP.py I-machine: move and invert a segment
  ``recombinationMachine(m1, m2)``  MolecularTSP.py R-machine: graft a segment of
                               the first tour into the second, keep the best two

  (the machines act on tours: data molecules that are permutations of
  0..cities-1 on the ring instance of MolecularTSP.py, with the reference
  code's fitness, including its ``%g`` rounding of stored fitness values);

- or an **expression rule** in a small language, parsed and interpreted here::

      [name:] x, y -> expr, expr, ... [if condition]

  The variables on the left are the binding sites (arity = their number). The
  right side lists the products (it may be empty: the educts are destroyed).
  Expressions: integers, variables, ``true``/``false``, ``+ - *``, ``/``
  (floor division), ``%``, comparisons ``== != < <= > >=``, ``and or not``,
  parentheses and the functions ``min(a, b, ...)``, ``max(a, b, ...)``,
  ``abs(a)``. Arithmetic and ordering need integers. If the condition is
  false, or an expression is ill-typed (a list in an arithmetic operation,
  division by zero), the rule returns its educts unchanged: an elastic
  collision, exactly as ``divrule`` does for a non-divisible pair.

Data molecules are JSON integers or lists (nested lists allowed). Species ids:
``n<k>`` for the integer k, the compact JSON text (``[0,3,1,2]``) for a list;
``structure`` holds the canonical value. Rule species are ``rule:<name>``
(``rule:expr<k>`` for the k-th unnamed expression rule) with the rule text as
structure.

evolve runs the algorithm for `iterations` iterations, a frame every
generation (as many iterations as initial data molecules), and returns the
observed network (effective reactions with firing counts). generate is a
Chemart addition: every reaction reachable from the distinct data molecules
under the (deterministic) rules, via chemart.expand.expand.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from typing import Any, Callable

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame, ticks


# ---------------------------------------------------------------------------
# Data molecules
# ---------------------------------------------------------------------------
def molecule(value: Any, where: str = "data"):
    """JSON value -> internal molecule (int, or tuple of molecules for a list)."""
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, (list, tuple)):
        return tuple(molecule(v, where) for v in value)
    raise ValueError(f"{where} molecules must be integers or lists of them, got {value!r}")


def text(m) -> str:
    """Canonical compact text: 12, [0,3,1,2]."""
    return json.dumps(to_json(m), separators=(",", ":"))


def sid(m) -> str:
    return f"n{m}" if isinstance(m, int) else text(m)


def to_json(m):
    return [to_json(e) for e in m] if isinstance(m, tuple) else m


def order(m):
    """Sort key: integers by value, then lists by length and content."""
    return (0, m) if isinstance(m, int) else (1, len(m), tuple(order(e) for e in m))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class Rule:
    id: str                  # species id, rule:<name>
    text: str                # the rule molecule, e.g. "divrule(m1, m2)"
    arity: int               # number of binding sites
    apply: Callable          # apply(ctx, rng, *molecules) -> list of products
    stochastic: bool = False
    tours: bool = False      # acts on MolecularTSP tours (needs the tour graph)


#: PyCellChemistry's rule parser (HighOrderChem.iterate, line 10 of figure 3).
_RULE_CALL = re.compile(r"([\.\w]+)\((.*)\)")


def binding_sites(rule_text: str) -> int:
    """nsites = len(bsite.split(',')), as in HighOrderChem.iterate (figure 3, line 14)."""
    match = _RULE_CALL.match(rule_text)
    if match is None:
        raise ValueError(f"malformed rule molecule {rule_text!r}: expected function(m1, ...)")
    return len(match.group(2).split(","))


def divrule(ctx, rng, m1, m2):
    """NumberChemHO.divrule: divide the larger number by the smaller if it is a multiple.

    Molecules outside the rule's domain (non-integers, zero) collide elastically.
    """
    p = [m1, m2]
    if not (isinstance(m1, int) and isinstance(m2, int)) or m1 == 0 or m2 == 0:
        return p
    if m1 > m2 and m1 % m2 == 0:
        p[0] = m1 // m2
    elif m2 > m1 and m2 % m1 == 0:
        p[1] = m2 // m1
    return p


# --- MolecularTSP.py machines ----------------------------------------------
class TourGraph:
    """MolecularTSP.TSPgraph with ring=True: n cities on a ring, fully meshed roads."""

    def __init__(self, n: int):
        self.n = n
        gs = 2 * n
        radius = gs // 2
        angle = 2 * math.pi / n
        theta = 0.0
        self.coord = []
        for _ in range(n):
            self.coord.append((math.cos(theta) * radius + radius, math.sin(theta) * radius + radius))
            theta += angle
        self.toofar = self.distance((0.0, 0.0), (float(gs), float(gs)))

    @staticmethod
    def distance(p1, p2) -> float:
        return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

    def is_tour(self, m) -> bool:
        return isinstance(m, tuple) and sorted(m) == list(range(self.n)) \
            and all(isinstance(c, int) for c in m)

    def fitness(self, tour) -> float:
        """MolecularTSP.fitness: closed tour length, with penalties (lower is better)."""
        visited = set()
        fit = 0
        for i in range(len(tour)):
            j = (i + 1) % len(tour)
            c1, c2 = tour[i], tour[j]
            if c2 != c1:  # c2 in roads[c1]: the mesh connects every pair of cities
                fit += self.distance(self.coord[c1], self.coord[c2])
            else:
                fit += 2 * self.toofar
            if c1 in visited:
                fit += self.toofar
            visited.add(c1)
        return fit

    def stored_fitness(self, tour) -> float:
        """The fitness a molecule carries: newMolecule writes it with %g, parseMolecule reads it back."""
        return float("%g" % self.fitness(tour))

    def random_tour(self, rng) -> tuple:
        """TSPgraph.randomTour(ncities)."""
        c0 = int(rng.integers(0, self.n))
        tour = [c0]
        c1 = c0
        while len(tour) < self.n:
            roads = [c for c in range(self.n) if c != c1]
            while c1 in tour and len(roads) > 0:
                c1 = roads[int(rng.integers(0, len(roads)))]
                roads.remove(c1)
            if c1 in tour:
                return tuple(tour)
            tour.append(c1)
        return tuple(tour)


def _random_pair(rng, n: int) -> tuple[int, int]:
    p1 = int(rng.integers(n))
    p2 = p1
    while p1 == p2:
        p2 = int(rng.integers(n))
    return p1, p2


def split_tour(tour, p1: int, p2: int):
    """MolecularTSP.splitTour: the circular segments [p1, p2) and [p2, p1)."""
    if p1 % len(tour) == p2 % len(tour):
        return list(tour), []
    seg1, i = [], p1
    while i != p2:
        seg1.append(tour[i])
        i = (i + 1) % len(tour)
    seg2, i = [], p2
    while i != p1:
        seg2.append(tour[i])
        i = (i + 1) % len(tour)
    return seg1, seg2


def _exchange(tour, rng):
    p1, p2 = _random_pair(rng, len(tour))
    new = list(tour)
    new[p1], new[p2] = tour[p2], tour[p1]
    return new


def _cut_or_invert(tour, rng, invert: bool):
    tour = list(tour)
    new = list(tour)
    if len(tour) < 3:
        return new
    while new == tour:
        p1, p2 = _random_pair(rng, len(tour))
        seg1, seg2 = split_tour(tour, p1, p2)
        if len(seg2) > 1:
            p3 = int(rng.integers(len(seg2) - 1)) + 1
            seg3, seg4 = split_tour(seg2, 0, p3)
        else:
            seg3, seg4 = seg2, []
        if invert:
            seg1.reverse()
        new = seg3 + seg1 + seg4
    return new


def _recombine(tour1, tour2, rng):
    p1, p2 = _random_pair(rng, len(tour1))
    seg1, _ = split_tour(tour1, p1, p2)
    city = seg1.pop(0)
    new = list(tour2)
    for c in seg1:
        new.remove(c)
    p3 = new.index(city)
    seg2, seg3 = split_tour(new, 0, p3 + 1)
    return seg2 + seg1 + seg3


def _single_machine(operator):
    def machine(graph, rng, mol):
        if graph is None or not graph.is_tour(mol):
            return [mol]
        fit = graph.stored_fitness(mol)
        new = operator(mol, rng)
        newfit = graph.fitness(new)
        return [tuple(new)] if newfit < fit else [mol]
    return machine


exchangeMachine = _single_machine(_exchange)
cutMachine = _single_machine(lambda t, rng: _cut_or_invert(t, rng, False))
invertMachine = _single_machine(lambda t, rng: _cut_or_invert(t, rng, True))


def recombinationMachine(graph, rng, mol1, mol2):
    """MolecularTSP.recombinationMachine: release the 2 fitter of parents and child."""
    if graph is None or not (graph.is_tour(mol1) and graph.is_tour(mol2)):
        return [mol1, mol2]
    fit1, fit2 = graph.stored_fitness(mol1), graph.stored_fitness(mol2)
    tour3 = tuple(_recombine(mol1, mol2, rng))
    fit3 = graph.fitness(tour3)
    if fit1 < fit3 and fit2 < fit3:
        return [mol1, mol2]
    if fit1 < fit2 and fit3 < fit2:
        return [mol1, tour3]
    return [mol2, tour3]


#: name -> (rule molecule text, function, stochastic, needs the tour graph)
RULES: dict[str, tuple[str, Callable, bool, bool]] = {
    "divrule": ("divrule(m1, m2)", divrule, False, False),
    "exchangeMachine": ("exchangeMachine(m)", exchangeMachine, True, True),
    "cutMachine": ("cutMachine(m)", cutMachine, True, True),
    "invertMachine": ("invertMachine(m)", invertMachine, True, True),
    "recombinationMachine": ("recombinationMachine(m1, m2)", recombinationMachine, True, True),
}
TOUR_RULES = {name for name, spec in RULES.items() if spec[3]}


# ---------------------------------------------------------------------------
# Expression rules: tokenizer, parser, interpreter (never Python eval)
# ---------------------------------------------------------------------------
_TOKEN = re.compile(r"\s*(?:(?P<int>\d+)|(?P<name>[A-Za-z_]\w*)|(?P<op>->|<=|>=|==|!=|[-+*/%(),:<>]))")
_KEYWORDS = {"if", "and", "or", "not", "true", "false"}
_FUNCTIONS = {"min": (2, None), "max": (2, None), "abs": (1, 1)}   # (min args, max args)
_COMPARE = ("==", "!=", "<", "<=", ">", ">=")


class _Mismatch(Exception):
    """Ill-typed evaluation: the rule returns its educts (elastic collision)."""


def _tokenize(src: str) -> list[tuple[str, Any]]:
    out, pos, src = [], 0, src.rstrip()
    while pos < len(src):
        m = _TOKEN.match(src, pos)
        if not m or m.end() == pos:
            raise ValueError(f"expression rule {src!r}: cannot read {src[pos:].strip()!r}")
        pos = m.end()
        if m.group("int") is not None:
            out.append(("int", int(m.group("int"))))
        elif m.group("name") is not None:
            name = m.group("name")
            out.append(("kw", name) if name in _KEYWORDS else ("name", name))
        else:
            out.append(("op", m.group("op")))
    return out


class _Parser:
    def __init__(self, src: str):
        self.src = src
        self.toks = _tokenize(src)
        self.i = 0

    def error(self, msg: str):
        return ValueError(f"expression rule {self.src!r}: {msg}")

    def peek(self, kind=None, value=None) -> bool:
        if self.i >= len(self.toks):
            return False
        k, v = self.toks[self.i]
        return (kind is None or k == kind) and (value is None or v == value)

    def take(self, kind, value=None):
        if not self.peek(kind, value):
            found = self.toks[self.i][1] if self.i < len(self.toks) else "end of rule"
            raise self.error(f"expected {value or kind}, found {found!r}")
        self.i += 1
        return self.toks[self.i - 1][1]

    def rule(self):
        name = None
        if len(self.toks) > 1 and self.peek("name") and self.toks[1] == ("op", ":"):
            name = self.take("name")
            self.take("op", ":")
        variables = [self.take("name")]
        while self.peek("op", ","):
            self.take("op", ",")
            variables.append(self.take("name"))
        if len(set(variables)) != len(variables):
            raise self.error("binding-site variables must be distinct")
        self.bound = set(variables)
        self.take("op", "->")
        products = []
        if self.i < len(self.toks) and not self.peek("kw", "if"):
            products.append(self.expr())
            while self.peek("op", ","):
                self.take("op", ",")
                products.append(self.expr())
        condition = None
        if self.peek("kw", "if"):
            self.take("kw", "if")
            condition = self.expr()
        if self.i != len(self.toks):
            raise self.error(f"unexpected {self.toks[self.i][1]!r}")
        return name, variables, products, condition

    def expr(self):
        node = self.conjunction()
        while self.peek("kw", "or"):
            self.take("kw", "or")
            node = ("or", node, self.conjunction())
        return node

    def conjunction(self):
        node = self.negation()
        while self.peek("kw", "and"):
            self.take("kw", "and")
            node = ("and", node, self.negation())
        return node

    def negation(self):
        if self.peek("kw", "not"):
            self.take("kw", "not")
            return ("not", self.negation())
        return self.comparison()

    def comparison(self):
        node = self.sum()
        if self.i < len(self.toks) and self.toks[self.i][0] == "op" and self.toks[self.i][1] in _COMPARE:
            op = self.take("op")
            node = ("bin", op, node, self.sum())
        return node

    def sum(self):
        node = self.term()
        while self.peek("op", "+") or self.peek("op", "-"):
            op = self.take("op")
            node = ("bin", op, node, self.term())
        return node

    def term(self):
        node = self.unary()
        while self.peek("op", "*") or self.peek("op", "/") or self.peek("op", "%"):
            op = self.take("op")
            node = ("bin", op, node, self.unary())
        return node

    def unary(self):
        if self.peek("op", "-"):
            self.take("op", "-")
            return ("neg", self.unary())
        return self.atom()

    def atom(self):
        if self.peek("int"):
            return ("const", self.take("int"))
        if self.peek("kw", "true") or self.peek("kw", "false"):
            return ("const", self.take("kw") == "true")
        if self.peek("op", "("):
            self.take("op", "(")
            node = self.expr()
            self.take("op", ")")
            return node
        if self.peek("name"):
            name = self.take("name")
            if self.peek("op", "("):
                if name not in _FUNCTIONS:
                    raise self.error(f"unknown function {name!r}; available: {sorted(_FUNCTIONS)}")
                self.take("op", "(")
                args = [self.expr()]
                while self.peek("op", ","):
                    self.take("op", ",")
                    args.append(self.expr())
                self.take("op", ")")
                low, high = _FUNCTIONS[name]
                if len(args) < low or (high is not None and len(args) > high):
                    raise self.error(f"wrong number of arguments for {name}")
                return ("call", name, args)
            if name not in self.bound:
                raise self.error(f"unknown variable {name!r}; binding sites are {sorted(self.bound)}")
            return ("var", name)
        found = self.toks[self.i][1] if self.i < len(self.toks) else "end of rule"
        raise self.error(f"expected an expression, found {found!r}")


def _int(v):
    if type(v) is not int:
        raise _Mismatch
    return v


def _bool(v):
    if type(v) is not bool:
        raise _Mismatch
    return v


def _floordiv(a, b):
    if _int(b) == 0:
        raise _Mismatch
    return _int(a) // b


def _mod(a, b):
    if _int(b) == 0:
        raise _Mismatch
    return _int(a) % b


_BINARY = {
    "+": lambda a, b: _int(a) + _int(b),
    "-": lambda a, b: _int(a) - _int(b),
    "*": lambda a, b: _int(a) * _int(b),
    "/": _floordiv,
    "%": _mod,
    "==": lambda a, b: type(a) is type(b) and a == b,
    "!=": lambda a, b: not (type(a) is type(b) and a == b),
    "<": lambda a, b: _int(a) < _int(b),
    "<=": lambda a, b: _int(a) <= _int(b),
    ">": lambda a, b: _int(a) > _int(b),
    ">=": lambda a, b: _int(a) >= _int(b),
}
_CALLS = {
    "min": lambda args: min(_int(a) for a in args),
    "max": lambda args: max(_int(a) for a in args),
    "abs": lambda args: abs(_int(args[0])),
}


def _interpret(node, env):
    kind = node[0]
    if kind == "const":
        return node[1]
    if kind == "var":
        return env[node[1]]
    if kind == "neg":
        return -_int(_interpret(node[1], env))
    if kind == "not":
        return not _bool(_interpret(node[1], env))
    if kind == "and":
        return _bool(_interpret(node[1], env)) and _bool(_interpret(node[2], env))
    if kind == "or":
        return _bool(_interpret(node[1], env)) or _bool(_interpret(node[2], env))
    if kind == "bin":
        return _BINARY[node[1]](_interpret(node[2], env), _interpret(node[3], env))
    if kind == "call":
        return _CALLS[node[1]]([_interpret(a, env) for a in node[2]])
    raise AssertionError(kind)


def parse_expression_rule(src: str):
    """Return (name or None, arity, apply) for an expression rule."""
    name, variables, products, condition = _Parser(src).rule()

    def apply(ctx, rng, *mols):
        env = dict(zip(variables, mols))
        try:
            if condition is not None and not _bool(_interpret(condition, env)):
                return list(mols)
            out = [_interpret(e, env) for e in products]
        except _Mismatch:
            return list(mols)
        if any(type(v) is bool for v in out):
            return list(mols)
        return out

    return name, len(variables), apply


def parse_rules(rules: Any) -> list[tuple[Rule, int]]:
    """The rule multiset {rule: multiplicity} -> [(Rule, multiplicity)], merged by rule id."""
    if not isinstance(rules, dict) or not rules:
        raise ValueError(f"rules must be a non-empty object {{rule: multiplicity}}, got {rules!r}")
    parsed: dict[str, list] = {}
    unnamed = 0
    for key, mult in rules.items():
        if not isinstance(mult, int) or isinstance(mult, bool) or mult < 1:
            raise ValueError(f"rules[{key!r}]: multiplicity must be a positive integer, got {mult!r}")
        if "->" in key:
            name, arity, apply = parse_expression_rule(key)
            if name is None:
                unnamed += 1
                name = f"expr{unnamed}"
            rule = Rule(f"rule:{name}", " ".join(key.split()), arity, apply)
            if rule.id in parsed:  # also catches a clash with a named rule given earlier
                raise ValueError(f"rules: duplicate rule name {name!r}")
            parsed[rule.id] = [rule, mult]
            continue
        match = re.fullmatch(r"\s*(?:self\.)?(\w+)\s*(?:\((.*)\))?\s*", key)
        if match is None or match.group(1) not in RULES:
            raise ValueError(
                f"unknown rule {key!r}: named rules are {sorted(RULES)} (optionally written "
                f"as their rule molecule, e.g. 'divrule(m1, m2)'); an expression rule needs '->'. "
                "Rule strings are never executed as code.")
        name = match.group(1)
        rule_text, fn, stochastic, tours = RULES[name]
        arity = binding_sites(rule_text)
        if match.group(2) is not None and len(match.group(2).split(",")) != arity:
            raise ValueError(f"rule {key!r} has {len(match.group(2).split(','))} binding sites, "
                             f"but {name} takes {arity}: write {rule_text!r}")
        rid = f"rule:{name}"
        if rid in parsed and parsed[rid][0].apply is not fn:
            raise ValueError(f"rules: duplicate rule name {name!r}")
        if rid in parsed:
            parsed[rid][1] += mult
        else:
            parsed[rid] = [Rule(rid, rule_text, arity, fn, stochastic, tours), mult]
    return [(r, m) for r, m in parsed.values()]


# ---------------------------------------------------------------------------
# The reactor: HighOrderChem.iterate
# ---------------------------------------------------------------------------
class HighOrderChem:
    """Port of PyCellChemistry HighOrderChem (book appendix, figure 3).

    rset and mset are lists used as multisets; expelrnd removes a uniformly
    random molecule (Multiset.expelrnd picks each molecule with probability
    1/total).
    """

    def __init__(self, rules: list[tuple[Rule, int]], data, rng, ctx=None):
        self.rset = [rule for rule, mult in rules for _ in range(mult)]
        self.mset = list(data)
        self.rng = rng
        self.ctx = ctx

    def _expelrnd(self, pool):
        i = int(self.rng.integers(len(pool)))
        mol = pool[i]
        pool[i] = pool[-1]
        pool.pop()
        return mol

    def iterate(self):
        """One iteration; returns (rule, educts, products) like the reference code."""
        if not self.rset:
            return None, [], []
        rule = self._expelrnd(self.rset)
        educts, products = [], []
        if len(self.mset) >= rule.arity:
            for _ in range(rule.arity):
                educts.append(self._expelrnd(self.mset))
            products = list(rule.apply(self.ctx, self.rng, *educts))
            self.mset.extend(products)
        self.rset.append(rule)  # reinject rule for reuse
        return rule, educts, products


def is_effective(educts, products) -> bool:
    """HighOrderChem.is_effective: the product multiset differs from the educt multiset."""
    if not educts and not products:
        return False
    return Counter(educts) != Counter(products)


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------
def _data(p, rng, graph):
    if p.data:
        return [molecule(v) for v in p.data]
    if p.init == "tours":
        return [graph.random_tour(rng) for _ in range(p.M)]
    if p.maxn < p.minn:
        raise ValueError(f"maxn must be >= minn, got minn={p.minn}, maxn={p.maxn}")
    # NumberChemHO: np.random.randint(minn, maxn+1), both ends inclusive.
    return [int(v) for v in rng.integers(p.minn, p.maxn + 1, size=p.M)]


def _setup(p, rng):
    rules = parse_rules(p.rules)
    uses_tours = any(r.tours for r, _ in rules)
    graph = TourGraph(p.cities) if uses_tours or (p.init == "tours" and not p.data) else None
    data = _data(p, rng, graph)
    if uses_tours:
        bad = [to_json(m) for m in data if isinstance(m, tuple) and not graph.is_tour(m)]
        if bad:
            raise ValueError(f"tour rules need list molecules that are permutations of "
                             f"0..{p.cities - 1} (cities={p.cities}), got {bad[0]!r}")
    return rules, data, graph


def generate(p, rng) -> Network:
    """Every reaction reachable from the distinct data molecules under deterministic rules."""
    rules, data, graph = _setup(p, rng)
    return _closure(p, rules, data, graph)


def evolve(p, rng):
    """The book's HighOrderChem.iterate loop, a frame every generation of len(data) iterations."""
    rules, data, graph = _setup(p, rng)
    if not data:
        raise ValueError("the data multiset is empty")
    return (yield from _soup(p, rng, rules, data, graph))


def _rule_species(rules):
    return [Species(r.id, structure=r.text) for r, _ in rules]


def _data_species(molecules):
    return [Species(sid(m), structure=text(m)) for m in sorted(molecules, key=order)]


def _counts(molecules) -> dict[str, int]:
    return {sid(m): c for m, c in sorted(Counter(molecules).items(), key=lambda kv: order(kv[0]))}


def _soup(p, rng, rules, data, graph):
    chem = HighOrderChem(rules, data, rng, graph)
    size = len(data)
    track_primes = all(r.id == "rule:divrule" for r, _ in rules) and all(isinstance(m, int) for m in data)
    seen = dict.fromkeys(data)
    tally = Tally()
    draws: Counter = Counter()
    idle = 0

    def frame(t):
        rule_counts = Counter(r.id for r in chem.rset)
        state = {r.id: float(rule_counts[r.id]) for r, _ in rules if rule_counts[r.id]}
        state.update({s: float(c) for s, c in _counts(chem.mset).items()})
        observables = {}
        if track_primes and chem.mset:
            observables["prime_fraction"] = sum(map(is_prime, chem.mset)) / len(chem.mset)
        return Frame(t=float(t), state=state, fired=tally.flush(), observables=observables)

    yield frame(0)
    for t in ticks(p.iterations):
        i = t - 1
        rule, educts, products = chem.iterate()
        draws[rule.id] += 1
        if len(educts) < rule.arity:
            idle += 1
        elif is_effective(educts, products):
            tally.add([rule.id, *map(sid, educts)], [rule.id, *map(sid, products)])
            seen.update(dict.fromkeys(products))
        if t % size == 0 or t == p.iterations:
            yield frame(t)

    reactions = [Reaction.of(list(lhs), list(rhs), count=count) for lhs, rhs, count in tally.reactions()]
    start = set(data)
    analysis = {
        "iterations": p.iterations,
        "generation_size": size,
        "effective_collisions": sum(r.count for r in reactions),
        "idle_draws": idle,
        "rule_draws": {r.id: draws.get(r.id, 0) for r, _ in rules},
        "new_molecules": [sid(m) for m in sorted(seen, key=order) if m not in start],
    }
    final_rules = Counter(r.id for r in chem.rset)
    return Network(
        species=_rule_species(rules) + _data_species(seen),
        reactions=reactions,
        status="observed",
        initial_state={**{r.id: m for r, m in rules}, **_counts(data)},
        extras={
            "analysis": analysis,
            "final_state": {**{r.id: final_rules[r.id] for r, _ in rules}, **_counts(chem.mset)},
            "rules": {r.id: r.text for r, _ in rules},
        },
    )


def _closure(p, rules, data, graph):
    stochastic = [r.id for r, _ in rules if r.stochastic]
    if stochastic:
        raise ValueError(f"the closure (generate_network) needs deterministic rules; {stochastic} draw "
                         "random numbers, use chemart.evolve")
    by_arity: dict[int, list[Rule]] = {}
    for r, _ in rules:
        by_arity.setdefault(r.arity, []).append(r)
    recorded: dict[tuple, tuple] = {}

    def react(*mols):
        outs = []
        for rule in by_arity[len(mols)]:
            out = tuple(rule.apply(graph, None, *mols))
            if Counter(out) == Counter(mols):
                continue
            key = (rule.id, frozenset(Counter(mols).items()), frozenset(Counter(out).items()))
            recorded.setdefault(key, (rule, mols, out))
            outs.append(out)
        return outs or None

    seed = sorted(set(data), key=order)
    molecules, _, status = expand(react, seed, arity=sorted(by_arity), max_species=p.max_species,
                                  ordered=True, alternatives=True)
    known = set(molecules)
    reactions = [
        Reaction.of([rule.id, *map(sid, lhs)], [rule.id, *map(sid, rhs)])
        for rule, lhs, rhs in recorded.values() if all(m in known for m in rhs)
    ]
    return Network(
        species=_rule_species(rules) + _data_species(molecules),
        reactions=reactions,
        status=status,
        extras={
            "seed": [sid(m) for m in seed],
            "rules": {r.id: r.text for r, _ in rules},
        },
    )
