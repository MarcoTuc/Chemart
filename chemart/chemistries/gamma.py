"""Gamma: programming by multiset transformation (Banatre & Le Metayer).

Catalog id: gamma. Book 9.2; programs from Banatre, Fradet & Le Metayer,
"Gamma and the chemical reaction model: fifteen years after" (LNCS 2235, 2001)
and Banatre, Fradet & Radenac, "Principles of chemical programming" (RULE'04).

The only data structure is a multiset. A reaction ``x1, ..., xn -> A <= R``
removes n elements satisfying the reaction condition R and adds the elements
produced by the action A. The Gamma operator applies reactions, in any order
and to any tuple, until none applies: the multiset is then stable (inert) and
is the result. Composed programs ``P2 o P1`` run P1 to stability, then P2.

Elements are integers or tuples of elements. Species ids are ``n<k>`` for the
integer k and the canonical text ``(1,5)`` for a tuple; ``structure`` holds the
canonical text.

Program text (never evaluated as Python). One reaction per line::

    # comment
    name: pattern, ..., pattern -> expr, ..., expr if condition
    then                      # sequential composition: next stage
    other: x -> x - 1, x - 2 if x > 1

- name and ``if condition`` are optional; the book's ``→``, ``⇐``, ``≤``, ``≥``,
  ``≠`` and ``=`` are accepted. An empty action (or ``{}``) removes the elements.
- patterns: a variable, an integer, or a tuple ``(p, p, ...)``; a repeated
  variable must bind equal elements.
- expressions: integers, variables, tuples ``(e, e, ...)``, ``+ - *``, ``/``
  (floor division; ``[e]`` is plain grouping, for the papers' integer part),
  ``mod`` or ``%``, comparisons ``== != < <= > >=``, ``and or not``,
  ``true false``, and the functions ``multiple(x, y)`` (x is a multiple of y),
  ``min``, ``max``, ``abs``.
- an ill-typed expression (a tuple compared with ``<``, division by zero) makes
  the reaction inapplicable to that tuple.

The network is the closure (chemart.expand.expand) of the initial elements
under each stage in turn; extras.analysis holds the stable multiset of one
sampled run and, within a state budget, every stable multiset reachable.
"""

from __future__ import annotations

import re
from collections import Counter
from itertools import product
from math import prod

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

# --- the published programs ------------------------------------------------------
PROGRAMS = {
    # book eq. 9.7; Gamma15 sec. 1
    "max": ("max: x, y -> y if x <= y", [3, 8, 1, 8, 4]),
    # Gamma15 sec. 1 (exchange sort of (index, value) pairs)
    "sort": ("sort: (i, x), (j, y) -> (i, y), (j, x) if i > j and x < y",
             [[1, 8], [2, 3], [3, 6], [4, 1], [5, 9]]),
    # Gamma15 sec. 2.1: primes(n) = rem(iota({(2, n)}))
    "primes": ("iota_split: (x, y) -> (x, [(x + y) / 2]), ([(x + y) / 2] + 1, y) if x != y\n"
               "iota_leaf: (x, y) -> x if x == y\n"
               "then\n"
               "rem: x, y -> y if multiple(x, y)",
               [[2, 20]]),
    # Gamma15 sec. 2.1: fib(n) = add(dec1(n))
    "fibonacci": ("dec1_split: x -> x - 1, x - 2 if x > 1\n"
                  "dec1_zero: x -> 1 if x == 0\n"
                  "then\n"
                  "add: x, y -> x + y if true",
                  [5]),
    # Gamma15 sec. 2.1: maxss(M) = maxg(maxl(M)), elements (i, x, s) with s = x initially
    "max-segment-sum": ("maxl: (i, x, s), (j, y, t) -> (i, x, s), (j, y, s + y) if j == i + 1 and s + y > t\n"
                        "then\n"
                        "maxg: (i, x, s), (j, y, t) -> (j, y, t) if t > s",
                        [[1, 3, 3], [2, -4, -4], [3, 5, 5], [4, -1, -1], [5, 2, 2]]),
    # RULE'04 fig. 1 and the York abstract: majority element
    "majority": ("maj: x, y -> if x != y", [1, 2, 1, 3, 1, 1, 2]),
    # York abstract / RULE'04: <<2..10, prime>, replace <prime, w> by <w, max>>, as primes then max
    "largest-prime": ("primes: x, y -> y if multiple(x, y)\n"
                      "then\n"
                      "max: x, y -> x if x >= y",
                      [2, 3, 4, 5, 6, 7, 8, 9, 10]),
}


# --- elements --------------------------------------------------------------------
class _Mismatch(Exception):
    """An ill-typed evaluation: the reaction does not apply to this tuple."""


def text(v) -> str:
    """Canonical text of an element: 5, (1,5), ((1,2),3)."""
    if isinstance(v, tuple):
        return "(" + ",".join(text(e) for e in v) + ")"
    return str(v)


def sid(v) -> str:
    return f"n{v}" if isinstance(v, int) else text(v)


def order(v):
    """Sort key: integers by value, then tuples by length and components."""
    return (0, v) if isinstance(v, int) else (1, len(v), tuple(order(e) for e in v))


def to_json(v):
    return [to_json(e) for e in v] if isinstance(v, tuple) else v


def element(value, where: str = "multiset"):
    """JSON value (int or list of >= 2 elements) -> element (int or tuple)."""
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, (list, tuple)) and len(value) >= 2:
        return tuple(element(e, where) for e in value)
    raise ValueError(f"{where} elements must be integers or lists of at least two elements, got {value!r}")


# --- the reaction language -------------------------------------------------------
_TOKEN = re.compile(r"\s*(?:(?P<int>\d+)|(?P<name>[A-Za-z_][A-Za-z0-9_]*'*)"
                    r"|(?P<op>->|→|⇐|<=|>=|==|!=|≤|≥|≠|∅|[-+*/%(),:<>=\[\]{}]))")
_OP_ALIASES = {"→": "->", "≤": "<=", "≥": ">=", "≠": "!=", "=": "==", "%": "mod"}
_KEYWORDS = {"if", "and", "or", "not", "true", "false", "mod", "then"}
_CMP = ("==", "!=", "<", "<=", ">", ">=")


def _tokens(line: str) -> list[tuple[str, object]]:
    out, pos, line = [], 0, line.rstrip()
    while pos < len(line):
        m = _TOKEN.match(line, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot read {line[pos:].strip()!r} in reaction {line.strip()!r}")
        pos = m.end()
        if m.group("int") is not None:
            out.append(("int", int(m.group("int"))))
        elif m.group("name") is not None:
            name = m.group("name")
            name = name.lower() if name in ("True", "False") else name
            out.append(("kw", name) if name in _KEYWORDS else ("name", name))
        else:
            op = m.group("op")
            out.append(("kw", "if") if op == "⇐" else ("kw", "mod") if op == "%" else
                       ("op", _OP_ALIASES.get(op, op)))
    return out


def _int(v):
    if type(v) is not int:
        raise _Mismatch
    return v


def _bool(v):
    if type(v) is not bool:
        raise _Mismatch
    return v


def _element_value(v):
    if type(v) is int or type(v) is tuple:
        return v
    raise _Mismatch


def _div(a, b):
    if _int(b) == 0:
        raise _Mismatch
    return _int(a) // b


def _mod(a, b):
    if _int(b) == 0:
        raise _Mismatch
    return _int(a) % b


def _eq(a, b):
    return type(a) is type(b) and a == b


_OPS = {
    "+": lambda a, b: _int(a) + _int(b),
    "-": lambda a, b: _int(a) - _int(b),
    "*": lambda a, b: _int(a) * _int(b),
    "/": _div,
    "mod": _mod,
    "==": _eq,
    "!=": lambda a, b: not _eq(a, b),
    "<": lambda a, b: _int(a) < _int(b),
    "<=": lambda a, b: _int(a) <= _int(b),
    ">": lambda a, b: _int(a) > _int(b),
    ">=": lambda a, b: _int(a) >= _int(b),
}
_FUNCS = {
    "multiple": (2, lambda x, y: _int(y) != 0 and _int(x) % y == 0),
    "min": (2, lambda x, y: min(_int(x), _int(y))),
    "max": (2, lambda x, y: max(_int(x), _int(y))),
    "abs": (1, lambda x: abs(_int(x))),
}


def _binary(op, left, right):
    if op == "or":
        return lambda env: _bool(left(env)) or _bool(right(env))
    if op == "and":
        return lambda env: _bool(left(env)) and _bool(right(env))
    f = _OPS[op]
    return lambda env: f(left(env), right(env))


class _Parser:
    def __init__(self, tokens, line):
        self.tokens, self.i, self.line, self.variables = tokens, 0, line, set()

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else (None, None)

    def next(self):
        tok = self.peek()
        if tok[0] is None:
            self.error("unexpected end")
        self.i += 1
        return tok

    def accept(self, kind, value):
        if self.peek() == (kind, value):
            self.i += 1
            return True
        return False

    def expect(self, kind, value):
        if not self.accept(kind, value):
            self.error(f"expected {value!r} but found {self.peek()[1]!r}")

    def error(self, message):
        raise ValueError(f"{message} in reaction {self.line.strip()!r}")

    # patterns
    def pattern(self):
        if self.accept("op", "("):
            items = [self.pattern()]
            while self.accept("op", ","):
                items.append(self.pattern())
            self.expect("op", ")")
            if len(items) < 2:
                self.error("a tuple pattern needs at least two components")
            return ("tuple", tuple(items))
        negative = self.accept("op", "-")
        kind, value = self.next()
        if kind == "int":
            return ("const", -value if negative else value)
        if kind == "name" and not negative:
            self.variables.add(value)
            return ("var", value)
        self.error(f"a pattern is a variable, an integer or a tuple, not {value!r}")

    # expressions
    def expr(self):
        left = self.conj()
        while self.accept("kw", "or"):
            left = _binary("or", left, self.conj())
        return left

    def conj(self):
        left = self.neg()
        while self.accept("kw", "and"):
            left = _binary("and", left, self.neg())
        return left

    def neg(self):
        if self.accept("kw", "not"):
            inner = self.neg()
            return lambda env: not _bool(inner(env))
        return self.comparison()

    def comparison(self):
        left = self.additive()
        kind, value = self.peek()
        if kind == "op" and value in _CMP:
            self.i += 1
            return _binary(value, left, self.additive())
        return left

    def additive(self):
        left = self.term()
        while self.peek() in (("op", "+"), ("op", "-")):
            op = self.next()[1]
            left = _binary(op, left, self.term())
        return left

    def term(self):
        left = self.unary()
        while self.peek() in (("op", "*"), ("op", "/"), ("kw", "mod")):
            op = self.next()[1]
            left = _binary(op, left, self.unary())
        return left

    def unary(self):
        if self.accept("op", "-"):
            inner = self.unary()
            return lambda env: -_int(inner(env))
        return self.primary()

    def primary(self):
        kind, value = self.next()
        if kind == "int":
            return lambda env: value
        if kind == "kw" and value in ("true", "false"):
            truth = value == "true"
            return lambda env: truth
        if kind == "name":
            if self.accept("op", "("):
                args = [self.expr()]
                while self.accept("op", ","):
                    args.append(self.expr())
                self.expect("op", ")")
                if value not in _FUNCS:
                    self.error(f"unknown function {value!r} (available: {', '.join(sorted(_FUNCS))})")
                n, f = _FUNCS[value]
                if len(args) != n:
                    self.error(f"{value} takes {n} arguments")
                return lambda env: f(*(a(env) for a in args))
            if value not in self.variables:
                self.error(f"variable {value!r} does not appear in the patterns")
            return lambda env: env[value]
        if kind == "op" and value in ("(", "["):
            close = ")" if value == "(" else "]"
            items = [self.expr()]
            while close == ")" and self.accept("op", ","):
                items.append(self.expr())
            self.expect("op", close)
            if len(items) == 1:
                return items[0]
            return lambda env: tuple(_element_value(f(env)) for f in items)
        self.error(f"unexpected {value!r}")


def _match(pattern, value, env) -> bool:
    kind, p = pattern
    if kind == "var":
        if p in env:
            return _eq(env[p], value)
        env[p] = value
        return True
    if kind == "const":
        return type(value) is int and value == p
    return (type(value) is tuple and len(value) == len(p)
            and all(_match(q, v, env) for q, v in zip(p, value)))


class Rule:
    """One Gamma reaction, e.g. Rule('max: x, y -> y if x <= y')."""

    def __init__(self, line: str, name: str = "r"):
        tokens = _tokens(line)
        body = line.strip()
        if len(tokens) >= 2 and tokens[0][0] == "name" and tokens[1] == ("op", ":"):
            name = str(tokens[0][1])
            tokens = tokens[2:]
            body = body.split(":", 1)[1].strip()
        self.name, self.text = name, f"{name}: {' '.join(body.split())}"
        p = _Parser(tokens, line)
        self.patterns = [p.pattern()]
        while p.accept("op", ","):
            self.patterns.append(p.pattern())
        p.expect("op", "->")
        self.actions = []
        if p.accept("op", "{"):
            p.expect("op", "}")
        elif not p.accept("op", "∅") and p.peek()[0] is not None and p.peek() != ("kw", "if"):
            self.actions.append(p.expr())
            while p.accept("op", ","):
                self.actions.append(p.expr())
        self.condition = p.expr() if p.accept("kw", "if") else None
        if p.peek()[0] is not None:
            p.error(f"unexpected {p.peek()[1]!r}")

    @property
    def arity(self) -> int:
        return len(self.patterns)

    def apply(self, *lhs) -> tuple | None:
        """The action's elements for this ordered tuple, or None if the reaction does not apply."""
        if len(lhs) != self.arity:
            return None
        env: dict = {}
        if not all(_match(q, v, env) for q, v in zip(self.patterns, lhs)):
            return None
        try:
            if self.condition is not None and not _bool(self.condition(env)):
                return None
            return tuple(_element_value(a(env)) for a in self.actions)
        except _Mismatch:
            return None


class Stage:
    """A parallel composition of reactions (one stage of a sequential composition)."""

    def __init__(self, rules: list[Rule]):
        self.rules = rules
        self.arities = sorted({r.arity for r in rules})
        self._memo: dict[tuple, list] = {}

    def outcomes(self, lhs: tuple) -> list[tuple[tuple, str]]:
        """Distinct (products, rule names) of every reaction applicable to the ordered tuple."""
        if lhs not in self._memo:
            found: dict[frozenset, list] = {}
            left = Counter(lhs)
            for r in self.rules:
                rhs = r.apply(*lhs)
                if rhs is None or Counter(rhs) == left:
                    continue
                entry = found.setdefault(frozenset(Counter(rhs).items()), [rhs, []])
                if r.name not in entry[1]:
                    entry[1].append(r.name)
            self._memo[lhs] = [(rhs, "|".join(names)) for rhs, names in found.values()]
        return self._memo[lhs]


def parse_program(program: str) -> list[Stage]:
    """Program text -> stages; ``then`` lines separate sequentially composed stages."""
    stages: list[list[Rule]] = [[]]
    count = 0
    for raw in program.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line == "then":
            if not stages[-1]:
                raise ValueError("'then' must separate two non-empty groups of reactions")
            stages.append([])
            continue
        count += 1
        stages[-1].append(Rule(line, f"r{count}"))
    if not stages[-1]:
        raise ValueError("a program needs at least one reaction (and no trailing 'then')")
    return [Stage(rules) for rules in stages]


# --- Gamma semantics -----------------------------------------------------------------
def _falling(n: int, k: int) -> int:
    return prod(range(n - k + 1, n + 1))


def applicable(stage: Stage, counts: Counter) -> list[tuple[tuple, tuple, str, int]]:
    """Every (reactants, products, rule, weight) applicable to the multiset.

    weight = number of ordered tuples of distinct element occurrences with those values.
    """
    values = sorted(counts, key=order)
    size = sum(counts.values())
    out = []
    for k in stage.arities:
        if size < k:
            continue
        for lhs in product(values, repeat=k):
            need = Counter(lhs)
            if any(counts[v] < c for v, c in need.items()):
                continue
            found = stage.outcomes(lhs)
            if found:
                weight = prod(_falling(counts[v], c) for v, c in need.items())
                out.extend((lhs, rhs, names, weight) for rhs, names in found)
    return out


def _fire(counts: Counter, lhs, rhs) -> Counter:
    new = counts.copy()
    for v in lhs:
        new[v] -= 1
        if not new[v]:
            del new[v]
    new.update(rhs)
    return new


def run(stages: list[Stage], initial: Counter, rng, max_steps: int):
    """One execution of the Gamma operator on every stage in turn.

    Each step fires one applicable reaction, chosen with probability proportional
    to the number of element tuples it applies to. Returns (final multiset,
    reactions fired, whether a stable state was reached, multiset after each
    completed stage).
    """
    counts, steps, stable, after = Counter(initial), 0, True, []
    for stage in stages:
        while True:
            candidates = applicable(stage, counts)
            if not candidates:
                break
            if steps >= max_steps:
                stable = False
                break
            x = rng.random() * sum(c[3] for c in candidates)
            for lhs, rhs, _, weight in candidates:
                if x < weight:
                    break
                x -= weight
            counts = _fire(counts, lhs, rhs)
            steps += 1
        if not stable:
            break
        after.append(counts)
    return counts, steps, stable, after


def _key(counts: Counter) -> tuple:
    return tuple(sorted(counts.items(), key=lambda vc: order(vc[0])))


def results(stages: list[Stage], initial: Counter, max_states: int) -> list[Counter] | None:
    """Every stable multiset reachable under any order of reactions (the Gamma relation).

    None when more than max_states multisets would have to be visited.
    """
    frontier, visited = {_key(Counter(initial))}, 0
    for stage in stages:
        finals, seen, stack = set(), set(frontier), list(frontier)
        while stack:
            state = stack.pop()
            visited += 1
            if visited > max_states:
                return None
            counts = Counter(dict(state))
            candidates = applicable(stage, counts)
            if not candidates:
                finals.add(state)
            for lhs, rhs, _, _ in candidates:
                new = _key(_fire(counts, lhs, rhs))
                if new not in seen:
                    seen.add(new)
                    stack.append(new)
        frontier = finals
    return [Counter(dict(s)) for s in sorted(frontier, key=lambda s: [(order(v), c) for v, c in s])]


def closure(stages: list[Stage], seed: list, max_species: int):
    """Closure of each stage in turn, seeded with every species found so far.

    Returns (species, [(reactants, products, rule names)], status).
    """
    species, found, index, status = list(seed), [], {}, "complete"
    for stage in stages:
        def react(*lhs, stage=stage):
            return [rhs for rhs, _ in stage.outcomes(lhs)]

        species, reactions, st = expand(react, species, arity=stage.arities, max_species=max_species,
                                        ordered=True, alternatives=True)
        if st == "truncated":
            status = "truncated"
        for lhs, rhs in reactions:
            right = Counter(rhs)
            names = next(n for r, n in stage.outcomes(lhs) if Counter(r) == right)
            key = (frozenset(Counter(lhs).items()), frozenset(right.items()))
            if key in index:
                entry = found[index[key]]
                entry[2] = "|".join(dict.fromkeys(entry[2].split("|") + names.split("|")))
                continue
            index[key] = len(found)
            found.append([lhs, rhs, names])
    return species, [tuple(f) for f in found], status


# --- generator -------------------------------------------------------------------------
def _inputs(p):
    if p.program == "custom":
        if not isinstance(p.rules, str) or not p.rules.strip():
            raise ValueError("program 'custom' needs rules: the reaction text, e.g. 'max: x, y -> y if x <= y'")
        if not p.multiset:
            raise ValueError("program 'custom' needs a non-empty multiset")
        text_, default = p.rules, None
    else:
        if p.rules:
            raise ValueError("rules is only used with program='custom'")
        text_, default = PROGRAMS[p.program]
    raw = p.multiset if p.multiset else default
    if not isinstance(raw, list):
        raise ValueError(f"multiset must be a list, got {raw!r}")
    return text_, [element(v) for v in raw]


def _multiset(counts: Counter) -> dict[str, int]:
    return {sid(v): int(c) for v, c in sorted(counts.items(), key=lambda vc: order(vc[0]))}


def generate(p, rng):
    program, elements = _inputs(p)
    stages = parse_program(program)
    initial = Counter(elements)
    seed = sorted(initial, key=order)

    species, found, status = closure(stages, seed, p.max_species)
    final, steps, stable, after = run(stages, initial, rng, p.max_steps)
    every = results(stages, initial, p.max_states)

    analysis = {
        "final_multiset": _multiset(final),
        "final_values": [to_json(v) for v in sorted(final.elements(), key=order)],
        "stable": stable,
        "steps": steps,
        "stage_multisets": [_multiset(c) for c in after],
        "results": None if every is None else [_multiset(c) for c in every],
        "deterministic": None if every is None else len(every) == 1,
    }
    return Network(
        species=[Species(sid(v), structure=text(v)) for v in species],
        reactions=[Reaction.of([sid(v) for v in lhs], [sid(v) for v in rhs]) for lhs, rhs, _ in found],
        status=status,
        initial_state={sid(v): float(initial[v]) for v in seed},
        extras={
            "program": p.program,
            "stages": [[r.text for r in s.rules] for s in stages],
            "reaction_rules": [names for _, _, names in found],
            "analysis": analysis,
        },
    )
