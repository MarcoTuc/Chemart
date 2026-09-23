"""Artificial chemistry as a proof search system: RESAC (book 16.5; Busch 2004).

Catalog id: proof-ac.

Molecules are clauses of first-order predicate logic; the only reaction is
binary resolution (Busch 2004, Def. 1.2.12 and eq. 3.1). Two clauses with a
complementary pair of unifiable literals react and produce their resolvent;
otherwise the collision is elastic. The theory's clauses plus the negated
theorem are the start clauses; a refutation is the construction pathway of
the empty clause [], read off from which reactions built it (book 16.5).

Clause syntax (the `clauses`, `goal` and `target` parameters): clauses are
separated by ';', literals by ',' (or '|'), negation is '~' (also '-' or
'¬'). Variables start with an upper-case letter or '_' (Prolog convention);
predicates, functions and constants start with a lower-case letter or a
digit. `[]` is the empty clause. Example: "~dog(X), howls(X); has(john,s)".

Two faces. `generate` is the level-saturation closure of the start clauses
(thesis 1.3.3 and 6.1): every resolvent of every pair, level by level,
truncated by `max_species`. Premises survive (a theory only grows), so a
reaction is a + b -> a + b + r. `evolve` is RESAC's reactor (thesis
algorithms 3.1-3.3): `multiplicity` copies of each start clause, random
collisions, educt or free replacement, optional inflow of start clauses,
stopped when the target appears; a frame per generation (one collision per
molecule of the reactor).
"""

from __future__ import annotations

import re
from collections import Counter
from itertools import permutations, product
from math import factorial

from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame, ticks

EMPTY = "[]"

# --- built-in problems (thesis tables 6.1 and 6.2) ------------------------------
PROBLEMS = {
    # Gini (1995) via Busch (2004) table 6.1: s1..s5 axioms, s6..s8 the negated goal.
    "gini-1995": (
        "~dog(X), howls(X);"
        "~has(X,Y), ~cat(Y), ~has(X,Z), ~mice(Z);"
        "~lightsleep(X), ~has(X,Y), ~howls(Y);"
        "has(john,s);"
        "cat(s), dog(s)",
        "lightsleep(john); has(john,m); mice(m)",
    ),
    # Loveland (1978) via Busch (2004) table 6.2: existence of right inverses in a group.
    "group-right-inverse": (
        "p(X,Y,f(X,Y));"
        "~p(X,Y,U), ~p(Y,Z,V), ~p(X,V,W), p(U,Z,W);"
        "~p(X,Y,U), ~p(Y,Z,V), ~p(U,Z,W), p(X,V,W);"
        "p(e,Y,Y);"
        "p(g(Y),Y,e)",
        "~p(X,h(X),h(X)), ~p(k(X),Z,X)",
    ),
}

STRATEGIES = ("unrestricted", "set-of-support", "negative", "negative-set-of-support")

# Canonicalisation tries every ordering of literals with the same skeleton up to
# this many orderings; beyond it a single deterministic ordering is used.
_PERMUTATION_LIMIT = 720


# --- terms, literals, clauses ---------------------------------------------------
# A variable is a str; a constant or compound term is (name, args-tuple).
# A literal is (positive, predicate, args-tuple); a clause is a frozenset of literals.

_TOKEN = re.compile(r"\s*(?:(\[\])|([A-Za-z0-9_]+)|(.))")


def _is_var(name: str) -> bool:
    return name[0].isupper() or name[0] == "_"


class _Parser:
    def __init__(self, text: str):
        self.tokens = []
        pos = 0
        text = text.strip()
        while pos < len(text):
            m = _TOKEN.match(text, pos)
            if m.group(0).strip() == "":
                break
            self.tokens.append(m.group(1) or m.group(2) or m.group(3))
            pos = m.end()
        self.i = 0
        self.text = text

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self, expected=None):
        tok = self.peek()
        if tok is None or (expected is not None and tok != expected):
            raise ValueError(f"cannot parse clause text {self.text!r}: expected {expected or 'a token'}, got {tok!r}")
        self.i += 1
        return tok

    def term(self):
        name = self.take()
        if not re.fullmatch(r"[A-Za-z0-9_]+", name):
            raise ValueError(f"cannot parse clause text {self.text!r}: unexpected {name!r}")
        if _is_var(name):
            return name
        return (name, self.args())

    def args(self):
        if self.peek() != "(":
            return ()
        self.take("(")
        out = [self.term()]
        while self.peek() == ",":
            self.take(",")
            out.append(self.term())
        self.take(")")
        return tuple(out)

    def literal(self):
        positive = True
        while self.peek() in ("~", "-", "¬"):
            self.take()
            positive = not positive
        name = self.take()
        if not re.fullmatch(r"[a-z0-9][A-Za-z0-9_]*", name):
            raise ValueError(f"cannot parse clause text {self.text!r}: predicate must start lower-case, got {name!r}")
        return (positive, name, self.args())

    def clause(self):
        if self.peek() == "[]":
            self.take()
            return frozenset()
        lits = [self.literal()]
        while self.peek() in (",", "|"):
            self.take()
            lits.append(self.literal())
        return frozenset(lits)


def parse_clauses(text: str) -> list[frozenset]:
    """Parse ';'-separated clauses; empty pieces are ignored."""
    out = []
    for piece in text.replace("\n", ";").split(";"):
        if piece.strip():
            p = _Parser(piece)
            c = p.clause()
            if p.peek() is not None:
                raise ValueError(f"cannot parse clause {piece.strip()!r}: unexpected {p.peek()!r}")
            out.append(c)
    return out


def term_str(t, ren=None) -> str:
    if isinstance(t, str):
        return ren(t) if ren else t
    name, args = t
    return name + ("(" + ",".join(term_str(a, ren) for a in args) + ")" if args else "")


def literal_str(lit, ren=None) -> str:
    positive, pred, args = lit
    body = pred + ("(" + ",".join(term_str(a, ren) for a in args) + ")" if args else "")
    return body if positive else "~" + body


def _term_vars(t, out: list):
    if isinstance(t, str):
        if t not in out:
            out.append(t)
    else:
        for a in t[1]:
            _term_vars(a, out)


def clause_vars(lits) -> list[str]:
    out: list[str] = []
    for _, _, args in lits:
        for a in args:
            _term_vars(a, out)
    return out


def _rename_term(t, mapping):
    if isinstance(t, str):
        return mapping.get(t, t)
    return (t[0], tuple(_rename_term(a, mapping) for a in t[1]))


def _rename(clause, mapping):
    return frozenset((p, n, tuple(_rename_term(a, mapping) for a in args)) for p, n, args in clause)


def canonical(clause: frozenset) -> tuple[str, frozenset]:
    """Canonical id (identical for alphabetic variants) and the renamed clause."""
    if not clause:
        return EMPTY, clause
    blank = lambda v: "_"  # noqa: E731
    lits = sorted(clause, key=lambda l: (literal_str(l, blank), literal_str(l)))
    groups: list[list] = []
    for lit in lits:
        if groups and literal_str(groups[-1][0], blank) == literal_str(lit, blank):
            groups[-1].append(lit)
        else:
            groups.append([lit])
    count = 1
    for g in groups:
        count *= factorial(len(g))
    orderings = product(*(permutations(g) for g in groups)) if count <= _PERMUTATION_LIMIT else [groups]
    best = None
    for choice in orderings:
        order = [l for g in choice for l in g]
        mapping = {v: f"X{i + 1}" for i, v in enumerate(clause_vars(order))}
        text = ",".join(sorted(literal_str(l, mapping.get) for l in order))
        if best is None or text < best[0]:
            best = (text, mapping)
    return best[0], _rename(clause, best[1])


def symbols(clause) -> int:
    """Number of non-logical symbols: predicates, functions, constants, variables."""
    def count(t):
        return 1 if isinstance(t, str) else 1 + sum(count(a) for a in t[1])
    return sum(1 + sum(count(a) for a in args) for _, _, args in clause)


# --- unification (Robinson, with occurs check) ----------------------------------
def _walk(t, s):
    while isinstance(t, str) and t in s:
        t = s[t]
    return t


def _occurs(v, t, s) -> bool:
    t = _walk(t, s)
    if isinstance(t, str):
        return t == v
    return any(_occurs(v, a, s) for a in t[1])


def unify(a, b, s: dict | None = None) -> dict | None:
    """Most general unifier extending s (triangular form), or None."""
    s = {} if s is None else dict(s)
    stack = [(a, b)]
    while stack:
        x, y = stack.pop()
        x, y = _walk(x, s), _walk(y, s)
        if x == y:
            continue
        if isinstance(x, str):
            if _occurs(x, y, s):
                return None
            s[x] = y
        elif isinstance(y, str):
            if _occurs(y, x, s):
                return None
            s[y] = x
        else:
            if x[0] != y[0] or len(x[1]) != len(y[1]):
                return None
            stack.extend(zip(x[1], y[1]))
    return s


def _subst_term(t, s):
    t = _walk(t, s)
    if isinstance(t, str):
        return t
    return (t[0], tuple(_subst_term(a, s) for a in t[1]))


def substitute(clause, s) -> frozenset:
    return frozenset((p, n, tuple(_subst_term(a, s) for a in args)) for p, n, args in clause)


def is_tautology(clause) -> bool:
    return any((not p, n, args) in clause for p, n, args in clause)


def factors(clause: frozenset) -> list[frozenset]:
    """All proper factors (Def. 1.2.13): unify literals of equal sign and predicate."""
    seen = {canonical(clause)[0]}
    out, todo = [], [clause]
    while todo:
        c = todo.pop()
        lits = sorted(c, key=literal_str)
        for i, a in enumerate(lits):
            for b in lits[i + 1:]:
                if a[0] != b[0] or a[1] != b[1] or len(a[2]) != len(b[2]):
                    continue
                s = unify(("", a[2]), ("", b[2]))
                if s is None:
                    continue
                f = substitute(c, s)
                key = canonical(f)[0]
                if key not in seen:
                    seen.add(key)
                    out.append(f)
                    todo.append(f)
    return out


def binary_resolvents(a: frozenset, b: frozenset) -> list[frozenset]:
    """Binary resolvents of a and b (Def. 1.2.12), one per resolvable literal pair.

    b is standardised apart from a first. Identical literals merge (set union).
    """
    va = set(clause_vars(a))
    mapping = {v: f"{v}_" for v in clause_vars(b)}
    while any(m in va for m in mapping.values()):
        mapping = {v: f"{m}_" for v, m in mapping.items()}
    b = _rename(b, mapping)
    out = []
    for la in sorted(a, key=literal_str):
        for lb in sorted(b, key=literal_str):
            if la[0] == lb[0] or la[1] != lb[1] or len(la[2]) != len(lb[2]):
                continue
            s = unify(("", la[2]), ("", lb[2]))
            if s is None:
                continue
            out.append(substitute((a - {la}) | (b - {lb}), s))
    return out


# --- the reaction ----------------------------------------------------------------
class Chemistry:
    """Res2 of the thesis, with its prover settings, over canonical clause ids."""

    def __init__(self, factoring: bool, max_length: int, strategy: str):
        self.factoring = factoring
        self.max_length = max_length
        self.strategy = strategy
        self.clauses: dict[str, frozenset] = {}
        self._factors: dict[str, list[frozenset]] = {}
        self._res: dict[tuple[str, str], list[str]] = {}

    def add(self, clause: frozenset) -> str:
        cid, renamed = canonical(clause)
        self.clauses.setdefault(cid, renamed)
        return cid

    def variants(self, cid: str) -> list[frozenset]:
        """The clause and, with factoring, its proper factors."""
        if cid not in self._factors:
            c = self.clauses[cid]
            self._factors[cid] = [c] + (factors(c) if self.factoring else [])
        return self._factors[cid]

    def outcomes(self, a: str, b: str) -> list[str]:
        """Resolvents of every resolvable literal pairing (of factors), with repetition.

        Tautologies and clauses longer than max_length are unstable and dropped.
        """
        key = (a, b) if a <= b else (b, a)
        if key not in self._res:
            out = []
            for ca in self.variants(key[0]):
                for cb in self.variants(key[1]):
                    for r in binary_resolvents(ca, cb):
                        if is_tautology(r) or (self.max_length and symbols(r) > self.max_length):
                            continue
                        out.append(self.add(r))
            self._res[key] = out
        return self._res[key]

    def allowed(self, a: str, b: str, support: dict[str, bool]) -> bool:
        if "negative" in self.strategy:
            if not any(all(not p for p, _, _ in self.clauses[c]) for c in (a, b)):
                return False
        if "set-of-support" in self.strategy:
            if not (support[a] or support[b]):
                return False
        return True


def _setup(p):
    """The chemistry, the start clauses (axioms + goal), the target and the support flags."""
    axioms_text, goal_text = _problem(p)
    chem = Chemistry(p.factoring, p.max_length, p.strategy)
    axioms = list(dict.fromkeys(chem.add(c) for c in parse_clauses(axioms_text)))
    goal = list(dict.fromkeys(chem.add(c) for c in parse_clauses(goal_text)))
    start = list(dict.fromkeys(axioms + goal))
    if not start:
        raise ValueError("no start clauses: give clauses (and goal) when problem is custom")
    targets = parse_clauses(p.target)
    if len(targets) != 1:
        raise ValueError(f"target must be exactly one clause (use [] for the empty clause), got {p.target!r}")
    target = chem.add(targets[0])
    if "set-of-support" in p.strategy and not goal:
        raise ValueError("the set-of-support strategies need goal clauses (the negated theorem)")
    support = {c: c in goal for c in start}
    return chem, start, axioms, goal, target, support


def generate(p, rng):
    """The level-saturation closure of the start clauses, cut off by max_species."""
    return _closure(p, *_setup(p))


def _problem(p) -> tuple[str, str]:
    if p.problem == "custom":
        if not p.clauses.strip() and not p.goal.strip():
            raise ValueError("problem custom needs clauses (and optionally goal)")
        return p.clauses, p.goal
    if p.clauses.strip() or p.goal.strip():
        raise ValueError("clauses and goal are only read when problem is custom")
    return PROBLEMS[p.problem]


def _key(lhs, rhs):
    return (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))


def proof_of(target: str, producer: dict[str, tuple[str, str]]) -> list[dict]:
    """The reactions that build target from the start clauses, premises first."""
    steps: list[dict] = []
    done: set[str] = set()
    stack = [(target, False)]
    while stack:
        c, expanded = stack.pop()
        if c in done or c not in producer:
            continue
        if expanded:
            done.add(c)
            a, b = producer[c]
            steps.append({"reactants": [a, b], "product": c})
        else:
            stack.append((c, True))
            a, b = producer[c]
            stack.extend([(b, False), (a, False)])
    return steps


def _analysis(chem, target, producer, start, extra):
    proved = target in producer or target in start
    proof = proof_of(target, producer) if proved else []
    return {
        "target": target,
        "proved": proved,
        "proof": proof,
        "proof_length": len(proof),
        "longest_clause_in_proof": max((len(chem.clauses[s["product"]]) for s in proof), default=0),
        **extra,
    }


def _species(chem, ids):
    return [Species(c, structure=" | ".join(sorted(literal_str(l) for l in chem.clauses[c])) or EMPTY)
            for c in ids]


def _closure(p, chem, start, axioms, goal, target, support):
    species = list(start)
    known = set(species)
    level = {c: 0 for c in species}
    producer: dict[str, tuple[str, str]] = {}
    reactions: dict[tuple, Reaction] = {}
    tried: set[tuple[str, str]] = set()
    truncated = False
    frontier = list(species)
    # Once the budget is full nothing new can be added, so the saturation stops
    # (the remaining pairs of first-order clauses can be very expensive).
    while frontier and not truncated:
        snapshot = list(species)
        index = {c: i for i, c in enumerate(snapshot)}
        front = set(frontier)
        nxt: list[str] = []
        for x in frontier:
            if truncated:
                break
            for y in snapshot:
                if y in front and index[y] < index[x]:
                    continue
                pair = (x, y) if x <= y else (y, x)
                if pair in tried or not chem.allowed(x, y, support):
                    continue
                tried.add(pair)
                supported = support[x] or support[y]
                for r in dict.fromkeys(chem.outcomes(x, y)):
                    if r not in known:
                        if len(species) >= p.max_species:
                            truncated = True
                            continue
                        species.append(r)
                        known.add(r)
                        level[r] = max(level[x], level[y]) + 1
                        producer[r] = pair
                        support[r] = supported
                        nxt.append(r)
                    elif supported and not support[r]:
                        support[r] = True
                        nxt.append(r)
                    lhs = (pair[0], pair[1])
                    rhs = (pair[0], pair[1], r)
                    reactions.setdefault(_key(lhs, rhs), Reaction.of(lhs, rhs))
        frontier = list(dict.fromkeys(nxt))
    return Network(
        species=_species(chem, species),
        reactions=list(reactions.values()),
        status="truncated" if truncated else "complete",
        extras={
            "axioms": axioms,
            "goal": goal,
            "analysis": _analysis(chem, target, producer, start, {
                "target_level": level.get(target),
                "levels": {c: level[c] for c in species},
            }),
        },
    )


def evolve(p, rng):
    """RESAC's reactor (thesis algorithms 3.1-3.3): a frame per generation, stopped at the target.

    The reactor is its own loop rather than chemart.soup.stir: the resolvent
    overwrites one reactant or a random position in place, and start clauses
    flow in between collisions.
    """
    if p.inflow_rate > 0 and p.elastic_inflow:
        raise ValueError("choose one inflow: set elastic_inflow=false to use inflow_rate > 0")
    chem, start, axioms, goal, target, support = _setup(p)
    if len(start) * p.multiplicity < 2:
        raise ValueError("the soup needs at least 2 molecules: raise multiplicity or add clauses")
    size = len(start) * p.multiplicity
    pop = [start[int(i)] for i in rng.permutation(size) % len(start)]
    seen = dict.fromkeys(start)
    producer: dict[str, tuple[str, str]] = {}
    tally = Tally()
    period = round(1 / p.inflow_rate) if p.inflow_rate > 0 else 0
    collisions = productive = inflows = 0
    proved_at = 0 if target in seen else None

    def frame() -> Frame:
        return Frame(t=float(collisions), state={c: float(n) for c, n in Counter(pop).items()},
                     fired=tally.flush())

    def running() -> bool:
        return (not p.max_collisions or collisions < p.max_collisions) and proved_at is None

    yield frame()
    while running():
        i = int(rng.integers(size))
        j = int(rng.integers(size))
        while j == i:
            j = int(rng.integers(size))
        # Algorithm 3.1 draws k as a boolean, algorithm 3.2 as a position.
        k = bool(rng.integers(2)) if p.replacement == "educt" else int(rng.integers(size))
        a, b = pop[i], pop[j]
        options = chem.outcomes(a, b) if chem.allowed(a, b, support) else []
        collisions += 1
        if options:
            productive += 1
            r = options[int(rng.integers(len(options)))]
            if p.replacement == "educt":
                lhs, rhs = (a, b), ((b, r) if k else (a, r))
                pop[i if k else j] = r
            else:
                lost = {i: a, j: b}.get(k)
                lhs = (a, b)
                rhs = ((b if k == i else a), r) if lost is not None else (a, b, r)
                pop[k] = r
            if r not in seen:
                seen[r] = None
                producer[r] = (a, b) if a <= b else (b, a)
            if support[a] or support[b]:
                support[r] = True
            support.setdefault(r, False)
            tally.add(lhs, rhs)
            if r == target:
                proved_at = collisions
        if (p.elastic_inflow and not options) or (period and collisions % period == 0):
            inflows += 1
            pop[int(rng.integers(size))] = start[int(rng.integers(len(start)))]
        if collisions % size == 0 and running():
            yield frame()
    if collisions:
        yield frame()

    reactions = [Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in tally.reactions()]
    final = Counter(pop)
    return Network(
        species=_species(chem, seen),
        reactions=reactions,
        status="observed",
        initial_state={c: float(p.multiplicity) for c in start},
        extras={
            "axioms": axioms,
            "goal": goal,
            "final_state": {c: n for c, n in sorted(final.items())},
            "analysis": _analysis(chem, target, producer, start, {
                "collisions": collisions,
                "collisions_to_proof": proved_at,
                "productive_collisions": productive,
                "inflows": inflows,
                "generations": collisions / size,
            }),
        },
    )
