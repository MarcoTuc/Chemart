"""AlChemy: Fontana's lambda-calculus chemistry (book 9.1). Catalog id: alchemy.

Molecules are closed normal forms of the untyped lambda-calculus, identified
modulo renaming of bound variables. A collision of an ordered pair applies
the first object to the second and reduces the result to normal form:

    s1 + s2 -> s1 + s2 + NF((s1) s2)          (book eq. 9.4)

Reduction is "pragmatic" (Fontana & Buss 1994, 4.4 and 5.3): it is cut off
after a number of reduction steps or when the term grows beyond a size limit,
and then the collision is elastic. Products larger than a normal-form size
limit, and products rejected by the boundary conditions (no copy actions,
forbidden syntactic patterns), are elastic too.

Terms are stored with de Bruijn indices, so alpha-equivalent terms are equal
and the species id is the canonical de Bruijn string: ``^`` is an abstraction
(its body extends as far right as possible), ``(M)N`` is the application of M
to N, and a number is a variable (1 = the nearest enclosing binder). The
species structure is the paper's standardized form, with bound variables
named x1, x2, ... in the order their binders occur, e.g. ``λx1.λx2.(x1)x2``.

method "soup" runs the paper's stochastic flow reactor (chemart.soup with
constant dilution) and returns the observed reactions; method "closure"
returns the reaction closure of a seed set (chemart.expand), complete or cut
off by max_species.
"""

from __future__ import annotations

import re
import sys
from collections import Counter

from chemart.expand import expand
from chemart.helpers.params import apportion
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import soup

# --- terms --------------------------------------------------------------------
# A term is a tuple (tag, a, b, size, fv):
#   VAR: a = de Bruijn index (0 = nearest binder)
#   LAM: a = body
#   APP: a = operator, b = operand
# size is the length in characters of the paper's notation, counting every
# variable name as 2 characters ("x1"), "λx1." as 4 and the parentheses of an
# application as 2; fv is 1 + the largest free index (0 for a closed term).
VAR, LAM, APP = 0, 1, 2


def var(i: int) -> tuple:
    return (VAR, i, None, 2, i + 1)


def lam(body: tuple) -> tuple:
    return (LAM, body, None, body[3] + 4, max(body[4] - 1, 0))


def app(f: tuple, x: tuple) -> tuple:
    return (APP, f, x, f[3] + x[3] + 2, max(f[4], x[4]))


def size(t: tuple) -> int:
    return t[3]


def is_closed(t: tuple) -> bool:
    return t[4] == 0


def _shift(t: tuple, d: int, cutoff: int) -> tuple:
    """Add d to every free index >= cutoff."""
    if t[4] <= cutoff:
        return t
    tag = t[0]
    if tag == VAR:
        return var(t[1] + d)
    if tag == LAM:
        return lam(_shift(t[1], d, cutoff + 1))
    return app(_shift(t[1], d, cutoff), _shift(t[2], d, cutoff))


def _subst(t: tuple, s: tuple, depth: int, shifted: dict) -> tuple:
    """t[depth := s] for the binder being removed; indices above depth drop by one."""
    if t[4] <= depth:
        return t
    tag = t[0]
    if tag == VAR:
        i = t[1]
        if i == depth:
            if depth not in shifted:
                shifted[depth] = _shift(s, depth, 0)
            return shifted[depth]
        return var(i - 1)
    if tag == LAM:
        return lam(_subst(t[1], s, depth + 1, shifted))
    return app(_subst(t[1], s, depth, shifted), _subst(t[2], s, depth, shifted))


def beta(abstraction: tuple, argument: tuple) -> tuple:
    """Contract the redex (abstraction) argument, capture-free."""
    return _subst(abstraction[1], argument, 0, {})


class Diverged(Exception):
    """Pragmatic reduction gave up: step or size budget exceeded."""


def normal_form(t: tuple, max_steps: int = 10000, max_size: int = 4000) -> tuple[tuple, int]:
    """Normal-order (leftmost-outermost) reduction to normal form.

    Returns (normal form, number of beta steps). Raises Diverged when more
    than max_steps contractions are needed or when the whole term exceeds
    max_size characters after a contraction.
    """
    steps = 0

    def tick(total: int) -> None:
        nonlocal steps
        steps += 1
        if steps > max_steps:
            raise Diverged(f"no normal form within {max_steps} steps")
        if total > max_size:
            raise Diverged(f"term grew beyond {max_size} characters")

    # ctx is the number of characters of the whole term outside the subterm.
    def whnf(t, ctx):
        while t[0] == APP:
            x = t[2]
            f = whnf(t[1], ctx + 2 + x[3])
            if f[0] != LAM:
                return app(f, x) if f is not t[1] else t
            t = beta(f, x)
            tick(ctx + t[3])
        return t

    def nf(t, ctx):
        tag = t[0]
        if tag == VAR:
            return t
        if tag == LAM:
            body = nf(t[1], ctx + 4)
            return t if body is t[1] else lam(body)
        h = whnf(t, ctx)
        if h[0] != APP:
            return nf(h, ctx)
        f = nf(h[1], ctx + 2 + h[2][3])
        x = nf(h[2], ctx + 2 + f[3])
        return app(f, x)

    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(limit, 4 * max_size + 1000))
    try:
        return nf(t, 0), steps
    finally:
        sys.setrecursionlimit(limit)


# --- notation -----------------------------------------------------------------
def to_id(t: tuple) -> str:
    """Canonical de Bruijn string: ^ abstraction, (M)N application, 1-based indices."""
    out: list[str] = []

    def walk(t):
        while True:
            tag = t[0]
            if tag == VAR:
                out.append(str(t[1] + 1))
                return
            if tag == LAM:
                out.append("^")
                t = t[1]
                continue
            out.append("(")
            walk(t[1])
            out.append(")")
            t = t[2]

    walk(t)
    return "".join(out)


def to_text(t: tuple) -> str:
    """The paper's standardized notation: binders named x1, x2, ... in order of occurrence."""
    out: list[str] = []
    names: list[str] = []
    counter = 0

    def walk(t):
        nonlocal counter
        pushed = 0
        while True:
            tag = t[0]
            if tag == VAR:
                i = t[1]
                out.append(names[-1 - i] if i < len(names) else f"free{i - len(names) + 1}")
                break
            if tag == LAM:
                counter += 1
                names.append(f"x{counter}")
                pushed += 1
                out.append(f"λ{names[-1]}.")
                t = t[1]
                continue
            out.append("(")
            walk(t[1])
            out.append(")")
            t = t[2]
        del names[len(names) - pushed:]

    walk(t)
    return "".join(out)


_TOKEN = re.compile(r"\s*(?:(?P<lam>[λ\\^])|(?P<name>[A-Za-z_][A-Za-z0-9_']*)|(?P<num>\d+)|(?P<sym>[().]))")


def parse(text: str) -> tuple:
    """Parse the paper's notation or a species id into a term.

    Grammar: term := λx.term | (term)term | x, where λ may be written λ or \\,
    and application is always written (M)N as in Fontana & Buss. De Bruijn
    ids are accepted too: ^term for an abstraction and 1-based indices for
    variables. Raises ValueError for syntax errors and free variables.
    """
    if not isinstance(text, str):
        raise ValueError(f"a lambda term must be a string, got {text!r}")
    tokens: list[tuple[str, str]] = []
    pos = 0
    stripped = text.rstrip()
    while pos < len(stripped):
        m = _TOKEN.match(stripped, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot parse lambda term {text!r} at position {pos}")
        kind = m.lastgroup
        tokens.append((kind, m.group(kind)))
        pos = m.end()
    i = 0
    binders: list[str | None] = []

    def expect(sym):
        nonlocal i
        if i >= len(tokens) or tokens[i] != ("sym", sym):
            raise ValueError(f"cannot parse lambda term {text!r}: expected {sym!r}")
        i += 1

    def term():
        nonlocal i
        if i >= len(tokens):
            raise ValueError(f"cannot parse lambda term {text!r}: unexpected end")
        kind, value = tokens[i]
        if kind == "lam":
            i += 1
            name = None
            if value != "^":
                if i >= len(tokens) or tokens[i][0] != "name":
                    raise ValueError(f"cannot parse lambda term {text!r}: λ needs a variable name")
                name = tokens[i][1]
                i += 1
                expect(".")
            binders.append(name)
            body = term()
            binders.pop()
            return lam(body)
        if kind == "sym" and value == "(":
            i += 1
            f = term()
            expect(")")
            return app(f, term())
        if kind == "name":
            i += 1
            for depth, b in enumerate(reversed(binders)):
                if b == value:
                    return var(depth)
            raise ValueError(f"lambda term {text!r} has the free variable {value!r}; "
                             "molecules must be closed (bind it with a leading λ)")
        if kind == "num":
            i += 1
            k = int(value)
            if not 1 <= k <= len(binders):
                raise ValueError(f"lambda term {text!r}: de Bruijn index {k} has no binder")
            return var(k - 1)
        raise ValueError(f"cannot parse lambda term {text!r}: unexpected {value!r}")

    t = term()
    if i != len(tokens):
        raise ValueError(f"cannot parse lambda term {text!r}: trailing {tokens[i][1]!r}")
    return t


# --- random terms (paper 5.3, step 2) ----------------------------------------
def random_term(rng, p_variable: float, p_abstraction: float, max_depth: int,
                p_bound: float, n_free: int) -> tuple:
    """A random closed term, before reduction.

    Recursively: a variable with probability p_variable, an abstraction with
    p_abstraction, an application otherwise; at max_depth a variable. A
    variable is one of the enclosing binders (uniformly) with probability
    p_bound, else one of n_free free names. Free names are then bound by
    abstractions prefixed to the term (standardization).
    """
    free_order: list[int] = []

    def gen(depth, nbinders):
        r = float(rng.random())
        if depth >= max_depth or r < p_variable:
            if nbinders and float(rng.random()) < p_bound:
                return ("b", int(rng.integers(nbinders)))
            k = int(rng.integers(n_free))
            if k not in free_order:
                free_order.append(k)
            return ("f", k)
        if r < p_variable + p_abstraction:
            return ("l", gen(depth + 1, nbinders + 1))
        return ("a", gen(depth + 1, nbinders), gen(depth + 1, nbinders))

    raw = gen(0, 0)
    # the first free name gets the outermost binder
    order = {k: j for j, k in enumerate(free_order)}
    n = len(free_order)

    def build(node, nbinders):
        kind = node[0]
        if kind == "b":
            return var(node[1])
        if kind == "f":
            return var(nbinders + n - 1 - order[node[1]])
        if kind == "l":
            return lam(build(node[1], nbinders + 1))
        return app(build(node[1], nbinders), build(node[2], nbinders))

    t = build(raw, 0)
    for _ in range(n):
        t = lam(t)
    return t


# --- the reaction -------------------------------------------------------------
class Chemistry:
    """Collision rule plus boundary conditions, with a cache of outcomes by species id."""

    def __init__(self, p):
        self.max_steps = p.max_steps
        self.max_size = p.max_size
        self.max_nf_size = p.max_nf_size
        self.no_copy = p.filter == "no-copy"
        self.patterns = []
        for pat in p.forbidden_patterns:
            if not isinstance(pat, str):
                raise ValueError(f"forbidden_patterns must be regular expressions (strings), got {pat!r}")
            try:
                self.patterns.append(re.compile(pat))
            except re.error as err:
                raise ValueError(f"forbidden_patterns: {pat!r} is not a valid regular expression ({err})") from None
        self.mediator = None
        if p.mediator:
            self.mediator = parse(p.mediator)
        self.terms: dict[str, tuple] = {}
        self.cache: dict[tuple[str, str], str | None] = {}
        self.stats = Counter()

    def add(self, t: tuple) -> str:
        sid = to_id(t)
        self.terms.setdefault(sid, t)
        return sid

    def forbidden(self, sid: str) -> bool:
        return any(p.search(sid) for p in self.patterns)

    def reduce(self, t: tuple) -> tuple | None:
        """Pragmatic reduction: the normal form, or None if a budget is exceeded."""
        try:
            nf, _ = normal_form(t, self.max_steps, self.max_size)
        except Diverged:
            return None
        return nf if nf[3] <= self.max_nf_size else None

    def product(self, a: str, b: str) -> str | None:
        """Species id of a o b, or None for an elastic collision."""
        key = (a, b)
        if key in self.cache:
            return self.cache[key]
        ta, tb = self.terms[a], self.terms[b]
        expr = app(ta, tb) if self.mediator is None else app(app(self.mediator, ta), tb)
        nf = self.reduce(expr)
        out = None
        if nf is None:
            self.stats["no_normal_form"] += 1
        else:
            sid = to_id(nf)
            if self.no_copy and sid in (a, b):
                self.stats["copy"] += 1
            elif self.forbidden(sid):
                self.stats["forbidden"] += 1
            else:
                self.terms.setdefault(sid, nf)
                out = sid
        self.cache[key] = out
        return out

    def react(self, a: str, b: str):
        c = self.product(a, b)
        return None if c is None else (a, b, c)

    def rate(self, lhs, rhs) -> dict:
        """Unit rate per ordered collision (paper eq. 18): k counts the orders giving the product."""
        a, b = lhs
        (c,) = (Counter(rhs) - Counter(lhs)).elements()
        return {"law": "mass-action", "k": float(sum(1 for x, y in {(a, b), (b, a)} if self.product(x, y) == c))}

    def reaction(self, lhs, rhs, count=None) -> Reaction:
        return Reaction.of(lhs, rhs, rate=self.rate(lhs, rhs), count=count)

    def species(self, ids) -> list[Species]:
        return [Species(s, structure=to_text(self.terms[s])) for s in ids]


def _seed_terms(p, chem: Chemistry, rng) -> list[str]:
    if p.terms:
        if not isinstance(p.terms, list):
            raise ValueError(f"terms must be a list of lambda terms, got {p.terms!r}")
        out = []
        for text in p.terms:
            t = parse(text)
            nf = chem.reduce(t)
            if nf is None:
                raise ValueError(f"terms: {text!r} has no normal form within the reduction budget")
            out.append(chem.add(nf))
        return list(dict.fromkeys(out))
    if p.p_variable + p.p_abstraction > 1:
        raise ValueError(f"p_variable + p_abstraction must be <= 1 (the rest is p_application), "
                         f"got {p.p_variable} + {p.p_abstraction}")
    out: dict[str, None] = {}
    attempts = 0
    while len(out) < p.M:
        attempts += 1
        if attempts > 50 * p.M:
            raise ValueError(f"could only draw {len(out)} distinct random normal forms for M={p.M}; "
                             "increase max_depth or p_application, or lower M")
        t = random_term(rng, p.p_variable, p.p_abstraction, p.max_depth, p.p_bound, p.n_free)
        nf = chem.reduce(t)
        if nf is None:
            continue
        sid = to_id(nf)
        if sid in out or chem.forbidden(sid):
            continue
        chem.add(nf)
        out[sid] = None
    return list(out)


def generate(p, rng):
    chem = Chemistry(p)
    seeds = _seed_terms(p, chem, rng)
    extras = {
        "notation": "species id: de Bruijn string (^ abstraction, (M)N application, 1-based indices); "
                    "structure: Fontana-Buss standardized form",
    }
    if p.method == "closure":
        found, pairs, status = expand(chem.react, seeds, arity=2, max_species=p.max_species, ordered=True)
        reactions = [chem.reaction(lhs, rhs) for lhs, rhs in pairs]
        extras["seed"] = seeds
        extras["analysis"] = {
            "elastic": dict(chem.stats),
            "self_maintaining": _self_maintaining(chem, found),
        }
        return Network(species=chem.species(found), reactions=reactions, status=status,
                       outflow=CONSTANT_TOTAL, extras=extras)
    return _soup(p, chem, seeds, rng, extras)


def _self_maintaining(chem: Chemistry, ids) -> bool:
    """Paper eq. 38: every object is produced by a collision within the set."""
    known = set(ids)
    made = {chem.product(a, b) for a in ids for b in ids} & known
    return made == known


def _soup(p, chem, seeds, rng, extras):
    if p.terms:
        if p.M < len(seeds):
            raise ValueError(f"M={p.M} is smaller than the {len(seeds)} distinct terms given")
        counts = apportion(p.M, [1.0] * len(seeds))
        pop = [s for s, m in zip(seeds, counts) for _ in range(m)]
    else:
        pop = list(seeds)
    start = Counter(pop)
    size_ = len(pop)
    fired: dict[tuple, list] = {}
    diversity = [len(start)]
    done = 0
    while done < p.collisions:
        n = min(size_, p.collisions - done)
        chunk, pop = soup(chem.react, pop, n, rng, arity=2, dilution="constant")
        done += n
        for lhs, rhs, count in chunk:
            key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
            fired.setdefault(key, [lhs, rhs, 0])[2] += count
        diversity.append(len(set(pop)))
    reactions = [chem.reaction(lhs, rhs, count) for lhs, rhs, count in fired.values()]
    ids = list(dict.fromkeys([*start, *(c for _, rhs, _ in fired.values() for c in rhs)]))
    final = Counter(pop)
    extras["analysis"] = {
        "collisions_per_sample": size_,
        "distinct_species": diversity,
        "elastic_pairs": dict(chem.stats),
    }
    extras["final_state"] = {s: n for s, n in final.most_common()}
    return Network(species=chem.species(ids), reactions=reactions, status="observed",
                   initial_state={s: n for s, n in start.items()}, outflow=CONSTANT_TOTAL,
                   extras=extras)
