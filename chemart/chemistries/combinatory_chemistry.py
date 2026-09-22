"""Combinatory Chemistry: Kruszewski & Mikolov (ALIFE 2020).

Catalog id: combinatory-chemistry.

Molecules are expressions of combinatory logic over S, K and I, written with
as few parentheses as possible (application is left-associative, so ``SII`` is
``((SI)I)``). The reactor starts from free atoms only and conserves every
atom: the logic's reduction rules are rewritten as reactions that release
what plain combinatory logic would destroy and take from the multiset what it
would copy (paper eqs. 1-3):

    α(I f)β          ->  α f β          + I
    α(K f g)β        ->  α f β          + g + K
    α(S f g x)β + x  ->  α(f x (g x))β  + S

An S-redex counts only when a copy of its third argument x, the reactant, is
in the multiset. Irreducible expressions are cleaved at their top split,
xy -> x + y, or condensed with the previous irreducible expression
(eq. 4). Reductions take precedence, which is why the paper calls them
auto-catalysed.

Algorithm 1: sample an expression e with probability proportional to its
count. If e has a redex, apply one drawn uniformly from the first
``max_reductions`` in outer-to-inner order. Otherwise flip a coin: cleave e,
or condense it with the remembered expression e_LEFT (or remember e if there is
none). ``F`` > 1 turns on reactant assemblage (Algorithm 2): a missing reactant
of at most F atoms is built on the spot from free atoms.

The network is what one run did: the distinct reactions that fired, with
counts, and their type in ``extras["reaction_kinds"]``.
"""

from __future__ import annotations

from collections import Counter

from chemart.network import Network, Reaction, Species

ATOMS = "SKI"
ARITY = {"I": 1, "K": 2, "S": 3}

# A molecule is a tuple of items whose first item is an atom; an item is an
# atom (one-letter str) or a sub-expression (a tuple of at least two items).
# ("S",) is the atom S as a molecule. Lowercase letters are inert variables,
# used to state the rules and in tests.


# --- expressions -------------------------------------------------------------------
def _seq(items) -> tuple:
    """Normalise: a leading sub-expression is spliced in, ((SI)I) is SII."""
    items = tuple(items)
    while items and isinstance(items[0], tuple):
        items = items[0] + items[1:]
    return items


def _item(seq):
    """A molecule as an item of a larger expression."""
    return seq[0] if len(seq) == 1 else seq


def parse(text: str, variables: bool = False) -> tuple:
    """Parse an expression; parentheses are simplified (``(S(I)((K)K))`` is ``SI(KK)``)."""
    stack: list[list] = [[]]
    for ch in text:
        if ch.isspace():
            continue
        if ch == "(":
            stack.append([])
        elif ch == ")":
            if len(stack) == 1:
                raise ValueError(f"unbalanced ')' in expression {text!r}")
            group = stack.pop()
            if not group:
                raise ValueError(f"empty parentheses in expression {text!r}")
            stack[-1].append(_item(_seq(group)))
        elif ch in ARITY or (variables and ch.isalpha() and ch.islower()):
            stack[-1].append(ch)
        else:
            raise ValueError(f"unknown symbol {ch!r} in expression {text!r}; atoms are {ATOMS}")
    if len(stack) != 1:
        raise ValueError(f"unbalanced '(' in expression {text!r}")
    seq = _seq(stack[0])
    if not seq:
        raise ValueError(f"empty expression {text!r}")
    return seq


def show(seq) -> str:
    """Canonical string with as few parentheses as possible."""
    if isinstance(seq, str):
        return seq
    return "".join(x if isinstance(x, str) else f"({show(x)})" for x in seq)


def bracketed(seq) -> str:
    """Fully bracketed binary application tree: ``SII`` -> ``((SI)I)``."""
    def item(x):
        return x if isinstance(x, str) else bracketed(x)
    out = item(seq[0])
    for x in seq[1:]:
        out = f"({out}{item(x)})"
    return out


def size(seq) -> int:
    """Number of atoms."""
    if isinstance(seq, str):
        return 1
    return sum(size(x) for x in seq)


def atom_counts(text: str) -> Counter:
    return Counter(ch for ch in text if ch in ARITY)


def apply(a: tuple, b: tuple) -> tuple:
    """Condensation: a applied to b."""
    return a + (_item(b),)


def cleave(seq: tuple) -> tuple[tuple, tuple] | None:
    """The top split: e = (x y) -> x, y. None for an atom."""
    if len(seq) == 1:
        return None
    last = seq[-1]
    return seq[:-1], (last,) if isinstance(last, str) else last


# --- reductions ----------------------------------------------------------------------
def redexes(seq, path=()) -> list[tuple[tuple, str, str | None]]:
    """Every redex of plain combinatory logic, outer to inner (pre-order).

    Returns (path, combinator, reactant): path indexes nested sub-expressions,
    reactant is the S-redex's third argument as a molecule string, else None.
    The order is that of the authors' ``all-reductions`` (a redex before those
    inside its arguments, arguments left to right).
    """
    out = []
    head = seq[0]
    if head in ARITY and len(seq) > ARITY[head]:
        out.append((path, head, show(_mol(seq[3])) if head == "S" else None))
    for i, x in enumerate(seq[1:], start=1):
        if isinstance(x, tuple):
            out.extend(redexes(x, path + (i,)))
    return out


def _mol(item) -> tuple:
    return (item,) if isinstance(item, str) else item


def _rewrite(seq) -> tuple[tuple, list[tuple]]:
    """Reduce the redex at the head of seq: (result, by-products)."""
    head = seq[0]
    if head == "I":
        return _seq((seq[1],) + seq[2:]), [("I",)]
    if head == "K":
        return _seq((seq[1],) + seq[3:]), [_mol(seq[2]), ("K",)]
    f, g, x = seq[1], seq[2], seq[3]
    gx = _item(_seq((g, x)))
    return _seq((f, x, gx) + seq[4:]), [("S",)]


def reduce_at(seq: tuple, path: tuple) -> tuple[tuple, list[tuple]]:
    """Reduce the redex at `path`: (new expression, by-products)."""
    if not path:
        return _rewrite(seq)
    i = path[0]
    inner, by = reduce_at(seq[i], path[1:])
    return seq[:i] + (_item(inner),) + seq[i + 1:], by


# --- the reactor ------------------------------------------------------------------
class Pool:
    """The multiset P, with sampling proportional to counts."""

    def __init__(self, molecules):
        self.items: list[str] = []
        self.where: dict[str, set[int]] = {}
        self.count: Counter = Counter()
        for m in molecules:
            self.add(m)

    def add(self, m: str) -> None:
        self.where.setdefault(m, set()).add(len(self.items))
        self.items.append(m)
        self.count[m] += 1

    def remove(self, m: str) -> None:
        i = self.where[m].pop()
        last = self.items.pop()
        if i < len(self.items):
            self.items[i] = last
            slots = self.where[last]
            slots.remove(len(self.items))
            slots.add(i)
        self.count[m] -= 1
        if not self.count[m]:
            del self.count[m]
            del self.where[m]

    def sample(self, rng) -> str:
        return self.items[int(rng.integers(len(self.items)))]


class Reactor:
    def __init__(self, pool: Pool, F: int, max_reductions: int):
        self.pool = pool
        self.F = F
        self.max_reductions = max_reductions
        self.left: str | None = None
        self._parsed: dict[str, tuple] = {}
        self._redexes: dict[str, list] = {}

    def parsed(self, m: str) -> tuple:
        seq = self._parsed.get(m)
        if seq is None:
            seq = self._parsed[m] = parse(m)
        return seq

    def available(self, x: str) -> bool:
        """Is the reactant x in P, or can reactant assemblage build it?"""
        if self.pool.count[x]:
            return True
        if len(x) > self.F:       # atoms <= characters, so this is a cheap pre-check
            if size(self.parsed(x)) > self.F:
                return False
        need = atom_counts(x)
        return all(self.pool.count[a] >= n for a, n in need.items())

    def reductions(self, m: str) -> list:
        """The CC-redexes of m (S-redexes need their reactant), at most max_reductions."""
        found = self._redexes.get(m)
        if found is None:
            found = self._redexes[m] = redexes(self.parsed(m))
        out = []
        for r in found:
            if r[1] != "S" or self.available(r[2]):
                out.append(r)
                if len(out) == self.max_reductions:
                    break
        return out

    def step(self, rng):
        """One iteration of Algorithm 1; returns a list of (kind, lhs, rhs) events."""
        return self.react(self.pool.sample(rng), rng)

    def react(self, e: str, rng):
        """Algorithm 1 applied to the expression e (which must be in the pool)."""
        pool = self.pool
        options = self.reductions(e)
        if options:
            path, kind, x = options[int(rng.integers(len(options)))]
            events = []
            if kind == "S" and not pool.count[x]:            # reactant assemblage
                atoms = sorted(atom_counts(x).elements())
                for a in atoms:
                    pool.remove(a)
                pool.add(x)
                events.append(("assemblage", tuple(atoms), (x,)))
            new, by = reduce_at(self.parsed(e), path)
            products = (show(new), *(show(b) for b in by))
            lhs = (e, x) if kind == "S" else (e,)
            for m in lhs:
                pool.remove(m)
            for m in products:
                pool.add(m)
            events.append((kind, lhs, products))
            return events
        if rng.random() < 0.5:
            parts = cleave(self.parsed(e))
            if parts is None:
                return []
            x, y = show(parts[0]), show(parts[1])
            pool.remove(e)
            pool.add(x)
            pool.add(y)
            return [("cleave", (e,), (x, y))]
        left = self.left
        if left is None or pool.count[left] < (2 if left == e else 1):
            self.left = e
            return []
        product = show(apply(self.parsed(left), self.parsed(e)))
        pool.remove(left)
        pool.remove(e)
        pool.add(product)
        self.left = None
        return [("condense", (left, e), (product,))]


def generate(p, rng):
    initial = {a: n for a, n in (("I", p.n_I), ("K", p.n_K), ("S", p.n_S)) if n}
    if sum(initial.values()) < 2:
        raise ValueError("n_I + n_K + n_S must be at least 2")
    pool = Pool([a for a, n in initial.items() for _ in range(n)])
    reactor = Reactor(pool, p.F, p.max_reductions)

    fired: dict[tuple, list] = {}
    kinds: dict[tuple, str] = {}
    consumed: Counter = Counter()           # reactants of S-reactions, per window
    window = Counter()                      # event kinds in the current window
    analysis = {"iteration": [], "diversity": [], "mean_length": [], "reductions": [],
                "free_atoms": {a: [] for a in ATOMS}, "top_reactants": []}

    def measure(t):
        n = len(pool.items)
        analysis["iteration"].append(t)
        analysis["diversity"].append(len(pool.count))
        analysis["mean_length"].append(round(sum(size(reactor.parsed(m)) * c for m, c in pool.count.items()) / n, 3))
        steps = sum(window.values())
        reduced = window["I"] + window["K"] + window["S"]
        analysis["reductions"].append(round(reduced / steps, 4) if steps else 0.0)
        for a in ATOMS:
            analysis["free_atoms"][a].append(pool.count[a])
        analysis["top_reactants"].append(dict(consumed.most_common(5)))
        consumed.clear()
        window.clear()

    measure(0)
    for t in range(1, p.iterations + 1):
        events = reactor.step(rng)
        window["idle" if not events else events[-1][0]] += 1
        for kind, lhs, rhs in events:
            key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
            if key not in fired:
                fired[key] = [lhs, rhs, 0]
                kinds[key] = kind
            fired[key][2] += 1
            if kind == "S":
                consumed[lhs[1]] += 1
        if t % p.record_every == 0 or t == p.iterations:
            measure(t)

    names = dict.fromkeys(initial)
    for lhs, rhs, _ in fired.values():
        names.update(dict.fromkeys((*lhs, *rhs)))
    species = [Species(m, structure=bracketed(reactor.parsed(m))) for m in names]
    conservation = [
        {"name": f"atom {a}", "vector": {s.id: atom_counts(s.id)[a] for s in species}}
        for a in ATOMS
    ]
    return Network(
        species=species,
        reactions=[Reaction.of(lhs, rhs, count=c) for lhs, rhs, c in fired.values()],
        status="observed",
        initial_state={a: float(n) for a, n in initial.items()},
        extras={
            "reaction_kinds": [kinds[k] for k in fired],
            "conservation": conservation,
            "final_state": dict(pool.count.most_common()),
            "analysis": analysis,
        },
    )
