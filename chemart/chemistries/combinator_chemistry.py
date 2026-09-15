"""Combinator chemistry: Speroni di Fenizio (2000; thesis 2007, ch. 6-7).

Catalog id: combinator-chemistry.

Molecules are combinators in normal form: strings of atoms (B, C, I, K, R, S,
W) and balanced parentheses, written with as few parentheses as possible
(``SI(KK)``). A collision of a with b applies a to b, giving ``a(b)``, and
reduces it to normal form with the atom rules of thesis table 6.1:

    B x y z -> x(yz)     C x y z -> x z y     I x -> x
    K x y   -> x         S x y z -> x z(yz)   W x y -> x y y
    R x y   -> x, and y is released as a separate molecule

(in the 2000 paper K itself releases y: ``k_action="release"``). Only the
first atom of a (sub-)combinator acts, when it is followed by at least as many
items as its arity. The order of reduction matters once sizes, atoms and
releases are bounded, so it follows thesis table 6.3: all K redexes first,
then all B, C, R and I redexes, and only then the outermost S or W redex.

A reaction is elastic when the normal forms are not reached within
``max_reductions`` atom reductions, when a combinator grows beyond
``max_react_size`` atoms during the reduction, or when a product has more than
``max_size`` atoms or ``max_depth`` nested parentheses.

reaction "reactive" (thesis level 3, 2000 paper): a + b -> c1 + ... + cn, the
reactants are used up, and the atoms come from and go back to a pool of free
atoms (species ``free:X``), so every reaction conserves atoms. reaction
"catalytic" (thesis level 1, AlChemy): a + b -> a + b + c1 + ... + cn with an
unlimited supply of atoms and a dilution that keeps the population constant.

method "closure" returns every reaction reachable from the seed molecules
(chemart.expand.expand); method "soup" runs the well-stirred flow reactor and
returns the reactions that fired.
"""

from __future__ import annotations

from collections import Counter

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import soup

ATOMS = "BCIKRSW"
ARITY = {"B": 3, "C": 3, "I": 1, "K": 2, "R": 2, "S": 3, "W": 2}
FREE = "free:"
MAX_TRY = 10000   # ECAL 2001 MaxTry: attempts to assemble one random molecule

# A term is a tuple of items whose first item is an atom (a one-letter str);
# an item is an atom or a sub-term with at least two items. Lowercase letters
# are inert variables (no rule), used to state the rules and in tests.


class Elastic(Exception):
    """The reduction was cut off by one of the limits."""


# --- terms ---------------------------------------------------------------------
def _seq(items) -> tuple:
    """Normalise a sequence: a leading parenthesis is dropped, (x) is x."""
    items = tuple(items)
    while items and isinstance(items[0], tuple):
        items = items[0] + items[1:]
    return items


def _item(items):
    seq = _seq(items)
    return seq[0] if len(seq) == 1 else seq


def parse(text: str, variables: bool = False) -> tuple:
    """Parse a combinator; parentheses are simplified (``(S(I)((K)K))`` is ``SI(KK)``)."""
    stack: list[list] = [[]]
    for ch in text:
        if ch.isspace():
            continue
        if ch == "(":
            stack.append([])
        elif ch == ")":
            if len(stack) == 1:
                raise ValueError(f"unbalanced ')' in combinator {text!r}")
            group = stack.pop()
            if group:
                stack[-1].append(_item(group))
        elif ch in ARITY or (variables and ch.isalpha() and ch.islower()):
            stack[-1].append(ch)
        else:
            raise ValueError(f"unknown symbol {ch!r} in combinator {text!r}; atoms are {ATOMS}")
    if len(stack) != 1:
        raise ValueError(f"unbalanced '(' in combinator {text!r}")
    term = _seq(stack[0])
    if not term:
        raise ValueError(f"empty combinator {text!r}")
    return term


def show(term) -> str:
    """Canonical string with as few parentheses as possible."""
    if isinstance(term, str):
        return term
    return "".join(x if isinstance(x, str) else f"({show(x)})" for x in term)


def bracketed(term) -> str:
    """Fully bracketed binary application tree: ``SI(KK)`` -> ``((SI)(KK))``."""
    def item(x):
        return x if isinstance(x, str) else bracketed(x)
    out = item(term[0])
    for x in term[1:]:
        out = f"({out}{item(x)})"
    return out


def size(term) -> int:
    return sum(1 if isinstance(x, str) else size(x) for x in term)


def depth(term) -> int:
    return max((1 + depth(x) for x in term if isinstance(x, tuple)), default=0)


def atom_counts(text: str) -> Counter:
    return Counter(ch for ch in text if ch in ARITY)


def apply(a: tuple, b: tuple) -> tuple:
    """a applied to b: b becomes a sub-combinator of a."""
    return a + (_item(b),)


# --- reduction ------------------------------------------------------------------
def _redex(seq) -> bool:
    arity = ARITY.get(seq[0])
    return arity is not None and len(seq) - 1 >= arity


def _heads(seq, out: set) -> set:
    if _redex(seq):
        out.add(seq[0])
    for x in seq[1:]:
        if isinstance(x, tuple):
            _heads(x, out)
    return out


def _rule(seq, k_release: bool, released: list) -> tuple:
    head, args = seq[0], seq[1:]
    n = ARITY[head]
    x, rest = args[:n], args[n:]
    if head == "B":
        new = (x[0], _item((x[1], x[2])))
    elif head == "C":
        new = (x[0], x[2], x[1])
    elif head == "I":
        new = (x[0],)
    elif head == "K":
        new = (x[0],)
        if k_release:
            released.append(_seq((x[1],)))
    elif head == "R":
        new = (x[0],)
        released.append(_seq((x[1],)))
    elif head == "S":
        new = (x[0], x[2], _item((x[1], x[2])))
    else:  # W
        new = (x[0], x[1], x[1])
    return _seq(new + rest)


def _parallel_pass(seq, heads: frozenset, k_release: bool, released: list, count: list) -> tuple:
    """Reduce every redex whose head is in `heads`, innermost first."""
    items = [x if isinstance(x, str) else _item(_parallel_pass(x, heads, k_release, released, count))
             for x in seq[1:]]
    seq = (seq[0], *items)
    if seq[0] in heads and _redex(seq):
        count[0] += 1
        seq = _rule(seq, k_release, released)
    return seq


def _outermost(seq, heads: frozenset, path=(), level=0, best=None):
    """(level, path) of the outermost, then leftmost, redex with a head in `heads`."""
    if seq[0] in heads and _redex(seq):
        cand = (level, path)
        if best is None or cand < best:
            best = cand
    for i, x in enumerate(seq):
        if isinstance(x, tuple):
            best = _outermost(x, heads, path + (i,), level + 1, best)
    return best


def _replace(seq, path, k_release, released):
    if not path:
        return _rule(seq, k_release, released)
    i = path[0]
    items = list(seq)
    items[i] = _item(_replace(seq[i], path[1:], k_release, released))
    return _seq(items)


_SW = frozenset("SW")
_BCRI = frozenset("BCRI")
_K = frozenset("K")


def reduce_terms(terms, k_release: bool = False, max_reductions: int = 100,
                 max_react_size: int = 50):
    """Normal forms of `terms` and of every molecule released on the way (thesis table 6.3).

    Returns (normal forms, reductions used, atom counts held after each pass).
    Raises Elastic when a limit is exceeded. The reduction budget is shared by
    all the molecules of one reaction.
    """
    queue = list(terms)
    done: list[tuple] = []
    count = [0]
    snapshots = [sum((atom_counts(show(t)) for t in queue), Counter())]
    while queue:
        term = queue.pop(0)
        while True:
            heads = _heads(term, set())
            released: list[tuple] = []
            if "K" in heads:
                term = _parallel_pass(term, _K, k_release, released, count)
            elif heads & _BCRI:
                term = _parallel_pass(term, _BCRI, k_release, released, count)
            elif heads & _SW:
                _, path = _outermost(term, _SW)
                term = _replace(term, path, k_release, released)
                count[0] += 1
            else:
                done.append(term)
                break
            if count[0] > max_reductions:
                raise Elastic(f"no normal form within {max_reductions} reductions")
            if size(term) > max_react_size or any(size(r) > max_react_size for r in released):
                raise Elastic(f"combinator longer than {max_react_size} atoms during the reduction")
            queue.extend(released)
            snapshots.append(sum((atom_counts(show(t)) for t in [term, *queue, *done]), Counter()))
    return done, count[0], snapshots


def normal_form(text: str, k_release: bool = False, max_reductions: int = 100,
                max_react_size: int = 50) -> list[str] | None:
    """Normal forms of a combinator (first the main one, then released ones), or None."""
    try:
        done, _, _ = reduce_terms([parse(text, variables=True)], k_release, max_reductions, max_react_size)
    except Elastic:
        return None
    return [show(t) for t in done]


class Chemistry:
    """The reaction a + b of one parameter set, cached."""

    def __init__(self, k_release=False, max_reductions=100, max_react_size=50, max_size=20,
                 max_depth=5, filter_reproduction=False):
        self.k_release = k_release
        self.max_reductions = max_reductions
        self.max_react_size = max_react_size
        self.max_size = max_size
        self.max_depth = max_depth
        self.filter_reproduction = filter_reproduction
        self._cache: dict = {}

    def normalise(self, term: tuple):
        """(products, peak atoms held beyond the start) or None if elastic."""
        try:
            done, _, snaps = reduce_terms([term], self.k_release, self.max_reductions, self.max_react_size)
        except Elastic:
            return None
        if any(size(t) > self.max_size or depth(t) > self.max_depth for t in done):
            return None
        start = snaps[0]
        peak = Counter()
        for s in snaps:
            for x, n in s.items():
                peak[x] = max(peak[x], n - start[x])
        return tuple(show(t) for t in done), +peak

    def react(self, a: str, b: str):
        """(products, atoms needed from the pool at the peak) for a applied to b, or None."""
        key = (a, b)
        if key not in self._cache:
            out = None
            if size(parse(a)) + size(parse(b)) <= self.max_react_size:
                out = self.normalise(apply(parse(a), parse(b)))
                if out is not None and self.filter_reproduction and (a in out[0] or b in out[0]):
                    out = None
            self._cache[key] = out
        return self._cache[key]


# --- parameters -----------------------------------------------------------------
def _basis(text: str) -> str:
    if not text or any(ch not in ARITY for ch in text) or len(set(text)) != len(text):
        raise ValueError(f"atoms must be distinct letters from {ATOMS!r}, got {text!r}")
    return "".join(ch for ch in ATOMS if ch in text)


def _molecules(p, chem: Chemistry, basis: str) -> list[str]:
    out = []
    for m in p.molecules:
        if not isinstance(m, str):
            raise ValueError(f"molecules must be combinator strings such as 'SI(KK)', got {m!r}")
        term = parse(m)
        text = show(term)
        extra = set(atom_counts(text)) - set(basis)
        if extra:
            raise ValueError(f"molecule {m!r} uses atoms {sorted(extra)} outside atoms={basis!r}")
        if _heads(term, set()):
            nf = normal_form(text, chem.k_release, chem.max_reductions, chem.max_react_size)
            raise ValueError(f"molecule {m!r} is not in normal form (it reduces to {nf}); "
                             "molecules are stored in normal form")
        if size(term) > chem.max_size or depth(term) > chem.max_depth:
            raise ValueError(f"molecule {m!r} exceeds max_size={chem.max_size} or max_depth={chem.max_depth}")
        out.append(text)
    return out


def generate(p, rng):
    basis = _basis(p.atoms)
    if p.max_size > p.max_react_size:
        raise ValueError(f"max_size ({p.max_size}) cannot exceed max_react_size ({p.max_react_size})")
    chem = Chemistry(p.k_action == "release", p.max_reductions, p.max_react_size, p.max_size,
                     p.max_depth, p.filter_reproduction)
    molecules = _molecules(p, chem, basis)
    if p.method == "closure":
        return _closure(p, chem, basis, molecules or list(basis))
    if p.reaction == "catalytic":
        return _catalytic_soup(p, chem, basis, molecules, rng)
    return _reactive_soup(p, chem, basis, molecules, rng)


# --- network building -----------------------------------------------------------
def _free_balance(lhs, rhs) -> tuple[list[str], list[str]]:
    """Free atoms taken from (left) and returned to (right) the pool by a reaction."""
    before = sum((atom_counts(m) for m in lhs), Counter())
    after = sum((atom_counts(m) for m in rhs), Counter())
    taken = after - before
    returned = before - after
    return (list((FREE + x for x in sorted(taken.elements()))),
            list((FREE + x for x in sorted(returned.elements()))))


def _species(molecules, basis, reactive) -> list[Species]:
    out = [Species(m, structure=bracketed(parse(m))) for m in dict.fromkeys(molecules)]
    if reactive:
        out += [Species(FREE + x) for x in basis]
    return out


def _conservation(species: list[Species], basis: str) -> list[dict]:
    out = []
    for x in basis:
        vector = {}
        for s in species:
            vector[s.id] = 1 if s.id == FREE + x else (0 if s.id.startswith(FREE) else atom_counts(s.id)[x])
        out.append({"name": f"atom {x}", "vector": vector})
    return out


def _reaction(lhs, rhs, reactive, count=None) -> Reaction:
    if reactive:
        taken, returned = _free_balance(lhs, rhs)
        return Reaction.of([*lhs, *taken], [*rhs, *returned], count=count)
    return Reaction.of(lhs, rhs, count=count)


def _closure(p, chem: Chemistry, basis: str, seed: list[str]) -> Network:
    reactive = p.reaction == "reactive"

    def react(a, b):
        out = chem.react(a, b)
        if out is None:
            return None
        return out[0] if reactive else (a, b, *out[0])

    seed = list(dict.fromkeys(seed))
    found, reactions, status = expand(react, seed, arity=2, max_species=p.max_species, ordered=True)
    species = _species(found, basis, reactive)
    extras = {"seed": seed}
    if reactive:
        extras["conservation"] = _conservation(species, basis)
    return Network(
        species=species,
        reactions=[_reaction(lhs, rhs, reactive) for lhs, rhs in reactions],
        status=status,
        outflow=None if reactive else CONSTANT_TOTAL,
        extras=extras,
    )


# --- soups ----------------------------------------------------------------------
def _random_molecule(rng, chem: Chemistry, basis: str, pool: Counter | None):
    """ECAL 2001: each atom, '(' and ')' with equal probability; an unmatched ')' ends it.

    Returns the normal forms (a molecule may release others while it is
    normalised) or None after MAX_TRY failed attempts. pool None = unlimited.
    """
    symbols = basis + "()"
    for _ in range(MAX_TRY):
        stack: list[list] = [[]]
        used = Counter()
        ok = True
        while ok:
            ch = symbols[int(rng.integers(len(symbols)))]
            if ch == "(":
                stack.append([])
            elif ch == ")":
                if len(stack) == 1:
                    break
                group = stack.pop()
                if group:
                    stack[-1].append(_item(group))
            else:
                used[ch] += 1
                if (pool is not None and used[ch] > pool[ch]) or sum(used.values()) > chem.max_react_size:
                    ok = False
                stack[-1].append(ch)
        if not ok:
            continue
        while len(stack) > 1:
            group = stack.pop()
            if group:
                stack[-1].append(_item(group))
        term = _seq(stack[0])
        if not term:
            continue
        out = chem.normalise(term)
        if out is None:
            continue
        products, peak = out
        if pool is not None and any(used[x] + peak[x] > pool[x] for x in peak):
            continue
        return products
    return None


def _catalytic_soup(p, chem: Chemistry, basis: str, molecules: list[str], rng) -> Network:
    start: list[str] = list(molecules)
    if not start:
        while len(start) < p.M:
            made = _random_molecule(rng, chem, basis, None)
            if made is None:
                raise ValueError("could not assemble random molecules within the limits; raise max_size or max_reductions")
            start.extend(made[: p.M - len(start)])
    if len(start) < 2:
        raise ValueError(f"the soup needs at least 2 molecules, got {len(start)}")

    def react(a, b):
        out = chem.react(a, b)
        return None if out is None else (a, b, *out[0])

    fired, pop = soup(react, start, p.generations * len(start), rng, arity=2, dilution="constant")
    found = list(start) + [m for _, rhs, _ in fired for m in rhs]
    species = _species(sorted(set(found)), basis, False)
    final = Counter(pop)
    return Network(
        species=species,
        reactions=[_reaction(lhs, rhs, False, count) for lhs, rhs, count in fired],
        status="observed",
        initial_state={m: float(n) for m, n in sorted(Counter(start).items())},
        outflow=CONSTANT_TOTAL,
        extras={"final_state": {m: n for m, n in final.most_common()},
                "analysis": {"final_diversity": len(final)}},
    )


def _reactive_soup(p, chem: Chemistry, basis: str, molecules: list[str], rng) -> Network:
    pool = Counter({x: p.atoms_per_type for x in basis})
    pop: list[str] = []
    for m in molecules:
        need = atom_counts(m)
        if any(need[x] > pool[x] for x in need):
            raise ValueError(f"atoms_per_type={p.atoms_per_type} is too small for the given molecules")
        pool -= need
        pop.append(m)
    if not molecules:
        while len(pop) < p.M:
            made = _random_molecule(rng, chem, basis, pool)
            if made is None:
                break
            for m in made:
                pool -= atom_counts(m)
            pop.extend(made)
    pop.sort()
    initial = Counter(pop)
    initial_pool = Counter(pool)

    fired: dict[tuple, list] = {}

    def record(lhs, rhs):
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        fired.setdefault(key, [tuple(lhs), tuple(rhs), 0])[2] += 1

    def remove(i):
        pop[i] = pop[-1]
        pop.pop()

    def p_add(n):
        return 1.0 if n <= p.min_molecules else 0.5 ** ((n - p.min_molecules) / p.half_add_prob)

    def insert():
        made = _random_molecule(rng, chem, basis, pool)
        if made is not None:
            for m in made:
                pool.subtract(atom_counts(m))
            pop.extend(made)
            record((), made)

    seen = set(pop)
    analysis = {"population": [len(pop)], "diversity": [len(set(pop))],
                "free_atoms": {x: [pool[x]] for x in basis}}
    t, generation = 0.0, 0
    while t < p.generations:
        n = len(pop)
        if n < 2:
            if rng.random() < p_add(n):
                insert()
            t += 1.0
        else:
            i, j = (int(v) for v in rng.choice(n, size=2, replace=False))
            a, b = pop[i], pop[j]
            out = chem.react(a, b)
            if out is not None:
                products, peak = out
                if all(pool[x] >= peak[x] for x in peak):
                    for k in sorted((i, j), reverse=True):
                        remove(k)
                    pool.update(atom_counts(a) + atom_counts(b))
                    for m in products:
                        pool.subtract(atom_counts(m))
                    pop.extend(products)
                    if Counter((a, b)) != Counter(products):
                        record((a, b), products)
            if pop and rng.random() < p.prob_destroy:
                k = int(rng.integers(len(pop)))
                m = pop[k]
                remove(k)
                pool.update(atom_counts(m))
                record((m,), ())
            if rng.random() < p_add(len(pop)) / max(n, 1):
                insert()
            t += 1.0 / n
        while generation + 1 <= t and generation < p.generations:
            generation += 1
            analysis["population"].append(len(pop))
            analysis["diversity"].append(len(set(pop)))
            for x in basis:
                analysis["free_atoms"][x].append(pool[x])

    names = sorted(seen | set(initial) | {m for lhs, rhs, _ in fired.values() for m in (*lhs, *rhs)})
    species = _species(names, basis, True)
    initial_state = {m: float(c) for m, c in sorted(initial.items())}
    initial_state.update({FREE + x: float(initial_pool[x]) for x in basis})
    final = Counter(pop)
    return Network(
        species=species,
        reactions=[_reaction(lhs, rhs, True, count) for lhs, rhs, count in fired.values()],
        status="observed",
        initial_state=initial_state,
        extras={
            "conservation": _conservation(species, basis),
            "final_state": {m: c for m, c in final.most_common()},
            "final_free_atoms": {x: int(pool[x]) for x in basis},
            "analysis": analysis,
        },
    )
