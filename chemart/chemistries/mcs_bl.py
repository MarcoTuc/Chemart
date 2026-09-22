"""MCS.bl: the Molecular Classifier System based on Holland's Broadcast Language.

Catalog id: mcs-bl. Book 11.1.1; Decraene & McMullin (thesis 2009, ACS 2011).

Molecules are broadcast devices, strings over the alphabet {0, 1, *, :, #, $,
%, '} (thesis glyphs: # = diamond, $ = reversed triangle, % = triangle,
' = prime). A device holds zero or more broadcast units "*condition:action".
In a collision the first molecule acts as the enzyme and the second as the
substrate: if a unit's condition matches the whole substrate string, its action
is written out as a product and both reactants survive:

    enzyme + substrate -> enzyme + substrate + product

Matching and action follow the thesis (4.3) and the authors' C++ Broadcast
Language implementation (technical report ALL-06-01, class BUnit):

- '0', '1' and quoted symbols are literals; an unquoted quote only marks the
  next symbol, and is dropped from both condition and action.
- '#' matches any one symbol; as the last symbol of the condition it matches
  any non-empty suffix.
- '$' matches a non-empty prefix (first position) or suffix (last position);
  anywhere else it is ignored, and at both ends only the first counts. In the
  action every '$' is replaced by the matched string (transposition).
- '%' matches one symbol; only the leftmost '%' counts. In the action every
  '%' is replaced by the matched symbol.
- '$' and '%' in the action without a counterpart in the condition, and
  unquoted '#' in the action, are ignored.

Two faces: `generate` gives every reaction reachable from the seed strings
(chemart.expand); `evolve` runs the thesis' single-reactor model with
per-symbol mutation, a frame per generation of the initial population, and
returns the reactions that fired.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame

ALPHABET = "01*:#$%'"
# Glyphs used in the thesis and papers -> ASCII symbols used by the book's text.
GLYPHS = {"∗": "*", "♦": "#", "⋄": "#", "◊": "#", "▽": "$", "△": "%", "′": "'"}
_ID = {"0": "0", "1": "1", "*": "S", ":": "C", "#": "H", "$": "D", "%": "P", "'": "Q"}

# Book table 11.2 / thesis table 7.1: the seed cell type c0.
C0 = ("*$0:$1", "*$0:$0", "*$1:$0", "*$1:$1")

_LIT, _ANY, _TAIL, _VAR, _BIND = range(5)


def to_ascii(text: str) -> str:
    """Translate thesis glyphs (▽ ♦ △ ′ ∗) to the ASCII alphabet and drop spaces."""
    return "".join(GLYPHS.get(c, c) for c in text if not c.isspace())


def species_id(device: str) -> str:
    """Safe species id: m_ + one letter per symbol (* S, : C, # H, $ D, % P, ' Q)."""
    return "m_" + "".join(_ID[c] for c in device)


def _quoted(s: str, i: int) -> bool:
    # BUnit.h isQuoted: a symbol is quoted if the symbol before it is a quote.
    return i > 0 and s[i - 1] == "'"


def units(device: str) -> list[tuple[str, str]]:
    """Broadcast units (condition, action) of a device, left to right (thesis 4.3.2).

    Symbols before the first unquoted '*' are junk; a unit runs to the next
    unquoted '*'; a second unquoted ':' and everything after it are ignored;
    condition and action must both be non-empty.
    """
    stars = [i for i, c in enumerate(device) if c == "*" and not _quoted(device, i)]
    out = []
    for k, start in enumerate(stars):
        end = stars[k + 1] if k + 1 < len(stars) else len(device)
        body = device[start + 1:end]
        colons = [j for j, c in enumerate(body) if c == ":" and not _quoted(body, j)]
        if not colons:
            continue
        cond = body[:colons[0]]
        act = body[colons[0] + 1:colons[1] if len(colons) > 1 else len(body)]
        if cond and act:
            out.append((cond, act))
    return out


@lru_cache(maxsize=None)
def _compile(cond: str, act: str) -> tuple[tuple, tuple]:
    n = len(cond)
    var_at = 0 if cond[0] == "$" else (n - 1 if cond[-1] == "$" and not _quoted(cond, n - 1) else None)
    bind_at = next((i for i, c in enumerate(cond) if c == "%" and not _quoted(cond, i)), None)
    pattern = []
    for i, c in enumerate(cond):
        if _quoted(cond, i) or c in "01":
            pattern.append((_LIT, c))
        elif c == "$":
            if i == var_at:
                pattern.append((_VAR, None))
        elif c == "%":
            if i == bind_at:
                pattern.append((_BIND, None))
        elif c == "#":
            pattern.append((_TAIL if i == n - 1 else _ANY, None))
    out: list = []
    for i, c in enumerate(act):
        if _quoted(act, i) or c in "01":
            out.append(c)
        elif c == "$" and var_at is not None:
            out.append(_VAR)
        elif c == "%" and bind_at is not None:
            out.append(_BIND)
    return tuple(pattern), tuple(out)


def fire(cond: str, act: str, substrate: str) -> str | None:
    """The product of one broadcast unit on a substrate, or None if the condition fails."""
    pattern, out = _compile(cond, act)
    m = len(pattern)
    if m == 0:
        return None
    value = None
    if pattern[0][0] == _VAR:
        rest = pattern[1:]
        offset = len(substrate) - len(rest)
        if offset < 1:
            return None
        value = substrate[:offset]
    elif pattern[-1][0] == _VAR:
        rest = pattern[:-1]
        offset = 0
        if len(substrate) < len(rest) + 1:
            return None
        value = substrate[len(rest):]
    elif pattern[-1][0] == _TAIL:
        rest, offset = pattern, 0
        if len(substrate) < m:
            return None
    else:
        rest, offset = pattern, 0
        if len(substrate) != m:
            return None
    bind = None
    for j, (kind, c) in enumerate(rest):
        x = substrate[offset + j]
        if kind == _LIT and x != c:
            return None
        if kind == _BIND:
            bind = x
    return "".join(value if t == _VAR else bind if t == _BIND else t
                   for t in out) if out else ""


@lru_cache(maxsize=None)
def candidates(enzyme: str, substrate: str) -> tuple[str, ...]:
    """One product per unit of the enzyme whose condition the substrate satisfies."""
    return tuple(p for cond, act in units(enzyme)
                 if (p := fire(cond, act, substrate)) is not None)


class Rules:
    """Which products a collision can give once the reactor's filters apply."""

    def __init__(self, max_length: int, self_replication: bool):
        self.max_length = max_length
        self.self_replication = self_replication
        self._cache: dict[tuple[str, str], tuple[int, Counter]] = {}

    def allowed(self, enzyme: str, substrate: str, product: str) -> bool:
        if not product or len(product) > self.max_length:
            return False
        return self.self_replication or not (product == enzyme == substrate)

    def outcomes(self, enzyme: str, substrate: str) -> tuple[int, Counter]:
        """(number of satisfied units, Counter of allowed products over those units)."""
        key = (enzyme, substrate)
        if key not in self._cache:
            cands = candidates(enzyme, substrate)
            ok = Counter(p for p in cands if self.allowed(enzyme, substrate, p))
            self._cache[key] = (len(cands), ok)
        return self._cache[key]

    def products(self, enzyme: str, substrate: str) -> list[str]:
        return list(self.outcomes(enzyme, substrate)[1])

    def k(self, a: str, b: str, product: str) -> float:
        """Collision probability weight of {a, b} -> {a, b, product}, summed over both orders."""
        total = 0.0
        for e, s in {(a, b), (b, a)}:
            n, ok = self.outcomes(e, s)
            if ok.get(product):
                total += ok[product] / n
        return total


def mutate(device: str, p_s: float, rng) -> str:
    """Each symbol mutates with probability p_s: flip, insert after, or delete (equal odds)."""
    n = len(device)
    hits = int(rng.binomial(n, p_s)) if n else 0
    if not hits:
        return device
    where = {int(i) for i in rng.choice(n, size=hits, replace=False)}
    out = []
    for i, c in enumerate(device):
        if i not in where:
            out.append(c)
            continue
        kind = int(rng.integers(3))
        if kind == 0:
            others = ALPHABET.replace(c, "")
            out.append(others[int(rng.integers(len(others)))])
        elif kind == 1:
            out.append(c)
            out.append(ALPHABET[int(rng.integers(len(ALPHABET)))])
    return "".join(out)


def _strings(p, rng) -> list[str]:
    seeds = []
    for s in p.strings:
        if not isinstance(s, str):
            raise ValueError(f"strings must be broadcast-device strings, got {s!r}")
        s = to_ascii(s)
        if not s or any(c not in ALPHABET for c in s):
            raise ValueError(f"strings must be non-empty strings over {ALPHABET!r}, got {s!r}")
        if len(s) > p.max_length:
            raise ValueError(f"string {s!r} is longer than max_length={p.max_length}")
        seeds.append(s)
    if p.n_random:
        if p.random_length > p.max_length:
            raise ValueError(f"random_length={p.random_length} exceeds max_length={p.max_length}")
        draws = rng.integers(len(ALPHABET), size=(p.n_random, p.random_length))
        seeds += ["".join(ALPHABET[int(i)] for i in row) for row in draws]
    if not seeds:
        raise ValueError("give seed strings or n_random > 0")
    return seeds


def generate(p, rng):
    """Every reaction reachable from the seed strings, cut off by max_species."""
    seeds = _strings(p, rng)
    return _closure(p, Rules(p.max_length, p.self_replication), seeds)


def _reaction(rules: Rules, a: str, b: str, product: str, count: int | None = None) -> Reaction:
    k = rules.k(a, b, product)
    rate = {"law": "mass-action", "k": k} if k > 0 else None
    return Reaction.of([species_id(a), species_id(b)],
                       [species_id(a), species_id(b), species_id(product)], rate=rate, count=count)


def _network(species, reactions, status, seeds, p, extras) -> Network:
    start = Counter(seeds)
    return Network(
        species=[Species(species_id(s), structure=s) for s in species],
        reactions=reactions,
        status=status,
        initial_state={species_id(s): float(n * p.initial_copies) for s, n in start.items()},
        outflow=CONSTANT_TOTAL,
        extras={"seed": [species_id(s) for s in start], **extras},
    )


def _closure(p, rules: Rules, seeds: list[str]) -> Network:
    def react(a, b):
        prods = rules.products(a, b)
        return (a, b, *prods) if prods else None

    found, _, status = expand(react, list(dict.fromkeys(seeds)), arity=2,
                              max_species=p.max_species, ordered=True)
    known = set(found)
    seen: dict[tuple, tuple[str, str, str]] = {}
    for a in found:
        for b in found:
            for prod in rules.products(a, b):
                if prod in known:
                    key = (min(a, b), max(a, b), prod)
                    seen.setdefault(key, (a, b, prod))
    reactions = [_reaction(rules, a, b, prod) for a, b, prod in seen.values()]
    return _network(found, reactions, status, seeds, p, {})


def evolve(p, rng):
    """The thesis' single reactor (4.2.1): a frame per generation of the initial population.

    The reactor is its own loop rather than chemart.soup.stir: below n_max the
    product is added, at n_max it overwrites a random molecule other than the
    two reactants in place, and the product is mutated before it enters.
    """
    seeds = _strings(p, rng)
    rules = Rules(p.max_length, p.self_replication)
    pop = [s for s in seeds for _ in range(p.initial_copies)]
    if len(pop) < 2:
        raise ValueError("the soup needs at least 2 molecules (strings x initial_copies)")
    if len(pop) > p.n_max:
        raise ValueError(f"{len(pop)} initial molecules exceed the capacity n_max={p.n_max}")
    size = len(pop)
    order = dict.fromkeys(pop)
    tally = Tally()
    mutants = 0

    def collide() -> None:
        nonlocal mutants
        i = j = 0
        while i == j:
            i, j = (int(x) for x in rng.integers(len(pop), size=2))
        enzyme, substrate = pop[i], pop[j]
        cands = candidates(enzyme, substrate)
        if not cands:
            return
        product = cands[int(rng.integers(len(cands)))] if len(cands) > 1 else cands[0]
        if not rules.allowed(enzyme, substrate, product):
            return
        if p.p_s > 0:
            mutated = mutate(product, p.p_s, rng)
            mutants += mutated != product
            product = mutated
            if not product or len(product) > p.max_length:
                return
        if len(pop) < p.n_max:
            pop.append(product)
        else:   # displace a random molecule other than the two reactants
            while (x := int(rng.integers(len(pop)))) in (i, j):
                pass
            pop[x] = product
        order.setdefault(product)
        tally.add((enzyme, substrate), (enzyme, substrate, product))

    def frame(step: int) -> Frame:
        fired = [[[species_id(s) for s in lhs], [species_id(s) for s in rhs], n] for lhs, rhs, n in tally.flush()]
        return Frame(t=float(step), state={species_id(s): float(n) for s, n in Counter(pop).items()}, fired=fired)

    yield frame(0)
    for done in range(1, p.steps + 1):
        collide()
        if done % size == 0 and done < p.steps:
            yield frame(done)
    if p.steps:
        yield frame(p.steps)

    fired = tally.reactions()
    reactions = [_reaction(rules, a, b, rhs[2], count) for (a, b), rhs, count in fired]
    final = Counter(pop)
    return _network(list(order), reactions, "observed", seeds, p, {
        "final_state": {species_id(s): n for s, n in final.most_common()},
        "analysis": {"collisions": p.steps, "productive": sum(n for _, _, n in fired),
                     "mutant_products": mutants},
    })
