"""Prime number (number-division) chemistry (book 2.5.2; appendix NumberChem.py).

Catalog id: prime-number-chemistry.

Molecules are integers >= 2. Two molecules collide; if the smaller divides
the larger (and is strictly smaller), the larger is replaced by the quotient
and the divisor acts as a catalyst:

    s1 + s2 -> s1 + s2/s1      (s1 < s2, s1 | s2; otherwise elastic)

Two faces. `generate` is the reaction closure of the distinct seed numbers
(always finite: every product divides an existing number), truncated only by
`max_species`. `evolve` is the book's run (NumberChem.py): M integers drawn
uniformly from [minn, maxn], `iterations` collisions of two distinct random
molecules, a frame per generation (M collisions) with the prime fraction, and
the observed reactions with firing counts at the end.
"""

from collections import Counter

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame


def divide(a: int, b: int):
    """react(a, b): order the pair, divide the larger by the smaller if possible."""
    small, large = (a, b) if a <= b else (b, a)
    if large > small and large % small == 0:
        return (small, large // small)
    return None


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


def sid(n: int) -> str:
    return f"n{n}"


def _initial(p, rng) -> list[int]:
    if p.numbers:
        bad = [v for v in p.numbers if not isinstance(v, int) or isinstance(v, bool) or v < 2]
        if bad:
            raise ValueError(f"numbers must be integers >= 2 (0 and 1 are excluded), got {bad!r}")
        return list(p.numbers)
    if p.maxn < p.minn:
        raise ValueError(f"maxn must be >= minn, got minn={p.minn}, maxn={p.maxn}")
    # NumberChem.py: np.random.randint(minn, maxn+1), both ends inclusive.
    return [int(v) for v in rng.integers(p.minn, p.maxn + 1, size=p.M)]


def _reaction(lhs, rhs, count=None) -> Reaction:
    return Reaction.of([sid(n) for n in lhs], [sid(n) for n in rhs], count=count)


def _prime_fraction(pop) -> float:
    return sum(1 for n in pop if is_prime(n)) / len(pop)


def _state(pop) -> dict[str, float]:
    return {sid(n): float(c) for n, c in sorted(Counter(pop).items())}


def generate(p, rng):
    """The closure of the distinct seed numbers (a Chemart addition), cut off by max_species."""
    start = _initial(p, rng)
    seed = sorted(set(start))
    numbers, reactions, status = expand(divide, seed, arity=2, max_species=p.max_species, ordered=False)
    numbers = sorted(numbers)
    return Network(
        species=[Species(sid(n), structure=str(n)) for n in numbers],
        reactions=[_reaction(lhs, rhs) for lhs, rhs in reactions],
        status=status,
        extras={
            "seed": [sid(n) for n in seed],
            "primes": [sid(n) for n in numbers if is_prime(n)],
        },
    )


def evolve(p, rng):
    """The book's soup (NumberChem.py): a frame per generation of M collisions."""
    start = _initial(p, rng)
    if len(start) < 2:
        raise ValueError(f"the soup needs at least 2 molecules, got {len(start)}")
    size = len(start)
    tally = Tally()
    for step, pop, tally in stir(divide, start, p.iterations, rng, arity=2, tally=tally):
        fired = [[[sid(n) for n in lhs], [sid(n) for n in rhs], count] for lhs, rhs, count in tally.flush()]
        yield Frame(t=float(step), state=_state(pop), fired=fired,
                    observables={"prime_fraction": _prime_fraction(pop)})

    events = tally.reactions()
    numbers = sorted({*start, *(n for _, rhs, _ in events for n in rhs)})
    reactions = [_reaction(lhs, rhs, count) for lhs, rhs, count in events]
    final = Counter(pop)
    return Network(
        species=[Species(sid(n), structure=str(n)) for n in numbers],
        reactions=reactions,
        status="observed",
        initial_state={sid(n): c for n, c in sorted(Counter(start).items())},
        extras={
            "analysis": {
                "generation_size": size,
                "effective_collisions": sum(r.count for r in reactions),
                "new_numbers": [n for n in numbers if n not in set(start)],
            },
            "final_state": {sid(n): c for n, c in sorted(final.items())},
            "primes": [sid(n) for n in numbers if is_prime(n)],
        },
    )
