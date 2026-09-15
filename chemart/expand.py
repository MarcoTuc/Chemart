"""Closure of a constructive chemistry: the generating operator G_C.

Book §12.6: start from a set of species, react every combination, add any
novel product, and repeat until nothing new appears. Each round only tries
combinations that involve at least one species discovered in the previous
round, so no combination is ever evaluated twice.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Callable, Hashable, Iterable, Iterator

Molecule = Hashable
ReactFn = Callable[..., Iterable[Molecule] | None]


def expand(
    react: ReactFn,
    seed: Iterable[Molecule],
    arity: int | Iterable[int] = 2,
    max_species: int = 1000,
    ordered: bool = True,
) -> tuple[list[Molecule], list[tuple[tuple, tuple]], str]:
    """Compute the reaction network reachable from `seed`.

    react(*molecules) returns the complete right-hand side of the reaction
    (reactants that survive, e.g. catalysts, must be listed) or None for an
    elastic collision. A right-hand side equal to the left-hand side is also
    treated as elastic.

    arity: reaction order(s) to try.
    ordered: True if react(a, b) and react(b, a) can differ (operator and
        operand); False to try each multiset of reactants once.
    max_species: budget. A reaction whose novel products would exceed it is
        dropped, and the status becomes "truncated".

    Returns (species in discovery order, reactions as (reactants, products)
    tuples deduplicated as multisets, status "complete" | "truncated").
    """
    arities = (arity,) if isinstance(arity, int) else tuple(arity)
    species = list(dict.fromkeys(seed))
    truncated = len(species) > max_species
    species = species[:max_species]
    known = set(species)
    reactions: dict[tuple, tuple[tuple, tuple]] = {}

    done = 0
    while done < len(species):
        old, n = done, len(species)
        done = n
        for k in arities:
            for combo in _combinations(old, n, k, ordered):
                lhs = tuple(species[i] for i in combo)
                out = react(*lhs)
                if out is None:
                    continue
                rhs = tuple(out)
                left, right = Counter(lhs), Counter(rhs)
                if left == right:
                    continue
                novel = [m for m in dict.fromkeys(rhs) if m not in known]
                if len(species) + len(novel) > max_species:
                    truncated = True
                    continue
                species.extend(novel)
                known.update(novel)
                key = (frozenset(left.items()), frozenset(right.items()))
                reactions.setdefault(key, (lhs, rhs))
    return species, list(reactions.values()), "truncated" if truncated else "complete"


def _combinations(old: int, n: int, k: int, ordered: bool) -> Iterator[tuple[int, ...]]:
    """Index tuples over range(n) of length k with at least one index >= old."""
    if ordered:
        # Split on the first position that holds a new index.
        for first in range(k):
            yield from product(
                *([range(old)] * first), range(old, n), *([range(n)] * (k - first - 1))
            )
    else:
        # Non-decreasing tuples whose last (largest) index is new.
        for last in range(old, n):
            for prefix in combinations_with_replacement(range(last + 1), k - 1):
                yield (*prefix, last)
