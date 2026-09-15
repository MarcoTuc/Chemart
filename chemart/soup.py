"""Well-stirred multiset reactor that records which reactions fire.

The generic algorithm of the book's constructive examples: draw `arity`
molecules at random, react them, replace reactants with products, optionally
dilute back to the initial population size, and repeat. The observed network
is the set of distinct reactions that fired, with firing counts.
"""

from __future__ import annotations

from collections import Counter
from typing import Callable, Hashable, Iterable

import numpy as np

Molecule = Hashable
DILUTION = ("none", "constant")


def soup(
    react: Callable[..., Iterable[Molecule] | None],
    population: Iterable[Molecule],
    steps: int,
    rng: np.random.Generator,
    arity: int = 2,
    dilution: str = "none",
) -> tuple[list[tuple[tuple, tuple, int]], list[Molecule]]:
    """Run `steps` collisions.

    react(*molecules) returns the complete right-hand side or None (elastic),
    as in `chemart.expand.expand`.

    dilution:
      "none"      reactants are replaced by products; the population may grow.
      "constant"  after each reaction, random molecules are removed until the
                  population is back to its initial size (e.g. the matrix
                  chemistry's removal of a random string s4).

    Returns (reactions, final population), where reactions are
    (reactants, products, count), deduplicated as multisets, in order of first
    firing.
    """
    if dilution not in DILUTION:
        raise ValueError(f"dilution must be one of {DILUTION}, got {dilution!r}")
    pop = list(population)
    size = len(pop)
    fired: dict[tuple, list] = {}

    for _ in range(steps):
        if len(pop) < arity:
            break
        idx = _draw(rng, len(pop), arity)
        lhs = tuple(pop[i] for i in idx)
        out = react(*lhs)
        if out is None:
            continue
        rhs = tuple(out)
        left, right = Counter(lhs), Counter(rhs)
        if left == right:
            continue
        for i in sorted(idx, reverse=True):
            _swap_remove(pop, i)
        pop.extend(rhs)
        if dilution == "constant":
            while len(pop) > size:
                _swap_remove(pop, int(rng.integers(len(pop))))
        key = (frozenset(left.items()), frozenset(right.items()))
        entry = fired.setdefault(key, [lhs, rhs, 0])
        entry[2] += 1

    return [(lhs, rhs, count) for lhs, rhs, count in fired.values()], pop


def _draw(rng: np.random.Generator, n: int, k: int) -> list[int]:
    """k distinct indices in range(n), without an O(n) permutation."""
    while True:
        idx = [int(i) for i in rng.integers(n, size=k)]
        if len(set(idx)) == k:
            return idx


def _swap_remove(items: list, i: int) -> None:
    items[i] = items[-1]
    items.pop()
