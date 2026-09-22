"""Well-stirred multiset reactor that records which reactions fire.

The generic algorithm of the book's constructive examples: draw `arity`
molecules at random, react them, replace reactants with products, optionally
dilute back to the initial population size, and repeat. The observed network
is the set of distinct reactions that fired, with firing counts.

`stir` runs the reactor as a generator of windows, so a chemistry's `evolve`
can yield a frame per window; `soup` runs it to the end. `Tally` is the
bookkeeping both use: firing counts keyed by reactant and product multisets.
"""

from __future__ import annotations

from collections import Counter
from typing import Callable, Hashable, Iterable

import numpy as np

Molecule = Hashable
DILUTION = ("none", "constant")


class Tally:
    """Firing counts of reactions, overall and since the last `flush`.

    Reactions are keyed by their reactant and product multisets, so the order
    of molecules does not matter; each keeps the sides it was first seen with.
    """

    def __init__(self) -> None:
        self.total: dict[tuple, list] = {}
        self.window: dict[tuple, list] = {}

    def add(self, lhs, rhs, n: int = 1) -> None:
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        for book in (self.total, self.window):
            entry = book.get(key)
            if entry is None:
                book[key] = [tuple(lhs), tuple(rhs), n]
            else:
                entry[2] += n

    def flush(self) -> list[list]:
        """The window as a frame's `fired`: [[reactants], [products], count]; then reset it."""
        out = [[list(lhs), list(rhs), n] for lhs, rhs, n in self.window.values()]
        self.window = {}
        return out

    def reactions(self) -> list[tuple[tuple, tuple, int]]:
        """Every reaction so far as (reactants, products, count), in order of first firing."""
        return [(lhs, rhs, n) for lhs, rhs, n in self.total.values()]


def stir(
    react: Callable[..., Iterable[Molecule] | None],
    population: Iterable[Molecule],
    steps: int,
    rng: np.random.Generator,
    *,
    every: int | None = None,
    arity: int = 2,
    dilution: str = "none",
    alternatives: bool = False,
    tally: Tally | None = None,
):
    """Run `steps` collisions, yielding (step, population, tally) at step 0,
    every `every` collisions (default: the initial population size, one
    generation) and at the end.

    The loop is `soup`'s, draw for draw. The population yielded is the live
    list, and the tally holds the reactions fired since the previous yield
    (`tally.flush()` gives them as a frame's `fired`). With dilution
    "constant" the population is diluted back to its initial size.
    """
    if dilution not in DILUTION:
        raise ValueError(f"dilution must be one of {DILUTION}, got {dilution!r}")
    pop = list(population)
    size = len(pop)
    every = every or max(size, 1)
    tally = tally if tally is not None else Tally()
    yield 0, pop, tally
    done = 0
    for done in range(1, steps + 1):
        if len(pop) < arity:
            break
        idx = _draw(rng, len(pop), arity)
        lhs = tuple(pop[i] for i in idx)
        out = react(*lhs)
        if out is not None:
            if alternatives:
                options = [tuple(o) for o in out]
                rhs = options[int(rng.integers(len(options)))] if options else None
            else:
                rhs = tuple(out)
            if rhs is not None and Counter(lhs) != Counter(rhs):
                for i in sorted(idx, reverse=True):
                    _swap_remove(pop, i)
                pop.extend(rhs)
                if dilution == "constant":
                    while len(pop) > size:
                        _swap_remove(pop, int(rng.integers(len(pop))))
                tally.add(lhs, rhs)
        if done % every == 0 and done < steps:
            yield done, pop, tally
    if steps:
        yield done, pop, tally


def soup(
    react: Callable[..., Iterable[Molecule] | None],
    population: Iterable[Molecule],
    steps: int,
    rng: np.random.Generator,
    arity: int = 2,
    dilution: str = "none",
    alternatives: bool = False,
) -> tuple[list[tuple[tuple, tuple, int]], list[Molecule]]:
    """Run `steps` collisions.

    react(*molecules) returns the complete right-hand side or None (elastic),
    as in `chemart.expand.expand`.

    dilution:
      "none"      reactants are replaced by products; the population may grow.
      "constant"  after each reaction, random molecules are removed until the
                  population is back to its initial size (e.g. the matrix
                  chemistry's removal of a random string s4).

    alternatives: if True, react returns an iterable of alternative
        right-hand sides (or None), and one of them is drawn uniformly with
        `rng` -- the stochastic counterpart of `expand(alternatives=True)`,
        which instead records every outcome as its own reaction. An empty
        iterable is treated as elastic.

    Note that with alternatives=False a returned list is taken as one
    right-hand side, whatever it contains. Molecules only have to be hashable,
    so tuples are legal molecules and a list of alternatives cannot be told
    apart from a right-hand side of tuple-valued molecules -- passing a
    multi-outcome rule without setting the flag silently injects the
    alternatives themselves into the population. Set the flag deliberately.

    Returns (reactions, final population), where reactions are
    (reactants, products, count), deduplicated as multisets, in order of first
    firing.
    """
    pop, tally = list(population), Tally()
    for _, pop, tally in stir(react, pop, steps, rng, arity=arity, dilution=dilution,
                              alternatives=alternatives, tally=tally):
        pass
    return tally.reactions(), pop


def _draw(rng: np.random.Generator, n: int, k: int) -> list[int]:
    """k distinct indices in range(n), without an O(n) permutation."""
    while True:
        idx = [int(i) for i in rng.integers(n, size=k)]
        if len(set(idx)) == k:
            return idx


def _swap_remove(items: list, i: int) -> None:
    items[i] = items[-1]
    items.pop()
