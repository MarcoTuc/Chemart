"""Checks for list/dict parameters, whose inner shape the catalog cannot express."""

from __future__ import annotations

from typing import Any


def _is_number(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def vector(name: str, value: Any, length: int | None = None, *,
           minimum: float | None = None) -> list[float]:
    if not isinstance(value, list) or not all(_is_number(v) for v in value):
        raise ValueError(f"{name} must be a list of numbers, got {value!r}")
    if length is not None and len(value) != length:
        raise ValueError(f"{name} must have {length} entries, got {len(value)}")
    if minimum is not None and any(v < minimum for v in value):
        raise ValueError(f"{name} entries must be >= {minimum}, got {value!r}")
    return [float(v) for v in value]


def square_matrix(name: str, value: Any, n: int | None = None) -> list[list[float]]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a non-empty list of rows, got {value!r}")
    size = len(value) if n is None else n
    if len(value) != size:
        raise ValueError(f"{name} must be {size}x{size}, got {len(value)} rows")
    return [vector(f"{name}[{i}]", row, size) for i, row in enumerate(value)]


def edges(name: str, value: Any) -> list[tuple[int, int]]:
    """Undirected edge list [[u, v], ...] of integer vertex ids, no self-loops."""
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a non-empty list of [u, v] pairs, got {value!r}")
    out = []
    for e in value:
        if (not isinstance(e, list) or len(e) != 2
                or not all(isinstance(v, int) and not isinstance(v, bool) for v in e)):
            raise ValueError(f"{name}: every edge must be a pair of integers, got {e!r}")
        if e[0] == e[1]:
            raise ValueError(f"{name}: self-loop {e!r} is not allowed")
        out.append((e[0], e[1]))
    return out


def apportion(total: int, weights: list[float]) -> list[int]:
    """Non-negative integers summing to `total`, proportional to `weights`."""
    s = sum(weights)
    if s <= 0:
        raise ValueError(f"weights must have a positive sum, got {weights!r}")
    raw = [total * w / s for w in weights]
    base = [int(r) for r in raw]
    by_remainder = sorted(range(len(raw)), key=lambda i: raw[i] - base[i], reverse=True)
    for i in by_remainder[: total - sum(base)]:
        base[i] += 1
    return base
