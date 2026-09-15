"""Build a Network from written-down reactions.

    network([
        ("A + B -> C", 1.0),                                  # mass action, k = 1.0
        ("2 X + Y -> 3 X", {"law": "mass-action", "k": 2.0}),
        ("X -> ", None),                                      # no rate; empty side = ∅
    ], initial_state={"A": 2.0})

Terms are separated by " + " (with spaces), so species names may contain
"+" ("e+"). A stoichiometric coefficient is separated from the species by a
space ("2 X"), so names may start with digits ("12C", "4He").
"""

from __future__ import annotations

import re
from typing import Any, Iterable

from chemart.network import Network, Reaction, Species

_TERM = re.compile(r"^(\d+)\s+(\S+)$")


def term(n: int, species: str) -> str:
    """Format one side term: term(1, "P") -> "P", term(3, "P") -> "3 P"."""
    return species if n == 1 else f"{n} {species}"


def parse(text: str) -> tuple[dict[str, int], dict[str, int]]:
    if text.count("->") != 1:
        raise ValueError(f"reaction {text!r} needs exactly one '->'")
    left, right = text.split("->")
    return _side(left, text), _side(right, text)


def _side(side: str, text: str) -> dict[str, int]:
    out: dict[str, int] = {}
    side = side.strip()
    if side in ("", "∅"):
        return out
    for raw in side.split(" + "):
        item = raw.strip()
        match = _TERM.match(item)
        n, name = (int(match.group(1)), match.group(2)) if match else (1, item)
        if not name or " " in name or n < 1:
            raise ValueError(f"cannot parse term {item!r} in reaction {text!r}")
        out[name] = out.get(name, 0) + n
    return out


def rate(value: Any) -> dict[str, Any] | None:
    """None -> no rate; a number -> mass action with that k; a dict -> as given."""
    if value is None:
        return None
    if isinstance(value, dict):
        return dict(value)
    return {"law": "mass-action", "k": float(value)}


def network(
    reactions: Iterable[tuple[str, Any]],
    species: Iterable[str | Species] = (),
    **kwargs: Any,
) -> Network:
    """Species come first in the order given, then in order of appearance."""
    order: dict[str, Species] = {}
    for s in species:
        sp = s if isinstance(s, Species) else Species(s)
        order[sp.id] = sp
    built = []
    for text, r in reactions:
        lhs, rhs = parse(text)
        for name in (*lhs, *rhs):
            order.setdefault(name, Species(name))
        built.append(Reaction(lhs, rhs, rate(r)))
    for key in ("initial_state", "inflow"):
        for name in kwargs.get(key) or {}:
            order.setdefault(name, Species(name))
    return Network(species=list(order.values()), reactions=built, **kwargs)
