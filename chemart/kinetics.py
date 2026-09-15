"""Rate-law vocabulary and unit conversions.

A reaction's rate is a plain dict, ``{"law": <name>, <param>: <number>, ...}``.
The law must be one of `RATE_LAWS`, and the listed parameters must be
present; extra keys are allowed (e.g. ``"mode": "repression"`` for Hill).
"""

from __future__ import annotations

from math import factorial, prod
from typing import Any

#: law name -> required parameters
#: ("power" is k [X]^order for a reaction with a single reactant species, e.g.
#: the generalised selection equation x' = f x^c.)
RATE_LAWS: dict[str, tuple[str, ...]] = {
    "mass-action": ("k",),
    "power": ("k", "order"),
    "michaelis-menten": ("vmax", "km"),
    "hill": ("vmax", "K", "n"),
    "saturating": ("k", "K"),
    "arrhenius": ("A", "Ea"),
}

AVOGADRO = 6.02214076e23


def rate_problem(rate: Any) -> str | None:
    """Describe what is wrong with a rate dict, or return None if it is valid."""
    if not isinstance(rate, dict):
        return f"rate must be a dict, got {type(rate).__name__}"
    law = rate.get("law")
    if law not in RATE_LAWS:
        return f"rate law {law!r} not in {sorted(RATE_LAWS)}"
    missing = [k for k in RATE_LAWS[law] if k not in rate]
    if missing:
        return f"{law} rate is missing {missing}"
    for key, value in rate.items():
        if key == "law":
            continue
        if not isinstance(value, (int, float, str, bool)):
            return f"rate parameter {key!r} must be a JSON scalar, got {type(value).__name__}"
        if key in RATE_LAWS[law] and (isinstance(value, bool) or not isinstance(value, (int, float))):
            return f"rate parameter {key!r} must be a number, got {value!r}"
    return None


def k_to_c(k: float, reactants: dict[str, int], volume: float,
           avogadro: float = AVOGADRO) -> float:
    """Macroscopic rate constant k -> mesoscopic (stochastic) constant c.

    c = k / (N_A V)^(m-1) * prod_i l_i!, where m is the total number of
    reactant molecules and l_i the multiplicity of each reactant species.
    """
    m = sum(reactants.values())
    return k / (avogadro * volume) ** (m - 1) * prod(factorial(n) for n in reactants.values())
