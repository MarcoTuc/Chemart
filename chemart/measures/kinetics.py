"""F. Kinetics: measures that need rate constants."""

from __future__ import annotations

import numpy as np

from chemart.measures import register


def _mass_action_k(ctx) -> list[float]:
    return [float(r.rate["k"]) for r in ctx.net.reactions
            if r.rate["law"] == "mass-action" and float(r.rate["k"]) > 0]


@register("rate_spread", "F", needs="K")
def rate_spread(ctx) -> float | None:
    """Orders of magnitude spanned by the mass-action rate constants,
    log10(max k / min k). None without mass-action rates."""
    k = _mass_action_k(ctx)
    return float(np.log10(max(k) / min(k))) if k else None


def reversible_pairs(net) -> list[tuple[int, int]]:
    """Index pairs (forward, backward) of reactions that are each other's reverse."""
    seen: dict[tuple, int] = {}
    pairs = []
    for j, r in enumerate(net.reactions):
        key = (frozenset(r.reactants.items()), frozenset(r.products.items()))
        back = (key[1], key[0])
        if back in seen:
            pairs.append((seen.pop(back), j))
        else:
            seen[key] = j
    return pairs


@register("wegscheider_residual", "F", needs="K S")
def wegscheider_residual(ctx) -> float | None:
    """How far the reversible mass-action pairs are from allowing detailed
    balance (Wegscheider's conditions): the least-squares residual of
    ln(k+/k-) against the reactions' stoichiometry, 0 when some chemical
    potentials make every pair balance. None without reversible
    mass-action pairs."""
    rows, target = [], []
    for f, b in reversible_pairs(ctx.net):
        rf, rb = ctx.net.reactions[f].rate, ctx.net.reactions[b].rate
        if rf["law"] == rb["law"] == "mass-action" and rf["k"] > 0 and rb["k"] > 0:
            rows.append(ctx.S[:, f])
            target.append(np.log(rf["k"] / rb["k"]))
    if not rows:
        return None
    A, y = np.array(rows, dtype=float), np.array(target)
    # detailed balance at amounts x*: k+ x*^R = k- x*^P, i.e. ln(k+/k-) = S_jᵀ ln x*
    mu, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.linalg.norm(A @ mu - y))
