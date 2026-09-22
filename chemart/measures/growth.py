"""E. Constructiveness and growth: whether a chemistry keeps making new species.

Growth exponents across sizes are `chemart.measures.scaling` over a `sweep`;
these measures read one network or one run.
"""

from __future__ import annotations

import numpy as np

from chemart.measures import register


@register("mean_structure_length", "E", needs="str")
def mean_structure_length(ctx) -> float:
    """Mean length of the species' structure strings (a λ-term, a sequence, a
    fold): a crude size of the molecules in their own representation."""
    lengths = [len(s.structure) for s in ctx.net.species if s.structure]
    return float(np.mean(lengths))


@register("novelty_rate", "E", input="trajectory", needs="D")
def novelty_rate(traj) -> float | None:
    """Species never seen before, per unit of the trajectory's clock."""
    seen = set(traj.frames[0].state)
    new = 0
    for f in traj.frames[1:]:
        fresh = set(f.state) - seen
        new += len(fresh)
        seen |= fresh
    span = traj.frames[-1].t - traj.frames[0].t
    return new / span if span > 0 else None


@register("complexity_drift", "E", input="trajectory", needs="D")
def complexity_drift(traj) -> float | None:
    """Slope over time of the abundance-weighted mean length of species ids: do
    molecules get bigger as the run goes on?"""
    t, y = [], []
    for f in traj.frames:
        total = sum(f.state.values())
        if total > 0:
            t.append(f.t)
            y.append(sum(len(s) * x for s, x in f.state.items()) / total)
    if len(t) < 2 or len(set(t)) < 2:
        return None
    return float(np.polyfit(t, y, 1)[0])
