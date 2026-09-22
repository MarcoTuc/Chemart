"""G. Dynamics: measures of a population state, and of a whole trajectory.

Population measures take one frame's state ({species: amount}); `track` and
`over` apply them frame by frame. Trajectory measures summarise a run.
"""

from __future__ import annotations

import numpy as np

from chemart.measures import register


def _amounts(state: dict) -> np.ndarray:
    x = np.array([v for v in state.values() if v > 0], dtype=float)
    return x


@register("richness", "G", input="state")
def richness(state) -> int:
    """Number of species present."""
    return int(sum(1 for v in state.values() if v > 0))


@register("shannon", "G", input="state")
def shannon(state) -> float:
    """Shannon diversity (nats) of the population's composition."""
    x = _amounts(state)
    if not len(x):
        return 0.0
    p = x / x.sum()
    return float(-(p * np.log(p)).sum())


@register("dominance", "G", input="state")
def dominance(state) -> float:
    """Share of the most abundant species (Berger-Parker index)."""
    x = _amounts(state)
    return float(x.max() / x.sum()) if len(x) else 0.0


@register("population", "G", input="state")
def population(state) -> float:
    """Total amount of all species."""
    return float(sum(state.values()))


@register("turnover", "G", input="trajectory", needs="D")
def turnover(traj) -> float | None:
    """Mean Jaccard distance between the species sets of consecutive frames:
    how fast the population's composition changes. None with one frame."""
    sets = [{s for s, v in f.state.items() if v > 0} for f in traj.frames]
    d = [1 - len(a & b) / len(a | b) for a, b in zip(sets, sets[1:]) if a | b]
    return float(np.mean(d)) if d else None


@register("collapse_time", "G", input="trajectory", needs="D")
def collapse_time(traj, fraction: float = 0.1) -> float | None:
    """First time, after its peak, at which richness falls to `fraction` (10%)
    of the peak. None if it never does."""
    rich = [richness(f.state) for f in traj.frames]
    peak = int(np.argmax(rich))
    for f, r in zip(traj.frames[peak:], rich[peak:]):
        if r <= fraction * rich[peak]:
            return float(f.t)
    return None


@register("final_richness_ratio", "G", input="trajectory", needs="D")
def final_richness_ratio(traj) -> float:
    """Richness at the end over the peak richness of the run."""
    rich = [richness(f.state) for f in traj.frames]
    return rich[-1] / max(rich) if max(rich) else 0.0
