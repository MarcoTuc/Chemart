"""French flag model: threshold differentiation along a morphogen gradient (book 18.6).

Catalog id: french-flag. Each cell position i is an undifferentiated species
U_i that differentiates into the fate selected by the morphogen
concentration at i; the gradient itself is imposed (maternal).
"""

import math

from chemart.helpers.explicit import network
from chemart.helpers.params import vector


def generate(p, rng):
    thresholds = sorted(vector("thresholds", p.thresholds, minimum=0.0), reverse=True)
    if not thresholds:
        raise ValueError("thresholds must list at least one concentration")
    fates = ["blue", "white", "red"] if len(thresholds) == 2 else [f"type{k}" for k in range(len(thresholds) + 1)]
    width = p.n_cells
    height = p.n_cells if p.dimensions == 2 else 1

    reactions, cells, morphogen = [], [], {}
    for y in range(height):
        for x in range(width):
            suffix = f"_{x}" if p.dimensions == 1 else f"_{x}_{y}"
            position = x / (width - 1)
            c = 1.0 - position if p.gradient == "linear" else math.exp(-position / p.decay_length)
            fate = next((k for k, t in enumerate(thresholds) if c >= t), len(thresholds))
            cells.append(f"U{suffix}")
            morphogen[f"U{suffix}"] = c
            reactions.append((f"U{suffix} -> {fates[fate]}{suffix}", None))
    return network(
        reactions,
        species=cells,
        extras={"space": {"dimensions": p.dimensions, "shape": [width, height], "morphogen": morphogen}},
    )
