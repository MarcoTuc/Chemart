"""Regenerate the measure tables of docs/guide/measures.md from the registry.

    uv run python tools/gen_measures_page.py           # rewrite the tables
    uv run python tools/gen_measures_page.py --check   # exit 1 if they are stale

Each section of the page has a block between `<!-- measures X -->` and
`<!-- /measures -->`; the block is rebuilt from `chemart.measures.REGISTRY`,
so the page cannot drift from the code. Measures that are specified on the
page but not implemented yet are listed from PLANNED, marked as planned.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from chemart.measures import REGISTRY  # noqa: E402

PAGE = ROOT / "docs" / "guide" / "measures.md"

#: (name, meaning, needs, input, cost) of the measures still to implement.
PLANNED: dict[str, list[tuple[str, str, str, str, str]]] = {
    "B": [
        ("p_invariants", "the semi-positive conservation laws (conserved moieties): extreme rays of {m ≥ 0, mᵀS = 0}", "S", "network", "exponential"),
        ("elementary_flux_modes", "number and mean length of the minimal pathways through the network with its food boundary", "S F", "network", "exponential"),
        ("blocked_fraction", "share of reactions that can never carry flux at steady state (flux variability analysis)", "S F", "network", "moderate"),
    ],
    "C": [
        ("modularity", "Louvain modularity Q of the bipartite graph, relative to the null model", "T", "network", "moderate"),
        ("motif_profile", "significance profile of the triads of the substrate -> product graph (Milo et al. 2004)", "T", "network", "moderate"),
    ],
    "D": [
        ("irreducible_rafs", "a lower bound on the number of irreducible RAFs, by sampling", "C F", "network", "exponential"),
        ("autocatalytic_cores", "minimal stoichiometric autocatalytic cores (Blokhuis, Lacoste & Nghe 2020)", "S", "network", "exponential"),
        ("organisations", "number of chemical organisations and the size of the largest (Dittrich & Speroni di Fenizio 2007)", "S", "network", "exponential"),
    ],
    "F": [
        ("steady_states", "number of steady states found by multi-start root finding", "K", "network", "moderate"),
        ("stability", "largest real part of the Jacobian's eigenvalues at steady state, and the stiffness ratio", "K", "network", "moderate"),
        ("oscillation", "whether a simulation settles into sustained oscillation, and its period", "K", "network", "moderate"),
        ("flux_concentration", "Gini coefficient of the steady-state fluxes", "K", "network", "moderate"),
        ("sloppiness", "eigenvalue spread of the Fisher information of the rate constants", "K", "network", "moderate"),
        ("entropy_production", "entropy production at steady state over the reversible pairs", "K", "network", "moderate"),
    ],
    "G": [
        ("attractor_type", "fixed point, cycle or chaos, from the recurrence of the trajectory", "D", "trajectory", "moderate"),
    ],
    "H": [
        ("degeneracy", "structurally different pathways to each species from the food set", "S F", "network", "moderate"),
        ("knockout_tolerance", "share of reactions whose removal leaves the scope and the maxRAF unchanged", "S F", "network", "moderate"),
        ("synthetic_lethal_pairs", "pairs of reactions that back each other up (sampled on large networks)", "S F", "network", "moderate"),
    ],
    "I": [
        ("structure_function_mi", "mutual information between reactant and product structure features", "str", "network", "moderate"),
    ],
}


def table(section: str) -> str:
    rows = ["| measure | meaning | needs | input | cost |", "|---|---|---|---|---|"]
    for m in REGISTRY.values():
        if m.section == section:
            limit = f", up to {m.limit} nodes" if m.limit else ""
            rows.append(f"| `{m.name}` | {m.doc} | {' '.join(sorted(m.needs))} | {m.input} | {m.cost}{limit} |")
    for name, meaning, needs, kind, cost in PLANNED.get(section, []):
        rows.append(f"| *{name}* (planned) | {meaning} | {needs} | {kind} | {cost} |")
    return "\n".join(rows)


BLOCK = re.compile(r"<!-- measures ([A-J]) -->\n.*?<!-- /measures -->", re.S)


def render(text: str) -> str:
    return BLOCK.sub(lambda m: f"<!-- measures {m.group(1)} -->\n{table(m.group(1))}\n<!-- /measures -->", text)


def main(argv: list[str]) -> int:
    text = PAGE.read_text()
    new = render(text)
    if "--check" in argv:
        if new != text:
            print("docs/guide/measures.md is stale: run tools/gen_measures_page.py")
            return 1
        return 0
    PAGE.write_text(new)
    print(f"wrote {PAGE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
