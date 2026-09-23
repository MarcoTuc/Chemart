#!/usr/bin/env python3
"""Compare chemistries structurally — the point of a common record.

    uv run python .claude/skills/chemart/examples/survey_and_compare.py
"""

import numpy as np

import chemart
from chemart.catalog import active, load

# ------------------------------------------------------- catalog-level survey
entries = active(load())          # the catalog, without the archive
with_kinetics = [c for c in entries if "rate-constants" in c.provides]
with_energy = [c for c in entries if "energies" in c.provides]
both = [c for c in entries if {"rate-constants", "energies"} <= set(c.provides)]

print(f"{len(entries)} entries: {len(with_kinetics)} carry kinetics, "
      f"{len(with_energy)} carry energetics, {len(both)} carry both")
print("  both:", ", ".join(sorted(c.id for c in both)), "\n")

# --------------------------------------------------- structural profiles
# chemart.measure computes registered measures; one that does not apply is
# left out (shown as "-"), and measures.applicable says why.
names = ["n_species", "n_reactions", "stoichiometric_rank", "conservation_laws",
         "deficiency", "catalysed_fraction"]
print(f"{'chemistry':<24} {'type':<10} {'S':>5} {'R':>6} {'rank':>5} {'laws':>5} "
      f"{'defic':>6} {'cat':>5}  status")
print("-" * 84)
for cid in ["brusselator", "oregonator", "michaelis-menten",
            "matrix-chemistry", "prime-number-chemistry", "gard"]:
    net = chemart.generate_network(cid, seed=1)
    m = chemart.measure(net, names)
    cells = [m.get(n, "-") for n in names]
    cells[-1] = f"{cells[-1]:.2f}" if isinstance(cells[-1], float) else cells[-1]
    kind = chemart.describe_chemistry(cid)["type"]
    print(f"{cid:<24} {kind:<10} {cells[0]:>5} {cells[1]:>6} {cells[2]:>5} {cells[3]:>5} "
          f"{cells[4]:>6} {cells[5]:>5}  {net.status}")

print("\nNote: these are only comparable within a status. An 'observed' network is")
print("a sample of one run, not the chemistry's definition. 'laws' is the number")
print("of independent conservation laws of S (species - rank), not the number a")
print("chemistry declares: most entries declare none, and a few declare fewer.")

# ------------------------------------------------------- conservation laws
# Use a chemistry that actually declares them: many do not, and checking a
# declared law against S is the fastest way to catch a stoichiometry bug.
net = chemart.generate_network("chameleon", seed=1)
ids, R, P = net.matrices()
S = (P - R).toarray()

declared = net.extras.get("conservation", [])
print(f"\nchameleon declares {len(declared)} conservation law(s)")
for i, law in enumerate(declared):
    # Read the payload defensively: `vector` is a dict for most entries but a
    # list for some; `name` is not guaranteed; and a `modulus` makes the law
    # modular, so a plain zero test would flag a correct law as violated.
    v = law["vector"]
    m = np.array([v.get(s, 0) for s in ids]) if isinstance(v, dict) else np.array(v)
    residual = (S.astype(int)).T @ m
    modulus = law.get("modulus")
    holds = (residual % modulus == 0).all() if modulus else not residual.any()
    name = law.get("name", f"law {i}")
    extra = f" (mod {modulus})" if modulus else ""
    print(f"  {name:<24}{extra:<9} {'holds' if holds else f'VIOLATED {residual}'}")

# --------------------------------------------------------- catalysis
net = chemart.generate_network("brusselator", seed=1)
for r in net.reactions:
    if r.catalysts:
        print(f"\ncatalytic reaction: {r.to_text()}")
        print(f"  catalysts (survive the reaction): {r.catalysts}")

autocat = [r for r in net.reactions
           if any(r.products.get(s, 0) > n for s, n in r.reactants.items())]
for r in autocat:
    print(f"autocatalytic: {r.to_text()}  (a reactant comes out amplified)")
