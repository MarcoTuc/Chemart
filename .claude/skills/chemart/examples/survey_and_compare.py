#!/usr/bin/env python3
"""Compare chemistries structurally — the point of a common record.

    uv run python .claude/skills/chemart/examples/survey_and_compare.py
"""

import numpy as np

import chemart
from chemart.catalog import load

# ------------------------------------------------------- catalog-level survey
entries = load()
with_kinetics = [c for c in entries if "rate-constants" in c.provides]
with_energy = [c for c in entries if "energies" in c.provides]
both = [c for c in entries if {"rate-constants", "energies"} <= set(c.provides)]

print(f"{len(entries)} entries: {len(with_kinetics)} carry kinetics, "
      f"{len(with_energy)} carry energetics, {len(both)} carry both")
print("  both:", ", ".join(sorted(c.id for c in both)), "\n")

# --------------------------------------------------- structural profiles
def profile(cid, seed=1):
    net = chemart.generate_network(cid, seed=seed)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    rank = int(np.linalg.matrix_rank(S)) if S.size else 0
    return {
        "id": cid,
        "species": len(ids),
        "reactions": len(net.reactions),
        "status": net.status,
        "rank": rank,
        "cons_laws": len(ids) - rank,          # upper bound on independent laws
        "catalysed": sum(1 for r in net.reactions if r.catalysts),
    }

print(f"{'chemistry':<26} {'S':>5} {'R':>6} {'rank':>5} {'laws':>5} {'cat':>5}  status")
print("-" * 70)
for cid in ["dimerization", "brusselator", "oregonator", "michaelis-menten",
            "matrix-chemistry", "prime-number-chemistry", "gard"]:
    p = profile(cid)
    print(f"{p['id']:<26} {p['species']:>5} {p['reactions']:>6} {p['rank']:>5} "
          f"{p['cons_laws']:>5} {p['catalysed']:>5}  {p['status']}")

print("\nNote: rank and 'laws' are only comparable within a status. An 'observed'")
print("network is a sample of one run, not the chemistry's definition.")
print("'laws' is the upper bound species - rank(S), not the number a chemistry")
print("declares: most entries declare none, and a few declare fewer than the bound.")

# ------------------------------------------------------- conservation laws
# Use a chemistry that actually declares them: 31 of 98 do, and checking a
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
