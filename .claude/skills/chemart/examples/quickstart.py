#!/usr/bin/env python3
"""Chemart in five minutes.

    uv run python .claude/skills/chemart/examples/quickstart.py
"""

import json

import chemart

# ---------------------------------------------------------------- discovery
listing = chemart.list_chemistries()
print(f"{len(listing)} chemistries in the catalog\n")

for row in listing[:3]:
    print(f"  {row['id']:<22} {row['fidelity']:<16} {row['summary'][:60]}…")

# ------------------------------------------------------------- one chemistry
described = chemart.describe_chemistry("brusselator")
print("\nbrusselator")
print("  origin :", described["origin"])
print("  book   :", described["book"])
print("  params :", list(described["params"]["properties"]))
print("  provides:", ", ".join(described["provides"]))

# Parameters come with a JSON Schema, so an LLM or a UI can drive them safely.
schema = described["params"]["properties"]["b"]
print("  'b' ->", json.dumps(schema))

# ----------------------------------------------------------------- generate
net = chemart.generate_network("brusselator", seed=1)
print("\n" + net.summary())
print()
print(net.to_text())

# Buffered species are held constant by the model — the Brusselator's A and B
# are reservoirs, so flat trajectories are the model, not a bug.
print("\nbuffered:", net.extras.get("buffered"))

# ------------------------------------------------------- parameters and seeds
wide = chemart.generate_network("brusselator", seed=1, b=5.0)
print(f"\nb=5.0 -> same shape, different kinetics: "
      f"{len(wide.species)} species, {len(wide.reactions)} reactions")

again = chemart.generate_network("brusselator", seed=1)
assert again.to_dict() == net.to_dict(), "same seed must give the same network"
print("same seed reproduces the network exactly")

# Errors name the parameter and its meaning rather than failing obscurely.
try:
    chemart.generate_network("brusselator", nonsense=1)
except ValueError as err:
    print("\nunknown parameter ->", err)

# ------------------------------------------------------------------- structure
ids, R, P = net.matrices()
S = (P - R).toarray()
print(f"\nstoichiometric matrix S = P - R: {S.shape[0]} species x {S.shape[1]} reactions")
print("net change of X per reaction:", S[ids.index("X")].tolist())

# ------------------------------------------------------------------ JSON round trip
data = net.to_dict()
restored = chemart.Network.from_dict(json.loads(json.dumps(data)))
assert restored == net
print("\nrecord round-trips through JSON exactly")

# ------------------------------------------------ a constructive chemistry
# Species are not fixed up front; they are discovered as reactions produce them.
prime = chemart.generate_network("prime-number-chemistry", seed=1)
print(f"\nprime-number-chemistry: {len(prime.species)} species, "
      f"{len(prime.reactions)} reactions, status={prime.status}")

# Always read `status` before interpreting the reaction list. This default is
# an *observed* run: these are the reactions that actually fired in one
# simulated soup, each with a firing count. A reaction's absence says something
# about this run, not about the chemistry.
assert prime.status == "observed"
busiest = max(prime.reactions, key=lambda r: r.count or 0)
print(f"  busiest reaction: {busiest.to_text()}")
print(f"  'complete' would mean the whole defined network or a finished closure;")
print(f"  'truncated' would mean a size budget stopped a closure early.")
print("  See references/comparing-and-analysing.md for computing closures with expand().")
