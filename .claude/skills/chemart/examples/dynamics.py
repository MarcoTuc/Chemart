#!/usr/bin/env python3
"""Simulate a catalogued chemistry and check a published behaviour.

    uv run python .claude/skills/chemart/examples/dynamics.py

Shows the two things that trip people up: finding which chemistries are
integrable at all, and remembering that buffered species are meant to be flat.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import chemart                                    # noqa: E402
from simulate import NotIntegrable, integrate     # noqa: E402

# ------------------------------------------------- which ones can be simulated
integrable = [row["id"] for row in chemart.list_chemistries()
              if "rate-constants" in chemart.describe_chemistry(row["id"])["provides"]]
print(f"{len(integrable)} of 98 entries claim rate constants; e.g. {integrable[:6]}\n")

# --------------------------------------------------------- the Brusselator
net = chemart.generate_network("brusselator", seed=1)
traj, t = integrate(net, t_end=40.0, points=800)

buffered = set(net.extras.get("buffered", []))
print(f"brusselator: buffered = {sorted(buffered)} (reservoirs, held constant)")

x, y = traj["X"], traj["Y"]

# The published behaviour: with b > 1 + a^2 the steady state (a, b/a) is
# unstable and the system settles into a limit cycle. Check that X keeps
# oscillating instead of relaxing, by comparing the swing late in the run.
late = slice(len(t) // 2, None)
swing = x[late].max() - x[late].min()
print(f"  X still swings by {swing:.3f} over the second half -> sustained oscillation")
assert swing > 0.1, "expected a limit cycle at these parameters"

# Buffered species really are constant.
for s in buffered:
    assert np.allclose(traj[s], traj[s][0]), s
print(f"  buffered species stayed flat, as the model intends")

# ------------------------------------------- a chemistry that cannot be integrated
try:
    integrate(chemart.generate_network("tierra", seed=1), t_end=1.0)
except NotIntegrable as err:
    print(f"\ntierra: {err}")

# ------------------------------------------------- stochastic conversion
# Mass-action k is a deterministic rate; Gillespie needs a stochastic c.
from chemart.kinetics import k_to_c                # noqa: E402

k = 1.0
print(f"\nk -> c at volume 1e-15, heterodimer A + B : {k_to_c(k, {'A': 1, 'B': 1}, 1e-15):.4g}")
print(f"k -> c at volume 1e-15, homodimer  2 A     : {k_to_c(k, {'A': 2}, 1e-15):.4g}"
      "   (the factor 2 is the combinatorial correction)")
