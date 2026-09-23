#!/usr/bin/env python3
"""Simulate a catalogued chemistry, evolve a gas, and check published behaviours.

    uv run python .claude/skills/chemart/examples/dynamics.py

Shows what trips people up: finding which chemistries can be simulated as they
come, that buffered species are meant to be flat, how to give rates to a
network that has none, and that a Turing gas is evolved, not simulated.
"""

import numpy as np

import chemart
from chemart import measures, simulate

# ------------------------------------------------- which ones carry kinetics
rated = [row["id"] for row in chemart.list_chemistries()
         if "rate-constants" in chemart.describe_chemistry(row["id"])["provides"]]
print(f"{len(rated)} entries carry rate constants; e.g. {rated[:6]}\n")

# --------------------------------------------------------- the Brusselator
net = chemart.generate_network("brusselator", seed=1)
traj = simulate.ode(net, t_end=40.0, points=800)
ids, t, X = traj.array()
x = X[:, ids.index("X")]

buffered = set(net.extras.get("buffered", []))
print(f"brusselator: buffered = {sorted(buffered)} (reservoirs, held constant)")

# The published behaviour: with b > 1 + a^2 the steady state (a, b/a) is
# unstable and the system settles into a limit cycle. Check that X keeps
# oscillating instead of relaxing, by comparing the swing late in the run.
late = t > t[-1] / 2
swing = x[late].max() - x[late].min()
print(f"  X still swings by {swing:.3f} over the second half -> sustained oscillation")
assert swing > 0.1, "expected a limit cycle at these parameters"

for s in buffered:
    column = X[:, ids.index(s)]
    assert np.allclose(column, column[0]), s
print("  buffered species stayed flat, as the model intends")

# One stochastic path of the same network: counts = amount x volume.
path = simulate.ssa(net, t_end=40.0, volume=100, seed=0)
print(f"  one SSA path at volume 100: {path.settings['events']} events")

# ------------------------------------------------- a network without kinetics
kauffman = chemart.generate_network("kauffman-autocatalytic-sets", seed=1)
try:
    simulate.ode(kauffman, t_end=1.0, x0=1.0)
except simulate.NotSimulable as err:
    print(f"\nkauffman-autocatalytic-sets: {err}")
traj = simulate.ode(kauffman, t_end=10.0, x0=1.0, seed=0,
                    rates={"dist": "lognormal", "mean": 0, "sigma": 1})
print(f"  with lognormal rates drawn for its {len(kauffman.reactions)} reactions it runs: "
      f"{len(traj.frames)} frames")

# ------------------------------------------------------- a gas is evolved
# The prime-number chemistry's published behaviour: composite numbers keep
# being divided away, so the soup turns prime. Its process runs in collisions.
run = chemart.evolve("prime-number-chemistry", seed=1)
fraction = run.series("prime_fraction")
print(f"\nprime-number-chemistry, {len(run.frames)} frames in {run.clock}: "
      f"prime fraction {fraction[0]:.2f} -> {fraction[-1]:.2f}")
assert fraction[-1] > 0.9, "expected the soup to turn prime"

# Measures in evolutionary time: population measures read each frame's soup.
diversity = measures.over(run, ["richness", "shannon"])
print(f"  richness {diversity['richness'][0]} -> {diversity['richness'][-1]}, "
      f"Shannon {diversity['shannon'][0]:.2f} -> {diversity['shannon'][-1]:.2f} nats")
