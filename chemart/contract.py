"""The contract every generator must satisfy, as reusable checks.

The test suite runs these over the built-in catalog; `chemart check` and
`chemart.hub.push_generator` run them on a chemistry before it is shared, so
what reaches the hub meets the same bar as what ships with Chemart.

    problems(entry, generate)   -> list of human-readable problems; [] is a pass
"""

from __future__ import annotations

import json
import time
from typing import Callable

from chemart import api
from chemart.catalog import Chemistry
from chemart.network import Network

#: Default parameters must run in seconds: paper-scale values go in `range`.
TIME_LIMIT = 5.0


def problems(entry: Chemistry, generate: Callable, *, seed: int = 1,
             time_limit: float = TIME_LIMIT) -> list[str]:
    """Check `generate(p, rng)` against its catalog `entry`.

    Runs the generator twice at default parameters: it must finish within
    `time_limit`, return a network that is plain JSON data (it round-trips
    exactly), give the same network for the same seed, and contain nothing
    the entry does not claim in `provides`.
    """
    out: list[str] = []
    start = time.perf_counter()
    try:
        net = api.run_generator(entry, generate, seed, {})
    except Exception as err:  # the point is to report, not to crash
        return [f"generation with default parameters failed: {type(err).__name__}: {err}"]
    elapsed = time.perf_counter() - start
    if elapsed > time_limit:
        out.append(f"default parameters took {elapsed:.1f} s; they must run in under "
                   f"{time_limit:g} s (put paper-scale values in the param's range)")

    try:
        restored = Network.from_dict(json.loads(json.dumps(net.to_dict())))
    except (TypeError, ValueError) as err:
        out.append(f"network is not plain JSON data: {err}")
    else:
        if restored != net:
            out.append("network does not round-trip through JSON exactly "
                       "(numpy scalars or tuples in the data?)")

    try:
        again = api.run_generator(entry, generate, seed, {})
    except Exception as err:
        out.append(f"second run with the same seed failed: {type(err).__name__}: {err}")
    else:
        if again.to_dict() != net.to_dict():
            out.append(f"seed={seed} gave two different networks; draw all randomness from rng")

    unclaimed = sorted(set(net.provides) - set(entry.provides))
    if unclaimed:
        out.append(f"network contains {unclaimed} but the entry does not claim it in provides")
    return out
