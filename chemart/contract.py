"""The contract every generator must satisfy, as reusable checks.

The test suite runs these over the built-in catalog; `chemart check` and
`chemart.hub.push_generator` run them on a chemistry before it is shared, so
what reaches the hub meets the same bar as what ships with Chemart.

    problems(entry, generate)              -> list of problems; [] is a pass
    problems(entry, generate, evolve=f)    also checks the evolve face
"""

from __future__ import annotations

import json
import time
from collections import Counter
from dataclasses import asdict
from typing import Callable

from chemart import api
from chemart.catalog import Chemistry
from chemart.network import Network

#: Default parameters must run in seconds: paper-scale values go in `range`.
TIME_LIMIT = 5.0


def problems(entry: Chemistry, generate: Callable, *, evolve: Callable | None = None,
             seed: int = 1, time_limit: float = TIME_LIMIT) -> list[str]:
    """Check `generate(p, rng)`, and `evolve(p, rng)` if given, against the
    catalog `entry`.

    Runs each face twice at default parameters: it must finish within
    `time_limit`, return a network that is plain JSON data (it round-trips
    exactly), give the same result for the same seed, and contain nothing the
    entry does not claim in `provides`. A given network must not depend on the
    seed. The evolve face must also yield frames that start at the initial
    state, never go back in time, round-trip through JSON, and add up to the
    firing counts of the network it returns.
    """
    out = _generate_problems(entry, generate, seed, time_limit)
    if evolve is not None:
        out += _evolve_problems(entry, evolve, seed, time_limit)
    return out


def _generate_problems(entry: Chemistry, generate: Callable, seed: int, time_limit: float) -> list[str]:
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

    if entry.type == "given":
        try:
            other = api.run_generator(entry, generate, seed + 1, {})
        except Exception as err:
            out.append(f"run with seed={seed + 1} failed: {type(err).__name__}: {err}")
        else:
            if _topology(other) != _topology(net):
                out.append("a given network changed its topology with the seed; "
                           "a network drawn at random is a generator (type: generator)")
    return out


def _topology(net: Network) -> set:
    return {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in net.reactions}


def _run(entry, evolve, seed):
    run = api.run_evolver(entry, evolve, seed, {})
    frames = []
    while True:
        try:
            frames.append(next(run))
        except StopIteration as stop:
            return frames, stop.value


def _evolve_problems(entry: Chemistry, evolve: Callable, seed: int, time_limit: float) -> list[str]:
    out: list[str] = []
    if not entry.clock:
        out.append("the module has an evolve face but the entry has no clock "
                   "(the unit of its time: collisions, epochs, ...)")
    start = time.perf_counter()
    try:
        frames, net = _run(entry, evolve, seed)
    except Exception as err:
        return out + [f"evolve with default parameters failed: {type(err).__name__}: {err}"]
    elapsed = time.perf_counter() - start
    if elapsed > time_limit:
        out.append(f"evolve at default parameters took {elapsed:.1f} s; it must run in under "
                   f"{time_limit:g} s")
    if not frames:
        return out + ["evolve yielded no frames; the first frame is the initial state"]
    if frames[0].fired:
        out.append("the first frame must be the initial state, with nothing fired yet")
    times = [f.t for f in frames]
    if any(b < a for a, b in zip(times, times[1:])):
        out.append("frame times go backwards")
    data = [asdict(f) for f in frames]
    if json.loads(json.dumps(data)) != data:
        out.append("frames are not plain JSON data (numpy scalars or tuples?)")
    try:
        again, _ = _run(entry, evolve, seed)
    except Exception as err:
        out.append(f"second evolve with the same seed failed: {type(err).__name__}: {err}")
    else:
        if [asdict(f) for f in again] != data:
            out.append(f"seed={seed} gave two different runs; draw all randomness from rng")
    fired = Counter()
    for f in frames:
        for lhs, rhs, n in f.fired:
            fired[(frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))] += n
    counted = Counter({(frozenset(r.reactants.items()), frozenset(r.products.items())): r.count
                       for r in net.reactions if r.count})
    if fired != counted:
        out.append("the reactions fired in the frames do not add up to the returned network's counts")
    known = {s.id for s in net.species}
    missing = sorted({s for f in frames for s in f.state} - known)[:5]
    if missing:
        out.append(f"frame states name species the network does not have: {missing}")
    unclaimed = sorted(set(net.provides) - set(entry.provides))
    if unclaimed:
        out.append(f"the evolved network contains {unclaimed} but the entry does not claim it")
    return out
