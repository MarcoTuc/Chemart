"""Machine-tape chemistry: Ikegami & Hashimoto (1995, 1996). Catalog id: ikegami-hashimoto.

Tapes are circular 7-bit strings, stored as the bits read clockwise from the
tape's source (species ``T<2 hex digits>``, site 0 = most significant bit).
Machines are 16-bit words (species ``M<4 hex digits>``): the hex digits are the
T' column, the M' column, the head and the tail. The transition table maps
(tape bit t, machine state m), in the order 00, 01, 10, 11 of the index 2t + m,
to (t', m'), so T' and M' are its two output columns (paper fig. 1).

A machine M reads a tape T if its head matches the tape at some site h (first
match clockwise from the source) and its tail at a different site (first
match clockwise after h). It rewrites the reading frame h .. tail-1 with its
table, starting in state 0 for one half of the population and in state 1 for
the other. The rewritten tape T' is translated into a machine M' by reading 16
bits clockwise from the source (7 bits reused): even bits give the T' column,
odd bits the M' column, then even bits the head and odd bits the tail (fig. 1,
panel 3). So

    M + T -> M + T + M' + T'

method "dynamics" runs the papers' population dynamics (capacity N per
population, a fraction c replaced by products each generation, external noise
flipping bits in the reading frame) and returns the reactions that fired;
method "closure" returns the noise-free reaction network reachable from the
seed machines and tapes.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import product as cartesian

import numpy as np

from chemart.expand import expand
from chemart.helpers.params import apportion
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species

TAPE_LEN = 7
TAPE_MASK = (1 << TAPE_LEN) - 1

# fig. 1 panel 3: tape bit i (from the source) -> machine bit, as machine bit order
# T'(00) T'(01) T'(10) T'(11) M'(00) .. M'(11) head(4 bits) tail(4 bits)
_TRANSLATION = (0, 2, 4, 6, 1, 3, 5, 7, 8, 10, 12, 14, 9, 11, 13, 15)


# --- molecules ------------------------------------------------------------------
def machine_id(m: int) -> str:
    return f"M{m:04x}"


def tape_id(t: int) -> str:
    return f"T{t:02x}"


def tape_bit(t: int, site: int) -> int:
    return (t >> (TAPE_LEN - 1 - site % TAPE_LEN)) & 1


def rotate(t: int, r: int) -> int:
    """The same circular tape read from site r (its source moved r sites clockwise)."""
    r %= TAPE_LEN
    return ((t << r) | (t >> (TAPE_LEN - r))) & TAPE_MASK if r else t


_NECKLACE = [min(rotate(t, r) for r in range(TAPE_LEN)) for t in range(1 << TAPE_LEN)]


def necklace(t: int) -> int:
    """Translation-invariant name: the smallest rotation (the papers' T1, T3, T1b, ...)."""
    return _NECKLACE[t]


def fields(m: int) -> tuple[list[int], list[int], int, int]:
    """(T' outputs, M' outputs, head, tail); outputs indexed by 2 t + m."""
    tp = [(m >> (15 - i)) & 1 for i in range(4)]
    mp = [(m >> (11 - i)) & 1 for i in range(4)]
    return tp, mp, (m >> 4) & 0xF, m & 0xF


def _window(t: int, site: int) -> int:
    return sum(tape_bit(t, site + i) << (3 - i) for i in range(4))


def bind(m: int, t: int) -> tuple[int, int] | None:
    """(head site h, reading-frame length L) or None if the machine cannot read the tape."""
    _, _, head, tail = fields(m)
    h = next((s for s in range(TAPE_LEN) if _window(t, s) == head), None)
    if h is None:
        return None
    for k in range(1, TAPE_LEN):
        if _window(t, h + k) == tail:
            return h, k
    return None


def rewrite_trace(m: int, t: int, state: int) -> list[tuple[int, int]]:
    """[(tape, machine state)] before each step of the reading frame, then the final pair."""
    found = bind(m, t)
    if found is None:
        return []
    h, length = found
    tp, mp, _, _ = fields(m)
    trace = [(t, state)]
    for i in range(length):
        shift = TAPE_LEN - 1 - (h + i) % TAPE_LEN
        k = 2 * ((t >> shift) & 1) + state
        t = (t & ~(1 << shift)) | (tp[k] << shift)
        state = mp[k]
        trace.append((t, state))
    return trace


def _translate(t: int) -> int:
    bits = [tape_bit(t, i) for i in range(16)]
    out = 0
    for i in _TRANSLATION:
        out = (out << 1) | bits[i]
    return out


_MACHINE_OF = [_translate(t) for t in range(1 << TAPE_LEN)]


def translate(t: int) -> int:
    """The machine encoded by a tape, read 16 bits clockwise from its source."""
    return _MACHINE_OF[t]


@lru_cache(maxsize=None)
def react(m: int, t: int, state: int) -> tuple[int, int, int, int, int] | None:
    """(M', T', head site, frame length L, rewritten bits w), or None if M cannot read T."""
    found = bind(m, t)
    if found is None:
        return None
    h, length = found
    new, _ = rewrite_trace(m, t, state)[-1]
    return translate(new), new, h, length, bin(new ^ t).count("1")


def frame_mask(h: int, flips: int, length: int) -> int:
    """Tape mask flipping the frame positions i (0 = head site) whose bit i is set in `flips`."""
    out = 0
    for i in range(length):
        if flips >> i & 1:
            out |= 1 << (TAPE_LEN - 1 - (h + i) % TAPE_LEN)
    return out


@lru_cache(maxsize=None)
def _flip_distribution(length: int, mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Non-empty flip sets of a frame and their probabilities given at least one error (eq. 7)."""
    sets = np.arange(1, 1 << length)
    ones = np.array([bin(s).count("1") for s in sets])
    p = mu ** ones * (1.0 - mu) ** (length - ones)
    return sets, p / p.sum()


# --- parameters -----------------------------------------------------------------
def _parse(values, bits: int, what: str, example: str) -> list[int]:
    out = []
    for v in values:
        try:
            x = int(v, 16) if isinstance(v, str) else (v if isinstance(v, int) and not isinstance(v, bool) else -1)
        except ValueError:
            x = -1
        if not 0 <= x < 1 << bits:
            raise ValueError(f"{what} must be {bits}-bit values as hex strings (e.g. {example!r}) or ints, got {v!r}")
        out.append(x)
    return out


def _seeds(p, rng) -> tuple[list[int], list[int]]:
    machines = _parse(p.machines, 16, "machines", "1002")
    tapes = _parse(p.tapes, TAPE_LEN, "tapes", "01")
    if not machines:
        machines = [int(x) for x in rng.integers(0, 1 << 16, size=p.n_machines)]
    if not tapes:
        tapes = [int(x) for x in rng.integers(0, 1 << TAPE_LEN, size=p.n_tapes)]
    return list(dict.fromkeys(machines)), list(dict.fromkeys(tapes))


def generate(p, rng):
    machines, tapes = _seeds(p, rng)
    if p.method == "closure":
        return _closure(p, machines, tapes)
    if p.noise_off > p.generations:
        raise ValueError(f"noise_off ({p.noise_off}) must be -1 (never) or at most generations ({p.generations})")
    return _dynamics(p, rng, machines, tapes)


# --- network pieces -------------------------------------------------------------
def _rate(c: float, states: int, length: int) -> dict:
    """Eq. 4 numerator: c f_M f_T, half of it per initial machine state."""
    return {"law": "mass-action", "k": c * states / 2, "frame_length": length}


def _reaction(m, t, m2, t2, rate, count=None) -> Reaction:
    return Reaction.of([machine_id(m), tape_id(t)], [machine_id(m), tape_id(t), machine_id(m2), tape_id(t2)],
                       rate=rate, count=count)


def _species(machines, tapes) -> list[Species]:
    return ([Species(machine_id(m), f"{m:016b}") for m in sorted(set(machines))]
            + [Species(tape_id(t), f"{t:07b}") for t in sorted(set(tapes))])


def _paper_names(tapes) -> dict[str, str]:
    return {tape_id(t): f"T{necklace(t):x}" for t in sorted(set(tapes))}


# --- closure --------------------------------------------------------------------
def _closure(p, machines, tapes):
    def full(a, b):
        if not (isinstance(a, tuple) and a[0] == "M" and b[0] == "T"):
            return None
        out = [react(a[1], b[1], s) for s in (0, 1)]
        if out[0] is None:
            return None
        return (a, b, *[x for r in out for x in (("M", r[0]), ("T", r[1]))])

    seed = [("M", m) for m in machines] + [("T", t) for t in tapes]
    found, _, status = expand(full, seed, arity=2, max_species=p.max_species, ordered=True)
    ms = [x for k, x in found if k == "M"]
    ts = [x for k, x in found if k == "T"]
    known_m, known_t = set(ms), set(ts)
    reactions = []
    for m, t in cartesian(ms, ts):
        products = Counter()
        length = 0
        for s in (0, 1):
            r = react(m, t, s)
            if r is not None and r[0] in known_m and r[1] in known_t:
                products[r[0], r[1]] += 1
                length = r[3]
        for (m2, t2), n in products.items():
            reactions.append(_reaction(m, t, m2, t2, _rate(p.c, n, length)))
    return Network(
        species=_species(ms, ts),
        reactions=reactions,
        status=status,
        outflow=CONSTANT_TOTAL,
        extras={"seed": [machine_id(m) for m in machines] + [tape_id(t) for t in tapes],
                "paper_names": _paper_names(ts)},
    )


# --- population dynamics ---------------------------------------------------------
def _stochastic_round(x: float, rng) -> int:
    n = int(x)
    return n + int(rng.random() < x - n)


def _dynamics(p, rng, machines, tapes):
    N, c = p.N, p.c
    mpop = {m: n for m, n in zip(machines, apportion(N, [1.0] * len(machines)))} if machines else {}
    tpop: dict[int, int] = {}
    for t, n in (zip(tapes, apportion(N, [1.0] * len(tapes))) if tapes else ()):
        # rotations of one circular tape share a source: merge them into the first
        key = next((u for u in tpop if necklace(u) == necklace(t)), t)
        tpop[key] = tpop.get(key, 0) + n
    mpop = {m: n for m, n in mpop.items() if n > 0}
    tpop = {t: n for t, n in tpop.items() if n > 0}
    initial = {**{machine_id(m): n for m, n in mpop.items()}, **{tape_id(t): n for t, n in tpop.items()}}

    seen_m, seen_t = set(mpop), set(tpop)
    fired: dict[tuple, dict] = {}
    analysis = {k: [] for k in ("distinct_machines", "distinct_tapes", "active_mutation", "reading_length")}

    def record(key, gen, exact, states=0, length=0):
        entry = fired.setdefault(key, {"gens": set(), "states": 0, "length": length, "exact": False})
        entry["gens"].add(gen)
        if exact:
            entry["exact"] = True
            entry["states"] = max(entry["states"], states)
            entry["length"] = length

    for gen in range(p.generations):
        analysis["distinct_machines"].append(len(mpop))
        analysis["distinct_tapes"].append(len(tpop))
        mu = p.noise if p.noise_off < 0 or gen < p.noise_off else 0.0
        registry = {necklace(t): t for t in tpop}

        def place(t2):
            k = necklace(t2)
            if k not in registry:
                registry[k] = rotate(t2, int(rng.integers(TAPE_LEN))) if p.source == "random" else t2
            return registry[k]

        pairs = []
        weight = 0.0
        am_num = len_num = 0.0
        for m, nm in mpop.items():
            for t, nt in tpop.items():
                r0 = react(m, t, 0)
                if r0 is None:
                    continue
                w = nm * nt
                pairs.append((m, t, w, (r0, react(m, t, 1))))
                weight += w
                am_num += w * sum(r[4] for r in (r0, react(m, t, 1))) / (2 * r0[3])
                len_num += w * r0[3]
        analysis["active_mutation"].append(am_num / weight if weight else 0.0)
        analysis["reading_length"].append(len_num / weight if weight else 0.0)

        new_m: Counter = Counter()
        new_t: Counter = Counter()
        for m, t, w, outs in pairs:
            flux = c * N * w / weight / 2          # new machines (and tapes) per initial state
            exact = Counter()
            for r in outs:
                m2, t2, h, length, _ = r
                eps = 1.0 - (1.0 - mu) ** length   # eq. 7
                t2p = place(t2)
                new_m[m2] += flux * (1.0 - eps)
                new_t[t2p] += flux * (1.0 - eps)
                exact[m2, t2p] += 1
                n_mut = _stochastic_round(flux * eps, rng) if eps > 0 else 0
                if n_mut:
                    sets, probs = _flip_distribution(length, mu)
                    for flips in rng.choice(sets, size=n_mut, p=probs):
                        tm = t2 ^ frame_mask(h, int(flips), length)
                        mm, tmp = translate(tm), place(tm)
                        new_m[mm] += 1
                        new_t[tmp] += 1
                        record((m, t, mm, tmp), gen, exact=False)
            for (m2, t2p), n in exact.items():
                record((m, t, m2, t2p), gen, exact=True, states=n, length=outs[0][3])

        def update(pop, new):
            out = {}
            for x in set(pop) | set(new):
                n = int((1.0 - c) * pop.get(x, 0) + new.get(x, 0.0) + 1e-9)
                if n > 0:
                    out[x] = n
            return out

        mpop, tpop = update(mpop, new_m), update(tpop, new_t)
        seen_m.update(mpop)
        seen_t.update(tpop)

    # products that never reached one copy still took part in a reaction
    seen_m.update(k[2] for k in fired)
    seen_t.update(k[3] for k in fired)
    reactions, noise_only = [], []
    for (m, t, m2, t2), e in fired.items():
        rate = _rate(c, e["states"], e["length"]) if e["exact"] else None
        if not e["exact"]:
            noise_only.append(len(reactions))
        reactions.append(_reaction(m, t, m2, t2, rate, count=len(e["gens"])))
    analysis["distinct_machines"].append(len(mpop))
    analysis["distinct_tapes"].append(len(tpop))
    return Network(
        species=_species(seen_m, seen_t),
        reactions=reactions,
        status="observed",
        initial_state=initial,
        outflow=CONSTANT_TOTAL,
        extras={
            "analysis": analysis,
            "final_state": {**{machine_id(m): n for m, n in sorted(mpop.items(), key=lambda kv: -kv[1])},
                            **{tape_id(t): n for t, n in sorted(tpop.items(), key=lambda kv: -kv[1])}},
            "noise_induced": noise_only,
            "paper_names": _paper_names(seen_t),
        },
    )
