"""Coreworld: Rasmussen, Knudsen, Feldberg & Hindsholm (1990), book 10.6.3.

Catalog id: coreworld.

The simulator VENUS (LANL report LA-UR-89-2618, the conference version of
Physica D 42:111): a one-dimensional cyclic core of `core_size` words, each
one of the ten Redcode instructions of Table 1 (DAT JMP JMZ JMN DJN ADD SUB
MOV CMP SPL) with two operands in the addressing modes # $ @ <.

- Execution pointers sit in a queue of at most `queue_len` places. A pointer
  executes the word it points to and moves to the next address unless the
  instruction sends it elsewhere; DAT kills it, SPL splits it.
- Computational resources (eqs. 1-3): every address holds r <= r_max < 1
  exec. A pointer executes only if the 2 R_res + 1 addresses around it hold
  at least one exec together; the execution then removes exactly one exec
  from them, proportionally (factor (S - 1)/S). After the whole queue, every
  address is refilled by Delta r, capped at r_max.
- Locality: every relative address is folded into [-R_opr, R_opr] around
  the executing word (the operation radius).
- Noise: a MOV writes a random word with probability P_mut (uniform over the
  ten instructions, random operands); with probability P_point per core
  update a new pointer appears at a random address.
- Venus I (the paper's simulator): all pointers of an update read the same
  core and their writes are applied afterwards; on a conflict the pointer
  later in the queue wins. Venus II (book): pointers are executed one after
  the other on the live core.

Organisms are not defined in Coreworld (book), so the network is the set of
observed *write events*: the executing word, the word it copies or adds
from (MOV, ADD, SUB with a non-immediate A operand) and every word it
changes are the reactants; the same cells after the write are the products.
With species = instruction (default), MOV $0,$1 overwriting a DAT is
`MOV + DAT -> 2 MOV`, and a copy loop duplicating an SPL over an ADD is
`MOV + SPL + ADD -> MOV + 2 SPL`. Events that leave the species multiset
unchanged (e.g. ADD changing a number) are elastic and counted only in
extras.analysis. Every event conserves the number of cells.
"""

from __future__ import annotations

import re
from collections import Counter

from chemart.network import Network, Reaction, Species

OPS = ("DAT", "JMP", "JMZ", "JMN", "DJN", "ADD", "SUB", "MOV", "CMP", "SPL")
DAT, JMP, JMZ, JMN, DJN, ADD, SUB, MOV, CMP, SPL = range(10)
MODES = "#$@<"
IMM, DIR, IND, DEC = range(4)
EPS = 1e-9

#: Appendix: MICE (Chip Wendell), the 8-instruction self-replicator used as seed; pointer on the MOV.
MICE = ["DAT #7", "MOV #7, $-1", "MOV @-2, <5", "DJN $-1, $-3", "SPL @3",
        "ADD #417, $2", "JMZ $-5, $-6", "DAT #714"]
SEEDS = {
    "mice": (MICE, 1),
    "jmp": (["JMP $0"], 0),                                     # "an active JMP($0)"
    "imp": (["MOV $0, $1"], 0),                                 # fig. 3(a)
    "copy-loop": (["MOV #4, $10", "MOV $-1, @3", "ADD #1, $2", "JMP $-2", "DAT #2"], 1),   # fig. 3(b)
    "none": ([], None),
}

_WORD = re.compile(r"^\s*([A-Za-z]{3})\s*(.*?)\s*$")
_OPERAND = re.compile(r"^([#$@<]?)\s*([+-]?\d+)$")


# ----------------------------------------------------------------------------
# words

def parse(text: str, core_size: int) -> tuple:
    """'MOV @-2, <5' -> (op, amode, a, bmode, b), values modulo core_size.

    A missing mode is $. One operand goes to B for DAT and SPL (Table 1:
    "DAT B", "SPL B"), with A = #0; to A for the others, with B = $0.
    """
    m = _WORD.match(text)
    if not m or m.group(1).upper() not in OPS:
        raise ValueError(f"cannot parse Redcode word {text!r}: opcode must be one of {OPS}")
    op = OPS.index(m.group(1).upper())
    parts = [s.strip() for s in m.group(2).split(",")] if m.group(2) else []
    if len(parts) > 2:
        raise ValueError(f"cannot parse Redcode word {text!r}: at most two operands")
    fields = []
    for s in parts:
        o = _OPERAND.match(s)
        if not o:
            raise ValueError(f"cannot parse operand {s!r} in {text!r}")
        fields.append((MODES.index(o.group(1) or "$"), int(o.group(2)) % core_size))
    if not fields:
        fields = [(IMM, 0), (IMM, 0)] if op == DAT else [(DIR, 0), (DIR, 0)]
    elif len(fields) == 1:
        fields = [(IMM, 0), fields[0]] if op in (DAT, SPL) else [fields[0], (DIR, 0)]
    (am, a), (bm, b) = fields
    return (op, am, a, bm, b)


def signed(v: int, core_size: int) -> int:
    return v - core_size if v > core_size // 2 else v


def text(word: tuple, core_size: int) -> str:
    op, am, a, bm, b = word
    return f"{OPS[op]} {MODES[am]}{signed(a, core_size)}, {MODES[bm]}{signed(b, core_size)}"


def word_id(word: tuple, core_size: int) -> str:
    op, am, a, bm, b = word
    return f"{OPS[op]}_{MODES[am]}{signed(a, core_size)}_{MODES[bm]}{signed(b, core_size)}"


# ----------------------------------------------------------------------------
# the machine

class Venus:
    """The core, the pointer queue, the resources and the recorder of write events."""

    def __init__(self, words, rand, *, mode="venus-i", queue_len=220, operation_radius=None,
                 resource_radius=3, resource_influx=0.5, resource_max=0.5, mutation_rate=0.0,
                 injection_rate=0.0, resources=True, record=True):
        self.core = list(words)
        self.M = M = len(self.core)
        self.rand = rand
        self.parallel = mode == "venus-i"
        self.L = queue_len
        self.R = operation_radius if operation_radius is not None and 2 * operation_radius + 1 < M else None
        self.resources = resources
        self.rmax, self.dr = resource_max, resource_influx
        self.r = [resource_max] * M
        self.low: set[int] = set()
        span = min(2 * resource_radius + 1, M)
        self.offsets = list(range(-resource_radius, -resource_radius + span))
        self.pmut, self.ppoint = mutation_rate, injection_rate
        self.record = record
        self.queue: list[int] = []
        self.time = 0
        self.events: dict[tuple, list] = {}
        self.stats = Counter()
        self.executed = Counter()

    # -- randomness -----------------------------------------------------------
    def random_word(self) -> tuple:
        u = self.rand
        M = self.M
        return (int(u() * 10), int(u() * 4), int(u() * M), int(u() * 4), int(u() * M))

    # -- addressing -------------------------------------------------------------
    def fold(self, off: int) -> int:
        """A relative address as a signed offset, folded into the operation radius."""
        M = self.M
        off %= M
        if off > M // 2:
            off -= M
        R = self.R
        if R is not None and not -R <= off <= R:
            off = (off + R) % (2 * R + 1) - R
        return off

    def operand(self, x: int, mode: int, val: int, w: dict) -> int:
        """Absolute address of an operand (ICWS'88 evaluation; # means the executing word)."""
        if mode == IMM:
            return x
        off = self.fold(val)
        if mode == DIR:
            return (x + off) % self.M
        p = (x + off) % self.M
        word = w.get(p) or self.core[p]
        b = word[4]
        if mode == DEC:
            b = (b - 1) % self.M
            w[p] = (word[0], word[1], word[2], word[3], b)
        return (x + self.fold(off + b)) % self.M

    def execute(self, x: int, w: dict, mutate: bool = True):
        """Run the word at x against the core plus the overlay w (which receives the writes).

        Returns (successor addresses, source address or None).
        """
        M = self.M
        core = self.core
        op, am, a, bm, b = w.get(x) or core[x]
        if op == DAT:
            return (), None
        pa = self.operand(x, am, a, w)
        A = w.get(pa) or core[pa]
        pb = self.operand(x, bm, b, w)
        B = w.get(pb) or core[pb]
        nxt = (x + 1) % M
        if op == MOV:
            if am == IMM:
                new, src = (B[0], B[1], B[2], B[3], a), None
            else:
                new, src = A, pa
            if mutate and self.pmut and self.rand() < self.pmut:
                new = self.random_word()
                self.stats["mutations"] += 1
            w[pb] = new
            return (nxt,), src
        if op == ADD or op == SUB:
            s = 1 if op == ADD else -1
            if am == IMM:
                w[pb] = (B[0], B[1], B[2], B[3], (B[4] + s * a) % M)
                return (nxt,), None
            w[pb] = (B[0], B[1], (B[2] + s * A[2]) % M, B[3], (B[4] + s * A[4]) % M)
            return (nxt,), pa
        if op == JMP:
            return (pa,), None
        if op == JMZ:
            return ((pa,) if B[4] == 0 else (nxt,)), None
        if op == JMN:
            return ((pa,) if B[4] != 0 else (nxt,)), None
        if op == DJN:
            v = (B[4] - 1) % M
            w[pb] = (B[0], B[1], B[2], B[3], v)
            return ((pa,) if v != 0 else (nxt,)), None
        if op == CMP:
            equal = (a == B[4]) if am == IMM else (A == B)
            return ((nxt,) if equal else ((x + 2) % M,)), None     # Table 1: skip if unequal
        return (nxt, pb), None                                     # SPL: next statement, then B

    # -- resources (eqs. 1-3) ---------------------------------------------------
    def consume(self, x: int) -> bool:
        if not self.resources:
            return True
        r, M = self.r, self.M
        idx = [(x + k) % M for k in self.offsets]
        s = 0.0
        for i in idx:
            s += r[i]
        if s < 1.0 - EPS:
            return False
        f = (s - 1.0) / s if s > 1.0 else 0.0
        for i in idx:
            r[i] *= f
        self.low.update(idx)
        return True

    def renew(self) -> None:
        if not self.resources:
            return
        r, rmax, dr = self.r, self.rmax, self.dr
        still = set()
        for i in self.low:
            v = r[i] + dr
            if v < rmax:
                r[i] = v
                still.add(i)
            else:
                r[i] = rmax
        self.low = still

    # -- one core update --------------------------------------------------------
    def update(self) -> None:
        core, L, M = self.core, self.L, self.M
        queue = self.queue
        n = len(queue)
        new: list[int] = []
        pending = []
        stats, executed = self.stats, self.executed
        for i, x in enumerate(queue):
            if not self.consume(x):
                new.append(x)
                stats["waits"] += 1
                continue
            executed[core[x][0]] += 1
            w: dict = {}
            succ, src = self.execute(x, w)
            if not succ:
                stats["deaths"] += 1
            elif len(succ) == 2:
                new.append(succ[0])
                if len(new) + (n - i - 1) < L:
                    new.append(succ[1])
                    stats["splits"] += 1
                else:
                    stats["splits_refused"] += 1
            else:
                new.append(succ[0])
            if not w:
                continue
            if self.parallel:
                pending.append((x, src, w))
            else:
                if self.record:
                    self._event(x, src, w)
                for addr, word in w.items():
                    core[addr] = word
        if pending:
            last = {}
            for k, (_, _, w) in enumerate(pending):
                for addr in w:
                    last[addr] = k
            effective = []
            for k, (x, src, w) in enumerate(pending):
                eff = {addr: word for addr, word in w.items() if last[addr] == k}
                if len(eff) < len(w):
                    stats["overwritten_writes"] += len(w) - len(eff)
                if eff:
                    if self.record:
                        self._event(x, src, eff)
                    effective.append(eff)
            for eff in effective:
                for addr, word in eff.items():
                    core[addr] = word
        if self.ppoint and self.rand() < self.ppoint:
            if len(new) < L:
                new.append(int(self.rand() * M))
                stats["injections"] += 1
            else:
                stats["injections_refused"] += 1
        self.queue = new
        self.renew()
        self.time += 1

    def _event(self, x: int, src, w: dict) -> None:
        core = self.core
        cells = set(w)
        cells.add(x)
        if src is not None:
            cells.add(src)
        before = sorted(core[c] for c in cells)
        after = sorted(w.get(c) or core[c] for c in cells)
        self.stats["writes"] += 1
        if self.key is not None:
            before = sorted(self.key(t) for t in before)
            after = sorted(self.key(t) for t in after)
        if before == after:
            self.stats["elastic_writes"] += 1
            return
        k = (tuple(before), tuple(after))
        entry = self.events.get(k)
        if entry is None:
            self.events[k] = [1]
        else:
            entry[0] += 1

    key = None   # set to a word -> species-key function, or None for whole words

    # -- inspection -------------------------------------------------------------
    def self_loop(self, x: int) -> bool:
        """True if the pointer at x is at a fixed point: a jump that sends it back to x."""
        word = self.core[x]
        if word[0] not in (JMP, JMZ, JMN, DJN):
            return False
        succ, _ = self.execute(x, {}, mutate=False)
        return succ == (x,)


# ----------------------------------------------------------------------------
# set-up and analysis

def _buffered(rng, size=8192):
    buf: list[float] = []

    def rand() -> float:
        if not buf:
            buf.extend(rng.random(size).tolist()[::-1])
        return buf.pop()

    return rand


def random_core(rand, core_size: int) -> list[tuple]:
    return [(int(rand() * 10), int(rand() * 4), int(rand() * core_size), int(rand() * 4),
             int(rand() * core_size)) for _ in range(core_size)]


def run_lengths(core: list[tuple]) -> dict[str, dict[int, int]]:
    """Fig. 4 statistic: how often each instruction occurs alone, in pairs, triples, ... (cyclic core)."""
    ops = [w[0] for w in core]
    n = len(ops)
    out = {name: Counter() for name in OPS}
    if len(set(ops)) == 1:
        out[OPS[ops[0]]][n] += 1
    else:
        start = next(i for i in range(n) if ops[i] != ops[i - 1])
        length = 1
        for k in range(1, n + 1):
            i = (start + k) % n
            if k < n and ops[i] == ops[(i - 1) % n]:
                length += 1
            else:
                out[OPS[ops[(i - 1) % n]]][length] += 1
                length = 1
    return {name: {str(k): v for k, v in sorted(c.items())} for name, c in out.items() if c}


def composition(core: list[tuple]) -> dict[str, int]:
    c = Counter(w[0] for w in core)
    return {name: c.get(i, 0) for i, name in enumerate(OPS)}


def simulate(machine: Venus, updates: int, samples: int = 100) -> dict:
    every = max(1, updates // samples) if updates else 1
    series = {"time": [], "pointers": [], "executions": [], "mean_resource": [],
              "composition": {name: [] for name in OPS}}
    done = 0
    for t in range(updates):
        before = sum(machine.executed.values())
        machine.update()
        done += 1
        if t % every == 0 or t == updates - 1:
            series["time"].append(machine.time)
            series["pointers"].append(len(machine.queue))
            series["executions"].append(sum(machine.executed.values()) - before)
            series["mean_resource"].append(round(sum(machine.r) / machine.M, 6))
            for name, v in composition(machine.core).items():
                series["composition"][name].append(v)
    return series


# ----------------------------------------------------------------------------

def generate(p, rng):
    M = p.core_size
    program, start = SEEDS[p.seed_program]
    if len(program) > M:
        raise ValueError(f"seed_program {p.seed_program!r} has {len(program)} words, more than core_size = {M}")
    rand = _buffered(rng)
    core = random_core(rand, M)
    origin = int(rand() * M)
    for k, line in enumerate(program):
        core[(origin + k) % M] = parse(line, M)
    machine = Venus(core, rand, mode=p.mode, queue_len=p.queue_len,
                    operation_radius=p.operation_radius, resource_radius=p.resource_radius,
                    resource_influx=p.resource_influx, resource_max=p.resource_max,
                    mutation_rate=p.mutation_rate, injection_rate=p.pointer_injection_rate)
    by_word = p.species == "word"
    if not by_word:
        machine.key = lambda w: w[0]
    if start is not None:
        machine.queue.append((origin + start) % M)
    initial_core = list(core)
    series = simulate(machine, p.updates)

    def sid(k) -> str:
        return word_id(k, M) if by_word else OPS[k]

    reactions = []
    for (lhs, rhs), (count,) in machine.events.items():
        reactions.append(Reaction(dict(Counter(sid(k) for k in lhs)), dict(Counter(sid(k) for k in rhs)),
                                  count=count))
    if by_word:
        seen = dict.fromkeys(initial_core)
        for lhs, rhs in machine.events:
            seen.update(dict.fromkeys(lhs))
            seen.update(dict.fromkeys(rhs))
        species = [Species(word_id(w, M), structure=text(w, M)) for w in seen]
        initial = Counter(word_id(w, M) for w in initial_core)
    else:
        species = [Species(name) for name in OPS]
        initial = Counter(OPS[w[0]] for w in initial_core)
    stats = machine.stats
    final_queue = machine.queue
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={s.id: initial[s.id] for s in species if initial[s.id]},
        extras={
            "conservation": [{"name": "core cells", "vector": {s.id: 1 for s in species}}],
            "space": {"geometry": "cyclic one-dimensional core", "size": M,
                      "operation_radius": p.operation_radius, "resource_radius": p.resource_radius},
            "analysis": {
                "updates": machine.time,
                "series": series,
                "executed": {OPS[k]: v for k, v in sorted(machine.executed.items())},
                "events": {k: int(stats[k]) for k in (
                    "writes", "elastic_writes", "mutations", "splits", "splits_refused", "deaths",
                    "waits", "injections", "injections_refused", "overwritten_writes")},
                "final_composition": composition(machine.core),
                "final_run_lengths": run_lengths(machine.core),
                "final_pointers": len(final_queue),
                "self_loop_pointers": sum(machine.self_loop(x) for x in final_queue),
            },
            "coreworld": {
                "mode": p.mode,
                "seed_program": p.seed_program,
                "seed_origin": origin,
                "event": ("reactants = the executing word, the word it copies/adds from (MOV, ADD, SUB "
                          "with non-immediate A) and the words it changes; products = the same cells after "
                          "the write; elastic events are not listed"),
                "final_core": [text(w, M) for w in machine.core],
                "final_queue": list(final_queue),
            },
        },
    )
