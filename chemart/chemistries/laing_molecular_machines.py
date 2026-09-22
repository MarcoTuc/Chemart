"""Laing's artificial molecular machines (Laing 1972-1977). Catalog id: laing-molecular-machines.

Molecules are strings of constituents in one of two forms (book 10.5.1):

    passive tape     t:0110            binary constituents in state 0 or 1
    active machine   m:CT1.W1.H.TT1    a folded chain of instruction constituents

The instruction set is the one of book fig. 10.5 (from Laing 1975): W0, W1, L,
R, NOP and the branch pair CT/TT, plus H (halt) and D (detach) from the list of
Laing's 1975 constituents in Freitas & Merkle (2004, 4.8). The fold that brings
a conditional transfer CT into contact with its transfer target TT is written
as a shared label: CT1 ... TT1.

A reaction is one machine attached to one tape at a single position of each:
the machine starts at its first instruction and runs the instruction in
contact with the tape, then moves its contact to the next instruction:

    W0 / W1   put the contacted tape constituent in state 0 / 1
    L / R     slide the contact one constituent left / right; sliding off the end
              recruits a new 0 constituent there (Freitas & Merkle 4.8)
    CT<k>     if the contacted constituent is 1, continue at TT<k> (Wang's
              conditional transfer); otherwise continue with the next instruction
    TT<k>     transfer target, no operation when reached in sequence
    NOP       no operation
    D         sever the contacted constituent from its predecessors: the part to
              its left is released as a separate tape
    H         halt; running past the last instruction also halts

The machine is a catalyst; the tape is replaced by its pieces:

    m + t -> m + t_1 + ... + t_n

generate returns every reaction reachable from the seed machines and tapes
(chemart.expand.expand); evolve draws random pairs from a population
(chemart.soup.stir), a frame every generation, and returns the reactions that
fired.
turing_program compiles a 2-symbol Turing machine into a Laing machine, which is
how the tests exercise universal computation (book 10.5.1, Laing [485]).
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from itertools import product as cartesian

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

SIMPLE = ("W0", "W1", "L", "R", "H", "NOP", "D")
_TOKEN = re.compile(r"^(?:W0|W1|L|R|H|NOP|D|(CT|TT)([A-Za-z0-9_]+))$")
_TAPE = re.compile(r"^[01]+$")
BINDINGS = ("leftmost", "all", "random")


# ---------------------------------------------------------------------------
# molecules
def machine_id(program: list[str]) -> str:
    return "m:" + ".".join(program)


def tape_id(bits: str) -> str:
    return "t:" + bits


def parse_program(text: str) -> list[str]:
    """Instruction tokens of a machine written 'CT1.W1.H' (an 'm:' prefix is allowed)."""
    if not isinstance(text, str):
        raise ValueError(f"machines must be strings like 'CT1.W1.H.TT1.W0', got {text!r}")
    body = text[2:] if text.startswith("m:") else text
    tokens = body.split(".") if body else []
    bad = [t for t in tokens if not _TOKEN.match(t)]
    if not tokens or bad:
        raise ValueError(f"machine {text!r}: instructions are {', '.join(SIMPLE)}, CT<label> and "
                         f"TT<label> separated by '.', got {bad or 'nothing'}")
    targets = Counter(t[2:] for t in tokens if t.startswith("TT"))
    doubled = sorted(k for k, n in targets.items() if n > 1)
    missing = sorted({t[2:] for t in tokens if t.startswith("CT")} - set(targets))
    if doubled or missing:
        raise ValueError(f"machine {text!r}: every CT<label> needs exactly one TT<label> "
                         f"(the fold joins them); duplicate targets {doubled}, missing targets {missing}")
    return tokens


def parse_tape(text: str) -> str:
    if not isinstance(text, str):
        raise ValueError(f"tapes must be strings over 0 and 1, got {text!r}")
    body = text[2:] if text.startswith("t:") else text
    if not _TAPE.match(body):
        raise ValueError(f"tapes must be non-empty strings over 0 and 1, got {text!r}")
    return body


def fold(program: list[str]) -> str:
    """The machine's shape: each instruction with its index, CTs pointing at their TT."""
    where = {t[2:]: i for i, t in enumerate(program) if t.startswith("TT")}
    return " ".join(f"{i}:{t}" + (f"->{where[t[2:]]}" if t.startswith("CT") else "")
                    for i, t in enumerate(program))


# ---------------------------------------------------------------------------
# the machine
@dataclass
class Run:
    pieces: list[str]              # released tapes, left to right, then the one in contact
    outcome: str                   # halt | end | steps | length
    steps: int = 0
    counts: Counter = field(default_factory=Counter)

    @property
    def halted(self) -> bool:
        return self.outcome in ("halt", "end")


def run(program: list[str], tape: str, site: int = 0, max_steps: int = 1000,
        max_length: int | None = None) -> Run:
    """Attach `program` at its first instruction to unit `site` of `tape` and run it."""
    if not 0 <= site < len(tape):
        raise ValueError(f"binding site {site} is outside the tape {tape!r}")
    target = {t[2:]: i for i, t in enumerate(program) if t.startswith("TT")}
    cells, pos, released = list(tape), site, []
    counts: Counter = Counter()
    pc = steps = 0
    outcome = "end"
    while pc < len(program):
        if steps >= max_steps:
            outcome = "steps"
            break
        op = program[pc]
        steps += 1
        counts[op[:2] if op[:2] in ("CT", "TT") else op] += 1
        pc += 1
        if op == "W0" or op == "W1":
            cells[pos] = op[1]
        elif op == "R":
            pos += 1
            if pos == len(cells):
                cells.append("0")
        elif op == "L":
            if pos == 0:
                cells.insert(0, "0")
            else:
                pos -= 1
        elif op.startswith("CT"):
            if cells[pos] == "1":
                pc = target[op[2:]] + 1       # the target itself is a no-op
        elif op == "D":
            if pos > 0:
                released.append("".join(cells[:pos]))
                cells, pos = cells[pos:], 0
        elif op == "H":
            outcome = "halt"
            break
        if max_length is not None and len(cells) > max_length:
            outcome = "length"
            break
    return Run(released + ["".join(cells)], outcome, steps, counts)


# ---------------------------------------------------------------------------
# Turing machines as Laing machines
def turing_program(table: dict, start: str = "A", halt: str = "H") -> list[str]:
    """Compile a 2-symbol Turing machine into Laing instructions.

    table maps (state, symbol) to (write, move, next state), symbols 0/1, move
    'L' or 'R', next state `halt` to stop. CT only branches on 1, so jumping to
    state q is 'CTq1 W1 CTq0': a 0 cell is turned into 1 to force the transfer
    and restored to 0 at TTq0 (Chemart construction, not Laing's).
    """
    states = list(dict.fromkeys(q for q, _ in table))

    def goto(q):
        return ["H"] if q == halt else [f"CT{q}1", "W1", f"CT{q}0"]

    program = goto(start)
    for q in states:
        for symbol in (0, 1):
            if (q, symbol) not in table:
                raise ValueError(f"turing table misses ({q!r}, {symbol})")
            write, move, nxt = table[(q, symbol)]
            if move not in ("L", "R"):
                raise ValueError(f"move must be 'L' or 'R', got {move!r}")
            if nxt != halt and nxt not in states:
                raise ValueError(f"unknown next state {nxt!r}")
            program += [f"TT{q}{symbol}"] + (["W0"] if symbol == 0 else [])
            program += [f"W{int(write)}", move] + goto(nxt)
    return program


# ---------------------------------------------------------------------------
class _Chemistry:
    def __init__(self, p, rng, either_order: bool = False):
        self.p, self.rng = p, rng
        self.either_order = either_order    # the soup reacts machine + tape drawn in either order
        self.truncated = False
        self.cache: dict[tuple, list[tuple]] = {}

    def outcomes(self, machine: str, tape: str, fresh: bool = False) -> list[tuple[str, ...]]:
        """Right-hand sides (machine included) of machine + tape, one per binding choice."""
        key = (machine, tape)
        if not fresh and key in self.cache:
            return self.cache[key]
        program, bits = parse_program(machine), tape[2:]
        if self.p.binding == "leftmost":
            sites = [0]
        elif self.p.binding == "all":
            sites = list(range(len(bits)))
        else:
            sites = [int(self.rng.integers(len(bits)))]
        out = []
        for site in sites:
            r = run(program, bits, site, self.p.max_steps, self.p.max_length)
            if not r.halted:
                self.truncated = True
                continue
            rhs = (machine, *sorted(tape_id(b) for b in r.pieces))
            if Counter(rhs) != Counter((machine, tape)) and rhs not in out:
                out.append(rhs)
        if not fresh:
            self.cache[key] = out
        return out

    def results(self, lhs, fresh: bool = False) -> list[tuple[str, ...]]:
        a, b = lhs
        if a.startswith("m:") and b.startswith("t:"):
            return self.outcomes(a, b, fresh)
        if self.either_order and a.startswith("t:") and b.startswith("m:"):
            return self.outcomes(b, a, fresh)
        return []                      # machine + machine, tape + tape: no reaction


def _species(ids) -> list[Species]:
    return [Species(s, structure="active: " + fold(parse_program(s)) if s.startswith("m:") else "passive")
            for s in ids]


def _seeds(p) -> tuple[list[str], list[str]]:
    if not isinstance(p.machines, list) or not p.machines:
        raise ValueError("machines must be a non-empty list of instruction strings like 'CT1.W1.H.TT1.W0'")
    if not isinstance(p.tapes, list) or not p.tapes:
        raise ValueError("tapes must be a non-empty list of strings over 0 and 1")
    machines = list(dict.fromkeys(machine_id(parse_program(m)) for m in p.machines))
    tapes = list(dict.fromkeys(tape_id(parse_tape(t)) for t in p.tapes))
    if any(len(t) - 2 > p.max_length for t in tapes):
        raise ValueError(f"seed tapes must not be longer than max_length={p.max_length}")
    return machines, tapes


def generate(p, rng) -> Network:
    """Every reaction reachable from the seed machines and tapes, cut off by max_species."""
    machines, tapes = _seeds(p)
    return _closure(p, _Chemistry(p, rng), machines + tapes)


def evolve(p, rng):
    """A well-stirred soup of `copies` of each seed (Chemart addition), a frame every generation."""
    machines, tapes = _seeds(p)
    if p.binding == "all":
        raise ValueError("evolve needs one outcome per collision: use binding leftmost or "
                         "random (generate_network takes 'all')")
    return (yield from _soup(p, _Chemistry(p, rng, either_order=True), machines, tapes, rng))


def _closure(p, chem, seed) -> Network:
    def discover(*lhs):
        found = [s for rhs in chem.results(lhs) for s in rhs]
        return tuple(dict.fromkeys(found)) if found else None

    ids, _, status = expand(discover, seed, arity=2, max_species=p.max_species, ordered=True)
    known = set(ids)
    reactions = []
    for lhs in cartesian(ids, repeat=2):
        for rhs in chem.results(lhs):
            if set(rhs) <= known:
                reactions.append(Reaction.of(lhs, rhs))
            else:
                status = "truncated"
    if chem.truncated:
        status = "truncated"
    return Network(species=_species(ids), reactions=reactions, status=status, extras={"seed": seed})


def _soup(p, chem, machines, tapes, rng):
    def react(*lhs):
        found = chem.results(lhs, fresh=p.binding == "random")
        return found[0] if found else None

    start = [s for s in machines + tapes for _ in range(p.copies)]
    tally = Tally()
    final = start
    for step, final, tally in stir(react, start, p.steps, rng, arity=2, dilution="constant", tally=tally):
        yield Frame(t=float(step), state={s: float(n) for s, n in Counter(final).items()},
                    fired=[[list(lhs), list(rhs), n] for lhs, rhs, n in tally.flush()])
    fired = tally.reactions()
    ids = list(dict.fromkeys([*start, *(s for _, rhs, _ in fired for s in rhs)]))
    return Network(
        species=_species(ids),
        reactions=[Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in fired],
        status="observed",
        initial_state={s: c for s, c in Counter(start).items()},
        outflow=CONSTANT_TOTAL,
        extras={"final_state": dict(Counter(final).most_common())},
    )
