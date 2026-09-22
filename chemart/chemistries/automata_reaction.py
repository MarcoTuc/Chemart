"""Automata reaction: Dittrich & Banzhaf (1998), 32-bit binary string chemistry.

Catalog id: automata-reaction.

Molecules are 32-bit words. In a collision s1 + s2 => s3, s1 is "folded" into
a program of eight 4-bit instructions that runs on a small register machine:
the operator register holds s1 (read only), the IO register starts as s2 and
holds the product s3 when the program halts. Both reactants are catalysts:

    s1 + s2 -> s1 + s2 + s3        (and a random molecule is diluted out)

The machine is a line-by-line port of the authors' ANSI C source
(autoreac-1.0, autoReac.c, cited as [8] in the paper), which the paper names
as the formal specification. It reproduces the paper's Table 2 and reaction
tables (Figs. 3, 4, 9, 10); see the tests.

evolve runs the paper's reactor algorithm (chemart.soup.stir), a frame per
generation of M collisions, and returns the observed reactions; generate
returns the reaction closure of the seed words (e.g. a published
organization).
"""

from __future__ import annotations

from collections import Counter

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

MASK = 0xFFFFFFFF

# autoReac.c: instructionTable01 and instructionTable02 (paper Fig. 1, right).
# They differ only at 1010: NOT (table 1) vs EQ (table 2).
_COMMON = ["ID", "MOV", "SETP", "TMM", "TDIR", "UNSETP", "CPON", "CPOFF",
           "STOP", "OR", None, "EXOR", "EXOR", "NOP", "ID", "AND"]
CODE_TABLES = {
    1: [op if op else "NOT" for op in _COMMON],
    2: [op if op else "EQ" for op in _COMMON],
}
LOGIC = frozenset({"ID", "NOT", "AND", "OR", "EXOR", "EQ"})

_LEFT, _RIGHT = 0, 1
_BOTH, _OPERATOR_PTR, _IO_PTR = 0, 1, 2   # autoReac.c: BOTH, OPERAND_PT (opReg), RESULT_PT


def word_id(w: int) -> str:
    return f"w{w:08x}"


def disassemble(s1: int, table: int = 1) -> list[str]:
    """The program folded from s1: nibbles from the least significant one, SETP takes the next."""
    codes = CODE_TABLES[table]
    out, i = [], 0
    while i < 8:
        op = codes[(s1 >> (4 * i)) & 0xF]
        if op == "SETP":
            i += 1
            arg = (s1 >> (4 * i)) & 0xF if i < 8 else 0
            op = f"SETP {arg:04b}"
        out.append(op)
        i += 1
    return out


def _logic(op: str, io: int, opr: int) -> int:
    """ALU on bit 0 of both registers (the pointers are kept at bit 0 by rotation)."""
    b = opr & 1
    if op == "ID":
        return (io & ~1 & MASK) | b
    if op == "NOT":
        return (io & ~1 & MASK) | (b ^ 1)
    if op == "AND":
        return io & (b | (MASK ^ 1))
    if op == "OR":
        return io | b
    if op == "EXOR":
        return io ^ b
    # EQ: invert bit 0, then EXOR (mEq)
    return (io ^ 1) ^ b


def automata(s1: int, s2: int, table: int = 1) -> int:
    """s1 + s2 => s3 by the automata reaction with the given code table (autoReac.c runN)."""
    codes = CODE_TABLES[table]
    opr, io, io_pos = s1, s2, 0
    pattern = None
    direction, mode = _LEFT, _BOTH
    copy, copy_fn = False, "ID"

    def step():
        nonlocal opr, io, io_pos
        if copy:
            io = _logic(copy_fn, io, opr)
        if direction == _LEFT:   # pointer moves to the next higher bit: rotate right
            if mode != _IO_PTR:
                opr = ((opr >> 1) | ((opr & 1) << 31)) & MASK
            if mode != _OPERATOR_PTR:
                io = ((io >> 1) | ((io & 1) << 31)) & MASK
                io_pos += 1
        else:
            if mode != _IO_PTR:
                opr = ((opr << 1) | (opr >> 31)) & MASK
            if mode != _OPERATOR_PTR:
                io = ((io << 1) | (io >> 31)) & MASK
                io_pos -= 1

    i = 0
    while i < 8:
        op = codes[(s1 >> (4 * i)) & 0xF]
        i += 1
        if op in LOGIC:
            copy_fn = op
            io = _logic(op, io, opr)
            step()
        elif op == "MOV":
            if pattern is None:
                step()
            else:
                for _ in range(32):
                    step()
                    reg = io if mode == _IO_PTR else opr
                    if reg & 0xF == pattern:
                        break
        elif op == "SETP":
            pattern = (s1 >> (4 * i)) & 0xF if i < 8 else 0
            i += 1
        elif op == "UNSETP":
            pattern = None
        elif op == "TMM":
            mode = (mode + 1) % 3
        elif op == "TDIR":
            direction = _RIGHT if direction == _LEFT else _LEFT
        elif op == "CPON":
            copy = True
        elif op == "CPOFF":
            copy = False
        elif op == "STOP":
            break
        # NOP: nothing
    # correctResult: undo the rotation of the IO register
    k = io_pos % 32
    return ((io << k) | (io >> (32 - k))) & MASK if k else io


def and_reaction(s1: int, s2: int) -> int:
    """The paper's reference AND reaction (section 3.1)."""
    return s1 & s2


def make_react(mechanism: str, table: int, forbid_exact_replication: bool):
    """react(s1, s2) -> s3, or None for an elastic collision (filter f1, paper eq. 3)."""
    cache: dict[tuple[int, int], int | None] = {}

    def react(s1: int, s2: int):
        key = (s1, s2)
        if key not in cache:
            s3 = automata(s1, s2, table) if mechanism == "automata" else s1 & s2
            if forbid_exact_replication and (s3 == s1 or s3 == s2):
                s3 = None
            cache[key] = s3
        return cache[key]

    return react


def _catalytic(react):
    def full(s1, s2):
        s3 = react(s1, s2)
        return None if s3 is None else (s1, s2, s3)
    return full


def _rate(react, lhs, rhs) -> dict:
    """Paper eq. 2: k = 1 per ordered pair (i, j) whose product is s3, summed over both orders."""
    a, b = lhs
    product = Counter(rhs) - Counter(lhs)
    (s3,) = product.elements()
    orders = {(a, b), (b, a)}
    return {"law": "mass-action", "k": float(sum(1 for x, y in orders if react(x, y) == s3))}


def _parse_words(words) -> list[int]:
    out = []
    for w in words:
        if isinstance(w, str):
            text = w[1:] if w.startswith("w") else w
            try:
                value = int(text, 16)
            except ValueError:
                value = -1
        elif isinstance(w, int) and not isinstance(w, bool):
            value = w
        else:
            value = -1
        if not 0 <= value <= MASK:
            raise ValueError(f"words must be 32-bit words as hex strings ('1e1ca260') or ints, got {w!r}")
        out.append(value)
    return out


def _reaction(react, lhs, rhs, count=None) -> Reaction:
    return Reaction.of([word_id(w) for w in lhs], [word_id(w) for w in rhs],
                       rate=_rate(react, lhs, rhs), count=count)


def _species(words) -> list[Species]:
    return [Species(word_id(w), structure=f"{w:032b}") for w in sorted(set(words))]


def _check(p) -> None:
    if p.mechanism == "and" and p.code_table != 1:
        raise ValueError("code_table only applies to mechanism 'automata'; leave it at 1 for 'and'")


def _random_words(rng, n: int) -> list[int]:
    return [int(w) for w in rng.integers(0, 1 << 32, size=n, dtype="uint64")]


def generate(p, rng):
    """Every reaction reachable from the distinct seed words, cut off by max_species."""
    _check(p)
    start = _parse_words(p.words) if p.words else _random_words(rng, p.n_seeds)
    react = make_react(p.mechanism, p.code_table, p.forbid_exact_replication)
    seed = sorted(set(start))
    words, found, status = expand(_catalytic(react), seed, arity=2, max_species=p.max_species, ordered=True)
    return Network(
        species=_species(words),
        reactions=[_reaction(react, lhs, rhs) for lhs, rhs in found],
        status=status,
        outflow=CONSTANT_TOTAL,
        extras={"seed": [word_id(w) for w in seed]},
    )


def _ids(fired):
    return [[[word_id(w) for w in lhs], [word_id(w) for w in rhs], n] for lhs, rhs, n in fired]


def evolve(p, rng):
    """The paper's reactor: `generations` x M collisions, a frame per generation.

    Each frame after the first reports the paper's productivity (inserted
    products / M) and innovativity (never-seen products / M) of that generation.
    """
    _check(p)
    start = _parse_words(p.words) if p.words else _random_words(rng, p.M)
    react = make_react(p.mechanism, p.code_table, p.forbid_exact_replication)
    if len(start) < 2:
        raise ValueError(f"the soup needs at least 2 molecules, got {len(start)}")
    size = len(start)
    seen = set(start)
    tally = Tally()
    pop = start
    for step, pop, tally in stir(_catalytic(react), start, p.generations * size, rng, every=size, arity=2,
                                 dilution="constant", tally=tally):
        fired = tally.flush()
        observables = {}
        if step:
            new = 0
            for lhs, rhs, _ in fired:
                (s3,) = (Counter(rhs) - Counter(lhs)).elements()
                if s3 not in seen:
                    seen.add(s3)
                    new += 1
            observables = {"productivity": sum(n for _, _, n in fired) / size, "innovativity": new / size}
        yield Frame(t=float(step // size), state={word_id(w): float(n) for w, n in sorted(Counter(pop).items())},
                    fired=_ids(fired), observables=observables)

    fired = tally.reactions()
    reactions = [_reaction(react, lhs, rhs, count) for lhs, rhs, count in fired]
    words = set(start).union(*(rhs for _, rhs, _ in fired))
    final = Counter(pop)
    return Network(
        species=_species(words),
        reactions=reactions,
        status="observed",
        initial_state={word_id(w): c for w, c in sorted(Counter(start).items())},
        outflow=CONSTANT_TOTAL,
        extras={
            "analysis": {"generation_size": size},
            "final_state": {word_id(w): c for w, c in final.most_common()},
        },
    )
