"""BFF, self-modifying Brainfuck programs: Agüera y Arcas et al. (2024).

Catalog id: bff.

Molecules are 64-byte programs ("tapes"). A reaction draws two of them, A and
B, concatenates them into one 128-byte tape AB and runs it as a BFF program for
at most ``max_steps`` steps; the tape is then split back into A' and B':

    A + B -> split(exec(AB)) = A' + B'

Code and data share the tape, so a program can rewrite itself and its partner.
The instruction set (paper section 2; cubff ``bff_noheads``) acts on two heads,
head0 and head1, which start at 0 and wrap around the 128 bytes:

    <  head0 -= 1        >  head0 += 1        {  head1 -= 1     }  head1 += 1
    -  tape[head0] -= 1  +  tape[head0] += 1
    .  tape[head1] = tape[head0]              ,  tape[head0] = tape[head1]
    [  if tape[head0] == 0: jump forward to the matching ]
    ]  if tape[head0] != 0: jump back to the matching [

Every other byte is a no-op. A jump without a match ends the program, and so
does the instruction pointer leaving the tape. Every byte read counts as a
step. Arithmetic wraps modulo 256.

The reactor is the paper's "primordial soup": each epoch shuffles the
programs, pairs them up and runs every pair once, after background mutation
has replaced each byte with a random one with probability ``mutation_rate``.
``space="grid"`` is the paper's 2D variant (section 2.2): programs sit on a
grid and pair only with a neighbour at most two cells away along each axis.

The chemistry is a gas with one face, ``evolve``: a frame per epoch, whose
observables are the paper's high-order entropy and three simpler traces. The
network is the record of one run: the distinct reactions A + B -> A' + B'
that changed a tape, and the mutations A -> A', with counts. A self-replicator
S shows up as S + F -> S + S, the paper's equation (5).
"""

from __future__ import annotations

from collections import Counter

import brotli
import numpy as np

from chemart.network import Network, Reaction, Species
from chemart.soup import Tally
from chemart.trajectory import Frame, ticks

TAPE = 64
COMMANDS = "<>{}-+.,[]"

# The hand-written self-replicator of the paper's Figure 4: a palindrome of
# BFF code padded with spaces (non-coding bytes). Concatenated with any
# program B, it copies itself into B, in reverse, which for a palindrome is
# an exact copy.
FIG4_REPLICATOR = "[[{.>]-]" + " " * 48 + "]-]>.{[["

# Opcodes, indexed by byte value.
_NOOP, _DEC0, _INC0, _DEC1, _INC1, _MINUS, _PLUS, _COPY01, _COPY10, _OPEN, _CLOSE = range(11)
_OP = [_NOOP] * 256
for _code, _ch in enumerate(COMMANDS, start=1):
    _OP[ord(_ch)] = _code
_OP = tuple(_OP)


# --- programs -------------------------------------------------------------------
def encode(text: str) -> bytes:
    """A program written as text: each character is one byte, '0' is the zero byte.

    Other characters must be Latin-1; cubff's symbols for non-coding bytes
    (U+0100 + byte) are accepted too, so ``encode(show(t)) == t``.
    """
    out = bytearray()
    for ch in text:
        c = ord(ch)
        if ch == "0":
            out.append(0)
        elif c < 256:
            out.append(c)
        elif 0x100 <= c < 0x200:
            out.append(c - 0x100)
        else:
            raise ValueError(f"character {ch!r} is not a byte")
    return bytes(out)


def show(tape: bytes) -> str:
    """Printable form, one character per byte, after cubff: BFF commands as
    themselves, the zero byte as '0', any other byte b as the character U+0100 + b.
    """
    return "".join(
        ch if (ch := chr(b)) in COMMANDS else "0" if b == 0 else chr(0x100 + b)
        for b in tape
    )


def run(tape: bytearray, max_steps: int = 8192) -> int:
    """Execute a BFF tape in place; return the number of non-no-op instructions.

    Follows cubff's ``Bff::Evaluate``: the heads start at 0 and wrap modulo the
    tape length; every byte read is one step.
    """
    n = len(tape)
    op = _OP
    pc = head0 = head1 = 0
    executed = 0
    for _ in range(max_steps):
        code = op[tape[pc]]
        if code:
            executed += 1
            if code == _DEC0:
                head0 = (head0 - 1) % n
            elif code == _INC0:
                head0 = (head0 + 1) % n
            elif code == _DEC1:
                head1 = (head1 - 1) % n
            elif code == _INC1:
                head1 = (head1 + 1) % n
            elif code == _MINUS:
                tape[head0] = (tape[head0] - 1) & 0xFF
            elif code == _PLUS:
                tape[head0] = (tape[head0] + 1) & 0xFF
            elif code == _COPY01:
                tape[head1] = tape[head0]
            elif code == _COPY10:
                tape[head0] = tape[head1]
            elif code == _OPEN:
                if tape[head0] == 0:
                    depth = 1
                    pc += 1
                    while pc < n:
                        c = op[tape[pc]]
                        if c == _CLOSE:
                            depth -= 1
                            if depth == 0:
                                break
                        elif c == _OPEN:
                            depth += 1
                        pc += 1
                    if depth:
                        return executed
            else:  # _CLOSE
                if tape[head0] != 0:
                    depth = 1
                    pc -= 1
                    while pc >= 0:
                        c = op[tape[pc]]
                        if c == _OPEN:
                            depth -= 1
                            if depth == 0:
                                break
                        elif c == _CLOSE:
                            depth += 1
                        pc -= 1
                    if depth:
                        return executed
        pc += 1
        if pc >= n:
            break
    return executed


def high_order_entropy(soup: bytes) -> float:
    """Shannon entropy of the bytes minus the Brotli-compressed bits per byte.

    The paper's complexity measure (section 2.1), computed as cubff does:
    Brotli at quality 2 with a 2^24 window stands in for Kolmogorov complexity.
    About 0 for random bytes; high when the soup holds many copies of a string.
    """
    counts = np.bincount(np.frombuffer(soup, dtype=np.uint8), minlength=256)
    frac = counts[counts > 0] / len(soup)
    h0 = float(-(frac * np.log2(frac)).sum())
    compressed = brotli.compress(soup, quality=2, lgwin=24)
    return h0 - len(compressed) * 8.0 / len(soup)


# --- the reactor ----------------------------------------------------------------
def _neighbours(width: int, height: int, radius: int = 2) -> list[list[int]]:
    """Grid cells within `radius` along each axis, no wrap-around (cubff make_2d_pattern.py)."""
    out = []
    for y in range(height):
        for x in range(width):
            out.append([
                yy * width + xx
                for yy in range(max(0, y - radius), min(height, y + radius + 1))
                for xx in range(max(0, x - radius), min(width, x + radius + 1))
                if (yy, xx) != (y, x)
            ])
    return out


def _pairs(p, rng, neighbours) -> list[tuple[int, int]]:
    order = rng.permutation(p.tapes)
    if neighbours is None:
        return [(int(order[i]), int(order[i + 1])) for i in range(0, p.tapes - 1, 2)]
    used = np.zeros(p.tapes, dtype=bool)
    out = []
    for i in order:
        i = int(i)
        options = neighbours[i]
        j = options[int(rng.integers(len(options)))]
        if used[i] or used[j]:
            continue
        used[i] = used[j] = True
        out.append((i, j))
    return out


def _key(lhs, rhs) -> tuple:
    return (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))


def evolve(p, rng):
    """The primordial soup: a frame per epoch (frame 0 is the initial soup)."""
    if p.replicators > p.tapes:
        raise ValueError(f"replicators ({p.replicators}) cannot exceed tapes ({p.tapes})")
    neighbours = None
    if p.space == "grid":
        if p.tapes % p.width:
            raise ValueError(f"space='grid' needs tapes ({p.tapes}) to be a multiple of width ({p.width})")
        height = p.tapes // p.width
        if p.width * height < 2:
            raise ValueError("space='grid' needs at least two cells")
        neighbours = _neighbours(p.width, height)

    soup = rng.integers(0, 256, size=(p.tapes, TAPE), dtype=np.uint8)
    if p.replicators:
        seeded = np.frombuffer(encode(FIG4_REPLICATOR), dtype=np.uint8)
        for i in rng.choice(p.tapes, size=p.replicators, replace=False):
            soup[int(i)] = seeded

    initial = Counter(show(bytes(row)) for row in soup)
    tally = Tally()
    kinds: dict[tuple, str] = {}

    def record(lhs, rhs, kind):
        kinds.setdefault(_key(lhs, rhs), kind)
        tally.add(lhs, rhs)

    def frame(epoch, ops, runs):
        flat = soup.tobytes()
        tapes = Counter(show(soup[i].tobytes()) for i in range(p.tapes))
        return Frame(t=float(epoch), state={s: float(n) for s, n in tapes.items()}, fired=tally.flush(),
                     observables={"high_order_entropy": round(high_order_entropy(flat), 4),
                                  "top_tape_count": tapes.most_common(1)[0][1],
                                  "ops_per_run": round(ops / runs, 2) if runs else 0.0,
                                  "zero_bytes": int(flat.count(0))})

    yield frame(0, 0, 0)
    for epoch in ticks(p.epochs):
        # Background mutation, before execution, of every byte of every program.
        if p.mutation_rate > 0:
            hit = rng.random(soup.shape) < p.mutation_rate
            rows = np.flatnonzero(hit.any(axis=1))
            if rows.size:
                before = [show(soup[i].tobytes()) for i in rows]
                soup[hit] = rng.integers(0, 256, size=int(hit.sum()), dtype=np.uint8)
                for old, i in zip(before, rows):
                    new = show(soup[i].tobytes())
                    if new != old:
                        record((old,), (new,), "mutation")
        ops = 0
        pairs = _pairs(p, rng, neighbours)
        for i, j in pairs:
            tape = bytearray(soup[i].tobytes() + soup[j].tobytes())
            a, b = show(tape[:TAPE]), show(tape[TAPE:])
            ops += run(tape, p.max_steps)
            soup[i] = np.frombuffer(tape, dtype=np.uint8, count=TAPE)
            soup[j] = np.frombuffer(tape, dtype=np.uint8, count=TAPE, offset=TAPE)
            a2, b2 = show(tape[:TAPE]), show(tape[TAPE:])
            if Counter((a, b)) != Counter((a2, b2)):
                record((a, b), (a2, b2), "execution")
        yield frame(epoch, ops, len(pairs))

    fired = tally.reactions()
    final = Counter(show(soup[i].tobytes()) for i in range(p.tapes))
    names = dict.fromkeys(initial)
    for lhs, rhs, _ in fired:
        names.update(dict.fromkeys((*lhs, *rhs)))
    names.update(dict.fromkeys(final))

    extras = {
        "reaction_kinds": [kinds[_key(lhs, rhs)] for lhs, rhs, _ in fired],
        "final_state": dict(final.most_common()),
        "notation": ("species are 64-byte programs, one character per byte: BFF commands "
                     "<>{}-+.,[] as themselves, the zero byte as '0', any other byte b as U+0100+b"),
    }
    if neighbours is not None:
        extras["space"] = {"type": "grid", "width": p.width, "height": p.tapes // p.width,
                           "radius": 2, "final_grid": [show(soup[i].tobytes()) for i in range(p.tapes)]}
    return Network(
        species=[Species(s) for s in names],
        reactions=[Reaction.of(lhs, rhs, count=count) for lhs, rhs, count in fired],
        status="observed",
        initial_state={s: float(c) for s, c in initial.items()},
        extras=extras,
    )
