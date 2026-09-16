"""Core War / Redcode: Dewdney (1984), book 10.6.2 (assembler automata).

Catalog id: corewar.

Warriors written in Redcode share a cyclic core of M instructions. Each
warrior owns a FIFO queue of tasks (program counters); every cycle MARS
executes one instruction of each living warrior, in turn. DAT (and division by
zero) removes the task; SPL adds one, up to the task limit. A warrior with no
task left is dead; a round ends when at most one warrior lives or after
max_cycles cycles (a tie).

The machine follows the ICWS'94 draft standard (v3.3, 1995) as implemented by
pMARS 0.9.2 (the reference MARS of the KotH hills, built with EXT94: SEQ, SNE,
NOP, the A-number modes * { } and P-space LDP/STP), ported instruction by
instruction from pMARS sim.c: in-register operand evaluation, the pMARS cycle
counter (a dead warrior's share of the remaining cycles is removed), warrior 1
loaded at address 0, the round starter rotating over the warriors, the last
round result in P-space cell 0. The assembler covers the ICWS'94 assembly
format (labels, EQU text substitution, ORG, END, expressions, predefined
labels such as CORESIZE) with pMARS's defaults for missing operands and
modifiers; FOR/ROF and PIN are not supported.

Chemical reading. Molecules are the instructions occupying core cells (species
id = the instruction in load-file form without spaces, e.g. "MOV.I$0,$1") and
the tasks of each warrior (species id = warrior name + a hash of its assembled
code; its count is the number of tasks). One executed instruction is one
reaction event: the reactants are the executing task and the instructions in
every cell the step reads or writes (the executing cell, pointer cells of
indirect operands, the A-cell when its value is used, the B-cell when it is
read or written), the products are the task(s) that go back into the queue (0
for DAT, 2 for a successful SPL) and the same cells after the step. So the Imp
replicates, P_imp + MOV.I$0,$1 + DAT.F$0,$0 -> P_imp + 2 MOV.I$0,$1; a bomb
destroys code; a task executing foreign code is infected; the number of cells
is conserved. The network aggregates the events of a battle (status observed).
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter, deque

from chemart.network import Network, Reaction, Species

OPCODES = ("DAT", "MOV", "ADD", "SUB", "MUL", "DIV", "MOD", "JMP", "JMZ", "JMN", "DJN",
           "CMP", "SEQ", "SNE", "SLT", "SPL", "NOP", "LDP", "STP")
(DAT, MOV, ADD, SUB, MUL, DIV, MOD, JMP, JMZ, JMN, DJN,
 CMP, SEQ, SNE, SLT, SPL, NOP, LDP, STP) = range(len(OPCODES))
MODIFIERS = ("A", "B", "AB", "BA", "F", "X", "I")
mA, mB, mAB, mBA, mF, mX, mI = range(7)
#: Addressing modes: immediate, direct, B-number indirect / predecrement /
#: postincrement, A-number indirect / predecrement / postincrement.
MODES = "#$@<>*{}"
IMM, DIR, BIND, BDEC, BINC, AIND, ADEC, AINC = range(8)

#: The initial instruction of the KOTH variable set (ICWS'94 sec. 4.3, pMARS INITIALINST).
EMPTY = (DAT, mF, DIR, 0, DIR, 0)
_READS_A = frozenset({MOV, ADD, SUB, MUL, DIV, MOD, CMP, SEQ, SNE, SLT, LDP, STP})
_USES_B = frozenset(range(len(OPCODES))) - {DAT, JMP, SPL, NOP}

#: Published warriors, in ICWS'94 assembly.
WARRIORS = {
    # Dewdney (1984): MOV 0 1. Karonen's guide writes it MOV.I 0, 1.
    "imp": ";name Imp\n;author A. K. Dewdney\n        MOV.I  $0, $1\n",
    # Dewdney (1984), in the ICWS'88/'94 form of Karonen's guide (DAT bomb every 4th cell).
    "dwarf": (";name Dwarf\n;author A. K. Dewdney\n"
              "        ADD.AB #4, 3\n        MOV.I  2, @2\n        JMP    -2\n        DAT    #0, #0\n"),
    # Book sec. 10.6.2: SPL 2 / JMP -1 / MOV 0 1, "an avalanche of self-replicating MOV 0 1 programs".
    "imp-avalanche": ";name Imp avalanche\n;author book 10.6.2\n        SPL 2\n        JMP -1\n        MOV 0, 1\n",
}


# --- instructions ----------------------------------------------------------------
def signed(v: int, M: int) -> int:
    return v - M if v > M // 2 else v


def instruction_text(cell, M: int) -> str:
    op, mod, am, av, bm, bv = cell
    return f"{OPCODES[op]}.{MODIFIERS[mod]} {MODES[am]}{signed(av, M)}, {MODES[bm]}{signed(bv, M)}"


def instruction_id(cell, M: int) -> str:
    return instruction_text(cell, M).replace(" ", "")


class Warrior:
    def __init__(self, name: str, author: str, code: list, start: int, M: int):
        self.name, self.author, self.code, self.start = name, author, code, start
        self.load_file = "\n".join([f"ORG {start}"] + [instruction_text(c, M) for c in code]) + "\n"
        self.hash = hashlib.sha1(self.load_file.encode()).hexdigest()[:8]
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_") or "warrior"
        self.id = f"{safe}_{self.hash}"


# --- assembler ---------------------------------------------------------------------
_EXPR_TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z_]\w*)|(&&|\|\||==|!=|<=|>=|[-+*/%()<>!]))")
_IDENT = re.compile(r"[A-Za-z_]\w*")
_BINARY = [("||",), ("&&",), ("==", "!="), ("<", "<=", ">", ">="), ("+", "-"), ("*", "/", "%")]


def _tokens(text: str) -> list:
    out, pos, text = [], 0, text.rstrip()
    while pos < len(text):
        m = _EXPR_TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot parse expression {text!r}")
        out.append(int(m.group(1)) if m.group(1) else m.group(2) or m.group(3))
        pos = m.end()
    return out


def evaluate(text: str, names) -> int:
    """Integer expression with C operators and truncating division; names(identifier) -> int."""
    toks = _tokens(text)
    if not toks:
        raise ValueError("empty expression")
    pos = 0

    def unary():
        nonlocal pos
        if pos >= len(toks):
            raise ValueError(f"incomplete expression {text!r}")
        t = toks[pos]
        pos += 1
        if t == "-":
            return -unary()
        if t == "+":
            return unary()
        if t == "!":
            return int(not unary())
        if t == "(":
            v = binary(0)
            if pos >= len(toks) or toks[pos] != ")":
                raise ValueError(f"unbalanced parentheses in {text!r}")
            pos += 1
            return v
        if isinstance(t, int):
            return t
        if isinstance(t, str) and _IDENT.fullmatch(t):
            return names(t)
        raise ValueError(f"unexpected {t!r} in expression {text!r}")

    def binary(level):
        nonlocal pos
        if level == len(_BINARY):
            return unary()
        v = binary(level + 1)
        while pos < len(toks) and toks[pos] in _BINARY[level]:
            o = toks[pos]
            pos += 1
            r = binary(level + 1)
            if o in ("/", "%"):
                if r == 0:
                    raise ValueError(f"division by zero in {text!r}")
                q = abs(v) // abs(r) * (1 if (v >= 0) == (r >= 0) else -1)
                v = q if o == "/" else v - q * r
            else:
                v = {"+": v + r, "-": v - r, "*": v * r, "==": int(v == r), "!=": int(v != r),
                     "<": int(v < r), "<=": int(v <= r), ">": int(v > r), ">=": int(v >= r),
                     "&&": int(bool(v) and bool(r)), "||": int(bool(v) or bool(r))}[o]
        return v

    v = binary(0)
    if pos != len(toks):
        raise ValueError(f"trailing tokens in expression {text!r}")
    return v


def _split_operands(text: str) -> list[str]:
    depth, parts, cur = 0, [], []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts]


def _default_modifier(op: int, am: int, bm: int) -> int:
    """pMARS asm.c (= ICWS'94 app. A.2.1.2, the ICWS'88 conversion)."""
    if op in (DAT, NOP):
        return mF
    if op in (MOV, CMP, SEQ, SNE):
        return mAB if am == IMM else mB if bm == IMM else mI
    if op in (ADD, SUB, MUL, DIV, MOD):
        return mAB if am == IMM else mB if bm == IMM else mF
    if op in (SLT, LDP, STP):
        return mAB if am == IMM else mB
    return mB


def assemble(source: str, core_size: int = 8000, max_length: int = 100, constants: dict | None = None) -> Warrior:
    """Assemble ICWS'94 Redcode (or a load file) into a Warrior; numbers are taken modulo core_size."""
    M = core_size
    name, author = "", ""
    labels: dict[str, int] = {}
    equs: dict[str, str] = {}
    lines: list = []
    pending: list[str] = []
    org = None
    for lineno, raw in enumerate(source.splitlines(), 1):
        text = raw.strip()
        if text.startswith(";"):
            m = re.match(r";\s*(name|author)\b\s*(.*)", text, re.I)
            if m:
                if m.group(1).lower() == "name":
                    name = m.group(2).strip()
                else:
                    author = m.group(2).strip()
            continue
        text = text.split(";", 1)[0].strip()
        ended = False
        while text:
            head = re.split(r"\s", text, maxsplit=1)[0]
            rest = text[len(head):].strip()
            word = head.rstrip(":")
            base, _, modifier = word.upper().partition(".")
            if base in OPCODES and (not modifier or modifier in MODIFIERS):
                for label in pending:
                    labels[label] = len(lines)
                pending = []
                lines.append((OPCODES.index(base), MODIFIERS.index(modifier) if modifier else None, rest, lineno))
                break
            if base in ("EQU", "ORG", "END") and not modifier:
                if base == "EQU":
                    if not pending:
                        raise ValueError(f"line {lineno}: EQU without a label")
                    for label in pending:
                        equs[label] = rest
                    pending = []
                elif rest:
                    org = (rest, lineno)
                ended = base == "END"
                break
            if base in ("FOR", "ROF", "PIN"):
                raise ValueError(f"line {lineno}: {base} is not supported by this assembler")
            if not _IDENT.fullmatch(word):
                raise ValueError(f"line {lineno}: cannot parse {raw.strip()!r}")
            if word in labels or word in equs or word in pending:
                raise ValueError(f"line {lineno}: duplicate label {word!r}")
            pending.append(word)
            text = rest
        if ended:
            break
    if not lines:
        raise ValueError("the warrior has no instructions")
    if len(lines) > max_length:
        raise ValueError(f"the warrior has {len(lines)} instructions, more than max_length = {max_length}")
    predefined = {"CORESIZE": M, **(constants or {})}

    def substitute(text: str) -> str:
        for _ in range(100):
            new = _IDENT.sub(lambda m: f"{equs[m.group(0)]}" if m.group(0) in equs else m.group(0), text)
            if new == text:
                return text
            text = new
        raise ValueError(f"EQU substitution does not terminate in {text!r}")

    def names_at(index: int):
        def names(ident: str) -> int:
            if ident in labels:
                return labels[ident] - index
            if ident == "CURLINE":
                return index
            if ident in predefined:
                return predefined[ident]
            raise ValueError(f"unknown label {ident!r}")
        return names

    code = []
    for index, (op, mod, rest, lineno) in enumerate(lines):
        operands = [o for o in _split_operands(substitute(rest))] if rest.strip() else []
        if len(operands) > 2 or any(not o for o in operands):
            raise ValueError(f"line {lineno}: expected one or two operands, got {rest!r}")
        fields = []
        for o in operands:
            mode = MODES.index(o[0]) if o[0] in MODES else DIR
            expr = o[1:] if o[0] in MODES else o
            try:
                fields.append((mode, evaluate(expr, names_at(index)) % M))
            except ValueError as e:
                raise ValueError(f"line {lineno}: {e}") from None
        if not fields:
            raise ValueError(f"line {lineno}: {OPCODES[op]} needs an operand")
        if len(fields) == 1:
            if op == DAT:                       # pMARS: DAT x -> DAT #0, x
                fields = [(IMM, 0), fields[0]]
            elif op in (JMP, SPL, NOP):         # pMARS: JMP x -> JMP x, $0
                fields.append((DIR, 0))
            else:
                raise ValueError(f"line {lineno}: {OPCODES[op]} needs two operands")
        (am, av), (bm, bv) = fields
        code.append((op, _default_modifier(op, am, bm) if mod is None else mod, am, av, bm, bv))
    start = 0
    if org is not None:
        start = evaluate(substitute(org[0]), names_at(0)) % M
    return Warrior(name, author, code, start, M)


# --- MARS --------------------------------------------------------------------------------
class Mars:
    """A pMARS-compatible battle: core, task queues, P-space, per-round loop."""

    def __init__(self, warriors: list[Warrior], core_size=8000, max_cycles=80000, max_processes=8000,
                 pspace_size=None):
        self.warriors, self.M = warriors, core_size
        self.max_cycles, self.max_processes = max_cycles, max_processes
        if pspace_size is None:                 # clparse.c: coreSize / i for the largest i <= 16 dividing it
            pspace_size = next(core_size // i for i in range(16, 0, -1) if core_size % i == 0)
        self.pspace_size = pspace_size
        n = len(warriors)
        self.pspace = [[0] * pspace_size for _ in range(n)]
        self.last_result = [core_size - 1] * n  # pmars.c pspace_init: "unlikely number"
        self.core: list = []
        self.owner: list = []

    def load(self, positions):
        M = self.M
        self.core = core = [EMPTY] * M
        self.owner = owner = [-1] * M
        for i, w in enumerate(self.warriors):
            for k, cell in enumerate(w.code):
                a = (positions[i] + k) % M
                core[a], owner[a] = cell, i

    def play(self, positions, starter=0, events=None, stats=None, series=None, every=1) -> dict:
        """One round (pMARS simulator1 loop). events: Counter of (warrior, tasks out, cells before, cells after)."""
        M, n = self.M, len(self.warriors)
        self.load(positions)
        core, owner = self.core, self.owner
        queues = [deque([(positions[i] + w.start) % M]) for i, w in enumerate(self.warriors)]
        nxt = [(i + 1) % n for i in range(n)]
        w, oldw = starter, (starter - 1) % n
        cycle, left, steps = n * self.max_cycles, n, 0
        maxp, S, pspace, lastres = self.max_processes, self.pspace_size, self.pspace, self.last_result
        death = [None] * n
        record = events is not None
        while True:
            q = queues[w]
            pc = q.popleft()
            steps += 1
            op, mod, am, av, bm, bv = core[pc]
            old = {}
            touched = [pc]
            # A-operand
            if am == IMM:
                addrA, AA, AB = pc, av, bv
            else:
                addrA = pc + av
                if addrA >= M:
                    addrA -= M
                if am == DIR:
                    c = core[addrA]
                    AA, AB = c[3], c[5]
                else:
                    p = addrA
                    touched.append(p)
                    c = core[p]
                    f = 3 if am >= AIND else 5
                    t = c[f]
                    if am == BDEC or am == ADEC:
                        t = t - 1 if t else M - 1
                        old[p] = c
                        c = core[p] = (c[0], c[1], c[2], t, c[4], c[5]) if f == 3 else c[:5] + (t,)
                    addrA = p + t
                    if addrA >= M:
                        addrA -= M
                    d = core[addrA]
                    AA, AB = d[3], d[5]
                    if am == BINC or am == AINC:
                        t = t + 1 if t + 1 < M else 0
                        c = core[p]
                        old.setdefault(p, c)
                        core[p] = (c[0], c[1], c[2], t, c[4], c[5]) if f == 3 else c[:5] + (t,)
            # B-operand
            if bm == IMM:
                addrB, BA, BB = pc, av, bv
            else:
                addrB = pc + bv
                if addrB >= M:
                    addrB -= M
                if bm == DIR:
                    c = core[addrB]
                    BA, BB = c[3], c[5]
                else:
                    p = addrB
                    touched.append(p)
                    c = core[p]
                    f = 3 if bm >= AIND else 5
                    t = c[f]
                    if bm == BDEC or bm == ADEC:
                        t = t - 1 if t else M - 1
                        old.setdefault(p, c)
                        c = core[p] = (c[0], c[1], c[2], t, c[4], c[5]) if f == 3 else c[:5] + (t,)
                    addrB = p + t
                    if addrB >= M:
                        addrB -= M
                    d = core[addrB]
                    BA, BB = d[3], d[5]
                    if bm == BINC or bm == AINC:
                        t = t + 1 if t + 1 < M else 0
                        c = core[p]
                        old.setdefault(p, c)
                        core[p] = (c[0], c[1], c[2], t, c[4], c[5]) if f == 3 else c[:5] + (t,)
            if op in _READS_A:
                touched.append(addrA)
            if op in _USES_B:
                touched.append(addrB)
            # execute
            nxtpc = pc + 1 if pc + 1 < M else 0
            die, out = False, 1
            newA = newB = None
            if op == MOV:
                if mod == mI:
                    s = core[addrA]
                    old.setdefault(addrB, core[addrB])
                    core[addrB] = (s[0], s[1], s[2], AA, s[4], AB)
                elif mod == mA:
                    newA = AA
                elif mod == mB:
                    newB = AB
                elif mod == mAB:
                    newB = AA
                elif mod == mBA:
                    newA = AB
                elif mod == mF:
                    newA, newB = AA, AB
                else:
                    newA, newB = AB, AA
                q.append(nxtpc)
            elif op == ADD or op == SUB or op == MUL:
                if op == ADD:
                    def calc(x, y):
                        return (x + y) % M
                elif op == SUB:
                    def calc(x, y):
                        return (x - y) % M
                else:
                    def calc(x, y):
                        return (x * y) % M
                if mod == mA:
                    newA = calc(BA, AA)
                elif mod == mB:
                    newB = calc(BB, AB)
                elif mod == mAB:
                    newB = calc(BB, AA)
                elif mod == mBA:
                    newA = calc(BA, AB)
                elif mod == mX:
                    newB, newA = calc(BB, AA), calc(BA, AB)
                else:                                   # F and I
                    newA, newB = calc(BA, AA), calc(BB, AB)
                q.append(nxtpc)
            elif op == DIV or op == MOD:
                if op == DIV:
                    def calc(x, y):
                        return x // y
                else:
                    def calc(x, y):
                        return x % y
                if mod == mA:
                    if AA:
                        newA = calc(BA, AA)
                    else:
                        die = True
                elif mod == mB:
                    if AB:
                        newB = calc(BB, AB)
                    else:
                        die = True
                elif mod == mAB:
                    if AA:
                        newB = calc(BB, AA)
                    else:
                        die = True
                elif mod == mBA:
                    if AB:
                        newA = calc(BA, AB)
                    else:
                        die = True
                elif mod == mX:
                    if AB:
                        newA = calc(BA, AB)
                    if AA:
                        newB = calc(BB, AA)
                    die = not (AA and AB)
                else:                                   # F and I
                    if AA:
                        newA = calc(BA, AA)
                    if AB:
                        newB = calc(BB, AB)
                    die = not (AA and AB)
                if not die:
                    q.append(nxtpc)
            elif op == JMP:
                q.append(addrA)
            elif op == JMZ or op == JMN:
                if mod == mA or mod == mBA:
                    zero = not BA
                elif mod == mB or mod == mAB:
                    zero = not BB
                else:
                    zero = not BA and not BB
                q.append(addrA if zero == (op == JMZ) else nxtpc)
            elif op == DJN:
                c = core[addrB]
                if mod == mA or mod == mBA:
                    newA = c[3] - 1 if c[3] else M - 1
                    jump = BA != 1
                elif mod == mB or mod == mAB:
                    newB = c[5] - 1 if c[5] else M - 1
                    jump = BB != 1
                else:
                    newA = c[3] - 1 if c[3] else M - 1
                    newB = c[5] - 1 if c[5] else M - 1
                    jump = not (BA == 1 and BB == 1)
                q.append(addrA if jump else nxtpc)
            elif op == CMP or op == SEQ or op == SNE:
                if mod == mA:
                    equal = BA == AA
                elif mod == mB:
                    equal = BB == AB
                elif mod == mAB:
                    equal = BB == AA
                elif mod == mBA:
                    equal = BA == AB
                elif mod == mF:
                    equal = BA == AA and BB == AB
                elif mod == mX:
                    equal = BB == AA and BA == AB
                else:
                    x, y = core[addrA], core[addrB]
                    equal = x[0] == y[0] and x[1] == y[1] and x[2] == y[2] and x[4] == y[4] \
                        and BA == AA and BB == AB
                skip = equal != (op == SNE)
                q.append((pc + 2) % M if skip else nxtpc)
            elif op == SLT:
                if mod == mA:
                    less = AA < BA
                elif mod == mB:
                    less = AB < BB
                elif mod == mAB:
                    less = AA < BB
                elif mod == mBA:
                    less = AB < BA
                elif mod == mX:
                    less = AA < BB and AB < BA
                else:
                    less = AA < BA and AB < BB
                q.append((pc + 2) % M if less else nxtpc)
            elif op == SPL:
                q.append(nxtpc)
                if len(q) < maxp:
                    q.append(addrA)
                    out = 2
            elif op == NOP:
                q.append(nxtpc)
            elif op == LDP or op == STP:
                if op == LDP:
                    if mod == mA:
                        idx = AA
                    elif mod == mAB:
                        idx = AA
                    else:
                        idx = AB
                    i = idx % S
                    value = pspace[w][i] if i else lastres[w]
                    if mod == mA or mod == mBA:
                        newA = value
                    else:
                        newB = value
                else:
                    if mod == mA:
                        idx, value = BA, AA
                    elif mod == mAB:
                        idx, value = BB, AA
                    elif mod == mBA:
                        idx, value = BA, AB
                    else:
                        idx, value = BB, AB
                    i = idx % S
                    if i:
                        pspace[w][i] = value
                    else:
                        lastres[w] = value
                q.append(nxtpc)
            else:                                       # DAT
                die = True
            if newA is not None or newB is not None:
                c = core[addrB]
                old.setdefault(addrB, c)
                core[addrB] = (c[0], c[1], c[2], c[3] if newA is None else newA, c[4],
                               c[5] if newB is None else newB)
            if die:
                out = 0
            # observe
            if record:
                cells = set(touched)
                cells.update(old)
                before = tuple(sorted(old.get(a, core[a]) for a in cells))
                after = tuple(sorted(core[a] for a in cells))
                events[(w, out, before, after)] += 1
            if stats is not None:
                st = stats[w]
                st["executed"] += 1
                if owner[pc] >= 0 and owner[pc] != w:
                    st["foreign_executions"] += 1
                if out == 2:
                    st["splits"] += 1
                elif out == 0:
                    st["task_deaths"] += 1
                for a, c in old.items():
                    if core[a] != c:
                        st["writes"] += 1
                        if owner[a] >= 0 and owner[a] != w:
                            st["foreign_writes"] += 1
                        owner[a] = w
                if op == MOV and mod == mI and old.get(addrB, core[addrB]) != core[addrB]:
                    st["copies"] += 1
            if series is not None and steps % every == 0:
                series["step"].append(steps)
                for i in range(n):
                    series["tasks"][i].append(len(queues[i]))
            # schedule (pMARS: nopush / die / continue while (--cycle))
            if die and not q:
                death[w] = steps
                cycle = cycle - 1 - (cycle - 1) // left
                left -= 1
                if left < 2:
                    break
                nxt[oldw] = nxt[w]
                w = nxt[w]
            else:
                oldw, w = w, nxt[w]
            cycle -= 1
            if cycle == 0:
                break
        alive = [i for i in range(n) if queues[i]]
        for i in range(n):
            lastres[i] = len(alive) if queues[i] else 0
        return {"steps": steps, "survivors": alive, "tasks": [len(qq) for qq in queues], "death_step": death}


def place(n: int, core_size: int, separation: int, rng) -> list[int]:
    """Random load addresses; warrior 1 at 0 (pMARS simulator1, posit and npos with rng for pMARS's generator)."""
    if n == 1:
        return [0]
    if n == 2:
        return [0, separation + int(rng.integers(core_size + 1 - 2 * separation))]
    pos = [0] * n
    k, retries1, retries2 = 1, 20, 4
    while k < n:
        pos[k] = int(rng.integers(core_size - 2 * separation + 1)) + separation
        clash = next((i for i in range(1, k) if abs(pos[k] - pos[i]) < separation), None)
        if clash is None:
            k += 1
        elif not retries2:
            break
        elif not retries1:
            k, retries2, retries1 = clash, retries2 - 1, 20
        else:
            retries1 -= 1
    else:
        return pos
    room = core_size - separation * n + 1                # npos
    for i in range(1, n):
        temp = int(rng.integers(room))
        j = i - 1
        while j > 0 and not temp > pos[j]:
            pos[j + 1] = pos[j]
            j -= 1
        pos[j + 1] = temp
    for i in range(1, n):
        pos[i] += separation * i
    for i in range(1, n):
        j = int(rng.integers(n - i)) + i
        pos[i], pos[j] = pos[j], pos[i]
    return pos


# --- generate ------------------------------------------------------------------------------
_STATS = ("executed", "splits", "task_deaths", "foreign_executions", "writes", "foreign_writes", "copies")


def generate(p, rng):
    M = p.core_size
    if not isinstance(p.warriors, list) or not p.warriors or not all(isinstance(s, str) for s in p.warriors):
        raise ValueError("warriors must be a non-empty list of strings: built-in names "
                         f"{sorted(WARRIORS)} or Redcode sources")
    n = len(p.warriors)
    if p.min_distance < p.max_length:
        raise ValueError(f"min_distance ({p.min_distance}) must be >= max_length ({p.max_length}) as in pMARS")
    if M < n * p.min_distance:
        raise ValueError(f"core_size ({M}) must be >= number of warriors x min_distance ({n * p.min_distance})")
    if p.pspace_size > M:
        raise ValueError(f"pspace_size ({p.pspace_size}) must not exceed core_size ({M})")
    constants = {"MAXPROCESSES": p.max_processes, "MAXCYCLES": p.max_cycles, "MAXLENGTH": p.max_length,
                 "MINDISTANCE": p.min_distance, "ROUNDS": p.rounds, "WARRIORS": n}
    pspace = p.pspace_size or None
    if pspace:
        constants["PSPACESIZE"] = pspace
    warriors = []
    for k, src in enumerate(p.warriors):
        text = WARRIORS.get(src.strip().lower(), src) if "\n" not in src else src
        if text is src and "\n" not in src and not re.search(r"[\s,]", src.strip()):
            raise ValueError(f"warriors[{k}]: {src!r} is neither a built-in warrior {sorted(WARRIORS)} "
                             "nor Redcode source")
        try:
            warriors.append(assemble(text, M, p.max_length, constants))
        except ValueError as e:
            raise ValueError(f"warriors[{k}]: {e}") from None
    ids = []
    for k, wr in enumerate(warriors):
        wid = wr.id if wr.id not in ids else f"{wr.id}_{k + 1}"
        ids.append(wid)
    fixed = p.positions
    if fixed:
        if (not isinstance(fixed, list) or len(fixed) != n - 1
                or not all(isinstance(x, int) and not isinstance(x, bool) for x in fixed)):
            raise ValueError(f"positions must be empty (random) or a list of {n - 1} integer load addresses "
                             "for warriors 2..n (warrior 1 is loaded at 0)")
        pos = [0, *fixed]
        for i in range(1, n):
            if not p.min_distance <= pos[i] <= M - p.min_distance:
                raise ValueError(f"positions: {pos[i]} must lie in [min_distance, core_size - min_distance]")
            if any(abs(pos[i] - pos[j]) < p.min_distance for j in range(1, i)):
                raise ValueError("positions: load addresses must be at least min_distance apart")

    mars = Mars(warriors, M, p.max_cycles, p.max_processes, pspace)
    events: Counter = Counter()
    stats = [Counter({k: 0 for k in _STATS}) for _ in range(n)]
    every = n * max(1, p.max_cycles // 200)
    series = {"step": [], "tasks": [[] for _ in range(n)]}
    rounds, final = [], Counter()
    score, wins, ties, losses = [0] * n, [0] * n, [0] * n, [0] * n
    initial = None
    for r in range(p.rounds):
        pos = [0, *fixed] if fixed else place(n, M, p.min_distance, rng)
        starter = r % n
        res = mars.play(pos, starter, events, stats, series if r == 0 else None, every)
        if initial is None:
            initial = Counter(mars.core)
        alive = res["survivors"]
        for i in range(n):
            if i in alive:
                score[i] += (n * n - 1) // len(alive)
                if len(alive) == 1 and n > 1:
                    wins[i] += 1
                else:
                    ties[i] += 1
            else:
                losses[i] += 1
        final.update(mars.core)
        rounds.append({"round": r + 1, "positions": pos, "starter": ids[starter], "steps": res["steps"],
                       "survivors": [ids[i] for i in alive], "tasks": dict(zip(ids, res["tasks"])),
                       "death_step": {ids[i]: s for i, s in enumerate(res["death_step"]) if s is not None}})
        final.update({(-1, i): res["tasks"][i] for i in range(n)})
    initial = Counter({cell: 0 for cell in initial})
    for i, wr in enumerate(warriors):
        for cell in wr.code:
            initial[cell] += 1
    initial[EMPTY] += M - sum(len(wr.code) for wr in warriors)
    initial = +initial

    names: dict = {}

    def sid(cell) -> str:
        s = names.get(cell)
        if s is None:
            s = names[cell] = instruction_id(cell, M)
        return s

    species = [Species(ids[i], structure=warriors[i].load_file) for i in range(n)]
    seen = set()

    def add(cell):
        if cell not in seen:
            seen.add(cell)
            species.append(Species(sid(cell), structure=instruction_text(cell, M)))

    for cell in initial:
        add(cell)
    reactions = []
    for (w, out, before, after), count in events.items():
        for cell in (*before, *after):
            add(cell)
        lhs = [ids[w], *map(sid, before)]
        rhs = [ids[w]] * out + [sid(c) for c in after]
        reactions.append(Reaction.of(lhs, rhs, count=count))
    final_state = {}
    for key, v in final.items():
        if isinstance(key, tuple) and len(key) == 2 and key[0] == -1:
            if v:
                final_state[ids[key[1]]] = v
        else:
            add(key)
            final_state[sid(key)] = final_state.get(sid(key), 0) + v
    initial_state = {ids[i]: 1 for i in range(n)}
    initial_state.update({sid(c): k for c, k in initial.items()})
    analysis = {
        "rounds": rounds,
        "score": dict(zip(ids, score)), "wins": dict(zip(ids, wins)),
        "ties": dict(zip(ids, ties)), "losses": dict(zip(ids, losses)),
        "events": {ids[i]: dict(stats[i]) for i in range(n)},
        "series_round_1": {"step": series["step"], "tasks": dict(zip(ids, series["tasks"]))},
    }
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={k: float(v) for k, v in initial_state.items()},
        extras={
            "analysis": analysis,
            "warriors": [{"id": ids[i], "name": wr.name, "author": wr.author, "length": len(wr.code),
                          "start": wr.start, "load_file": wr.load_file} for i, wr in enumerate(warriors)],
            "final_state": final_state,
            "conservation": [{"name": "core cells", "vector": {s.id: (0 if k < n else 1)
                                                               for k, s in enumerate(species)}}],
            "space": {"dimensions": 1, "shape": [M], "boundary": "periodic",
                      "note": "reactions are position-free; the core is the cyclic memory the events act on"},
            "rules": {"standard": "ICWS'94 draft v3.3 as implemented by pMARS 0.9.2", "core_size": M,
                      "max_cycles": p.max_cycles, "max_processes": p.max_processes,
                      "min_distance": p.min_distance, "max_length": p.max_length,
                      "pspace_size": mars.pspace_size, "rounds": p.rounds},
        },
    )
