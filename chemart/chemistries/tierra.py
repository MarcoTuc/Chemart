"""Tierra (Ray 1991): self-replicating machine code in a shared soup. Catalog id: tierra.

A block of memory (the soup) holds creatures: programs over the 32
instructions of instruction set 0 (Ray 1991, appendix B; book table 10.3),
five bits each and without operands. Every creature has a virtual CPU (ax,
bx, cx, dx, a ten-word circular stack and an instruction pointer).

- Addressing by template: jmp, jmpb, call, adr, adrb and adrf read the run of
  nop_0/nop_1 that follows them and look for the complementary pattern,
  outward in both directions (forward wins ties), backward or forward, up to
  search_limit instructions away (ctemplate in instruct.c). adr* put the
  address just after the found pattern in ax and the template size in cx.
- Memory allocation and protection: mal allocates a daughter block of cx
  instructions (first fit) and puts its address in ax. A creature may write
  only in its own two blocks and in free memory; it can read and execute
  anything. divide makes the daughter block an independent creature with its
  own CPU (instruction pointer at its start), provided at least 70% of it was
  written by mov_iab (divide in instruct.c).
- The slicer: a circular queue; each creature executes size^slice_power
  instructions per turn (or a constant slice); a newborn enters just ahead of
  its mother, so it runs last in the round.
- The reaper: a linear queue; newborns enter at the bottom. An instruction
  that fails (template not found, illegal mal, write to protected memory,
  premature divide) moves the creature one place up if it has at least as
  many errors as the one above; a successful mal or divide moves it one place
  down if it has no more errors than the one below. Whenever allocated
  memory exceeds reaper_threshold of the soup, the creature at the top dies:
  its blocks are freed but its code stays in the soup.
- Mutations: cosmic rays flip a random bit anywhere in the soup, on average
  once per cosmic_ray_interval instructions executed; copy errors flip a bit
  of a copied instruction once per copy_error_interval moves; flaws make an
  instruction's result off by +-1 (one more or one fewer shift, the
  neighbouring bit flipped, a template found one place off...) once per
  flaw_interval instructions. Intervals vary uniformly around their mean,
  as in Tierra, to avoid periodic effects.

The run is returned as the observed network of events between genotypes
(named by size and a base-26 label in order of appearance, e.g. 0080aaa, the
genome as structure):

- replication ``G -> G + D`` (a creature divides; D = the daughter's genome),
- replication using another creature's code ``G + H -> G + H + D`` when the
  dividing creature executed instructions inside the genome of a living
  creature of genotype H since its last division (parasitism, or mutualism
  when both partners do it),
- death by the reaper ``G -> ∅``,
- mutation of a living genome ``G -> G'`` (a cosmic ray or a write into its
  own genome).

Starting from initial_state and applying every event count reproduces the
final population exactly.
"""

from __future__ import annotations

from bisect import bisect_left, insort
from collections import Counter

from chemart.network import Network, Reaction, Species

N_INST = 32
STACK_SIZE = 10          # configur.h
MIN_CELL_SIZE = 12       # soup_in: MinCellSize, MinGenMemSiz
MAX_MAL_MULT = 3         # soup_in.h: MaxMalMult
MOV_PROP_THR_DIV = 0.7   # soup_in: MovPropThrDiv

#: Instruction set 0 in opcode order, with the names of Ray (1991), appendix B.
MNEMONICS = (
    "nop_0", "nop_1", "or1", "shl", "zero", "if_cz", "sub_ab", "sub_ac",
    "inc_a", "inc_b", "dec_c", "inc_c", "push_ax", "push_bx", "push_cx", "push_dx",
    "pop_ax", "pop_bx", "pop_cx", "pop_dx", "jmp", "jmpb", "call", "ret",
    "mov_cd", "mov_ab", "mov_iab", "adr", "adrb", "adrf", "mal", "divide",
)
#: The same opcodes with the names of the Tierra 6.02 genebank (gb0/opcode.map).
MNEMONICS_602 = (
    "nop0", "nop1", "not0", "shl", "zero", "ifz", "subCAB", "subAAC",
    "incA", "incB", "decC", "incC", "pushA", "pushB", "pushC", "pushD",
    "popA", "popB", "popC", "popD", "jmpo", "jmpb", "call", "ret",
    "movDC", "movBA", "movii", "adro", "adrb", "adrf", "mal", "divide",
)
OPCODES = {name: i for names in (MNEMONICS, MNEMONICS_602) for i, name in enumerate(names)}
OPCODES["jmp_b"] = 21    # spelling of book table 10.3

TEMPLATE_OPS = {20: ("o", True), 21: ("b", True), 22: ("o", True),
                27: ("o", False), 28: ("b", False), 29: ("f", False)}


def assemble(listing: str) -> bytes:
    """Genome from whitespace- or comma-separated mnemonics (1991 or 6.02 names)."""
    words = listing.replace(",", " ").split()
    unknown = [w for w in words if w not in OPCODES]
    if unknown:
        raise ValueError(f"unknown Tierra instruction(s) {unknown[:5]}; known: {', '.join(MNEMONICS)}")
    return bytes(OPCODES[w] for w in words)


def disassemble(genome: bytes) -> str:
    return " ".join(MNEMONICS[b] for b in genome)


#: The ancestor 0080aaa (Ray 1991, appendix C; gb0/0080aaa.tie).
ANCESTOR = assemble("""
    nop_1 nop_1 nop_1 nop_1 zero or1 shl shl mov_cd adrb nop_0 nop_0 nop_0 nop_0
    sub_ac mov_ab adrf nop_0 nop_0 nop_0 nop_1 inc_a sub_ab nop_1 nop_1 nop_0 nop_1
    mal call nop_0 nop_0 nop_1 nop_1 divide jmp nop_0 nop_0 nop_1 nop_0 if_cz
    nop_1 nop_1 nop_0 nop_0 push_ax push_bx push_cx nop_1 nop_0 nop_1 nop_0
    mov_iab dec_c if_cz jmp nop_0 nop_1 nop_0 nop_0 inc_a inc_b jmp nop_0 nop_1 nop_0 nop_1
    if_cz nop_1 nop_0 nop_1 nop_1 pop_cx pop_bx pop_ax ret nop_1 nop_1 nop_1 nop_0 if_cz
""")
#: The parasite 0045aaa: the low bit of instruction 42 of the ancestor flipped
#: (copy-procedure template 1100 -> 1110) makes it measure and copy only
#: instructions 0-44 (Ray 1991, "parasites"; gb0/0045aaa.tie).
PARASITE = ANCESTOR[:42] + bytes([1]) + ANCESTOR[43:45]

GENOMES = {
    "0080aaa": ANCESTOR,
    "0045aaa": PARASITE,
    # immune to 0045aaa (Ray 1991, "immunity to parasites")
    "0079aab": assemble("""
        nop1 nop1 nop1 nop1 zero shl shl shl movDC adrb nop0 nop0 nop0 nop0 subAAC movBA
        adrf nop0 nop0 nop0 nop1 movDC subCAB nop1 nop1 nop0 nop1 mal call nop0 nop0 nop1
        nop1 divide call nop0 nop0 nop1 nop0 pushB nop1 nop1 nop0 nop0 pushA pushB pushC
        nop1 nop0 nop1 nop0 movii decC ifz jmpo nop0 nop1 nop0 nop0 incA incB jmpb nop0
        nop1 nop0 nop1 jmpb nop1 nop0 nop1 nop1 popC popB popA ret nop1 nop1 nop1 nop0
    """),
    # the ancestor dissected by hand into two mutualists (Ray 1991, "ECOLOGY", fig. 3)
    "0046aaa": assemble("""
        nop1 nop1 nop1 nop1 adrb nop0 nop0 nop0 nop0 subAAC movBA adrf nop0 nop0 nop0 nop1
        incA subCAB jmpo nop0 nop0 nop1 nop0 ifz nop1 nop1 nop0 nop1 mal call nop0 nop0
        nop1 nop1 divide jmpb nop0 nop0 nop1 nop0 ifz nop1 nop1 nop1 nop0 ifz
    """),
    "0064aaa": assemble("""
        nop1 nop1 nop1 nop1 adrb nop0 nop0 nop0 nop0 subAAC movBA adrf nop0 nop0 nop0 nop1
        incA subCAB jmpo nop0 nop0 nop1 nop0 ifz nop1 nop1 nop0 nop0 pushA pushB pushC
        nop1 nop0 nop1 nop0 movii decC ifz jmpo nop0 nop1 nop0 nop0 incA incB jmpb nop0
        nop1 nop0 nop1 ifz nop1 nop0 nop1 nop1 popC popB popA ret nop1 nop1 nop1 nop0 ifz
    """),
}


def label(index: int) -> str:
    """Base-26 genotype label: 0 -> aaa, 1 -> aab, 26 -> aba."""
    out = ""
    for _ in range(3):
        index, r = divmod(index, 26)
        out = chr(97 + r) + out
    return out if index == 0 else f"{out}{index}"


def _wrap32(x: int) -> int:
    return ((x + 0x80000000) & 0xFFFFFFFF) - 0x80000000


class Cell:
    __slots__ = ("id", "ax", "bx", "cx", "dx", "sp", "st", "ip", "mp", "ms", "dp", "ds",
                 "errors", "genotype", "mov", "mov_min", "mov_max", "inst", "divisions",
                 "last_div_inst", "hosts", "rpos", "alive")

    def __init__(self, cid: int, start: int, size: int, genotype: str):
        self.id = cid
        self.ax = self.bx = self.cx = self.dx = self.sp = 0
        self.st = [0] * STACK_SIZE
        self.ip = start
        self.mp, self.ms = start, size       # mother block (the creature)
        self.dp = self.ds = 0                # daughter block
        self.errors = 0
        self.genotype = genotype
        self.mov = self.mov_min = self.mov_max = 0
        self.inst = 0
        self.divisions = 0
        self.last_div_inst = 0
        self.hosts: Counter = Counter()
        self.rpos = 0
        self.alive = True


class Soup:
    """The Tierra virtual computer with its operating system (slicer, reaper, memory manager)."""

    def __init__(self, rng, soup_size: int, *, search_limit: int = 400, reaper_threshold: float = 0.8,
                 slicer: str = "size-dependent", slice_power: float = 1.0, slice_size: int = 25,
                 cosmic_ray_interval: int = 0, copy_error_interval: int = 0, flaw_interval: int = 0,
                 semantics: str = "ray1991", record: bool = True):
        self.rng = rng
        self.size = soup_size
        self.soup = bytearray(soup_size)         # a new soup is all nop_0
        self.owner = [0] * soup_size             # cell id owning each address (0 = free)
        self.search_limit = search_limit
        self.threshold = reaper_threshold * soup_size
        self.constant_slice = slice_size if slicer == "constant" else 0
        self.slice_power = slice_power
        self.skip_templates = semantics == "tierra602"
        self.intervals = (cosmic_ray_interval, copy_error_interval, flaw_interval)
        self.cr_left = self._draw(cosmic_ray_interval)
        self.cp_left = self._draw(copy_error_interval)
        self.fl_left = self._draw(flaw_interval)
        self.record = record

        self.cells: dict[int, Cell] = {}
        self.next_id = 1
        self.starts: list[int] = []              # allocated blocks, sorted by start
        self.block_size: dict[int, int] = {}
        self.allocated = 0
        self.slicer: list[Cell] = []
        self.cur = 0
        self.reap: list[Cell] = []

        self.genotypes: dict[bytes, str] = {}
        self.genomes: dict[str, bytes] = {}
        self.labels: dict[int, int] = {}
        self.events: dict[tuple, list] = {}
        self.metabolism: dict[str, list[int]] = {}   # instructions for 1st and 2nd replication
        self.metabolism_cell: dict[str, int] = {}    # the creature those were measured on
        self.executed = 0
        self.stats = Counter()

    # ------------------------------------------------------------------ bookkeeping
    def _draw(self, mean: int) -> int:
        """Instructions until the next event: uniform on 1..2*mean (0 = never)."""
        return int(self.rng.integers(1, 2 * mean + 1)) if mean else 0

    def name(self, genome: bytes, preferred: str | None = None) -> str:
        known = self.genotypes.get(genome)
        if known:
            return known
        size = len(genome)
        if preferred and preferred not in self.genomes and preferred[:4] == f"{size:04d}":
            gname = preferred
        else:
            while True:
                i = self.labels.get(size, 0)
                self.labels[size] = i + 1
                gname = f"{size:04d}{label(i)}"
                if gname not in self.genomes:
                    break
        self.genotypes[genome] = gname
        self.genomes[gname] = genome
        return gname

    def event(self, reactants: list[str], products: list[str]) -> None:
        if not self.record:
            return
        lhs, rhs = Counter(reactants), Counter(products)
        key = (frozenset(lhs.items()), frozenset(rhs.items()))
        entry = self.events.get(key)
        if entry is None:
            self.events[key] = entry = [dict(lhs), dict(rhs), 0]
        entry[2] += 1

    def reidentify(self, cell: Cell) -> None:
        new = self.name(bytes(self.soup[cell.mp:cell.mp + cell.ms]))
        if new != cell.genotype:
            self.event([cell.genotype], [new])
            cell.genotype = new

    # ------------------------------------------------------------------ memory
    def first_fit(self, size: int) -> int:
        end = 0
        for s in self.starts:
            if s - end >= size:
                return end
            end = s + self.block_size[s]
        return end if self.size - end >= size else -1

    def claim(self, start: int, size: int, cid: int) -> None:
        insort(self.starts, start)
        self.block_size[start] = size
        self.allocated += size
        self.owner[start:start + size] = [cid] * size

    def release(self, start: int, size: int) -> None:
        del self.starts[bisect_left(self.starts, start)]
        del self.block_size[start]
        self.allocated -= size
        self.owner[start:start + size] = [0] * size

    def inoculate(self, genome: bytes, start: int, preferred: str | None = None) -> Cell:
        if start < 0 or start + len(genome) > self.size or any(self.owner[start:start + len(genome)]):
            raise ValueError(f"cannot place a genome of size {len(genome)} at address {start}")
        cell = Cell(self.next_id, start, len(genome), self.name(genome, preferred))
        self.next_id += 1
        self.soup[start:start + len(genome)] = genome
        self.claim(start, len(genome), cell.id)
        self.cells[cell.id] = cell
        self.slicer.append(cell)
        cell.rpos = len(self.reap)
        self.reap.append(cell)
        return cell

    # ------------------------------------------------------------------ queues
    def up_reaper(self, c: Cell) -> None:
        i = c.rpos
        if i > 0 and c.errors >= self.reap[i - 1].errors:
            other = self.reap[i - 1]
            self.reap[i - 1], self.reap[i] = c, other
            c.rpos, other.rpos = i - 1, i

    def down_reaper(self, c: Cell) -> None:
        i = c.rpos
        if i < len(self.reap) - 1 and c.errors <= self.reap[i + 1].errors:
            other = self.reap[i + 1]
            self.reap[i + 1], self.reap[i] = c, other
            c.rpos, other.rpos = i + 1, i

    def kill(self, cell: Cell) -> None:
        self.event([cell.genotype], [])
        self.stats["deaths"] += 1
        self.release(cell.mp, cell.ms)
        if cell.ds:
            self.release(cell.dp, cell.ds)
        i = self.slicer.index(cell)
        del self.slicer[i]
        if i < self.cur:
            self.cur -= 1
        if self.cur >= len(self.slicer):
            self.cur = 0
        del self.reap[cell.rpos]
        for j in range(cell.rpos, len(self.reap)):
            self.reap[j].rpos = j
        del self.cells[cell.id]
        cell.alive = False

    def reap_for(self, active: Cell) -> bool:
        """Kill the top of the reaper queue (the next one if that is the active cell)."""
        if len(self.cells) <= 1:
            return False
        victim = self.reap[0] if self.reap[0] is not active else self.reap[1]
        self.kill(victim)
        return True

    # ------------------------------------------------------------------ OS calls
    def mal(self, c: Cell, request: int, flaw: int) -> int:
        """malchm/mal: allocate and protect a daughter block; the address, or -1 on error."""
        if request < MIN_CELL_SIZE or request >= self.size:
            return -1
        if request == c.ds or request > MAX_MAL_MULT * c.ms:
            return -1
        size = request + flaw
        if c.ds:
            self.release(c.dp, c.ds)
            c.dp = c.ds = c.mov = c.mov_min = c.mov_max = 0
        while True:
            start = self.first_fit(size)
            if start >= 0:
                break
            if not self.reap_for(c):
                return -1
        self.claim(start, size, c.id)
        c.dp, c.ds = start, size
        self.down_reaper(c)
        return start

    def divide(self, c: Cell, inst: int) -> bool:
        ds = c.ds
        thresh = int(ds * MOV_PROP_THR_DIV)
        span = c.mov_max - c.mov_min + 1 if c.mov else 0
        if ds < MIN_CELL_SIZE or span < MIN_CELL_SIZE or span < thresh or c.mov < thresh:
            return False
        genome = bytes(self.soup[c.dp:c.dp + ds])
        d = Cell(self.next_id, c.dp, ds, self.name(genome))
        self.next_id += 1
        self.owner[c.dp:c.dp + ds] = [d.id] * ds
        self.cells[d.id] = d
        if c.hosts:
            host = c.hosts.most_common(1)[0][0]
            self.event([c.genotype, host], [c.genotype, host, d.genotype])
            self.stats["parasitic_births"] += 1
        else:
            self.event([c.genotype], [c.genotype, d.genotype])
        self.stats["births"] += 1
        c.divisions += 1
        if c.divisions <= 2 and not c.hosts:
            # metabolic data of the first creature of each genotype to replicate on its own
            if self.metabolism_cell.setdefault(c.genotype, c.id) == c.id:
                self.metabolism.setdefault(c.genotype, []).append(inst - c.last_div_inst)
        c.last_div_inst = inst
        c.hosts = Counter()
        c.dp = c.ds = c.mov = c.mov_min = c.mov_max = 0
        # slicer: just ahead of the mother; reaper: bottom
        i = self.slicer.index(c)
        self.slicer.insert(i, d)
        if i <= self.cur:
            self.cur += 1
        d.rpos = len(self.reap)
        self.reap.append(d)
        self.down_reaper(c)
        return True

    def cosmic_ray(self) -> None:
        addr = int(self.rng.integers(self.size))
        self.soup[addr] ^= 1 << int(self.rng.integers(5))
        self.stats["cosmic_rays"] += 1
        o = self.owner[addr]
        if o:
            cell = self.cells[o]
            if cell.mp <= addr < cell.mp + cell.ms:
                self.reidentify(cell)

    # ------------------------------------------------------------------ templates
    def template_size(self, ip: int) -> int:
        soup, size = self.soup, self.size
        s = 0
        while s < size and soup[(ip + 1 + s) % size] < 2:
            s += 1
        return s

    def _window(self, start: int, length: int) -> bytes:
        size = self.size
        start %= size
        if start + length <= size:
            return self.soup[start:start + length]
        return self.soup[start:] + self.soup[:start + length - size]

    def search(self, ip: int, s: int, direction: str, flaw: int = 0) -> int:
        """ctemplate: address after the nearest complement of the template after ip, or -1.

        Forward candidates start at ip + s + 2, backward ones at ip - s; both
        directions advance together for search_limit steps; forward wins ties.
        """
        soup, size = self.soup, self.size
        pattern = bytes(1 - soup[(ip + 1 + i) % size] for i in range(s))
        limit = min(self.search_limit, size - s + 1)
        if limit < 1:
            return -1
        best = None
        if direction != "b":
            f0 = ip + s + 2
            k = self._window(f0, limit - 1 + s).find(pattern)
            if k >= 0:
                best = (k, f0 + k)
        if direction != "f":
            b0 = ip - s
            idx = self._window(b0 - limit + 1, limit - 1 + s).rfind(pattern)
            if idx >= 0:
                k = limit - 1 - idx
                if best is None or k < best[0]:
                    best = (k, b0 - k)
        if best is None:
            return -1
        return (best[1] + flaw + s) % size

    # ------------------------------------------------------------------ the CPU
    def run_slice(self, c: Cell, n: int) -> None:
        soup, size, owner, rng = self.soup, self.size, self.owner, self.rng
        skip = self.skip_templates
        cr_mean, cp_mean, fl_mean = self.intervals
        ax, bx, cx, dx, sp, ip, st = c.ax, c.bx, c.cx, c.dx, c.sp, c.ip, c.st
        cid, mp, me = c.id, c.mp, c.mp + c.ms
        errors_before = c.errors
        done = 0
        while done < n:
            op = soup[ip]
            flaw = 0
            if self.fl_left:
                self.fl_left -= 1
                if not self.fl_left:
                    flaw = 1 if rng.integers(2) else -1
                    self.fl_left = self._draw(fl_mean)
                    self.stats["flaws"] += 1
            err = False
            nxt = ip + 1
            if not mp <= ip < me:                            # executing someone else's code?
                o = owner[ip]
                if o and o != cid:
                    host = self.cells[o]
                    if host.mp <= ip < host.mp + host.ms:
                        c.hosts[host.genotype] += 1
            if op < 2:
                pass
            elif op == 26:                                   # mov_iab: [ax] = [bx]
                dst, src = ax % size, (bx + flaw) % size
                o = owner[dst]
                if dst != src and (o == 0 or o == cid):
                    soup[dst] = soup[src]
                    if self.cp_left:
                        self.cp_left -= 1
                        if not self.cp_left:
                            soup[dst] ^= 1 << int(rng.integers(5))
                            self.cp_left = self._draw(cp_mean)
                            self.stats["copy_errors"] += 1
                    if c.ds and c.dp <= dst < c.dp + c.ds:
                        off = dst - c.dp
                        if not c.mov:
                            c.mov_min = c.mov_max = off
                        elif off < c.mov_min:
                            c.mov_min = off
                        elif off > c.mov_max:
                            c.mov_max = off
                        c.mov += 1
                    elif mp <= dst < me:
                        self.reidentify(c)
                else:
                    err = True
            elif op == 10:                                   # dec_c
                cx = cx + flaw - 1
            elif op == 5:                                    # if_cz
                if cx + flaw != 0:
                    nxt = ip + 2
            elif op == 8:                                    # inc_a
                ax = ax + flaw + 1
            elif op == 9:                                    # inc_b
                bx = bx + flaw + 1
            elif op in TEMPLATE_OPS:
                direction, is_jump = TEMPLATE_OPS[op]
                s = self.template_size(ip)
                if s == 0:
                    if not skip:
                        err = True                           # no template: error, ignored
                    elif op == 22:                           # 6.02 tcall: push ip + 1
                        sp = (sp + 1) % STACK_SIZE
                        st[sp] = ip + 1 + flaw
                    elif is_jump:                            # 6.02 jmpo/jmpb: ip = bx
                        nxt = bx + flaw
                else:
                    target = self.search(ip, s, direction, flaw)
                    after = ip + s + 1 if skip else ip + 1
                    if target < 0:
                        err = True
                        nxt = after
                    elif op == 22:                           # call
                        sp = (sp + 1) % STACK_SIZE
                        st[sp] = after
                        nxt = target
                    elif is_jump:
                        nxt = target
                    else:                                    # adr, adrb, adrf
                        ax, cx = target, s
                        nxt = after
            elif op == 12:
                sp = (sp + 1) % STACK_SIZE
                st[sp] = ax + flaw
            elif op == 13:
                sp = (sp + 1) % STACK_SIZE
                st[sp] = bx + flaw
            elif op == 14:
                sp = (sp + 1) % STACK_SIZE
                st[sp] = cx + flaw
            elif op == 15:
                sp = (sp + 1) % STACK_SIZE
                st[sp] = dx + flaw
            elif 16 <= op <= 19 or op == 23:                 # pop_*, ret
                value = st[sp] + flaw
                sp = STACK_SIZE - 1 if sp == 0 else sp - 1
                if op == 16:
                    ax = value
                elif op == 17:
                    bx = value
                elif op == 18:
                    cx = value
                elif op == 19:
                    dx = value
                else:
                    nxt = value
            elif op == 2:                                    # or1
                cx ^= 1 + flaw
            elif op == 3:                                    # shl
                cx = _wrap32(cx << (1 + flaw))
            elif op == 4:                                    # zero
                cx = flaw
            elif op == 6:                                    # sub_ab: cx = ax - bx
                cx = ax + flaw - bx
            elif op == 7:                                    # sub_ac: ax = ax - cx
                ax = ax + flaw - cx
            elif op == 11:                                   # inc_c
                cx = cx + flaw + 1
            elif op == 24:                                   # mov_cd
                dx = cx + flaw
            elif op == 25:                                   # mov_ab
                bx = ax + flaw
            elif op == 30:                                   # mal
                start = self.mal(c, cx, flaw)
                if start < 0:
                    err = True
                else:
                    ax = start
            else:                                            # divide
                if not self.divide(c, c.inst + done + 1):
                    err = True
            ip = nxt % size
            done += 1
            if err:
                c.errors += 1
                self.up_reaper(c)
            if self.cr_left:
                self.cr_left -= 1
                if not self.cr_left:
                    self.cosmic_ray()
                    self.cr_left = self._draw(cr_mean)
        c.ax, c.bx, c.cx, c.dx, c.sp, c.ip = ax, bx, cx, dx, sp, ip
        c.inst += done
        self.executed += done
        self.stats["errors"] += c.errors - errors_before

    def slice_for(self, c: Cell) -> int:
        if self.constant_slice:
            return self.constant_slice
        return max(1, int(c.ms ** self.slice_power))

    def run(self, instructions: int, samples: int = 0) -> list[dict]:
        history = []
        every = max(1, instructions // samples) if samples else 0
        mark = 0
        while self.executed < instructions and self.slicer:
            if every and self.executed >= mark:
                history.append(self.sample())
                mark += every
            c = self.slicer[self.cur]
            self.run_slice(c, min(self.slice_for(c), instructions - self.executed))
            self.cur = (self.cur + 1) % len(self.slicer)
            while self.allocated > self.threshold and self.reap_for(None):
                pass
        if every:
            history.append(self.sample())
        return history

    def sample(self) -> dict:
        sizes = Counter(c.ms for c in self.cells.values())
        return {"instructions": self.executed, "creatures": len(self.cells),
                "memory_fill": round(self.allocated / self.size, 4),
                "size_classes": {str(k): sizes[k] for k in sorted(sizes)}}


# ---------------------------------------------------------------------------- the chemistry

def build(p, rng, record: bool = True) -> tuple[Soup, list[tuple[str, int]]]:
    if not p.inoculum:
        raise ValueError("inoculum must name at least one genome")
    genomes = []
    for key, count in p.inoculum.items():
        if not isinstance(count, int) or isinstance(count, bool) or count < 1:
            raise ValueError(f"inoculum counts must be positive integers, got {key!r}: {count!r}")
        genome = GENOMES.get(key) or assemble(key)
        if len(genome) < MIN_CELL_SIZE:
            raise ValueError(f"inoculum genome {key!r} has {len(genome)} instructions; the minimum cell size is {MIN_CELL_SIZE}")
        genomes += [(key if key in GENOMES else None, genome)] * count
    need = sum(len(g) for _, g in genomes) + p.gap * len(genomes)
    if need > p.soup_size:
        raise ValueError(f"soup_size {p.soup_size} cannot hold the inoculum and its gaps ({need} instructions)")
    soup = Soup(rng, p.soup_size, search_limit=p.search_limit, reaper_threshold=p.reaper_threshold,
                slicer=p.slicer, slice_power=p.slice_power, slice_size=p.slice_size,
                cosmic_ray_interval=p.cosmic_ray_interval, copy_error_interval=p.copy_error_interval,
                flaw_interval=p.flaw_interval, semantics=p.semantics, record=record)
    start, placed = 0, []
    for preferred, genome in genomes:
        cell = soup.inoculate(genome, start, preferred)
        placed.append((cell.genotype, cell.mp))
        start += len(genome) + p.gap
    return soup, placed


def generate(p, rng):
    soup, placed = build(p, rng)
    initial = Counter(g for g, _ in placed)
    history = soup.run(p.instructions, samples=50)
    final = Counter(c.genotype for c in soup.cells.values())

    names = list(dict.fromkeys([*initial, *(s for lhs, rhs, _ in soup.events.values() for s in (*lhs, *rhs)), *final]))
    species = [Species(n, structure=disassemble(soup.genomes[n])) for n in names]
    reactions = [Reaction(lhs, rhs, count=n) for lhs, rhs, n in soup.events.values()]
    order = lambda g: (int(g[:4]), g[4:])
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={g: float(n) for g, n in initial.items()},
        extras={
            "tierra": {
                "instruction_set": list(MNEMONICS),
                "genomes": "species structure: the genome as instruction-set-0 mnemonics",
                "events": "G -> G + D replication; G + H -> G + H + D replication by executing H's copy code; G -> ∅ death by the reaper; G -> G' mutation of a living genome",
                "semantics": p.semantics,
                "inoculum_addresses": [[g, a] for g, a in placed],
            },
            "analysis": {
                "instructions_executed": soup.executed,
                "births": soup.stats["births"],
                "parasitic_births": soup.stats["parasitic_births"],
                "deaths": soup.stats["deaths"],
                "cosmic_rays": soup.stats["cosmic_rays"],
                "copy_errors": soup.stats["copy_errors"],
                "flaws": soup.stats["flaws"],
                "errors": soup.stats["errors"],
                "genotypes_seen": len(soup.genomes),
                "final_population": {g: final[g] for g in sorted(final, key=order)},
                "replication_instructions": {g: m for g, m in soup.metabolism.items() if m},
                "history": history,
            },
        },
    )
