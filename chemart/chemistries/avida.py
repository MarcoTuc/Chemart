"""Avida (Adami & Brown 1994; Ofria & Wilke 2004), book 10.7.1. Catalog id: avida.

Digital organisms on a 2D toroidal lattice. Each cell holds one organism: a
circular genome of instructions from the 26-instruction ``heads_default`` set
run on a virtual CPU with three registers (AX, BX, CX), two stacks, four heads
(instruction, read, write, flow) and input/output buffers. An organism
allocates memory (h-alloc), copies itself instruction by instruction (h-copy,
with copy mutations), and divides (h-divide, with insertion and deletion
mutations). The offspring goes to an empty neighbour if there is one,
otherwise to a random neighbour or the parent's own cell, killing the
occupant. CPU time is handed out by a probabilistic scheduler, proportional
to merit; merit is the organism's size (the smallest of genome length,
copied and executed lines) times the bonus it earned for logic tasks (NOT,
NAND, ... EQU) performed on its inputs during its last gestation. Organisms
die of old age after executing age_limit times their genome length.

The machine is a port of the Avida 2.14 source (devosoft/avida:
cHardwareCPU.cc, cHardwareBase.cc, cHeadCPU.h, cCPUStack.h, cTaskLib.cc,
cEnvironment.cc, cPhenotype.cc, cPopulation.cc, Apto Scheduler::Probabilistic)
with the default avida.cfg and environment.cfg. The random number generator
is numpy's, so runs are not bit-identical to upstream.

The run is returned as the observed network of events between genotypes
(species id ``<length>-<8 hex of sha1>``, structure = the genome as letters
a..z, the instruction set's own letter code):

- ``P -> P + O``      birth into an empty cell (O may be P itself: ``P -> 2 P``)
- ``P + V -> P + O``  birth that overwrites the organism V in a neighbouring cell
- ``P -> O``          birth into the parent's own cell (the parent dies)
- ``P -> ∅``          death of old age

Task completions are recorded in extras (task_events), not as species.
"""

from __future__ import annotations

import hashlib
from collections import Counter

from chemart.network import Network, Reaction, Species

# ----------------------------------------------------------------------------
# instruction set heads_default (instset-heads.cfg), letters a..z

INSTRUCTIONS = (
    "nop-A", "nop-B", "nop-C", "if-n-equ", "if-less", "if-label", "mov-head",
    "jmp-head", "get-head", "set-flow", "shift-r", "shift-l", "inc", "dec",
    "push", "pop", "swap-stk", "swap", "add", "sub", "nand", "h-copy",
    "h-alloc", "h-divide", "IO", "h-search",
)
LETTERS = "abcdefghijklmnopqrstuvwxyz"
NUM_INST = len(INSTRUCTIONS)
(NOP_A, NOP_B, NOP_C, IF_N_EQU, IF_LESS, IF_LABEL, MOV_HEAD, JMP_HEAD, GET_HEAD,
 SET_FLOW, SHIFT_R, SHIFT_L, INC, DEC, PUSH, POP, SWAP_STK, SWAP, ADD, SUB, NAND,
 H_COPY, H_ALLOC, H_DIVIDE, IO, H_SEARCH) = range(NUM_INST)
NUM_NOPS = 3
INST_ERROR = 255                      # cInstSet::GetInstError, returned past the end

#: support/config/default-heads.org: 100 instructions, 14 of them not nop-C.
DEFAULT_ANCESTOR = "wzcag" + "c" * 86 + "zvfcaxgab"

AX, BX, CX = 0, 1, 2
HEAD_IP, HEAD_READ, HEAD_WRITE, HEAD_FLOW = 0, 1, 2, 3
STACK_SIZE = 10                       # nHardware::STACK_SIZE
MAX_LABEL = 10                        # cCodeLabel::MAX_LENGTH
MIN_GENOME_LENGTH = 8                 # avida/core/Definitions.h
MAX_GENOME_LENGTH = 2048
INPUT_SIZE = 3                        # INPUT_SIZE_DEFAULT
OFFSPRING_SIZE_RANGE = 2.0            # avida.cfg defaults
MIN_COPIED_LINES = 0.5
MIN_EXE_LINES = 0.5
MAX_LABEL_EXE_SIZE = 1
TEST_CPU_TIME_MOD = 20

#: cTaskLib.cc: the logic ids (8-bit truth tables over inputs A, B, C) of each task.
TASKS = {
    "NOT": (15, 51, 85),
    "NAND": (63, 95, 119),
    "AND": (136, 160, 192),
    "ORN": (175, 187, 207, 221, 243, 245),
    "OR": (238, 250, 252),
    "ANDN": (10, 12, 34, 48, 68, 80),
    "NOR": (3, 5, 17),
    "XOR": (60, 90, 102),
    "EQU": (153, 165, 195),
}
#: support/config/environment.cfg: process type=pow, value v -> bonus x 2^v, max_count=1.
DEFAULT_REWARDS = {"NOT": 1.0, "NAND": 1.0, "AND": 2.0, "ORN": 2.0, "OR": 3.0,
                   "ANDN": 3.0, "NOR": 4.0, "XOR": 4.0, "EQU": 5.0}
#: cEnvironment::SetupInputs(random=false), used by test CPUs.
TEST_INPUTS = (0x0F13149F, 0x3308E53E, 0x556241EB)


def to_ops(genome: str) -> list[int]:
    bad = sorted(set(genome) - set(LETTERS))
    if bad or not genome:
        raise ValueError(f"a genome is a non-empty string of instruction letters a..z, got bad letters {bad}")
    return [ord(c) - 97 for c in genome]


def to_letters(ops) -> str:
    return "".join(LETTERS[o] for o in ops)


def genotype_id(genome: str) -> str:
    return f"{len(genome)}-{hashlib.sha1(genome.encode()).hexdigest()[:8]}"


def _i32(x: int) -> int:
    return ((x + 0x80000000) & 0xFFFFFFFF) - 0x80000000


def logic_id(inputs: list[int], output: int) -> int:
    """cTaskLib::SetupTests: truth table of the output against up to three inputs (most recent first), -1 if inconsistent."""
    n = len(inputs)
    test = [(inputs[i] if i < n else 0) & 0xFFFFFFFF for i in range(3)]
    out = output & 0xFFFFFFFF
    table = [-1] * 8
    for bit in range(32):
        pos = ((test[0] >> bit) & 1) | (((test[1] >> bit) & 1) << 1) | (((test[2] >> bit) & 1) << 2)
        o = (out >> bit) & 1
        if table[pos] != -1 and table[pos] != o:
            return -1
        table[pos] = o
    if n < 1:
        table[1] = table[0]
    if n < 2:
        table[2], table[3] = table[0], table[1]
    if n < 3:
        table[4:8] = table[0:4]
    return sum(v << i for i, v in enumerate(table))


class _Uniform:
    """Buffered uniform [0, 1) draws from a numpy Generator."""

    def __init__(self, rng, block: int = 4096):
        self.rng, self.block = rng, block
        self.buf, self.i = rng.random(block).tolist(), 0

    def __call__(self) -> float:
        if self.i == self.block:
            self.buf, self.i = self.rng.random(self.block).tolist(), 0
        u = self.buf[self.i]
        self.i += 1
        return u


def _adjust(pos: int, size: int) -> int:
    """cHeadCPU::fullAdjust: wrap forward once (or modulo), clamp negatives to 0."""
    if 0 <= pos < size:
        return pos
    if size == 0 or pos < 0:
        return 0
    return pos - size if pos < 2 * size else pos % size


# ----------------------------------------------------------------------------
# the organism: virtual hardware (cHardwareCPU) and phenotype (cPhenotype)

class Organism:
    __slots__ = (
        "genome", "gid", "mem", "executed", "copied", "reg", "heads", "stacks",
        "sp", "cur_stack", "label", "read_label", "mal_active", "advance_ip",
        "merit", "cur_bonus", "copied_size", "executed_size", "child_copied_size",
        "time_used", "gestation_start", "gestation_time", "reactions", "last_tasks",
        "input_buf", "output", "input_pointer", "max_executed", "generation",
        "cell", "dead",
    )

    def __init__(self, genome: str, age_limit: int):
        self.genome = genome
        self.gid = genotype_id(genome)
        self.mem = to_ops(genome)
        n = len(self.mem)
        self.executed = bytearray(n)
        self.copied = bytearray(n)
        self.reset_hardware()
        # phenotype as after cPhenotype::SetupInject (overwritten by SetupOffspring)
        self.merit = float(n)
        self.cur_bonus = 1.0
        self.copied_size = n
        self.executed_size = n
        self.child_copied_size = n
        self.time_used = 0
        self.gestation_start = 0
        self.gestation_time = 0
        self.reactions: dict[str, int] = {}
        self.last_tasks: dict[str, int] = {}
        self.input_buf: list[int] = []       # most recent first, at most INPUT_SIZE
        self.output: int | None = None
        self.input_pointer = 0
        self.max_executed = age_limit * n if age_limit > 0 else -1   # DEATH_METHOD 2
        self.generation = 0
        self.cell = -1
        self.dead = False

    def reset_hardware(self) -> None:
        """cHardwareCPU::internalReset."""
        self.reg = [0, 0, 0]
        self.heads = [0, 0, 0, 0]
        self.stacks = [[0] * STACK_SIZE, [0] * STACK_SIZE]
        self.sp = [0, 0]
        self.cur_stack = 0
        self.label: list[int] = []
        self.read_label: list[int] = []
        self.mal_active = False
        self.advance_ip = True

    # --- stacks (cCPUStack) --------------------------------------------------
    def push(self, value: int) -> None:
        s = self.cur_stack
        self.sp[s] = STACK_SIZE - 1 if self.sp[s] == 0 else self.sp[s] - 1
        self.stacks[s][self.sp[s]] = value

    def pop(self) -> int:
        s, p = self.cur_stack, self.sp[self.cur_stack]
        value = self.stacks[s][p]
        self.stacks[s][p] = 0
        self.sp[s] = 0 if p + 1 == STACK_SIZE else p + 1
        return value

    # --- nop modifiers ---------------------------------------------------------
    def next_inst(self) -> int:
        ip, mem = self.heads[HEAD_IP], self.mem
        return INST_ERROR if ip + 1 == len(mem) else mem[ip + 1]

    def modifier(self, default: int) -> int:
        """FindModifiedRegister / FindModifiedHead: a following nop selects A, B or C."""
        if self.next_inst() < NUM_NOPS:
            ip = _adjust(self.heads[HEAD_IP] + 1, len(self.mem))
            self.heads[HEAD_IP] = ip
            self.executed[ip] = 1
            return self.mem[ip]
        return default

    def read_label_after_ip(self) -> None:
        """cHardwareCPU::ReadLabel: the nops after the IP; the IP ends on the last one."""
        label = []
        while self.next_inst() < NUM_NOPS and len(label) < MAX_LABEL:
            ip = _adjust(self.heads[HEAD_IP] + 1, len(self.mem))
            self.heads[HEAD_IP] = ip
            label.append(self.mem[ip])
            if len(label) <= MAX_LABEL_EXE_SIZE:
                self.executed[ip] = 1
        self.label = label

    def find_label_forward(self, label: list[int]) -> int:
        """cHardwareCPU::FindLabel_Forward from position 0: the first line after the match, or -1."""
        mem, size, n = self.mem, len(self.mem), len(label)
        pos = n
        while pos < size:
            if mem[pos] < NUM_NOPS:
                start, end = pos, pos + 1
                while start > 0 and mem[start - 1] < NUM_NOPS:
                    start -= 1
                while end < size and mem[end] < NUM_NOPS:
                    end += 1
                for offset in range(start, start + (end - start) - n + 1):
                    if mem[offset:offset + n] == label:
                        return n + offset
                pos = end
            pos += n
        return -1


# ----------------------------------------------------------------------------
# the world: population, scheduler, environment

class World:
    def __init__(self, rng, *, width=10, height=10, rewards=None, copy_mut_prob=0.0075,
                 divide_ins_prob=0.05, divide_del_prob=0.05, ave_time_slice=30,
                 age_limit=20, birth_method="neighborhood", prefer_empty=True,
                 record=True):
        self.rng = rng
        self.u = _Uniform(rng)
        self.width, self.height = width, height
        self.size = width * height
        self.rewards = dict(DEFAULT_REWARDS if rewards is None else rewards)
        self.task_of_id = {i: t for t in self.rewards for i in TASKS[t]}
        self.copy_mut_prob = copy_mut_prob
        self.divide_ins_prob = divide_ins_prob
        self.divide_del_prob = divide_del_prob
        self.ave_time_slice = ave_time_slice
        self.age_limit = age_limit
        self.birth_method = birth_method
        self.prefer_empty = prefer_empty
        self.cells: list[Organism | None] = [None] * self.size
        self.inputs: list[list[int]] = [list(TEST_INPUTS) for _ in range(self.size)]
        self.neighbors = [self._torus(c) for c in range(self.size)]
        self.weight = [0.0] * self.size        # Apto WeightedIndex (binary heap layout)
        self.subtree = [0.0] * self.size
        self.num_orgs = 0
        self.update = 0
        self.test_cpu = False
        self.record = record
        self.events: dict[tuple, list] = {}
        self.task_events: dict[str, Counter] = {}
        self.genotypes: dict[str, dict] = {}
        self.stats = Counter()

    def _torus(self, cell: int) -> list[int]:
        """cTopology build_torus: 8 neighbours in the order NW, N, NE, E, SE, S, SW, W."""
        x, y = cell % self.width, cell // self.width
        out = []
        for dx, dy in ((-1, -1), (0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0)):
            out.append((y + dy) % self.height * self.width + (x + dx) % self.width)
        return out

    # --- scheduler (Apto::Scheduler::Probabilistic) --------------------------
    def set_merit(self, cell: int, merit: float) -> None:
        w, sub, n = self.weight, self.subtree, self.size
        w[cell] = merit
        i = cell
        while True:
            left, right = 2 * i + 1, 2 * i + 2
            sub[i] = w[i] + (sub[left] if left < n else 0.0) + (sub[right] if right < n else 0.0)
            if i == 0:
                break
            i = (i - 1) // 2

    def next_cell(self) -> int:
        total = self.subtree[0]
        if total <= 0.0:
            return -1
        pos, i = self.u() * total, 0
        w, sub = self.weight, self.subtree
        while True:
            if pos < w[i]:
                return i
            pos -= w[i]
            left = 2 * i + 1
            if left >= self.size:           # floating-point round-off at the last leaf
                return i
            if pos < sub[left] or left + 1 >= self.size:
                i = left
            else:
                pos -= sub[left]
                i = left + 1

    # --- population -----------------------------------------------------------
    def setup_inputs(self, cell: int) -> None:
        """cEnvironment::SetupInputs(random=true): top bytes 0x0F, 0x33, 0x55 plus 24 random bits."""
        u = self.u
        self.inputs[cell] = [(15 << 24) + int(u() * (1 << 24)), (51 << 24) + int(u() * (1 << 24)),
                             (85 << 24) + int(u() * (1 << 24))]

    def kill(self, cell: int) -> None:
        org = self.cells[cell]
        if org is not None:
            org.dead = True
            self.cells[cell] = None
            self.set_merit(cell, 0.0)
            self.num_orgs -= 1

    def activate(self, org: Organism, cell: int) -> None:
        """cPopulation::ActivateOrganism: kill the occupant, new inputs, schedule."""
        self.kill(cell)
        self.setup_inputs(cell)
        org.cell = cell
        self.cells[cell] = org
        self.num_orgs += 1
        self.set_merit(cell, org.merit)

    def inject(self, genome: str, cell: int = 0) -> Organism:
        org = Organism(genome, self.age_limit)
        self.activate(org, cell)
        self._genotype(org.gid, genome, None)
        return org

    def position_offspring(self, parent_cell: int) -> int:
        """cPopulation::PositionOffspring with ALLOW_PARENT 1 (BIRTH_METHOD 0 or 4)."""
        u, cells = self.u, self.cells
        if self.birth_method == "mass-action":
            if self.prefer_empty:
                empty = [c for c in range(self.size) if cells[c] is None]
                if empty:
                    return empty[int(u() * len(empty))]
            return int(u() * self.size)
        conn = self.neighbors[parent_cell]
        found = [c for c in conn if cells[c] is None] if self.prefer_empty else []
        if not found:
            found = conn + [parent_cell]
        return found[int(u() * len(found))]

    def _genotype(self, gid: str, genome: str, parent: str | None) -> None:
        if self.record and gid not in self.genotypes:
            self.genotypes[gid] = {"genome": genome, "parent": parent, "first_update": self.update}

    def _event(self, reactants: list[str], products: list[str]) -> None:
        if not self.record:
            return
        lhs, rhs = Counter(reactants), Counter(products)
        key = (frozenset(lhs.items()), frozenset(rhs.items()))
        entry = self.events.get(key)
        if entry is None:
            self.events[key] = entry = [dict(lhs), dict(rhs), 0]
        entry[2] += 1

    # --- environment ----------------------------------------------------------
    def do_output(self, org: Organism, value: int) -> None:
        """cPhenotype::TestOutput with the pow reactions of environment.cfg (max_count=1)."""
        org.output = value
        task = self.task_of_id.get(logic_id(org.input_buf, value))
        if task is None or org.reactions.get(task, 0) >= 1:
            return
        org.reactions[task] = org.reactions.get(task, 0) + 1
        org.cur_bonus *= 2.0 ** self.rewards[task]
        if self.record and not self.test_cpu:
            self.task_events.setdefault(org.gid, Counter())[task] += 1
            self.stats["task_" + task] += 1

    # --- divide ---------------------------------------------------------------
    def divide(self, org: Organism) -> bool:
        """Inst_HeadDivide -> Divide_Main -> Divide_CheckViable, Divide_DoMutations, ActivateOffspring."""
        mem = org.mem
        size = len(mem)
        for h in range(4):
            org.heads[h] = _adjust(org.heads[h], size)
        div_point = org.heads[HEAD_READ]
        child_end = org.heads[HEAD_WRITE] or size
        child_size = child_end - div_point
        genome_size = len(org.genome)
        min_size = max(MIN_GENOME_LENGTH, int(genome_size / OFFSPRING_SIZE_RANGE))
        max_size = min(MAX_GENOME_LENGTH, int(genome_size * OFFSPRING_SIZE_RANGE))
        if not (min_size <= child_size <= max_size and min_size <= div_point <= max_size):
            return False
        executed_size = sum(org.executed[:div_point])
        if executed_size < int(div_point * MIN_EXE_LINES):
            return False
        copied_size = sum(org.copied[div_point:div_point + child_size])
        if copied_size < int(child_size * MIN_COPIED_LINES):
            return False
        org.executed_size = executed_size
        org.child_copied_size = copied_size

        child = mem[div_point:div_point + child_size]
        del mem[div_point:]
        del org.executed[div_point:]
        del org.copied[div_point:]
        u = self.u
        if self.divide_ins_prob > 0 and u() < self.divide_ins_prob and len(child) < MAX_GENOME_LENGTH:
            child.insert(int(u() * (len(child) + 1)), int(u() * NUM_INST))
        if self.divide_del_prob > 0 and u() < self.divide_del_prob and len(child) > MIN_GENOME_LENGTH:
            del child[int(u() * len(child))]
        org.mal_active = False
        org.advance_ip = False
        self.stats["divides"] += 1
        parent_alive = self.activate_offspring(org, to_letters(child))
        if parent_alive:
            org.reset_hardware()
            org.advance_ip = False
            org.executed = bytearray(len(mem))
            org.copied = bytearray(len(mem))
        return True

    def gestation_reset(self, org: Organism) -> None:
        """cPhenotype::DivideReset: lock in merit = size merit x bonus, reset the bonus and counts."""
        size_merit = min(len(org.genome), org.copied_size, org.executed_size)   # BASE_MERIT_METHOD 4
        org.merit = float(size_merit * org.cur_bonus)
        org.gestation_time = org.time_used - org.gestation_start
        org.gestation_start = org.time_used
        org.last_tasks = dict(org.reactions)
        org.cur_bonus = 1.0
        org.reactions = {}

    def activate_offspring(self, parent: Organism, genome: str) -> bool:
        self.gestation_reset(parent)
        if self.test_cpu:
            self.offspring = genome
            return True
        child = Organism(genome, self.age_limit)
        child.merit = parent.merit                          # INHERIT_MERIT 1
        child.copied_size = parent.child_copied_size
        child.executed_size = parent.executed_size
        child.gestation_time = parent.gestation_time
        parent.generation += 1                              # GENERATION_INC_METHOD 1
        child.generation = parent.generation
        pcell = parent.cell
        target = self.position_offspring(pcell)
        victim = self.cells[target]
        self._genotype(child.gid, genome, parent.gid)
        self.stats["births"] += 1
        if target == pcell:
            self._event([parent.gid], [child.gid])
            self.stats["parent_replaced"] += 1
            parent_alive = False
        else:
            if victim is None:
                self._event([parent.gid], [parent.gid, child.gid])
            else:
                self._event([parent.gid, victim.gid], [parent.gid, child.gid])
                self.stats["overwrites"] += 1
            self.set_merit(pcell, parent.merit)
            parent_alive = True
        self.activate(child, target)
        return parent_alive

    # --- the CPU cycle (cHardwareCPU::SingleProcess) ---------------------------
    def process(self, org: Organism) -> None:
        mem, heads, reg = org.mem, org.heads, org.reg
        org.time_used += 1
        org.advance_ip = True
        ip = _adjust(heads[HEAD_IP], len(mem))
        heads[HEAD_IP] = ip
        op = mem[ip]
        org.executed[ip] = 1

        if op < NUM_NOPS:
            pass
        elif op == H_COPY:
            size = len(mem)
            r = heads[HEAD_READ] = _adjust(heads[HEAD_READ], size)
            w = heads[HEAD_WRITE] = _adjust(heads[HEAD_WRITE], size)
            inst = mem[r]
            if inst < NUM_NOPS:
                if len(org.read_label) < MAX_LABEL:
                    org.read_label.append(inst)
            else:
                org.read_label = []
            if self.copy_mut_prob > 0 and self.u() < self.copy_mut_prob:
                inst = int(self.u() * NUM_INST)
                self.stats["copy_mutations"] += 1
            mem[w] = inst
            org.copied[w] = 1
            heads[HEAD_READ] = _adjust(r + 1, size)
            heads[HEAD_WRITE] = _adjust(w + 1, size)
        elif op == IF_LABEL:
            org.read_label_after_ip()
            label = [(n + 1) % NUM_NOPS for n in org.label]
            if label != org.read_label:
                heads[HEAD_IP] = _adjust(heads[HEAD_IP] + 1, len(mem))
        elif op == MOV_HEAD:
            h = org.modifier(HEAD_IP)
            heads[h] = heads[HEAD_FLOW]
            if h == HEAD_IP:
                org.advance_ip = False
        elif op == H_SEARCH:
            org.read_label_after_ip()
            label = [(n + 1) % NUM_NOPS for n in org.label]
            found = heads[HEAD_IP]
            if label:
                pos = org.find_label_forward(label)
                if pos >= 0:
                    found = _adjust(pos - 1, len(mem))
            reg[BX] = _i32(found - heads[HEAD_IP])
            reg[CX] = len(label)
            heads[HEAD_FLOW] = _adjust(found + 1, len(mem))
        elif op == H_ALLOC:
            cur = len(mem)
            alloc = min(int(OFFSPRING_SIZE_RANGE * cur), MAX_GENOME_LENGTH - cur)
            new = cur + alloc
            if (not org.mal_active and alloc >= 1 and MIN_GENOME_LENGTH <= new <= MAX_GENOME_LENGTH
                    and alloc <= int(cur * OFFSPRING_SIZE_RANGE) and cur <= int(alloc * OFFSPRING_SIZE_RANGE)):
                mem.extend([NOP_A] * alloc)                  # ALLOC_METHOD 0: default instruction
                org.executed.extend(bytes(alloc))
                org.copied.extend(bytes(alloc))
                org.mal_active = True
                reg[AX] = cur
        elif op == H_DIVIDE:
            if self.divide(org) and (org.dead or self.test_cpu):
                return
        elif op == IO:
            r = org.modifier(BX)
            self.do_output(org, reg[r])
            cell_inputs = self.inputs[org.cell] if org.cell >= 0 else list(TEST_INPUTS)
            org.input_pointer %= len(cell_inputs)
            value = cell_inputs[org.input_pointer]
            org.input_pointer += 1
            reg[r] = value
            org.input_buf.insert(0, value)
            del org.input_buf[INPUT_SIZE:]
        elif op == IF_N_EQU:
            a = org.modifier(BX)
            if reg[a] == reg[(a + 1) % 3]:
                heads[HEAD_IP] = _adjust(heads[HEAD_IP] + 1, len(mem))
        elif op == IF_LESS:
            a = org.modifier(BX)
            if reg[a] >= reg[(a + 1) % 3]:
                heads[HEAD_IP] = _adjust(heads[HEAD_IP] + 1, len(mem))
        elif op == JMP_HEAD:
            h = org.modifier(HEAD_IP)
            heads[h] = _adjust(heads[h] + reg[CX], len(mem))
        elif op == GET_HEAD:
            h = org.modifier(HEAD_IP)
            reg[CX] = heads[h]
        elif op == SET_FLOW:
            r = org.modifier(CX)
            heads[HEAD_FLOW] = _adjust(reg[r], len(mem))
        elif op == SHIFT_R:
            r = org.modifier(BX)
            reg[r] >>= 1
        elif op == SHIFT_L:
            r = org.modifier(BX)
            reg[r] = _i32(reg[r] << 1)
        elif op == INC:
            r = org.modifier(BX)
            reg[r] = _i32(reg[r] + 1)
        elif op == DEC:
            r = org.modifier(BX)
            reg[r] = _i32(reg[r] - 1)
        elif op == PUSH:
            r = org.modifier(BX)
            org.push(reg[r])
        elif op == POP:
            r = org.modifier(BX)
            reg[r] = org.pop()
        elif op == SWAP_STK:
            org.cur_stack = 1 - org.cur_stack
        elif op == SWAP:
            a = org.modifier(BX)
            b = (a + 1) % 3
            reg[a], reg[b] = reg[b], reg[a]
        elif op == ADD:
            r = org.modifier(BX)
            reg[r] = _i32(reg[BX] + reg[CX])
        elif op == SUB:
            r = org.modifier(BX)
            reg[r] = _i32(reg[BX] - reg[CX])
        elif op == NAND:
            r = org.modifier(BX)
            reg[r] = ~(reg[BX] & reg[CX])

        if org.advance_ip:
            heads[HEAD_IP] = _adjust(heads[HEAD_IP] + 1, len(org.mem))
        if org.max_executed > 0 and org.time_used >= org.max_executed and not self.test_cpu:
            self._event([org.gid], [])
            self.stats["age_deaths"] += 1
            self.kill(org.cell)

    # --- updates --------------------------------------------------------------
    def run(self, updates: int, history: dict | None = None) -> None:
        cells = self.cells
        for _ in range(updates):
            self.update += 1
            steps = self.ave_time_slice * self.num_orgs        # cWorld::CalculateUpdateSize
            for _ in range(steps):
                if self.num_orgs == 0:
                    break
                self.process(cells[self.next_cell()])
            self.stats["instructions"] += steps
            if history is not None:
                self._sample(history)

    def _sample(self, h: dict) -> None:
        orgs = [o for o in self.cells if o is not None]
        h["organisms"].append(len(orgs))
        h["genotypes"].append(len({o.gid for o in orgs}))
        for key in ("births", "age_deaths", "overwrites"):      # events during this update
            h[key].append(self.stats[key] - self.stats["sampled_" + key])
            self.stats["sampled_" + key] = self.stats[key]
        h["ave_merit"].append(sum(o.merit for o in orgs) / len(orgs) if orgs else 0.0)
        gest = [o.gestation_time for o in orgs if o.gestation_time > 0]
        h["ave_gestation"].append(sum(gest) / len(gest) if gest else 0.0)
        h["ave_generation"].append(sum(o.generation for o in orgs) / len(orgs) if orgs else 0.0)
        for task in TASKS:
            h["tasks"][task].append(sum(1 for o in orgs if o.last_tasks.get(task)))


def test_cpu(genome: str, rewards: dict | None = None, inputs=TEST_INPUTS, rng=None) -> dict:
    """Run one organism alone until its first divide (cTestCPU): gestation, merit, sizes, tasks, offspring.

    Mutations are off. Gives up after TEST_CPU_TIME_MOD x genome length cycles.
    """
    import numpy as np

    world = World(rng if rng is not None else np.random.default_rng(0), width=1, height=1,
                  rewards=rewards, copy_mut_prob=0.0, divide_ins_prob=0.0, divide_del_prob=0.0,
                  age_limit=0, record=False)
    world.test_cpu = True
    org = Organism(genome, 0)
    world.inputs = [list(inputs)]
    org.cell = 0
    world.offspring = None
    for _ in range(TEST_CPU_TIME_MOD * len(genome)):
        world.process(org)
        if world.offspring is not None:
            break
    return {
        "divided": world.offspring is not None,
        "gestation": org.gestation_time if world.offspring is not None else None,
        "merit": org.merit,
        "executed_size": org.executed_size,
        "copied_size": org.child_copied_size,
        "tasks": sorted(org.last_tasks),
        "offspring": world.offspring,
    }


# ----------------------------------------------------------------------------

def check_rewards(rewards) -> dict:
    if not isinstance(rewards, dict):
        raise ValueError(f"rewards must be a dict {{task: exponent}}, got {rewards!r}")
    bad = sorted(set(rewards) - set(TASKS))
    if bad:
        raise ValueError(f"rewards: unknown tasks {bad}; tasks are {list(TASKS)}")
    for task, v in rewards.items():
        if not isinstance(v, (int, float)) or isinstance(v, bool):
            raise ValueError(f"rewards[{task!r}] must be a number (bonus multiplier 2^value), got {v!r}")
    return {t: float(rewards[t]) for t in TASKS if t in rewards}


def generate(p, rng):
    ops = to_ops(p.ancestor)
    if not MIN_GENOME_LENGTH <= len(ops) <= MAX_GENOME_LENGTH:
        raise ValueError(f"ancestor length must be in [{MIN_GENOME_LENGTH}, {MAX_GENOME_LENGTH}], got {len(ops)}")
    rewards = check_rewards(p.rewards)
    world = World(rng, width=p.world_x, height=p.world_y, rewards=rewards,
                  copy_mut_prob=p.copy_mut_prob, divide_ins_prob=p.divide_ins_prob,
                  divide_del_prob=p.divide_del_prob, ave_time_slice=p.ave_time_slice,
                  age_limit=p.age_limit, birth_method=p.birth_method, prefer_empty=p.prefer_empty)
    ancestor = world.inject(p.ancestor, 0)
    history = {k: [] for k in ("organisms", "genotypes", "births", "age_deaths", "overwrites",
                                "ave_merit", "ave_gestation", "ave_generation")}
    history["tasks"] = {t: [] for t in TASKS}
    world.run(p.updates, history)

    final = Counter(o.gid for o in world.cells if o is not None)
    used = {ancestor.gid} | set(final)
    reactions = []
    for lhs, rhs, n in world.events.values():
        used.update(lhs)
        used.update(rhs)
        reactions.append(Reaction(lhs, rhs, count=n))
    species = [Species(g, structure=world.genotypes[g]["genome"])
               for g in sorted(used, key=lambda g: (world.genotypes[g]["first_update"], g))]
    dominant = max(final, key=lambda g: (final[g], g)) if final else None
    s = world.stats
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={ancestor.gid: 1.0},
        extras={
            "space": {
                "lattice": "torus", "width": p.world_x, "height": p.world_y, "neighbourhood": "moore",
                "final_grid": [o.gid if o is not None else None for o in world.cells],
            },
            "final_state": dict(sorted(final.items())),
            "task_events": {g: dict(sorted(c.items())) for g, c in sorted(world.task_events.items())},
            "genotype_parent": {g: world.genotypes[g]["parent"] for g in sorted(used)},
            "environment": {t: {"logic_ids": list(TASKS[t]), "bonus_multiplier": 2.0 ** v}
                            for t, v in rewards.items()},
            "instruction_set": {LETTERS[i]: name for i, name in enumerate(INSTRUCTIONS)},
            "analysis": {
                "updates": p.updates,
                "instructions_executed": s["instructions"],
                "births": s["births"],
                "overwrites": s["overwrites"],
                "parent_replaced": s["parent_replaced"],
                "age_deaths": s["age_deaths"],
                "copy_mutations": s["copy_mutations"],
                "tasks_performed": {t: s["task_" + t] for t in TASKS},
                "dominant": dominant,
                "per_update": history,
            },
        },
    )
