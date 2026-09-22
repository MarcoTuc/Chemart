"""Stringmol: Hickinbotham, Clark, Stepney et al. (York, 2009-2011), book 11.1.2.

Catalog id: stringmol.

Molecules are strings over 33 symbols: 26 template codes A..Z and 7 function
codes ($ search, > move, ^ toggle, ? if, = copy, % cleave, } end). Two
molecules bind with a probability derived from a Smith-Waterman alignment of
the complement of one string (template codes shifted by 13, function codes
self-complementary) against the other. The molecule whose aligned segment
starts further from the beginning of its string becomes active and runs its
string as a program over four pointers (instruction, flow, read, write), each
of which can sit on either string. Copy ('=') may mutate. The reaction ends
at '}' (or when the instruction pointer runs off the string) and the pair
dissociates; molecules cleaved off ('%') on the way are the products.

The machine is a port of the authors' C++ source (franticspider/stringmol,
stringPM.cpp, instructions.cpp, alignment.cpp, stringmanip.cpp), including
its single-precision float arithmetic, and follows the technical report
YCS-2010-458 (spec v0.2) where the two agree; see the catalog decisions for
where they don't.

Faces:
- generate: chemart.expand.expand of the seed set, with exact copying (the
  closure).
- evolve, reactor "container": the authors' time-stepped container
  (stringPM::make_next with the ALife XII main loop): every time step each
  molecule, in random order, may decay, try to bind, or execute one
  instruction of its complex, and energy is added per step. Observed network:
  one reaction per complete bind-execute-dissociate event (reactants at bind,
  products at dissociation).
- evolve, reactor "soup": chemart.soup.stir with instantaneous pair reactions
  (bind test, then the program run to completion) and constant population
  size.
"""

from __future__ import annotations

import math
from array import array
from collections import Counter

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species
from chemart.soup import Tally, stir
from chemart.trajectory import Frame

#: Symbol order of the substitution matrix and of the mutation loop
#: (alignment.cpp default_table / config/ALXII.mtx).
KEY = "ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ"
FUNCTION_CODES = "$>^?=%}"
#: Seed replicase of spec v0.2 appendix B.1 (proper Smith-Waterman traceback).
SEED_REPLICASE = "OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"
#: Seed replicase of the upstream config files (e.g. config/test1.conf).
CONFIG_REPLICASE = "WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"
#: ALife XII species 9 (paper fig. 7, run 277).
ALIFE12_SPECIES_9 = "OBEQBXUUUDYGRHBBOSEOLHHHRLUEUOBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB"

#: stringPM::get_bprob and HSearch: min(score, l - m) / (l - m) with m = 1.124.
MISMATCH = 1.124

UNBOUND, PASSIVE, ACTIVE = 1, 2, 3
_A, _B, _C = ord("A"), ord("B"), ord("C")
_KEY_BYTES = KEY.encode()
_ID = str.maketrans({"$": "s", ">": "m", "^": "t", "?": "i", "=": "c", "%": "x", "}": "e"})

_F = array("f", [0.0])


def f32(x: float) -> float:
    """Round to single precision (the C source stores most reals in float)."""
    _F[0] = x
    return _F[0]


def _complement_table() -> bytes:
    """stringmanip.cpp AlphaComp: template codes +13 mod 26, everything else unchanged."""
    t = bytearray(range(256))
    for c in range(65, 91):
        t[c] = (c - 65 + 13) % 26 + 65
    return bytes(t)


COMPLEMENT = _complement_table()


def _substitution_table():
    """alignment.cpp default_table ("ALXII values"): spec eqs. 7-9 over KEY, printed to 3 decimals."""
    n = len(KEY)
    total = sum(min(k, n - k) for k in range(n))        # 272
    sub = [None] * 256
    for i, a in enumerate(KEY):
        row = [0.0] * 256
        for j, b in enumerate(KEY):
            if i == j:
                v = 1.0 if a.isalpha() else 0.5
            else:
                d = abs(i - j)
                v = -min(d, n - d) * n / total
            row[ord(b)] = f32(float(f"{v:.3f}"))
        sub[ord(a)] = row
    return sub, f32(-1.333)


SUBSTITUTION, INDEL = _substitution_table()


def complement(seq: str) -> str:
    return seq.encode().translate(COMPLEMENT).decode()


def species_id(seq: bytes | str) -> str:
    """Safe id: template codes as they are, function codes as s m t i c x e."""
    text = seq.decode("ascii") if isinstance(seq, bytes) else seq
    return text.translate(_ID)


# --- alignment -----------------------------------------------------------------
_SW_CACHE: dict = {}


def smith_waterman(s1: bytes, s2: bytes, traceback: str = "matrix", table=None):
    """Local alignment of s1 (rows) against s2 (columns): (score, start1, end1, start2, end2).

    Starts are 0-based, ends exclusive. traceback "matrix" is SmithWatermanV2
    (spec v0.2); "highest-neighbour" is the v0.1 SmithWaterman that walks back
    through the largest neighbouring cell (spec appendix B.1). Returns None if
    either string is empty (the C code then leaves the result undefined).
    table: optional (substitution rows by byte, indel score); default the source table.
    """
    key = (s1, s2, traceback)
    if table is None:
        hit = _SW_CACHE.get(key)
        if hit is not None or key in _SW_CACHE:
            return hit
    l1, l2 = len(s1), len(s2)
    if not l1 or not l2:
        return None
    sub, g = (SUBSTITUTION, INDEL) if table is None else table
    f = _F
    prev = [0.0] * (l2 + 1)
    rows = [prev]
    trace = [bytes(l2 + 1)]
    emax, endi, endj = 0.0, 0, 0
    for i in range(1, l1 + 1):
        w = sub[s1[i - 1]]
        cur = [0.0] * (l2 + 1)
        tr = bytearray(l2 + 1)
        left = 0.0
        for j in range(1, l2 + 1):
            h, t = 0.0, 0
            f[0] = prev[j - 1] + w[s2[j - 1]]
            if f[0] > h:
                h, t = f[0], 1
            f[0] = prev[j] + g
            if f[0] > h:
                h, t = f[0], 2
            f[0] = left + g
            if f[0] > h:
                h, t = f[0], 3
            cur[j] = left = h
            tr[j] = t
            if h > emax:
                emax, endi, endj = h, i, j
        rows.append(cur)
        trace.append(tr)
        prev = cur
    i, j = endi, endj
    if traceback == "matrix":
        while trace[i][j]:
            t = trace[i][j]
            if t == 1:
                i -= 1
                j -= 1
            elif t == 2:
                i -= 1
            else:
                j -= 1
    else:
        while rows[i][j] > 0:
            up, diag, left = rows[i - 1][j], rows[i - 1][j - 1], rows[i][j - 1]
            action = -1
            if up >= diag and up >= left:
                action = 0
            if left >= diag and left >= up:
                action = 2
            if diag >= up and diag >= left:
                action = 1
            if action == 0:
                i -= 1
            elif action == 1:
                i -= 1
                j -= 1
            else:
                j -= 1
    out = (emax, i, endi, j, endj)
    if table is None:
        if len(_SW_CACHE) > 100_000:
            _SW_CACHE.clear()
        _SW_CACHE[key] = out
    return out


def _soft_probability(sw) -> float:
    """min(score, l - 1.124) / (l - 1.124) with l the shorter aligned length, in C float arithmetic."""
    score, s1, e1, s2, e2 = sw
    lim = min(e1 - s1, e2 - s2) - MISMATCH
    s = f32(score if score < lim else lim)
    return f32(s / lim)


def bind_probability(a: str | bytes, b: str | bytes, traceback: str = "matrix") -> float:
    """P(bind) of molecule a (initiator) with b: stringPM::get_sw + get_bprob."""
    a = a.encode() if isinstance(a, str) else a
    b = b.encode() if isinstance(b, str) else b
    sw = smith_waterman(a.translate(COMPLEMENT), b, traceback)
    if sw is None or min(sw[2] - sw[1], sw[4] - sw[3]) <= 2:
        return 0.0
    return _soft_probability(sw)


# --- molecules and the machine -------------------------------------------------
def _strlen(buf: bytearray) -> int:
    n = buf.find(0)
    return len(buf) if n < 0 else n


def _char(buf: bytearray, k: int) -> int:
    return buf[k] if 0 <= k < len(buf) else 0


class Molecule:
    """An s_ag: a zero-filled string buffer of max_length + 1 bytes, pointers and toggles."""

    __slots__ = ("S", "len", "status", "exec", "pas", "i", "r", "w", "f",
                 "it", "rt", "wt", "ft", "rec")

    def __init__(self, seq: bytes, size: int):
        self.S = bytearray(size)
        self.S[:len(seq)] = seq
        self.len = len(seq)
        self.status = UNBOUND
        self.exec = self.pas = self.rec = None
        self.i, self.r, self.w, self.f = [0, 0], [0, 0], [0, 0], [0, 0]
        self.it = self.rt = self.wt = self.ft = 0

    def seq(self) -> bytes:
        return bytes(self.S[:_strlen(self.S)])


class Machine:
    """Binding and microprogram execution (stringPM, instructions.cpp), with a pluggable rand()."""

    def __init__(self, rand, *, max_length=2000, substitution_rate=0.0, indel_rate=0.0,
                 traceback="matrix", on_finish=None):
        self.rand = rand
        self.maxl = max_length
        self.size = max_length + 1
        self.indel = f32(indel_rate)
        self.subindel = f32(f32(substitution_rate) + f32(indel_rate))
        self.traceback = traceback
        self.on_finish = on_finish

    # binding ------------------------------------------------------------
    def get_sw(self, a: Molecule, b: Molecule):
        return smith_waterman(a.seq().translate(COMPLEMENT), b.seq(), self.traceback)

    @staticmethod
    def bprob(sw) -> float:
        if min(sw[2] - sw[1], sw[4] - sw[3]) <= 2:
            return 0.0
        return _soft_probability(sw)

    @staticmethod
    def set_exec(a: Molecule, b: Molecule, sw):
        _, s1, _, s2, _ = sw
        if s1 >= s2:
            act, pas, aidx, pidx = a, b, s1, s2
        else:
            act, pas, aidx, pidx = b, a, s2, s1
        act.status, act.exec, act.pas = ACTIVE, None, pas
        pas.status, pas.exec, pas.pas = PASSIVE, act, None
        act.i, act.r, act.w, act.f = [pidx, aidx], [pidx, aidx], [pidx, aidx], [pidx, aidx]
        act.it = act.rt = act.wt = act.ft = 1
        pas.i, pas.r, pas.w, pas.f = [0, 0], [0, 0], [0, 0], [0, 0]
        pas.it = pas.rt = pas.wt = pas.ft = 0
        act.rec = {"reactants": (act.seq(), pas.seq()), "released": []}
        return act, pas

    @staticmethod
    def unbind(m: Molecule):
        m.status, m.exec, m.pas, m.rec = UNBOUND, None, None, None
        m.i, m.r, m.w, m.f = [0, 0], [0, 0], [0, 0], [0, 0]
        m.it = m.rt = m.wt = m.ft = 0

    def finish(self, act: Molecule, survivors):
        rec = act.rec
        rec["products"] = tuple(m.seq() for m in survivors) + tuple(rec["released"])
        if self.on_finish:
            self.on_finish(rec)

    # instructions -------------------------------------------------------
    @staticmethod
    def lab_length(buf: bytearray, pos: int) -> int:
        n, k = 0, pos + 1
        while 65 <= _char(buf, k) <= 90:
            n += 1
            k += 1
        return n

    def hsearch(self, ibuf, ipos, cs):
        """'$' (HSearch with SOFT_SEARCH): (new flow position, True if it is on cs)."""
        n = self.lab_length(ibuf, ipos)
        if not n:
            return ipos, False
        tmp = bytes(ibuf[ipos + 1:ipos + 1 + n]).translate(COMPLEMENT)
        sp = bytes(cs[:_strlen(cs)])
        sw = smith_waterman(tmp, sp, self.traceback)
        if sw is None:                       # empty string: undefined in C
            sw = (0.0, 0, 0, 0, 0)
        bprob = _soft_probability(sw)
        if f32(self.rand()) < bprob:
            return sw[4], True
        return ipos + n, False

    def iflabel(self, ibuf, ipos, rbuf, rpos):
        """'?' (IfLabel): the new instruction position."""
        n = self.lab_length(ibuf, ipos)
        ip = ipos + 1
        if n <= 1:
            return ip + n + 1 if not _char(rbuf, rpos) else ip + n
        tmp = bytes(ibuf[ip:ip + n]).translate(COMPLEMENT)
        tmp2 = bytes(rbuf[rpos:rpos + n]) if 0 <= rpos < len(rbuf) else b""
        z = tmp2.find(0)
        if z >= 0:
            tmp2 = tmp2[:z]
        sw = smith_waterman(tmp, tmp2, self.traceback)
        score = sw[0] if sw is not None else 0.0
        prob = f32(f32(score / n) ** n)
        if f32(self.rand()) < prob:
            return ip + n + 1
        return ip + n

    def sym_from_adj(self, c: int) -> int:
        idx = _KEY_BYTES.index(c)
        n = len(KEY)
        neighbours = sorted({(idx - 1) % n, (idx + 1) % n})
        k = math.ceil(f32(2.0 * f32(self.rand())))
        return _KEY_BYTES[neighbours[k - 1]] if k else _KEY_BYTES[-1]

    def hcopy(self, act: Molecule) -> int:
        pas = act.pas
        act.len, pas.len = _strlen(act.S), _strlen(pas.S)
        wt, rt = act.wt, act.rt
        wbuf = act.S if wt else pas.S
        wpos = act.w[wt]
        if wpos >= self.maxl:
            wbuf[self.maxl] = 0
            act.i[act.it] += 1
            return -1
        if act.r[rt] >= self.maxl:
            act.i[act.it] += 1
            return -2
        rbuf = act.S if rt else pas.S
        c = _char(rbuf, act.r[rt])
        if c:
            rno = f32(self.rand())
            if rno < self.indel:
                if f32(self.rand()) < 0.5:           # insert: the copy plus a random symbol
                    wbuf[wpos] = c
                    act.w[wt] += 1
                    cidx = int(f32(f32(self.rand()) * len(KEY)))
                    if act.w[wt] < len(wbuf):
                        wbuf[act.w[wt]] = _KEY_BYTES[cidx]
                    act.w[wt] += 1
                else:                                  # delete: skip the symbol (and I moves on)
                    act.i[act.it] += 1
                act.r[rt] += 1
            else:
                wbuf[wpos] = self.sym_from_adj(c) if rno < self.subindel else c
                act.w[wt] += 1
                act.r[rt] += 1
        act.len, pas.len = _strlen(act.S), _strlen(pas.S)
        act.i[act.it] += 1
        return 0

    @staticmethod
    def rewind_bad_ptrs(act: Molecule):
        plen = _strlen(act.pas.S)
        arrays = (act.i, act.r, act.w, act.f)
        if plen:
            for a in arrays:
                if a[0] > plen or a[0] < 0:
                    a[0] = plen
        else:
            for a in arrays:
                a[0] = 0
            act.it = act.rt = act.wt = act.ft = 1
        alen = _strlen(act.S)
        if alen:
            for a in arrays:
                if a[1] > alen or a[1] < 0:
                    a[1] = alen
        else:
            for a in arrays:
                a[1] = 0
            if plen:
                act.it = act.rt = act.wt = act.ft = 0

    def cleave(self, act: Molecule, out: list) -> int:
        pas = act.pas
        csite = act if act.ft else pas
        fpos = act.f[act.ft]
        dac = 0
        if fpos < csite.len:
            cs = csite.S
            n = _strlen(cs)
            child = Molecule(bytes(cs[fpos:n]), self.size)
            act.rec["released"].append(child.seq())
            out.append(child)
            cs[fpos:n] = bytes(n - fpos)
            csite.len = _strlen(cs)
            self.rewind_bad_ptrs(act)
            if not _strlen(act.S):
                dac = 1
            elif not _strlen(pas.S):
                dac = 2
            if dac == 1:                               # active string used up
                self.finish(act, [pas])
                self.unbind(pas)
                out.append(pas)
            elif dac == 2:                             # passive string used up
                self.finish(act, [act])
                self.unbind(act)
                out.append(act)
        if not dac:
            act.i[act.it] += 1
        return dac

    def exec_step(self, act: Molecule, pas: Molecule, out: list):
        """One instruction of the complex (stringPM::exec_step); appends survivors to out."""
        ibuf = act.S if act.it else pas.S
        ipos = act.i[act.it]
        op = _char(ibuf, ipos)
        safe_append = True
        if op == 36:                                   # $ search
            cs = act.S if act.ft else pas.S
            pos, on_cs = self.hsearch(ibuf, ipos, cs)
            if not on_cs:
                act.ft = act.it
            act.f[act.ft] = pos
            act.i[act.it] += 1
        elif op == 62:                                 # > move
            mod = _char(ibuf, ipos + 1)
            if mod == _B:
                act.rt = act.ft
                act.r[act.rt] = act.f[act.ft]
            elif mod == _C:
                act.wt = act.ft
                act.w[act.wt] = act.f[act.ft]
            else:
                act.it = act.ft
                act.i[act.it] = act.f[act.ft]
            act.i[act.it] += 1
        elif op == 61:                                 # = copy
            if self.hcopy(act) < 0:
                self.finish(act, [act, pas])
                self.unbind(act)
                self.unbind(pas)
        elif op == 94:                                 # ^ toggle
            mod = _char(ibuf, ipos + 1)
            if mod == _A:
                act.it = 1 - act.it
            elif mod == _B:
                act.rt = 1 - act.rt
            elif mod == _C:
                act.wt = 1 - act.wt
            else:
                act.ft = 1 - act.ft
            act.i[act.it] += 1
        elif op == 63:                                 # ? if
            rbuf = act.S if act.rt else pas.S
            act.i[act.it] = self.iflabel(ibuf, ipos, rbuf, act.r[act.rt])
        elif op == 37:                                 # % cleave
            if self.cleave(act, out):
                safe_append = False
        elif op == 0 or op == 125:                     # } end (or I off the string)
            self.finish(act, [act, pas])
            self.unbind(act)
            self.unbind(pas)
        else:                                          # template codes: no-op
            act.i[act.it] += 1
        if safe_append:
            out.append(act)
            out.append(pas)

    # one isolated reaction ----------------------------------------------
    def react(self, a: bytes, b: bytes, bind: str, budget: int):
        """a initiates a bind with b; returns the products, None (no bind) or False (no end within budget)."""
        ma, mb = Molecule(a, self.size), Molecule(b, self.size)
        sw = self.get_sw(ma, mb)
        if sw is None:
            return None
        p = self.bprob(sw)
        if bind == "sample" and not f32(self.rand()) < p:
            return None
        if bind == "possible" and p <= 0:
            return None
        act, pas = self.set_exec(ma, mb, sw)
        rec = act.rec
        sink: list = []
        for _ in range(budget):
            self.exec_step(act, pas, sink)
            if "products" in rec:
                return rec["products"]
            sink.clear()
        return False


# --- the container -------------------------------------------------------------
class Container:
    """stringPM::make_next + update, with the SmPm_AlifeXII loop (energy added after each step)."""

    def __init__(self, machine: Machine, molecules, *, energy_per_step, cell_radius, agent_radius, decay):
        self.m = machine
        machine.on_finish = self._finished
        self.now = [Molecule(s, machine.size) for s in molecules]
        self.energy = 0
        self.estep = energy_per_step
        self.decay = f32(decay)
        agarea = f32(f32(math.pi) * float(f32(agent_radius)) ** 2)
        cellarea = f32(math.pi * float(f32(cell_radius)) ** 2)
        self.arearatio = f32(agarea / cellarea)
        self.seen: dict[bytes, None] = dict.fromkeys(molecules)
        self.fired: dict = {}
        self.tally = Tally()
        self.epochs: list = []           # [steps done, dominant sequence] when the most common changes
        self.steps_run = 0
        self.aborted: dict = {}
        self.decayed: Counter = Counter()

    def _key(self, lhs, rhs):
        return frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items())

    def _finished(self, rec):
        lhs, rhs = rec["reactants"], rec["products"]
        self.seen.update(dict.fromkeys(rhs))
        entry = self.fired.setdefault(self._key(lhs, rhs), [lhs, rhs, 0, Counter()])
        entry[2] += 1
        entry[3][lhs[0]] += 1
        self.tally.add(_ids(lhs), _ids(rhs))

    def eqn_prop(self, n: int) -> bool:
        if not n:
            return False
        cov = f32(1.0 - (1.0 - self.arearatio) ** n)
        return f32(self.m.rand()) < cov

    def step(self):
        m, rand = self.m, self.m.rand
        now, nxt = self.now, []
        unbound = sum(1 for a in now if a.status == UNBOUND)
        while now:
            pag = now.pop(int(len(now) * rand()))
            bag = None
            if pag.status == UNBOUND:
                unbound -= 1
            elif pag.status == ACTIVE:
                bag = pag.pas
                now.remove(bag)
            else:
                bag = pag.exec
                now.remove(bag)
            if f32(rand()) < self.decay:
                if bag is None:
                    self.decayed[pag.seq()] += 1
                else:
                    act = pag if pag.status == ACTIVE else bag
                    rec = act.rec
                    key = (rec["reactants"], tuple(rec["released"]))
                    self.aborted[key] = self.aborted.get(key, 0) + 1
                continue
            changed = False
            if self.energy > 0:
                if pag.status == UNBOUND:
                    changed = found = self.eqn_prop(unbound)
                    if found:
                        pool = [a for a in now if a.status == UNBOUND]
                        pos = unbound
                        while pos == unbound:
                            pos = int(unbound * rand())
                        bag = pool[pos]
                        sw = m.get_sw(pag, bag)
                        if sw is not None and f32(rand()) < m.bprob(sw):
                            m.set_exec(pag, bag, sw)
                            nxt.append(pag)
                            now.remove(bag)
                            unbound -= 1
                            nxt.append(bag)
                            self.energy -= 1
                        else:
                            changed = False
                        bag = None
                else:
                    act, pas = (pag, pag.pas) if pag.status == ACTIVE else (pag.exec, pag)
                    m.exec_step(act, pas, nxt)
                    self.energy -= 1
                    changed = True
            if not changed:
                nxt.append(pag)
                if bag is not None:
                    nxt.append(bag)
        self.now = nxt

    def state(self) -> dict[str, float]:
        """Every molecule by species; a bound one as the sequence it had when it bound."""
        counts = Counter(a.seq() for a in self.now if a.status == UNBOUND)
        for a in self.now:
            if a.status == ACTIVE:
                counts.update(a.rec["reactants"])
        return {species_id(s): float(n) for s, n in counts.items()}

    def run(self, steps: int, samples: int = 500):
        """Run `steps` time steps, yielding the number of steps done: 0, then after steps
        1, 1 + every, 1 + 2 every, ... with every = max(1, steps // samples), and at the
        end (or at extinction). A time step ends with the energy influx."""
        every = max(1, steps // samples) if steps else 1
        dominant = None
        yield 0
        t = 0
        for t in range(steps):
            self.step()
            if not self.now:
                break
            sample = t % every == 0
            if sample:
                counts = Counter(a.seq() for a in self.now)
                top = counts.most_common(1)[0][0]
                if top != dominant:
                    dominant = top
                    self.seen.setdefault(top, None)
                    self.epochs.append([t + 1, top])
            self.energy += self.estep
            if sample:
                self.steps_run = t + 1
                yield t + 1
        done = t + 1 if steps else 0
        if done != self.steps_run:
            self.steps_run = done
            yield done


# --- generate ------------------------------------------------------------------
def _buffered(rng, size=8192):
    buf: list[float] = []

    def rand() -> float:
        if not buf:
            buf.extend(rng.random(size).tolist()[::-1])
        return buf.pop()

    return rand


def _parse_molecules(value, max_length) -> dict[bytes, int]:
    if not isinstance(value, dict) or not value:
        raise ValueError("molecules must be a non-empty mapping {sequence: count}")
    out = {}
    for seq, n in value.items():
        if not isinstance(seq, str) or not seq or any(c not in KEY for c in seq):
            raise ValueError(f"molecules: {seq!r} is not a string over the 33 Stringmol symbols {KEY}")
        if len(seq) > max_length:
            raise ValueError(f"molecules: {seq!r} is longer than max_length = {max_length}")
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise ValueError(f"molecules: count of {seq!r} must be a positive integer, got {n!r}")
        out[seq.encode()] = n
    return out


def _species(seqs) -> list[Species]:
    return [Species(species_id(s), structure=s.decode()) for s in dict.fromkeys(seqs)]


def _ids(seqs) -> list[str]:
    return [species_id(s) for s in seqs]


def _machine(p, rng, exact: bool) -> Machine:
    return Machine(_buffered(rng), max_length=p.max_length, traceback=p.traceback,
                   substitution_rate=0.0 if exact else p.substitution_rate,
                   indel_rate=0.0 if exact else p.indel_rate)


def generate(p, rng):
    """The closure of the distinct seed sequences, with exact copying."""
    molecules = _parse_molecules(p.molecules, p.max_length)
    machine = _machine(p, rng, exact=True)
    stuck: list = []

    def react(a, b):
        out = machine.react(a, b, "possible", p.max_exec_steps)
        if out is False:
            stuck.append([species_id(a), species_id(b)])
            return None
        return out

    seqs, found, status = expand(react, list(molecules), arity=2, max_species=p.max_species, ordered=True)
    return Network(
        species=_species(seqs),
        reactions=[Reaction.of(_ids(lhs), _ids(rhs)) for lhs, rhs in found],
        status=status,
        extras={"seed": _ids(molecules), "non_terminating": stuck},
    )


def evolve(p, rng):
    """The container (a frame about every steps/500 time steps) or the soup (a frame per
    generation of collisions)."""
    if p.agent_radius > p.cell_radius:
        raise ValueError(f"agent_radius ({p.agent_radius}) must not exceed cell_radius ({p.cell_radius})")
    molecules = _parse_molecules(p.molecules, p.max_length)
    machine = _machine(p, rng, exact=False)
    initial = [s for s, n in molecules.items() for _ in range(n)]
    if p.reactor == "container":
        return (yield from _container(p, machine, molecules, initial))
    if len(initial) < 2:
        raise ValueError("the soup needs at least 2 molecules in total")
    return (yield from _soup(p, machine, molecules, initial, rng))


def _container(p, machine, molecules, initial):
    box = Container(machine, initial, energy_per_step=p.energy_per_step, cell_radius=p.cell_radius,
                    agent_radius=p.agent_radius, decay=p.decay)
    for t in box.run(p.steps):
        yield Frame(t=float(t), state=box.state(), fired=box.tally.flush(),
                    observables={"energy": box.energy,
                                 "complexes": sum(1 for a in box.now if a.status == ACTIVE)})
    final = Counter(a.seq() for a in box.now if a.status == UNBOUND)
    in_progress = [a.rec for a in box.now if a.status == ACTIVE]
    seqs = list(molecules)
    for lhs, rhs, _, _ in box.fired.values():
        seqs += [*lhs, *rhs]
    for (lhs, rel) in box.aborted:
        seqs += [*lhs, *rel]
    seqs += list(box.decayed) + list(final) + [s for _, s in box.epochs]
    for rec in in_progress:
        seqs += [*rec["reactants"], *rec["released"]]
    reactions, active = [], []
    for lhs, rhs, count, enzymes in box.fired.values():
        reactions.append(Reaction.of(_ids(lhs), _ids(rhs), count=count))
        active.append({species_id(s): n for s, n in enzymes.items()})
    return Network(
        species=_species(seqs),
        reactions=reactions,
        status="observed",
        initial_state={species_id(s): n for s, n in molecules.items()},
        outflow=float(p.decay) if p.decay > 0 else None,
        extras={
            "epochs": [[t, species_id(s)] for t, s in box.epochs],
            "time_steps": box.steps_run,
            "extinct": not box.now,
            "active_counts": active,
            "final_state": {species_id(s): n for s, n in final.most_common()},
            "decayed": {species_id(s): n for s, n in box.decayed.most_common()},
            "aborted": [{"reactants": _ids(lhs), "released": _ids(rel), "count": n}
                        for (lhs, rel), n in box.aborted.items()],
            "in_progress": [{"reactants": _ids(r["reactants"]), "released": _ids(r["released"])}
                            for r in in_progress],
        },
    )


def _soup(p, machine, molecules, initial, rng):
    stuck = Counter()

    def react(a, b):
        out = machine.react(a, b, "sample", p.max_exec_steps)
        if out is False:
            stuck[(a, b)] += 1
            return None
        return out

    tally = Tally()
    pop = initial
    for step, pop, tally in stir(react, initial, p.steps, rng, arity=2, dilution="constant", tally=tally):
        yield Frame(t=float(step), state={species_id(s): float(n) for s, n in Counter(pop).items()},
                    fired=[[_ids(lhs), _ids(rhs), n] for lhs, rhs, n in tally.flush()],
                    observables={"non_terminating": sum(stuck.values())})
    fired = tally.reactions()
    seqs = list(molecules) + [s for lhs, rhs, _ in fired for s in (*lhs, *rhs)] + pop
    return Network(
        species=_species(seqs),
        reactions=[Reaction.of(_ids(lhs), _ids(rhs), count=n) for lhs, rhs, n in fired],
        status="observed",
        initial_state={species_id(s): n for s, n in molecules.items()},
        outflow=CONSTANT_TOTAL,
        extras={
            "final_state": {species_id(s): n for s, n in Counter(pop).most_common()},
            "non_terminating": sum(stuck.values()),
        },
    )
