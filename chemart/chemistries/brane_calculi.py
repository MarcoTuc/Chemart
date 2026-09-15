"""Brane calculi: Cardelli (2004), book 9.7 [157]. Catalog id: brane-calculi.

A system is a nest of membranes; computation happens *on* membranes. Syntax
(Cardelli 2004, sections 2.1, 3.1, 3.2, 4.1), in the ASCII notation read and
written here:

    systems  P, Q ::= 0 | P, Q | !P | s[P] | m          (0 is the empty system)
    branes   s, t ::= 0 | s|t | !s | a.s
    actions  a    ::= phago_n | cophago_n(r) | exo_n | coexo_n | pino(r)
                    | mate_n | comate_n | bud_n | cobud_n(r) | drip(r)
                    | p1(p2)=>q1(q2)                     bind & release
                    | X                                  inert placeholder

``,`` is the paper's system composition ∘, ``[ ]`` its membrane brackets,
``co...`` its co-actions (⊥), ``=>`` its bind&release arrow ⇉. The subscript
``_n`` pairs an action with its co-action (omitted names match omitted
names). ``a.s|t`` is ``(a.s)|t``; ``!s|t`` is ``(!s)|t``; ``s|t[P]`` is
``(s|t)[P]``; a replicated membrane is written ``!(s[P])``. Molecule lists are
read greedily, so after a molecule a brane that starts with a bind&release
must be parenthesised: ``m, (a=>b)[P]`` (``m, a=>b[P]`` is one membrane). ``p1(p2)=>q1(q2)``
binds molecules p1 outside and p2 inside and releases q1 outside and q2
inside; ``p=>q`` omits empty inner parts. Molecule multisets are comma lists.
Any other identifier in a brane is an inert placeholder (the paper's
metavariables σ0, X, ρ); in a system it is a molecule. The Unicode forms
∘ ⟨ ⟩ ⦇ ⦈ ⇉ ◇ ≜ are accepted.

Reactions (Phago, Exo, Pino in 3.1; Mate, Bud, Drip as primitives, 3.2;
B&R, 4.1), closed under P→Q ⇒ P,R→Q,R and s[P]→s[Q] and structural
congruence (commutative monoids with replication, 0[0] ≡ 0):

    phago_n.s|s0[P], cophago_n(r).t|t0[Q]  ->  t|t0[r[s|s0[P]], Q]
    coexo_n.t|t0[exo_n.s|s0[P], Q]          ->  P, s|s0|t|t0[Q]
    pino(r).s|s0[P]                         ->  s|s0[r[], P]
    mate_n.s|s0[P], comate_n.t|t0[Q]        ->  s|s0|t|t0[P, Q]
    cobud_n(r).t|t0[bud_n.s|s0[P], Q]       ->  r[s|s0[P]], t|t0[Q]
    drip(r).s|s0[P]                         ->  r[], s|s0[P]
    p1, (p1(p2)=>q1(q2)).s|s0[p2, P]        ->  q1, s|s0[q2, P]

A program is a list of definitions ``name := term`` (the paper's ≜; a line
not containing ``:=`` continues the previous definition, ``#`` starts a
comment) and must define ``main``, the initial system. Definitions are
expanded by name in brane or system position; recursion is rejected.

A species is a whole configuration: the canonical text of the normalised
system (composition sorted, units dropped, ``!`` distributed over
composition and absorbing plain copies). Each reduction step is a
unimolecular reaction configuration -> configuration, and the network is the
closure of ``main`` under the steps (chemart.expand.expand with
alternatives=True), truncated by ``max_species``.
"""

from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

SYSTEMS = (
    "viral-infection", "viral-reproduction", "viral-replication",
    "eat-me", "seek-and-store", "plant-vacuole", "custom",
)

NAMED = ("phago", "cophago", "exo", "coexo", "mate", "comate", "bud", "cobud", "pino", "drip")
WITH_MEMBRANE = ("cophago", "cobud", "pino", "drip")
_KEYWORD = re.compile(r"(cophago|phago|coexo|exo|comate|mate|cobud|bud|pino|drip)(?:_([A-Za-z0-9_'+\-]+))?")
_IDENT = re.compile(r"[A-Za-z0-9_'+\-]+")
_TOKEN = re.compile(r"\s+|:=|=>|[()\[\],|.!]|[A-Za-z0-9_'+\-]+")
_UNICODE = {"∘": ",", "⟨": "[", "⦇": "[", "⟩": "]", "⦈": "]", "⇉": "=>", "◇": "0", "≜": ":="}


# --- terms ------------------------------------------------------------------------
class _Term:
    """A normalised term; equality and hashing go through its canonical text."""

    __slots__ = ("text",)

    def __eq__(self, other):
        return type(self) is type(other) and self.text == other.text

    def __hash__(self):
        return hash(self.text)

    def __repr__(self):
        return f"{type(self).__name__}({self.text!r})"


class Mol(_Term):
    __slots__ = ("name",)

    def __init__(self, name: str):
        self.name = name
        self.text = name


class Item(_Term):
    """One action with its continuation, a.s."""

    __slots__ = ("action", "cont")

    def __init__(self, action: tuple, cont: Brane):
        self.action, self.cont = action, cont
        head = action_text(action)
        if cont.size == 0:
            self.text = head
        elif cont.size == 1:
            self.text = f"{head}.{cont.text}"
        else:
            self.text = f"{head}.({cont.text})"


class Brane(_Term):
    """Multiset of plain items plus the set of replicated items."""

    __slots__ = ("plain", "rep", "size")

    def __init__(self, plain: tuple, rep: tuple):
        self.plain, self.rep = plain, rep          # ((Item, count), ...), (Item, ...)
        self.size = sum(c for _, c in plain) + len(rep)
        parts = [it.text for it, c in plain for _ in range(c)] + ["!" + it.text for it in rep]
        self.text = "|".join(parts) if parts else "0"


class Mem(_Term):
    __slots__ = ("brane", "contents")

    def __init__(self, brane: Brane, contents: System):
        self.brane, self.contents = brane, contents
        head = brane.text if brane.size else ""
        first = brane.plain[0][0] if brane.plain else (brane.rep[0] if brane.rep else None)
        if first is not None and first.action[0] == "br" and first.action[1]:
            head = f"({head})"      # 'm,(a=>b)[P]': a leading p1 list would swallow a preceding molecule
        self.text = f"{head}[{contents.text if contents.size else ''}]"


class System(_Term):
    """Multiset of plain elements (molecules, membranes) plus the replicated ones."""

    __slots__ = ("plain", "rep", "size")

    def __init__(self, plain: tuple, rep: tuple):
        self.plain, self.rep = plain, rep
        self.size = sum(c for _, c in plain) + len(rep)
        parts = [e.text for e, c in plain for _ in range(c)]
        parts += ["!" + (e.text if isinstance(e, Mol) else f"({e.text})") for e in rep]
        self.text = ",".join(parts) if parts else "0"


def action_text(a: tuple) -> str:
    kind = a[0]
    if kind == "const":
        return a[1]
    if kind == "br":
        p1, p2, q1, q2 = (",".join(x) for x in a[1:])
        return f"{p1}({p2})=>{q1}({q2})"
    head = kind + (f"_{a[1]}" if a[1] else "")
    return f"{head}({a[2].text})" if kind in WITH_MEMBRANE else head


def _sorted_counts(counter: Counter) -> tuple:
    return tuple(sorted(counter.items(), key=lambda kv: kv[0].text))


def mk_brane(parts) -> Brane:
    """Compose (Brane | Item, replicated) parts into a normalised brane."""
    plain, rep = Counter(), set()
    for x, replicated in parts:
        if isinstance(x, Brane):
            for it, c in x.plain:
                if replicated:
                    rep.add(it)
                else:
                    plain[it] += c
            rep.update(x.rep)
        elif replicated:
            rep.add(x)
        else:
            plain[x] += 1
    for it in rep:
        plain.pop(it, None)
    return Brane(_sorted_counts(plain), tuple(sorted(rep, key=lambda t: t.text)))


def mk_system(parts) -> System:
    """Compose (System | Mol | Mem, replicated) parts into a normalised system."""
    plain, rep = Counter(), set()
    for x, replicated in parts:
        if isinstance(x, System):
            for e, c in x.plain:
                if replicated:
                    rep.add(e)
                else:
                    plain[e] += c
            rep.update(x.rep)
        elif isinstance(x, Mem) and x.brane.size == 0 and x.contents.size == 0:
            continue                                      # 0[0] ≡ 0
        elif replicated:
            rep.add(x)
        else:
            plain[x] += 1
    for e in rep:
        plain.pop(e, None)
    return System(_sorted_counts(plain), tuple(sorted(rep, key=lambda t: t.text)))


EMPTY_BRANE = Brane((), ())
EMPTY = System((), ())


def compose(*branes: Brane) -> Brane:
    return mk_brane([(b, False) for b in branes])


def add(*parts) -> System:
    return mk_system([(x, False) for x in parts])


def molecules(names) -> System:
    return add(*(Mol(n) for n in names))


# --- parsing --------------------------------------------------------------------------
class _Fail(Exception):
    def __init__(self, pos: int, message: str):
        super().__init__(message)
        self.pos, self.message = pos, message


def tokenize(text: str) -> list[str]:
    for u, a in _UNICODE.items():
        text = text.replace(u, a)
    out, pos = [], 0
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m:
            raise ValueError(f"unexpected character {text[pos]!r} in {text!r}")
        if not m.group().isspace():
            out.append(m.group())
        pos = m.end()
    return out


def _is_ident(t) -> bool:
    return isinstance(t, str) and bool(_IDENT.fullmatch(t))


def _is_molecule_name(t) -> bool:
    return _is_ident(t) and t != "0" and not _KEYWORD.fullmatch(t)


class _Parser:
    def __init__(self, tokens, defs, cache, stack, where):
        self.toks, self.defs, self.cache, self.stack, self.where = tokens, defs, cache, stack, where
        self.pos = 0
        self.far = _Fail(0, "empty input")

    def peek(self, k=0):
        i = self.pos + k
        return self.toks[i] if i < len(self.toks) else None

    def fail(self, message):
        f = _Fail(self.pos, message)
        if f.pos >= self.far.pos:
            self.far = f
        raise f

    def expect(self, tok):
        if self.peek() != tok:
            self.fail(f"expected {tok!r}, found {self.peek()!r}")
        self.pos += 1

    def parse_all(self, mode):
        try:
            result = self.system() if mode == "system" else self.brane()
            if self.pos != len(self.toks):
                self.fail(f"unexpected {self.peek()!r}")
            return result
        except _Fail:
            f = self.far
            near = " ".join(self.toks[max(0, f.pos - 3): f.pos + 3])
            raise _Fail(f.pos, f"cannot parse {self.where} as a {mode}: {f.message} near '{near}'") from None

    def expand(self, name, mode):
        key = (name, mode)
        if key in self.cache:
            hit = self.cache[key]
            if isinstance(hit, _Fail):
                raise hit
            return hit
        if name in self.stack:
            raise ValueError(f"recursive definition: {' -> '.join(self.stack + [name])}")
        sub = _Parser(self.defs[name], self.defs, self.cache, self.stack + [name], f"definition {name!r}")
        try:
            result = sub.parse_all(mode)
        except _Fail as f:
            self.cache[key] = f
            raise
        self.cache[key] = result
        return result

    # branes
    def brane(self) -> Brane:
        parts = [self.item()]
        while self.peek() == "|":
            self.pos += 1
            parts.append(self.item())
        return compose(*parts)

    def _is_br(self) -> bool:
        t = self.peek()
        if t == "=>":
            return True
        if t == "(":
            depth, i = 0, self.pos
            while i < len(self.toks):
                depth += {"(": 1, ")": -1}.get(self.toks[i], 0)
                if depth == 0:
                    return i + 1 < len(self.toks) and self.toks[i + 1] == "=>"
                i += 1
            return False
        return _is_molecule_name(t) and self.peek(1) in ("(", ",", "=>")

    def mols(self) -> tuple:
        names = []
        if _is_molecule_name(self.peek()):
            names.append(self.peek())
            self.pos += 1
            while self.peek() == "," and _is_molecule_name(self.peek(1)):
                names.append(self.peek(1))
                self.pos += 2
        return tuple(sorted(names))

    def bind_release(self) -> tuple:
        p1 = self.mols()
        p2 = ()
        if self.peek() == "(":
            self.pos += 1
            p2 = self.mols()
            self.expect(")")
        self.expect("=>")
        q1 = self.mols()
        q2 = ()
        if self.peek() == "(":
            self.pos += 1
            q2 = self.mols()
            self.expect(")")
        return ("br", p1, p2, q1, q2)

    def item(self) -> Brane:
        t = self.peek()
        if t == "!":
            self.pos += 1
            return mk_brane([(self.item(), True)])
        if t == "0":
            self.pos += 1
            return EMPTY_BRANE
        if self._is_br():
            action = self.bind_release()
        elif t == "(":
            self.pos += 1
            b = self.brane()
            self.expect(")")
            return b
        elif _is_ident(t):
            if t in self.defs:
                self.pos += 1
                return self.expand(t, "brane")
            m = _KEYWORD.fullmatch(t)
            self.pos += 1
            if m:
                kind, name = m.group(1), m.group(2) or ""
                if kind in WITH_MEMBRANE:
                    self.expect("(")
                    rho = self.brane()
                    self.expect(")")
                    action = (kind, name, rho)
                else:
                    action = (kind, name)
            else:
                action = ("const", t)
        else:
            self.fail(f"expected an action, found {t!r}")
        cont = EMPTY_BRANE
        if self.peek() == ".":
            self.pos += 1
            cont = self.item()
        return mk_brane([(Item(action, cont), False)])

    # systems
    def system(self) -> System:
        if self.peek() in (None, "]", ")"):
            return EMPTY
        parts = [self.elem()]
        while self.peek() == ",":
            self.pos += 1
            parts.append(self.elem())
        return add(*parts)

    def elem(self) -> System:
        save = self.pos
        brane = EMPTY_BRANE
        if self.peek() != "[":
            try:
                brane = self.brane()
                if self.peek() != "[":
                    self.fail("a brane must be followed by '['")
            except _Fail:
                self.pos = save
                return self.plain_elem()
        self.expect("[")
        contents = self.system()
        self.expect("]")
        return add(Mem(brane, contents))

    def plain_elem(self) -> System:
        t = self.peek()
        if t == "!":
            self.pos += 1
            return mk_system([(self.elem(), True)])
        if t == "0":
            self.pos += 1
            return EMPTY
        if t == "(":
            self.pos += 1
            s = self.system()
            self.expect(")")
            return s
        if _is_ident(t):
            if t in self.defs:
                self.pos += 1
                return self.expand(t, "system")
            if _KEYWORD.fullmatch(t):     # no other reading exists: report it, not the brane attempt
                raise ValueError(f"cannot parse {self.where}: action {t!r} must be on a membrane, s[P]")
            self.pos += 1
            return add(Mol(t))
        self.fail(f"expected a molecule or a membrane, found {t!r}")


def _parse(text: str, defs: dict, mode: str, where: str):
    parser = _Parser(tokenize(text), defs, {}, [], where)
    try:
        return parser.parse_all(mode)
    except _Fail as f:
        raise ValueError(f.message) from None


def parse_brane(text: str, definitions: dict[str, str] | None = None) -> Brane:
    return _parse(text, _def_tokens(definitions or {}), "brane", repr(text))


def parse_system(text: str, definitions: dict[str, str] | None = None) -> System:
    """Parse a system term, optionally with definitions {name: term text}."""
    return _parse(text, _def_tokens(definitions or {}), "system", repr(text))


def _def_tokens(definitions: dict[str, str]) -> dict[str, list[str]]:
    out = {}
    for name, body in definitions.items():
        if not _is_molecule_name(name):
            raise ValueError(f"definition name {name!r} must be an identifier and not an action keyword")
        out[name] = tokenize(body)
    return out


def parse_definitions(program: str) -> dict[str, str]:
    """``name := term`` lines (continuation lines without ':=', '#' comments)."""
    defs: dict[str, str] = {}
    current = None
    for raw in program.replace("≜", ":=").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if ":=" in line:
            name, body = (s.strip() for s in line.split(":=", 1))
            if name in defs:
                raise ValueError(f"{name!r} is defined twice")
            defs[name] = body
            current = name
        elif current is None:
            raise ValueError(f"program line {raw.strip()!r} is not a definition 'name := term'")
        else:
            defs[current] += " " + line
    return defs


def parse_program(program: str) -> tuple[dict[str, str], System]:
    defs = parse_definitions(program)
    if "main" not in defs:
        raise ValueError("the program must define the initial system: 'main := ...'")
    tokens = _def_tokens(defs)
    parser = _Parser(tokens["main"], tokens, {}, ["main"], "definition 'main'")
    try:
        return defs, parser.parse_all("system")
    except _Fail as f:
        raise ValueError(f.message) from None


# --- reduction ------------------------------------------------------------------------
def _picks(s: System):
    """(element, rest of the system): a plain copy, or an unfolding of a replicated one."""
    for e, c in s.plain:
        rest = tuple((x, n - (x == e)) for x, n in s.plain if x != e or n > 1)
        yield e, System(rest, s.rep)
    for e in s.rep:
        yield e, s


def _actions(b: Brane):
    """(action, continuation, rest of the brane) for every action ready to fire."""
    for it, c in b.plain:
        rest = tuple((x, n - (x == it)) for x, n in b.plain if x != it or n > 1)
        yield it.action, it.cont, Brane(rest, b.rep)
    for it in b.rep:
        yield it.action, it.cont, b


def _remove(s: System, names) -> System | None:
    counts = Counter(names)
    plain = dict(s.plain)
    rep = set(s.rep)
    for name, k in counts.items():
        m = Mol(name)
        if m in rep:
            continue
        if plain.get(m, 0) < k:
            return None
        plain[m] -= k
    return System(tuple((e, c) for e, c in plain.items() if c), s.rep)


@lru_cache(maxsize=200_000)
def steps(s: System) -> tuple[tuple[System, str], ...]:
    """Every configuration reachable in one reduction, with the rule that fired."""
    out: dict[System, str] = {}

    def emit(t: System, rule: str):
        if t != s:
            out.setdefault(t, rule)

    for e, rest in _picks(s):
        if not isinstance(e, Mem):
            continue
        for a, cont, others in _actions(e.brane):
            kind, here = a[0], compose(cont, others)                     # s|s0
            if kind in ("phago", "mate"):
                partner = "cophago" if kind == "phago" else "comate"
                for e2, rest2 in _picks(rest):
                    if not isinstance(e2, Mem):
                        continue
                    for a2, cont2, others2 in _actions(e2.brane):
                        if a2[0] != partner or a2[1] != a[1]:
                            continue
                        there = compose(cont2, others2)                  # t|t0
                        if kind == "phago":
                            inner = Mem(a2[2], add(Mem(here, e.contents)))
                            emit(add(rest2, Mem(there, add(inner, e2.contents))), "phago")
                        else:
                            emit(add(rest2, Mem(compose(here, there), add(e.contents, e2.contents))), "mate")
            elif kind in ("coexo", "cobud"):
                partner = "exo" if kind == "coexo" else "bud"
                for e2, q in _picks(e.contents):
                    if not isinstance(e2, Mem):
                        continue
                    for a2, cont2, others2 in _actions(e2.brane):
                        if a2[0] != partner or a2[1] != a[1]:
                            continue
                        inner_brane = compose(cont2, others2)            # s|s0 of the inner membrane
                        if kind == "coexo":
                            emit(add(rest, e2.contents, Mem(compose(inner_brane, here), q)), "exo")
                        else:
                            bud = Mem(a[2], add(Mem(inner_brane, e2.contents)))
                            emit(add(rest, bud, Mem(here, q)), "bud")
            elif kind == "pino":
                emit(add(rest, Mem(here, add(Mem(a[2], EMPTY), e.contents))), "pino")
            elif kind == "drip":
                emit(add(rest, Mem(a[2], EMPTY), Mem(here, e.contents)), "drip")
            elif kind == "br":
                _, p1, p2, q1, q2 = a
                outside, inside = _remove(rest, p1), _remove(e.contents, p2)
                if outside is not None and inside is not None:
                    emit(add(outside, molecules(q1), Mem(here, add(inside, molecules(q2)))), "bind&release")
        for t, rule in steps(e.contents):
            emit(add(rest, Mem(e.brane, t)), rule)
    return tuple(out.items())


def reachable(start: System, max_species: int = 1000):
    """Closure of `start`: (configurations, [(lhs, rhs, rule)], status)."""
    def react(s):
        nxt = steps(s)
        return [(t,) for t, _ in nxt] if nxt else None

    species, found, status = expand(react, [start], arity=1, max_species=max_species, alternatives=True)
    rules = []
    for (lhs,), (rhs,) in found:
        rules.append((lhs, rhs, dict(steps(lhs))[rhs]))
    return species, rules, status


# --- inspection -------------------------------------------------------------------------
def layout(s: System, indent: int = 0) -> str:
    """Indented tree: one line per molecule or membrane (brane text), contents below."""
    pad = "  " * indent
    lines = []
    for e, c, rep in [(e, c, False) for e, c in s.plain] + [(e, 1, True) for e in s.rep]:
        mult = "!" if rep else (f"{c} x " if c > 1 else "")
        if isinstance(e, Mol):
            lines.append(f"{pad}{mult}{e.name}")
        else:
            lines.append(f"{pad}{mult}membrane {e.brane.text}")
            if e.contents.size:
                lines.append(layout(e.contents, indent + 1))
    return "\n".join(lines) if lines else pad + "0"


def tree(s: System) -> list[dict]:
    """JSON nesting of a configuration."""
    out = []
    for e, c, rep in [(e, c, False) for e, c in s.plain] + [(e, 1, True) for e in s.rep]:
        node = {"molecule": e.name} if isinstance(e, Mol) else {
            "membrane": e.brane.text, "contents": tree(e.contents)}
        node["copies"] = "replicated" if rep else c
        out.append(node)
    return out


def membranes(s: System) -> int:
    """Number of membranes, nested ones included (a replicated membrane counts once)."""
    elems = [(e, c) for e, c in s.plain] + [(e, 1) for e in s.rep]
    return sum(c * (1 + membranes(e.contents)) for e, c in elems if isinstance(e, Mem))


def molecule_parities(s: System, depth: int = 0) -> Counter:
    """Counter of (molecule, nesting depth mod 2): invariant under bitonal reactions."""
    out = Counter()
    for e, c in s.plain:
        if isinstance(e, Mol):
            out[(e.name, depth % 2)] += c
        else:
            for key, n in molecule_parities(e.contents, depth + 1).items():
                out[key] += n * c
    for e in s.rep:
        if isinstance(e, Mol):
            out[("!" + e.name, depth % 2)] = 1
        else:
            for (name, parity), _ in molecule_parities(e.contents, depth + 1).items():
                out[("!" + name.lstrip("!"), parity)] = 1
    return out


def copies(s: System, e) -> int:
    """Copies of element `e` at the top level (-1 when replicated)."""
    if e in s.rep:
        return -1
    return dict(s.plain).get(e, 0)


# --- the published systems -----------------------------------------------------------------
VIRAL_PART1 = """\
# Cardelli 2004, section 3.3 (Semliki Forest virus)
virus            := phago.exo[nucap]
nucap            := !bud|X[vRNA]
cell             := membrane[cytosol]
membrane         := !cophago(mate)|!coexo
cytosol          := endosome, Z
endosome         := !comate|!coexo[]
viral-envelope   := cobud(phago.exo)
envelope-vesicle := exo.viral-envelope[]
"""

PROGRAMS = {
    "viral-infection": VIRAL_PART1 + "main := virus, cell\n",
    "viral-reproduction": VIRAL_PART1 + "main := membrane[nucap, envelope-vesicle, Z']\n",
    "viral-replication": """\
# Cardelli 2004, section 4.6 (nucleocapsid replication, figure 11)
nucap            := capsid[vRNA]
capsid           := !bud|disasm
disasm           := disasm-trigger(vRNA)=>vRNA
vRNA-repl        := !vRNA=>vRNA,vRNA[]            # vRNA -> vRNA, vRNA (section 4.2)
capsomer-tran    := !vRNA=>vRNA.drip(capsomers)[]
capsomers        := vRNA=>(vRNA).capsid
ER               := !vRNA=>vRNA.drip(exo.viral-envelope)[Nucleus]
cytosol          := endosome, !disasm-trigger, vRNA-repl, capsomer-tran, ER
endosome         := !comate|!coexo[]
viral-envelope   := cobud(phago.exo)
envelope-vesicle := exo.viral-envelope[]
main             := nucap, cytosol
""",
    "eat-me": """\
# Cardelli 2004, section 4.5: A releases n, which makes B eat A
A    := =>n.phago[P]
B    := n=>.cophago(rho)[Q]
main := A, B
""",
    "seek-and-store": """\
# Cardelli 2004, section 4.5: cell C takes in nutrient n and stores it in a vesicle
seek_n := !n=>.pino(=>(n).mate_store)
store  := !comate_store
C      := seek_n[store[]]
main   := n, C
""",
    "plant-vacuole": """\
# Cardelli 2004, section 4.4 (figure 10); the outside molecules are a Chemart choice
ProtonPump       := !ATP=>ADP,Pi(H+,H+)
IonChannel       := !Cl-(H+)=>(H+,Cl-)
ProtonAntiporter := !Na+(H+)=>H+(Na+)
PlantVacuole     := ProtonPump|IonChannel|ProtonAntiporter[]
main             := ATP, Cl-, Na+, PlantVacuole
""",
}

#: Configurations the paper derives from `main` (terms over the program's definitions).
PUBLISHED = {
    "viral-infection": {"infection (fig. 7)": "membrane[nucap, cytosol]"},
    "viral-reproduction": {"reproduction (fig. 8)": "membrane[Z'], virus"},
    "viral-replication": {"disassembly (sec. 4.6)": "cytosol, vRNA, !bud[]"},
    "eat-me": {"A eaten by B (sec. 4.5)": "[rho[[P]], Q]"},
    "seek-and-store": {"n stored (sec. 4.5)": "seek_n[store[n]]"},
    "plant-vacuole": {},
}


def program_for(p) -> str:
    if p.system == "custom":
        if not p.program.strip():
            raise ValueError("system='custom' needs program: definitions 'name := term' including 'main := ...'")
        return p.program
    if p.program.strip():
        raise ValueError("program is only used with system='custom'")
    return PROGRAMS[p.system]


# --- generator ---------------------------------------------------------------------------------
def generate(p, rng):
    program = program_for(p)
    defs, start = parse_program(program)
    configs, found, status = reachable(start, p.max_species)

    ids = {s: s.text for s in configs}
    species = [Species(s.text, structure=layout(s)) for s in configs]
    reactions = [Reaction.of([ids[lhs]], [ids[rhs]]) for lhs, rhs, _ in found]

    analysis: dict = {"membranes_initial": membranes(start)}
    known = set(configs)
    for label, text in PUBLISHED.get(p.system, {}).items():
        target = parse_system(text, defs)
        analysis[label] = {"configuration": target.text, "reached": target in known}
    if status == "complete":
        sources = {lhs for lhs, _, _ in found}
        analysis["terminal"] = [s.text for s in configs if s not in sources]
    if p.system == "viral-replication":
        nucap, vesicle = parse_system("nucap", defs), parse_system("envelope-vesicle", defs)
        (n_el,), (v_el,) = [e for e, _ in nucap.plain], [e for e, _ in vesicle.plain]
        analysis["max_nucap_copies"] = max(copies(s, n_el) for s in configs)
        analysis["max_envelope_vesicle_copies"] = max(copies(s, v_el) for s in configs)

    extras = {
        "system": p.system,
        "program": program,
        "definitions": defs,
        "main": start.text,
        "reaction_rules": [rule for _, _, rule in found],
        "compartments": {
            "representation": "each species is a whole configuration (canonical nested-membrane text); "
                              "this is the membrane tree of the initial one",
            "initial": tree(start),
        },
        "analysis": analysis,
    }
    return Network(species=species, reactions=reactions, status=status,
                   initial_state={start.text: 1.0}, extras=extras)
