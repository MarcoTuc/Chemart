"""The Chemical Abstract Machine (CHAM). Catalog id: cham.

Book 9.3; G. Berry & G. Boudol, "The chemical abstract machine", POPL 1990,
pp. 81-94 (journal version: Theoretical Computer Science 96(1):217-248, 1992,
book [113]).

A CHAM is given by molecules (terms of an algebra), solutions (finite
multisets of molecules, written {| m1, ..., mk |}) and transformation rules
m1, ..., mk -> m1', ..., ml'. A solution is itself a molecule (the membrane
operator), so solutions nest. Rules are classified as heating (a molecule is
broken into simpler ones), cooling (the inverse) and reaction rules
(irreversible, they change the information content). They are applied by
general laws (sec. 3.1):

- reaction law: an instance of a rule's right-hand side replaces the
  corresponding instance of its left-hand side;
- chemical law: reactions can be performed freely within any solution
  (S -> S' gives S + S'' -> S' + S'');
- membrane law: a subsolution can evolve freely in any context
  ({| C[S] |} -> {| C[S'] |});
- airlock law: {| m |} + S <-> {| m <| S |}, a reversible structural
  transition that isolates one molecule next to the rest of the solution.

The molecules are those of the paper's process-calculus CHAMs (sec. 2 and
4.3), in ASCII:

    0                   inaction            P, Q, X    process constants / fix variables
    a.m  ~a.m  tau.m    prefix (ion)        m | m      parallel (agents only)
    m\\a                restriction         m[b/a,~c/d] relabelling (a -> b, d -> ~c)
    m (+) m             internal sum        m [] m     external sum (agents only)
    {| m, ... |}        solution            m <| S     airlock
    <m, m>  l:<m, m>  r:<m, m>              pairs
    fix_X(X=a.Y;Y=b.X)  fixpoint, selecting the definition of X

Postfix restriction and relabelling bind tightest, then prefix, airlock,
parallel and the sums. Species ids are the canonical text without spaces;
the members of a solution are sorted.

Rules are text, one per line, ``name: lhs ARROW rhs [if guard]``, with the
arrow surrounded by spaces:

    <=> (or ⇌)  heating from left to right, paired with its inverse cooling
    =>  (or ⇀)  heating only          ~>  (or ⇁)  cooling only
    ->  (or →)  reaction

Each side is a comma-separated list of molecule patterns (the right side may
be empty). A pattern is a molecule with metavariables ``?name`` whose kind is
fixed by the first letter: ?p ?q agents, ?m ?n any molecule, ?S ?T solutions,
?a ?b ?c ?d names, ?l labels (a, ~a or tau), ?f relabellings, ?F fixpoint
terms. In a prefix, ``?a.?m`` matches a name, ``~?a.?m`` a co-name and
``?l.?m`` any label; ``?f(?l)`` is the label ?l relabelled by ?f, and
``unfold(?F)`` is the one-step unfolding p_i[fix(x=p)/x]. Following the
paper, a solution pattern is ?S, {||} or {| pattern |} (no multiset
matching). The only guard is ``if <label> notin {<label>, ...}`` (or ``in``).
Comments start with '#'. Nothing is evaluated as Python.

The network is the closure of the initial solution (chemart.expand with
alternatives): species are the molecules floating in the top solution, a
reaction is one rule instance, one airlock transition, or a transition inside
a membrane lifted to the enclosing molecule (membrane law). Structural
transitions are recorded in both directions. To keep the closure finite the
exploration follows heating, reaction and cooling-only rules, the cooling
direction of reversible rules that decompose or rearrange a single molecule,
and the airlock law inside membranes on a whole subsolution (one molecule
isolated from all the others); heating that puts a membrane directly around a
membrane ({|m|}\\a with m already a solution) is skipped, and cooling that
joins several molecules (p, q -> p|q) and airlocks built in the top solution
only appear as inverses of observed heatings. ``transitions(..., full=True)`` gives the unrestricted
one-step relation of a finite solution.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

MACHINES = ("tccs", "ccs-minus", "custom")
PROGRAMS = ("heavy-ion-communication", "ccs-minus-execution", "nondeterministic-choice",
            "external-sum", "custom")

#: Berry & Boudol sec. 2.2-2.3: parallel, reaction, restriction membrane and ion, heavy ion.
CCS_MINUS_RULES = r"""
parallel: ?p | ?q <=> ?p, ?q
reaction: ?a.?m, ~?a.?n -> ?m, ?n
restriction membrane: ?m\?a <=> {|?m|}\?a
restriction ion: (?l.?m)\?a <=> ?l.(?m\?a) if ?l notin {?a, ~?a}
heavy ion: ?l.?m <| ?S <=> ?l.(?m <| ?S)
"""

#: Berry & Boudol sec. 4.3, the complete TCCS cham (plus the heavy ion rule of sec. 2.3).
TCCS_RULES = CCS_MINUS_RULES + r"""
relabelling membrane: ?m[?f] <=> {|?m|}[?f]
relabelling ion: (?l.?m)[?f] <=> ?f(?l).(?m[?f])
internal sum left: ?p (+) ?q -> ?p
internal sum right: ?p (+) ?q -> ?q
external sum expansion: ?p [] ?q <=> <{|?p|}, {|?q|}>
left external sum ion: <{|?l.?m|}, ?S> <=> ?l.l:<?m, ?S>
right external sum ion: <?S, {|?l.?m|}> <=> ?l.r:<?S, ?m>
left projection: l:<?m, ?n> ~> ?m
right projection: r:<?m, ?n> ~> ?n
fixpoint: ?F => unfold(?F)
"""

#: "Additional cleanup rules" (sec. 2.2 and 4.3).
CCS_MINUS_CLEANUP = "inaction cleanup: 0 =>\n"
TCCS_CLEANUP = CCS_MINUS_CLEANUP + "restriction cleanup: {||}\\?a =>\nrelabelling cleanup: {||}[?f] =>\n"

#: The paper's worked examples (p, q of the paper are the process constants P, Q).
PROGRAM_SOLUTIONS = {
    "heavy-ion-communication": r"a.0|(~a.P|Q)\b",          # sec. 2.3
    "ccs-minus-execution": "a.b.0|~a.0|~b.0",               # sec. 2.2
    "nondeterministic-choice": "a.0, ~a.b.0, ~a.c.0",       # sec. 2.2
    "external-sum": "a.~b.0|((~a.0|b.0)[]Q)",               # sec. 4.2
}

KINDS = {"p": "agent", "q": "agent", "m": "mol", "n": "mol", "S": "sol", "T": "sol",
         "a": "name", "b": "name", "c": "name", "d": "name", "l": "label", "f": "relab",
         "F": "fix"}

NIL = ("0",)

# --- terms ------------------------------------------------------------------------
# ("0",) ("K", X) (".", label, m) ("|", p, q) ("[]", p, q) ("(+)", p, q) ("\\", m, a)
# ("rel", m, ((old, new), ...)) ("sol", (m, ...)) ("<|", m, S) ("<>", tag, m, n)
# ("fix", X, ((X, p), ...)); patterns add ("?", kind, name), ("~?", name) and
# ("app", f, label) in label positions, and ("unfold", name).

_PREC = {"[]": 0, "(+)": 0, "|": 1, "<|": 2, ".": 3, "\\": 4, "rel": 4}
_CHILDREN = {".": (2,), "|": (1, 2), "[]": (1, 2), "(+)": (1, 2), "\\": (1,), "rel": (1,),
             "<|": (1, 2), "<>": (2, 3)}


def _sol(items) -> tuple:
    return ("sol", tuple(sorted(items, key=text)))


@lru_cache(maxsize=None)
def text(t) -> str:
    """Canonical text of a molecule or pattern, without spaces: the species id."""
    return _txt(t, False)


def pretty(t) -> str:
    """Readable text with spaces, as in the paper."""
    return _txt(t, True)


def _wrap(t, level: int, sp: bool) -> str:
    s = _txt(t, sp)
    return f"({s})" if _PREC.get(t[0], 5) < level else s


def _label_txt(lab) -> str:
    if isinstance(lab, str):
        return lab
    if lab[0] == "?":
        return "?" + lab[2]
    if lab[0] == "~?":
        return "~?" + lab[1]
    return f"?{lab[1]}({_label_txt(lab[2])})"


def _txt(t, sp: bool) -> str:
    k = t[0]
    comma = ", " if sp else ","
    if k == "0":
        return "0"
    if k == "K":
        return t[1]
    if k == "?":
        return "?" + t[2]
    if k == "unfold":
        return f"unfold(?{t[1]})"
    if k == ".":
        return f"{_label_txt(t[1])}.{_wrap(t[2], 3, sp)}"
    if k in ("|", "[]", "(+)"):
        level = _PREC[k]
        op = f" {k} " if sp else k
        return _wrap(t[1], level, sp) + op + _wrap(t[2], level + 1, sp)
    if k == "<|":
        return _wrap(t[1], 2, sp) + (" <| " if sp else "<|") + _wrap(t[2], 5, sp)
    if k == "\\":
        name = t[2] if isinstance(t[2], str) else "?" + t[2][2]
        return f"{_wrap(t[1], 4, sp)}\\{name}"
    if k == "rel":
        phi = t[2]
        body = "?" + phi[2] if phi and phi[0] == "?" else ",".join(f"{new}/{old}" for old, new in phi)
        return f"{_wrap(t[1], 4, sp)}[{body}]"
    if k == "sol":
        inner = comma.join(_txt(m, sp) for m in t[1])
        return f"{{| {inner} |}}" if sp and inner else f"{{|{inner}|}}"
    if k == "<>":
        tag = f"{t[1]}:" if t[1] else ""
        return f"{tag}<{_txt(t[2], sp)}{comma}{_txt(t[3], sp)}>"
    if k == "fix":
        defs = ";".join(f"{v}={_txt(p, sp)}" for v, p in t[2])
        return f"fix_{t[1]}({defs})"
    raise ValueError(f"not a term: {t!r}")


def is_agent(t) -> bool:
    """Agents of the TCCS syntax: 0, constants, prefix, |, \\, relabelling, sums, fix."""
    k = t[0]
    if k in ("0", "K", "fix"):
        return True
    if k == "?":
        return t[1] in ("agent", "fix")
    if k == ".":
        return is_agent(t[2])
    if k in ("|", "[]", "(+)"):
        return is_agent(t[1]) and is_agent(t[2])
    if k in ("\\", "rel"):
        return is_agent(t[1])
    return False


# --- parser -------------------------------------------------------------------------
_TOKEN = re.compile(r"""\s*(?:
    (?P<fix>fix_[A-Z][A-Za-z0-9_']*)
  | (?P<meta>\?[A-Za-z][A-Za-z0-9_']*)
  | (?P<sym>\{\||\|\}|<\||\(\+\)|\[\]|[|.\\\[\]/,;=()<>:~{}])
  | (?P<zero>0)
  | (?P<name>[a-z][a-z0-9_]*)
  | (?P<const>[A-Z][A-Za-z0-9_']*)
)""", re.X)


def _tokens(src: str) -> list[tuple[str, str]]:
    out, pos, src = [], 0, src.rstrip()
    while pos < len(src):
        m = _TOKEN.match(src, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot read {src[pos:]!r} in {src!r}")
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
        pos = m.end()
    return out


class _Parser:
    def __init__(self, src: str, patterns: bool):
        self.src, self.patterns = src, patterns
        self.toks = _tokens(src)
        self.i = 0

    def error(self, what: str):
        where = " ".join(v for _, v in self.toks[self.i:self.i + 4]) or "end of input"
        return ValueError(f"{what} at '{where}' in {self.src!r}")

    def peek(self, k: int = 0):
        j = self.i + k
        return self.toks[j] if j < len(self.toks) else (None, None)

    def at(self, value: str, k: int = 0) -> bool:
        kind, v = self.peek(k)
        return kind in ("sym", "name") and v == value

    def take(self, value: str | None = None):
        kind, v = self.peek()
        if kind is None or (value is not None and v != value):
            raise self.error(f"expected {value!r}" if value else "unexpected end")
        self.i += 1
        return kind, v

    def done(self) -> bool:
        return self.i >= len(self.toks)

    def meta(self, allowed: tuple[str, ...]):
        kind, v = self.take()
        if kind != "meta":
            raise self.error("expected a metavariable")
        if not self.patterns:
            raise self.error("metavariables are only allowed in rules")
        mk = KINDS.get(v[1])
        if mk is None:
            raise self.error(f"metavariable {v} must start with one of {''.join(KINDS)}")
        if mk not in allowed:
            raise self.error(f"metavariable {v} ({mk}) cannot stand here")
        return ("?", mk, v[1:])

    # molecules
    def molecules(self, stop=()) -> list:
        items = []
        if self.done() or self.peek()[1] in stop:
            return items
        items.append(self.mol())
        while self.at(","):
            self.take(",")
            items.append(self.mol())
        return items

    def mol(self):
        left = self.par()
        while self.at("[]") or self.at("(+)"):
            op = self.take()[1]
            left = self.agents(op, left, self.par())
        return left

    def agents(self, op, left, right):
        if not (is_agent(left) and is_agent(right)):
            raise self.error(f"'{op}' composes agents, not molecules ({text(left)}, {text(right)})")
        return (op, left, right)

    def par(self):
        left = self.air()
        while self.at("|"):
            self.take("|")
            left = self.agents("|", left, self.air())
        return left

    def air(self):
        left = self.pre()
        while self.at("<|"):
            self.take("<|")
            right = self.post()
            if not (right[0] == "sol" or (right[0] == "?" and right[1] == "sol")):
                raise self.error("an airlock m <| S needs a solution on its right")
            left = ("<|", left, right)
        return left

    def label_ahead(self) -> bool:
        kind, v = self.peek()
        if kind == "sym":
            return v == "~"
        if kind == "name":
            return not (v in ("l", "r") and self.at(":", 1)) and not (v == "unfold" and self.at("(", 1))
        if kind == "meta" and self.patterns:
            mk = KINDS.get(v[1])
            return (mk in ("label", "name") and self.at(".", 1)) or (mk == "relab" and self.at("(", 1))
        return False

    def label(self):
        kind, v = self.peek()
        if kind == "sym" and v == "~":
            self.take("~")
            kind, v = self.peek()
            if kind == "name" and v != "tau":
                self.take()
                return "~" + v
            return ("~?", self.meta(("name",))[2])
        if kind == "name":
            self.take()
            return v
        if kind == "meta" and KINDS.get(v[1]) == "relab":
            f = self.meta(("relab",))[2]
            self.take("(")
            inner = self.label()
            self.take(")")
            return ("app", f, inner)
        return self.meta(("label", "name"))

    def pre(self):
        if self.label_ahead():
            lab = self.label()
            self.take(".")
            return (".", lab, self.pre())
        return self.post()

    def name(self):
        kind, v = self.peek()
        if kind == "meta":
            return self.meta(("name",))
        if kind != "name" or v == "tau":
            raise self.error("expected a name")
        self.take()
        return v

    def post(self):
        m = self.atom()
        while self.at("\\") or self.at("["):
            if self.at("\\"):
                self.take("\\")
                m = ("\\", m, self.name())
                continue
            self.take("[")
            if self.peek()[0] == "meta":
                phi = self.meta(("relab",))
            else:
                pairs = {}
                while True:
                    new = self.label()
                    if not isinstance(new, str) or new == "tau":
                        raise self.error("a relabelling maps names to names or co-names")
                    self.take("/")
                    old = self.name()
                    if not isinstance(old, str) or old in pairs:
                        raise self.error(f"bad or repeated name {old!r} in relabelling")
                    pairs[old] = new
                    if not self.at(","):
                        break
                    self.take(",")
                phi = tuple(sorted(pairs.items()))
            self.take("]")
            m = ("rel", m, phi)
        return m

    def atom(self):
        kind, v = self.peek()
        if kind == "zero":
            self.take()
            return NIL
        if kind == "const":
            self.take()
            return ("K", v)
        if kind == "meta":
            return self.meta(("agent", "mol", "sol", "fix"))
        if kind == "fix":
            self.take()
            return self.fix(v[4:])
        if kind == "name" and v == "unfold" and self.at("(", 1):
            if not self.patterns:
                raise self.error("unfold(...) is only allowed in rules")
            self.take()
            self.take("(")
            fvar = self.meta(("fix",))
            self.take(")")
            return ("unfold", fvar[2])
        if kind == "name" and v in ("l", "r") and self.at(":", 1):
            self.take()
            self.take(":")
            return self.pair(v)
        if v == "(":
            self.take("(")
            m = self.mol()
            self.take(")")
            return m
        if v == "<":
            return self.pair("")
        if v == "{|":
            self.take("{|")
            items = self.molecules(stop=("|}",))
            self.take("|}")
            if self.patterns and len(items) > 1:
                raise self.error("a solution pattern is ?S, {||} or {|m|} (no multiset matching)")
            return _sol(items)
        raise self.error("expected a molecule")

    def pair(self, tag: str):
        self.take("<")
        m = self.mol()
        self.take(",")
        n = self.mol()
        self.take(">")
        return ("<>", tag, m, n)

    def fix(self, var: str):
        self.take("(")
        defs = {}
        while True:
            kind, v = self.take()
            if kind != "const" or v in defs:
                raise self.error("expected a new variable X=agent in fix")
            self.take("=")
            body = self.mol()
            if not is_agent(body):
                raise self.error("fix defines agents")
            defs[v] = body
            if not self.at(";"):
                break
            self.take(";")
        self.take(")")
        if var not in defs:
            raise self.error(f"fix_{var} selects an undefined variable")
        return ("fix", var, tuple(sorted(defs.items())))


def parse_molecule(src: str):
    p = _Parser(src, patterns=False)
    m = p.mol()
    if not p.done():
        raise p.error("unexpected text")
    return m


def parse_solution(src: str) -> tuple:
    """'{|a.0, b.0|}' or 'a.0, b.0' -> the sorted tuple of molecules."""
    p = _Parser(src, patterns=False)
    items = p.molecules()
    if not p.done():
        raise p.error("unexpected text")
    if len(items) == 1 and items[0][0] == "sol" and src.strip().startswith("{|"):
        return items[0][1]
    return _sol(items)[1]


def solution_text(items) -> str:
    return text(_sol(items))


# --- rules --------------------------------------------------------------------------
_ARROWS = {"<=>": "reversible", "⇌": "reversible", "=>": "heating", "⇀": "heating",
           "~>": "cooling", "⇁": "cooling", "->": "reaction", "→": "reaction"}
_ARROW = re.compile(r"(?:^|\s)(<=>|=>|~>|->|⇌|⇀|⇁|→)(?=\s|$)")
_FLIP = {"heating": "cooling", "cooling": "heating"}


@dataclass(frozen=True)
class Rule:
    name: str
    kind: str            # reversible | heating | cooling | reaction
    lhs: tuple
    rhs: tuple
    guard: tuple | None
    text: str


@dataclass(frozen=True)
class Direction:
    cls: str             # heating | cooling | reaction
    lhs: tuple
    rhs: tuple
    guard: tuple | None
    rule: str
    reversible: bool


def _metas(t, out: set) -> set:
    if isinstance(t, str):
        return out
    if isinstance(t, tuple) and t:
        if not isinstance(t[0], str):          # a tuple of terms (solution members)
            for x in t:
                _metas(x, out)
            return out
        if t[0] == "?":
            out.add(t[2])
            return out
        if t[0] == "~?":
            out.add(t[1])
            return out
        if t[0] == "app":
            out.add(t[1])
            return _metas(t[2], out)
        if t[0] == "unfold":
            out.add(t[1])
            return out
        for x in t[1:]:
            _metas(x, out)
    return out


def _has(t, head: str) -> bool:
    if not isinstance(t, tuple) or not t:
        return False
    if not isinstance(t[0], str):
        return any(_has(x, head) for x in t)
    return t[0] == head or any(_has(x, head) for x in t[1:])


def _parse_rule(line: str) -> Rule:
    if ": " not in line:
        raise ValueError(f"rule {line!r} must start with 'name: '")
    name, body = line.split(": ", 1)
    name = name.strip()
    found = list(_ARROW.finditer(body))
    if len(found) != 1:
        raise ValueError(f"rule {line!r} needs exactly one arrow (<=> => ~> -> or ⇌ ⇀ ⇁ →) between spaces")
    arrow = found[0]
    kind = _ARROWS[arrow.group(1)]
    left, right = body[:arrow.start(1)], body[arrow.end(1):]
    guard_src = None
    g = re.search(r"(?:^|\s)if\s", right)
    if g:
        right, guard_src = right[:g.start()], right[g.end():]
    lp = _Parser(left, patterns=True)
    lhs = tuple(lp.molecules())
    if not lp.done():
        raise lp.error("unexpected text")
    rp = _Parser(right, patterns=True)
    rhs = tuple(rp.molecules())
    if not rp.done():
        raise rp.error("unexpected text")
    guard = _parse_guard(guard_src) if guard_src is not None else None

    if not lhs:
        raise ValueError(f"rule {name!r} has an empty left-hand side")
    left_vars = set().union(*(_metas(t, set()) for t in lhs))
    right_vars = set().union(set(), *(_metas(t, set()) for t in rhs))
    if not right_vars <= left_vars:
        raise ValueError(f"rule {name!r}: {sorted(right_vars - left_vars)} appear only on the right")
    if guard and not _metas(guard, set()) <= left_vars:
        raise ValueError(f"rule {name!r}: the guard uses metavariables not bound on the left")
    if any(_has(t, "unfold") for t in lhs):
        raise ValueError(f"rule {name!r}: unfold(...) is only allowed on the right")
    if kind == "reversible":
        if not rhs or left_vars != right_vars or any(_has(t, "unfold") for t in rhs):
            raise ValueError(f"reversible rule {name!r} needs both sides with the same metavariables and no unfold")
    canon = f"{name}: {', '.join(pretty(t) for t in lhs)} {arrow.group(1)}"
    if rhs:
        canon += " " + ", ".join(pretty(t) for t in rhs)
    if guard_src is not None:
        canon += " if " + " ".join(guard_src.split())
    return Rule(name, kind, lhs, rhs, guard, canon)


def _parse_guard(src: str):
    p = _Parser(src, patterns=True)
    lab = p.label()
    kind, op = p.take()
    if op not in ("in", "notin"):
        raise p.error("a guard is 'if <label> notin {<label>, ...}' or 'in'")
    p.take("{")
    labels = [p.label()]
    while p.at(","):
        p.take(",")
        labels.append(p.label())
    p.take("}")
    if not p.done():
        raise p.error("unexpected text")
    return (op, lab, tuple(labels))


def parse_rules(src: str) -> list[Rule]:
    rules = []
    for raw in src.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            rules.append(_parse_rule(line))
    names = [r.name for r in rules]
    dup = [n for n, c in Counter(names).items() if c > 1]
    if dup:
        raise ValueError(f"rule names must be unique, repeated: {dup}")
    return rules


def machine_rules(machine: str, cleanup: bool = True, rules: str = "") -> list[Rule]:
    if machine == "custom":
        if not rules.strip():
            raise ValueError("machine=custom needs rules (one 'name: lhs -> rhs' per line)")
        return parse_rules(rules)
    if rules.strip():
        raise ValueError("rules is only used with machine=custom")
    if machine == "ccs-minus":
        return parse_rules(CCS_MINUS_RULES + (CCS_MINUS_CLEANUP if cleanup else ""))
    if machine == "tccs":
        return parse_rules(TCCS_RULES + (TCCS_CLEANUP if cleanup else ""))
    raise ValueError(f"machine must be one of {MACHINES}, got {machine!r}")


# --- matching -----------------------------------------------------------------------
def _co(lab: str) -> str:
    return lab[1:] if lab.startswith("~") else "~" + lab


def relabel(phi, lab: str) -> str:
    if lab == "tau":
        return lab
    name = lab.lstrip("~")
    new = dict(phi).get(name, name)
    return _co(new) if lab.startswith("~") else new


def _preimages(phi, lab: str) -> list[str]:
    if lab == "tau":
        return ["tau"]
    names = {lab.lstrip("~")} | {old for old, _ in phi}
    return [c for n in sorted(names) for c in (n, "~" + n) if relabel(phi, c) == lab]


def _bind(b: dict, name: str, value):
    if name in b:
        return b if b[name] == value else None
    return {**b, name: value}


def _match(pat, term, b):
    k = pat[0]
    if k == "?":
        kind = pat[1]
        ok = (kind == "mol" or (kind == "agent" and is_agent(term)) or
              (kind == "sol" and term[0] == "sol") or (kind == "fix" and term[0] == "fix"))
        if ok:
            b2 = _bind(b, pat[2], term)
            if b2 is not None:
                yield b2
        return
    if term[0] != k:
        return
    if k in ("0", "K", "fix"):
        if pat == term:
            yield b
    elif k == ".":
        for b1 in _match(pat[2], term[2], b):
            yield from _match_label(pat[1], term[1], b1)
    elif k in ("|", "[]", "(+)", "<|"):
        for b1 in _match(pat[1], term[1], b):
            yield from _match(pat[2], term[2], b1)
    elif k == "\\":
        for b1 in _match(pat[1], term[1], b):
            yield from _match_atom(pat[2], term[2], b1)
    elif k == "rel":
        for b1 in _match(pat[1], term[1], b):
            yield from _match_atom(pat[2], term[2], b1)
    elif k == "sol":
        if len(pat[1]) == len(term[1]) == 0:
            yield b
        elif len(pat[1]) == len(term[1]) == 1:
            yield from _match(pat[1][0], term[1][0], b)
    elif k == "<>":
        if pat[1] == term[1]:
            for b1 in _match(pat[2], term[2], b):
                yield from _match(pat[3], term[3], b1)


def _match_atom(pat, value, b):
    """Names (restriction) and relabellings: literal or a metavariable."""
    if isinstance(pat, tuple) and pat and pat[0] == "?":
        b2 = _bind(b, pat[2], value)
        if b2 is not None:
            yield b2
    elif pat == value:
        yield b


def _match_label(pat, lab: str, b):
    if isinstance(pat, str):
        if pat == lab:
            yield b
        return
    k = pat[0]
    if k == "?":
        if pat[1] == "name" and (lab == "tau" or lab.startswith("~")):
            return
        b2 = _bind(b, pat[2], lab)
        if b2 is not None:
            yield b2
    elif k == "~?":
        if lab.startswith("~"):
            b2 = _bind(b, pat[1], lab[1:])
            if b2 is not None:
                yield b2
    elif k == "app":
        if pat[1] not in b:
            raise ValueError(f"?{pat[1]}(...) is matched before ?{pat[1]} is bound; put ?{pat[1]} in the prefix body")
        for alpha in _preimages(b[pat[1]], lab):
            yield from _match_label(pat[2], alpha, b)


def _build_label(pat, b) -> str:
    if isinstance(pat, str):
        return pat
    if pat[0] == "?":
        return b[pat[2]]
    if pat[0] == "~?":
        return "~" + b[pat[1]]
    return relabel(b[pat[1]], _build_label(pat[2], b))


def _build(pat, b):
    k = pat[0]
    if k == "?":
        return b[pat[2]]
    if k == "unfold":
        return unfold(b[pat[1]])
    if k in ("0", "K", "fix"):
        return pat
    if k == ".":
        return (".", _build_label(pat[1], b), _build(pat[2], b))
    if k in ("|", "[]", "(+)", "<|"):
        return (k, _build(pat[1], b), _build(pat[2], b))
    if k in ("\\", "rel"):
        atom = pat[2]
        value = b[atom[2]] if isinstance(atom, tuple) and atom and atom[0] == "?" else atom
        return (k, _build(pat[1], b), value)
    if k == "sol":
        return _sol(_build(x, b) for x in pat[1])
    if k == "<>":
        return ("<>", pat[1], _build(pat[2], b), _build(pat[3], b))
    raise ValueError(f"cannot build {pat!r}")


def _guard_ok(guard, b) -> bool:
    op, lab, labels = guard
    inside = _build_label(lab, b) in {_build_label(x, b) for x in labels}
    return inside if op == "in" else not inside


def _substitute(t, env: dict):
    k = t[0]
    if k == "K":
        return env.get(t[1], t)
    if k == "fix":
        inner = {v: s for v, s in env.items() if v not in dict(t[2])}
        return ("fix", t[1], tuple((v, _substitute(p, inner)) for v, p in t[2]))
    if k == ".":
        return (".", t[1], _substitute(t[2], env))
    if k in ("|", "[]", "(+)"):
        return (k, _substitute(t[1], env), _substitute(t[2], env))
    if k in ("\\", "rel"):
        return (k, _substitute(t[1], env), t[2])
    return t


def unfold(fix):
    """fix_i(x = p) -> p_i[fix(x = p)/x]."""
    defs = dict(fix[2])
    env = {v: ("fix", v, fix[2]) for v in defs}
    return _substitute(defs[fix[1]], env)


def _match_all(d: Direction, mols):
    def rec(i, b):
        if i == len(d.lhs):
            if d.guard is None or _guard_ok(d.guard, b):
                yield b
            return
        for b1 in _match(d.lhs[i], mols[i], b):
            yield from rec(i + 1, b1)
    yield from rec(0, {})


# --- the machine ----------------------------------------------------------------------
class Machine:
    """One-step transitions of a CHAM under its rules and the four general laws.

    full=False is the exploration policy of the network closure (see module
    docstring); full=True is the unrestricted relation.
    """

    def __init__(self, rules: list[Rule], airlock: bool = True, full: bool = False):
        self.rules, self.airlock, self.full = rules, airlock, full
        self.dirs = []
        for r in rules:
            if r.kind == "reversible":
                self.dirs.append(Direction("heating", r.lhs, r.rhs, r.guard, r.name, True))
                if full or len(r.rhs) == 1:
                    self.dirs.append(Direction("cooling", r.rhs, r.lhs, r.guard, r.name, True))
            else:
                cls = r.kind
                self.dirs.append(Direction(cls, r.lhs, r.rhs, r.guard, r.name, False))
        self._cache: dict = {}

    def _wraps(self, d: Direction, mols, products) -> bool:
        """Closure policy: skip a heating that puts a membrane directly around a membrane."""
        if self.full or d.cls != "heating":
            return False
        return sum(map(_wrapped_membranes, products)) > sum(map(_wrapped_membranes, mols))

    @staticmethod
    def _label(cls, rule, reversible, depth=0):
        return {"class": cls, "rule": rule, "reversible": reversible, "depth": depth}

    def react(self, mols) -> list[tuple[tuple, dict]]:
        """Rule, airlock and membrane transitions whose reactants are exactly `mols`."""
        out = []
        for d in self.dirs:
            if len(d.lhs) == len(mols):
                for b in _match_all(d, mols):
                    products = tuple(_build(x, b) for x in d.rhs)
                    if self._wraps(d, mols, products):
                        continue
                    out.append((products, self._label(d.cls, d.rule, d.reversible)))
        if len(mols) == 1:
            m = mols[0]
            if self.airlock and m[0] == "<|":
                out.append(((m[1],) + m[2][1], self._label("cooling", "airlock", True)))
            for m2, lab in self.molecule_steps(m):
                out.append(((m2,), lab))
        return out

    def molecule_steps(self, m) -> list[tuple[tuple, dict]]:
        """Membrane law: transitions of the outermost solutions inside m."""
        out = []
        for path, s in _outer_solutions(m):
            for items, lab in self.solution_steps(s[1], top=False):
                out.append((_replace(m, path, ("sol", items)), {**lab, "depth": lab["depth"] + 1}))
        return out

    def solution_steps(self, items: tuple, top: bool = True) -> list[tuple[tuple, dict]]:
        key = (items, top)
        if key in self._cache:
            return self._cache[key]
        n = len(items)
        found: dict[tuple, dict] = {}

        def add(new, lab):
            new = _sol(new)[1]
            if new != items:
                found.setdefault(new, lab)

        for d in self.dirs:
            k = len(d.lhs)
            for idx in permutations(range(n), k):
                mols = [items[i] for i in idx]
                for b in _match_all(d, mols):
                    products = [_build(x, b) for x in d.rhs]
                    if self._wraps(d, mols, products):
                        continue
                    rest = [items[i] for i in range(n) if i not in idx]
                    add(rest + products, self._label(d.cls, d.rule, d.reversible))
        if self.airlock:
            for i, m in enumerate(items):
                others = items[:i] + items[i + 1:]
                if m[0] == "<|":
                    add(list(others) + [m[1], *m[2][1]], self._label("cooling", "airlock", True))
                if self.full:
                    for r in range(len(others) + 1):
                        for chosen in combinations(range(len(others)), r):
                            kept = [others[j] for j in range(len(others)) if j not in chosen]
                            add(kept + [("<|", m, _sol(others[j] for j in chosen))],
                                self._label("heating", "airlock", True))
                elif not top and n >= 2:
                    add([("<|", m, _sol(others))], self._label("heating", "airlock", True))
        for i, m in enumerate(items):
            for m2, lab in self.molecule_steps(m):
                add(list(items[:i] + items[i + 1:]) + [m2], lab)
        result = list(found.items())
        self._cache[key] = result
        return result


def _wrapped_membranes(t) -> int:
    """Number of solutions whose only member is itself a solution ({|{|...|}|})."""
    if not isinstance(t, tuple) or not t or not isinstance(t[0], str):
        return 0
    if t[0] == "sol":
        own = int(len(t[1]) == 1 and t[1][0][0] == "sol")
        return own + sum(_wrapped_membranes(x) for x in t[1])
    return sum(_wrapped_membranes(t[i]) for i in _CHILDREN.get(t[0], ()))


def _outer_solutions(t, path=()):
    if t[0] == "sol":
        yield path, t
        return
    for i in _CHILDREN.get(t[0], ()):
        yield from _outer_solutions(t[i], path + (i,))


def _replace(t, path, new):
    if not path:
        return new
    i = path[0]
    return t[:i] + (_replace(t[i], path[1:], new),) + t[i + 1:]


def transitions(solution, rules: list[Rule], airlock: bool = True, full: bool = True):
    """One-step transitions S -> S' of a whole solution: [(solution text, label)]."""
    items = parse_solution(solution) if isinstance(solution, str) else _sol(solution)[1]
    machine = Machine(rules, airlock, full)
    return [(solution_text(new), lab) for new, lab in machine.solution_steps(items, top=True)]


def closure(rules: list[Rule], seed, airlock: bool = True, max_species: int = 500):
    """Reaction closure of the top solution; returns (species, [(lhs, rhs, label)], status)."""
    machine = Machine(rules, airlock, full=False)
    arities = sorted({1} | {len(d.lhs) for d in machine.dirs})

    def react(*mols):
        out = [rhs for rhs, _ in machine.react(mols)]
        return out or None

    species, pairs, status = expand(react, seed, arity=arities, max_species=max_species,
                                    alternatives=True)
    found, seen = [], set()

    def add(lhs, rhs, lab):
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        if key not in seen:
            seen.add(key)
            found.append((tuple(lhs), tuple(rhs), lab))

    for lhs, rhs in pairs:
        target = Counter(rhs)
        lab = next(l for r, l in machine.react(lhs) if Counter(r) == target)
        add(lhs, rhs, lab)
        if lab["reversible"]:
            add(rhs, lhs, {**lab, "class": _FLIP[lab["class"]]})
    return species, found, status


# --- compartments ---------------------------------------------------------------------
def _context(parent, i) -> str:
    k = parent[0]
    if k == "\\":
        return f"restriction \\{parent[2]}"
    if k == "rel":
        return "relabelling [" + ",".join(f"{new}/{old}" for old, new in parent[2]) + "]"
    if k == "<|":
        return "airlock"
    if k == "<>":
        side = "left" if i == 2 else "right"
        return f"{parent[1] + ':' if parent[1] else ''}pair {side}"
    if k == ".":
        return f"prefix {parent[1]}"
    return "molecule"


def membranes(m, depth: int = 1) -> list[dict]:
    """The solutions nested in a molecule, in preorder, with their enclosing operator."""
    out = []

    def walk(t, parent, i, d):
        if t[0] == "sol":
            ctx = "solution" if parent is None else _context(parent, i)
            out.append({"depth": d, "context": ctx, "molecules": [text(x) for x in t[1]]})
            for x in t[1]:
                walk(x, None, None, d + 1)
            return
        for j in _CHILDREN.get(t[0], ()):
            walk(t[j], t, j, d)

    walk(m, None, None, depth)
    return out


# --- generator ------------------------------------------------------------------------
def initial_solution(program: str, solution: str) -> tuple:
    if program == "custom":
        if not solution.strip():
            raise ValueError("program=custom needs solution, e.g. 'a.0|(~a.P|Q)\\b'")
        return parse_solution(solution)
    if solution.strip():
        raise ValueError("solution is only used with program=custom")
    if program not in PROGRAM_SOLUTIONS:
        raise ValueError(f"program must be one of {PROGRAMS}, got {program!r}")
    return parse_solution(PROGRAM_SOLUTIONS[program])


def generate(p, rng):
    rules = machine_rules(p.machine, p.cleanup, p.rules)
    items = initial_solution(p.program, p.solution)
    species, found, status = closure(rules, list(dict.fromkeys(items)), p.airlock, p.max_species)

    compartments = {}
    for m in species:
        nested = membranes(m)
        if nested:
            compartments[text(m)] = nested
    extras = {
        "machine": p.machine,
        "rules": [r.text for r in rules],
        "airlock_law": bool(p.airlock),
        "reaction_rules": [lab for _, _, lab in found],
        "compartments": compartments,
    }
    return Network(
        species=[Species(text(m), pretty(m)) for m in species],
        reactions=[Reaction.of([text(x) for x in lhs], [text(x) for x in rhs]) for lhs, rhs, _ in found],
        status=status,
        initial_state={s: float(c) for s, c in Counter(text(m) for m in items).items()},
        extras=extras,
    )
