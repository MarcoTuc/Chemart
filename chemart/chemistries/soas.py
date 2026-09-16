"""Self-Organising Assembly Systems (SOAS). Catalog id: soas.

Book 20.1 (book [296], [297]); R. Frei, G. Di Marzo Serugendo & T. F. Serbanuta, "Ambient
intelligence in self-organising assembly systems using the chemical reaction model", J. Ambient
Intell. Humanized Comput. 1:163-184 (2010); R. Frei, T. F. Serbanuta & G. Di Marzo Serugendo,
"Self-organising assembly systems formally specified in Maude", J. Ambient Intell. Humanized
Comput. 5:491-510 (online 2012); and the authors' Maude specification soas-maude 1.0
(code.google.com/p/soas-maude), whose rules this module reimplements.

An assembly system designs itself for a product order. The chemical solution holds the
manufacturing resource agents (MRAs: robots, axes, grippers, feeders, positioning devices, the
human operator), the order (a generic assembly plan, GAP, made of tasks) and the parts. The
reactions are the rule types of the papers:

- initiate (type 4, task-coalition matching): an MRA offering a skill that a task requires starts
  a coalition molecule for that task;
- join (types 1-3: interface compatibility, composition patterns, composite skills): an MRA joins
  a coalition when it offers a skill the coalition still requires, a composition pattern links
  its skills with the coalition's, and both have a free interface of the same type and opposite
  sign (the two interfaces become occupied);
- assign (coalition selection): a complete coalition (nothing required any more) is assigned to
  its open task in a partial assembly line when it shares no MRA with the line and every task
  still open keeps a complete coalition sharing no MRA with the new line. Coalitions of the
  generic human operator are assigned first, all at once, and are never exclusive
  (assign-generic);
- layout (type 5): each complete line is laid out as a serpentine of conveyors on the shop floor
  (extras.analysis; it is a function of the line, not a reaction).

Before any reaction the order and the grippers are completed from the parts (2012 extension of
rule 4): a task whose first object part needs a gripper subtype G also requires
grip(gripper-type=G), with the part's grip width as range for 2-finger grippers; a task on one
part of type P starting at the feeder also requires feed(subtype=feeds(P)); a gripper of subtype
G offers grip(gripper-type=G), absorbing its open-close range. Composite skills are decomposed
(type 3): pick&place => move, move(direction=vertical).

Text notation (one item per line, '#' starts a comment):

    MRA      name type[/subtype] | skills | required skills | interfaces
             g1 gripper/2finger | open-close(range=70.0) | | circular+
             a4 axis | move(subtype=linear,direction=horizontal,range=250.0), hold | base-hold | square+, square-
    skill    type or type(key=value,...); keys subtype, direction, gripper-type, part-type (words),
             range, move-range, open-range (numbers), count, accuracy (integers)
    interface  type then + (held, passive) or - (holding, active): circular-, triangular+
    part     name type | gripper subtype | grip positions     p1 body-case | 2finger | (25,0,10) (25,60,10)
    GAP      first line "name number" (products to assemble), then one task per line:
             t(N) | type | object parts joined by + | start point | operations
             t(2) | pick&place | p1 | feeder | pick&place
             operations are separated by ';', the skills of one operation by ','
    pattern  skills & skills                                   grip, move & store(subtype=grippers)

An offered skill matches a required one when both have the same type and every item of the
required skill is an item of the offered one, except range, where the offered range must be at
least the required range. A pattern side matches a skill set when some skill of the set matches
some skill of the side. Nothing is evaluated as Python.

Species ids: an MRA is its name (r1); the order is the GAP name (gap2), which is also the empty
assembly line; a coalition is gap2->t(2)[f1,g1,r1]; a (partial) assembly line is
gap2{t(1)=unassigned-human;t(2)=f1+g1+r1;...}. The network is the closure of the initial solution
(the order and the MRAs); reactions have no rates.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

SYSTEMS = ("tape-roller-gap2", "tape-roller-gap1", "custom")

#: The MRA standing for any number of human operators (soas-maude generic-human).
GENERIC_HUMAN = "unassigned-human"

#: Layout constants of soas-maude 1.0 (syntax/mra-types.maude).
CONVEYOR_LENGTH, CONVEYOR_WIDTH, MRA_WIDTH = 300, 50, 50

#: soas-maude 1.0 mras/*.maude, the repository listed in mras/mras.maude (JAmI 2012 sec. 5.1).
TAPE_ROLLER_MRAS = """\
r3 robot | move(subtype=rotational,direction=vertical,range=220.0), move(subtype=rotational,direction=horizontal,range=280.0), move(subtype=rotational,direction=horizontal,range=60.0) | | circular-, triangular-, straight-, diamond-
r1 robot | move(subtype=linear,direction=horizontal,range=225.0), move(subtype=linear,direction=horizontal,range=220.0), move(subtype=linear,direction=vertical,range=125.0) | | circular-, triangular-, straight-, diamond-
r2 robot | move(subtype=linear,direction=vertical,range=300.0), move(subtype=rotational,direction=horizontal,range=180.0), move(subtype=rotational,direction=horizontal,range=90.0) | | circular-, triangular-, straight-, diamond-
a1 axis | move(subtype=linear,direction=horizontal,range=150.0), hold, base-hold | | square-, triangular-, straight-
a4 axis | move(subtype=linear,direction=horizontal,range=250.0), move(subtype=rotational,direction=vertical,range=180.0), hold | base-hold | square+, square-
a7 axis | move(subtype=linear,direction=horizontal,range=150.0) | hold | square+, circular-
f1 feeder | feed(subtype=feeds(body-case)) | | triangular+
f2 feeder | feed(subtype=feeds(tape-roll)) | | triangular+
f3 feeder | feed(subtype=feeds(body-case)) | | triangular+
f4 feeder | feed(subtype=feeds(screw)) | | triangular+
g1 gripper/2finger | open-close(range=70.0) | | circular+
g2 gripper/2finger | open-close(range=55.0) | | circular+
g3 gripper/vacuum | open-close | | circular+
g4 gripper/screw-driver | open-close | | circular+
pd1 positioning-device | position-carrier(accuracy=1) | | straight+, oval-
pd2 positioning-device | position-carrier(accuracy=2) | | straight+, oval-
pd3 positioning-device | position-carrier(accuracy=3) | | straight+, oval-
pd4 positioning-device | position-carrier(accuracy=3) | | straight+, oval-
unassigned-human human | load, unload | |
"""

#: soas-maude 1.0 parts/tape-roller.maude (only the items the rules read).
TAPE_ROLLER_PARTS = """\
p1 body-case | 2finger | (25,0,10) (25,60,10)
p2 tape-roll | 2finger | (20,0,10) (20,40,10)
p3 body-case | vacuum | (25,0,9) (25,50,9)
p4 screw | screw-driver | (2,2,15)
wpca carrier | |
"""

#: soas-maude 1.0 assembly/gap2.maude: the four-part tape roller of JAmI 2012 (Fig. 2).
GAP2 = """\
gap2 10
t(1) | other | wpca | (-10,50,-20) | load
t(2) | pick&place | p1 | feeder | pick&place
t(3) | pick&place | p2 | feeder | pick&place
t(4) | pick&place | p3 | feeder | pick&place
t(5) | pick&place | p4 | feeder | pick&place
t(6) | other | wpca+p1+p2+p3+p4 | | unload
"""

#: soas-maude 1.0 assembly/gap1.maude: the two-part tape roller of JAmI 2010 (Figs. 4, 18, 20),
#: whose pick&place tasks start at named feeders f1 and f2.
GAP1 = """\
gap1 10
t(1) | other | wpca | (-10,50,-20) | load
t(2) | pick&place | p1 | f1 | pick&place
t(3) | pick&place | p2 | f2 | pick&place
t(4) | other | wpca+p1+p2 | wpca | unload
"""

#: soas-maude 1.0 semantics/patterns.maude (rule type 2).
PATTERNS = """\
move & move
move & grip
grip, move & store(subtype=grippers)
transport(subtype=linear) & transport(subtype=linear)
move & position-carrier
feed & move
"""

NAMED = {
    "tape-roller-gap2": (TAPE_ROLLER_MRAS, GAP2, TAPE_ROLLER_PARTS),
    "tape-roller-gap1": (TAPE_ROLLER_MRAS, GAP1, TAPE_ROLLER_PARTS),
}

#: The reaction schemata, as reported in extras.rules.
RULES = (
    "initiate: M, O -> M, O, [T|M]  if MRA M offers a skill that task T of order O requires",
    "join: M, C -> M, C, C+M  if M is not in coalition C, M offers a skill that C still requires, "
    "a composition pattern links the skills of M and C, and M and C have free interfaces of one "
    "type and opposite signs (both become occupied)",
    "assign-generic: O, [T1|unassigned-human], ..., [Tk|unassigned-human] -> O{T1..Tk}  "
    "(every task completed by the generic human operator, at once)",
    "assign: L, C -> L+C  if C is complete, its task is open in line L, C shares no MRA with L, "
    "and every task still open keeps a complete coalition sharing no MRA with L+C",
)
COMPOSITE_SKILLS = (
    "pick&place => move(direction, range=move-range), move(direction=vertical)",
    "load-unload-gripper => pick&place, store(subtype=grippers, count)",
    "index-carrier-for-task => move(direction, range=move-range), position-carrier",
)
PART_MATCHING = (
    "a task whose first object part needs gripper subtype G also requires grip(gripper-type=G); "
    "for G = 2finger its range is the distance between the part's first two grip positions",
    "a task on a single part of type P whose start point is the feeder also requires "
    "feed(subtype=feeds(P))",
    "a gripper of subtype G offers grip(gripper-type=G); its open-close(range=F) becomes "
    "grip(gripper-type=G,range=F)",
)

KEYS = ("subtype", "direction", "gripper-type", "part-type", "range", "move-range", "open-range",
        "count", "accuracy")
FLOAT_KEYS = {"range", "move-range", "open-range"}
INT_KEYS = {"count", "accuracy"}

_WORD = r"[A-Za-z0-9][A-Za-z0-9&_.-]*"
_SKILL = re.compile(rf"^({_WORD})(?:\((.*)\))?$")
_VALUE = re.compile(rf"^{_WORD}(?:\({_WORD}\))?$")
_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
_INTERFACE = re.compile(r"^([a-z][a-z-]*)([+-])$")
_TASK = re.compile(r"^t\((\d+)\)$")
_POINT = re.compile(r"\((-?\d+),(-?\d+),(-?\d+)\)")


# --- molecules --------------------------------------------------------------------------
@dataclass(frozen=True)
class Skill:
    type: str
    items: tuple = ()        # ((key, value), ...) in KEYS order

    def text(self) -> str:
        if not self.items:
            return self.type
        return f"{self.type}({','.join(f'{k}={v}' for k, v in self.items)})"


@dataclass(frozen=True)
class Interface:
    type: str
    sign: str

    def text(self) -> str:
        return self.type + self.sign


@dataclass(frozen=True)
class Mra:
    name: str
    type: str
    subtype: str | None
    skills: tuple
    requires: tuple
    interfaces: tuple

    def text(self) -> str:
        kind = self.type + (f"/{self.subtype}" if self.subtype else "")
        return " | ".join([f"{self.name} {kind}", _skills_text(self.skills),
                           _skills_text(self.requires),
                           ", ".join(i.text() for i in self.interfaces)]).rstrip()


@dataclass(frozen=True)
class Part:
    name: str
    type: str
    gripper: str | None
    grips: tuple


@dataclass(frozen=True)
class Task:
    number: int
    type: str
    object: tuple
    start: str
    operations: tuple        # (tuple of skills, ...)

    @property
    def name(self) -> str:
        return f"t({self.number})"

    def text(self) -> str:
        ops = "; ".join(_skills_text(op) for op in self.operations)
        return f"{self.name} | {self.type} | {'+'.join(self.object)} | {self.start} | {ops}"


@dataclass(frozen=True)
class Order:
    name: str
    number: int
    tasks: tuple

    def text(self) -> str:
        return "\n".join([f"{self.name} {self.number}", *(t.text() for t in self.tasks)])


@dataclass(frozen=True)
class Coalition:
    order: str
    task: int
    members: tuple           # sorted MRA names
    provided: tuple          # skills offered by the members
    required: tuple          # operations still required (tuples of skills)
    interfaces: tuple        # free interfaces

    @property
    def id(self) -> str:
        return f"{self.order}->t({self.task})[{','.join(self.members)}]"

    @property
    def complete(self) -> bool:
        return not self.required

    def text(self) -> str:
        required = " | ".join(_skills_text(op) for op in self.required) or "none"
        return (f"task={self.order}->t({self.task}); coalition={', '.join(self.members)}; "
                f"provided-skills={_skills_text(self.provided)}; required-ops={required}; "
                f"open-interfaces={', '.join(i.text() for i in self.interfaces)}")


def _skills_text(skills) -> str:
    return ", ".join(s.text() for s in skills)


# --- parsing --------------------------------------------------------------------------------
def _split(text: str, sep: str = ",") -> list[str]:
    out, depth, cur = [], 0, []
    for ch in text:
        depth += (ch == "(") - (ch == ")")
        if ch == sep and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [s.strip() for s in out if s.strip()]


def _lines(text: str) -> list[str]:
    return [line for raw in text.splitlines() if (line := raw.split("#", 1)[0].strip())]


def parse_skill(src: str) -> Skill:
    """'move(subtype=linear,range=150.0)' -> Skill('move', (('subtype', 'linear'), ('range', 150.0)))."""
    m = _SKILL.match(src.strip())
    if not m:
        raise ValueError(f"cannot read skill {src!r}: expected type or type(key=value,...)")
    items = {}
    for part in _split(m.group(2) or ""):
        key, eq, value = (s.strip() for s in part.partition("="))
        if not eq or key not in KEYS:
            raise ValueError(f"skill {src!r}: items are key=value with key in {KEYS}, got {part!r}")
        if key in items:
            raise ValueError(f"skill {src!r}: repeated item {key!r}")
        try:
            if key in FLOAT_KEYS:
                items[key] = float(value)
            elif key in INT_KEYS:
                items[key] = int(value)
            elif _VALUE.match(value):
                items[key] = value
            else:
                raise ValueError
        except ValueError:
            raise ValueError(f"skill {src!r}: bad value {value!r} for {key!r}") from None
    return Skill(m.group(1), tuple((k, items[k]) for k in KEYS if k in items))


def _skill_list(src: str) -> tuple:
    return tuple(parse_skill(s) for s in _split(src))


def parse_mras(text: str) -> list[Mra]:
    out, names = [], set()
    for line in _lines(text):
        fields = [f.strip() for f in line.split("|")]
        if len(fields) != 4:
            raise ValueError(f"MRA {line!r} needs 4 fields: name type[/subtype] | skills | "
                             "required skills | interfaces")
        head = fields[0].split()
        if len(head) != 2 or not _NAME.match(head[0]):
            raise ValueError(f"MRA {line!r} must start with 'name type[/subtype]'")
        name, (kind, _, subtype) = head[0], head[1].partition("/")
        if name in names:
            raise ValueError(f"MRA names must be unique, repeated: {name!r}")
        names.add(name)
        interfaces = []
        for src in _split(fields[3]):
            m = _INTERFACE.match(src)
            if not m:
                raise ValueError(f"MRA {name}: interface {src!r} must be a type followed by + or -")
            interfaces.append(Interface(m.group(1), m.group(2)))
        mra = Mra(name, kind, subtype or None, _skill_list(fields[1]), _skill_list(fields[2]),
                  tuple(interfaces))
        out.append(Mra(mra.name, mra.type, mra.subtype, _gripper_skills(mra), mra.requires,
                       mra.interfaces))
    return out


def parse_parts(text: str) -> list[Part]:
    out, names = [], set()
    for line in _lines(text):
        fields = [f.strip() for f in line.split("|")]
        head = fields[0].split()
        if len(fields) != 3 or len(head) != 2 or not _NAME.match(head[0]):
            raise ValueError(f"part {line!r} must be 'name type | gripper subtype | grip positions'")
        if head[0] in names:
            raise ValueError(f"part names must be unique, repeated: {head[0]!r}")
        names.add(head[0])
        compact = fields[2].replace(" ", "")
        points = tuple(tuple(int(v) for v in m) for m in _POINT.findall(compact))
        if _POINT.sub("", compact):
            raise ValueError(f"part {head[0]}: grip positions must be (x,y,z) integer triples, "
                             f"got {fields[2]!r}")
        out.append(Part(head[0], head[1], fields[1] or None, points))
    return out


def parse_gap(text: str) -> Order:
    lines = _lines(text)
    if not lines:
        raise ValueError("the GAP needs a first line 'name number' and one task per line")
    head = lines[0].split()
    if len(head) != 2 or not _NAME.match(head[0]) or not head[1].isdigit():
        raise ValueError(f"the GAP must start with 'name number', got {lines[0]!r}")
    tasks, numbers = [], set()
    for line in lines[1:]:
        fields = [f.strip() for f in line.split("|")]
        if len(fields) != 5:
            raise ValueError(f"task {line!r} needs 5 fields: t(N) | type | object | start | operations")
        m = _TASK.match(fields[0])
        if not m or int(m.group(1)) < 1:
            raise ValueError(f"task names are t(N) with N >= 1, got {fields[0]!r}")
        number = int(m.group(1))
        if number in numbers:
            raise ValueError(f"task t({number}) is repeated")
        numbers.add(number)
        parts = tuple(p.strip() for p in fields[2].split("+") if p.strip())
        if any(not _NAME.match(p) for p in parts):
            raise ValueError(f"task t({number}): the object is part names joined by '+'")
        operations = tuple(_skill_list(op) for op in fields[4].split(";") if op.strip())
        if not operations:
            raise ValueError(f"task t({number}) needs at least one operation")
        tasks.append(Task(number, fields[1], parts, fields[3], operations))
    if not tasks:
        raise ValueError("the GAP has no tasks")
    return Order(head[0], int(head[1]), tuple(sorted(tasks, key=lambda t: t.number)))


def parse_patterns(text: str) -> list[tuple[tuple, tuple]]:
    out = []
    for line in _lines(text):
        sides = line.split("&")
        if len(sides) != 2 or not sides[0].strip() or not sides[1].strip():
            raise ValueError(f"composition pattern {line!r} must be 'skills & skills'")
        out.append((_skill_list(sides[0]), _skill_list(sides[1])))
    return out


# --- skills ---------------------------------------------------------------------------------
def skill_match(offered: Skill, required: Skill) -> bool:
    """soas-maude MATCH-SKILLS: same type, required items present, offered range >= required."""
    if offered.type != required.type:
        return False
    for key, value in required.items:
        if key == "range":
            if not any(k == "range" and v >= value for k, v in offered.items):
                return False
        elif (key, value) not in offered.items:
            return False
    return True


def matches(offered, required) -> bool:
    return any(skill_match(o, r) for o in offered for r in required)


def matches_ops(offered, operations) -> bool:
    return any(matches(offered, op) for op in operations)


def subtract(operations, offered) -> tuple:
    """Remove from every operation the skills some offered skill matches; drop empty ones."""
    out = []
    for op in operations:
        rest = tuple(r for r in op if not any(skill_match(o, r) for o in offered))
        if rest:
            out.append(rest)
    return tuple(out)


def _ordered(items) -> tuple:
    d = dict(items)
    return tuple((k, d[k]) for k in KEYS if k in d)


def _move_items(items) -> tuple:
    return _ordered(("range" if k == "move-range" else k, v) for k, v in items
                    if k in ("move-range", "direction"))


def decompose(skills) -> tuple:
    """Composite skills into simple ones (soas-maude DECOMPOSE-SKILLS)."""
    out = []
    for s in skills:
        if s.type == "pick&place":
            out += [Skill("move", _move_items(s.items)), Skill("move", (("direction", "vertical"),))]
        elif s.type == "load-unload-gripper":
            out += decompose([Skill("pick&place", s.items)])
            out.append(Skill("store", _ordered([("subtype", "grippers")]
                                               + [(k, v) for k, v in s.items if k == "count"])))
        elif s.type == "index-carrier-for-task":
            out += [Skill("move", _move_items(s.items)), Skill("position-carrier")]
        else:
            out.append(s)
    return tuple(out)


def _gripper_skills(mra: Mra) -> tuple:
    """PART-MATCH-GRIPPER on MRAs: a gripper of subtype G offers grip(gripper-type=G)."""
    if mra.type != "gripper" or not mra.subtype:
        return mra.skills
    skills = list(mra.skills)
    grip = Skill("grip", (("gripper-type", mra.subtype),))
    probe = Skill("grip", (("gripper-type", mra.subtype), ("range", math.inf)))
    if not any(skill_match(probe, s) for s in skills):
        skills.append(grip)
    oc = next((s for s in skills if s.type == "open-close" and len(s.items) == 1
               and s.items[0][0] == "range"), None)
    if oc is not None and grip in skills:
        skills.remove(oc)
        skills.remove(grip)
        skills.append(Skill("grip", (("gripper-type", mra.subtype), ("range", oc.items[0][1]))))
    return tuple(skills)


def complete_order(order: Order, parts: list[Part]) -> Order:
    """PART-MATCH-GRIPPER and PART-MATCH-FEEDER on the tasks of the order."""
    by_name = {p.name: p for p in parts}
    tasks = []
    for task in order.tasks:
        ops = list(task.operations)
        first = by_name.get(task.object[0]) if task.object else None
        if first is not None and first.gripper:
            probe = Skill("grip", (("gripper-type", first.gripper), ("range", math.inf)))
            if not matches_ops((probe,), ops):
                ops.append((Skill("grip", (("gripper-type", first.gripper),)),))
            if first.gripper == "2finger" and len(first.grips) >= 2:
                bare = (Skill("grip", (("gripper-type", "2finger"),)),)
                width = math.dist(first.grips[0], first.grips[1])
                wide = (Skill("grip", (("gripper-type", "2finger"), ("range", width))),)
                ops = [wide if op == bare else op for op in ops]
        if len(task.object) == 1 and first is not None and task.start == "feeder":
            feed = Skill("feed", (("subtype", f"feeds({first.type})"),))
            if not matches_ops((feed,), ops):
                ops.append((feed,))
        tasks.append(Task(task.number, task.type, task.object, task.start, tuple(ops)))
    return Order(order.name, order.number, tuple(tasks))


# --- reactions ------------------------------------------------------------------------------
def _canon_skills(skills) -> tuple:
    return tuple(sorted(skills, key=Skill.text))


def _canon_ops(ops) -> tuple:
    return tuple(sorted((_canon_skills(op) for op in ops), key=_skills_text))


def initiate(mra: Mra, order: Order, task: Task) -> Coalition | None:
    """Rule type 4: MRA M starts a coalition for task T if it offers a required skill."""
    offered, required = decompose(mra.skills), tuple(decompose(op) for op in task.operations)
    if not matches_ops(offered, required):
        return None
    ops = subtract(required, offered) + ((mra.requires,) if mra.requires else ())
    return Coalition(order.name, task.number, (mra.name,), _canon_skills(mra.skills),
                     _canon_ops(ops), tuple(sorted(mra.interfaces, key=Interface.text)))


def join(mra: Mra, coalition: Coalition, patterns) -> list[tuple[Coalition, str]]:
    """Rule types 1-3: every way MRA M can join coalition C, with the interface pair used."""
    if mra.name in coalition.members:
        return []
    offered, provided = decompose(mra.skills), decompose(coalition.provided)
    if not matches_ops(offered, coalition.required):
        return []
    if not any((matches(offered, a) and matches(provided, b))
               or (matches(offered, b) and matches(provided, a)) for a, b in patterns):
        return []
    requires = tuple(s for s in mra.requires if s not in provided)
    ops = subtract(coalition.required, offered) + ((requires,) if requires else ())
    out, seen = [], set()
    for i, own in enumerate(mra.interfaces):
        for j, other in enumerate(coalition.interfaces):
            if own.type != other.type or own.sign == other.sign or (own, other) in seen:
                continue
            seen.add((own, other))
            free = mra.interfaces[:i] + mra.interfaces[i + 1:] + coalition.interfaces[:j] \
                + coalition.interfaces[j + 1:]
            new = Coalition(coalition.order, coalition.task,
                            tuple(sorted(coalition.members + (mra.name,))),
                            _canon_skills(coalition.provided + mra.skills), _canon_ops(ops),
                            tuple(sorted(free, key=Interface.text)))
            out.append((new, f"{own.text()}/{other.text()}"))
    return out


def coalition_closure(order: Order, mras: list[Mra], patterns, max_species: int):
    """Initiate and join reactions from the order and the MRAs (chemart.expand, arity 2).

    A coalition is identified by its task and members, as in soas-maude, where a coalition is
    only created if none with the same members exists for the task; the first derivation found
    is kept and later derivations of the same members with a different state are counted in
    `conflicts`.
    """
    registry: dict[tuple, Coalition] = {}
    conflicts: set[tuple] = set()
    labels: dict[tuple, dict] = {}

    def register(c: Coalition) -> Coalition:
        key = (c.task, c.members)
        known = registry.setdefault(key, c)
        if known != c:
            conflicts.add(key)
        return known

    def react(a, b):
        mra, other = (a, b) if isinstance(a, Mra) else (b, a)
        if not isinstance(mra, Mra) or isinstance(other, Mra):
            return None
        out = []
        if isinstance(other, Order):
            for task in other.tasks:
                c = initiate(mra, other, task)
                if c is not None:
                    c = register(c)
                    labels.setdefault((mra.name, other.name, c.id), {"rule": "initiate", "task": task.name})
                    out.append((mra, other, c))
        elif isinstance(other, Coalition):
            for c, pair in join(mra, other, patterns):
                c = register(c)
                labels.setdefault((mra.name, other.id, c.id),
                                  {"rule": "join", "task": f"t({c.task})", "interfaces": pair})
                out.append((mra, other, c))
        return out or None

    species, pairs, status = expand(react, [order, *mras], arity=2, max_species=max_species,
                                    ordered=False, alternatives=True)
    reactions = []
    for lhs, rhs in pairs:
        mra = next(x for x in lhs if isinstance(x, Mra))
        other = next(x for x in lhs if not isinstance(x, Mra))
        other_id = other.name if isinstance(other, Order) else other.id
        reactions.append((lhs, rhs, labels[(mra.name, other_id, rhs[2].id)]))
    return species, reactions, status, len(conflicts)


def line_closure(order: Order, coalitions: list[Coalition], mra_names: list[str], max_species: int):
    """Assign reactions from the order (the empty line), arity 1 over lines.

    A line is a tuple of (task, index into the complete coalitions). Returns the lines in
    discovery order (the first is the order itself), the transitions (line, new line) and the
    status.
    """
    complete = [c for c in coalitions if c.complete]
    bit = {name: 1 << i for i, name in enumerate(mra_names)}
    masks = [sum(bit[m] for m in c.members) for c in complete]
    generic = [c.members == (GENERIC_HUMAN,) for c in complete]
    tasks = [t.number for t in order.tasks]
    by_task = {t: [i for i, c in enumerate(complete) if c.task == t] for t in tasks}
    human = tuple((t, next(i for i in by_task[t] if generic[i]))
                  for t in tasks if any(generic[i] for i in by_task[t]))

    def react(x):
        if isinstance(x, Order):
            if human:
                return [(("line", human),)]
            assign = ()
        else:
            assign = x[1]
        done = {t for t, _ in assign}
        open_tasks = [t for t in tasks if t not in done]
        used = 0
        for _, i in assign:
            if not generic[i]:
                used |= masks[i]
        out = []
        for t in open_tasks:
            rest = [s for s in open_tasks if s != t]
            for i in by_task[t]:
                if masks[i] & used:
                    continue
                now = used | masks[i]
                if all(any(not masks[j] & now for j in by_task[s]) for s in rest):
                    out.append((("line", tuple(sorted(assign + ((t, i),)))),))
        return out or None

    lines, pairs, status = expand(react, [order], arity=1, max_species=max_species,
                                  ordered=True, alternatives=True)
    return complete, lines, [(lhs[0], rhs[0]) for lhs, rhs in pairs], status, bool(human)


# --- layout (rule type 5) -------------------------------------------------------------------
def _task_id(n: Fraction) -> str:
    return f"t({n.numerator})" if n.denominator == 1 else f"t({n.numerator}/{n.denominator})"


def layout(assignment: dict[int, str], floor=(1000, 10000)) -> list[dict] | None:
    """soas-maude 1.0 semantics/layout.maude: a serpentine line starting at (0,0,0) eastwards.

    `assignment` maps task numbers to coalitions (members joined by '+'). A station is placed
    with a linear conveyor in front of it while the floor allows, otherwise a corner, a
    northward conveyor and a corner turn the line. The loading operator of t(1) is not placed.
    Returns None when the floor is too small (the specification then has no solution).
    """
    xe, ye = floor
    length, width = CONVEYOR_LENGTH, CONVEYOR_WIDTH
    order = sorted(assignment)
    if order and order[0] == 1 and assignment[1] == GENERIC_HUMAN:
        order = order[1:]
    stations: list[dict] = []

    def put(task, kind, subtype, x, y, angle, coalition=None):
        station = {"task": _task_id(Fraction(task)), "kind": kind, "subtype": subtype,
                   "position": [x, y, 0], "angle": [angle, 0, 0]}
        if coalition is not None:
            station["coalition"] = coalition
        stations.append(station)

    x = y = 0
    east = True
    current = Fraction(order[0]) if order else Fraction(0)
    i = 0
    while i < len(order):
        t = order[i]
        if east and x + length + width <= xe:
            put(t, "conveyor", "linear", x, y, 0)
            put(t, "coalition", None, x + length - 2 * MRA_WIDTH, y + width, 0, assignment[t])
            x, current, i = x + length, Fraction(t), i + 1
        elif not east and x - length - width >= 0:
            put(t, "conveyor", "linear", x - length, y - width, 180)
            put(t, "coalition", None, x - length, y - width - MRA_WIDTH, 180, assignment[t])
            x, current, i = x - length, Fraction(t), i + 1
        elif east and y + length + 2 * width <= ye:
            put(current + Fraction(1, 4), "conveyor", "corner", x, y, 0)
            put(current + Fraction(1, 2), "conveyor", "linear", x, y + width, -90)
            put(current + Fraction(3, 4), "conveyor", "corner", x, y + width + length, -90)
            y, east, current = y + length + 2 * width, False, Fraction(t)
        elif not east and y + length + width <= ye:
            put(current + Fraction(1, 4), "conveyor", "corner", x - width, y - width, 90)
            put(current + Fraction(1, 2), "conveyor", "linear", x - width, y, -90)
            put(current + Fraction(3, 4), "conveyor", "corner", x - width, y + length, 180)
            y, east, current = y + length, True, Fraction(t)
        else:
            return None
    return stations


# --- generator ------------------------------------------------------------------------------
def load_system(system: str, mras: str = "", gap: str = "", parts: str = "", patterns: str = ""):
    if system == "custom":
        if not mras.strip() or not gap.strip():
            raise ValueError("system=custom needs mras and gap (and parts when tasks handle parts), "
                             "in the notation of the module docstring")
        texts = (mras, gap, parts)
    elif system in NAMED:
        if mras.strip() or gap.strip() or parts.strip():
            raise ValueError("mras, gap and parts are only used with system=custom")
        texts = NAMED[system]
    else:
        raise ValueError(f"system must be one of {SYSTEMS}, got {system!r}")
    modules = parse_mras(texts[0])
    part_list = parse_parts(texts[2])
    order = complete_order(parse_gap(texts[1]), part_list)
    rules = parse_patterns(patterns if patterns.strip() else PATTERNS)
    clash = {m.name for m in modules} & {order.name}
    if clash:
        raise ValueError(f"the order and the MRAs need different names, both use {sorted(clash)}")
    return modules, order, rules


def _check_floor(floor) -> tuple[int, int]:
    if (not isinstance(floor, list) or len(floor) != 2
            or not all(isinstance(v, int) and not isinstance(v, bool) and v > 0 for v in floor)):
        raise ValueError(f"floor must be [length, depth] with two positive integers, got {floor!r}")
    return floor[0], floor[1]


def generate(p, rng):
    modules, order, patterns = load_system(p.system, p.mras, p.gap, p.parts, p.patterns)
    floor = _check_floor(p.floor)
    names = {m.name for m in modules}
    if not isinstance(p.unavailable, list) or not all(isinstance(n, str) for n in p.unavailable):
        raise ValueError(f"unavailable must be a list of MRA names, got {p.unavailable!r}")
    unknown = sorted(set(p.unavailable) - names)
    if unknown:
        raise ValueError(f"unavailable names unknown MRAs {unknown}; known: {sorted(names)}")
    modules = [m for m in modules if m.name not in set(p.unavailable)]

    species1, reactions1, status1, conflicts = coalition_closure(order, modules, patterns,
                                                                 p.max_species)
    coalitions = [x for x in species1 if isinstance(x, Coalition)]
    budget = p.max_species - len(species1)
    if budget >= 1:
        complete, lines, transitions, status2, generic_first = line_closure(
            order, coalitions, [m.name for m in modules], budget + 1)
    else:
        complete = [c for c in coalitions if c.complete]
        lines, transitions, status2, generic_first = [order], [], "truncated", False

    def line_id(assign) -> str:
        body = ";".join(f"t({t})={'+'.join(complete[i].members)}" for t, i in assign)
        return f"{order.name}{{{body}}}"

    species = []
    for x in species1:
        if isinstance(x, Order):
            species.append(Species(x.name, x.text()))
        elif isinstance(x, Mra):
            species.append(Species(x.name, x.text()))
        else:
            species.append(Species(x.id, x.text()))
    for x in lines[1:]:
        species.append(Species(line_id(x[1]), "\n".join(
            f"t({t}): {', '.join(complete[i].members)}" for t, i in x[1])))

    def mol_id(x) -> str:
        return x.name if isinstance(x, (Order, Mra)) else x.id

    reactions, labels = [], []
    for lhs, rhs, label in reactions1:
        reactions.append(Reaction.of([mol_id(x) for x in lhs], [mol_id(x) for x in rhs]))
        labels.append(label)
    outgoing = set()
    for before, after in transitions:
        old = set() if isinstance(before, Order) else set(before[1])
        new = [pair for pair in after[1] if pair not in old]
        left = order.name if isinstance(before, Order) else line_id(before[1])
        outgoing.add(left)
        reactions.append(Reaction.of([left, *(complete[i].id for _, i in new)], [line_id(after[1])]))
        generic_step = generic_first and isinstance(before, Order)
        labels.append({"rule": "assign-generic" if generic_step else "assign",
                       "task": ",".join(f"t({t})" for t, _ in new)})

    tasks = [t.number for t in order.tasks]
    finals = [x for x in lines[1:] if len(x[1]) == len(tasks)]
    assembly_lines = []
    for x in finals:
        assignment = {t: "+".join(complete[i].members) for t, i in x[1]}
        assembly_lines.append({
            "line": line_id(x[1]),
            "assignment": {f"t({t})": v for t, v in assignment.items()},
            "layout": layout(assignment, floor),
        })
    counts = Counter(c.task for c in coalitions)
    analysis = {
        "order": order.name,
        "tasks": [t.name for t in order.tasks],
        "coalitions": {f"t({t})": counts.get(t, 0) for t in tasks},
        "complete_coalitions": {f"t({t})": ["+".join(c.members) for c in complete if c.task == t]
                                for t in tasks},
        "unassignable_tasks": [f"t({t})" for t in tasks if not any(c.task == t for c in complete)],
        "coalition_conflicts": conflicts,
        "n_assembly_lines": len(assembly_lines),
        "assembly_lines": assembly_lines,
        "dead_end_lines": sum(1 for x in lines[1:]
                              if len(x[1]) < len(tasks) and line_id(x[1]) not in outgoing),
    }
    status = "truncated" if "truncated" in (status1, status2) else "complete"
    return Network(
        species=species,
        reactions=reactions,
        status=status,
        initial_state={order.name: 1.0, **{m.name: 1.0 for m in modules}},
        extras={
            "system": p.system,
            "rules": list(RULES),
            "composite_skills": list(COMPOSITE_SKILLS),
            "part_matching": list(PART_MATCHING),
            "patterns": [f"{_skills_text(a)} & {_skills_text(b)}" for a, b in patterns],
            "reaction_rules": labels,
            "analysis": analysis,
        },
    )
