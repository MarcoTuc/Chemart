"""P systems (membrane computing): Păun (1998/2000). Catalog id: p-systems.

Book 9.5; formal definition and semantics from Păun, "Introduction to Membrane
Computing" (2006). A transition P system

    Π = (O, μ, w_1..w_m, (R_1, ρ_1)..(R_m, ρ_m), i_o)

has a tree of nested membranes μ written with labelled brackets
``[1 [2 ]2 [3 ]3 ]1`` (1 is the skin), a multiset w_i in every region, a set
of multiset-rewriting rules R_i per region with a priority relation ρ_i, and an
output region i_o (a membrane label, or ``env``/``0`` for the environment).

Rule syntax (the explicit reaction syntax plus target suffixes):

    r1: a + c -> c'                 here (default)
    2 a + b -> a + b@out + c@in     out = the surrounding region, in = a random inner membrane
    d + c + c' -> a@in3             in3 = the inner membrane labelled 3
    d -> d + delta                  delta (or δ) dissolves the membrane of the region

The label ``r1:`` is optional (default r<position> within the region);
priorities are ``"r1 > r3"`` strings per region (chains ``r1 > r2 > r3``
allowed; the relation is closed transitively).

Network: every rule instantiated over compartment-qualified species
``object@region``. A target resolves to every region it can reach in some
structure reachable by dissolution (the parent, or a further ancestor when the
intermediate membranes can dissolve; likewise for inner membranes), so the
network is complete. Reactions carry no rate: P systems have none. Dissolution
and priorities have no formula and are recorded in extras.

Analysis: one computation in the non-deterministic maximally parallel mode
(strong or weak priorities), up to halting or max_steps, with the output
read in the halting configuration.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from itertools import combinations_with_replacement, product
from math import comb, isqrt, prod

from chemart.helpers.explicit import parse as parse_reaction
from chemart.network import Network, Reaction, Species

SYSTEMS = ("divisibility", "n-squared", "custom")
PRIORITY_MODES = ("strong", "weak")
ENV = "env"
DELTA = frozenset({"delta", "δ"})
MAX_REACTIONS = 100_000
TRACE_STEPS = 100
CUSTOM_FIELDS = ("membranes", "objects", "rules", "priorities", "output")


def divisibility(n: int, k: int) -> dict:
    """Book eq. 9.21 and fig. 9.3: is n a multiple of k? (Calude & Păun 2000, Păun 2002)."""
    return {
        "membranes": "[1 [2 ]2 [3 ]3 ]1",
        "objects": {"2": {"a": n, "c": k, "d": 1}},
        "rules": {
            "1": ["r4: d + c + c' -> a@in3"],
            "2": ["r1: a + c -> c'", "r2: a + c' -> c", "r3: d -> d + delta"],
        },
        "priorities": {"2": ["r1 > r3", "r2 > r3"]},
        "output": "3",
    }


#: Păun (2006), Introduction to Membrane Computing, fig. 3: N(Π) = {n^2 | n >= 1}.
N_SQUARED = {
    "membranes": "[1 [2 [3 ]3 ]2 ]1",
    "objects": {"3": "a + f + c"},
    "rules": {
        "3": ["a -> a + b", "a -> b + delta", "f -> 2 f"],
        "2": ["b -> d", "d -> d + e", "2 f -> f", "c + f -> c + d + delta"],
        "1": ["e -> e@out", "f -> f"],
    },
    "priorities": {},
    "output": "env",
}


# --- parsing ------------------------------------------------------------------------
@dataclass
class Rule:
    region: str
    label: str
    text: str
    lhs: dict[str, int]
    products: list[tuple[str, int, str]]    # (object, count, "here" | "out" | "in" | "in:<label>")
    dissolve: bool

    @property
    def id(self) -> str:
        return f"{self.region}/{self.label}"


@dataclass
class System:
    skin: str
    membranes: list[str]                     # document order
    parent: dict[str, str]                   # every membrane but the skin -> its parent
    objects: dict[str, Counter]
    rules: dict[str, list[Rule]]
    higher: dict[str, dict[str, set[str]]]   # region -> label -> labels with priority over it
    output: str

    def children(self, m: str) -> list[str]:
        return [c for c in self.membranes if self.parent.get(c) == m]

    def dissolvable(self, m: str) -> bool:
        return m != self.skin and any(r.dissolve for r in self.rules.get(m, ()))

    def reach_up(self, m: str) -> list[str]:
        """Regions an object sent out of m can land in, over all dissolution histories."""
        if m == self.skin:
            return [ENV]
        out, p = [], self.parent[m]
        while True:
            out.append(p)
            if not self.dissolvable(p):
                return out
            p = self.parent[p]

    def reach_down(self, m: str) -> list[str]:
        """Membranes that are, or can become (by dissolution), directly inside m."""
        out, queue = [], self.children(m)
        while queue:
            c = queue.pop(0)
            out.append(c)
            if self.dissolvable(c):
                queue.extend(self.children(c))
        return out


_BRACKET = re.compile(r"\s*(?:\[\s*(\w+)|\]\s*(\w+))")
_LABEL = re.compile(r"^\s*(\w+)\s*:\s*(.*)$", re.S)


def parse_membranes(text: str) -> tuple[str, list[str], dict[str, str]]:
    """'[1 [2 ]2 [3 ]3 ]1' -> (skin, labels in order, parent map)."""
    text = text.strip()
    if not text:
        raise ValueError("membranes is empty: give a structure such as '[1 [2 ]2 [3 ]3 ]1'")
    pos, stack, order, parent, skin = 0, [], [], {}, None
    while pos < len(text):
        m = _BRACKET.match(text, pos)
        if not m:
            raise ValueError(f"cannot read membranes {text!r} at {text[pos:]!r}: expected '[label' or ']label'")
        pos = m.end()
        if m.group(1) is not None:
            label = m.group(1)
            if label in (ENV, "0"):
                raise ValueError(f"membrane label {label!r} is reserved for the environment")
            if label in order:
                raise ValueError(f"membrane label {label!r} is used twice in {text!r}")
            if stack:
                parent[label] = stack[-1]
            elif skin is None:
                skin = label
            else:
                raise ValueError(f"membranes {text!r} must have a single skin membrane enclosing all others")
            order.append(label)
            stack.append(label)
        else:
            if not stack or stack[-1] != m.group(2):
                raise ValueError(f"membranes {text!r}: ']{m.group(2)}' does not close the innermost open membrane")
            stack.pop()
    if stack:
        raise ValueError(f"membranes {text!r}: membrane {stack[-1]!r} is never closed")
    return skin, order, parent


def render(skin: str, membranes: list[str], parent: dict[str, str]) -> str:
    def one(m):
        inner = "".join(one(c) + " " for c in membranes if parent.get(c) == m)
        return f"[{m} {inner}]{m}"
    return one(skin)


def _object(name: str, where: str) -> str:
    if not name or "@" in name or ":" in name or name in DELTA:
        raise ValueError(f"{where}: invalid object name {name!r} (no '@' or ':', and delta/δ is reserved)")
    return name


def parse_multiset(value, where: str) -> Counter:
    if isinstance(value, dict):
        out = Counter()
        for obj, n in value.items():
            if not isinstance(n, int) or isinstance(n, bool) or n < 0:
                raise ValueError(f"{where}: multiplicity of {obj!r} must be a non-negative integer, got {n!r}")
            out[_object(obj, where)] += n
        return +out
    if isinstance(value, str):
        if "->" in value:
            raise ValueError(f"{where}: a multiset cannot contain '->', got {value!r}")
        lhs, _ = parse_reaction(value + " -> ")
        return Counter({_object(o, where): n for o, n in lhs.items()})
    raise ValueError(f"{where}: a multiset is a string like '3 a + c' or an object {{'a': 3}}, got {value!r}")


def parse_rule(region: str, index: int, text: str, membranes: list[str]) -> Rule:
    if not isinstance(text, str):
        raise ValueError(f"rules[{region!r}]: every rule must be a string, got {text!r}")
    m = _LABEL.match(text)
    label, body = (m.group(1), m.group(2)) if m else (f"r{index}", text)
    where = f"rule {text!r} in region {region}"
    lhs_raw, rhs_raw = parse_reaction(body)
    if not lhs_raw:
        raise ValueError(f"{where}: the left-hand side must not be empty")
    lhs = {_object(o, where): n for o, n in lhs_raw.items()}
    dissolve, products = False, Counter()
    for name, n in rhs_raw.items():
        if name in DELTA:
            if n != 1:
                raise ValueError(f"{where}: write delta once")
            dissolve = True
            continue
        obj, _, target = name.partition("@")
        target = target or "here"
        if target not in ("here", "out", "in"):
            if not (target.startswith("in") and target[2:] in membranes):
                raise ValueError(f"{where}: target {target!r} must be here, out, in or in<label> with an existing membrane label")
            target = "in:" + target[2:]
        products[(_object(obj, where), target)] += n
    return Rule(region, label, text, lhs, [(o, n, t) for (o, t), n in products.items()], dissolve)


def _as_list(value, where: str) -> list:
    if isinstance(value, str):
        return [line for line in value.splitlines() if line.strip()]
    if isinstance(value, list):
        return value
    raise ValueError(f"{where} must be a list of strings, got {value!r}")


def build_system(spec: dict) -> System:
    skin, membranes, parent = parse_membranes(spec.get("membranes", ""))

    def region(key, field):
        key = str(key)
        if key not in membranes:
            raise ValueError(f"{field}: region {key!r} is not a membrane of {spec['membranes']!r}")
        return key

    objects = {region(k, "objects"): parse_multiset(v, f"objects[{k!r}]")
               for k, v in (spec.get("objects") or {}).items()}
    rules: dict[str, list[Rule]] = {}
    for k, texts in (spec.get("rules") or {}).items():
        r = region(k, "rules")
        rules[r] = [parse_rule(r, i + 1, t, membranes) for i, t in enumerate(_as_list(texts, f"rules[{k!r}]"))]
        labels = [x.label for x in rules[r]]
        if len(set(labels)) != len(labels):
            raise ValueError(f"rules[{k!r}]: rule labels must be unique within a region, got {labels}")
        for x in rules[r]:
            if x.dissolve and r == skin:
                raise ValueError(f"rule {x.text!r}: the skin membrane {skin} can never be dissolved")
    higher: dict[str, dict[str, set[str]]] = {}
    for k, entries in (spec.get("priorities") or {}).items():
        r = region(k, "priorities")
        known = {x.label for x in rules.get(r, [])}
        above: dict[str, set[str]] = {}
        for entry in _as_list(entries, f"priorities[{k!r}]"):
            chain = [s.strip() for s in str(entry).split(">")]
            if len(chain) < 2 or any(s not in known for s in chain):
                raise ValueError(f"priorities[{k!r}]: {entry!r} must be 'r1 > r2' over rule labels {sorted(known)} of region {r}")
            for hi, lo in zip(chain, chain[1:]):
                above.setdefault(lo, set()).add(hi)
        changed = True
        while changed:                                   # transitive closure
            changed = False
            for lo in list(above):
                extra = set().union(*(above.get(h, set()) for h in above[lo])) - above[lo]
                if extra:
                    above[lo] |= extra
                    changed = True
        if any(lo in his for lo, his in above.items()):
            raise ValueError(f"priorities[{k!r}] contain a cycle; a priority relation is a strict partial order")
        higher[r] = above
    output = str(spec.get("output", "")).strip()
    if output in ("0", ENV):
        output = ENV
    elif output not in membranes:
        raise ValueError(f"output must be a membrane label of {spec['membranes']!r} or 'env' (0), got {output!r}")
    return System(skin, membranes, parent, objects, rules, higher, output)


# --- the network ----------------------------------------------------------------------
def instantiate(system: System) -> tuple[list[Reaction], list[str]]:
    """Reactions of every rule over every reachable target; returns (reactions, rule id of each)."""
    reactions: list[Reaction] = []
    rule_of: list[str] = []
    for m in system.membranes:
        for rule in system.rules.get(m, []):
            lhs = {f"{o}@{m}": n for o, n in rule.lhs.items()}
            options = []
            for obj, n, target in rule.products:
                if target == "here":
                    regions = [m]
                elif target == "out":
                    regions = system.reach_up(m)
                elif target == "in":
                    regions = system.reach_down(m)
                else:
                    regions = [target[3:]] if target[3:] in system.reach_down(m) else []
                options.append((obj, n, regions))
            count = prod(comb(len(rs) + n - 1, n) for _, n, rs in options)
            if len(reactions) + count > MAX_REACTIONS:
                raise ValueError(f"rule {rule.text!r} expands to {count} reactions; more than {MAX_REACTIONS} in total")
            seen = set()
            for combo in product(*(_distributions(o, n, rs) for o, n, rs in options)):
                rhs = Counter()
                for part in combo:
                    rhs.update(part)
                key = tuple(sorted(rhs.items()))
                if key not in seen:
                    seen.add(key)
                    reactions.append(Reaction(dict(lhs), dict(rhs)))
                    rule_of.append(rule.id)
    return reactions, rule_of


def _distributions(obj: str, n: int, regions: list[str]):
    if len(regions) == 1:
        return [{f"{obj}@{regions[0]}": n}]
    return [dict(Counter(f"{obj}@{r}" for r in c)) for c in combinations_with_replacement(regions, n)]


# --- the maximally parallel computation ---------------------------------------------------
def _uniform(rng, n: int) -> int:
    """Uniform integer in [0, n), exact also beyond int64."""
    if n <= 2**62:
        return int(rng.integers(n))
    bits = n.bit_length()
    while True:
        v, got = 0, 0
        while got < bits:
            v = (v << 62) | int(rng.integers(2**62))
            got += 62
        v >>= got - bits
        if v < n:
            return v


def _max_apps(rule: Rule, rem: Counter, children: list[str]) -> int:
    for _, _, target in rule.products:
        if (target == "in" and not children) or (target.startswith("in:") and target[3:] not in children):
            return 0          # Păun 2006 sec. 7: an 'in' with no such inner membrane cannot be followed
    return min(rem[o] // n for o, n in rule.lhs.items())


def select(rules: list[Rule], rem: Counter, children: list[str], higher: dict[str, set[str]],
           priority: str, rng) -> Counter:
    """A maximal multiset of rules for one region; consumes its reactants from `rem`.

    Rules and objects are assigned at random until no further assignment is
    possible: a rule is drawn with weight equal to the number of times it
    could still be applied, and applied a uniform number of times in [1, that].
    strong: a rule is not used if a rule with priority over it is applicable in
    the configuration; weak: it is used only on objects that the rules with
    priority over it cannot take any more (Păun 2006 sec. 11).
    """
    by_label = {r.label: r for r in rules}
    start = {r.label: _max_apps(r, rem, children) > 0 for r in rules}
    picks = Counter()
    while True:
        cands = []
        for r in rules:
            m = _max_apps(r, rem, children)
            if m == 0:
                continue
            above = higher.get(r.label, ())
            if priority == "strong":
                blocked = any(start[h] for h in above)
            else:
                blocked = any(_max_apps(by_label[h], rem, children) > 0 for h in above)
            if not blocked:
                cands.append((r, m))
        if not cands:
            return picks
        x = _uniform(rng, sum(m for _, m in cands))
        for r, m in cands:
            if x < m:
                break
            x -= m
        t = 1 + _uniform(rng, m)
        for o, n in r.lhs.items():
            rem[o] -= n * t
        picks[r.label] += t


def _split(rng, total: int, k: int) -> list[int]:
    if total < 2**62:
        return [int(x) for x in rng.multinomial(total, [1.0 / k] * k)]
    base = [total // k] * k                               # beyond int64: even split, random remainder
    for i in rng.choice(k, size=total % k, replace=False):
        base[int(i)] += 1
    return base


def _config(membranes, contents, env) -> dict:
    out = {m: dict(sorted((+contents[m]).items())) for m in membranes}
    out[ENV] = dict(sorted((+env).items()))
    return out


def run(system: System, rng, priority: str = "strong", max_steps: int = 1000,
        trace_steps: int = TRACE_STEPS) -> dict:
    """One computation; returns the halting (or last) configuration and the output."""
    if priority not in PRIORITY_MODES:
        raise ValueError(f"priority must be one of {PRIORITY_MODES}, got {priority!r}")
    parent = dict(system.parent)
    alive = list(system.membranes)
    contents = {m: Counter(system.objects.get(m, {})) for m in alive}
    env: Counter = Counter()
    applications: Counter = Counter()
    dissolved, trace = [], []
    steps, halted = 0, False
    while True:
        children = {m: [c for c in alive if parent.get(c) == m] for m in alive}
        chosen, remaining = {}, {}
        for m in alive:
            rules = system.rules.get(m, [])
            if not rules:
                continue
            rem = Counter(contents[m])
            picks = select(rules, rem, children[m], system.higher.get(m, {}), priority, rng)
            if picks:
                chosen[m], remaining[m] = picks, rem
        if not chosen:
            halted = True
            break
        if steps >= max_steps:
            break
        steps += 1
        new = {m: +remaining[m] if m in remaining else Counter(contents[m]) for m in alive}
        dissolve, applied = set(), {}
        for m, picks in chosen.items():
            by_label = {r.label: r for r in system.rules[m]}
            for label, t in picks.items():
                rule = by_label[label]
                applied[rule.id] = t
                applications[rule.id] += t
                if rule.dissolve:
                    dissolve.add(m)
                for obj, n, target in rule.products:
                    amount = n * t
                    if target == "here":
                        new[m][obj] += amount
                    elif target == "out":
                        (env if m == system.skin else new[parent[m]])[obj] += amount
                    elif target == "in":
                        for c, x in zip(children[m], _split(rng, amount, len(children[m]))):
                            new[c][obj] += x
                    else:
                        new[target[3:]][obj] += amount
        if dissolve:
            def survivor(x):
                while x in dissolve:
                    x = parent[x]
                return x
            for d in [m for m in alive if m in dissolve]:
                new[survivor(parent[d])].update(new.pop(d))
                dissolved.append({"membrane": d, "step": steps})
            for m in alive:
                if m not in dissolve and m in parent:
                    parent[m] = survivor(parent[m])
            for d in dissolve:
                del parent[d]
            alive = [m for m in alive if m not in dissolve]
        contents = new
        if steps <= trace_steps:
            trace.append({
                "step": steps, "applied": applied,
                "dissolved": [m for m in system.membranes if m in dissolve],
                "configuration": _config(alive, contents, env),
            })
    if system.output == ENV:
        alive_output, final = True, +env
    else:
        alive_output, final = system.output in alive, +contents.get(system.output, Counter())
    result = sum(final.values()) if halted and alive_output else None
    return {
        "semantics": f"non-deterministic maximally parallel, {priority} priorities",
        "priority": priority,
        "halted": halted,
        "steps": steps,
        "max_steps": max_steps,
        "output_region": system.output,
        "output_dissolved": not alive_output,
        "output": dict(sorted(final.items())) if result is not None else None,
        "result": result,
        "final_structure": render(system.skin, alive, parent),
        "final_configuration": _config(alive, contents, env),
        "dissolved": dissolved,
        "applications": dict(applications),
        "trace": trace,
        "trace_truncated": steps > trace_steps,
    }


# --- generator ---------------------------------------------------------------------------
def _spec(p) -> dict:
    custom = {f: getattr(p, f) for f in CUSTOM_FIELDS}
    if p.system == "custom":
        return custom
    given = [f for f, v in custom.items() if v]
    if given:
        raise ValueError(f"{', '.join(given)} only apply with system='custom' (system={p.system!r})")
    return divisibility(p.n, p.k) if p.system == "divisibility" else N_SQUARED


def generate(p, rng) -> Network:
    system = build_system(_spec(p))
    reactions, rule_of = instantiate(system)

    order: dict[str, None] = {}
    initial = {}
    for m in system.membranes:
        for obj, n in system.objects.get(m, {}).items():
            order[f"{obj}@{m}"] = None
            initial[f"{obj}@{m}"] = float(n)
    for r in reactions:
        for s in (*r.reactants, *r.products):
            order.setdefault(s, None)
    species = [Species(s) for s in order]

    analysis = run(system, rng, p.priority, p.max_steps)
    if p.system == "divisibility":
        analysis["n"], analysis["k"] = p.n, p.k
        analysis["n_multiple_of_k"] = None if analysis["result"] is None else analysis["result"] == 0
    elif p.system == "n-squared":
        res = analysis["result"]
        analysis["result_is_square"] = None if res is None else isqrt(res) ** 2 == res

    all_rules = [r for m in system.membranes for r in system.rules.get(m, [])]
    reactions_of = {r.id: [i for i, rid in enumerate(rule_of) if rid == r.id] for r in all_rules}
    compartments = {}
    for m in system.membranes:
        compartments[m] = {
            "parent": system.parent.get(m, ENV),
            "children": system.children(m),
            "skin": m == system.skin,
            "dissolvable": system.dissolvable(m),
            "rules": [r.id for r in system.rules.get(m, [])],
            "species": [s for s in order if s.rsplit("@", 1)[1] == m],
        }
    compartments[ENV] = {
        "parent": None, "children": [system.skin], "environment": True, "rules": [],
        "species": [s for s in order if s.rsplit("@", 1)[1] == ENV],
    }
    extras = {
        "system": p.system,
        "membranes": render(system.skin, system.membranes, system.parent),
        "compartments": compartments,
        "rules": [{"id": r.id, "region": r.region, "label": r.label, "rule": r.text,
                   "dissolves": r.dissolve, "reactions": reactions_of[r.id]} for r in all_rules],
        "reaction_rules": rule_of,
        "priorities": {m: sorted([hi, lo] for lo, his in above.items() for hi in his)
                       for m, above in system.higher.items() if above},
        "priority_semantics": p.priority,
        "dissolution": {
            m: {"rules": [r.id for r in system.rules[m] if r.dissolve],
                "contents_move_to": system.reach_up(m),
                "effect": "the membrane disappears with its rules; its objects and inner membranes "
                          "join the nearest upper region that is not dissolved in the same step"}
            for m in system.membranes if system.dissolvable(m)
        },
        "output_region": system.output,
        "analysis": analysis,
    }
    never = [r.id for r in all_rules if not reactions_of[r.id]]
    if never:
        extras["inapplicable_rules"] = never
    return Network(species=species, reactions=reactions, status="complete",
                   initial_state=initial, extras=extras)
