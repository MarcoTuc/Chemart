"""Kappa calculus: rule-based rewriting of site graphs, flattened into a reaction network.

Catalog id: kappa-calculus. Book 9.7; Danos & Laneve, "Formal molecular
biology", TCS 325 (2004) [208]; the Kappa language reference guide v4
(Fontana, Boutillier, Feret & Krivine); KaDE (Camporesi, Feret & Ly, CMSB 2017).

An agent (a protein) has named sites. A site may carry an internal state and a
link to one site of another agent. A complex (a species) is a connected site
graph; a rule rewrites a pattern, a partially specified site graph, so one rule
stands for many reactions. The network is the closure of the rules over the
complexes reachable from the initial ones, with complexes identified up to
isomorphism.

Notation (a subset of Kappa 4, the syntax of the reference guide):

    signature   C(x1{u p},x2{u p})           states of a site; the first is the default
    pattern     A(x[.],c[1]),C(x1{u}[1])     [.] free, [n] bond n (twice per pattern),
                                             [_] bound to anything, [s.T] bound to site s
                                             of an agent of type T, [#] and {#} anything;
                                             a site that is not written is not tested
    rule        'label' L -> R @ k           arrow notation; agents pair up by position
                'label' L <-> R @ k, k_op    reversible rule (the reverse is 'label_op')
                'label' A(x[./1]),B(x[./1]) @ k      edit notation (before/after)
                ... @ k {k_uni}               two-component rule: k when the components lie
                                             in two complexes, k_uni when they lie in one
    init        {"A(),B()": 1000}            copies of a (possibly disconnected) expression;
                                             unwritten sites take the default state, free

Rules change internal states and create or delete bonds. Agent creation and
deletion are not supported, so the agents of every type are conserved.

Kinetics (KaDE's differential semantics). The disconnected components of a
rule's left-hand side match distinct complexes, and match inside a single
complex only when a unimolecular rate {k_uni} is given. For a reaction R with
reactant multiplicities m_i, k_R = sum over the rule's embeddings into the
reactant copies that produce R of  gamma * Omega / prod m_i!,  with gamma the
rule rate and Omega the symmetry factor of the chosen convention (reference
guide 3.3.1): "rule" 1/|automorphisms of L preserved in R| (stance D, the
biochemist convention), "lhs" 1/|Aut(L)| (stance ND, the SSA convention),
"none" 1 (KaSim's default). The mass-action flux k_R prod x_i^m_i then equals
the rule activity of the reference guide, eq. (1), in molecule counts.
"""

from __future__ import annotations

import re
from collections import Counter
from itertools import permutations, product
from math import factorial, prod

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

MODELS = ("abc", "dimerization", "polymerization", "egfr", "custom")
CONVENTIONS = ("rule", "lhs", "none")

_AGENT = re.compile(r"([A-Za-z][A-Za-z0-9_+\-]*)\s*\((.*)\)$", re.S)
_SITE = re.compile(r"([A-Za-z0-9_][A-Za-z0-9_~+\-]*)((?:\{[^{}]*\}|\[[^\[\]]*\])*)$")
_STUB = re.compile(r"([A-Za-z0-9_][A-Za-z0-9_~+\-]*)\.([A-Za-z][A-Za-z0-9_+\-]*)$")
_NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_RATE = re.compile(rf"({_NUMBER})\s*(?:\{{\s*({_NUMBER})\s*\}})?$")


# ============================================================================
# Parsing
# ============================================================================
def _split_top(text: str, seps: str) -> list[str]:
    """Split at separator characters that are outside (), {} and []."""
    out, cur, depth = [], [], 0
    for ch in text:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
            if depth < 0:
                raise ValueError(f"unbalanced brackets in {text!r}")
        if depth == 0 and ch in seps:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if depth:
        raise ValueError(f"unbalanced brackets in {text!r}")
    out.append("".join(cur))
    return out


def _agent_tokens(text: str, where: str) -> tuple[str, dict[str, tuple[str | None, str | None]]]:
    """'C(x1{u}[1],x2)' -> ('C', {'x1': ('u', '1'), 'x2': (None, None)}), raw tokens."""
    text = text.strip()
    if text == ".":
        raise ValueError(f"empty slot '.' in {where!r}: agent creation and deletion are not supported")
    if re.search(r"\)\s*[+-]$", text):
        raise ValueError(f"'{text}' in {where!r}: agent creation (+) and deletion (-) are not supported")
    m = _AGENT.match(text)
    if not m:
        raise ValueError(f"cannot read agent {text!r} in {where!r}: expected Name(site, ...)")
    tokens: list[str] = []
    for tok in _split_top(m.group(2), ", \t\n"):
        if not tok:
            continue
        if tokens and tok[0] in "{[":
            tokens[-1] += tok
        else:
            tokens.append(tok)
    sites: dict[str, tuple[str | None, str | None]] = {}
    for tok in tokens:
        sm = _SITE.match(tok)
        if not sm:
            raise ValueError(f"cannot read site {tok!r} of agent {text!r} in {where!r}")
        name = sm.group(1)
        if name in sites:
            raise ValueError(f"site {name!r} appears twice in agent {text!r} in {where!r}")
        groups = re.findall(r"\{[^{}]*\}|\[[^\[\]]*\]", sm.group(2))
        states = [g[1:-1].strip() for g in groups if g[0] == "{"]
        links = [g[1:-1].strip() for g in groups if g[0] == "["]
        if len(states) > 1 or len(links) > 1:
            raise ValueError(f"site {name!r} of {text!r} in {where!r} has more than one internal or link state")
        sites[name] = (states[0] if states else None, links[0] if links else None)
    return m.group(1), sites


def _agents(text: str, where: str) -> list[tuple[str, dict]]:
    return [_agent_tokens(a, where) for a in _split_top(text, ",") if a.strip()]


def signatures(texts) -> dict[str, dict[str, tuple[str, ...]]]:
    """['A(x,c)', 'C(x1{u p},x2{u p})'] -> {agent: {site: states}} (empty tuple: no internal state)."""
    if not isinstance(texts, list) or not texts or not all(isinstance(t, str) for t in texts):
        raise ValueError(f"signatures must be a non-empty list of strings like 'C(x1{{u p}},x2)', got {texts!r}")
    sig: dict[str, dict[str, tuple[str, ...]]] = {}
    for text in texts:
        body = text.strip()
        if body.startswith("%agent:"):
            body = body[len("%agent:"):]
        name, raw = _agent_tokens(body, text)
        if name in sig:
            raise ValueError(f"agent {name!r} is declared twice in signatures")
        decl = {}
        for site, (st, _binding_types) in raw.items():
            if st is not None and ("=" in st or "/" in st):
                raise ValueError(f"signature {text!r}: counters are not supported")
            states = tuple(st.split()) if st else ()
            if len(set(states)) != len(states):
                raise ValueError(f"signature {text!r}: site {site!r} lists a state twice")
            decl[site] = states
        sig[name] = decl
    return sig


def _link(tok: str | None, where: str):
    """None (untested), '#', '.', '_', ('stub', site, agent) or ('bond', n)."""
    if tok is None or tok in ("#", ".", "_"):
        return tok
    if tok.isdigit():
        return ("bond", int(tok))
    m = _STUB.match(tok)
    if m:
        return ("stub", m.group(1), m.group(2))
    raise ValueError(f"cannot read link state [{tok}] in {where!r}")


class Pattern:
    """Agents with site constraints: sites[i][site] = (state, link)."""

    def __init__(self, agents: list[tuple[str, dict]], sig: dict, where: str):
        self.types: list[str] = []
        self.sites: list[dict] = []
        ends: dict[int, list[tuple[int, str]]] = {}
        for i, (name, raw) in enumerate(agents):
            if name not in sig:
                raise ValueError(f"agent {name!r} in {where!r} is not declared in signatures")
            decl = sig[name]
            sites = {}
            for s, (st, lk) in raw.items():
                if s not in decl:
                    raise ValueError(f"agent {name} has no site {s!r} (in {where!r}); declared: {list(decl)}")
                if st is not None and not decl[s]:
                    raise ValueError(f"site {name}.{s} has no internal states (in {where!r})")
                if st not in (None, "#") and st not in decl[s]:
                    raise ValueError(f"site {name}.{s} has no internal state {st!r} (in {where!r}); "
                                     f"declared: {list(decl[s])}")
                link = _link(lk, where)
                if isinstance(link, tuple) and link[0] == "stub":
                    if link[2] not in sig or link[1] not in sig[link[2]]:
                        raise ValueError(f"binding type {link[1]}.{link[2]} in {where!r} names an undeclared site")
                if isinstance(link, tuple) and link[0] == "bond":
                    ends.setdefault(link[1], []).append((i, s))
                sites[s] = (st, link)
            self.types.append(name)
            self.sites.append(sites)
        self.partner: dict[tuple[int, str], tuple[int, str]] = {}
        for n, e in ends.items():
            if len(e) != 2:
                raise ValueError(f"bond label {n} must appear exactly twice in {where!r}, found {len(e)}")
            self.partner[e[0]], self.partner[e[1]] = e[1], e[0]

    def components(self) -> list[list[int]]:
        root = list(range(len(self.types)))

        def find(i):
            while root[i] != i:
                root[i] = root[root[i]]
                i = root[i]
            return i

        for (a, _), (b, _) in self.partner.items():
            root[find(a)] = find(b)
        groups: dict[int, list[int]] = {}
        for i in range(len(self.types)):
            groups.setdefault(find(i), []).append(i)
        return list(groups.values())

    def _local(self, i: int):
        return tuple(sorted((s, repr(st), repr("bond" if isinstance(lk, tuple) and lk[0] == "bond" else lk))
                            for s, (st, lk) in self.sites[i].items()))

    def preserved_by(self, sigma: list[int]) -> bool:
        return all(self.types[sigma[i]] == self.types[i] and self._local(sigma[i]) == self._local(i)
                   for i in range(len(sigma))) and all(
            self.partner.get((sigma[a], s)) == (sigma[b], t) for (a, s), (b, t) in self.partner.items())


def automorphisms(p: Pattern, also: Pattern | None = None) -> int:
    """Number of agent permutations that map p onto itself (and `also` onto itself)."""
    n = len(p.types)
    keys = [(p.types[i], p._local(i), also._local(i) if also else None) for i in range(n)]
    count, sigma, used = 0, [0] * n, [False] * n

    def extend(i):
        nonlocal count
        if i == n:
            count += p.preserved_by(sigma) and (also is None or also.preserved_by(sigma))
            return
        for j in range(n):
            if not used[j] and keys[j] == keys[i]:
                used[j], sigma[i] = True, j
                extend(i + 1)
                used[j] = False

    extend(0)
    return count


class Rule:
    """A one-way rule L -> R: agents pair up by position."""

    def __init__(self, label: str, lhs: Pattern, rhs: Pattern, rate: float, uni: float | None, text: str):
        self.label, self.lhs, self.rhs, self.rate, self.uni, self.text = label, lhs, rhs, rate, uni, text
        if len(lhs.types) != len(rhs.types):
            raise ValueError(f"rule {text!r}: both sides need the same agents (creation and deletion are not supported)")
        if not lhs.types:
            raise ValueError(f"rule {text!r} has no agents")
        self.state_changes: list[tuple[int, str, str]] = []
        self.unbind: list[tuple[int, str]] = []
        self.bind: list[tuple[tuple[int, str], tuple[int, str]]] = []
        for a, (tl, tr) in enumerate(zip(lhs.types, rhs.types)):
            if tl != tr:
                raise ValueError(f"rule {text!r}: slot {a + 1} holds {tl} on the left and {tr} on the right")
            if set(lhs.sites[a]) != set(rhs.sites[a]):
                raise ValueError(f"rule {text!r}: agent {tl} (slot {a + 1}) must mention the same sites on both sides")
            for s, (sl, ll) in lhs.sites[a].items():
                sr, lr = rhs.sites[a][s]
                if (sl is None) != (sr is None) or (ll is None) != (lr is None):
                    raise ValueError(f"rule {text!r}: site {tl}.{s} must specify the same kind of state on both sides")
                if sr == "#" and sl != "#":
                    raise ValueError(f"rule {text!r}: site {tl}.{s} cannot become an unspecified state")
                if sr not in (None, "#") and sr != sl:
                    self.state_changes.append((a, s, sr))
                if lr in ("#", "_") or (isinstance(lr, tuple) and lr[0] == "stub"):
                    if lr != ll:
                        raise ValueError(f"rule {text!r}: site {tl}.{s} cannot become the unspecified link [{_show(lr)}]")
                elif lr == ".":
                    if ll != ".":
                        self.unbind.append((a, s))
                elif isinstance(lr, tuple):
                    other = rhs.partner[(a, s)]
                    if not (isinstance(ll, tuple) and ll[0] == "bond" and lhs.partner[(a, s)] == other):
                        if ll != ".":
                            self.unbind.append((a, s))
                        if (a, s) < other:
                            self.bind.append(((a, s), other))
        self.components = lhs.components()
        if uni is not None and len(self.components) != 2:
            raise ValueError(f"rule {text!r}: a unimolecular rate {{k}} needs a left-hand side of exactly two "
                             f"connected components, found {len(self.components)}")
        self.aut_lhs = automorphisms(lhs)
        self.aut_rule = automorphisms(lhs, rhs)

    def omega(self, convention: str) -> float:
        return {"rule": 1 / self.aut_rule, "lhs": 1 / self.aut_lhs, "none": 1.0}[convention]


def _show(link) -> str:
    if isinstance(link, tuple):
        return str(link[1]) if link[0] == "bond" else f"{link[1]}.{link[2]}"
    return link


def _rates(text: str, rule: str) -> list[tuple[float, float | None]]:
    out = []
    for part in _split_top(text, ","):
        m = _RATE.match(part.strip())
        if not m:
            raise ValueError(f"rule {rule!r}: cannot read rate {part.strip()!r}; rates are numbers, "
                             f"optionally followed by a unimolecular rate in braces (0.001 {{0.1}})")
        k, uni = float(m.group(1)), (float(m.group(2)) if m.group(2) is not None else None)
        if k < 0 or (uni is not None and uni < 0):
            raise ValueError(f"rule {rule!r}: rates must be non-negative")
        out.append((k, uni))
    return out


def parse_rules(texts, sig: dict) -> list[Rule]:
    """Kappa rule strings -> one-way Rules (a reversible rule gives 'label' and 'label_op')."""
    if not isinstance(texts, list) or not texts or not all(isinstance(t, str) for t in texts):
        raise ValueError(f"rules must be a non-empty list of rule strings, got {texts!r}")
    out: list[Rule] = []
    labels: set[str] = set()
    for i, text in enumerate(texts):
        body = " ".join(text.split())
        m = re.match(r"'([^']*)'\s*(.*)$", body)
        label, body = (m.group(1), m.group(2)) if m else (f"r{i + 1}", body)
        if label in labels:
            raise ValueError(f"rule label {label!r} is used twice")
        labels.add(label)
        if "@" not in body:
            raise ValueError(f"rule {text!r} needs its rate after '@'")
        expr, rate_text = body.rsplit("@", 1)
        rates = _rates(rate_text, text)
        if "<->" in expr:
            if expr.count("<->") != 1 or len(rates) != 2:
                raise ValueError(f"reversible rule {text!r} needs one '<->' and two rates 'k, k_op'")
            left, right = expr.split("<->")
            L, R = Pattern(_agents(left, text), sig, text), Pattern(_agents(right, text), sig, text)
            out.append(Rule(label, L, R, rates[0][0], rates[0][1], text))
            out.append(Rule(label + "_op", R, L, rates[1][0], rates[1][1], text))
            continue
        if len(rates) != 1:
            raise ValueError(f"rule {text!r} needs exactly one rate")
        if "->" in expr:
            if expr.count("->") != 1:
                raise ValueError(f"rule {text!r} needs exactly one '->'")
            left, right = expr.split("->")
            L, R = Pattern(_agents(left, text), sig, text), Pattern(_agents(right, text), sig, text)
        else:
            L, R = _edit(expr, sig, text)
        out.append(Rule(label, L, R, rates[0][0], rates[0][1], text))
    return out


def _edit(expr: str, sig: dict, text: str) -> tuple[Pattern, Pattern]:
    """Edit notation: before/after written inside the states, A(x{u/p}[./1])."""
    left, right = [], []
    for name, raw in _agents(expr, text):
        lsites, rsites = {}, {}
        for s, (st, lk) in raw.items():
            sl, sr = (st.split("/", 1) if st is not None and "/" in st else (st, st))
            kl, kr = (lk.split("/", 1) if lk is not None and "/" in lk else (lk, lk))
            lsites[s] = (sl.strip() if sl is not None else None, kl.strip() if kl is not None else None)
            rsites[s] = (sr.strip() if sr is not None else None, kr.strip() if kr is not None else None)
        left.append((name, lsites))
        right.append((name, rsites))
    return Pattern(left, sig, text), Pattern(right, sig, text)


# ============================================================================
# Complexes: canonical forms, embeddings and rule application
# ============================================================================
class Complex:
    """A connected, fully specified site graph with agents in canonical order."""

    __slots__ = ("types", "states", "links")

    def __init__(self, types, states, links):
        self.types, self.states, self.links = types, states, links


class Registry:
    """Canonical complexes, their embeddings cache, and the rule application."""

    def __init__(self, sig: dict, rules: list[Rule], max_size: int):
        self.sig, self.rules, self.max_size = sig, rules, max_size
        self.complexes: dict[str, Complex] = {}
        self._emb: dict[tuple[str, int, int], list[tuple[int, ...]]] = {}
        self.oversized = False

    # -- canonical form ------------------------------------------------------
    def text(self, types, states, links, order) -> str:
        bonds: dict[tuple[int, str], int] = {}
        parts = []
        for g in order:
            sites = []
            for s, allowed in self.sig[types[g]].items():
                p = links.get((g, s))
                if p is None:
                    lk = "."
                else:
                    if (g, s) not in bonds:
                        bonds[(g, s)] = bonds[p] = len(bonds) // 2 + 1
                    lk = str(bonds[(g, s)])
                sites.append(f"{s}{{{states[g][s]}}}[{lk}]" if allowed else f"{s}[{lk}]")
            parts.append(f"{types[g]}({','.join(sites)})")
        return ",".join(parts)

    def add(self, types, states, links, agents: list[int]) -> str:
        """Register the connected component `agents` of a mixture; return its canonical id."""
        head = min(types[g] + "(" for g in agents)
        best = None
        for root in agents:
            if types[root] + "(" != head:
                continue
            order, seen, k = [root], {root}, 0
            while k < len(order):
                g = order[k]
                k += 1
                for s in self.sig[types[g]]:
                    p = links.get((g, s))
                    if p is not None and p[0] not in seen:
                        seen.add(p[0])
                        order.append(p[0])
            code = self.text(types, states, links, order)
            if best is None or code < best[0]:
                best = (code, order)
        code, order = best
        if code not in self.complexes:
            new = {g: i for i, g in enumerate(order)}
            self.complexes[code] = Complex(
                tuple(types[g] for g in order),
                tuple(dict(states[g]) for g in order),
                {(new[g], s): (new[h], t) for (g, s), (h, t) in links.items() if g in new},
            )
        return code

    def split(self, types, states, links) -> list[list[int]]:
        adjacent: dict[int, list[int]] = {}
        for (g, _), (h, _) in links.items():
            adjacent.setdefault(g, []).append(h)
        seen, out = set(), []
        for g in range(len(types)):
            if g in seen:
                continue
            comp, stack = [], [g]
            seen.add(g)
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adjacent.get(x, ()):
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            out.append(sorted(comp))
        return out

    # -- matching ------------------------------------------------------------
    def embeddings(self, cid: str, ri: int, ci: int) -> list[tuple[int, ...]]:
        """Embeddings of component ci of rule ri into complex cid (targets in component order)."""
        key = (cid, ri, ci)
        if key not in self._emb:
            rule = self.rules[ri]
            self._emb[key] = embed(rule.lhs, rule.components[ci], self.complexes[cid])
        return self._emb[key]

    # -- rewriting -----------------------------------------------------------
    def apply(self, rule: Rule, ids, assign, combo) -> list[str] | None:
        types, states, links, offsets = [], [], {}, []
        for cid in ids:
            c = self.complexes[cid]
            off = len(types)
            offsets.append(off)
            types.extend(c.types)
            states.extend(dict(d) for d in c.states)
            links.update({(g + off, s): (h + off, t) for (g, s), (h, t) in c.links.items()})
        where = {}
        for ci, targets in enumerate(combo):
            off = offsets[assign[ci]]
            for a, g in zip(rule.components[ci], targets):
                where[a] = off + g
        for a, s, new in rule.state_changes:
            states[where[a]][s] = new
        for a, s in rule.unbind:
            p = links.pop((where[a], s), None)
            if p is not None:
                links.pop(p, None)
        for (a, s), (b, t) in rule.bind:
            x, y = (where[a], s), (where[b], t)
            links[x], links[y] = y, x
        comps = self.split(types, states, links)
        if any(len(c) > self.max_size for c in comps):
            self.oversized = True
            return None
        return [self.add(types, states, links, c) for c in comps]


def embed(p: Pattern, component: list[int], c: Complex) -> list[tuple[int, ...]]:
    """All injective maps of the connected pattern component into complex c (Kappa rigidity:
    the image of the first agent fixes the rest)."""
    out = []
    anchor = component[0]
    for g0 in range(len(c.types)):
        if c.types[g0] != p.types[anchor]:
            continue
        m, used, stack, ok = {anchor: g0}, {g0}, [anchor], True
        while stack and ok:
            a = stack.pop()
            g = m[a]
            for s, (st, lk) in p.sites[a].items():
                if st not in (None, "#") and c.states[g][s] != st:
                    ok = False
                    break
                actual = c.links.get((g, s))
                if lk is None or lk == "#":
                    continue
                if lk == ".":
                    ok = actual is None
                elif lk == "_":
                    ok = actual is not None
                elif lk[0] == "stub":
                    ok = actual is not None and actual[1] == lk[1] and c.types[actual[0]] == lk[2]
                else:
                    b, t = p.partner[(a, s)]
                    ok = actual is not None and actual[1] == t and c.types[actual[0]] == p.types[b]
                    if ok and b in m:
                        ok = m[b] == actual[0]
                    elif ok:
                        ok = actual[0] not in used
                        if ok:
                            m[b] = actual[0]
                            used.add(actual[0])
                            stack.append(b)
                if not ok:
                    break
        if ok:
            out.append(tuple(m[a] for a in component))
    return out


def initial_complexes(init: dict, sig: dict, registry: Registry) -> dict[str, float | None]:
    """{expression: copies} -> {complex id: copies}; unwritten sites take the default state, free."""
    out: dict[str, float | None] = {}
    for expr, copies in init.items():
        agents = _agents(expr, expr)
        pat = Pattern(agents, sig, expr)
        types, states, links = list(pat.types), [], {}
        for i, name in enumerate(pat.types):
            st = {}
            for s, allowed in sig[name].items():
                given, link = pat.sites[i].get(s, (None, None))
                if given == "#" or link in ("#", "_") or (isinstance(link, tuple) and link[0] == "stub"):
                    raise ValueError(f"initial complex {expr!r} must be fully specified (no #, _ or binding types)")
                if allowed:
                    st[s] = given if given is not None else allowed[0]
            states.append(st)
        links = {k: v for k, v in pat.partner.items()}
        for comp in registry.split(types, states, links):
            if len(comp) > registry.max_size:
                raise ValueError(f"initial complex {expr!r} has more than max_complex_size = {registry.max_size} agents")
            cid = registry.add(types, states, links, comp)
            if copies is None:
                out.setdefault(cid, None)
            else:
                out[cid] = (out.get(cid) or 0.0) + float(copies)
    return out


def canonical(expr: str, signature_texts: list[str]) -> list[str]:
    """Canonical ids of the connected complexes written in a Kappa expression."""
    sig = signatures(signature_texts)
    return list(initial_complexes({expr: None}, sig, Registry(sig, [], 10 ** 9)))


# ============================================================================
# Published models
# ============================================================================
# Kappa reference guide v4, section 1.3 ("Hello ABC"), rates from its %var lines.
ABC = {
    "signatures": ["A(x,c)", "B(x)", "C(x1{u p},x2{u p})"],
    "rules": [
        "'rule 1' A(x[.]), B(x[.]) <-> A(x[1]), B(x[1]) @ 1.0E-4, 0.1",
        "'rule 2' A(x[_],c[.]), C(x1{u}[.]) -> A(x[_],c[2]), C(x1{u}[2]) @ 1.0E-4",
        "'rule 3' C(x1{u}[1]), A(c[1]) -> C(x1{p}[.]), A(c[.]) @ 1",
        "'rule 4' A(x[.],c[.]), C(x1{p}[.],x2{u}[.]) -> A(x[.],c[1]), C(x1{p}[.],x2{u}[1]) @ 1.0E-4",
        "'rule 5' A(x[.],c[1]), C(x1{p}[.],x2{u}[1]) -> A(x[.],c[.]), C(x1{p}[.],x2{p}[.]) @ 1",
    ],
    "init": {"A(),B()": 1000, "C(x1{u},x2{u})": 10000},
}

# Reference guide 2.3.1.1: 'symmetric dimerization' (strict dimers), no initial amounts published.
DIMERIZATION = {
    "signatures": ["A(x)"],
    "rules": ["'symmetric dimerization' A(x[.]), A(x[.]) -> A(x[1]), A(x[1]) @ 0.001"],
    "init": {"A()": None},
}

# Reference guide 2.3.1.1: 'ambiguous molecularity' (chains grow bimolecularly, close into rings unimolecularly).
POLYMERIZATION = {
    "signatures": ["A(x,y)"],
    "rules": ["'ambiguous molecularity' A(x[.]), A(y[.]) -> A(x[1]), A(y[1]) @ 0.001 {0.1}"],
    "init": {"A()": None},
}

# Blinov, Faeder, Goldstein & Hlavacek (2006), as translated to Kappa 3 in KappaTools
# examples/benchmarks/CMSB2017-KaDe-tool-paper/BNGL/egfr_net.ka (verbatim rules, %var values substituted).
EGFR_VARS = {
    "egf_tot": 1200000, "egfr_tot": 180000, "Grb2_tot": 100000, "Shc_tot": 270000, "Sos_tot": 13000,
    "Grb2_Sos_tot": 49000, "kp1": 1.667e-06, "km1": 0.06, "kp2": 5.556e-06, "km2": 0.1, "kp3": 0.5,
    "km3": 4.505, "kp14": 3, "km14": 0.03, "km16": 0.005, "kp9": 8.333e-07, "km9": 0.05, "kp10": 5.556e-06,
    "km10": 0.06, "kp11": 1.25e-06, "km11": 0.03, "kp13": 2.5e-05, "km13": 0.6, "kp15": 2.5e-07, "km15": 0.3,
    "kp17": 1.667e-06, "km17": 0.1, "kp18": 2.5e-07, "km18": 0.3, "kp19": 5.556e-06, "km19": 0.0214,
    "kp20": 6.667e-08, "km20": 0.12, "kp24": 5e-06, "km24": 0.0429, "kp21": 1.667e-06, "km21": 0.01,
    "kp23": 1.167e-05, "km23": 0.1, "kp12": 5.556e-08, "km12": 0.0015, "kp22": 1.667e-05, "km22": 0.064,
}
EGFR_KAPPA3_INIT = [
    ("egf_tot", "egf(r)"), ("Grb2_tot", "Grb2(SH2,SH3)"), ("Shc_tot", "Shc(PTB,Y317~Y)"), ("Sos_tot", "Sos(dom)"),
    ("egfr_tot", "egfr(l,r,Y1068~Y,Y1148~Y)"), ("Grb2_Sos_tot", "Grb2(SH2,SH3!1),Sos(dom!1)"),
]
EGFR_KAPPA3_RULES = [
    "egfr(l,r) , egf(r) <-> egfr(l!1,r),egf(r!1) @'kp1'{0}, 'km1'",
    "egfr(l!_,r) , egfr(l!_,r) <-> egfr(l!_,r!3),egfr(l!_,r!3) @'kp2'{0}, 'km2'",
    "egfr(r!_,Y1068~Y) -> egfr(r!_,Y1068~pY) @'kp3'",
    "egfr(r!_,Y1148~Y) -> egfr(r!_,Y1148~pY) @'kp3'",
    "egfr(Y1068~pY) -> egfr(Y1068~Y) @'km3'",
    "egfr(Y1148~pY) -> egfr(Y1148~Y) @'km3'",
    "egfr(r!_,Y1148~pY!1),Shc(PTB!1,Y317~Y) -> egfr(r!_,Y1148~pY!1),Shc(PTB!1,Y317~pY) @'kp14'",
    "Shc(PTB!_,Y317~pY) -> Shc(PTB!_,Y317~Y) @'km14'",
    "egfr(Y1068~pY) , Grb2(SH2,SH3) <-> egfr(Y1068~pY!1),Grb2(SH2!1,SH3) @'kp9'{0}, 'km9'",
    "egfr(Y1068~pY) , Grb2(SH2,SH3!_) <-> egfr(Y1068~pY!1),Grb2(SH2!1,SH3!_) @'kp11'{0}, 'km11'",
    "egfr(Y1068~pY!1),Grb2(SH2!1,SH3) ,Sos(dom) <-> egfr(Y1068~pY!1),Grb2(SH2!1,SH3!2),Sos(dom!2) @'kp10'{0}, 'km10'",
    "egfr(Y1148~pY) , Shc(PTB,Y317~Y) <-> egfr(Y1148~pY!1),Shc(PTB!1,Y317~Y) @'kp13'{0}, 'km13'",
    "egfr(Y1148~pY) , Shc(PTB,Y317~pY) <-> egfr(Y1148~pY!1),Shc(PTB!1,Y317~pY) @'kp15'{0}, 'km15'",
    "egfr(Y1148~pY) , Shc(PTB,Y317~pY!1),Grb2(SH2!1,SH3) <-> "
    "egfr(Y1148~pY!2),Shc(PTB!2,Y317~pY!1),Grb2(SH2!1,SH3) @'kp18'{0}, 'km18'",
    "egfr(Y1148~pY) , Shc(PTB,Y317~pY!1),Grb2(SH2!1,SH3!3),Sos(dom!3) <-> "
    "egfr(Y1148~pY!2),Shc(PTB!2,Y317~pY!1),Grb2(SH2!1,SH3!3),Sos(dom!3) @'kp20'{0}, 'km20'",
    "egfr(Y1148~pY!1),Shc(PTB!1,Y317~pY) , Grb2(SH2,SH3) <-> "
    "egfr(Y1148~pY!1),Shc(PTB!1,Y317~pY!2),Grb2(SH2!2,SH3) @'kp17'{0}, 'km17'",
    "egfr(Y1148~pY!1),Shc(PTB!1,Y317~pY) , Grb2(SH2,SH3!3),Sos(dom!3) <-> "
    "egfr(Y1148~pY!1),Shc(PTB!1,Y317~pY!2),Grb2(SH2!2,SH3!3),Sos(dom!3) @'kp24'{0}, 'km24'",
    "Shc(PTB!_,Y317~pY!2),Grb2(SH2!2,SH3) , Sos(dom) <-> "
    "Shc(PTB!_,Y317~pY!2),Grb2(SH2!2,SH3!3),Sos(dom!3) @'kp19'{0}, 'km19'",
    "Shc(PTB,Y317~pY) , Grb2(SH2,SH3) <-> Shc(PTB,Y317~pY!1),Grb2(SH2!1,SH3) @'kp21'{0}, 'km21'",
    "Shc(PTB,Y317~pY) , Grb2(SH2,SH3!_) <-> Shc(PTB,Y317~pY!1),Grb2(SH2!1,SH3!_) @'kp23'{0}, 'km23'",
    "Shc(PTB,Y317~pY) -> Shc(PTB,Y317~Y) @'km16'",
    "Grb2(SH2,SH3) , Sos(dom) <-> Grb2(SH2,SH3!1),Sos(dom!1) @'kp12'{0}, 'km12'",
    "Shc(PTB,Y317~pY!2),Grb2(SH2!2,SH3) , Sos(dom) <-> "
    "Shc(PTB,Y317~pY!2),Grb2(SH2!2,SH3!3),Sos(dom!3) @'kp22', 'km22'",
]
EGFR_SIGNATURES = ["egf(r)", "egfr(l,r,Y1068{Y pY},Y1148{Y pY})", "Grb2(SH2,SH3)", "Shc(PTB,Y317{Y pY})", "Sos(dom)"]


def kappa3_to_kappa4(expr: str, variables: dict | None = None) -> str:
    """Translate Kappa 3 agents (site~state!link, a bare site is free, ? is untested) and quoted rate names."""
    def agent(m):
        sites = []
        for tok in (t.strip() for t in m.group(2).split(",") if t.strip()):
            sm = re.fullmatch(r"([A-Za-z0-9_]+)(?:~([A-Za-z0-9_]+))?(\?|!(_|\d+|[A-Za-z0-9_]+\.[A-Za-z0-9_]+))?", tok)
            if not sm:
                raise ValueError(f"cannot translate Kappa 3 site {tok!r}")
            state = f"{{{sm.group(2)}}}" if sm.group(2) else ""
            link = "" if sm.group(3) == "?" else f"[{sm.group(4) if sm.group(4) else '.'}]"
            sites.append(sm.group(1) + state + link)
        return f"{m.group(1)}({','.join(sites)})"

    out = re.sub(r"([A-Za-z][A-Za-z0-9_]*)\(([^()]*)\)", agent, expr)
    if variables is not None:
        out = re.sub(r"'([^']+)'", lambda m: repr(float(variables[m.group(1)])), out)
    return out


def _egfr() -> dict:
    rules = []
    for i, text in enumerate(EGFR_KAPPA3_RULES):
        expr, rate = text.rsplit("@", 1)
        rules.append(f"'R{i + 1}' {kappa3_to_kappa4(expr)} @ {kappa3_to_kappa4(rate, EGFR_VARS)}")
    init = {kappa3_to_kappa4(e): EGFR_VARS[v] for v, e in EGFR_KAPPA3_INIT}
    return {"signatures": EGFR_SIGNATURES, "rules": rules, "init": init}


def model(p) -> dict:
    """The model's signatures, rules and init (custom: the parameters)."""
    if p.model != "custom":
        if p.signatures or p.rules or p.init:
            raise ValueError("signatures, rules and init are only used with model='custom'")
        return {"abc": ABC, "dimerization": DIMERIZATION, "polymerization": POLYMERIZATION}.get(p.model) or _egfr()
    if not isinstance(p.init, dict) or not p.init or not all(
            isinstance(k, str) and isinstance(v, (int, float)) and not isinstance(v, bool) and v >= 0
            for k, v in p.init.items()):
        raise ValueError("model custom needs init: a non-empty object mapping Kappa expressions to "
                         f"non-negative copy numbers, e.g. {{\"A(),B()\": 1000}}; got {p.init!r}")
    return {"signatures": p.signatures, "rules": p.rules, "init": p.init}


# ============================================================================
# Generator
# ============================================================================
def generate(p, rng):
    spec = model(p)
    sig = signatures(list(spec["signatures"]))
    rules = parse_rules(list(spec["rules"]), sig)
    registry = Registry(sig, rules, p.max_complex_size)
    init = initial_complexes(spec["init"], sig, registry)

    found: dict[tuple, list] = {}

    def react(*ids):
        m = len(ids)
        outcomes: dict[tuple, list] = {}
        for ri, rule in enumerate(rules):
            q = len(rule.components)
            if q == m:
                rate, assigns = rule.rate, sorted(set(permutations(range(m))))
            elif q == 2 and m == 1 and rule.uni:
                rate, assigns = rule.uni, [(0, 0)]
            else:
                continue
            if not rate:
                continue
            weight = rate * rule.omega(p.symmetry)
            for assign in assigns:
                lists = [registry.embeddings(ids[assign[c]], ri, c) for c in range(q)]
                if not all(lists):
                    continue
                for combo in product(*lists):
                    if m < q and set(combo[0]) & set(combo[1]):
                        continue
                    products = registry.apply(rule, ids, assign, combo)
                    if products is None:
                        continue
                    entry = outcomes.setdefault(tuple(sorted(products)), [0.0, []])
                    entry[0] += weight
                    if rule.label not in entry[1]:
                        entry[1].append(rule.label)
        if not outcomes:
            return None
        symmetry = prod(factorial(n) for n in Counter(ids).values())
        lhs = tuple(sorted(ids))
        for rhs, (k, labels) in outcomes.items():
            found[(lhs, rhs)] = [k / symmetry, labels]
        return list(outcomes)

    arities = sorted({len(r.components) for r in rules} | ({1} if any(r.uni for r in rules) else set()))
    species, pairs, status = expand(react, list(init), arity=arities, max_species=p.max_species,
                                    ordered=False, alternatives=True)
    if registry.oversized:
        status = "truncated"

    reactions, reaction_rules = [], []
    for lhs, rhs in pairs:
        k, labels = found[(tuple(sorted(lhs)), tuple(sorted(rhs)))]
        reactions.append(Reaction(dict(Counter(lhs)), dict(Counter(rhs)), {"law": "mass-action", "k": k}))
        reaction_rules.append(labels)

    known = set(species)
    initial = {c: n for c, n in init.items() if n is not None and c in known}
    complexes = registry.complexes
    conservation = [{"name": f"agents of type {t}",
                     "vector": {s: sum(1 for x in complexes[s].types if x == t) for s in species}}
                    for t in sig]
    sizes = [len(complexes[s].types) for s in species]
    return Network(
        species=[Species(s, structure=s) for s in species],
        reactions=reactions,
        status=status,
        initial_state=initial or None,
        extras={
            "model": p.model,
            "signatures": list(spec["signatures"]),
            "rules": list(spec["rules"]),
            "init": {e: n for e, n in spec["init"].items()},
            "kappa": kappa_file(spec),
            "rate_convention": {
                "symmetry": p.symmetry,
                "rule_symmetries": {r.label: {"lhs": r.aut_lhs, "preserved": r.aut_rule} for r in rules},
                "k": "k_R = sum over embeddings producing R of gamma * Omega / prod(m_i!), molecule-count units "
                     "(volume 1), so k_R prod x_i^m_i is the Kappa rule activity",
                "molecularity": "components of a rule's left-hand side match distinct complexes; inside one "
                                "complex only with a unimolecular rate {k} (KaDE)",
            },
            "reaction_rules": reaction_rules,
            "conservation": conservation,
            "species_encoding": "id = structure = canonical Kappa expression: agents in breadth-first order "
                                "from the root whose text is least, sites in signature order, bonds numbered "
                                "in order of appearance",
            "analysis": {
                "largest_complex": max(sizes) if sizes else 0,
                "truncated_by_max_complex_size": registry.oversized,
            },
        },
    )


def kappa_file(spec: dict) -> str:
    """The unflattened model as a Kappa file (signatures, rules, initial conditions)."""
    lines = [f"%agent: {s}" for s in spec["signatures"]]
    lines += list(spec["rules"])
    lines += [f"%init: {n:g} {e}" for e, n in spec["init"].items() if n is not None]
    return "\n".join(lines) + "\n"
