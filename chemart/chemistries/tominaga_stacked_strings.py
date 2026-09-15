"""Tominaga's artificial chemistry of stacked strings (pattern matching and recombination).

Catalog id: tominaga-stacked-strings. Book 18.3.2; Tominaga et al., Artificial
Life 13(3) 2007 [856] and 15(1) 2009 [855].

A molecule (object, v-molecule) is a stack of lines. Each line is a sequence of
elements and has a displacement relative to the first line; the string
notation writes each line as ``disp#elements/``:

    0#GGATGTAC/0#CCTACATGCCGA/        a DNA duplex with a sticky end
    0#FFFFFFFFF/-5#GGATG.../-5#CCTAC.../   Fok I bound to it

Elements are a capital letter followed by lower-case letters (``Orc``, ``Rp``,
``Coa``); the 2007 paper's subscripted element names (U01) are written with an
underscore (``U_01``).

Patterns use the same notation plus wildcards: a digit ``n`` matches one
element, ``*n`` (left end of a line only) and ``n*`` (right end only) match a
sequence of zero or more elements. A line pattern is L M R with L, R optional
sequence wildcards; the displacement of a pattern line is the offset of its
middle part M (the length of L counts as zero). An object pattern matches an
object with the same number of lines at a displacement d when every line
pattern i matches line i with M starting at position pdisp_i - odisp_i + d
(2007, section 2.2.3).

A recombination rule ``lhs -> rhs`` (terms joined by ``+``) consumes objects
matched by the left-hand patterns and produces the objects obtained by filling
the right-hand patterns with the matched elements; every wildcard of the lhs is
unique and appears once in the rhs. Sources supply an object without limit;
drains remove objects matched by a pattern. The dynamics is nondeterministic
(apply a rule, operate a source, or operate a drain, in any order).

method "closure" returns every reaction reachable from the initial pool and the
sources (each rule applied to every tuple of known species); method "soup"
samples one run of the nondeterministic process and returns the reactions that
fired. Sources are reactions ∅ -> s, drains are reactions s -> ∅.
"""

from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache
from itertools import product

from chemart.network import Network, Reaction, Species

SYSTEMS = (
    "benenson-automaton", "transcription", "fatty-acid-oxidation",
    "adleman-hamiltonian-path", "ab-concatenation", "custom",
)

_ELEMENT = re.compile(r"[A-Z][a-z]*(?:_[0-9a-z]+)?")
_LINE = re.compile(r"(-?\d+)#([^/#]*)/")

Molecule = tuple[tuple[int, tuple[str, ...]], ...]


# --- notation -------------------------------------------------------------------
def _elements(body: str, text: str) -> tuple[str, ...]:
    out, i = [], 0
    while i < len(body):
        if body[i] == " ":                    # a space only separates tokens (U_31 4*)
            i += 1
            continue
        m = _ELEMENT.match(body, i)
        if not m:
            raise ValueError(f"cannot read element at {body[i:]!r} in {text!r}")
        out.append(m.group())
        i = m.end()
    return tuple(out)


def _lines(text: str) -> list[tuple[int, str]]:
    text = text.strip()
    pos, out = 0, []
    for m in _LINE.finditer(text):
        if m.start() != pos:
            break
        out.append((int(m.group(1)), m.group(2)))
        pos = m.end()
    if not out or pos != len(text):
        raise ValueError(f"{text!r} is not in the notation disp#line/disp#line/...")
    return out


def normalize(lines) -> Molecule:
    """Shift displacements so that the first line has displacement 0."""
    lines = [(int(d), tuple(e)) for d, e in lines]
    s0 = lines[0][0]
    return tuple((d - s0, e) for d, e in lines)


@lru_cache(maxsize=None)
def parse_molecule(text: str) -> Molecule:
    """'0#AB/3#CD/' -> ((0, ('A', 'B')), (3, ('C', 'D')))."""
    return normalize([(d, _elements(body, text)) for d, body in _lines(text)])


def molecule_text(mol: Molecule) -> str:
    """Canonical string notation; it is the species id."""
    return "".join(f"{d}#{''.join(e)}/" for d, e in mol)


def layout(mol: Molecule) -> str:
    """The stacked layout, one text row per line, one column per element."""
    width = max((len(e) for _, els in mol for e in els), default=1)
    low = min(d for d, _ in mol)
    return "\n".join(
        (" " * (width * (d - low)) + "".join(e.ljust(width) for e in els)).rstrip()
        for d, els in mol
    )


def _pattern_line(body: str, text: str):
    toks, i, n = [], 0, len(body)
    while i < n:
        c = body[i]
        if c == " ":
            i += 1
        elif c == "*" and i + 1 < n and body[i + 1].isdigit():
            toks.append(("L", "*" + body[i + 1]))
            i += 2
        elif c.isdigit() and i + 1 < n and body[i + 1] == "*":
            toks.append(("R", body[i] + "*"))
            i += 2
        elif c.isdigit():
            toks.append(("E", c))
            i += 1
        else:
            m = _ELEMENT.match(body, i)
            if not m:
                raise ValueError(f"cannot read pattern symbol at {body[i:]!r} in {text!r}")
            toks.append(("lit", m.group()))
            i = m.end()
    left = toks.pop(0)[1] if toks and toks[0][0] == "L" else None
    right = toks.pop()[1] if toks and toks[-1][0] == "R" else None
    if any(kind in ("L", "R") for kind, _ in toks):
        raise ValueError(f"in {text!r}: *n may only open a line and n* may only close it")
    return left, tuple(toks), right


@lru_cache(maxsize=None)
def parse_pattern(text: str):
    """Object pattern: tuple of (pdisp, L, M, R); L/R are wildcard names or None."""
    return tuple((d, *_pattern_line(body, text)) for d, body in _lines(text))


def _wildcards(pattern) -> list[str]:
    out = []
    for _, left, mid, right in pattern:
        out += [left] if left else []
        out += [v for kind, v in mid if kind == "E"]
        out += [right] if right else []
    return out


def match(pattern, mol: Molecule) -> list[dict]:
    """All wildcard bindings under which the object pattern matches the object."""
    if len(pattern) != len(mol):
        return []
    base = mol[0][0] - pattern[0][0]
    found = []
    for d in range(base, base + len(mol[0][1]) + 1):
        binding = {}
        for (pdisp, left, mid, right), (odisp, line) in zip(pattern, mol):
            k, n, m = pdisp - odisp + d, len(line), len(mid)
            if k < 0 or k + m > n or (left is None and k) or (right is None and k + m != n):
                break
            if any(kind == "lit" and v != e for (kind, v), e in zip(mid, line[k:k + m])):
                break
            binding.update({v: (e,) for (kind, v), e in zip(mid, line[k:k + m]) if kind == "E"})
            if left:
                binding[left] = line[:k]
            if right:
                binding[right] = line[k + m:]
        else:
            if binding not in found:
                found.append(binding)
    return found


def build(pattern, binding: dict) -> Molecule:
    """Fill a right-hand-side pattern with bound elements."""
    lines = []
    for pdisp, left, mid, right in pattern:
        lpart = binding[left] if left else ()
        mpart = tuple(v if kind == "lit" else binding[v][0] for kind, v in mid)
        rpart = binding[right] if right else ()
        lines.append((pdisp - len(lpart), lpart + mpart + rpart))
    return normalize(lines)


class Rule:
    """A recombination rule, e.g. Rule('0#*1AB/ + 0#CD/ -> 0#*1AB/1#CD/', '(1)')."""

    def __init__(self, text: str, label: str = ""):
        self.text, self.label = " ".join(text.split()), label
        sides = re.split(r"->|→", text)
        if len(sides) != 2:
            raise ValueError(f"rule {text!r} needs exactly one '->'")
        self.lhs, self.rhs = ([parse_pattern(t.strip()) for t in side.split("+")] for side in sides)
        left = [w for p in self.lhs for w in _wildcards(p)]
        right = [w for p in self.rhs for w in _wildcards(p)]
        dup = [w for w, c in Counter(left).items() if c > 1]
        if dup:
            raise ValueError(f"rule {text!r}: wildcard {dup[0]} appears twice on the left-hand side")
        if Counter(left) != Counter(right):
            raise ValueError(f"rule {text!r}: every left-hand wildcard must appear exactly once "
                             f"on the right-hand side (left {sorted(left)}, right {sorted(right)})")

    @property
    def arity(self) -> int:
        return len(self.lhs)

    def outcomes(self, bindings) -> list[tuple[Molecule, ...]]:
        out = []
        for choice in product(*bindings):
            merged = {}
            for b in choice:
                merged.update(b)
            rhs = tuple(build(p, merged) for p in self.rhs)
            if rhs not in out:
                out.append(rhs)
        return out

    def apply(self, *mols) -> list[tuple[Molecule, ...]]:
        """Every product tuple of this rule on the given objects (text or parsed)."""
        mols = [parse_molecule(m) if isinstance(m, str) else m for m in mols]
        if len(mols) != self.arity:
            return []
        bindings = [match(p, m) for p, m in zip(self.lhs, mols)]
        return self.outcomes(bindings) if all(bindings) else []


# --- the published systems --------------------------------------------------------
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "X": "X"}

BENENSON_RULES = [
    ("(3)", "0#F0*/ + 0#GGATG1*/0#2*/ -> 0#F0*/-5#GGATG1*/-5#2*/"),
    ("(4)", "0#*0F/1#*123*/1#*456789*/ -> 0#*0F/1#*1/1#*45678/ + 0#23*/4#9*/"),
    ("(5)", "0#F0*/-5#GGATG1*/-5#2*/ -> 0#F0*/ + 0#GGATG1*/0#2*/"),
    ("(6)", "0#*0/0#*1CCGA/ + 0#GGCT2*/4#3*/ -> 0#*0GGCT2*/0#*1CCGA3*/"),
    ("(7)", "0#*0/0#*1GTCG/ + 0#CAGC2*/4#3*/ -> 0#*0CAGC2*/0#*1GTCG3*/"),
    ("(8)", "0#*0/0#*1GACC/ + 0#CTGG2*/4#3*/ -> 0#*0CTGG2*/0#*1GACC3*/"),
    ("(9)", "0#*0/0#*1GCGT/ + 0#CGCA2*/4#3*/ -> 0#*0CGCA2*/0#*1GCGT3*/"),
    ("(10)", "0#*0/0#*1AGCG/ + 0#TCGC2*/4#3*/ -> 0#*0TCGC2*/0#*1AGCG3*/"),
    ("(11)", "0#*0/0#*1ACAG/ + 0#TGTC2*/4#3*/ -> 0#*0TGTC2*/0#*1ACAG3*/"),
]
FOK_I = "0#FFFFFFFFF/"
TRANSITIONS = {
    ("S0", "a"): "0#GGATGTAC/0#CCTACATGCCGA/",
    ("S0", "b"): "0#GGATGACGAC/0#CCTACTGCTGGTCG/",
    ("S1", "a"): "0#GGATGTCG/0#CCTACAGCGACC/",
    ("S1", "b"): "0#GGATGG/0#CCTACCGCGT/",
}
DETECTOR_S0 = "0#X/0#XAGCG/"
DETECTOR_S1 = "0#X/0#XACAG/"
SYMBOL = {"a": "CTGGCT", "b": "CGCAGC"}
TERMINATOR = "TGTCGC"


def benenson_input(word: str, trail: int) -> str:
    """Input molecule: Fok I site, spacer, symbols, terminator, trailing tag of `trail` X."""
    upper = "GGATG" + "X" * 7 + "".join(SYMBOL[c] for c in word) + TERMINATOR + "X" * trail
    lower = "".join(_COMPLEMENT[c] for c in upper)
    return f"0#{upper}/0#{lower}/"


def benenson_reporter(state: str, trail: int) -> str:
    """The reporter left by an input with `trail` X ending in S0 (rule 10) or S1 (rule 11)."""
    tail = "X" * trail
    if state == "S0":
        return f"0#XTCGC{tail}/0#XAGCG{tail}/"
    return f"0#XTGTCGC{tail}/0#XACAGCG{tail}/"


def _benenson(p):
    words = p.words
    if not isinstance(words, list) or not words or not all(isinstance(w, str) and set(w) <= {"a", "b"} for w in words):
        raise ValueError(f"words must be a non-empty list of strings over 'a' and 'b', got {words!r}")
    pool = Counter({FOK_I: 100})
    for t in TRANSITIONS.values():
        pool[t] += 20
    pool[DETECTOR_S0] += 20
    if p.s1_detector:
        pool[DETECTOR_S1] += 20
    for i, w in enumerate(words):
        pool[benenson_input(w, i + 1)] += 10
    return BENENSON_RULES, pool, [], []


ADLEMAN_NODES = ["0#/0#L_01L_02/"] + [f"0#/-1#L_{i}1L_{i}2/" for i in range(1, 6)] + ["0#/-2#L_61L_62/"]
ADLEMAN_EDGES = {
    (0, 1): "0#U_01U_02U_11/2#/", (0, 3): "0#U_01U_02U_31/2#/", (0, 6): "0#U_01U_02U_61/2#/",
    (1, 2): "0#U_12U_21/1#/", (2, 1): "0#U_22U_11/1#/", (2, 3): "0#U_22U_31/1#/",
    (3, 2): "0#U_32U_21/1#/", (3, 4): "0#U_32U_41/1#/", (4, 1): "0#U_42U_11/1#/",
    (4, 5): "0#U_42U_51/1#/", (5, 1): "0#U_52U_11/1#/", (5, 2): "0#U_52U_21/1#/",
    (5, 6): "0#U_52U_61U_62/1#/",
}
ADLEMAN_RULES = (
    [(f"(1) U_{i}1", f"0#*1U_{i}1/0#*2/ + 0#3*/-1#L_{i}1 4*/ -> 0#*1U_{i}1 3*/0#*2L_{i}1 4*/")
     for i in range(1, 7)]
    + [(f"(2) L_{i}2", f"0#*1/0#*2L_{i}2/ + 0#U_{i}2 3*/1#4*/ -> 0#*1U_{i}2 3*/0#*2L_{i}2 4*/")
       for i in range(0, 6)]
    + [("start", "0#U_01U_02 1*/2#2*/ + 0#/0#L_01L_02/ -> 0#U_01U_02 1*/0#L_01L_02 2*/"),
       ("end", "0#*1U_61U_62/0#*2/ + 0#/-2#L_61L_62/ -> 0#*1U_61U_62/0#*2L_61L_62/")]
)
ADLEMAN_ANSWER = ("0#" + "".join(f"U_{i}1U_{i}2" for i in range(7)) + "/0#"
                  + "".join(f"L_{i}1L_{i}2" for i in range(7)) + "/")

TRANSCRIPTION_RULES = [
    ("(12)", "0#*1TATATT2*/0#*3ATATAA45*/ + 0#Rp/ -> 0#*1TATATT2*/0#*3ATATAA45*/6#Rp/"),
    ("(13)", "0#*12*/0#*34*/0#Rp/ + 0#Cap/ -> 0#*12*/0#*34*/-1#CapRp/"),
    ("(14)", "0#*12*/0#*3A4*/-1#*56Rp/ + 0#U/ -> 0#*12*/0#*3A4*/-1#*56URp/"),
    ("(15)", "0#*12*/0#*3T4*/-1#*56Rp/ + 0#A/ -> 0#*12*/0#*3T4*/-1#*56ARp/"),
    ("(16)", "0#*12*/0#*3G4*/-1#*56Rp/ + 0#C/ -> 0#*12*/0#*3G4*/-1#*56CRp/"),
    ("(17)", "0#*12*/0#*3C4*/-1#*56Rp/ + 0#G/ -> 0#*12*/0#*3C4*/-1#*56GRp/"),
    ("(18)", "0#*12*/0#*34*/-4#*5UUUURp/ -> 0#*12*/0#*34*/ + 0#*5UUUU/ + 0#Rp/"),
]
TRANSCRIPTION_SOURCES = ["0#U/", "0#C/", "0#A/", "0#G/", "0#Cap/"]
DNA = "TATATTCGCAATGCTGAGCTAGTTTT"


def chromosome(dna: str) -> str:
    return f"0#Orc{dna}/0#Orc{''.join(_COMPLEMENT[c] for c in dna)}/"


def _transcription(p):
    if not isinstance(p.dna, str) or not p.dna or set(p.dna) - set("TCAG"):
        raise ValueError(f"dna must be a non-empty string over T, C, A, G, got {p.dna!r}")
    return TRANSCRIPTION_RULES, Counter({"0#Rp/": 1, chromosome(p.dna): 1}), TRANSCRIPTION_SOURCES, []


FATTY_ACID_RULES = [
    ("(25)", "0#*1HHHO/0#*2CCCC/0#*3HHHSCoa/ + 0#Fad/ -> 0#*1HHHO/0#*2CCCC/0#*3HXXSCoa/ + 0#HFadH/"),
    ("(26)", "0#*1HHHO/0#*2CCCC/0#*3HXXSCoa/ + 0#HOH/ -> 0#H/-1#*1HOHO/-1#*2CCCC/-1#*3HHHSCoa/"),
    ("(27)", "0#H/-1#*1HOHO/-1#*2CCCC/-1#*3HHHSCoa/ + 0#NadPo/ "
             "-> 0#*1HOHO/0#*2CCCC/0#*3HXHSCoa/ + 0#NadH/ + 0#HPo/"),
    ("(28)", "0#*1HOHO/0#*2CCCC/0#*3HXHSCoa/ + 0#CoaSH/ -> 0#HO/-1#HCC/0#HSCoa/ + 0#*1HO/0#*2CC/0#*3HSCoa/"),
]
ACETYL_COA = "0#HO/-1#HCC/0#HSCoa/"
CARRIERS = ["0#Fad/", "0#NadPo/", "0#CoaSH/", "0#HOH/"]


def acyl_coa(carbons: int) -> str:
    """Fatty acyl CoA with an even number of carbons (2 gives acetyl CoA)."""
    return f"0#{'H' * (carbons - 1)}O/-1#H{'C' * carbons}/0#{'H' * (carbons - 1)}SCoa/"


def _fatty_acid(p):
    if p.carbons % 2:
        raise ValueError(f"carbons must be even (the pathway oxidises even-chain acyl CoA), got {p.carbons}")
    pool = Counter({acyl_coa(p.carbons): 1})
    for c in CARRIERS:
        pool[c] += 1
    return FATTY_ACID_RULES, pool, [], []


AB_RULES = [
    ("(1)", "0#*1AB/ + 0#CD/ -> 0#*1AB/1#CD/"),
    ("(2)", "0#*1AB/1#CD/ + 0#AB2*/ -> 0#*1ABAB2*/1#CD/"),
    ("(3)", "0#*1ABAB2*/1#CD/ -> 0#*1ABAB2*/ + 0#CD/"),
]


def _custom(p):
    if not isinstance(p.rules, list) or not p.rules or not all(isinstance(r, str) for r in p.rules):
        raise ValueError("system custom needs rules: a non-empty list of rule strings 'lhs -> rhs'")
    if not isinstance(p.pool, dict) or not all(
            isinstance(k, str) and isinstance(v, int) and not isinstance(v, bool) and v >= 1
            for k, v in p.pool.items()):
        raise ValueError(f"pool must map molecule strings to positive integer counts, got {p.pool!r}")
    for name in ("sources", "drains"):
        value = getattr(p, name)
        if not isinstance(value, list) or not all(isinstance(s, str) for s in value):
            raise ValueError(f"{name} must be a list of strings, got {value!r}")
    if not p.pool and not p.sources:
        raise ValueError("system custom needs a non-empty pool or at least one source")
    rules = [(f"r{i + 1}", r) for i, r in enumerate(p.rules)]
    return rules, Counter(p.pool), list(p.sources), list(p.drains)


def system(p):
    """(rules [(label, text)], initial pool Counter of texts, sources, drains)."""
    if p.system != "custom" and (p.rules or p.pool or p.sources or p.drains):
        raise ValueError("rules, pool, sources and drains are only used with system='custom'")
    if p.system == "benenson-automaton":
        return _benenson(p)
    if p.system == "transcription":
        return _transcription(p)
    if p.system == "fatty-acid-oxidation":
        return _fatty_acid(p)
    if p.system == "adleman-hamiltonian-path":
        pool = Counter(ADLEMAN_NODES) + Counter(ADLEMAN_EDGES.values())
        return ADLEMAN_RULES, pool, [], []
    if p.system == "ab-concatenation":
        return AB_RULES, Counter({"0#CD/": 1}), ["0#AB/"], ["0#ABABABABAB1*/"]
    return _custom(p)


# --- closure and soup --------------------------------------------------------------
def closure(rules: list[Rule], seed: list[Molecule], max_species: int = 2000):
    """Level-saturation closure: every rule on every tuple with a newly found species.

    Returns (species, reactions [(lhs, rhs, rule label)], status).
    """
    species = list(dict.fromkeys(seed))
    truncated = len(species) > max_species
    species = species[:max_species]
    known = set(species)
    tables = [[[] for _ in r.lhs] for r in rules]      # per rule, per position: (index, bindings)
    reactions: dict[tuple, tuple] = {}
    done = 0
    while done < len(species):
        old, n = done, len(species)
        done = n
        for r, table in zip(rules, tables):
            for pos, pattern in enumerate(r.lhs):
                for i in range(old, n):
                    bindings = match(pattern, species[i])
                    if bindings:
                        table[pos].append((i, bindings))
        for r, table in zip(rules, tables):
            for first in range(r.arity):
                lists = ([[x for x in table[q] if x[0] < old] for q in range(first)]
                         + [[x for x in table[first] if x[0] >= old]]
                         + [table[q] for q in range(first + 1, r.arity)])
                for combo in product(*lists):
                    lhs = tuple(species[i] for i, _ in combo)
                    for rhs in r.outcomes([b for _, b in combo]):
                        if Counter(lhs) == Counter(rhs):
                            continue
                        novel = [m for m in dict.fromkeys(rhs) if m not in known]
                        if len(species) + len(novel) > max_species:
                            truncated = True
                            continue
                        species.extend(novel)
                        known.update(novel)
                        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
                        reactions.setdefault(key, (lhs, rhs, r.label))
    return species, list(reactions.values()), "truncated" if truncated else "complete"


def run(rules: list[Rule], pool: Counter, sources: list[Molecule], drains: list, steps: int, rng):
    """One sampled run of the nondeterministic process.

    Every step picks an event with probability proportional to its weight: a
    rule by the number of ordered reactant choices (product over its patterns
    of the copies matching that pattern), a source with weight 1, a drain by the
    number of copies it matches. Reactants are drawn in proportion to their
    copies and a match (displacement) uniformly; a draw that needs more copies
    than exist does nothing.
    """
    counts = Counter({m: c for m, c in pool.items() if c > 0})
    cache: dict[Molecule, list] = {}
    fired: dict[tuple, list] = {}
    seen = dict.fromkeys(list(counts) + list(sources))

    def matches(mol):
        if mol not in cache:
            cache[mol] = [[match(p, mol) for p in r.lhs] for r in rules]
            cache[mol].append([bool(match(d, mol)) for d in drains])
        return cache[mol]

    def record(lhs, rhs, label):
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        fired.setdefault(key, [lhs, rhs, label, 0])[3] += 1

    for _ in range(steps):
        present = [m for m, c in counts.items() if c > 0]
        events = []
        for j, r in enumerate(rules):
            per_pos = [[(m, counts[m]) for m in present if matches(m)[j][q]] for q in range(r.arity)]
            weight = 1
            for choices in per_pos:
                weight *= sum(c for _, c in choices)
            if weight:
                events.append((weight, "rule", j, per_pos))
        for s in sources:
            events.append((1, "source", s, None))
        for k in range(len(drains)):
            hit = [(m, counts[m]) for m in present if matches(m)[-1][k]]
            if hit:
                events.append((sum(c for _, c in hit), "drain", k, hit))
        total = sum(e[0] for e in events)
        if not total:
            break
        x = rng.random() * total
        for weight, kind, what, data in events:
            if x < weight:
                break
            x -= weight
        if kind == "source":
            counts[what] += 1
            record((), (what,), "source")
            continue
        if kind == "drain":
            mol = _pick(data, rng)
            counts[mol] -= 1
            record((mol,), (), "drain")
            continue
        r = rules[what]
        lhs = tuple(_pick(choices, rng) for choices in data)
        if any(counts[m] < c for m, c in Counter(lhs).items()):
            continue
        bindings = [matches(m)[what][q] for q, m in enumerate(lhs)]
        choice = [b[int(rng.integers(len(b)))] for b in bindings]
        rhs = r.outcomes([[b] for b in choice])[0]
        if Counter(lhs) == Counter(rhs):
            continue
        counts.subtract(Counter(lhs))
        counts.update(Counter(rhs))
        seen.update(dict.fromkeys(rhs))
        record(lhs, rhs, r.label)
    return list(fired.values()), +counts, list(seen)


def _pick(choices, rng):
    x = rng.random() * sum(c for _, c in choices)
    for mol, c in choices:
        if x < c:
            return mol
        x -= c
    return choices[-1][0]


# --- generator ----------------------------------------------------------------------
def _species(mols) -> list[Species]:
    return [Species(molecule_text(m), structure=layout(m)) for m in mols]


def _ids(mols) -> list[str]:
    return [molecule_text(m) for m in mols]


def generate(p, rng):
    texts, pool_texts, source_texts, drain_texts = system(p)
    rules = [Rule(text, label) for label, text in texts]
    pool = Counter()
    for text, c in pool_texts.items():
        pool[parse_molecule(text)] += c
    sources = [parse_molecule(s) for s in source_texts]
    drains = [parse_pattern(d) for d in drain_texts]

    extras = {
        "system": p.system,
        "rules": [f"{label}: {r.text}" for label, r in zip((t[0] for t in texts), rules)],
        "sources": _ids(sources),
        "drains": list(drain_texts),
    }
    if p.method == "closure":
        seed = list(pool) + sources
        species, found, status = closure(rules, seed, p.max_species)
        reactions = [Reaction.of(_ids(lhs), _ids(rhs)) for lhs, rhs, _ in found]
        labels = [label for _, _, label in found]
        for s in sources:
            reactions.append(Reaction({}, {molecule_text(s): 1}))
            labels.append("source")
        for m in species:
            if any(match(d, m) for d in drains):
                reactions.append(Reaction({molecule_text(m): 1}, {}))
                labels.append("drain")
        initial = {molecule_text(m): float(c) for m, c in pool.items()}
    else:
        if p.copies > 1:
            pool = Counter({m: c * p.copies for m, c in pool.items()})
        fired, final, species = run(rules, pool, sources, drains, p.steps, rng)
        reactions = [Reaction.of(_ids(lhs), _ids(rhs), count=c) for lhs, rhs, _, c in fired]
        labels = [label for _, _, label, _ in fired]
        status = "observed"
        initial = {molecule_text(m): float(c) for m, c in pool.items()}
        extras["final_state"] = {molecule_text(m): c for m, c in final.most_common()}
    extras["reaction_rules"] = labels
    extras["analysis"] = _analysis(p, {molecule_text(m) for m in species}, extras.get("final_state"))
    return Network(species=_species(species), reactions=reactions, status=status,
                   initial_state=initial, extras=extras)


def _analysis(p, ids: set[str], final: dict | None) -> dict:
    present = (lambda s: final.get(s, 0) > 0) if final is not None else (lambda s: s in ids)
    if p.system == "benenson-automaton":
        states = {}
        for i, w in enumerate(p.words):
            got = [s for s in ("S0", "S1") if present(benenson_reporter(s, i + 1))]
            states[w] = got
        return {"reporters": states,
                "accepted": [w for w, s in states.items() if "S0" in s]}
    if p.system == "fatty-acid-oxidation":
        return {"acyl_coa_carbons": [c for c in range(2, p.carbons + 1, 2) if present(acyl_coa(c))],
                "acetyl_coa": present(ACETYL_COA)}
    if p.system == "adleman-hamiltonian-path":
        return {"answer": ADLEMAN_ANSWER, "answer_found": present(ADLEMAN_ANSWER)}
    if p.system == "transcription":
        mrna = "0#Cap" + p.dna.replace("T", "U")[p.dna.find("TATATT") + 6:] + "/" if "TATATT" in p.dna else None
        return {"mrna": mrna, "mrna_found": bool(mrna) and present(mrna)}
    return {}
