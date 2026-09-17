"""Algorithmic music composition with an artificial chemistry.

Catalog id: music-ac. Book 16.4; Miura & Tominaga, GWAL-7 (2006) [592].

Music composition is written as a *nondeterministic algorithm* and executed by
Tominaga's artificial chemistry of pattern matching and recombination: every
choice point of the composition (which note, which chord, which bar comes next)
is a choice of colliding molecules in a well-stirred reactor.

Molecules (objects) are stacks of lines of elements, written ``disp#elements/``
per line, the displacement counting columns relative to the first line:

    0#C4D4E4/                                a melody fragment (3 eighth notes)
    0#C4D4E4F4G4A4B4C5/0#TCCCCCCC/           one bar under the tonic chord C
    0#StartE4D4.../1#TCCCCCCCDGGGGGGG.../    a finished four-bar phrase

Elements are the notes C4..C5, the triads C, Dm, Em, F, G, Am, Bm (Bm is the
paper's Bm-diminished), the chord functions T, S, D, D2 and the control
elements Degree, Chord, Avoid, Start, Stop and Dummy.

Patterns add wildcards to that notation: ``<1>`` matches one element, ``<*1>``
(left end of a line) and ``<1*>`` (right end) match a sequence of zero or more
elements, and ``<1..7>`` abbreviates ``<1><2>...<7>``. A recombination rule
``lhs -> rhs`` consumes the objects matched by its left-hand patterns and
produces the right-hand patterns filled with the matched elements; every
wildcard occurs once on each side, so elements are conserved exactly.

The 65 rules built here implement the five steps of the paper's algorithm
(section 3): grow a conjunct melody, cut off a bar and give it a chord that
contains its first note, replace avoid notes, concatenate bars along the
cadence rules, and cut out four-bar phrases that start on a tonic bar.

The reactor is a random-collision run (section 5); the observed network holds
the reactions that fired with their firing counts, and extras["analysis"]
carries the phrases that were composed, as notes, chords and functions.
"""

from __future__ import annotations

import re
from collections import Counter
from itertools import product

from chemart.network import Network, Reaction, Species

# --- the musical material (paper, section 3) -------------------------------------
#: the diatonic scale of C major in one octave, as eighth notes (section 3)
NOTES = ("C4", "D4", "E4", "F4", "G4", "A4", "B4", "C5")
#: pitch class of each note element
PITCH = {n: n[0] for n in NOTES}
#: the triads on the degrees of C major; Bm is the paper's Bm-diminished
CHORD_TONES = {
    "C": ("C", "E", "G"), "Dm": ("D", "F", "A"), "Em": ("E", "G", "B"),
    "F": ("F", "A", "C"), "G": ("G", "B", "D"), "Am": ("A", "C", "E"),
    "Bm": ("B", "D", "F"),
}
#: chords of each functional group (section 3.4)
FUNCTION_CHORDS = {
    "T": ("C", "Em", "Am"), "S": ("Dm", "F", "Am"),
    "D": ("Em", "G", "Bm"), "D2": ("Dm", "F"),
}
#: avoid pitch of each chord (section 3.3); F has none, "C" means C4 and C5
AVOID_PITCHES = {
    "C": ("F4",), "Dm": ("B4",), "Em": ("F4",), "G": ("C4", "C5"),
    "Am": ("F4",), "Bm": ("C4", "C5"),
}
#: cadence rules as implemented in section 4.3 ("T-D2 ... T-D, T-S, D-T, D2-D, D2-T")
CADENCES = (("T", "D2"), ("T", "D"), ("T", "S"), ("D", "T"), ("D2", "D"), ("D2", "T"))
#: the transition figure 3(1) needs and the typical cadence T-S-T implies
SUBDOMINANT_TO_TONIC = ("S", "T")
CADENCE_SETS = ("published", "with-s-t")
MODES = ("homophonic",)

CONTROL = ("Degree", "Chord", "Avoid", "Start", "Stop")

_ELEMENT = re.compile(r"[A-Z][a-z0-9]*")
_LINE = re.compile(r"(-?\d+)#([^/#]*)/")
_TOKEN = re.compile(r"<(\*?\d+\*?|\d+\.\.\d+)>")

Object = tuple[tuple[int, tuple[str, ...]], ...]


# --- notation ---------------------------------------------------------------------
def _split_lines(text: str) -> list[tuple[int, str]]:
    text = text.strip()
    pos, out = 0, []
    for m in _LINE.finditer(text):
        if m.start() != pos:
            break
        out.append((int(m.group(1)), m.group(2)))
        pos = m.end()
    if not out or pos != len(text):
        raise ValueError(f"{text!r} is not in the notation disp#elements/disp#elements/...")
    return out


def _elements(body: str, text: str) -> tuple[str, ...]:
    out, i = [], 0
    while i < len(body):
        if body[i] == " ":
            i += 1
            continue
        m = _ELEMENT.match(body, i)
        if not m:
            raise ValueError(f"cannot read an element at {body[i:]!r} in {text!r}")
        out.append(m.group())
        i = m.end()
    return tuple(out)


def normalize(lines) -> Object:
    """Shift displacements so that the first line sits at 0."""
    lines = [(int(d), tuple(e)) for d, e in lines]
    base = lines[0][0]
    return tuple((d - base, e) for d, e in lines)


def parse_object(text: str) -> Object:
    """'0#C4D4/0#TCC/' -> ((0, ('C4', 'D4')), (0, ('T', 'C', 'C')))."""
    return normalize([(d, _elements(body, text)) for d, body in _split_lines(text)])


def object_text(obj: Object) -> str:
    """Canonical string notation; this is the species id."""
    return "".join(f"{d}#{''.join(e)}/" for d, e in obj)


def _pattern_line(body: str, text: str):
    """(left, mid, right): sequence wildcards at the ends, middle as (kind, value)."""
    toks, i, n = [], 0, len(body)
    while i < n:
        if body[i] == " ":
            i += 1
            continue
        m = _TOKEN.match(body, i)
        if m:
            name = m.group(1)
            if ".." in name:
                lo, hi = (int(v) for v in name.split(".."))
                toks += [("E", str(k)) for k in range(lo, hi + 1)]
            elif name.startswith("*"):
                toks.append(("L", name))
            elif name.endswith("*"):
                toks.append(("R", name))
            else:
                toks.append(("E", name))
            i = m.end()
            continue
        m = _ELEMENT.match(body, i)
        if not m:
            raise ValueError(f"cannot read a pattern symbol at {body[i:]!r} in {text!r}")
        toks.append(("lit", m.group()))
        i = m.end()
    left = toks.pop(0)[1] if toks and toks[0][0] == "L" else None
    right = toks.pop()[1] if toks and toks[-1][0] == "R" else None
    if any(kind in ("L", "R") for kind, _ in toks):
        raise ValueError(f"in {text!r}: <*n> may only open a line and <n*> may only close it")
    return left, tuple(toks), right


def parse_pattern(text: str):
    """Object pattern: tuple of (displacement, left, middle, right) per line."""
    return tuple((d, *_pattern_line(body, text)) for d, body in _split_lines(text))


def _wildcards(pattern) -> list[str]:
    out = []
    for _, left, mid, right in pattern:
        out += [left] if left else []
        out += [v for kind, v in mid if kind == "E"]
        out += [right] if right else []
    return out


def match(pattern, obj: Object) -> list[dict]:
    """Every binding under which the pattern matches the object.

    A pattern line's displacement locates its middle part, the length of a
    sequence wildcard counting as zero (paper, section 2). The displacement d
    of the whole match is fixed by the first line unless that line has sequence
    wildcards on both ends.
    """
    if len(pattern) != len(obj):
        return []
    pdisp, left, mid, right = pattern[0]
    odisp, line = obj[0]
    n, m = len(line), len(mid)
    if m > n:
        return []
    if left is None and right is None:
        ds = [odisp - pdisp] if m == n else []
    elif left is None:
        ds = [odisp - pdisp]
    elif right is None:
        ds = [n - m + odisp - pdisp]
    else:
        ds = [k + odisp - pdisp for k in range(n - m + 1)]

    found = []
    for d in ds:
        binding: dict[str, tuple[str, ...]] = {}
        for (pd, lft, mdl, rgt), (od, ln) in zip(pattern, obj):
            k, size, width = pd - od + d, len(ln), len(mdl)
            if k < 0 or k + width > size or (lft is None and k) or (rgt is None and k + width != size):
                break
            if any(kind == "lit" and v != e for (kind, v), e in zip(mdl, ln[k:k + width])):
                break
            binding.update({v: (e,) for (kind, v), e in zip(mdl, ln[k:k + width]) if kind == "E"})
            if lft:
                binding[lft] = ln[:k]
            if rgt:
                binding[rgt] = ln[k + width:]
        else:
            if binding not in found:
                found.append(binding)
    return found


def build(pattern, binding: dict) -> Object:
    """Fill a right-hand-side pattern with the bound elements."""
    lines = []
    for pdisp, left, mid, right in pattern:
        lpart = binding[left] if left else ()
        mpart = tuple(v if kind == "lit" else binding[v][0] for kind, v in mid)
        rpart = binding[right] if right else ()
        lines.append((pdisp - len(lpart), lpart + mpart + rpart))
    return normalize(lines)


class Rule:
    """A recombination rule, e.g. Rule('0#<*1>D4/ + 0#Degree/0#D4<2>/ -> ...')."""

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
            raise ValueError(f"rule {text!r}: wildcard <{dup[0]}> appears twice on the left-hand side")
        if Counter(left) != Counter(right):
            raise ValueError(f"rule {text!r}: each wildcard must appear once on each side "
                             f"(left {sorted(left)}, right {sorted(right)})")

    @property
    def arity(self) -> int:
        return len(self.lhs)

    def outcomes(self, bindings) -> list[tuple[Object, ...]]:
        out = []
        for choice in product(*bindings):
            merged: dict = {}
            for b in choice:
                merged.update(b)
            rhs = tuple(build(p, merged) for p in self.rhs)
            if rhs not in out:
                out.append(rhs)
        return out

    def apply(self, *objects) -> list[tuple[Object, ...]]:
        """Every product tuple of this rule on the given objects (text or parsed)."""
        objs = [parse_object(o) if isinstance(o, str) else o for o in objects]
        if len(objs) != self.arity:
            return []
        bindings = [match(p, o) for p, o in zip(self.lhs, objs)]
        return self.outcomes(bindings) if all(bindings) else []


# --- the published system (paper, section 4) ---------------------------------------
def _seq(lo: int, hi: int) -> str:
    """Wildcard run <lo..hi>, or <lo> for a single one, or '' when empty."""
    if hi < lo:
        return ""
    return f"<{lo}>" if hi == lo else f"<{lo}..{hi}>"


def chord_line(function: str, chord: str, notes_per_bar: int) -> str:
    """Second line of a bar: the function element and then the chord, one per note."""
    return function + chord * (notes_per_bar - 1)


def elements(notes=NOTES) -> list[str]:
    """The element kinds of the description (section 4.1), plus Dummy."""
    return [*notes, *CHORD_TONES, *FUNCTION_CHORDS, *CONTROL, "Dummy"]


def initial_pool(notes=NOTES, notes_per_bar: int = 8) -> Counter:
    """The initial working multiset of section 4.2: 43 kinds, 1826 objects."""
    pool: Counter = Counter()
    for a, b in zip(notes, notes[1:]):                       # conjunct motion [100 each]
        pool[f"0#Degree/0#{a}{b}/"] = 100
        pool[f"0#Degree/0#{b}{a}/"] = 100
    for function, chords in FUNCTION_CHORDS.items():          # chord objects [20 each]
        for chord in chords:
            pool[f"0#Chord/0#{chord_line(function, chord, notes_per_bar)}/"] = 20
    pool["0#Avoid/0#Dummy/"] = 10                             # replaces an avoid note
    for note in notes:                                        # random notes for Dummy [10 each]
        pool[f"0#Avoid/0#{note}/"] = 10
    pool["0#StartStop/"] = 100                                # terminates a phrase
    for note in notes:                                        # seeds [2 each]
        pool[f"0#{note}/"] = 2
    return pool


def rule_set(notes=NOTES, notes_per_bar: int = 8, bars_per_phrase: int = 4,
             avoid_notes: bool = True, cadences: str = "published") -> list[Rule]:
    """The 65 recombination rules of section 4.3, in the order of the algorithm."""
    n, rules = notes_per_bar, []
    cadence_pairs = list(CADENCES)
    if cadences == "with-s-t":
        cadence_pairs.append(SUBDOMINANT_TO_TONIC)

    # Step 1, rule (1): add a note to a sequence, using a conjunct-motion object.
    for note in notes:
        rules.append(Rule(
            f"0#<*1>{note}/ + 0#Degree/0#{note}<2>/ -> 0#<*1>{note}<2>/ + 0#Degree/0#{note}/",
            f"(1) note after {note}"))

    # Step 2, rules (2)-(5): cut off a bar and give it a chord containing its first note.
    mid = _seq(1, n - 1)
    for function, chords in FUNCTION_CHORDS.items():
        for chord in chords:
            line = chord_line(function, chord, n)
            for note in notes:
                if PITCH[note] not in CHORD_TONES[chord]:
                    continue
                rules.append(Rule(
                    f"0#{note}{mid}<{n}><{n + 1}*>/ + 0#Chord/0#{line}/ "
                    f"-> 0#{note}{mid}/0#{line}/ + 0#<{n}><{n + 1}*>/ + 0#Chord/",
                    f"(2) {note} -> {function} {chord}"))

    # Step 3, rules (6)-(7): replace an avoid note, via Dummy.
    if avoid_notes:
        for chord, pitches in AVOID_PITCHES.items():
            for pitch in pitches:
                rules.append(Rule(
                    f"0#<*1>{pitch}<2*>/0#<*3>{chord}<4*>/ + 0#Avoid/0#Dummy/ "
                    f"-> 0#<*1>Dummy<2*>/0#<*3>{chord}<4*>/ + 0#Avoid/0#{pitch}/",
                    f"(6) avoid {pitch} of {chord}"))
        rules.append(Rule(
            "0#<*1>Dummy<2*>/0#<*3><4><5*>/ + 0#Avoid/0#<6>/ "
            "-> 0#<*1><6><2*>/0#<*3><4><5*>/ + 0#Avoid/0#Dummy/",
            "(7) Dummy -> a random note"))

    # Step 4, rule (8): concatenate two bar sequences along a cadence.
    # Wildcard numbers follow the paper's rule (8) for the default eight notes
    # per bar: <*1><2..9>/<*10>T<11..17>/ + <18><19><20*>/D2<21*>/.
    last = _seq(2, n + 1)                       # the notes of the last bar of the left object
    last_chords = _seq(n + 3, 2 * n + 1)        # its chords, after the function element
    right = 2 * n + 2                           # first free number for the right-hand object
    for first, second in cadence_pairs:
        rules.append(Rule(
            f"0#<*1>{last}/0#<*{n + 2}>{first}{last_chords}/ "
            f"+ 0#<{right}><{right + 1}><{right + 2}*>/0#{second}<{right + 3}*>/ "
            f"-> 0#<*1>{last}<{right}><{right + 1}><{right + 2}*>/"
            f"0#<*{n + 2}>{first}{last_chords}{second}<{right + 3}*>/",
            f"(8) {first}-{second}"))

    # Step 5, rule (9): cut out a phrase of bars_per_phrase bars starting on a tonic bar.
    # Paper: 0#<1..32><33><34*>/0#T<35..65><66><67*>/ + 0#StartStop/ -> ...
    span = n * bars_per_phrase
    head, tail = _seq(1, span), _seq(span + 3, 2 * span + 1)
    rules.append(Rule(
        f"0#{head}<{span + 1}><{span + 2}*>/0#T{tail}<{2 * span + 2}><{2 * span + 3}*>/ "
        f"+ 0#StartStop/ -> 0#Start{head}Stop/1#T{tail}/ "
        f"+ 0#<{span + 1}><{span + 2}*>/0#<{2 * span + 2}><{2 * span + 3}*>/",
        f"(9) {bars_per_phrase}-bar phrase"))

    # Step 5, rule (10): detach a leading non-tonic bar, which is reused.
    # Paper: 0#<1..8><9><10*>/0#D<11..17><18><19*>/ -> ...
    bar, chords = _seq(1, n), _seq(n + 3, 2 * n + 1)
    for function in ("D", "D2", "S"):
        rules.append(Rule(
            f"0#{bar}<{n + 1}><{n + 2}*>/0#{function}{chords}<{2 * n + 2}><{2 * n + 3}*>/ "
            f"-> 0#{bar}/0#{function}{chords}/ "
            f"+ 0#<{n + 1}><{n + 2}*>/0#<{2 * n + 2}><{2 * n + 3}*>/",
            f"(10) detach leading {function} bar"))
    return rules


# --- reading the music off an object ------------------------------------------------
def is_phrase(obj: Object) -> bool:
    """A finished phrase: terminated by Start and Stop, inert under every rule."""
    line = obj[0][1]
    return len(line) > 2 and line[0] == "Start" and line[-1] == "Stop"


def read_bars(obj: Object, notes_per_bar: int = 8) -> list[dict]:
    """Bars of a bar sequence or phrase: the notes, the chord and its function."""
    if len(obj) != 2:
        return []
    line, chords = list(obj[0][1]), list(obj[1][1])
    if is_phrase(obj):
        line = line[1:-1]
    if not chords or chords[0] not in FUNCTION_CHORDS:
        return []
    bars = []
    for start in range(0, len(chords), notes_per_bar):
        block = chords[start:start + notes_per_bar]
        tones = block[1:]
        if block[0] not in FUNCTION_CHORDS or not tones or len(set(tones)) != 1:
            return []
        bars.append({"function": block[0], "chord": tones[0],
                     "notes": line[start:start + notes_per_bar]})
    return bars


def describe(obj: Object, notes_per_bar: int = 8) -> str:
    """The musical content of an object, as the species structure."""
    bars = read_bars(obj, notes_per_bar)
    if bars:
        kind = "phrase" if is_phrase(obj) else "bars"
        return kind + ": " + " | ".join(
            f"{b['chord']}({b['function']}) " + " ".join(b["notes"]) for b in bars)
    if len(obj) == 1 and all(e in PITCH for e in obj[0][1]):
        return "notes: " + " ".join(obj[0][1])
    return object_text(obj)


# --- the reactor (paper, section 5: random collision) --------------------------------
class Reactor:
    """A well-stirred multiset of objects reacting by random collision.

    Each step draws one rule with a weight equal to its number of ordered
    reactant choices, then draws the reactants in proportion to their copies and
    one match of each uniformly, so that every step is a collision that can
    react. A draw that would need more copies than the multiset holds fails.
    """

    def __init__(self, rules: list[Rule], pool: Counter, rng):
        self.rules, self.rng = rules, rng
        self.counts: Counter = Counter({t: c for t, c in pool.items() if c > 0})
        self.parsed = {t: parse_object(t) for t in self.counts}
        self.slots: dict[str, list[tuple[int, int, list]]] = {}
        self.totals = [[0] * r.arity for r in rules]
        self.fired: dict[tuple, list] = {}
        self.seen = dict.fromkeys(self.counts)
        self.failed = 0
        for text, count in self.counts.items():
            self._index(text, count)

    def _index(self, text: str, delta: int) -> None:
        if text not in self.slots:
            obj = self.parsed.setdefault(text, parse_object(text))
            hits = []
            for i, rule in enumerate(self.rules):
                for pos, pattern in enumerate(rule.lhs):
                    bindings = match(pattern, obj)
                    if bindings:
                        hits.append((i, pos, bindings))
            self.slots[text] = hits
        for i, pos, _ in self.slots[text]:
            self.totals[i][pos] += delta

    def _change(self, text: str, delta: int) -> None:
        self.counts[text] += delta
        self._index(text, delta)

    def _pick(self, rule_index: int, pos: int) -> str:
        total = self.totals[rule_index][pos]
        x = self.rng.random() * total
        for text, count in self.counts.items():
            if count and any(i == rule_index and q == pos for i, q, _ in self.slots[text]):
                if x < count:
                    return text
                x -= count
        raise RuntimeError("no reactant found")           # pragma: no cover

    def step(self) -> bool:
        weights = []
        for i, rule in enumerate(self.rules):
            w = 1
            for pos in range(rule.arity):
                w *= self.totals[i][pos]
            weights.append(w)
        total = sum(weights)
        if not total:
            return False
        x = self.rng.random() * total
        index = len(weights) - 1
        for i, w in enumerate(weights):
            if x < w:
                index = i
                break
            x -= w
        rule = self.rules[index]
        lhs = tuple(self._pick(index, pos) for pos in range(rule.arity))
        if any(self.counts[t] < c for t, c in Counter(lhs).items()):
            self.failed += 1
            return True
        bindings = []
        for pos, text in enumerate(lhs):
            options = next(b for i, q, b in self.slots[text] if i == index and q == pos)
            bindings.append([options[int(self.rng.integers(len(options)))]])
        rhs_objects = rule.outcomes(bindings)[0]
        rhs = tuple(object_text(o) for o in rhs_objects)
        if Counter(lhs) == Counter(rhs):
            self.failed += 1
            return True
        for text in lhs:
            self._change(text, -1)
        for text, obj in zip(rhs, rhs_objects):
            self.parsed.setdefault(text, obj)
            self._change(text, 1)
            self.seen.setdefault(text, None)
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        self.fired.setdefault(key, [lhs, rhs, rule.label, 0])[3] += 1
        return True

    def phrases(self) -> list[str]:
        return [t for t, c in self.counts.items() if c > 0 and is_phrase(self.parsed[t])]

    def run(self, steps: int, stop_after: int = 0) -> None:
        for _ in range(steps):
            if stop_after and len(self.phrases()) >= stop_after:
                return
            if not self.step():
                return


# --- generator ------------------------------------------------------------------------
def generate(p, rng) -> Network:
    if p.mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {p.mode!r}")
    if p.cadences not in CADENCE_SETS:
        raise ValueError(f"cadences must be one of {CADENCE_SETS}, got {p.cadences!r}")
    if p.notes_per_bar < 2:
        raise ValueError(f"notes_per_bar must be at least 2, got {p.notes_per_bar}")

    rules = rule_set(NOTES, p.notes_per_bar, p.bars_per_phrase, p.avoid_notes, p.cadences)
    pool = initial_pool(NOTES, p.notes_per_bar)
    if p.copies != 1:
        pool = Counter({t: c * p.copies for t, c in pool.items()})

    reactor = Reactor(rules, pool, rng)
    reactor.run(p.steps, p.phrases)

    reactions = [Reaction.of(lhs, rhs, count=count) for lhs, rhs, _, count in reactor.fired.values()]
    labels = [label for _, _, label, _ in reactor.fired.values()]
    species_ids = list(reactor.seen)
    species = [Species(t, structure=describe(reactor.parsed[t], p.notes_per_bar))
               for t in species_ids]

    phrases = [{
        "object": t,
        "music": describe(reactor.parsed[t], p.notes_per_bar),
        "notes": [n for b in read_bars(reactor.parsed[t], p.notes_per_bar) for n in b["notes"]],
        "chords": [b["chord"] for b in read_bars(reactor.parsed[t], p.notes_per_bar)],
        "cadence": [b["function"] for b in read_bars(reactor.parsed[t], p.notes_per_bar)],
    } for t in reactor.phrases()]

    index = {t: i for i, t in enumerate(species_ids)}
    conservation = []
    for element in elements():
        vector = [0] * len(species_ids)
        for t in species_ids:
            vector[index[t]] = sum(line.count(element) for _, line in reactor.parsed[t])
        if any(vector):
            conservation.append({"name": f"element:{element}", "vector": vector})

    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={t: float(c) for t, c in pool.items()},
        extras={
            "analysis": {
                "phrases": phrases,
                "phrase_count": len(phrases),
                "elements": len(elements()),
                "initial_object_kinds": len(pool),
                "initial_objects": int(sum(pool.values())),
                "rules": len(rules),
                "failed_collisions": reactor.failed,
            },
            "conservation": conservation,
            "rules": [f"{r.label}: {r.text}" for r in rules],
            "reaction_rules": labels,
            "final_state": dict(Counter(
                {t: c for t, c in reactor.counts.items() if c > 0}).most_common()),
        },
    )
