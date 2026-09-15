"""L-systems (Lindenmayer systems): parallel string rewriting.

Catalog id: l-systems. Book 9.8 [512, 677]; Prusinkiewicz & Lindenmayer, The
Algorithmic Beauty of Plants (ABOP, 1990), chapter 1.

An L-system is an alphabet, a nonempty axiom word and a set of productions.
In one derivation step EVERY letter of the word is replaced simultaneously by
the successor of a production whose predecessor it is (letters without an
applicable production are replaced by themselves). Production syntax:

    a -> ab                     context-free (0L)
    b < a -> b                  left context (1L)
    F_a > F_b -> F_b            right context (1L)
    0 < 1 > 0 -> 1[+F1F1]       both contexts (2L); '*' stands for no context
    F -> F[+F]F : 0.33          stochastic production with its probability
    a -> ε                      empty successor (also 'a ->')

Symbols are single non-space characters, optionally followed by an
underscore subscript (a_r, F_b, a_12), so ABOP's subscripted letters can be
written directly; spaces inside words are ignored. Brackets delimit branches
and replace themselves. Context matching follows ABOP 1.8: symbols in
`ignore` are skipped, a left-context search moves out of branches ('[') and
over complete branches (']'), a right-context search skips complete lateral
branches, and ']' in a right context skips the rest of the current branch.
Context-sensitive productions take precedence over context-free ones; among
context-sensitive productions the first listed that matches applies.

Chemical readings (`reading`):

- words: species are the words of the derivation (structure: the word) and
  each derivation step is a unimolecular reaction w -> w'. A stochastic system
  gives one reaction per distinct successor word, its probability (the sum
  over the independent per-letter choices that produce it) in
  extras["probabilities"], aligned with the reactions. A word that derives
  itself is an elastic step (dropped, as in chemart.expand); its probability is
  kept in extras["stationary"].
- symbols: species are the letters and each production a -> chi is the
  reaction a -> multiset(chi). This loses order and context (context-sensitive
  systems are rejected) but its stoichiometry is exactly ABOP's growth matrix
  Q (section 1.9): letter counts evolve as n_{k+1} = n_k Q.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from itertools import product

from chemart.expand import expand
from chemart.network import Network, Reaction, Species

Word = tuple[str, ...]

BRACKETS = ("[", "]")
_TOKEN = re.compile(r"\S(?:_[A-Za-z0-9]+)?")
_PROBABILITY = re.compile(r"^(.*?)(?:(?:^|\s+):\s*(\d*\.?\d+(?:[eE][-+]?\d+)?))?$")

ABOP = "Prusinkiewicz & Lindenmayer, The Algorithmic Beauty of Plants (1990)"

#: The published systems. n is the derivation length of the published figure.
SYSTEMS: dict[str, dict] = {
    "algae": dict(
        axiom="a", productions=["a -> ab", "b -> a"],
        source="ABOP eq. (1.2), section 1.9 (letter counts follow the Fibonacci series); "
               "the same productions from axiom b are the derivation of section 1.2, figure 1.3"),
    "anabaena": dict(
        axiom="a_r", productions=["a_r -> a_l b_r", "a_l -> b_l a_r", "b_r -> a_r", "b_l -> a_l"],
        source="ABOP eq. (1.1), figure 1.4: filament of Anabaena catenula"),
    "koch-snowflake": dict(
        axiom="F--F--F", productions=["F -> F+F--F+F"], angle=60.0,
        source="the Koch construction of ABOP figure 1.1 (generator of N = 4 sides of length 1/3), "
               "coded as ABOP section 1.3 prescribes (initiator = axiom, generator = successor)"),
    "koch-island": dict(
        axiom="F-F-F-F", productions=["F -> F-F+F+FF-F-F+F"], angle=90.0, n=3,
        source="ABOP section 1.3, figure 1.6: quadratic Koch island (derivations 0 to 3)"),
    "quadratic-koch-island": dict(
        axiom="F-F-F-F", productions=["F -> F+FF-FF-F-F+F+FF-F-F+F+FF+FF-F"], angle=90.0, n=2,
        source="ABOP figure 1.7a: quadratic Koch island (Mandelbrot p. 52)"),
    "quadratic-snowflake": dict(
        axiom="-F", productions=["F -> F+F-F-F+F"], angle=90.0, n=4,
        source="ABOP figure 1.7b: quadratic modification of the snowflake curve"),
    "islands-and-lakes": dict(
        axiom="F+F+F+F", productions=["F -> F+f-FF+F+FF+Ff+FF-f+FF-F-FF-Ff-FFF", "f -> ffffff"],
        angle=90.0, n=2,
        source="ABOP figure 1.8: combination of islands and lakes (f moves without drawing)"),
    "koch-variant-a": dict(axiom="F-F-F-F", productions=["F -> FF-F-F-F-F-F+F"], angle=90.0, n=4,
                           source="ABOP figure 1.9a"),
    "koch-variant-b": dict(axiom="F-F-F-F", productions=["F -> FF-F-F-F-FF"], angle=90.0, n=4,
                           source="ABOP figure 1.9b"),
    "koch-variant-c": dict(axiom="F-F-F-F", productions=["F -> FF-F+F-F-FF"], angle=90.0, n=3,
                           source="ABOP figure 1.9c"),
    "koch-variant-d": dict(axiom="F-F-F-F", productions=["F -> FF-F--F-F"], angle=90.0, n=4,
                           source="ABOP figure 1.9d"),
    "koch-variant-e": dict(axiom="F-F-F-F", productions=["F -> F-FF--F-F"], angle=90.0, n=5,
                           source="ABOP figure 1.9e"),
    "koch-variant-f": dict(axiom="F-F-F-F", productions=["F -> F-F+F-F-F"], angle=90.0, n=4,
                           source="ABOP figure 1.9f"),
    "dragon-curve": dict(
        axiom="F_l", productions=["F_l -> F_l+F_r+", "F_r -> -F_l-F_r"], angle=90.0, n=10,
        source="ABOP figure 1.10a (edge rewriting; F_l and F_r both draw a segment)"),
    "sierpinski-gasket": dict(
        axiom="F_r", productions=["F_l -> F_r+F_l+F_r", "F_r -> F_l-F_r-F_l"], angle=60.0, n=6,
        source="ABOP figure 1.10b (edge rewriting)"),
    "plant-a": dict(axiom="F", productions=["F -> F[+F]F[-F]F"], angle=25.7, n=5,
                    source="ABOP figure 1.24a (bracketed 0L, edge rewriting)"),
    "plant-b": dict(axiom="F", productions=["F -> F[+F]F[-F][F]"], angle=20.0, n=5,
                    source="ABOP figure 1.24b"),
    "plant-c": dict(axiom="F", productions=["F -> FF-[-F+F+F]+[+F-F-F]"], angle=22.5, n=4,
                    source="ABOP figure 1.24c"),
    "plant-d": dict(axiom="X", productions=["X -> F[+X]F[-X]+X", "F -> FF"], angle=20.0, n=7,
                    source="ABOP figure 1.24d (node rewriting)"),
    "plant-e": dict(axiom="X", productions=["X -> F[+X][-X]FX", "F -> FF"], angle=25.7, n=7,
                    source="ABOP figure 1.24e"),
    "plant-f": dict(axiom="X", productions=["X -> F-[[X]+X]+F[+FX]-X", "F -> FF"], angle=22.5, n=5,
                    source="ABOP figure 1.24f"),
    "stochastic-plant": dict(
        axiom="F",
        productions=["F -> F[+F]F[-F]F : 0.33", "F -> F[+F]F : 0.33", "F -> F[-F]F : 0.34"],
        n=5, source="ABOP section 1.7, figure 1.27 (stochastic bracketed 0L; no angle given)"),
    "signal-propagation": dict(
        axiom="baaaaaaaa", productions=["b < a -> b", "b -> a"],
        source="ABOP section 1.8: 1L-system propagating a signal along a string"),
    "acropetal-signal": dict(
        axiom="F_b[+F_a]F_a[-F_a]F_a[+F_a]F_a", productions=["F_b < F_a -> F_b"], ignore="+-",
        source="ABOP section 1.8, figure 1.30a: acropetal signal in a bracketed string"),
    "basipetal-signal": dict(
        axiom="F_a[+F_a]F_a[-F_a]F_a[+F_a]F_b", productions=["F_a > F_b -> F_b"], ignore="+-",
        source="ABOP section 1.8, figure 1.30b: basipetal signal in a bracketed string"),
    "square-root-growth": dict(
        axiom="X F_u F_a X",
        productions=["F_u < F_a > F_a -> F_u", "F_u < F_a > X -> F_d F_a", "F_a < F_a > F_d -> F_d",
                     "X < F_a > F_d -> F_u", "F_u -> F_a", "F_d -> F_a"],
        source="ABOP section 1.9, eq. (1.5), figure 1.33: 2L-system with square-root growth"),
}

_HOGEWEG = {
    # (angle, n, axiom, successors of 000 001 010 011 100 101 110 111, successors of + and -)
    "a": (22.5, 30, "F1F1F1", ["0", "1[+F1F1]", "1", "1", "0", "1F1", "0", "0"], ("-", "+")),
    "b": (22.5, 30, "F1F1F1", ["1", "1[-F1F1]", "1", "1", "0", "1F1", "1", "0"], ("-", "+")),
    "c": (25.75, 26, "F1F1F1", ["0", "1", "0", "1[+F1F1]", "0", "1F1", "0", "0"], ("-", "+")),
    "d": (25.75, 24, "F0F1F1", ["1", "0", "0", "1F1", "1", "1[+F1F1]", "1", "0"], ("-", "+")),
    "e": (22.5, 26, "F1F1F1", ["0", "1[-F1F1]", "1", "1", "0", "1F1", "1", "0"], ("-", "+")),
}
for _label, (_angle, _n, _axiom, _succ, (_plus, _minus)) in _HOGEWEG.items():
    SYSTEMS[f"hogeweg-hesper-{_label}"] = dict(
        axiom=_axiom, ignore="+-F", angle=_angle, n=_n,
        productions=[f"{c[0]} < {c[1]} > {c[2]} -> {s}"
                     for c, s in zip(("000", "001", "010", "011", "100", "101", "110", "111"), _succ)]
        + [f"* < + > * -> {_plus}", f"* < - > * -> {_minus}"],
        source=f"ABOP section 1.8, figure 1.31{_label}: bracketed 2L-system after Hogeweg and Hesper (1974)")

CUSTOM_ONLY = ("axiom", "productions", "ignore", "angle")


# --- parsing ---------------------------------------------------------------------------------
def symbols(text: str) -> Word:
    """Split a word into symbols: one character plus an optional _subscript; spaces are ignored."""
    word = tuple(_TOKEN.findall(text))
    if "_" in word:
        raise ValueError(f"'_' is not a symbol on its own (it starts a subscript, as in a_r): {text!r}")
    return word


def _balanced(word: Word, what: str) -> None:
    depth = 0
    for s in word:
        depth += (s == "[") - (s == "]")
        if depth < 0:
            break
    if depth != 0:
        raise ValueError(f"unbalanced brackets in {what}: {''.join(word)!r}")


@dataclass(frozen=True)
class Production:
    left: Word
    pred: str
    right: Word
    succ: Word
    probability: float | None
    text: str

    @property
    def key(self) -> tuple[Word, str, Word]:
        return (self.left, self.pred, self.right)

    @property
    def context_free(self) -> bool:
        return not self.left and not self.right


def parse_production(text) -> Production:
    if not isinstance(text, str):
        raise ValueError(f"every production must be a string like 'a -> ab', got {text!r}")
    body = text.replace("→", "->")
    if body.count("->") != 1:
        raise ValueError(f"production {text!r} needs exactly one '->'")
    lhs, rhs = body.split("->")
    succ_text, prob_text = _PROBABILITY.match(rhs.strip()).groups()
    succ_text = succ_text.strip()
    succ = () if succ_text in ("", "ε") else symbols(succ_text)
    probability = None
    if prob_text is not None:
        probability = float(prob_text)
        if not 0 < probability <= 1:
            raise ValueError(f"production {text!r}: probability must be in (0, 1], got {probability}")

    lhs = f" {lhs.strip()} "
    parts = re.split(r"\s<\s", lhs)
    if len(parts) > 2:
        raise ValueError(f"production {text!r} has more than one '<'")
    left_text, rest = parts if len(parts) == 2 else ("*", parts[0])
    parts = re.split(r"\s>\s", f" {rest.strip()} ")
    if len(parts) > 2:
        raise ValueError(f"production {text!r} has more than one '>'")
    pred_text, right_text = parts if len(parts) == 2 else (parts[0], "*")

    def context(t: str) -> Word:
        return () if t.strip() == "*" else symbols(t)

    left, right, pred = context(left_text), context(right_text), symbols(pred_text)
    if len(pred) != 1:
        raise ValueError(
            f"production {text!r}: the strict predecessor must be one symbol (write contexts as 'l < a > r')")
    if pred[0] in BRACKETS:
        raise ValueError(f"production {text!r}: brackets replace themselves and cannot be rewritten")
    if any(s in BRACKETS for s in left):
        raise ValueError(f"production {text!r}: a left context is a path and cannot contain brackets")
    _balanced(right, f"the right context of {text!r}")
    _balanced(succ, f"the successor of {text!r}")
    return Production(left, pred[0], right, succ, probability, text)


class LSystem:
    """Axiom, productions and #ignore set, with the parallel derivation operator."""

    def __init__(self, axiom: str, productions: list, ignore: str = ""):
        if not isinstance(axiom, str) or not symbols(axiom):
            raise ValueError(f"the axiom must be a nonempty word, got {axiom!r}")
        self.axiom = symbols(axiom)
        _balanced(self.axiom, "the axiom")
        if not isinstance(productions, list) or not productions:
            raise ValueError("productions must be a non-empty list of strings like 'a -> ab'")
        self.productions = [parse_production(t) for t in productions]
        self.ignore = frozenset(symbols(ignore))
        if self.ignore & set(BRACKETS):
            raise ValueError("brackets cannot be ignored during context matching")

        groups: dict[tuple, list[Production]] = {}
        for prod in self.productions:
            groups.setdefault(prod.key, []).append(prod)
        self.context_free: dict[str, tuple] = {}
        self.context_sensitive: dict[str, list] = {}
        for (left, pred, right), prods in groups.items():
            given = [p.probability for p in prods]
            if len(prods) > 1 or given[0] is not None:
                if None in given:
                    name = prods[0].text.split("->")[0].strip()
                    raise ValueError(
                        f"productions with predecessor {name!r} have alternatives: give every one a "
                        f"probability (': 0.5')")
                if abs(sum(given) - 1) > 1e-6:
                    raise ValueError(
                        f"probabilities of the productions {[p.text for p in prods]} sum to {sum(given)}, not 1")
            alternatives = tuple((p.succ, 1.0 if p.probability is None else p.probability) for p in prods)
            if left or right:
                self.context_sensitive.setdefault(pred, []).append((left, right, alternatives))
            else:
                self.context_free[pred] = alternatives

    @property
    def stochastic(self) -> bool:
        return any(p.probability is not None and p.probability < 1 for p in self.productions)

    # --- context ------------------------------------------------------------------------------
    @staticmethod
    def _match_brackets(word: Word) -> list[int]:
        match = [-1] * len(word)
        stack = []
        for i, s in enumerate(word):
            if s == "[":
                stack.append(i)
            elif s == "]":
                j = stack.pop()
                match[i], match[j] = j, i
        return match

    def _left(self, word: Word, i: int, left: Word, match: list[int]) -> bool:
        j = i - 1
        for c in reversed(left):
            while j >= 0:
                s = word[j]
                if s == "]":
                    j = match[j] - 1          # skip a complete branch
                elif s == "[" or s in self.ignore:
                    j -= 1                    # leave the branch towards its parent
                else:
                    break
            if j < 0 or word[j] != c:
                return False
            j -= 1
        return True

    def _right(self, word: Word, i: int, right: Word, match: list[int]) -> bool:
        j, n = i + 1, len(word)
        for c in right:
            if c == "]":                      # skip the rest of the current branch
                while j < n and word[j] != "]":
                    j = match[j] + 1 if word[j] == "[" else j + 1
                if j >= n:
                    return False
                j += 1
                continue
            while j < n:
                s = word[j]
                if s == "[" and c != "[":
                    j = match[j] + 1          # skip a lateral branch
                elif s in self.ignore:
                    j += 1
                else:
                    break
            if j >= n or word[j] != c:
                return False
            j += 1
        return True

    def choices(self, word: Word, i: int, match: list[int]) -> tuple:
        """Alternative (successor, probability) pairs for the letter at position i."""
        s = word[i]
        for left, right, alternatives in self.context_sensitive.get(s, ()):
            if self._left(word, i, left, match) and self._right(word, i, right, match):
                return alternatives
        return self.context_free.get(s) or (((s,), 1.0),)

    # --- derivation ---------------------------------------------------------------------------
    def successors(self, word: Word, max_alternatives: int = 10_000,
                   max_length: int | None = None) -> tuple[dict[Word, float], bool]:
        """Every word directly derived from `word`, with its probability.

        Returns (successors, cut): cut is True when the word has more than
        `max_alternatives` combinations of per-letter choices (none is then
        returned) or when successors longer than `max_length` were left out.
        """
        match = self._match_brackets(word) if self.context_sensitive else []
        per = [self.choices(word, i, match) for i in range(len(word))]
        total = 1
        for alternatives in per:
            total *= len(alternatives)
            if total > max_alternatives:
                return {}, True
        fixed: list[list[str]] = [[]]
        variable = []
        for alternatives in per:
            if len(alternatives) == 1:
                fixed[-1].extend(alternatives[0][0])
            else:
                variable.append(alternatives)
                fixed.append([])
        out: dict[Word, float] = {}
        cut = False
        for combo in product(*variable):
            parts = list(fixed[0])
            probability = 1.0
            for (succ, pr), tail in zip(combo, fixed[1:]):
                parts.extend(succ)
                parts.extend(tail)
                probability *= pr
            if max_length is not None and len(parts) > max_length:
                cut = True
                continue
            key = tuple(parts)
            out[key] = out.get(key, 0.0) + probability
        return out, cut

    def derive(self, word: Word | None = None) -> Word:
        """One deterministic derivation step (ValueError for a stochastic system)."""
        if self.stochastic:
            raise ValueError("derive() needs a deterministic system; use successors()")
        (succ,), _ = self.successors(self.axiom if word is None else word)
        return succ

    def derivation(self, n: int) -> list[Word]:
        words = [self.axiom]
        for _ in range(n):
            words.append(self.derive(words[-1]))
        return words


# --- parameters --------------------------------------------------------------------------------
def system(p) -> tuple[LSystem, dict]:
    if p.system == "custom":
        lsys = LSystem(p.axiom, p.productions, p.ignore)
        info = {"source": "custom"}
        if p.angle:
            info["angle"] = float(p.angle)
        return lsys, info
    given = [name for name in CUSTOM_ONLY if getattr(p, name)]
    if given:
        raise ValueError(f"{', '.join(given)} are only used with system='custom' (got system={p.system!r})")
    spec = SYSTEMS[p.system]
    lsys = LSystem(spec["axiom"], spec["productions"], spec.get("ignore", ""))
    info = {"source": f"{ABOP}: {spec['source']}"}
    if "angle" in spec:
        info["angle"] = spec["angle"]
    if "n" in spec:
        info["published_iterations"] = spec["n"]
    return lsys, info


def word_text(word: Word) -> str:
    """Text of a word that symbols() reads back: compact, or space-separated if a symbol has a subscript."""
    return " ".join(word) if any(len(s) > 1 for s in word) else "".join(word)



def generate(p, rng):
    lsys, info = system(p)
    extras = {
        "system": p.system,
        "reading": p.reading,
        "axiom": word_text(lsys.axiom),
        "productions": [prod.text for prod in lsys.productions],
        "ignore": word_text(tuple(sorted(lsys.ignore))),
        "stochastic": lsys.stochastic,
        "source": info["source"],
    }
    turtle = {k: info[k] for k in ("angle", "published_iterations") if k in info}
    if turtle:
        extras["turtle"] = turtle
    if p.reading == "symbols":
        return _symbols(p, lsys, extras)
    return _words(p, lsys, extras)


def _words(p, lsys: LSystem, extras: dict) -> Network:
    depth = {lsys.axiom: 0}
    probabilities: dict[tuple[Word, Word], float] = {}
    stationary: dict[Word, float] = {}
    cut = []

    def react(word):
        d = depth[word]
        if d >= p.iterations:
            return None
        outcomes, was_cut = lsys.successors(word, p.max_species, p.max_length)
        if was_cut:
            cut.append(word)
        alternatives = []
        for succ, probability in outcomes.items():
            if succ == word:
                stationary[word] = probability
                continue
            depth.setdefault(succ, d + 1)
            probabilities[(word, succ)] = probability
            alternatives.append((succ,))
        return alternatives or None

    words, pairs, status = expand(react, [lsys.axiom], arity=1, max_species=p.max_species,
                                  alternatives=True)
    ids = {w: f"w{i}" for i, w in enumerate(words)}
    reactions = [Reaction({ids[lhs[0]]: 1}, {ids[rhs[0]]: 1}) for lhs, rhs in pairs]
    extras.update({
        "probabilities": [probabilities[(lhs[0], rhs[0])] for lhs, rhs in pairs],
        "stationary": {ids[w]: pr for w, pr in stationary.items()},
        "analysis": {
            "steps": [depth[w] for w in words],
            "lengths": [len(w) for w in words],
            "unexpanded": [ids[w] for w in cut if w in ids],
        },
    })
    return Network(
        species=[Species(ids[w], structure=word_text(w)) for w in words],
        reactions=reactions,
        status="truncated" if status == "truncated" or cut else "complete",
        initial_state={ids[lsys.axiom]: 1.0},
        extras=extras,
    )


def _symbols(p, lsys: LSystem, extras: dict) -> Network:
    if lsys.context_sensitive:
        raise ValueError(
            "reading='symbols' cannot express context-sensitive productions (letters lose their "
            "neighbours); use reading='words'")
    order = dict.fromkeys(lsys.axiom)
    for prod in lsys.productions:
        order.update(dict.fromkeys((prod.pred, *prod.succ)))
    reactions, probabilities, stationary = [], [], {}
    for prod in lsys.productions:
        probability = 1.0 if prod.probability is None else prod.probability
        if prod.succ == (prod.pred,):
            stationary[prod.pred] = stationary.get(prod.pred, 0.0) + probability
            continue
        reactions.append(Reaction({prod.pred: 1}, dict(Counter(prod.succ))))
        probabilities.append(probability)

    # Growth (ABOP 1.9): n_{k+1} = n_k Q, expected counts for a stochastic system.
    counts = Counter(lsys.axiom)
    history = [dict(counts)]
    for _ in range(p.iterations):
        new: Counter = Counter()
        for s, c in counts.items():
            for succ, pr in lsys.context_free.get(s, (((s,), 1.0),)):
                for t, m in Counter(succ).items():
                    new[t] += c * m * pr if lsys.stochastic else c * m
        counts = +new
        history.append({s: counts[s] for s in order if counts[s]})
    extras.update({
        "probabilities": probabilities,
        "stationary": stationary,
        "analysis": {
            "counts": history,
            "lengths": [sum(h.values()) for h in history],
        },
    })
    if not lsys.stochastic:
        extras["analysis"]["lengths"] = [int(v) for v in extras["analysis"]["lengths"]]
    return Network(
        species=[Species(s) for s in order],
        reactions=reactions,
        status="complete",
        initial_state={s: float(c) for s, c in Counter(lsys.axiom).items()},
        extras=extras,
    )
