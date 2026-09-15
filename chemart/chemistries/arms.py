"""ARMS, the Abstract Rewriting system on MultiSets (Suzuki & Tanaka). Catalog id: arms.

Book 9.4 (with 8.2.3 and 18.3.2); Suzuki & Tanaka, "Order parameter for a
symbolic chemical system", Artificial Life VI (1998) 130-139 (book [831]).

An ARMS is a pair Gamma = (A, R): symbols of an alphabet A (the species) and
multiset rewriting rules r_1..r_n (the reactions), written in the reaction
syntax of chemart.helpers.explicit ("3 a -> c", "a -> 2 a + 2 b"). The
reaction vessel is a multiset of symbols; input symbols may be injected, and
the multiset may be bounded by a maximal size. The explicit rule set is the
network; the rules carry no rates (kinetics is rule-application frequency).

The generator also runs the ARMS once from the initial multiset and records
the trajectory class in extras["analysis"]. Two rewriting disciplines:

- ordered (book 9.4, eqs. 9.15-9.18): at each step the rules are processed by
  index (ascending, or descending for the reverse order); each rule is applied
  as many times as its left side fits in what is left of the current multiset,
  and all products join the next multiset. This reproduces the book's halting
  sequence {a,b,a,a} -> {c,d} -> {e,f,f} -> {e,h,h} and its expanding reverse.
- random ([831], figs. 2 and 14): at each iteration the input symbols are
  injected (if the maximal size allows), one rule is selected, and it rewrites
  the multiset if its left side is present and the result does not exceed the
  maximal size; only rewrites count as steps. Selection is uniform over the
  rules, or (heating-probability) a heating rule (|rhs| > |lhs|) with
  probability p and a cooling rule (|rhs| < |lhs|) with probability 1 - p.

A run stops at the normal form ([831] def. 7): no selectable rule can be
applied and no symbol can be input without exceeding the maximal size.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations_with_replacement

from chemart.helpers.explicit import parse, term
from chemart.network import Network, Reaction, Species

SYSTEMS = ("book-example", "ru1", "brusselator", "two-symbol", "custom")
SELECTIONS = ("published", "order", "reverse-order", "random", "heating-probability")

#: Book eq. 9.15: a,a,a -> c : r1, b -> d : r2, c -> e : r3, d -> f,f : r4,
#: a -> a,b,b,a : r5, f -> h : r6; M0 = {a, b, a, a} (eq. 9.16).
BOOK_RULES = ["3 a -> c", "b -> d", "c -> e", "d -> 2 f", "a -> 2 a + 2 b", "f -> h"]
BOOK_ALPHABET = ["a", "b", "c", "d", "e", "f", "h"]

#: [831] p. 132: Ru1 = {aaa -> b : r1, b -> a : r2, b -> c : r3, a -> bb : r4}.
RU1_RULES = ["3 a -> b", "b -> a", "b -> c", "a -> 2 b"]

#: [831] fig. 8: the Brusselator as rewriting rules, inputs A and B.
BRUSSELATOR_RULES = ["A -> X", "B + X -> Y + D", "2 X + Y -> 3 X", "X -> E"]

#: name -> (rules, alphabet, initial multiset, inputs, maximal size or None, published selection)
PUBLISHED = {
    "book-example": (BOOK_RULES, BOOK_ALPHABET, {"a": 3, "b": 1}, [], None, "order"),
    "ru1": (RU1_RULES, ["a", "b", "c"], {"a": 4}, ["a"], 4, "random"),
    "brusselator": (BRUSSELATOR_RULES, ["A", "B", "X", "Y", "D", "E"], {}, ["A", "B"], 5000, "random"),
}


def _side(stoich: dict[str, int]) -> str:
    return " + ".join(term(n, s) for s, n in sorted(stoich.items()))


def two_symbol_rules() -> list[str]:
    """[831] 'rule set': both sides multisets over {a, b} of 1 to 5 symbols, left != right (380 rules)."""
    sides = [dict(Counter(c)) for k in range(1, 6) for c in combinations_with_replacement("ab", k)]
    return [f"{_side(l)} -> {_side(r)}" for l in sides for r in sides if l != r]


class Rule:
    def __init__(self, text: str, label: str):
        self.text, self.label = text, label
        self.lhs, self.rhs = parse(text)
        if not self.lhs and not self.rhs:
            raise ValueError(f"rule {text!r} is empty on both sides")
        self.delta = sum(self.rhs.values()) - sum(self.lhs.values())

    @property
    def kind(self) -> str:
        return "heating" if self.delta > 0 else "cooling" if self.delta < 0 else "neutral"

    def fits(self, state: dict[str, int], size: int, max_size: int | None) -> bool:
        return (all(state.get(s, 0) >= n for s, n in self.lhs.items())
                and (max_size is None or size + self.delta <= max_size))


def _apply(state: dict[str, int], rule: Rule, k: int = 1) -> None:
    for s, n in rule.lhs.items():
        state[s] -= k * n
        if not state[s]:
            del state[s]
    for s, n in rule.rhs.items():
        state[s] = state.get(s, 0) + k * n


def _key(state: dict[str, int]) -> tuple:
    return tuple(sorted(state.items()))


def run(rules: list[Rule], initial: dict[str, int], *, selection: str, steps: int, rng,
        inputs: list[str] = (), max_size: int | None = None, p: float = 0.5,
        input_probability: float = 1.0, sequence: list[int] | None = None,
        trace: bool = False) -> dict:
    """Run an ARMS; returns the analysis dict (plus 'trace' of states if asked).

    selection: 'order' / 'reverse-order' (book 9.4, one maximally parallel pass
    per step), 'random' or 'heating-probability' ([831] figs. 2 and 14).
    `sequence` (rule indices, cycled) replaces the random choice with a fixed
    rule order, one rule tried per iteration as in [831] fig. 3; the trace then
    has one state per iteration instead of per step.
    """
    state = {s: n for s, n in initial.items() if n > 0}
    size = sum(state.values())
    if max_size is not None and size > max_size:
        raise ValueError(f"the initial multiset has {size} symbols, more than max_size {max_size}")
    uses = [0] * len(rules)
    n_inputs = 0
    last_seen = {_key(state): 0}
    periods: Counter = Counter()
    states = [dict(state)]
    halted, period, step = False, None, 0
    deterministic = selection in ("order", "reverse-order") and input_probability in (0.0, 1.0)

    def can_input() -> bool:
        return bool(inputs) and input_probability > 0 and (max_size is None or size + 1 <= max_size)

    def do_inputs() -> None:
        nonlocal size, n_inputs
        for s in inputs:
            if (max_size is None or size + 1 <= max_size) and (
                    input_probability >= 1 or rng.random() < input_probability):
                state[s] = state.get(s, 0) + 1
                size += 1
                n_inputs += 1

    def record() -> bool:
        """Note the state after a step; True if a deterministic run became periodic."""
        nonlocal period
        key = _key(state)
        if key in last_seen:
            periods[step - last_seen[key]] += 1
            if deterministic:
                period = step - last_seen[key]
                return True
        last_seen[key] = step
        return False

    if selection in ("order", "reverse-order"):
        order = list(range(len(rules)))[:: 1 if selection == "order" else -1]
        while step < steps:
            before = size
            do_inputs()
            remaining, new, fired = dict(state), {}, False
            used = size
            for i in order:
                r = rules[i]
                if r.lhs:
                    k = min(remaining.get(s, 0) // n for s, n in r.lhs.items())
                else:
                    k = 1
                if max_size is not None and r.delta > 0:
                    k = min(k, (max_size - used) // r.delta)
                if k <= 0:
                    continue
                fired = True
                uses[i] += k
                used += k * r.delta
                for s, n in r.lhs.items():
                    remaining[s] -= k * n
                for s, n in r.rhs.items():
                    new[s] = new.get(s, 0) + k * n
            if not fired and size == before and not can_input():
                halted = True
                break
            for s, n in new.items():
                remaining[s] = remaining.get(s, 0) + n
            state.clear()
            state.update({s: n for s, n in remaining.items() if n > 0})
            size = sum(state.values())
            step += 1
            if trace:
                states.append(dict(state))
            if record():
                break
    elif selection in ("random", "heating-probability"):
        if selection == "heating-probability":
            heating = [i for i, r in enumerate(rules) if r.delta > 0]
            cooling = [i for i, r in enumerate(rules) if r.delta < 0]
            selectable = (heating if p > 0 else []) + (cooling if p < 1 else [])
        else:
            selectable = list(range(len(rules)))
        normal: dict[tuple, bool] = {}
        iteration, cap = 0, 1000 + 50 * steps * max(1, len(rules))
        while step < steps and iteration < cap:
            iteration += 1
            key = (_key(state), can_input())
            if key not in normal:
                normal[key] = not key[1] and not any(rules[i].fits(state, size, max_size) for i in selectable)
            if normal[key]:
                halted = True
                break
            do_inputs()
            if sequence is not None:
                i = sequence[(iteration - 1) % len(sequence)]
            elif selection == "random":
                i = int(rng.integers(len(rules))) if rules else None
            else:
                group = heating if rng.random() < p else cooling
                i = group[int(rng.integers(len(group)))] if group else None
            if i is not None and rules[i].fits(state, size, max_size):
                _apply(state, rules[i])
                size += rules[i].delta
                uses[i] += 1
                step += 1
                record()
            if trace and (sequence is not None or i is not None):
                states.append(dict(state))
    else:
        raise ValueError(f"unknown selection {selection!r}")

    heat = sum(u for u, r in zip(uses, rules) if r.delta > 0)
    cool = sum(u for u, r in zip(uses, rules) if r.delta < 0)
    if halted:
        cls = "halting"
    elif period is not None:
        cls = "periodic"
    elif periods:
        cls = "recurrent"
    else:
        cls = "non-recurrent"
    out = {
        "selection": selection,
        "steps": step,
        "class": cls,
        "halted": halted,
        "period": period,
        "final_state": {s: int(n) for s, n in sorted(state.items())},
        "final_size": int(size),
        "rule_uses": uses,
        "inputs": n_inputs,
        "heating_uses": heat,
        "cooling_uses": cool,
        # [831] eq. (3): lambda_e = sum r_heating / (1 + (sum r_cooling - 1)); undefined without cooling
        "lambda_e": float(heat / (1 + (cool - 1))) if cool else None,
        "revisits": int(sum(periods.values())),
        "kinds_of_periods": len(periods),
        "periods": sorted(int(t) for t in periods),
    }
    if trace:
        out["trace"] = states
    return out


def system(p, rng):
    """(rule texts, alphabet, initial multiset, inputs, max size or None, published selection)."""
    if p.system != "custom" and (p.rules or p.initial or p.inputs):
        raise ValueError("rules, initial and inputs are only used with system='custom'")
    if p.system != "two-symbol" and p.rule_count:
        raise ValueError("rule_count is only used with system='two-symbol'")
    if p.system in PUBLISHED:
        rules, alphabet, initial, inputs, max_size, sel = PUBLISHED[p.system]
        return list(rules), list(alphabet), dict(initial), list(inputs), max_size, sel
    if p.system == "two-symbol":
        rules = two_symbol_rules()
        if p.rule_count:
            chosen = rng.choice(len(rules), size=p.rule_count, replace=False)
            rules = [rules[int(i)] for i in sorted(chosen)]
        size = int(rng.integers(1, 11))
        initial = dict(Counter(str(s) for s in rng.choice(["a", "b"], size=size)))
        return rules, ["a", "b"], initial, [], 10, "heating-probability"
    if not isinstance(p.rules, list) or not p.rules or not all(isinstance(r, str) for r in p.rules):
        raise ValueError("system custom needs rules: a non-empty list of rules like '3 a -> c'")
    if not all(isinstance(k, str) and isinstance(v, int) and not isinstance(v, bool) and v >= 0
               for k, v in p.initial.items()):
        raise ValueError(f"initial must map symbols to non-negative integer counts, got {p.initial!r}")
    if not all(isinstance(s, str) and s and " " not in s for s in p.inputs):
        raise ValueError(f"inputs must be a list of symbols, got {p.inputs!r}")
    return list(p.rules), [], dict(p.initial), list(p.inputs), None, "random"


def generate(p, rng):
    texts, alphabet, initial, inputs, published_max, published_sel = system(p, rng)
    rules = [Rule(text, f"r{i + 1}") for i, text in enumerate(texts)]
    max_size = published_max if p.max_size == 0 else None if p.max_size < 0 else p.max_size
    selection = published_sel if p.selection == "published" else p.selection

    analysis = run(rules, initial, selection=selection, steps=p.steps, rng=rng, inputs=inputs,
                   max_size=max_size, p=p.p, input_probability=p.input_probability)

    order = dict.fromkeys(alphabet)
    for r in rules:
        order.update(dict.fromkeys([*r.lhs, *r.rhs]))
    order.update(dict.fromkeys([*initial, *inputs]))
    reactions = [Reaction(dict(r.lhs), dict(r.rhs)) for r in rules]
    labels = [r.label for r in rules]
    for s in dict.fromkeys(inputs):
        reactions.append(Reaction({}, {s: 1}))
        labels.append("input")
    extras = {
        "system": p.system,
        "rules": [f"{r.label}: {r.text}" for r in rules],
        "rule_kinds": [r.kind for r in rules],
        "reaction_rules": labels,
        "inputs": list(inputs),
        "max_size": max_size,
        "analysis": analysis,
    }
    return Network(species=[Species(s) for s in order], reactions=reactions, status="complete",
                   initial_state={s: float(n) for s, n in initial.items()}, extras=extras)
