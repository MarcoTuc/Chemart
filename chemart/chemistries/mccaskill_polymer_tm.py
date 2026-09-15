"""McCaskill's pattern processing chemistry: polymers as Turing machines.

Catalog id: mccaskill-polymer-tm. Book 10.5.3; reconstructed from McCaskill
(1988), "Polymer chemistry on tape: a computational model for emergent
genetics", MPI for Biophysical Chemistry internal report (sections 2-4).

Molecules are binary strings. Every string is read in two ways, both with
doublet codes (report section 3, items 4-5):

* recognizons: from each initiator triplet 111, doublets 00/11 -> #,
  01 -> 0, 10 -> 1 (Holland conditions over {0, 1, #}); a later 111 on a
  doublet boundary or the end of the string terminates, at most R symbols.
  Two molecules collide when they place the same specific pattern (each #
  resolved at random), i.e. when two recognizons of equal length agree
  wherever neither has a #.
* processor rules: 12-bit sextuples (READ, WRITE, STATE, NEXT STATE, READ
  MOTION, WRITE MOTION) marked by the initiator. The processor is a two-state
  head that reads the other molecule (the tape, not consumed) and writes new
  strings on a blank tape; it halts when no rule matches, and a later rule
  for the same (read symbol, state) overwrites an earlier one.

A collision processor + tape therefore gives s1 + s2 -> s1 + s2 + s3 (+ ...):
both reactants survive and the released strings are new molecules.

method "closure" returns every reaction reachable from the seed strings
(chemart.expand.expand over ordered processor/tape pairs); method "soup" runs
the report's collision algorithm (a pattern space filled by randomly drawn
molecules, products displacing random molecules) and returns the reactions
that fired.
"""

from __future__ import annotations

from collections import Counter

from chemart.expand import expand
from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species

INITIATOR = "111"
RULE_BITS = 12
BLANK = None
OPPOSITE = "opposite"
SAME = "same"

# The report's published strings (section 4).
REPLICATOR = "0001111100111010111"
PARASITE = "0000000100100010111"

# Report section 3, item 5: doublet codes of the six rule components.
_READ = {"01": (0,), "10": (1,), "00": (BLANK,), "11": (0, 1)}
_WRITE = {"01": 0, "10": 1, "00": BLANK, "11": OPPOSITE}
_STATE = {"01": (0,), "10": (1,), "00": (0, 1), "11": (0, 1)}
_NEXT = {"01": 0, "10": 1, "00": SAME, "11": OPPOSITE}
_MOTION = {"00": -1, "01": 0, "10": 1, "11": 1}
# Report section 3, item 4: recognizon doublet code.
_RECOGNIZON = {"01": "0", "10": "1", "00": "#", "11": "#"}


def complement(s: str) -> str:
    return s.translate(str.maketrans("01", "10"))


def decode_rules(s: str) -> list[tuple[int, dict]]:
    """Rules as (position of the initiator, rule): the 12 bits starting at every 111."""
    out = []
    for i in range(len(s) - RULE_BITS + 1):
        if s.startswith(INITIATOR, i):
            d = [s[i + 2 * k:i + 2 * k + 2] for k in range(6)]
            out.append((i, {
                "read": _READ[d[0]], "write": _WRITE[d[1]], "state": _STATE[d[2]],
                "next": _NEXT[d[3]], "read_motion": _MOTION[d[4]], "write_motion": _MOTION[d[5]],
            }))
    return out


def recognizons(s: str, R: int = 16) -> list[str]:
    """Condition strings over {0, 1, #}, one per initiator, in string order."""
    out = []
    n = len(s)
    for i in range(n - 2):
        if not s.startswith(INITIATOR, i):
            continue
        pattern, j = "", i
        while j + 1 < n and len(pattern) < R:
            if j > i + 2 and s.startswith(INITIATOR, j):
                break
            pattern += _RECOGNIZON[s[j:j + 2]]
            j += 2
        out.append(pattern)
    return out


def _compatible(a: str, b: str) -> bool:
    return len(a) == len(b) and all(x == y or x == "#" or y == "#" for x, y in zip(a, b))


def can_collide(a: str, b: str, R: int = 16) -> bool:
    """True if some specific pattern can be placed by both molecules."""
    rb = recognizons(b, R)
    return any(_compatible(x, y) for x in recognizons(a, R) for y in rb)


def _table(processor: str) -> dict:
    table = {}
    for _, rule in decode_rules(processor):   # later rules overwrite earlier ones
        for r in rule["read"]:
            for st in rule["state"]:
                table[(r, st)] = rule
    return table


def process(processor: str, tape: str, max_steps: int = 1000,
            error_rate: float = 0.0, rng=None, table: dict | None = None) -> list[str]:
    """Strings released when `processor` processes `tape` ([] if none or cut off)."""
    table = _table(processor) if table is None else table
    if not table:
        return []
    rp = wp = 0
    state = 0
    written: dict[int, int] = {}
    for _ in range(max_steps):
        r = int(tape[rp]) if 0 <= rp < len(tape) else BLANK
        rule = table.get((r, state))
        if rule is None:
            break
        w = rule["write"]
        if w == OPPOSITE:
            w = BLANK if r is BLANK else 1 - r
        if w is not BLANK and error_rate > 0 and rng.random() < error_rate:
            w = 1 - w
        if w is BLANK:
            written.pop(wp, None)
        else:
            written[wp] = w
        nxt = rule["next"]
        state = state if nxt == SAME else (1 - state if nxt == OPPOSITE else nxt)
        rp += rule["read_motion"]
        wp += rule["write_motion"]
    else:
        return []   # cut off: uncompleted processes have no effect
    strings, current = [], ""
    for pos in range(min(written, default=0), max(written, default=-1) + 2):
        if pos in written:
            current += str(written[pos])
        elif current:
            strings.append(current)
            current = ""
    return strings


def species_id(s: str) -> str:
    return f"p{s}"


def _parse_strings(value) -> list[str]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"strings must be a non-empty list of binary strings, got {value!r}")
    for s in value:
        if not isinstance(s, str) or not s or set(s) - {"0", "1"} or len(s) > 256:
            raise ValueError(f"strings entries must be binary strings of 1-256 characters '0'/'1', got {s!r}")
    return list(dict.fromkeys(value))


def generate(p, rng):
    strings = _parse_strings(p.strings)
    if p.method == "closure":
        if p.error_rate != 0:
            raise ValueError("error_rate only applies to method 'soup'; the closure is the error-free chemistry")
        return _closure(p, strings)
    return _soup(p, strings, rng)


def _reaction(lhs, rhs, count=None) -> Reaction:
    return Reaction.of([species_id(s) for s in lhs], [species_id(s) for s in rhs], count=count)


def _species(strings) -> list[Species]:
    return [Species(species_id(s), structure=s) for s in sorted(set(strings), key=lambda s: (len(s), s))]


def _closure(p, seeds):
    def react(a, b):
        if p.recognition == "patterns" and not can_collide(a, b, p.R):
            return None
        out = process(a, b, p.max_steps)
        return (a, b, *out) if out else None

    found, reactions, status = expand(react, seeds, arity=2, max_species=p.max_species, ordered=True)
    return Network(
        species=_species(found),
        reactions=[_reaction(lhs, rhs) for lhs, rhs in reactions],
        status=status,
        extras={"seed": [species_id(s) for s in seeds]},
    )


def _soup(p, inoculum, rng):
    n = p.population
    per = int(round(p.inoculum_fraction * n))
    if per * len(inoculum) > n:
        raise ValueError(f"inoculum_fraction x len(strings) = {per * len(inoculum)} molecules exceeds population {n}")
    pop = [s for s in inoculum for _ in range(per)]
    while len(pop) < n:
        pop.append("".join("1" if b else "0" for b in rng.integers(0, 2, size=p.l)))
    order = rng.permutation(n)
    pop = [pop[int(i)] for i in order]
    start = list(pop)

    space: dict[str, int] = {}            # specific pattern -> slot that placed it
    placed: dict[int, set[str]] = {}      # slot -> its patterns in the space
    rec_cache: dict[str, list[str]] = {}
    table_cache: dict[str, dict] = {}
    fired: dict[tuple, list] = {}
    collisions = 0

    def clear(slot):
        for pat in placed.pop(slot, ()):
            if space.get(pat) == slot:
                del space[pat]

    for _ in range(p.steps):
        i = int(rng.integers(n))
        s = pop[i]
        if p.recognition == "none":        # well-mixed limit: a random partner, no pattern space
            j = int(rng.integers(n - 1))
            j += j >= i
            collisions += 1
            proc, tape = pop[j], s
            table = table_cache.get(proc)
            if table is None:
                table = table_cache[proc] = _table(proc)
            products = process(proc, tape, p.max_steps, p.error_rate, rng, table=table)
            if products:
                for prod in products:
                    pop[int(rng.integers(n))] = prod
                lhs, rhs = (proc, tape), (proc, tape, *products)
                key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
                fired.setdefault(key, [lhs, rhs, 0])[2] += 1
            continue
        recs = rec_cache.get(s)
        if recs is None:
            recs = rec_cache[s] = recognizons(s, p.R)
        if not recs:
            continue
        cond = recs[int(rng.integers(len(recs)))]
        pat = "".join(c if c != "#" else ("1" if rng.random() < 0.5 else "0") for c in cond)
        j = space.get(pat)
        if j is None:
            space[pat] = i
            placed.setdefault(i, set()).add(pat)
            continue
        if j == i:
            continue
        collisions += 1
        proc, tape = pop[j], s
        clear(i)
        clear(j)
        table = table_cache.get(proc)
        if table is None:
            table = table_cache[proc] = _table(proc)
        products = process(proc, tape, p.max_steps, p.error_rate, rng, table=table)
        if not products:
            continue
        for prod in products:              # Moran displacement: overwrite a random molecule
            k = int(rng.integers(n))
            clear(k)
            pop[k] = prod
        lhs, rhs = (proc, tape), (proc, tape, *products)
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        entry = fired.setdefault(key, [lhs, rhs, 0])
        entry[2] += 1

    names = set(start).union(*(rhs for _, rhs, _ in fired.values()))
    return Network(
        species=_species(names),
        reactions=[_reaction(lhs, rhs, count) for lhs, rhs, count in fired.values()],
        status="observed",
        initial_state={species_id(s): c for s, c in sorted(Counter(start).items())},
        outflow=CONSTANT_TOTAL,
        extras={
            "analysis": {"steps": p.steps, "recognition_collisions": collisions},
            "final_state": {species_id(s): c for s, c in Counter(pop).most_common()},
        },
    )
