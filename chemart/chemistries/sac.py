"""SAC: string-based artificial chemistry with cells (Suzuki & Ono, 2002-2003).

Catalog id: sac.

Molecules are strings over 19 symbols (book table 11.4). Two strings react:
the first (operator) is decoded into rewriting rules and applied to the second
(operand); both are put back, so the operator is a catalyst:

    operator + operand -> operator + pieces of the rewritten operand

The operand is cut into pieces at every unsuppressed '.'. A collision whose
rules leave the operand unchanged is elastic.

Decoding (reconstructed from Suzuki's ECAL 2003 slides; see catalog decisions):

- ``\\x`` suppresses x: it matches the character x literally.
- ``'x`` deletes and ``"x`` creates the next symbol x (a literal, ``\\x``, or a
  wild card); every other symbol is kept (matched and left in place).
- ``!``, ``?``, ``%`` match one character; a symbol used again in the same
  rule must match the same character, and ``"%`` creates the character ``%``
  matched. ``*`` matches the shortest possible substring (the "next nearest"
  occurrence of what follows).
- ``/A/B/`` matches A followed by B and writes B followed by A.
- ``&`` and ``$`` separate rules, which run in order on the operand; after a
  rule, ``&`` goes on only if it matched, ``$`` only if it did not.
- A rule's pattern (kept and deleted symbols) is searched leftmost in the
  operand and only that occurrence is rewritten.

With these semantics the published ancestral cells of the three models work
as designed: copiers make exact copies of genes (``0000 genotype 0011``) and
constructors translate genes into phenotypes by removing the suppressors.

method "closure" is the reaction closure of the ancestral cell (or of given
strings) with chemart.expand.expand; method "soup" is a well-stirred run
inside one cell with chemart.soup.soup. Cell division, age counters and the
cell population live above a single network and are only recorded in extras.
"""

from __future__ import annotations

import re
from collections import Counter
from functools import lru_cache

from chemart.expand import expand
from chemart.network import Network, Reaction, Species
from chemart.soup import soup

ALPHABET = ".0123'\"!?%*/\\&$MELR"
_ONE = "!?%"
# Bijective, space-free transliteration used for species ids.
_ID = {".": "t", "'": "d", '"': "c", "!": "x", "?": "y", "%": "z", "*": "s",
       "/": "p", "\\": "b", "&": "a", "$": "n"}

# --- the ancestral cells (ECAL 2003 slides 16-18) -----------------------------
# Phenotypic enzymes; the genes are built from them (slide 15: a genotypic
# string inserts '\' before every character of its phenotype).
MODEL_I_COPIER = [
    r"""00'0'0"1"1"0"2"0"2*0011"0"2"0"2""",
    r"""/\0\2\0\2/\\%/*"\\"%0202""",
    r"""00'1'1"0"0*'0'2'0'20011"\."0"0"0"0*'0'2'0'2"0"0"1"1""",
]
MODEL_I_CONSTRUCTOR = [
    r"""00'0'0"1"1"0"2"2"0*0011"0"2"2"0""",
    r"""/\0\2\2\0/\\%/*"%0220""",
    r"""00'1'1"0"0*'0'2'2'00011"\.*'0'2'2'0""",
]
MODEL_II_COPIER = [
    r"""\0\0'0'0"1"1"0"2"0"2*3300"0"2"0"2""",
    r"""/\0\2\0\2/\\%/*"\\"%\0\2\0\2""",
    r"""/\0\2\0\2/!/*"!\0\2\0\2""",
    r"""\0\0'1'1"0"0*'0'2'0'23300"\."0"0"0"0*'0'2'0'2"3"3"0"0""",
]
MODEL_II_CONSTRUCTOR = [
    r"""\0\0'0'0"1"1"0"2"2"0*3300"0"2"2"0""",
    r"""/\0\2\2\0/\\%/*"%\0\2\2\0""",
    r"""/\0\2\2\0/0101*0110/*3300/"\.*/\0\2\2\0/""",
    r"""\0\0'1'1"0"0*'0'2'2'03300"\.*'0'2'2'0""",
]
MODEL_III_COPIER = [
    r"""\L\0\0$\R\0\0$\0\0"1"1\0"2\0"2*"0"2"0"2"\&\0\0\1\1"\."\'"\'"\\"\""\""\\"\'"\""\\"\'"\M"\\"\\"\\"\E"\\"\\"\\"\M""",
    r"""/\0\2\0\2/\\%/*"\\"%\0\0\1\1""",
    r"""\L\0\0$\R\0\0$\0\0*\0\0\1\1"\."\'"\'"\M"\""\\"\'"\""\\"\'"\M"\M"\M"\M""",
    r"""
"L*\0\0'1'1"0"0*'0'2'0'2\0'2\0'2"1"1'\&"\."R"0"0"0"0*\0\0\1\1"\."\E"\M"\."\""\M"\\"\E"\\"\M""".strip(),
]
MODEL_III_CONSTRUCTOR = [
    r"""\0\0'0'0"1"1"0"2"2"0*"0"2"2"0"\&\0\0\1\1""",
    r"""/\0\2\2\0/\\%/*"%\0\0\1\1""",
    r"""\0\0'1'1"0"0*'0'2'2'0\0'2'2\0"1"1'\&"\.*'0'0'1'1""",
]
CHROMOSOME_SEPARATOR = "01010110"
MEMBRANE_SEED = "EM"
MEMBRANE_BREEDER = '"M\\E\\M'                 # decodes to [EM → MEM]
MEMBRANE_DIVISION_M = 10

ORGANISATIONS = ("independent-genes", "single-chromosome", "spindle-membrane")


def escape(phenotype: str) -> str:
    """The genotypic form of a string: '\\' before every character (slide 15)."""
    return "".join("\\" + c for c in phenotype)


def translate(genotype: str) -> str:
    """Remove one level of suppression (what the constructor does to a gene body)."""
    out, i = [], 0
    while i < len(genotype):
        if genotype[i] == "\\" and i + 1 < len(genotype):
            out.append(genotype[i + 1])
            i += 2
        else:
            out.append(genotype[i])
            i += 1
    return "".join(out)


def ancestral_cell(organisation: str) -> dict[str, list[str]]:
    """The published string set of an ancestral cell, by role."""
    if organisation == "independent-genes":
        enzymes = MODEL_I_COPIER + MODEL_I_CONSTRUCTOR
        return {"genes": ["0000" + escape(p) + "0011" for p in enzymes],
                "copier": list(MODEL_I_COPIER), "constructor": list(MODEL_I_CONSTRUCTOR)}
    if organisation == "single-chromosome":
        enzymes = MODEL_II_COPIER + MODEL_II_CONSTRUCTOR
        chromosome = "0000" + CHROMOSOME_SEPARATOR.join(escape(p) for p in enzymes) + "3300"
        return {"genes": [chromosome], "copier": list(MODEL_II_COPIER),
                "constructor": list(MODEL_II_CONSTRUCTOR)}
    if organisation == "spindle-membrane":
        enzymes = MODEL_III_COPIER + MODEL_III_CONSTRUCTOR
        return {"genes": ["0000" + escape(p) + "0011" for p in enzymes],
                "copier": list(MODEL_III_COPIER), "constructor": list(MODEL_III_CONSTRUCTOR)}
    raise ValueError(f"organisation must be one of {ORGANISATIONS}, got {organisation!r}")


# --- decoding -------------------------------------------------------------------
def _tokenize(s: str) -> list[tuple]:
    """Tokens (kind, symbol, mode, start, end); kind lit|one|star|sep|slash, mode keep|del|new."""
    toks, i, n = [], 0, len(s)
    while i < n:
        start, c, mode = i, s[i], "keep"
        if c in "'\"":
            mode = "del" if c == "'" else "new"
            i += 1
            if i >= n:
                break                     # dangling prefix: ignored
            c = s[i]
        if c == "\\":
            if i + 1 >= n:
                break                     # dangling suppressor: ignored
            toks.append(("lit", s[i + 1], mode, start, i + 2))
            i += 2
            continue
        if c in _ONE:
            kind = "one"
        elif c == "*":
            kind = "star"
        elif mode == "keep" and c in "&$":
            kind = "sep"
        elif mode == "keep" and c == "/":
            kind = "slash"
        else:
            kind = "lit"                  # includes a prefixed ', ", /, & or $
        toks.append((kind, c, mode, start, i + 1))
        i += 1
    return toks


class Rule:
    """One guarded rewriting rule, compiled to a regular expression."""

    def __init__(self, toks: list[tuple], source: str):
        self.source = source
        slashes = [i for i, t in enumerate(toks) if t[0] == "slash"]
        usable = slashes[: len(slashes) - len(slashes) % 3]   # leftover slashes are ignored
        items, prev = [], 0
        for g in range(0, len(usable), 3):
            a, b, c = usable[g:g + 3]
            items += [("tok", t) for t in toks[prev:a] if t[0] != "slash"]
            items.append(("perm", toks[a + 1:b], toks[b + 1:c]))
            prev = c + 1
        items += [("tok", t) for t in toks[prev:] if t[0] != "slash"]

        seq, out_order, uid = [], [], 0          # (token, uid): match order and write order
        for item in items:
            if item[0] == "perm":
                A = [(t, uid + j) for j, t in enumerate(item[1])]
                uid += len(A)
                B = [(t, uid + j) for j, t in enumerate(item[2])]
                uid += len(B)
                seq += A + B
                out_order += B + A
            else:
                seq.append((item[1], uid))
                out_order.append((item[1], uid))
                uid += 1

        parts, group_of, bound, stars, g = [], {}, {}, [], 0
        for (kind, sym, mode, *_), u in seq:
            if mode == "new":
                continue
            g += 1
            group_of[u] = g
            if kind == "lit":
                parts.append("(" + re.escape(sym) + ")")
            elif kind == "star":
                parts.append("(.*?)")
                stars.append((u, g))
            elif sym in bound:
                parts.append(f"(\\{bound[sym]})")
            else:
                bound[sym] = g
                parts.append("(.)")
        self.regex = re.compile("".join(parts), re.S)

        self.valid, self.out = True, []
        for (kind, sym, mode, *_), u in out_order:
            if mode == "keep":
                self.out.append((True, group_of[u]))
            elif mode == "new":
                if kind == "lit":
                    self.out.append((False, sym))
                elif kind == "one" and sym in bound:
                    self.out.append((True, bound[sym]))
                elif kind == "star" and any(su < u for su, _ in stars):
                    self.out.append((True, [sg for su, sg in stars if su < u][-1]))
                else:
                    self.valid = False            # creates an unbound wild card
        self.rewrites = any(t[2] != "keep" for t, _ in seq) or [u for _, u in seq] != [u for _, u in out_order]

        def show(tok):
            return tok[1] if tok[0] == "lit" or tok[0] in ("one", "star") else ""
        self.left = "".join(show(t) for t, _ in seq if t[2] != "new")
        self.right = "".join(show(t) for t, _ in out_order if t[2] != "del")

    def apply(self, s: str) -> tuple[bool, str]:
        if not self.valid:
            return False, s
        m = self.regex.search(s)
        if m is None:
            return False, s
        body = "".join(m.group(v) if is_group else v for is_group, v in self.out)
        return True, s[:m.start()] + body + s[m.end():]

    def text(self) -> str:
        return f"[{self.left} → {self.right}]" if self.rewrites else f"[{self.source}]"


@lru_cache(maxsize=65536)
def program(operator: str) -> tuple[tuple[Rule, str | None], ...]:
    """Decode a string into its rules, each with the separator that follows it."""
    toks = _tokenize(operator)
    rules, cur = [], []
    for t in toks:
        if t[0] == "sep":
            rules.append((cur, t[1]))
            cur = []
        else:
            cur.append(t)
    rules.append((cur, None))
    out = []
    for body, sep in rules:
        source = operator[body[0][3]:body[-1][4]] if body else ""
        out.append((Rule(body, source), sep))
    return tuple(out)


def decode(operator: str) -> str:
    """Book notation, e.g. [\\0]&[\\!]$[00*01 → 01*02] (book eq. 11.2)."""
    return "".join(rule.text() + (sep or "") for rule, sep in program(operator))


def rewrite(operator: str, operand: str) -> str:
    """Apply the operator's rules to the operand (before cutting at '.')."""
    for rule, sep in program(operator):
        matched, operand = rule.apply(operand)
        if sep is None or (sep == "&" and not matched) or (sep == "$" and matched):
            break
    return operand


def split(s: str) -> list[str]:
    """Cut a string at every unsuppressed '.'; empty pieces vanish."""
    pieces, cur, i = [], [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s):
            cur.append(s[i:i + 2])
            i += 2
            continue
        if s[i] == ".":
            pieces.append("".join(cur))
            cur = []
        else:
            cur.append(s[i])
        i += 1
    pieces.append("".join(cur))
    return [p for p in pieces if p]


@lru_cache(maxsize=262144)
def react(operator: str, operand: str) -> tuple[str, ...] | None:
    """operator + operand -> (operator, *pieces), or None if elastic."""
    out = rewrite(operator, operand)
    if out == operand:
        return None
    return (operator, *split(out))


def sid(s: str) -> str:
    return "s_" + "".join(_ID.get(c, c) for c in s)


def membrane_m(s: str) -> int:
    """M symbols before the final EM of a membrane string M...MEM (slide: MMMMMMMMMMEM has ten), else 0."""
    return len(s) - 2 if re.fullmatch(r"M*EM", s) else 0


# --- generator ------------------------------------------------------------------
def _seed(p) -> tuple[list[str], dict[str, str]]:
    if p.strings:
        bad = [s for s in p.strings if not isinstance(s, str) or not s or set(s) - set(ALPHABET)]
        if bad:
            raise ValueError(f"strings must be non-empty SAC strings over {ALPHABET!r}, got {bad[:3]!r}")
        cut = [s for s in p.strings if split(s) != [s]]
        if cut:
            raise ValueError(f"strings must not contain an unsuppressed '.' (write '\\.'), got {cut[:3]!r}")
        return list(p.strings), {}
    cell = ancestral_cell(p.organisation)
    roles = {}
    for role, strings in cell.items():
        for s in strings:
            roles[s] = role.rstrip("s") if role == "genes" else role
    strings = cell["genes"] + cell["copier"] + cell["constructor"]
    if p.organisation == "spindle-membrane":
        strings.append(MEMBRANE_SEED)
        roles[MEMBRANE_SEED] = "membrane"
    return strings, roles


def _species(strings) -> list[Species]:
    return [Species(sid(s), structure=s) for s in strings]


def _extras(p, seed, roles) -> dict:
    return {
        "organisation": None if p.strings else p.organisation,
        "seed": [sid(s) for s in dict.fromkeys(seed)],
        "roles": {sid(s): r for s, r in roles.items()},
        "cell_level": {
            "note": "not part of the network: cells, division and selection act on whole string sets",
            "division": "a cell divides when its number of strings has doubled; its strings are shared "
                        "equally between the daughters (L/R strings forcibly left/right in model iii)",
            "membrane_condition": f"model iii: a cell cannot divide without a membrane string of "
                                  f"{MEMBRANE_DIVISION_M} M symbols (MMMMMMMMMMEM)",
            "population": "max 100 cells with max 50 strings; when a new cell is added an old one is removed; "
                          "cells carry an age counter reset at division (book 11.1.3)",
        },
    }


def generate(p, rng):
    seed, roles = _seed(p)
    if p.method == "closure":
        return _closure(p, seed, roles)
    return _soup(p, rng, seed, roles)


def _closure(p, seed, roles):
    strings, found, status = expand(react, list(dict.fromkeys(seed)), arity=2,
                                    max_species=p.max_species, ordered=True)
    reactions = [Reaction.of([sid(s) for s in lhs], [sid(s) for s in rhs]) for lhs, rhs in found]
    produced = set()
    for lhs, rhs in found:
        produced.update((Counter(rhs) - Counter(lhs)).elements())
    extras = _extras(p, seed, roles)
    extras["analysis"] = {
        "reproduced_seed": [sid(s) for s in dict.fromkeys(seed) if s in produced],
        "max_membrane_M": max((membrane_m(s) for s in strings), default=0),
    }
    return Network(
        species=_species(strings),
        reactions=reactions,
        status=status,
        initial_state={sid(s): float(c) for s, c in Counter(seed).items()},
        extras=extras,
    )


def _soup(p, rng, seed, roles):
    start = [s for s in seed for _ in range(p.copies)]
    if len(start) < 2:
        raise ValueError(f"the soup needs at least 2 strings, got {len(start)}")
    pop, fired, seen = list(start), {}, dict.fromkeys(start)
    size, trace, done = [len(pop)], [], 0
    chunk = max(1, len(start))
    while done < p.steps and len(pop) >= 2:
        steps = min(chunk, p.steps - done)
        part, pop = soup(react, pop, steps, rng, arity=2, dilution=p.dilution)
        for lhs, rhs, count in part:
            key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
            fired.setdefault(key, [lhs, rhs, 0])[2] += count
            seen.update(dict.fromkeys(rhs))
        done += steps
        size.append(len(pop))
        trace.append(max((membrane_m(s) for s in pop), default=0))
    final = Counter(pop)
    extras = _extras(p, seed, roles)
    extras["analysis"] = {
        "chunk_steps": chunk,
        "population_size": size,
        "max_membrane_M": trace,
        "seed_copies_final": {sid(s): final.get(s, 0) for s in dict.fromkeys(seed)},
    }
    extras["final_state"] = {sid(s): c for s, c in final.most_common()}
    return Network(
        species=_species(seen),
        reactions=[Reaction.of([sid(s) for s in lhs], [sid(s) for s in rhs], count=c)
                   for lhs, rhs, c in fired.values()],
        status="observed",
        initial_state={sid(s): float(c) for s, c in Counter(start).items()},
        extras=extras,
    )
