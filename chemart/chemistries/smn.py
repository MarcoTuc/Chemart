"""String Metabolic Network (SMN) of Ono, Fujiwara & Yuta (2005). Catalog id: smn.

Built from Banzhaf & Yamamoto (2015), section 18.3.1 (the original paper,
LNAI 3630 pp. 716-724, is not openly available):

- metabolites are strings over the alphabet {a, ..., h};
- an organism's genome is a list of enzymes, and each enzyme *is* one
  reversible reaction: a ligation A + B <-> AB (cleavage backwards, eq. 18.17)
  or a recombination AB + CD <-> AD + CB (eq. 18.18);
- the metabolism of one organism is the network of its enzymes' reactions over
  its compound set (initial metabolites plus everything the enzymes make);
- the book's mutation operator (duplicate an enzyme, replace one educt with a
  compound already in the compound set, choose a new recombination point, add
  the result to the genome) can be applied `mutations` times, without
  selection;
- a regulated inflow/outflow keeps the amount of substance in the cell
  constant (outflow = constant-total), so the letter mass, the fitness, grows
  only when longer compounds are made.

Enzymes are reactions, not molecules, so they are not species: each genome
entry in extras.genome names the reactions it encodes. The book gives no
kinetic coefficients (only that longer products are harder to make), so rates
are None. Selection across generations is not part of one network.
"""

from collections import Counter

from chemart.network import CONSTANT_TOTAL, Network, Reaction, Species

LETTERS = "abcdefgh"
ATTEMPTS = 1000


# ----------------------------------------------------------------------------
# Enzymes: {"kind": "ligation" | "recombination", "educts": [x, y], "cuts": [i, j] | None}
# ----------------------------------------------------------------------------
def products(enzyme: dict) -> list[str]:
    x, y = enzyme["educts"]
    if enzyme["kind"] == "ligation":
        return [x + y]
    i, j = enzyme["cuts"]
    return [x[:i] + y[j:], y[:j] + x[i:]]        # A|B + C|D -> AD + CB


def notation(enzyme: dict) -> str:
    """'A + B' for a ligation, 'A|B + C|D' for a recombination."""
    x, y = enzyme["educts"]
    if enzyme["kind"] == "ligation":
        return f"{x} + {y}"
    i, j = enzyme["cuts"]
    return f"{x[:i]}|{x[i:]} + {y[:j]}|{y[j:]}"


def _elastic(enzyme: dict) -> bool:
    return Counter(enzyme["educts"]) == Counter(products(enzyme))


def _fits(enzyme: dict, max_length: int) -> bool:
    return all(len(s) <= max_length for s in products(enzyme))


def _parse(text: str, letters: str) -> dict:
    if not isinstance(text, str):
        raise ValueError(f"genome entries must be strings like 'ab + cd' or 'a|b + cd|e', got {text!r}")
    parts = [t.strip() for t in text.split(" + ")]
    if len(parts) != 2:
        raise ValueError(f"genome: cannot parse {text!r}; expected 'A + B' (ligation) or 'A|B + C|D' (recombination)")
    bars = [t.count("|") for t in parts]
    if bars == [0, 0]:
        enzyme = {"kind": "ligation", "educts": parts, "cuts": None}
    elif bars == [1, 1]:
        halves = [t.split("|") for t in parts]
        if not all(h[0] and h[1] for h in halves):
            raise ValueError(f"genome: {text!r} has an empty fragment; both sides of each '|' must be non-empty")
        enzyme = {"kind": "recombination", "educts": [h[0] + h[1] for h in halves],
                  "cuts": [len(h[0]) for h in halves]}
    else:
        raise ValueError(f"genome: {text!r} must mark a recombination point in both educts ('A|B + C|D') or in neither ('A + B')")
    for s in enzyme["educts"]:
        if not s or set(s) - set(letters):
            raise ValueError(f"genome: {text!r} uses letters outside the alphabet {letters!r}")
    if _elastic(enzyme):
        raise ValueError(f"genome: {text!r} gives back its own educts (an elastic reaction)")
    return enzyme


# ----------------------------------------------------------------------------
def _draw_enzyme(rng, compounds, p):
    """A random enzyme acting on the current compound set."""
    long = [c for c in compounds if len(c) >= 2]
    for _ in range(ATTEMPTS):
        if rng.random() < p.ligation_fraction or not long:
            x, y = (compounds[int(k)] for k in rng.integers(len(compounds), size=2))
            enzyme = {"kind": "ligation", "educts": [x, y], "cuts": None}
        else:
            x, y = (long[int(k)] for k in rng.integers(len(long), size=2))
            cuts = [int(rng.integers(1, len(x))), int(rng.integers(1, len(y)))]
            enzyme = {"kind": "recombination", "educts": [x, y], "cuts": cuts}
        if _fits(enzyme, p.max_length) and not _elastic(enzyme):
            return enzyme
    raise ValueError(f"could not draw an enzyme with products of length <= max_length = {p.max_length}; increase max_length or use shorter metabolites")


def _mutate(rng, genome, compounds, p):
    """Book 18.3.1: duplicate a random enzyme, swap one educt for a known compound, re-cut."""
    for _ in range(ATTEMPTS):
        parent = int(rng.integers(len(genome)))
        old = genome[parent]
        slot = int(rng.integers(2))
        pool = [c for c in compounds if c != old["educts"][slot]
                and (old["kind"] == "ligation" or len(c) >= 2)]
        if not pool:
            continue
        new = pool[int(rng.integers(len(pool)))]
        educts = list(old["educts"])
        educts[slot] = new
        cuts = None
        if old["kind"] == "recombination":
            cuts = list(old["cuts"])
            cuts[slot] = int(rng.integers(1, len(new)))
        enzyme = {"kind": old["kind"], "educts": educts, "cuts": cuts, "parent": parent}
        if _fits(enzyme, p.max_length) and not _elastic(enzyme):
            return enzyme
    raise ValueError(f"could not find a valid mutation with products of length <= max_length = {p.max_length}; increase max_length")


def _random_metabolites(rng, p, letters):
    possible = sum(len(letters) ** n for n in range(1, p.max_metabolite_length + 1))
    if p.n_metabolites > possible:
        raise ValueError(f"n_metabolites = {p.n_metabolites} exceeds the {possible} distinct strings of length <= {p.max_metabolite_length} over {letters!r}")
    if p.max_metabolite_length > p.max_length:
        raise ValueError(f"max_metabolite_length = {p.max_metabolite_length} must be <= max_length = {p.max_length}")
    out: dict[str, None] = {}
    while len(out) < p.n_metabolites:
        n = int(rng.integers(1, p.max_metabolite_length + 1))
        out.setdefault("".join(letters[int(k)] for k in rng.integers(len(letters), size=n)))
    return list(out)


def generate(p, rng):
    letters = LETTERS[: p.alphabet_size]
    explicit = [_parse(t, letters) for t in p.genome]
    for e, text in zip(explicit, p.genome):
        if not _fits(e, p.max_length):
            raise ValueError(f"genome: {text!r} makes a product longer than max_length = {p.max_length}")

    if p.initial_metabolites:
        bad = [s for s in p.initial_metabolites
               if not isinstance(s, str) or not s or set(s) - set(letters) or len(s) > p.max_length]
        if bad:
            raise ValueError(f"initial_metabolites {bad} must be non-empty strings over {letters!r} of length <= max_length = {p.max_length}")
        initial = list(dict.fromkeys(p.initial_metabolites))
    elif explicit:
        initial = list(dict.fromkeys(s for e in explicit for s in e["educts"]))
    else:
        initial = _random_metabolites(rng, p, letters)

    compounds = list(initial)
    known = set(compounds)

    def add(enzyme):
        enzyme.setdefault("parent", None)
        genome.append(enzyme)
        for s in (*enzyme["educts"], *products(enzyme)):
            if s not in known:
                known.add(s)
                compounds.append(s)

    genome: list[dict] = []
    if explicit:
        for e in explicit:
            add(e)
    else:
        for _ in range(p.n_enzymes):
            add(_draw_enzyme(rng, compounds, p))
    if p.mutations and not genome:
        raise ValueError("mutations need a non-empty genome: set n_enzymes >= 1 or give a genome")
    for _ in range(p.mutations):
        add(_mutate(rng, genome, compounds, p))

    # Each enzyme encodes a reversible reaction: two directed reactions, shared
    # between enzymes that catalyse the same reaction.
    reactions: list[Reaction] = []
    index: dict[tuple, int] = {}

    def directed(lhs, rhs):
        key = (frozenset(Counter(lhs).items()), frozenset(Counter(rhs).items()))
        if key not in index:
            index[key] = len(reactions)
            reactions.append(Reaction.of(lhs, rhs))
        return index[key]

    records = []
    for n, e in enumerate(genome):
        made = products(e)
        records.append({
            "enzyme": f"E{n}",
            "kind": e["kind"],
            "educts": notation(e),
            "reaction": f"{' + '.join(e['educts'])} <-> {' + '.join(made)}",
            "reactions": [directed(e["educts"], made), directed(made, e["educts"])],
            "parent": None if e["parent"] is None else f"E{e['parent']}",
        })

    return Network(
        species=[Species(s) for s in compounds],
        reactions=reactions,
        status="complete",
        initial_state={s: p.initial_amount for s in initial},
        outflow=CONSTANT_TOTAL,
        extras={
            "genome": records,
            "initial_metabolites": initial,
            "conservation": [
                {"name": f"letter {x}", "vector": {s: s.count(x) for s in compounds if x in s}}
                for x in letters
            ],
            "flow_law": "regulated inflow/outflow keeps the total amount of substance in the cell constant (book 18.3.1): outflow = constant-total",
            "kinetics": "not published in the book: each enzyme catalyses its reaction in both directions, and synthesising longer molecules is harder; rates are None",
            "fitness": "mass of the organism in letters, sum_s len(s) [s], after a fixed number of timesteps (book 18.3.1); used for selection across generations, not computed here",
        },
    )


def letter_mass(state: dict[str, float]) -> float:
    """The book's fitness measure: total number of letters in the organism."""
    return sum(len(s) * v for s, v in state.items())


__all__ = ["generate", "letter_mass", "notation", "products"]
