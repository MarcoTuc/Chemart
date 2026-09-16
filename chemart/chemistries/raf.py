"""RAF sets: reflexively autocatalytic and F-generated (book 6.3.1). Catalog id: raf.

An ANALYSIS. `generate` builds a catalytic reaction system - Kauffman's binary
polymer model (all strings up to length n over an alphabet of size B, every
cleavage/ligation pair, each molecule catalysing each pair with probability
p = f / |R|) with the food set of all strings up to length t - and then runs the
RAF algorithm on it. `system` replaces the polymer model by a user-supplied CRS.

Definitions (Hordijk, Kauffman & Steel 2011, sec. 3.1; Steel, Hordijk & Smith
2013, sec. 2.1). Given a CRS Q = (X, R, C) and a food set F, a non-empty subset
R' of R is

- reflexively autocatalytic (RA) if every r in R' has a catalyst in cl_R'(F);
- F-generated if every reactant of every r in R' is in cl_R'(F);
- an RAF if it is both,

where cl_R'(F), the closure of F relative to R', is F plus everything that can
be built from F by repeatedly applying reactions of R'. The RAF algorithm
(ibid., sec. 3.2) iterates R_{i+1} = {r in R_i : all reactants and at least one
catalyst of r lie in cl_{R_i}(F)} to its fixed point, which is empty (no RAF) or
the unique maximal RAF. An irrRAF is an RAF no proper subset of which is an RAF;
one is found by trying to delete the reactions of the maxRAF in random order.

Cleavage and ligation are one reversible reaction, as in the model papers, so a
reaction's molecules are its reactants *and* its products.
"""

from itertools import product

from chemart.helpers.explicit import network, term
from chemart.network import Species

ALPHABET = "0123"
MAX_SPECIES = 40000
MAX_REACTIONS = 200000


# --------------------------------------------------------------------------
# The RAF algorithm, on reactions given as
# {"id", "reactants": [...], "products": [...], "catalysts": [...], "reversible": bool}
# --------------------------------------------------------------------------
def _directions(r: dict) -> list[tuple[list[str], list[str]]]:
    out = [(r["reactants"], r["products"])]
    if r["reversible"]:
        out.append((r["products"], r["reactants"]))
    return out


def molecules_of(r: dict) -> set[str]:
    """rho(r): the molecules that must be constructible for r to be usable."""
    return set(r["reactants"]) | (set(r["products"]) if r["reversible"] else set())


def closure(food, reactions) -> set[str]:
    """cl_R(F): F plus everything R can build from it, catalysts ignored."""
    needs: dict[str, list[int]] = {}
    missing: list[int] = []
    makes: list[list[str]] = []
    for r in reactions:
        for lhs, rhs in _directions(r):
            k = len(makes)
            makes.append(rhs)
            missing.append(0)
            for x in set(lhs):
                needs.setdefault(x, []).append(k)
                missing[k] += 1
    cl = set()
    queue = []
    for x in food:
        if x not in cl:
            cl.add(x)
            queue.append(x)
    while queue:
        x = queue.pop()
        for k in needs.get(x, ()):
            missing[k] -= 1
            if missing[k] == 0:
                for y in makes[k]:
                    if y not in cl:
                        cl.add(y)
                        queue.append(y)
    return cl


def max_raf(reactions, food) -> list[dict]:
    """The RAF algorithm: the unique maximal RAF within `reactions`, or []."""
    current = list(reactions)
    while current:
        cl = closure(food, current)
        keep = [r for r in current
                if molecules_of(r) <= cl and any(c in cl for c in r["catalysts"])]
        if len(keep) == len(current):
            return keep
        current = keep
    return []


def is_raf(reactions, food) -> bool:
    """Is this exact set of reactions an RAF (RA and F-generated)?"""
    if not reactions:
        return False
    cl = closure(food, reactions)
    return all(molecules_of(r) <= cl and any(c in cl for c in r["catalysts"])
               for r in reactions)


def irreducible_raf(reactions, food, rng) -> list[dict]:
    """One irrRAF: delete reactions of the maxRAF in random order, keeping every
    deletion whose remaining maxRAF is non-empty (Hordijk & Steel 2004)."""
    current = max_raf(reactions, food)
    if not current:
        return []
    order = [current[i] for i in rng.permutation(len(current))]
    for r in order:
        if r["id"] not in {x["id"] for x in current}:
            continue
        smaller = max_raf([x for x in current if x["id"] != r["id"]], food)
        if smaller:
            current = smaller
    return current


# --------------------------------------------------------------------------
# The two chemical reaction systems
# --------------------------------------------------------------------------
def binary_polymer_model(p, rng) -> tuple[list[str], list[dict], list[str], float]:
    """Kauffman's binary polymer model: species, cleavage/ligation reactions with
    catalysts drawn at random, food set, and the catalysis probability used."""
    letters = ALPHABET[: p.B]
    species = ["".join(s) for n in range(1, p.n + 1) for s in product(letters, repeat=n)]
    if len(species) > MAX_SPECIES:
        raise ValueError(
            f"B = {p.B}, n = {p.n} gives {len(species)} polymers, above the limit of {MAX_SPECIES}"
        )
    if p.t > p.n:
        raise ValueError(f"the food set length t = {p.t} must not exceed n = {p.n}")
    food = [s for s in species if len(s) <= p.t]

    reactions = [
        {"id": f"{c[:cut]}+{c[cut:]}<->{c}", "reactants": [c[:cut], c[cut:]],
         "products": [c], "catalysts": [], "reversible": True}
        for c in species for cut in range(1, len(c))
    ]
    prob = p.f / len(reactions)
    if prob > 1:
        raise ValueError(
            f"f = {p.f} needs a catalysis probability of {prob:.3g} > 1 with only "
            f"{len(reactions)} reactions; lower f or raise n"
        )

    # Bernoulli(prob) over all (molecule, reaction) pairs: draw the number of
    # catalysis events, then a uniform subset of the pairs of that size.
    total = len(species) * len(reactions)
    count = int(rng.binomial(total, prob))
    chosen: set[int] = set()
    while len(chosen) < count:
        chosen.update(int(v) for v in rng.integers(0, total, size=count - len(chosen)))
    for pair in chosen:
        reactions[pair % len(reactions)]["catalysts"].append(species[pair // len(reactions)])
    for r in reactions:
        r["catalysts"].sort()
    return species, reactions, food, prob


def given_system(system) -> tuple[list[str], list[dict], list[str]]:
    """Validate a user-supplied CRS: {"food": [...], "reactions": [{...}, ...]}."""
    if set(system) - {"food", "reactions"}:
        raise ValueError(f"system takes only 'food' and 'reactions', got {sorted(system)}")
    raw = system.get("reactions")
    if not isinstance(raw, list) or not raw:
        raise ValueError("system['reactions'] must be a non-empty list of reactions")
    food = system.get("food", [])
    if not isinstance(food, list) or not all(isinstance(x, str) and x and " " not in x for x in food):
        raise ValueError(f"system['food'] must be a list of species ids, got {food!r}")

    reactions, order = [], list(food)
    for i, item in enumerate(raw):
        if not isinstance(item, dict) or set(item) - {"id", "reactants", "products", "catalysts", "reversible"}:
            raise ValueError(
                f"system['reactions'][{i}] must be an object with keys id, reactants, "
                f"products, catalysts, reversible; got {item!r}"
            )
        r = {
            "id": str(item.get("id", f"r{i + 1}")),
            "reactants": list(item.get("reactants", [])),
            "products": list(item.get("products", [])),
            "catalysts": list(item.get("catalysts", [])),
            "reversible": bool(item.get("reversible", False)),
        }
        for key in ("reactants", "products", "catalysts"):
            bad = [x for x in r[key] if not isinstance(x, str) or not x or " " in x]
            if bad:
                raise ValueError(f"system['reactions'][{i}]['{key}'] has invalid species ids {bad}")
        if not r["reactants"] and not r["products"]:
            raise ValueError(f"system['reactions'][{i}] is empty on both sides")
        order += r["reactants"] + r["products"] + r["catalysts"]
        reactions.append(r)
    if len({r["id"] for r in reactions}) != len(reactions):
        raise ValueError("system reaction ids must be unique")

    seen: dict[str, None] = {}
    for x in order:
        seen.setdefault(x, None)
    return list(seen), reactions, list(food)


# --------------------------------------------------------------------------
def _side(names: list[str]) -> str:
    counts: dict[str, int] = {}
    for x in names:
        counts[x] = counts.get(x, 0) + 1
    return " + ".join(term(n, x) for x, n in counts.items())


def _texts(r: dict, catalysed_only: bool) -> list[str]:
    """The reaction(s) this CRS reaction contributes to the network: one per
    catalyst (catalyst on both sides) and, unless only catalysed reactions are
    asked for, the bare reaction itself."""
    out = []
    sides = _directions(r)
    if not catalysed_only:
        out += [f"{_side(lhs)} -> {_side(rhs)}" for lhs, rhs in sides]
    for e in r["catalysts"]:
        out += [f"{_side(lhs + [e])} -> {_side(rhs + [e])}" for lhs, rhs in sides]
    return out


def generate(p, rng):
    if p.system:
        species, reactions, food = given_system(p.system)
        prob = None
    else:
        species, reactions, food, prob = binary_polymer_model(p, rng)

    raf = max_raf(reactions, food)
    cl = closure(food, raf) if raf else set(food)
    irr = []
    for _ in range(p.irreducible_rafs if raf else 0):
        one = irreducible_raf(raf, food, rng)
        irr.append({"size": len(one), "reactions": sorted(r["id"] for r in one)})

    catalysed = [r for r in reactions if r["catalysts"]]
    texts = [t for r in reactions for t in _texts(r, p.reactions == "catalysed")]
    if len(texts) > MAX_REACTIONS:
        raise ValueError(
            f"this network would have {len(texts)} reactions, above the limit of "
            f"{MAX_REACTIONS}; lower n or f, or set reactions='catalysed'"
        )

    analysis = {
        "raf_exists": bool(raf),
        "max_raf_size": len(raf),
        "max_raf_reactions": sorted(r["id"] for r in raf),
        "max_raf_molecules": sorted({x for r in raf for x in molecules_of(r)}),
        "max_raf_closure_size": len(cl),
        "irreducible_rafs": irr,
        "reactions_total": len(reactions),
        "catalysed_reactions": len(catalysed),
        "catalysis_events": sum(len(r["catalysts"]) for r in reactions),
        "food": list(food),
    }
    if prob is not None:
        analysis["level_of_catalysis"] = float(p.f)
        analysis["catalysis_probability"] = float(prob)

    extras = {"analysis": analysis, "food": list(food), "buffered": list(food)}
    if not p.system:
        letters = ALPHABET[: p.B]
        extras["conservation"] = [
            {"name": f"monomer {x}", "vector": {s: s.count(x) for s in species if x in s}}
            for x in letters
        ]
    return network(
        [(t, None) for t in texts],
        species=[Species(s, structure=s) if not p.system else Species(s) for s in species],
        extras=extras,
    )
