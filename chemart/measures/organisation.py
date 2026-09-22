"""D. Catalysis, autocatalysis and organisation.

RAF sets (Hordijk & Steel 2004; Hordijk, Kauffman & Steel 2011): a set of
reactions is reflexively autocatalytic and food-generated when every reaction
has a catalyst, and all its reactants, in the closure of the food set under
the set itself. The RAF algorithm finds the unique maximal one.

Network expansion (Handorf, Ebenhöh & Heinrich 2005): the scope of a food set
is everything the network can make from it, reaction by reaction.

The RAF functions work on reactions given as
{"id", "reactants": [...], "products": [...], "catalysts": [...], "reversible": bool},
the form chemart.chemistries.raf uses; `raf_reactions` converts a network.
"""

from __future__ import annotations

from chemart.measures import register


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


def raf_reactions(net) -> list[dict]:
    """A network's reactions in RAF form: a species on both sides is a catalyst."""
    out = []
    for j, r in enumerate(net.reactions):
        cat = r.catalysts
        out.append({
            "id": j,
            "reactants": [s for s, n in r.reactants.items() if n > cat.get(s, 0)],
            "products": [s for s, n in r.products.items() if n > cat.get(s, 0)],
            "catalysts": sorted(cat),
            "reversible": False,
        })
    return out


@register("max_raf_fraction", "D", needs="C F")
def max_raf_fraction(ctx) -> float:
    """Share of the reactions that belong to the maximal RAF set."""
    return len(max_raf(raf_reactions(ctx.net), ctx.food)) / len(ctx.net.reactions)


def expansion(net, food) -> list[set[str]]:
    """Network expansion from a food set: the species available after each
    generation. As in the RAF closure, catalysts are not needed for a reaction
    to fire, only its consumed reactants (a species on both sides is a catalyst)."""
    reactions = raf_reactions(net)
    available = set(food)
    layers = [set(available)]
    while True:
        new = {s for r in reactions if set(r["reactants"]) <= available for s in r["products"]} - available
        if not new:
            return layers
        available |= new
        layers.append(set(available))


@register("scope_fraction", "D", needs="S F")
def scope_fraction(ctx) -> float:
    """Share of species the network can make from its food set (its scope),
    catalysts not required."""
    return len(expansion(ctx.net, ctx.food)[-1]) / len(ctx.ids)


@register("expansion_depth", "D", needs="S F")
def expansion_depth(ctx) -> int:
    """Generations network expansion takes to reach the scope of the food set."""
    return len(expansion(ctx.net, ctx.food)) - 1
