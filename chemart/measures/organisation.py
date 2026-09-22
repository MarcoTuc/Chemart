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


# --------------------------------------------------------------------------
# Exponential: irreducible RAFs, autocatalytic cores, organisations
# --------------------------------------------------------------------------

@register("irreducible_rafs", "D", needs="C F", cost="exponential", limit=3000)
def irreducible_rafs(ctx, samples: int = 20) -> int:
    """Distinct irreducible RAFs found in `samples` (20) random reduction orders
    of the maximal RAF: a lower bound on how many different ways the network
    can sustain itself (Hordijk & Steel 2004)."""
    import numpy as np

    reactions = raf_reactions(ctx.net)
    rng = np.random.default_rng(ctx.seed)
    found = {frozenset(r["id"] for r in irreducible_raf(reactions, ctx.food, rng)) for _ in range(samples)}
    found.discard(frozenset())
    return len(found)


def autocatalytic_core(net, exclude: list[set[int]] = (), cap: int = 10**6):
    """One smallest autocatalytic subnetwork (Blokhuis, Lacoste & Nghe 2020),
    as a set of reaction indices, or None: species M and reactions R such that
    every reaction of R has a reactant and a product in M, and some positive
    flux on R strictly increases every species of M. A mixed-integer programme
    minimising |R|; `exclude` rules out supersets of cores already found."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp

    ids = [s.id for s in net.species]
    index = {s: i for i, s in enumerate(ids)}
    n, r = len(ids), len(net.reactions)
    if not r:
        return None
    S = np.zeros((n, r))
    uses = np.zeros((n, r))
    makes = np.zeros((n, r))
    for j, rx in enumerate(net.reactions):
        for s, k in rx.reactants.items():
            S[index[s], j] -= k
            uses[index[s], j] = 1
        for s, k in rx.products.items():
            S[index[s], j] += k
            makes[index[s], j] = 1
    big = 1e3
    # variables: v (r, continuous 0..big), x (r, binary: reaction in R), y (n, binary: species in M)
    nv = 2 * r + n
    rows, lo, hi = [], [], []

    def add(coeffs, low, high):
        row = np.zeros(nv)
        for k, c in coeffs:
            row[k] += c
        rows.append(row)
        lo.append(low)
        hi.append(high)

    for i in range(n):          # productive: y_i = 1 -> sum_j S_ij v_j >= 1 (v = 0 off R)
        relax = big * np.abs(S[i]).sum() + 1.0      # large enough to switch the row off when y_i = 0
        add([(j, S[i, j]) for j in range(r)] + [(2 * r + i, -relax)], 1.0 - relax, np.inf)
    for j in range(r):          # v_j > 0 only on R; R's reactions touch M on both sides
        add([(j, 1.0), (r + j, -big)], -np.inf, 0.0)
        add([(r + j, 1.0)] + [(2 * r + i, -1.0) for i in range(n) if uses[i, j]], -np.inf, 0.0)
        add([(r + j, 1.0)] + [(2 * r + i, -1.0) for i in range(n) if makes[i, j]], -np.inf, 0.0)
    for i in range(n):          # a species of M must be consumed within R (else nothing is autocatalytic about it)
        add([(2 * r + i, 1.0)] + [(r + j, -1.0) for j in range(r) if uses[i, j]], -np.inf, 0.0)
    add([(2 * r + i, 1.0) for i in range(n)], 1.0, np.inf)
    for core in exclude:        # no-good cut: not every reaction of a known core
        add([(r + j, 1.0) for j in core], -np.inf, len(core) - 1)
    c = np.concatenate([np.zeros(r), np.ones(r), np.zeros(n)])
    integrality = np.concatenate([np.zeros(r), np.ones(r), np.ones(n)])
    bounds = Bounds(np.zeros(nv), np.concatenate([np.full(r, big), np.ones(r), np.ones(n)]))
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integrality,
               bounds=bounds, options={"time_limit": 30})
    if res.status != 0 or res.x is None:
        return None
    return {j for j in range(r) if res.x[r + j] > 0.5}


@register("autocatalytic_cores", "D", needs="S", cost="exponential", limit=400)
def autocatalytic_cores(ctx, cap: int = 20) -> int:
    """Number of minimal autocatalytic subnetworks, found one at a time by a
    mixed-integer programme with no-good cuts, up to `cap` (20): autocatalysis
    from stoichiometry alone, catalysts not labelled (Blokhuis, Lacoste & Nghe 2020)."""
    cores: list[set[int]] = []
    while len(cores) < cap:
        core = autocatalytic_core(ctx.net, cores)
        if core is None:
            break
        cores.append(core)
    return len(cores)


def closure_of(net, start) -> frozenset[str]:
    """The smallest closed set containing `start`: add the products of every
    reaction whose reactants are all present, until nothing changes."""
    have = set(start)
    changed = True
    while changed:
        changed = False
        for r in net.reactions:
            if set(r.reactants) <= have and not set(r.products) <= have:
                have |= set(r.products)
                changed = True
    return frozenset(have)


def self_maintaining(net, O: frozenset[str], S, index) -> bool:
    """Whether some strictly positive flux on the reactions that can run in O
    leaves no species of O decreasing (Dittrich & Speroni di Fenizio 2007)."""
    import numpy as np
    from scipy.optimize import linprog

    run = [j for j, r in enumerate(net.reactions) if set(r.reactants) <= O]
    if not run or not O:
        return True
    rows = [index[s] for s in O]
    A = S[np.ix_(rows, run)]
    res = linprog(np.zeros(len(run)), A_ub=-A, b_ub=np.zeros(len(rows)),
                  bounds=[(1.0, None)] * len(run), method="highs")
    return bool(res.status == 0)


@register("organisations", "D", needs="S", cost="exponential", limit=200)
def organisations(ctx, max_closed: int = 20000) -> dict[str, float] | None:
    """Chemical organisations (Dittrich & Speroni di Fenizio 2007): sets of
    species that are closed (make nothing outside themselves) and
    self-maintaining (can run all their reactions without depleting any member).
    Returns their number and the size of the largest as a share of all
    species; None when more than `max_closed` closed sets would need checking."""
    closed = {closure_of(ctx.net, ())}
    frontier = list(closed)
    while frontier:
        C = frontier.pop()
        for s in ctx.ids:
            if s in C:
                continue
            D = closure_of(ctx.net, C | {s})
            if D not in closed:
                closed.add(D)
                frontier.append(D)
                if len(closed) > max_closed:
                    return None
    orgs = [O for O in closed if self_maintaining(ctx.net, O, ctx.S, ctx.index)]
    return {"count": len(orgs), "largest": max(len(O) for O in orgs) / len(ctx.ids)}
