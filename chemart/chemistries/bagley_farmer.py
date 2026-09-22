"""Bagley & Farmer (1992) autocatalytic metabolism. Catalog id: bagley-farmer.

Reconstructed from "Spontaneous emergence of a metabolism" (Artificial Life II):

- polymers are strings over an alphabet of `alphabet_size` letters, up to
  `max_length`; every split gives a reversible reaction A + B <-> C + H (eq. 1)
  with rate constants kf (condensation) and kr (hydrolysis, water H buffered);
- a catalyst E multiplies both directions by (1 + nu E) (eq. 9), so catalysis
  never moves the equilibrium. Strong catalytic links are drawn at random with
  probability p per (reaction, catalyst) pair (section 3.1.1) or given
  explicitly;
- with `saturation`, catalysis goes through bound complexes (appendix,
  eqs. 22-29): A + B + E -> [CE] + H and C + H + E -> [ABE] at nu kf and nu kr,
  and each complex releases its members at ku. Complexes are tracked, as in
  the paper, by one bound pool `x_bound` per species (eq. 26);
- a chemostat supplies each food species at rate delta and removes everything
  at rate K (eqs. 10-12), with K fixed by the total monomer mass m0;
- `threshold` > 0 runs the paper's deterministic metadynamics (section 3.2):
  only species at or above the threshold react; the ODE is taken to its
  fixed point, species crossing the threshold are added or removed, and the
  process repeats until the graph stops changing.
"""

from collections import Counter
from itertools import product

import numpy as np

from chemart.network import Network, Reaction, Species

LETTERS = "abcd"
MAX_SPECIES = 512
WATER = "H"
MAX_METADYNAMICS_ROUNDS = 50


def bound(s: str) -> str:
    return f"{s}_bound"


def _parse_link(text: str, species: set[str]) -> tuple[str, str, str, str]:
    """'A + B <-> C | E' -> (A, B, C, E)."""
    try:
        reaction, catalyst = (part.strip() for part in text.split("|"))
        left, c = (part.strip() for part in reaction.split("<->"))
        a, b = (part.strip() for part in left.split(" + "))
    except ValueError:
        raise ValueError(f"links: cannot parse {text!r}; expected 'A + B <-> C | E', e.g. 'a + b <-> ab | bb'") from None
    if a + b != c:
        raise ValueError(f"links: {text!r} is not a condensation, {a!r} + {b!r} gives {a + b!r}, not {c!r}")
    unknown = [s for s in (a, b, c, catalyst) if s not in species]
    if unknown:
        raise ValueError(f"links: {text!r} uses {unknown}, which are not polymers of this chemistry")
    return a, b, c, catalyst


def generate(p, rng):
    letters = LETTERS[: p.alphabet_size]
    polymers = ["".join(s) for n in range(1, p.max_length + 1) for s in product(letters, repeat=n)]
    if len(polymers) > MAX_SPECIES:
        raise ValueError(f"alphabet_size = {p.alphabet_size}, max_length = {p.max_length} gives {len(polymers)} polymers, above the limit of {MAX_SPECIES}")
    known = set(polymers)
    if not p.food_set or not all(isinstance(f, str) for f in p.food_set):
        raise ValueError("food_set must be a non-empty list of polymer strings")
    unknown = [f for f in p.food_set if f not in known]
    if unknown:
        raise ValueError(f"food_set entries {unknown} are not polymers over {letters!r} of length <= {p.max_length}")
    food = list(dict.fromkeys(p.food_set))

    pairs = [(s[:cut], s[cut:], s) for s in polymers for cut in range(1, len(s))]
    if p.links:
        wanted = list(dict.fromkeys(_parse_link(t, known) for t in p.links))
        index = {pair: i for i, pair in enumerate(pairs)}
        chosen = [(index[(a, b, c)], e) for a, b, c, e in wanted]
    else:
        mask = rng.random((len(pairs), len(polymers))) < p.p
        chosen = [(int(i), polymers[j]) for i, j in zip(*np.nonzero(mask))]
    if p.nu_distribution == "constant":
        nus = [p.nu] * len(chosen)
    else:
        nus = [float(v) for v in rng.uniform(0.0, 2.0 * p.nu, len(chosen))]

    K = p.delta * sum(len(f) for f in food) / p.m0
    initial = {f: p.m0 / sum(len(f) for f in food) for f in food}

    # Every reaction as (reactants, products, rate, species that must be active).
    full = []
    for a, b, c in pairs:
        full.append(([a, b], [c, WATER], {"law": "mass-action", "k": p.kf}, {a, b}))
        full.append(([c, WATER], [a, b], {"law": "mass-action", "k": p.kr}, {c}))
    for (i, e), nu in zip(chosen, nus):
        a, b, c = pairs[i]
        fwd = {"law": "mass-action", "k": nu * p.kf, "catalyst": e, "catalytic_efficiency": nu}
        bwd = {"law": "mass-action", "k": nu * p.kr, "catalyst": e, "catalytic_efficiency": nu}
        if p.saturation:
            full.append(([a, b, e], [bound(c), bound(e), WATER], fwd, {a, b, e}))
            full.append(([c, WATER, e], [bound(a), bound(b), bound(e)], bwd, {c, e}))
        else:
            full.append(([a, b, e], [c, WATER, e], fwd, {a, b, e}))
            full.append(([c, WATER, e], [a, b, e], bwd, {c, e}))

    analysis = None
    if p.threshold > 0:
        active, fixed_point, rounds, outcome = _metadynamics(full, food, initial, K, p)
        kept = [r for r in full if r[3] <= active]
        analysis = {
            "metadynamical_fixed_point": fixed_point,
            "active_species": sorted(active, key=lambda s: (len(s), s)),
            "metadynamics_rounds": rounds,
            "outcome": outcome,
        }
    else:
        kept = full

    reactions = [Reaction.of(lhs, rhs, dict(rate)) for lhs, rhs, rate, _ in kept]
    present = set(food) | {s for r in reactions for s in (*r.reactants, *r.products)}
    for s in [s for s in polymers if bound(s) in present]:
        reactions.append(Reaction.of([bound(s)], [s], {"law": "mass-action", "k": p.ku}))
        present.add(s)
    ids = [s for s in polymers if s in present]
    ids += [WATER] if WATER in present else []
    ids += [bound(s) for s in polymers if bound(s) in present]

    extras = {
        "food": food,
        "buffered": [WATER] if WATER in present else [],
        "flow_law": "every food species enters at constant flux delta = inflow[s]; every species leaves at first-order rate K = outflow (eqs. 10-12)",
        "catalytic_links": [
            {"reaction": f"{pairs[i][0]} + {pairs[i][1]} <-> {pairs[i][2]}", "catalyst": e, "nu": nu}
            for (i, e), nu in zip(chosen, nus)
            if analysis is None or {pairs[i][0], pairs[i][1], e} <= set(analysis["active_species"])
            or {pairs[i][2], e} <= set(analysis["active_species"])
        ],
        "conservation": [
            {
                "name": f"monomer {x}",
                "vector": {s: s.replace("_bound", "").count(x) for s in ids if s != WATER and x in s.replace("_bound", "")},
            }
            for x in letters
        ],
    }
    if analysis is not None:
        extras["analysis"] = analysis

    state = dict(initial)
    if WATER in present:
        state[WATER] = p.H
    return Network(
        species=[Species(s) for s in ids],
        reactions=reactions,
        status="truncated" if p.threshold > 0 else "complete",
        initial_state=state,
        inflow={f: p.delta for f in food},
        outflow=K,
        extras=extras,
    )


# ----------------------------------------------------------------------------
# Deterministic metadynamics (section 3.2)
# ----------------------------------------------------------------------------
def _metadynamics(full, food, initial, K, p):
    """Return (active species, fixed point, rounds, outcome).

    outcome is "fixed" when the active set stops changing. If the updates
    cycle (the active set repeats an earlier one), the sets in the cycle are
    merged, the fixed point of the merged graph is returned, and outcome is
    "cycle-merged". Hitting the round limit gives "round-limit".
    """
    active = frozenset(food)
    history = [active]
    x = dict(initial)

    def run(species, start):
        kept = [r for r in full if r[3] <= species]
        point = _fixed_point(kept, food, start, K, p)
        totals = Counter()
        for s, v in point.items():
            totals[s.replace("_bound", "")] += v
        above = frozenset(food) | {s for s, v in totals.items() if v >= p.threshold}
        return point, above

    for rounds in range(1, MAX_METADYNAMICS_ROUNDS + 1):
        x, new_active = run(active, x)
        if new_active == active:
            return set(active), _rounded(x), rounds, "fixed"
        if new_active in history:
            merged = frozenset().union(*history[history.index(new_active):])
            x, _ = run(merged, x)
            return set(merged), _rounded(x), rounds + 1, "cycle-merged"
        history.append(new_active)
        active = new_active
    return set(active), _rounded(x), MAX_METADYNAMICS_ROUNDS, "round-limit"


def _rounded(x):
    return {s: float(v) for s, v in sorted(x.items())}


def _fixed_point(kept, food, start, K, p):
    """Integrate the chemostat ODE of the reactions `kept` to its fixed point."""
    from scipy.integrate import solve_ivp
    from scipy.sparse import coo_array, identity

    reactions = [(lhs, rhs, rate["k"]) for lhs, rhs, rate, _ in kept]
    names = set(food)
    for lhs, rhs, _ in reactions:
        names.update(lhs, rhs)
    for s in [s for s in names if s.endswith("_bound")]:
        reactions.append(([s], [s[: -len("_bound")]], p.ku))
        names.add(s[: -len("_bound")])
    names.discard(WATER)
    ids = sorted(names)
    index = {s: i for i, s in enumerate(ids)}
    n, m = len(ids), len(reactions)
    pad = n                                  # index of a constant 1.0 slot

    slots = np.full((m, 3), pad)
    orders = np.zeros((m, 3))
    rates = np.empty(m)
    rows, cols, vals = [], [], []
    for j, (lhs, rhs, k) in enumerate(reactions):
        water = lhs.count(WATER)
        rates[j] = k * p.H ** water
        for c, (s, order) in enumerate(Counter(s for s in lhs if s != WATER).items()):
            slots[j, c], orders[j, c] = index[s], order
        for s, order in Counter(s for s in lhs if s != WATER).items():
            rows.append(index[s]); cols.append(j); vals.append(-order)
        for s, order in Counter(s for s in rhs if s != WATER).items():
            rows.append(index[s]); cols.append(j); vals.append(order)
    S = coo_array((vals, (rows, cols)), shape=(n, m)).tocsr()
    inflow = np.zeros(n)
    for f in food:
        inflow[index[f]] = p.delta

    def powers(x):
        xp = np.append(np.maximum(x, 0.0), 1.0)
        return xp, xp[slots] ** orders

    def f(t, x):
        _, terms = powers(x)
        return S @ (rates * terms.prod(axis=1)) + inflow - K * x

    def jac(t, x):
        xp, terms = powers(x)
        drows, dcols, dvals = [], [], []
        for c in range(3):
            others = np.prod(np.delete(terms, c, axis=1), axis=1)
            use = slots[:, c] != pad
            deriv = rates * orders[:, c] * xp[slots[:, c]] ** np.maximum(orders[:, c] - 1, 0) * others
            drows.append(np.nonzero(use)[0]); dcols.append(slots[use, c]); dvals.append(deriv[use])
        dv = coo_array((np.concatenate(dvals), (np.concatenate(drows), np.concatenate(dcols))), shape=(m, n))
        return (S @ dv.tocsc() - K * identity(n, format="csc")).toarray()

    # LSODA with a dense Jacobian (at most ~1000 variables here) is an order of
    # magnitude faster than BDF on the stiff paper-scale constants.
    x0 = np.array([start.get(s, 0.0) for s in ids])
    sol = solve_ivp(f, (0.0, 200.0 / K), x0, method="LSODA", jac=jac, rtol=1e-6, atol=1e-6 * p.threshold)
    if not sol.success:
        raise ValueError(f"metadynamics: ODE integration failed ({sol.message}); try smaller rate constants")
    return {s: max(float(v), 0.0) for s, v in zip(ids, sol.y[:, -1])}
