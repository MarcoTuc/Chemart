"""Arrhenius-gated collision algorithm (book 18.3.3). Catalog id: energy-gated-collision.

An ANALYSIS of a *reactor algorithm*, not of a chemistry. Section 18.3.3 asks
whether a stochastic reaction algorithm can honour different kinetic rates
without computing the propensities of every possible reaction, and answers with
a three-line extension of the naive bimolecular collision algorithm of chapter 2
(sec. 2.3.3):

    Pick two random molecules for collision; if their kinetic energy is higher
    than the activation energy for the reaction, then perform the reaction
    (effective collision), otherwise return the molecules to the reaction
    unchanged (elastic collision).

The probability that a collision clears the barrier is the Arrhenius factor of
eq. 2.32, k = A exp(-Ea / RT), with R = 8.31451 J/K/mol (sec. 2.2.6). So the
algorithm is a *rejection* method: a reaction of barrier Ea fires on a fraction
exp(-Ea / RT) of the collisions of its reactants, and no propensity is ever
summed over the reaction set. The book's caveat is the acceptance rate: it only
pays off when elastic collisions are rare.

`generate` builds a small reaction system carrying activation energies, runs the
gated-collision algorithm on a well-stirred multiset (`chemart.soup`), and
returns the reactions that fired with their counts (status `observed`), their
`arrhenius` rates, the per-species potential energies in extras["energies"] and
the measured-versus-predicted validation in extras["analysis"].

Two energy models, both exact and both testable:

- ``bath``       the reactor is a thermostat: each collision draws its energy
                 fresh from the thermal distribution. The energy along the line
                 of centres of a Maxwell-Boltzmann gas is exponential with mean
                 RT (Atkins & de Paula, the book's own ref. [42] for eq. 2.32),
                 so the acceptance probability is exactly exp(-Ea / RT) and the
                 measured rate constant is exactly the Arrhenius one. Energy is
                 not conserved: the bath supplies or absorbs the reaction heat.
- ``conserved``  every molecule carries a kinetic energy, initialised i.i.d.
                 exponential with mean RT, and the gate compares Ea with the
                 *sum* of the two kinetic energies, which is the book's literal
                 wording in sec. 2.2.6 ("the sum of their kinetic energies").
                 The pooled energy minus the potential-energy change is handed
                 to the products, so E = Ek + Ep (eq. 2.28) is conserved
                 exactly. The sum of two exponentials is Gamma-2, so the
                 acceptance probability is (1 + Ea/RT) exp(-Ea / RT), *not* the
                 Arrhenius factor; see the catalog `decisions`.

Energies are in kJ/mol throughout; every rate dict carries its own T and R so
the units cannot be misread.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from chemart.helpers.params import apportion
from chemart.network import Network, Reaction, Species
from chemart.soup import soup

#: Gas constant in kJ K^-1 mol^-1. The book gives R = 8.31451 J K^-1 mol^-1
#: (sec. 2.2.6, under eq. 2.32); all energies here are in kJ/mol.
R_GAS = 8.31451e-3
ENERGY_UNITS = "kJ/mol"
ENERGY_MODELS = ("bath", "conserved")

#: Barriers, energy changes and detailed balance are compared with this
#: tolerance (kJ/mol); user-supplied numbers are floats, not exact values.
TOLERANCE = 1e-9


# --------------------------------------------------------------------------
# The Arrhenius gate
# --------------------------------------------------------------------------
def arrhenius_factor(Ea: float, temperature: float) -> float:
    """exp(-Ea / RT): the Arrhenius factor of eq. 2.32, Ea in kJ/mol."""
    return float(np.exp(-Ea / (R_GAS * temperature)))


def acceptance_probability(Ea: float, temperature: float,
                           energy_model: str = "bath") -> float:
    """The probability that one collision of the reactants clears the barrier.

    ``bath``: the line-of-centres energy is Exponential(RT), so the acceptance
    is exactly the Arrhenius factor exp(-Ea/RT). ``conserved``: the gate uses
    the sum of two Exponential(RT) kinetic energies, which is Gamma(2, RT), so
    the acceptance is (1 + Ea/RT) exp(-Ea/RT).
    """
    if energy_model not in ENERGY_MODELS:
        raise ValueError(f"energy_model must be one of {ENERGY_MODELS}, got {energy_model!r}")
    x = Ea / (R_GAS * temperature)
    return float(np.exp(-x)) if energy_model == "bath" else float((1.0 + x) * np.exp(-x))


# --------------------------------------------------------------------------
# The reaction system
# --------------------------------------------------------------------------
def default_system(Ea: float, delta_G: float) -> dict:
    """One reversible bimolecular reaction, X1 + X2 <-> Y1 + Y2.

    The product well sits at delta_G, so the reverse barrier is Ea - delta_G
    (book fig. 2.4 and the text under it: Ea(YX) = Ea(XY) - dG). Forward and
    reverse therefore satisfy detailed balance, kf/kr = exp(-dG/RT), eq. 2.33.
    """
    return {
        "energies": {"X1": 0.0, "X2": 0.0, "Y1": float(delta_G), "Y2": 0.0},
        "reactions": [
            {"id": "X1+X2->Y1+Y2", "reactants": ["X1", "X2"],
             "products": ["Y1", "Y2"], "Ea": float(Ea)},
            {"id": "Y1+Y2->X1+X2", "reactants": ["Y1", "Y2"],
             "products": ["X1", "X2"], "Ea": float(Ea) - float(delta_G)},
        ],
        "initial": {"X1": 1.0, "X2": 1.0},
    }


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _species_list(name: str, value) -> list[str]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a non-empty list of species ids, got {value!r}")
    bad = [x for x in value if not isinstance(x, str) or not x or " " in x]
    if bad:
        raise ValueError(f"{name} has invalid species ids {bad}")
    return list(value)


def _multiset(names) -> tuple:
    return tuple(sorted(Counter(names).items()))


def _pair(names) -> tuple:
    return tuple(sorted(names))


def validated_system(raw: dict) -> dict:
    """Check a reaction system and annotate every reaction with its delta_G.

    The algorithm is bimolecular (the book picks *two* molecules), every
    species needs a potential energy, and a barrier may not lie below the well
    it leads to: Ea >= max(0, delta_G), or the transition state of fig. 2.4
    would be under the products.
    """
    if not isinstance(raw, dict) or set(raw) - {"energies", "reactions", "initial"}:
        raise ValueError(
            "system must be an object with keys 'energies', 'reactions' and optionally "
            f"'initial', got {sorted(raw) if isinstance(raw, dict) else raw!r}"
        )
    energies = raw.get("energies")
    if not isinstance(energies, dict) or not energies:
        raise ValueError(f"system['energies'] must map every species id to a potential energy in {ENERGY_UNITS}")
    for s, e in energies.items():
        if not isinstance(s, str) or not s or " " in s or not _is_number(e):
            raise ValueError(f"system['energies'] entry {s!r}: {e!r} is not a species id and a number")
    potential = {s: float(e) for s, e in energies.items()}

    items = raw.get("reactions")
    if not isinstance(items, list) or not items:
        raise ValueError("system['reactions'] must be a non-empty list of reactions")

    reactions, seen_pairs = [], {}
    for i, item in enumerate(items):
        where = f"system['reactions'][{i}]"
        if not isinstance(item, dict) or set(item) - {"id", "reactants", "products", "Ea"}:
            raise ValueError(f"{where} must be an object with keys id, reactants, products, Ea; got {item!r}")
        reactants = _species_list(f"{where}['reactants']", item.get("reactants"))
        products = _species_list(f"{where}['products']", item.get("products"))
        if len(reactants) != 2:
            raise ValueError(
                f"{where} has {len(reactants)} reactants: the collision algorithm draws exactly "
                "two molecules, so every reaction must be bimolecular"
            )
        unknown = [x for x in reactants + products if x not in potential]
        if unknown:
            raise ValueError(f"{where} uses species {unknown} that have no entry in system['energies']")
        if _multiset(reactants) == _multiset(products):
            raise ValueError(f"{where} does not change the multiset, so it can never be an effective collision")
        Ea = item.get("Ea")
        if not _is_number(Ea) or Ea < 0:
            raise ValueError(f"{where}['Ea'] must be a non-negative activation energy in {ENERGY_UNITS}, got {Ea!r}")
        delta_G = sum(potential[x] for x in products) - sum(potential[x] for x in reactants)
        if Ea < delta_G - TOLERANCE:
            raise ValueError(
                f"{where}: Ea = {Ea} is below its energy change delta_G = {delta_G}; the "
                "activation barrier cannot lie under the product well (book fig. 2.4). Raise "
                "Ea or lower delta_G"
            )
        rid = str(item.get("id", f"r{i + 1}"))
        reactions.append({"id": rid, "reactants": reactants, "products": products,
                          "Ea": float(Ea), "delta_G": float(delta_G)})
        key = (_multiset(reactants), _multiset(products))
        if key in seen_pairs:
            raise ValueError(f"{where} repeats the reaction of {seen_pairs[key]!r}")
        seen_pairs[key] = rid
    if len({r["id"] for r in reactions}) != len(reactions):
        raise ValueError("system reaction ids must be unique")

    # Detailed balance: where both directions are present, fig. 2.4 fixes the
    # reverse barrier at Ea(YX) = Ea(XY) - dG.
    by_key = {(_multiset(r["reactants"]), _multiset(r["products"])): r for r in reactions}
    for r in reactions:
        back = by_key.get((_multiset(r["products"]), _multiset(r["reactants"])))
        if back is not None and abs(back["Ea"] - (r["Ea"] - r["delta_G"])) > 1e-6:
            raise ValueError(
                f"reactions {r['id']!r} and {back['id']!r} are each other's reverse but are not "
                f"thermodynamically consistent: Ea({back['id']}) must be Ea({r['id']}) - delta_G = "
                f"{r['Ea'] - r['delta_G']}, got {back['Ea']} (book fig. 2.4, eq. 2.33)"
            )

    initial = raw.get("initial")
    if initial is None:
        initial = {s: 1.0 for s in reactions[0]["reactants"]}
    if not isinstance(initial, dict) or not initial:
        raise ValueError("system['initial'] must map species ids to positive starting proportions")
    for s, w in initial.items():
        if s not in potential:
            raise ValueError(f"system['initial'] uses species {s!r} that has no entry in system['energies']")
        if not _is_number(w) or w <= 0:
            raise ValueError(f"system['initial'][{s!r}] must be a positive proportion, got {w!r}")

    return {"energies": potential, "reactions": reactions,
            "initial": {s: float(w) for s, w in initial.items()}}


# --------------------------------------------------------------------------
# The algorithm
# --------------------------------------------------------------------------
def _two_indices(rng, n: int) -> tuple[int, int]:
    """Two distinct indices in range(n), as `chemart.soup` draws molecules."""
    while True:
        i, j = (int(v) for v in rng.integers(n, size=2))
        if i != j:
            return i, j


def _split(rng, total: float, n: int) -> list[float]:
    """Share `total` among n products as uniform spacings, Dirichlet(1, ..., 1).

    For n = 2 this is the Beta(1,1) split that maps Gamma(2, RT) back onto two
    independent Exponential(RT) energies, so an equilibrated pool of kinetic
    energies stays Maxwell-Boltzmann under thermoneutral reactions. The last
    share is taken by subtraction, so the parts sum to `total`.
    """
    if n == 1:
        return [total]
    cuts = np.sort(rng.random(n - 1))
    parts, previous = [], 0.0
    for c in cuts:
        parts.append(total * (float(c) - previous))
        previous = float(c)
    parts.append(total - sum(parts))
    return parts


def _drop(pool: list, i: int) -> None:
    pool[i] = pool[-1]
    pool.pop()


def run(system: dict, rng, temperature: float = 300.0, population: int = 400,
        steps: int = 4000, energy_model: str = "bath") -> dict:
    """Run the gated-collision algorithm on a well-stirred multiset.

    One step is one collision: `chemart.soup` draws two distinct molecules, the
    rule for that unordered pair is looked up (no rule -> elastic), its barrier
    is tested against the collision energy (too low -> elastic), and only then
    are the reactants replaced by the products. Elastic collisions leave the
    vessel untouched, exactly as in sec. 2.3.3 and 18.3.3.

    Returns the fired reactions, the final population, the per-reaction attempt
    and acceptance tallies and the energy budget.
    """
    if energy_model not in ENERGY_MODELS:
        raise ValueError(f"energy_model must be one of {ENERGY_MODELS}, got {energy_model!r}")
    if temperature <= 0:
        raise ValueError(f"temperature must be positive, got {temperature}")

    potential = system["energies"]
    channels: dict[tuple, list[dict]] = {}
    for r in system["reactions"]:
        channels.setdefault(_pair(r["reactants"]), []).append(r)

    names = list(system["initial"])
    counts = apportion(population, [system["initial"][s] for s in names])
    start = [s for s, n in zip(names, counts) for _ in range(n)]
    if len(start) < 2:
        raise ValueError(f"population = {population} leaves fewer than two molecules in the vessel")

    kinetic = [float(v) for v in rng.exponential(R_GAS * temperature, size=len(start))] \
        if energy_model == "conserved" else []
    kinetic_initial = float(sum(kinetic))
    stats = {r["id"]: {"attempts": 0, "effective": 0} for r in system["reactions"]}

    def react(m1, m2):
        options = channels.get(_pair((m1, m2)))
        if not options:
            return None                                    # no rule: elastic
        r = options[0] if len(options) == 1 else options[int(rng.integers(len(options)))]
        tally = stats[r["id"]]
        tally["attempts"] += 1
        if energy_model == "bath":
            energy = float(rng.exponential(R_GAS * temperature))
            if energy < r["Ea"]:
                return None                                # barrier not cleared: elastic
        else:
            i, j = _two_indices(rng, len(kinetic))
            energy = kinetic[i] + kinetic[j]
            if energy < r["Ea"]:
                return None
            shares = _split(rng, energy - r["delta_G"], len(r["products"]))
            for k in sorted((i, j), reverse=True):
                _drop(kinetic, k)
            kinetic.extend(shares)
        tally["effective"] += 1
        return list(r["products"])

    fired, final = soup(react, start, steps, rng, arity=2)

    potential_initial = float(sum(potential[m] for m in start))
    potential_final = float(sum(potential[m] for m in final))
    budget = {"units": ENERGY_UNITS,
              "potential_initial": potential_initial,
              "potential_final": potential_final}
    if energy_model == "conserved":
        kinetic_final = float(sum(kinetic))
        budget.update({
            "kinetic_initial": kinetic_initial, "kinetic_final": kinetic_final,
            "total_initial": kinetic_initial + potential_initial,
            "total_final": kinetic_final + potential_final,
            "drift": (kinetic_final + potential_final) - (kinetic_initial + potential_initial),
            "law": "E = Ek + Ep is conserved exactly (eq. 2.28): the pooled kinetic energy of "
                   "the colliding pair, minus the potential energy change, goes to the products",
        })
    else:
        budget.update({
            "kinetic_tracked": False,
            "supplied_by_bath": potential_final - potential_initial,
            "law": "a thermostat holds T fixed and supplies (or absorbs) the reaction heat, so "
                   "Ek is not tracked and only Ep is booked",
        })

    return {"fired": fired, "final": final, "stats": stats, "budget": budget,
            "initial_counts": {s: int(n) for s, n in zip(names, counts) if n > 0}}


# --------------------------------------------------------------------------
def generate(p, rng):
    system = validated_system(p.system if p.system else default_system(p.Ea, p.delta_G))
    result = run(system, rng, temperature=p.temperature, population=p.population,
                 steps=p.steps, energy_model=p.energy_model)

    by_sides = {(_multiset(r["reactants"]), _multiset(r["products"])): r
                for r in system["reactions"]}
    reactions = []
    for lhs, rhs, count in result["fired"]:
        r = by_sides[(_multiset(lhs), _multiset(rhs))]
        rate = {"law": "arrhenius", "A": float(p.A), "Ea": r["Ea"],
                "T": float(p.temperature), "R": R_GAS, "units": ENERGY_UNITS}
        reactions.append(Reaction.of(lhs, rhs, rate=rate, count=count))

    effective = sum(t["effective"] for t in result["stats"].values())
    measured = {}
    for r in system["reactions"]:
        tally = result["stats"][r["id"]]
        acceptance = tally["effective"] / tally["attempts"] if tally["attempts"] else None
        factor = arrhenius_factor(r["Ea"], p.temperature)
        measured[r["id"]] = {
            "Ea": r["Ea"], "delta_G": r["delta_G"],
            "attempts": tally["attempts"], "effective": tally["effective"],
            "acceptance": acceptance,
            "predicted_acceptance": acceptance_probability(r["Ea"], p.temperature, p.energy_model),
            "arrhenius_factor": factor,
            "k_observed": float(p.A) * acceptance if acceptance is not None else None,
            "k_arrhenius": float(p.A) * factor,
        }

    analysis = {
        "energy_model": p.energy_model,
        "temperature": float(p.temperature),
        "gas_constant": R_GAS,
        "pre_exponential": float(p.A),
        "collisions": int(p.steps),
        "generations": float(p.steps) / len(result["final"]) if result["final"] else 0.0,
        "effective_collisions": effective,
        "elastic_collisions": int(p.steps) - effective,
        "elastic_fraction": (int(p.steps) - effective) / p.steps if p.steps else None,
        "reactions": measured,
        "final_population": {s: int(n) for s, n in sorted(Counter(result["final"]).items())},
        "energy_budget": result["budget"],
        "validation": "k_observed = A x (effective collisions / collisions of the reactants) "
                      "against k_arrhenius = A exp(-Ea/RT), eq. 2.32",
    }

    extras = {
        "energies": {
            "units": ENERGY_UNITS,
            "gas_constant": R_GAS,
            "temperature": float(p.temperature),
            "potential": dict(system["energies"]),
            "activation": {r["id"]: r["Ea"] for r in system["reactions"]},
            "delta_G": {r["id"]: r["delta_G"] for r in system["reactions"]},
            "law": "E = Ek + Ep (eq. 2.28); a collision is effective iff its kinetic energy "
                   "reaches Ea, which happens with probability exp(-Ea/RT) (eq. 2.32)",
        },
        "analysis": analysis,
    }
    return Network(
        species=[Species(s) for s in system["energies"]],
        reactions=reactions,
        status="observed",
        initial_state={s: float(n) for s, n in result["initial_counts"].items()},
        extras=extras,
    )
