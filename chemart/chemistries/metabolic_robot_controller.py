"""Artificial reaction graphs for metabolic robot control (Ziegler & Banzhaf 2001). Catalog id: metabolic-robot-controller.

One individual of the paper's genetic-programming system: a random weighted
bipartite reaction graph (sec. 3) that fulfils material balance (sec. 3.1.1).

- Reaction nodes have the types of eqs. 4-7 (s1 -> s1', s1 + s2 -> s1',
  s1 -> s1' + s2, s1 + s2 -> s1' + s2'), with participants drawn from all
  molecules, as in mutation case 4 (sec. 3.2.1). A node is kept only while a
  strictly positive species-weight vector m with m^T M = 0 still exists.
- Stoichiometric edges carry the default weight k_0 (spontaneous reaction). A
  node may carry one catalytic or inhibitory edge with k ~ U[k_min, k_max].
- Catalysis follows the mass-action form of Ziegler, Dittrich & Banzhaf (1998,
  appendix, example A): a parallel channel with the catalyst on both sides and
  constant k_0 * kappa, so the total rate is k_0 prod[X] (1 + kappa [C]).
  No source gives a functional form for inhibition, so the inhibitor and its
  constant are recorded as extra keys on the reaction's rate.
- Sensor substances receive the inflow of eq. 30; actuator substances are
  consumed at rate alpha above the threshold a_min (eq. 31).
"""

import math
from collections import Counter
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from chemart.helpers.explicit import network, term

#: reaction type (paper eqs. 4-7) -> (number of reactants, number of products)
TYPES = {4: (1, 1), 5: (2, 1), 6: (1, 2), 7: (2, 2)}

MAX_SENSOR = 1023           # maximum Khepera sensor value (eq. 30)
ACTIONS = {
    "rotate": "above the threshold, switch the rotation direction of one motor (Ziegler & Banzhaf 2001, sec. 5.2)",
    "left-motor": "concentration drives the left wheel motor (Ziegler, Dittrich & Banzhaf 1998, fig. 10)",
    "right-motor": "concentration drives the right wheel motor (Ziegler, Dittrich & Banzhaf 1998, fig. 10)",
}


def material_balance(M: np.ndarray) -> np.ndarray | None:
    """Positive species weights m >= 1 with M m = 0 (M: reactions x species), or None.

    This is the test of sec. 3.1.1: the network is consistent iff the
    system M(R) m = 0 has a solution with every weight positive.
    """
    M = np.atleast_2d(np.asarray(M, dtype=float))
    n = M.shape[1]
    result = linprog(np.zeros(n), A_eq=M, b_eq=np.zeros(M.shape[0]), bounds=[(1, None)] * n, method="highs")
    return result.x if result.status == 0 else None


def _integer_weights(m: np.ndarray, M: np.ndarray) -> list:
    scaled = [Fraction(v / m.min()).limit_denominator(1000) for v in m]
    lcm = math.lcm(*(f.denominator for f in scaled))
    ints = [int(f * lcm) for f in scaled]
    return ints if not M.size or not np.any(M @ np.array(ints)) else [float(v) for v in m]


def _check_name(what: str, name) -> str:
    if not isinstance(name, str) or not name or any(ch.isspace() for ch in name) or "->" in name:
        raise ValueError(f"{what}: substance names must be non-empty strings without spaces or '->', got {name!r}")
    return name


def _io(p):
    if not isinstance(p.sensor_map, dict) or not p.sensor_map:
        raise ValueError(f"sensor_map must map each input substance to a list of sensor names, got {p.sensor_map!r}")
    sensors = {}
    for substance, names in p.sensor_map.items():
        _check_name("sensor_map", substance)
        if not isinstance(names, list) or not names or not all(isinstance(s, str) and s for s in names):
            raise ValueError(f"sensor_map[{substance!r}] must be a non-empty list of sensor names, got {names!r}")
        sensors[substance] = names
    known = {s for names in sensors.values() for s in names}
    for sensor, value in p.sensor_readings.items():
        if sensor not in known:
            raise ValueError(f"sensor_readings: {sensor!r} is not a sensor named in sensor_map ({sorted(known)})")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= MAX_SENSOR:
            raise ValueError(f"sensor_readings[{sensor!r}] must be a number in [0, {MAX_SENSOR}], got {value!r}")
    if not isinstance(p.actuator_map, dict) or not p.actuator_map:
        raise ValueError(f"actuator_map must map each output substance to an action, got {p.actuator_map!r}")
    for substance, action in p.actuator_map.items():
        _check_name("actuator_map", substance)
        if action not in ACTIONS:
            raise ValueError(f"actuator_map[{substance!r}] must be one of {sorted(ACTIONS)}, got {action!r}")
    both = sorted(set(sensors) & set(p.actuator_map))
    if both:
        raise ValueError(f"input and output sets must be disjoint (I ∩ O = ∅); {both} appear in both sensor_map and actuator_map")
    return sensors


def generate(p, rng):
    sensors = _io(p)
    if p.k_min > p.k_max:
        raise ValueError(f"k_min must be <= k_max, got k_min={p.k_min}, k_max={p.k_max}")
    types = p.reaction_types
    if not isinstance(types, list) or not types or any(isinstance(t, bool) or t not in TYPES for t in types):
        raise ValueError(f"reaction_types must be a non-empty list drawn from {sorted(TYPES)} (paper eqs. 4-7), got {types!r}")
    types = sorted(set(types))

    io_names = [*p.actuator_map, *sensors]
    if p.n_substances < len(io_names):
        raise ValueError(f"n_substances={p.n_substances} must be at least the {len(io_names)} input and output substances {io_names}")
    names, i = list(io_names), 1
    while len(names) < p.n_substances:
        if f"s{i}" not in names:
            names.append(f"s{i}")
        i += 1
    m = len(names)

    rows, nodes, seen, mass, unique = [], [], set(), np.ones(m), False
    budget = 10000 + 100 * p.n_reactions
    for _ in range(budget):
        if len(nodes) >= p.n_reactions:
            break
        t = types[int(rng.integers(len(types)))]
        n_in, n_out = TYPES[t]
        lhs = tuple(sorted(int(x) for x in rng.integers(m, size=n_in)))
        rhs = tuple(sorted(int(x) for x in rng.integers(m, size=n_out)))
        if lhs == rhs or (lhs, rhs) in seen:
            continue                              # no-op or an existing reaction node
        row = np.zeros(m)
        for s in lhs:
            row[s] -= 1
        for s in rhs:
            row[s] += 1
        if unique:
            if abs(row @ mass) > 1e-9:
                continue
        elif abs(row @ mass) < 1e-9:
            # the current positive weights already balance it: keep them, no LP needed
            unique = m - np.linalg.matrix_rank(np.vstack(rows + [row])) == 1
        else:
            trial = np.vstack(rows + [row])
            w = material_balance(trial)
            if w is None:
                continue
            mass = w
            unique = m - np.linalg.matrix_rank(trial) == 1
        rows.append(row)
        seen.add((lhs, rhs))
        nodes.append((t, lhs, rhs))
    else:
        if len(nodes) < p.n_reactions:
            raise ValueError(
                f"only {len(nodes)} of n_reactions={p.n_reactions} materially balanced reactions were found "
                f"in {budget} draws; lower n_reactions or raise n_substances"
            )

    def side(indices):
        counts = Counter(names[s] for s in indices)
        return " + ".join(term(n, s) for s, n in counts.items())

    reactions, graph = [], []
    for j, (t, lhs, rhs) in enumerate(nodes):
        node_id = f"r{j + 1}"
        base = {"law": "mass-action", "k": float(p.k_0), "node": node_id}
        entry = {"id": node_id, "type": t, "reactants": [names[s] for s in lhs],
                 "products": [names[s] for s in rhs], "k": float(p.k_0), "network_reactions": [len(reactions)]}
        modifier = None
        if rng.random() < p.p_modifier:
            effect = "inhibitive" if rng.random() < p.p_inhibitor else "catalytic"
            species = names[int(rng.integers(m))]
            k = float(rng.uniform(p.k_min, p.k_max))
            modifier = {"species": species, "effect": effect, "k": k}
            entry["modifier"] = modifier
        text = f"{side(lhs)} -> {side(rhs)}"
        if modifier and modifier["effect"] == "inhibitive":
            base.update({"inhibitor": modifier["species"], "k_inhibition": modifier["k"]})
        reactions.append((text, base))
        if modifier and modifier["effect"] == "catalytic":
            c = modifier["species"]
            entry["network_reactions"].append(len(reactions))
            reactions.append((
                f"{side(lhs + (names.index(c),))} -> {side(rhs + (names.index(c),))}",
                {"law": "mass-action", "k": float(p.k_0 * modifier["k"]), "node": node_id,
                 "catalyst": c, "k_catalysis": modifier["k"]},
            ))
        graph.append(entry)

    inflow = {
        s: float(p.max_inflow * max(p.sensor_readings.get(x, 0) for x in sensors[s]) / MAX_SENSOR)
        for s in sensors
    }
    M = np.array(rows) if rows else np.zeros((0, m))
    weights = _integer_weights(mass, M) if rows else [1] * m
    return network(
        reactions,
        species=names,
        inflow=inflow,
        outflow={s: float(p.alpha) for s in p.actuator_map},
        extras={
            "conservation": [{"name": "material balance", "vector": dict(zip(names, weights)), "unique": bool(unique)}],
            "reaction_graph": graph,
            "input_set": list(sensors),
            "output_set": list(p.actuator_map),
            "sensors": {
                s: {"sensors": sensors[s],
                    "inflow": "max_inflow * max(sensor values) / 1023 per reactor iteration, as a fraction of the reactor volume (eq. 30)"}
                for s in sensors
            },
            "actuators": {
                s: {"action": a, "meaning": ACTIONS[a], "threshold": float(p.a_min), "consumption_rate": float(p.alpha),
                    "outflow": "alpha * [s] while [s] > a_min, otherwise 0 (eq. 31)"}
                for s, a in p.actuator_map.items()
            },
        },
    )
