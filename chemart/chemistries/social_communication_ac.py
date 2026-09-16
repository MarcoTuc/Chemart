"""Artificial chemistry of social communication (book 20.4; Dittrich, Kron & Banzhaf 2003).

Catalog id: social-communication-ac.

Agents ("Ego" and "Alter") show signs carrying one of N messages. An *activity*
is an agent changing the number on its own sign after observing another agent's
sign, so the population of displayed signs is a well-stirred multiset of
n_agents molecules and every communication event is the reaction

    a_i + a_k -> a_i + a_j

where a_i is the sign Ego reacted to (Alter's, a catalyst), a_k the sign Ego
replaced and a_j the activity Ego selected. The returned network is the set of
events that fired, with counts (status observed).

Activity selection (paper Sec. 2.5), for agent Ego reacting to message a:

    w_EE[i] = lookup(M_ego, a, i)                     expectation-expectation
    w_EC[i] = f_certainty(lookup(M_alter, i))         expectation-certainty
    f[i]    = (1 - alpha) w_EE[i] + alpha w_EC[i] + c_f / N          (eq. 2)
    w_AV    = normalize(f);  w_AP = normalize(w_AV ** gamma)         (eq. 3)

with f_certainty(p_1..p_N) = 1 + sum_i p_i log_N p_i (eq. 1), zero for the
uniform vector and one for a point mass. Activity i is drawn with probability
w_AP[i]. ee_memory="alter" replaces M_ego by M_alter in w_EE (Sec. 6.6): Ego
then expects of itself what *other* agents did when reacting to it, which is
the change that makes order scalable.

The memory is the paper's simple neuronal matrix memory (Sec. 2.12, memory
model 05): rows of an N x N matrix initialised to 1/N, and

    memorize(M, a, b):  m[a, b] += r_learn;  m += r_forget / N;
                        normalise every row.                     (eqs. 8-10)
"""

from collections import Counter

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from chemart.network import Network, Reaction, Species

#: paper's `intervalSize`: the window over which different activities are counted.
WINDOW = 50
#: paper's `startAverage`: measurement starts here, after the transient phase.
START_AVERAGE = 500

INTERACTIONS = ("dyadic", "population")
EE_MEMORIES = ("ego", "alter")


# --- the model ---------------------------------------------------------------
def certainty(p) -> float:
    """f_certainty(p_1, ..., p_N) = 1 + sum_i p_i log_N p_i (eq. 1).

    0.0 for the uniform vector (no information), 1.0 for a point mass.
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    if n < 2:
        return 1.0
    nz = p[p > 0]
    return float(1.0 + (nz * (np.log(nz) / np.log(n))).sum())


def new_memory(N: int) -> np.ndarray:
    """Memory matrix initialised with m[a, b] = 1/N (Sec. 2.12)."""
    return np.full((N, N), 1.0 / N)


def memorize(m: np.ndarray, a: int, b: int, r_learn: float, r_forget: float) -> None:
    """Store the event (a, b) in place: learn, forget globally, normalise rows."""
    m[a, b] += r_learn                       # eq. 8
    m += r_forget / m.shape[0]               # eq. 9
    m /= m.sum(axis=1, keepdims=True)        # eq. 10


class Agent:
    """One social agent: two memories and the message currently on its sign."""

    def __init__(self, N: int, sign: int):
        self.ego = new_memory(N)             # how I reacted to messages
        self.alter = new_memory(N)           # how others reacted to my messages
        self.sign = sign

    def behaviour(self, p) -> np.ndarray:
        """Activity probabilities w_AP as a matrix [received message, activity]."""
        ec = np.array([certainty(row) for row in self.alter])          # eq. 1
        ee = self.alter if p.ee_memory == "alter" else self.ego
        f = (1.0 - p.alpha) * ee + p.alpha * ec[None, :] + p.c_f / ee.shape[0]   # eq. 2
        total = f.sum(axis=1, keepdims=True)
        w_av = np.divide(f, total, out=np.full_like(f, 1.0 / f.shape[1]), where=total > 0)
        w = w_av**p.gamma                                              # eq. 3
        total = w.sum(axis=1, keepdims=True)
        w_ap = np.divide(w, total, out=np.full_like(w, 1.0 / w.shape[1]), where=total > 0)
        return w_av, w_ap


def _draw(rng, weights: np.ndarray) -> int:
    return int(min(np.searchsorted(np.cumsum(weights), rng.random()), weights.size - 1))


def _check(p) -> None:
    if p.interaction == "dyadic" and p.n_agents != 2:
        raise ValueError(
            "interaction='dyadic' is the paper's two-agent world: it needs "
            f"n_agents = 2, got n_agents={p.n_agents} (use interaction='population')"
        )
    if p.n_observers > p.n_agents - 2:
        raise ValueError(
            f"n_observers must leave room for Ego and Alter: n_observers={p.n_observers} "
            f"> n_agents - 2 = {p.n_agents - 2}"
        )


# --- the run -----------------------------------------------------------------
def generate(p, rng):
    _check(p)
    N, M = p.N, p.n_agents
    agents = [Agent(N, int(rng.integers(N))) for _ in range(M)]
    initial = Counter(a.sign for a in agents)

    start = min(START_AVERAGE, p.steps // 2)
    events: dict[tuple[int, int, int], list] = {}
    performed: list[int] = []
    certainties: list[float] = []

    for step in range(p.steps):
        if p.interaction == "dyadic":
            ego_i, alter_i = step % 2, (step + 1) % 2
        else:
            ego_i = int(rng.integers(M))
            alter_i = int(rng.integers(M - 1))
            alter_i += alter_i >= ego_i
        ego, alter = agents[ego_i], agents[alter_i]

        received = alter.sign
        replaced = ego.sign
        w_av, w_ap = ego.behaviour(p)
        selected = _draw(rng, w_ap[received])
        ego.sign = selected

        memorize(ego.ego, received, selected, p.r_learn, p.r_forget)
        memorize(alter.alter, received, selected, p.r_learn, p.r_forget)
        if p.n_observers:                    # Sec. 6.13: observers learn by watching
            pool = [i for i in range(M) if i not in (ego_i, alter_i)]
            for i in rng.choice(len(pool), size=p.n_observers, replace=False):
                memorize(agents[pool[int(i)]].alter, received, selected, p.r_learn, p.r_forget)

        key = (received, replaced, selected)
        entry = events.setdefault(key, [key, 0])
        entry[1] += 1
        if step >= start:
            performed.append(selected)
            certainties.append(certainty(w_av[received]))

    return _network(p, agents, initial, events, performed, certainties, start)


# --- the measures (Sec. 4) ---------------------------------------------------
def _behaviour_matrix(p, agents) -> np.ndarray:
    """p[i][j]: probability that a randomly drawn agent reacts with j to message i."""
    return np.mean([a.behaviour(p)[1] for a in agents], axis=0)


def _systems_level_order(agents, behaviour: np.ndarray) -> float:
    """O_P = (1/M) sum_i M_i f_certainty(p_i1, ..., p_iN) (eq. 15)."""
    shown = Counter(a.sign for a in agents)
    return float(sum(n * certainty(behaviour[i]) for i, n in shown.items()) / len(agents))


def _activity_systems(behaviour: np.ndarray, threshold: float) -> list[list[int]]:
    """Sets O that are closed and self-maintaining in the thresholded graph (Sec. 7.6).

    An edge (v1, v2) survives if w(v1, v2) > threshold. A set is closed iff no
    edge leaves it and self-maintaining iff every member has an incoming edge
    from inside, so the minimal ones are exactly the bottom (terminal) strongly
    connected components that feed themselves.
    """
    adj = behaviour > threshold
    n_comp, label = connected_components(csr_array(adj), directed=True, connection="strong")
    out = []
    for c in range(n_comp):
        members = np.flatnonzero(label == c)
        inside = np.isin(label, [c])
        if adj[np.ix_(members, ~inside)].any():          # an edge leaves: not closed
            continue
        if not adj[np.ix_(members, members)].any(axis=0).all():   # not self-maintaining
            continue
        out.append([int(v) + 1 for v in members])
    return sorted(out)


def _different_activities(performed: list[int]) -> float:
    """Average number of different activities per interval of WINDOW steps."""
    windows = [performed[i:i + WINDOW] for i in range(0, len(performed) - WINDOW + 1, WINDOW)]
    if not windows:
        return float(len(set(performed)))
    return float(np.mean([len(set(w)) for w in windows]))


def _network(p, agents, initial, events, performed, certainties, start) -> Network:
    def sid(i: int) -> str:
        return f"a{i + 1}"

    reactions = [
        Reaction.of([sid(received), sid(replaced)], [sid(received), sid(selected)], count=count)
        for (received, replaced, selected), count in events.values()
    ]
    behaviour = _behaviour_matrix(p, agents)
    extras = {
        "agents": [
            {
                "sign": int(a.sign) + 1,
                "ego_memory": np.round(a.ego, 6).tolist(),
                "alter_memory": np.round(a.alter, 6).tolist(),
            }
            for a in agents
        ],
        "interaction_law": {
            "unit": "an activity: one agent changing the message on its sign after "
                    "observing another agent's sign. Communications, not persons, are "
                    "the species; agents are the (stateful) context of a reaction.",
            "activity_value": "w_AV = normalize((1 - alpha) w_EE + alpha w_EC + c_f / N)  "
                              "(eqs. 1-2)",
            "activity_probability": "w_AP = normalize(w_AV ** gamma)  (eq. 3)",
            "expectation_expectation": f"w_EE[i] = lookup(M_{p.ee_memory}, received, i)",
            "expectation_certainty": "w_EC[i] = f_certainty(lookup(M_alter, i)), "
                                     "f_certainty(p) = 1 + sum_i p_i log_N p_i",
            "memory": "m[a, b] += r_learn; m += r_forget / N; rows normalised (eqs. 8-10)",
        },
        "analysis": {
            "steps": int(p.steps),
            "measurement_start": int(start),
            "window": WINDOW,
            "O_P": _systems_level_order(agents, behaviour),
            "O_AV": float(np.mean(certainties)) if certainties else 0.0,
            "different_activities": _different_activities(performed),
            "activities_used": sorted({int(i) + 1 for i in performed}),
            "behaviour_matrix": np.round(behaviour, 6).tolist(),
            "activity_graph": {
                "threshold": float(p.edge_threshold),
                "edges": [
                    [int(i) + 1, int(j) + 1, float(round(behaviour[i, j], 6))]
                    for i in range(p.N) for j in range(p.N)
                    if behaviour[i, j] > p.edge_threshold
                ],
            },
            "activity_systems": _activity_systems(behaviour, p.edge_threshold),
        },
    }
    return Network(
        species=[Species(sid(i)) for i in range(p.N)],
        reactions=reactions,
        status="observed",
        initial_state={sid(i): int(n) for i, n in sorted(initial.items())},
        extras=extras,
    )
