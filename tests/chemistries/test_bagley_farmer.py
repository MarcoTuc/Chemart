"""Bagley & Farmer (1992) "Spontaneous emergence of a metabolism": the paper's equations and results."""

import numpy as np
import pytest
from chemart.simulate import rhs
from scipy.integrate import solve_ivp

from chemart import generate_network

SMALL = dict(max_length=3, kf=2.0, kr=1.5, nu=3.0, ku=4.0, H=0.7, delta=1.3, m0=2.0)


def chemostat(net):
    """The chemostat: chemart.simulate.rhs applies the constant inflow, first-order outflow K and buffered H."""
    return rhs(net)


def closed(net):
    """The same vessel without flow."""
    net.inflow, net.outflow = None, None
    return net


def random_state(net, seed=0):
    rng = np.random.default_rng(seed)
    state = {s.id: float(v) for s, v in zip(net.species, rng.uniform(0.1, 1.0, len(net.species)))}
    state["H"] = net.initial_state["H"]
    return state


def evaluate(net, state):
    ids, g = chemostat(net)
    return dict(zip(ids, g(0.0, np.array([state[s] for s in ids]))))


def pairs_of(net):
    polymers = [s.id for s in net.species if s.id != "H" and not s.id.endswith("_bound")]
    return [(c[:cut], c[cut:], c) for c in polymers for cut in range(1, len(c))]


def links_of(net):
    out = {}
    for link in net.extras["catalytic_links"]:
        left, c = link["reaction"].split(" <-> ")
        a, b = left.split(" + ")
        out.setdefault((a, b, c), []).append((link["catalyst"], link["nu"]))
    return out


# --- eq. 9 and eqs. 10-11: unsaturated kinetics in the chemostat ------------------
def test_unsaturated_rates_are_eq_9_with_flow():
    net = generate_network("bagley-farmer", seed=3, p=0.05, saturation=False, **SMALL)
    assert net.extras["catalytic_links"], "the test needs some catalysis"
    x = random_state(net)
    paper = {s: 0.0 for s in x}
    links = links_of(net)
    for a, b, c in pairs_of(net):
        gamma = 1.0 + sum(nu * x[e] for e, nu in links.get((a, b, c), []))
        flux = gamma * (SMALL["kf"] * x[a] * x[b] - SMALL["kr"] * SMALL["H"] * x[c])   # eq. 9
        paper[c] += flux
        paper[a] -= flux
        paper[b] -= flux
    for s in paper:
        paper[s] += net.inflow.get(s, 0.0) - net.outflow * x[s]                          # eq. 10
    paper["H"] = 0.0
    ours = evaluate(net, x)
    for s in paper:
        assert ours[s] == pytest.approx(paper[s], rel=1e-10, abs=1e-12)


# --- appendix eqs. 22-29: saturation through bound pools -------------------------
def test_saturated_rates_are_appendix_eqs_24_to_29():
    net = generate_network("bagley-farmer", seed=4, p=0.05, **SMALL)
    x = random_state(net)
    kf, kr, H, ku = SMALL["kf"], SMALL["kr"], SMALL["H"], SMALL["ku"]
    paper = {s: 0.0 for s in x}
    b_ = lambda s: f"{s}_bound"   # noqa: E731
    links = links_of(net)
    for a, b, c in pairs_of(net):
        spont = kf * x[a] * x[b] - kr * H * x[c]
        paper[c] += spont
        paper[a] -= spont
        paper[b] -= spont
        for e, nu in links.get((a, b, c), []):
            bind_c = nu * kf * x[a] * x[b] * x[e]      # A + B + E -> CE + H   (eq. 22)
            bind_ab = nu * kr * x[c] * H * x[e]        # C + H + E -> ABE      (eq. 22)
            paper[a] -= bind_c
            paper[b] -= bind_c
            paper[e] -= bind_c + bind_ab
            paper[c] -= bind_ab
            paper[b_(c)] += bind_c                     # eq. 26: pools sum the complexes
            paper[b_(e)] += bind_c + bind_ab
            paper[b_(a)] += bind_ab
            paper[b_(b)] += bind_ab
    for s in list(paper):
        if s.endswith("_bound"):
            paper[s[:-6]] += ku * x[s]                 # eq. 27
            paper[s] -= ku * x[s]                      # eq. 29
    for s in paper:
        paper[s] += net.inflow.get(s, 0.0) - net.outflow * x[s]   # eqs. 10-11
    paper["H"] = 0.0
    ours = evaluate(net, x)
    assert set(ours) == set(paper)
    for s in paper:
        assert ours[s] == pytest.approx(paper[s], rel=1e-10, abs=1e-12)


# --- eq. 12: total mass relaxes to m0 ---------------------------------------------
@pytest.mark.parametrize("saturation", [True, False])
def test_total_mass_obeys_eq_12(saturation):
    net = generate_network("bagley-farmer", seed=5, p=0.05, saturation=saturation,
                           food_set=["a", "b", "ab", "ba"], **{**SMALL, "m0": 3.0})
    x = random_state(net, seed=2)
    dx = evaluate(net, x)
    length = {s.id: len(s.id.replace("_bound", "")) for s in net.species if s.id != "H"}
    m = sum(length[s] * x[s] for s in length)
    dm = sum(length[s] * dx[s] for s in length)
    food_mass = sum(len(f) for f in net.extras["food"])
    assert dm == pytest.approx(food_mass * SMALL["delta"] - net.outflow * m, rel=1e-10)
    assert sum(length[s] * v for s, v in net.initial_state.items() if s != "H") == pytest.approx(3.0)
    assert food_mass * SMALL["delta"] / net.outflow == pytest.approx(3.0)             # fixed point m0


# --- section 2.6, eqs. 13-15: catalytic focusing on the network of Figure 1 --------
def steady_state(net, t_end=60.0, x0=None):
    ids, g = chemostat(net)
    start = x0 or net.initial_state
    x = np.array([float(start.get(s, 0.0)) for s in ids])
    sol = solve_ivp(g, (0.0, t_end), x, method="LSODA", rtol=1e-10, atol=1e-12)
    assert sol.success
    return dict(zip(ids, sol.y[:, -1]))


@pytest.mark.parametrize("delta", [1e-3, 0.5, 5.0, 50.0])
def test_catalytic_focusing_eq_15(delta):
    kf, kr, H, nu, E = 1.0, 1.0, 1.0, 10.0, 1.0
    net = generate_network("bagley-farmer", max_length=2, links=["a + b <-> ab | bb"], saturation=False,
                           kf=kf, kr=kr, nu=nu, H=H, delta=delta, m0=2.0)
    net.extras["buffered"] = ["H", "bb"]                  # the enzyme is held fixed, as in Figure 1
    x = steady_state(net, t_end=60.0 / min(net.outflow, 1.0) if delta > 0.1 else 400.0,
                     x0={**net.initial_state, "bb": E})
    gamma = 1.0 + nu * E
    beta = net.outflow / (kr * H)
    assert x["ab"] / x["ba"] == pytest.approx((1 + beta) / (1 + beta / gamma), rel=1e-5)


def test_focusing_needs_driving_and_specific_catalysis():
    near_eq = generate_network("bagley-farmer", max_length=2, links=["a + b <-> ab | bb"], saturation=False,
                               kf=1.0, kr=1.0, nu=10.0, delta=1e-4)
    driven = generate_network("bagley-farmer", max_length=2, links=["a + b <-> ab | bb"], saturation=False,
                              kf=1.0, kr=1.0, nu=10.0, delta=100.0)
    ratio = {}
    for name, net in (("near", near_eq), ("driven", driven)):
        net.extras["buffered"] = ["H", "bb"]
        x = steady_state(net, t_end=2000.0 if name == "near" else 5.0, x0={**net.initial_state, "bb": 1.0})
        ratio[name] = x["ab"] / x["ba"]
    assert ratio["near"] == pytest.approx(1.0, abs=1e-3)   # beta -> 0: no focusing at equilibrium
    assert ratio["driven"] > 5.0                            # beta >> gamma: ratio -> gamma = 11


# --- section 2.3: catalysis does not change the equilibrium ------------------------
def test_detailed_balance_state_is_fixed_whatever_the_catalysis():
    for p_cat in (0.0, 0.3):
        net = generate_network("bagley-farmer", seed=6, p=p_cat, saturation=False, **SMALL)
        ids, f = rhs(closed(net))                          # closed vessel: no flow
        kappa = SMALL["kf"] / (SMALL["kr"] * SMALL["H"])
        monomer = {"a": 0.4, "b": 0.9}
        x = np.array([SMALL["H"] if s == "H" else
                      kappa ** (len(s) - 1) * np.prod([monomer[ch] for ch in s]) for s in ids])
        assert np.max(np.abs(f(0.0, x))) < 1e-12


def test_saturated_closed_vessel_relaxes_to_the_uncatalysed_equilibrium():
    net = generate_network("bagley-farmer", seed=7, p=0.2, **SMALL)
    ids, f = rhs(closed(net))
    x0 = np.array([{"a": 1.0, "b": 0.5, "H": SMALL["H"]}.get(s, 0.0) for s in ids])
    sol = solve_ivp(f, (0.0, 300.0), x0, method="LSODA", rtol=1e-10, atol=1e-13)
    x = dict(zip(ids, sol.y[:, -1]))
    kappa = SMALL["kf"] / (SMALL["kr"] * SMALL["H"])
    for a, b, c in pairs_of(net):
        assert x[c] == pytest.approx(kappa * x[a] * x[b], rel=1e-5)


# --- network structure and parameters ------------------------------------------------
def test_counts_catalysis_probability_and_conservation():
    none = generate_network("bagley-farmer", p=0.0)
    polymers = 2 + 4 + 8 + 16 + 32
    pairs = sum((n - 1) * 2 ** n for n in range(2, 6))
    assert len(none.species) == polymers + 1 and len(none.reactions) == 2 * pairs
    full = generate_network("bagley-farmer", p=1.0, max_length=3)
    n3, pairs3 = 14, 4 + 16
    assert len(full.extras["catalytic_links"]) == pairs3 * n3
    assert len(full.reactions) == 2 * pairs3 + 2 * pairs3 * n3 + n3          # + one unbinding per pool
    ids, R, P = full.matrices()
    S = (P - R).toarray()
    for law in full.extras["conservation"]:
        assert not np.any(np.array([law["vector"].get(s, 0) for s in ids]) @ S)


def test_paper_defaults_and_rate_records():
    net = generate_network("bagley-farmer", seed=1)
    assert net.status == "complete" and net.extras["buffered"] == ["H"]
    assert net.outflow == pytest.approx(17.9 * 2 / 2.0)
    catalysed = [r.rate for r in net.reactions if "catalyst" in (r.rate or {})]
    assert catalysed and all(r["k"] in (pytest.approx(8.97e5 * 649.0), pytest.approx(8.97e5 * 2.5)) for r in catalysed)


def test_links_are_validated():
    with pytest.raises(ValueError, match="not a condensation"):
        generate_network("bagley-farmer", links=["a + b <-> ba | a"])
    with pytest.raises(ValueError, match="not polymers"):
        generate_network("bagley-farmer", links=["a + b <-> ab | c"])
    with pytest.raises(ValueError, match="food_set"):
        generate_network("bagley-farmer", food_set=["x"])


# --- section 3.2: deterministic metadynamics ------------------------------------------
def test_metadynamics_reaches_a_fixed_graph_above_threshold():
    # Figure 3's constants (kf = 1e2, kr = 10, ku = 1e4 -> 1e3 here, H = 1, delta = 1e2, m0 = 3, food a, b, ab, ba)
    params = dict(max_length=4, p=0.08, kf=100.0, kr=10.0, nu=1e2, ku=1e3, delta=100.0, m0=3.0,
                  food_set=["a", "b", "ab", "ba"], threshold=1e-2)
    net = generate_network("bagley-farmer", seed=2, **params)
    analysis = net.extras["analysis"]
    assert net.status == "truncated" and analysis["outcome"] == "fixed"
    point = analysis["metadynamical_fixed_point"]
    total = {}
    for s, v in point.items():
        total[s.replace("_bound", "")] = total.get(s.replace("_bound", ""), 0.0) + v
    active = set(analysis["active_species"])
    assert set(net.extras["food"]) <= active
    assert {s for s, v in total.items() if v >= params["threshold"]} | set(net.extras["food"]) == active
    for r in net.reactions:                                 # only above-threshold species react
        assert all(s in active or s == "H" or s.endswith("_bound") for s in r.reactants)
    # the returned graph is at its dynamical fixed point
    ids, g = chemostat(net)
    x = np.array([point.get(s, net.initial_state["H"] if s == "H" else 0.0) for s in ids])
    scale = np.maximum(np.abs(x), params["threshold"])
    assert np.max(np.abs(g(0.0, x)) / scale) < 1e-3
    complete = generate_network("bagley-farmer", seed=2, **{**params, "threshold": 0.0})
    assert len(net.reactions) < len(complete.reactions)
    assert len(active) < sum(1 for s in complete.species if s.id != "H" and not s.id.endswith("_bound"))
