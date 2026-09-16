"""Isologous diversification: the equations, constants and behaviour of Kaneko & Yomo.

Sources: K. Kaneko & T. Yomo, "Isologous diversification: a theory of cell
differentiation", Bull. Math. Biol. 59(1):139-196 (1997), preprint
arXiv:adap-org/9606002 (eqs. 1-9, section 5, appendix 1, captions of Figs. 4,
6, 15 and 18), and "Emergence of rules in cell society", Bull. Math. Biol.
60:659-687 (1998), preprint arXiv:adap-org/9802002 (eqs. 1-6, section 2.4).
"""

import numpy as np
import pytest
from odes import rhs

from chemart import catalog, generate_network
from chemart.chemistries import isologous_diversification as I

ID = "isologous-diversification"
ENTRY = {c.id: c for c in catalog.load()}[ID]
DEFAULTS = {p.name: p.default for p in ENTRY.params}

# A small network pinned by hand, so that the tests do not depend on a draw:
# 1 -> 2 -> 3 -> 1 and 1 -> 3, each path catalysed by its own product.
RING = [[1, 2, 2], [2, 3, 3], [3, 1, 1], [1, 3, 3]]


def fixed(**kwargs):
    """A network with a pinned topology and no simulation (t_max = 0)."""
    params = dict(n_chemicals=3, paths=RING, n_division_paths=1,
                  max_cells=1, t_max=0.0)
    params.update(kwargs)
    return generate_network(ID, seed=0, **params)


def paper_derivative(model, x, X):
    """Eqs. (1)-(7) of the 1997 paper, written out path by path as printed."""
    k = model.k
    dx = np.zeros(k + 1)
    for m, l, j in model.paths:                                        # eq. (1)
        if model.enzyme == "quadratic":                                # 1998, eq. (2)
            flux = model.e1 * x[m] * x[j] ** 2
        else:
            flux = model.e1 * x[j] * x[m] / (1 + x[m] / model.x_M)
        dx[l] += flux
        dx[m] -= flux
    for l in model.source:                                             # S(l) = 1
        dx[l] += model.e0 * x[0] * x[l]
        dx[0] -= model.e0 * x[0] * x[l]                                # eq. (4)
    division_factor = 0.0
    for l in model.division:                                           # P(l) = 1
        dx[l] -= model.gamma * x[l]
        division_factor += model.gamma * x[l]
    activity = sum(x[1:])
    transport = np.array([model.transport * activity * X[m] for m in range(k + 1)])
    diffusion = np.array([model.diffusion * (X[m] - x[m]) for m in range(k + 1)])
    dx = dx + transport + diffusion                                    # eqs. (2), (3), (5)
    dX = -(transport + diffusion) / model.medium_volume                # eq. (6)
    dX[0] += model.flow * (model.nutrient - X[0])                      # eq. (7)
    dX[1:] -= model.washout * X[1:]
    return dx, division_factor, dX


def state_of(model, x, extra=0.0):
    cells = np.zeros((1, model.width))
    cells[0, : model.k + 1] = x
    cells[0, model.k + 1] = extra
    return cells


# --- the published constants ----------------------------------------------
def test_the_published_parameter_sets():
    # caption of Fig. 6 (network of Fig. 5a, the main simulation)
    assert I.PAPER == {
        "transport": 10.0, "e0": 1.0, "e1": 1.0, "nutrient": 40.0,
        "diffusion": 0.02, "gamma": 0.2, "x_M": 10.0,
        "division_threshold": 2000.0, "death_threshold": 0.05,
        "flow": 0.005, "washout": 0.005, "medium_volume": 1000.0,
        "split_noise": 0.001,
    }
    # caption of Fig. 18 (network of Fig. 5b, the simulations with cell death)
    assert I.PAPER_FIG18["nutrient"] == 10.0
    assert I.PAPER_FIG18["division_threshold"] == 100.0
    assert I.PAPER_FIG18["death_threshold"] == 0.01
    assert {k: v for k, v in I.PAPER_FIG18.items()
            if k not in ("nutrient", "division_threshold", "death_threshold")} == {
        k: v for k, v in I.PAPER.items()
        if k not in ("nutrient", "division_threshold", "death_threshold")}


def test_every_kinetic_default_is_a_published_value():
    for name in ("e0", "e1", "gamma", "x_M", "transport", "diffusion",
                 "flow", "washout", "nutrient", "split_noise"):
        assert DEFAULTS[name] == I.PAPER[name], name
    # the three deliberate deviations, all published values themselves
    assert DEFAULTS["division_threshold"] == I.PAPER_FIG18["division_threshold"] == 100.0
    assert DEFAULTS["death_threshold"] == I.PAPER_FIG18["death_threshold"] == 0.01
    assert DEFAULTS["medium_volume"] == 100.0        # V of the 1998 paper (its Fig. 1)
    assert DEFAULTS["n_chemicals"] == 8 and DEFAULTS["connections"] == 3   # Fig. 5a


# --- the network is the equations ------------------------------------------
def test_the_reactions_are_the_terms_of_equation_1():
    net = fixed()
    assert [s.id for s in net.species] == ["X0", "X1", "X2", "X3", "DF"]
    (division,) = net.extras["reaction_network"]["division_factor_paths"]
    text = net.to_text().splitlines()

    # one source path per chemical, S(l) = 1 for all of them (caption of Fig. 5)
    assert net.extras["reaction_network"]["source_paths"] == [1, 2, 3]
    for l in (1, 2, 3):
        assert f"X0 + X{l} -> 2 X{l}  [mass-action k=1.0]" in text
    # one reaction per catalytic path; here every path is autocatalytic (j = l)
    for m, l, j in RING:
        assert any(line.startswith(f"X{m} + X{j} -> 2 X{l}  [mass-action k=1.0") for line in text)
    # the chemicals with P(l) = 1 decay into the division factor at gamma
    assert f"X{division} -> DF  [mass-action k=0.2]" in text
    assert len(net.reactions) == 3 + len(RING) + 1


def test_the_saturation_is_recorded_and_not_invented():
    net = fixed()
    catalysed = [r for r in net.reactions if "DF" not in r.products and len(r.reactants) == 2
                 and "X0" not in r.reactants]
    assert len(catalysed) == len(RING)
    for r in catalysed:
        assert r.rate["law"] == "mass-action" and r.rate["k"] == 1.0
        assert r.rate["saturation"] == "michaelis-menten"
        assert r.rate["x_M"] == 10.0
        assert r.rate["saturated_species"] in r.reactants
    # the quadratic variant of the 1998 paper is exactly mass action, with no extras
    quadratic = fixed(enzyme="quadratic")
    for r in quadratic.reactions:
        assert set(r.rate) == {"law", "k"}
    assert "X1 + 2 X2 -> 3 X2  [mass-action k=1.0]" in quadratic.to_text()


def test_the_catalyst_survives_every_catalysed_path():
    # j = l: Xm + Xl -> 2 Xl, so the catalyst is a catalyst in Chemart's sense
    net = fixed()
    assert any(r.catalysts for r in net.reactions)
    other = fixed(paths=[[1, 2, 3], [2, 3, 1]])
    for r in other.reactions:
        if len(r.reactants) == 2 and "X0" not in r.reactants:
            (catalyst,) = list(r.catalysts)
            assert r.catalysts[catalyst] == 1


# --- the rate equations -----------------------------------------------------
@pytest.mark.parametrize("enzyme", I.ENZYMES)
def test_the_derivative_matches_the_equations_as_printed(enzyme):
    net = generate_network(ID, seed=4, enzyme=enzyme, n_chemicals=6, connections=2,
                           n_division_paths=3, max_cells=1, t_max=0.0,
                           network_attempts=1)
    model = I.model_from_network(net)
    rng = np.random.default_rng(0)
    x = rng.uniform(0.1, 5.0, model.k + 1)
    X = rng.uniform(0.1, 5.0, model.k + 1)

    d_cells, dX = model.derivative(state_of(model, x), X)
    want_x, want_df, want_X = paper_derivative(model, x, X)
    assert d_cells[0, : model.k + 1] == pytest.approx(want_x, rel=1e-12)
    assert d_cells[0, model.k + 1] == pytest.approx(want_df, rel=1e-12)
    assert dX == pytest.approx(want_X, rel=1e-12)


def test_the_quadratic_network_integrates_exactly_as_mass_action():
    # with the membrane closed the cell is its reaction network alone, and the
    # 1998 enzyme term e1 x(m) x(j)^2 is the mass action of Xm + 2 Xj -> Xl + 2 Xj
    net = fixed(enzyme="quadratic", transport=0.0, diffusion=0.0)
    model = I.model_from_network(net)
    ids, f = rhs(net)
    rng = np.random.default_rng(1)
    values = dict(zip(ids, rng.uniform(0.2, 3.0, len(ids))))
    x = np.array([values[f"X{i}"] for i in range(model.k + 1)])

    ours = dict(zip(ids, f(0.0, np.array([values[s] for s in ids]))))
    d_cells, _ = model.derivative(state_of(model, x), np.zeros(model.k + 1))
    for i in range(model.k + 1):
        assert ours[f"X{i}"] == pytest.approx(d_cells[0, i], rel=1e-12), f"X{i}"
    assert ours["DF"] == pytest.approx(d_cells[0, model.k + 1], rel=1e-12)


def test_mass_action_is_the_dilute_limit_of_the_recorded_saturation():
    # the emitted rate is the Michaelis-Menten form's x(m) << x_M limit
    net = fixed(transport=0.0, diffusion=0.0, x_M=1e9)
    model = I.model_from_network(net)
    ids, f = rhs(net)
    rng = np.random.default_rng(2)
    values = dict(zip(ids, rng.uniform(0.1, 2.0, len(ids))))
    x = np.array([values[f"X{i}"] for i in range(model.k + 1)])
    ours = dict(zip(ids, f(0.0, np.array([values[s] for s in ids]))))
    d_cells, _ = model.derivative(state_of(model, x), np.zeros(model.k + 1))
    for i in range(model.k + 1):
        assert ours[f"X{i}"] == pytest.approx(d_cells[0, i], rel=1e-7), f"X{i}"


def test_winner_takes_all_of_appendix_1():
    # "dx(1)/dt = x(2) x(1)/(1 + x(2)) - x(1) x(2)/(1 + x(1))" and its mirror:
    # "the difference x(1) - x(2) is amplified with time, and goes to a state
    # with either x(1) = 0 or x(2) = 0"
    net = fixed(n_chemicals=2, paths=[[1, 2, 2], [2, 1, 1]], n_division_paths=0,
                x_M=1.0, e0=0.0, gamma=0.0, transport=0.0, diffusion=0.0)
    model = I.model_from_network(net)
    x = np.array([0.0, 0.5, 0.6])
    d_cells, _ = model.derivative(state_of(model, x), np.zeros(3))
    x1, x2 = x[1], x[2]
    assert d_cells[0, 1] == pytest.approx(x2 * x1 / (1 + x2) - x1 * x2 / (1 + x1))
    assert d_cells[0, 2] == pytest.approx(x1 * x2 / (1 + x1) - x2 * x1 / (1 + x2))

    from scipy.integrate import solve_ivp
    sol = solve_ivp(model.rhs(1), (0.0, 200.0), model.pack(state_of(model, x), np.zeros(3)),
                    method="LSODA", rtol=1e-10, atol=1e-12)
    cells, _ = model.unpack(sol.y[:, -1], 1)
    end = cells[0, 1:3]
    assert min(end) == pytest.approx(0.0, abs=1e-9)        # one of them wins
    assert max(end) == pytest.approx(1.1, rel=1e-6)        # the sum is conserved
    assert abs(end[0] - end[1]) > abs(x1 - x2)             # the difference is amplified


# --- parameters -------------------------------------------------------------
def test_invalid_networks_are_refused():
    with pytest.raises(ValueError, match="connections"):
        generate_network(ID, n_chemicals=3, connections=4, t_max=0.0)
    with pytest.raises(ValueError, match="n_division_paths"):
        generate_network(ID, n_chemicals=3, connections=2, n_division_paths=4, t_max=0.0)
    with pytest.raises(ValueError, match="triple of integers"):
        generate_network(ID, paths=[[1, 2]], t_max=0.0)
    with pytest.raises(ValueError, match="outside 1..8"):
        generate_network(ID, paths=[[1, 9, 2]], t_max=0.0)
    with pytest.raises(ValueError, match="to itself"):
        generate_network(ID, paths=[[2, 2, 1]], t_max=0.0)


def test_the_drawn_network_has_the_right_shape():
    net = generate_network(ID, seed=3, n_chemicals=8, connections=3,
                           n_division_paths=4, max_cells=1, t_max=0.0,
                           network_attempts=1)
    drawn = net.extras["reaction_network"]
    assert len(drawn["paths"]) == 8 * 3
    assert len(drawn["division_factor_paths"]) == 4
    for m, l, j in drawn["paths"]:
        assert j == l and m != l                     # autocatalytic by default
    assert sorted({m for m, _, _ in drawn["paths"]}) == list(range(1, 9))
    random_catalysts = generate_network(ID, seed=3, autocatalytic=False, max_cells=1,
                                        t_max=0.0, network_attempts=1)
    assert any(j != l for _, l, j in random_catalysts.extras["reaction_network"]["paths"])


# --- the simulated cell society --------------------------------------------
def test_the_screened_network_oscillates_as_a_single_cell():
    # "Only for medium number of reaction paths, non-trivial oscillations of
    # chemicals appear as in Fig.1. We use such network for our simulation."
    net = generate_network(ID, seed=1, max_cells=1, t_max=0.0)
    a = net.extras["analysis"]
    assert a["screened"] and a["accepted"]
    assert a["single_cell_amplitude"] > I.SCREEN_AMPLITUDE
    assert a["division_factor_share"] >= I.SCREEN_DIVISION_SHARE
    assert a["stage"] == 0 and a["cells"] == 1

    # with the screening switched off the first draw is taken, whatever it does
    first = generate_network(ID, seed=1, max_cells=1, t_max=0.0, network_attempts=1)
    assert first.extras["analysis"]["screened"] is False
    assert first.extras["reaction_network"]["paths"] != net.extras["reaction_network"]["paths"]


def test_a_cell_divides_exactly_when_the_integral_of_equation_8_reaches_R():
    from scipy.integrate import solve_ivp

    net = generate_network(ID, seed=1, max_cells=1, t_max=0.0)
    model = I.model_from_network(net)
    divided = model.run(2, 400.0, np.random.default_rng(1))
    assert divided["divisions"] == 1
    when = divided["lineage"][1]["born"]

    # integrate the same cell without dividing it and find where DF crosses R
    cells, medium = model.initial(np.random.default_rng(1))
    sol = solve_ivp(model.rhs(1), (0.0, when * 1.5), model.pack(cells, medium),
                    method="LSODA", rtol=model.rtol, atol=model.atol, dense_output=True)
    factor = float(sol.sol(when)[model.k + 1])
    assert factor == pytest.approx(model.division_threshold, rel=1e-3)
    assert float(sol.sol(when * 0.9)[model.k + 1]) < model.division_threshold


def test_the_daughters_get_half_of_the_mother_and_the_split_conserves_her():
    net = generate_network(ID, seed=1, max_cells=1, t_max=0.0)
    model = I.model_from_network(net)
    when = model.run(2, 400.0, np.random.default_rng(1))["lineage"][1]["born"]

    just_after = model.run(2, when + 1e-9, np.random.default_rng(1))
    daughters = just_after["cells"][:, : model.k + 1]
    assert len(daughters) == 2
    mother = model.run(1, when, np.random.default_rng(1))["cells"][0, : model.k + 1]

    total = daughters.sum(0)
    assert total == pytest.approx(mother, rel=1e-4)                  # nothing is lost
    # each daughter gets (1/2 +- eps) of every chemical, so |x_A - x_B| = 2 eps x
    gap = np.abs(daughters[0] - daughters[1])
    assert np.all(gap <= 2 * DEFAULTS["split_noise"] * np.abs(total) + 1e-12 * np.abs(total).max())
    assert np.any(daughters[0] != daughters[1])                      # but not exactly half


def test_a_starved_cell_dies_and_hands_its_chemicals_to_the_medium():
    # eq. (9): the cell is removed when sum_l x(l) < S, and "the concentration
    # X(j) is added by x_i(j)/V at every cell death"
    net = generate_network(ID, seed=2, n_chemicals=2, paths=[[1, 2, 2], [2, 1, 1]],
                           n_division_paths=2, e0=0.0, gamma=0.5, transport=0.0,
                           diffusion=0.0, flow=0.0, washout=0.0, medium_volume=10.0,
                           death_threshold=0.5, division_threshold=1000.0,
                           max_cells=1, t_max=200.0)
    a = net.extras["analysis"]
    assert a["deaths"] == 1 and a["cells"] == 0
    released = net.extras["interaction_law"]["medium_state"]
    assert sum(released[1:]) == pytest.approx(0.5 / 10.0, rel=1e-3)


@pytest.mark.slow
def test_the_cells_divide_synchronously_and_their_number_doubles():
    # stage 1: "the cell number increases as 1, 2, 4, 8, ..." and up to a
    # threshold number "all cells have identical chemical concentrations"
    net = generate_network(ID, seed=1, max_cells=8, t_max=400.0)
    a = net.extras["analysis"]
    assert a["cells"] == 8 and a["divisions"] == 7 and a["deaths"] == 0
    assert a["n_types"] == 1 and a["stage"] == 1
    assert a["stage_name"].startswith("synchronous oscillation")
    assert a["snapshot_spread"] < a["type_tolerance"]


@pytest.mark.slow
def test_the_cells_stay_identical_only_while_the_division_is_exactly_equal():
    # the differentiation starts from the imbalance at division: "any tiny
    # difference is amplified to yield a macroscopic differentiation"
    equal = generate_network(ID, seed=1, split_noise=0.0, max_cells=8, t_max=400.0)
    assert equal.extras["analysis"]["snapshot_spread"] == 0.0
    noisy = generate_network(ID, seed=1, max_cells=8, t_max=400.0)
    assert noisy.extras["analysis"]["snapshot_spread"] > 0.0
