"""RBN (book 18.4.2; Drossel 2008). RBN World, built from the same machinery, is tested in test_rbn_world.py."""

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.rbn import attractors, in_degrees, random_rbn, successor_map


# --- classic RBN ----------------------------------------------------------------
def test_reactions_are_the_truth_tables():
    """Book 18.4.2: next state of node i = f_i(inputs). In every state, the enabled row reactions switch exactly the nodes whose output differs."""
    N = 6
    net = generate_network("rbn", seed=3, N=N, K=3)
    a = net.extras["analysis"]
    for x in range(2 ** N):
        present = {f"x{i}_{(x >> i) & 1}" for i in range(N)}
        switched = {}
        for r in net.reactions:
            if set(r.reactants) <= present:
                (old,) = set(r.reactants) - set(r.products)
                (new,) = set(r.products) - set(r.reactants)
                assert old not in switched, "one truth-table row applies per node"
                switched[old] = new
        for i in range(N):
            row = int("".join(str((x >> j) & 1) for j in a["inputs"][i]), 2)
            out, now = int(a["functions"][i][row]), (x >> i) & 1
            assert switched.get(f"x{i}_{now}") == (f"x{i}_{out}" if out != now else None)


def test_attractors_are_cycles_of_the_synchronous_map():
    N = 8
    net = generate_network("rbn", seed=5, N=N, K=2)
    a = net.extras["analysis"]
    f = successor_map(N, a["inputs"], [[int(c) for c in t] for t in a["functions"]])
    assert sum(att["basin_size"] for att in a["attractors"]) == 2 ** N
    for att in a["attractors"]:
        states = [int(s[::-1], 2) for s in att["cycle"]]
        assert len(states) == att["length"]
        assert all(f[s] == t for s, t in zip(states, states[1:] + states[:1]))
    x = sum(int(name.split("_")[1]) << int(name[1:].split("_")[0]) for name in net.initial_state)
    for _ in range(2 ** N):
        x = int(f[x])
    assert format(x, f"0{N}b")[::-1] in a["attractors"][a["initial_attractor"]]["cycle"]


def _mean_attractor_length(N, K, seed):
    inputs, tables = random_rbn(N, [K] * N, 0.5, np.random.default_rng(seed))
    found = attractors(successor_map(N, inputs, tables))
    return sum(x["length"] * x["basin_size"] for x in found) / 2 ** N


def test_order_to_chaos_as_K_grows():
    """Book 18.4.2: K = 1 frozen, K = 2 'edge of chaos', large K chaotic (attractors of order 2^(N/2), Drossel)."""
    N = 12
    median = {K: float(np.median([_mean_attractor_length(N, K, s) for s in range(15)])) for K in (1, 2, N)}
    assert median[1] <= median[2] < median[N]
    assert median[N] > 10 * median[1]
    assert median[N] > 2 ** (N / 2) / 8


@pytest.mark.parametrize("K, bias", [(2, 0.5), (3, 0.2)])
def test_sensitivity_matches_annealed_approximation(K, bias):
    """Drossel [243] eq. 10: lambda = 2 K p (1 - p); K = 2, p = 1/2 is critical."""
    values = [generate_network("rbn", seed=s, N=12, K=K, function_bias=bias).extras["analysis"]["average_sensitivity"]
              for s in range(30)]
    assert np.mean(values) == pytest.approx(2 * K * bias * (1 - bias), abs=0.08)


def test_k_distributions_have_mean_K():
    """Book 18.4.2: K replaced by probability distributions with mean K."""
    rng = np.random.default_rng(0)
    for dist in ("poisson", "power-law"):
        ks = [k for _ in range(300) for k in in_degrees(16, 3, dist, rng)[0]]
        assert np.mean(ks) == pytest.approx(3, abs=0.15)
        assert 1 <= min(ks) and max(ks) <= 16
    net = generate_network("rbn", seed=2, N=12, K=3, K_distribution="power-law")
    assert net.extras["analysis"]["power_law_exponent"] > 0


def test_invalid_combinations():
    with pytest.raises(ValueError, match="K must satisfy"):
        generate_network("rbn", N=4, K=5)
    with pytest.raises(ValueError, match="power law"):
        generate_network("rbn", N=6, K=4, K_distribution="power-law")
    # RBN World is its own entry now
    with pytest.raises(ValueError, match="'model' is gone"):
        generate_network("rbn", model="rbn-world")
