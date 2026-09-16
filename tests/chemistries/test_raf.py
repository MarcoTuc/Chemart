"""RAF theory (book 6.3.1): the published algorithm, examples and the level of
catalysis needed for autocatalytic sets to appear.

Sources, all applied to the binary polymer model with food set t = 2:

- Hordijk, Kauffman & Steel, Int. J. Mol. Sci. 12:3085-3101 (2011): the RAF
  definition and algorithm (secs. 3.1-3.2) and Table 1, f(n) = 1.0970 + 0.0189 n
  at P_n = 0.50;
- Hordijk, Hein & Steel, Entropy 12:1733-1742 (2010): "between 1 and 2 reactions
  per molecule";
- Steel, Hordijk & Smith, J. Theor. Biol. 332:96-107 (2013), Fig. 4: at n = 10 no
  RAFs below f = 1.20, maxRAF ~1222 reactions when they first appear and growing
  with f, irrRAF size roughly constant;
- Hordijk, Smith & Steel, Algorithms Mol. Biol. 10:15 (2015): p = 0.00041 (n = 8)
  and p = 0.0000792 (n = 10) give P_n ~ 0.5, average maxRAF 375 reactions at
  n = 8; the linear-chain and pseudo-RAF examples;
- Hordijk & Steel, arXiv:1206.1017 (2012), Fig. 1: a 4-reaction system whose
  maxRAF is the 2 mutually catalysed reactions.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.raf import (binary_polymer_model, closure, irreducible_raf,
                                     is_raf, max_raf, molecules_of)

# |R| for the binary polymer model, and the food set of all strings up to length 2.
REACTIONS = {n: (n - 2) * 2 ** (n + 1) + 4 for n in range(2, 15)}
FOOD2 = ["0", "1", "00", "01", "10", "11"]


def reaction(rid, reactants, products, catalysts=(), reversible=False):
    return {"id": rid, "reactants": list(reactants), "products": list(products),
            "catalysts": list(catalysts), "reversible": reversible}


def instance(n, f, seed, t=2, B=2):
    """One random instance of the binary polymer model."""
    p = SimpleNamespace(n=n, t=t, B=B, f=f)
    _, reactions, food, prob = binary_polymer_model(p, np.random.default_rng(seed))
    return reactions, food, prob


def probability_of_raf(n, f, trials, offset=0):
    """P_n: the fraction of model instances that contain an RAF, and the mean
    maxRAF size over those that do."""
    sizes = []
    for s in range(trials):
        reactions, food, _ = instance(n, f, 7717 * n + offset + s)
        raf = max_raf(reactions, food)
        if raf:
            sizes.append(len(raf))
    return len(sizes) / trials, float(np.mean(sizes)) if sizes else 0.0


# --- the algorithm on systems with a known answer ---------------------------------
def test_mutual_catalysis_example_of_hordijk_steel_2012_figure_1():
    """Four ligations from the food set {00, 01, 10, 11}: the two whose products
    catalyse each other are an RAF, the two uncatalysed ones are in no RAF."""
    system = [
        reaction("r1", ["00", "01"], ["0001"], ["1011"], reversible=True),
        reaction("r2", ["10", "11"], ["1011"], ["0001"], reversible=True),
        reaction("r3", ["00", "11"], ["0011"], reversible=True),
        reaction("r4", ["01", "10"], ["0110"], reversible=True),
    ]
    food = ["00", "01", "10", "11"]
    raf = max_raf(system, food)
    assert sorted(r["id"] for r in raf) == ["r1", "r2"]
    assert is_raf(raf, food)
    # Neither half survives on its own: each catalyst is the other's product.
    assert max_raf([system[0]], food) == []
    assert max_raf([system[1]], food) == []


def test_linear_chain_of_algorithms_paper_has_no_raf():
    """x_{i-1} -> x_i catalysed by x_{i+1}, with the last reaction uncatalysed:
    the algorithm peels the chain from the top and reports no RAF."""
    N = 12
    chain = [reaction(f"r{i}", [f"x{i - 1}"], [f"x{i}"], [f"x{i + 1}"] if i < N else [])
             for i in range(1, N + 1)]
    assert max_raf(chain, ["x0"]) == []
    # Catalyse that last reaction from inside the chain and the whole chain is an RAF.
    chain[-1]["catalysts"] = ["x1"]
    assert len(max_raf(chain, ["x0"])) == N


def test_pseudo_raf_is_not_f_generated_so_is_not_a_raf():
    """Figure 2 of the algorithms paper: every reaction needs a food molecule AND
    another product, so the system is reflexively autocatalytic in the F u pi(R')
    sense but cannot get started from F; it is a pseudo-RAF, not an RAF."""
    food = ["f1", "f2", "f3"]
    system = [
        reaction("r1", ["f1", "p3"], ["p1"], ["p2"]),
        reaction("r2", ["f2", "p1"], ["p2"], ["p3"]),
        reaction("r3", ["f3", "p2"], ["p3"], ["p1"]),
    ]
    products = {"p1", "p2", "p3"}
    assert all(set(r["reactants"]) <= set(food) | products for r in system)   # pseudo-RAF
    assert all(set(r["catalysts"]) <= set(food) | products for r in system)
    assert closure(food, system) == set(food)                                # nothing can fire
    assert max_raf(system, food) == []


def test_f_generated_without_catalysis_is_not_a_raf():
    """A chain that runs from the food set but catalyses nothing: F-generated,
    not reflexively autocatalytic."""
    chain = [reaction("r1", ["a"], ["b"]), reaction("r2", ["b"], ["c"])]
    assert closure(["a"], chain) == {"a", "b", "c"}          # F-generated
    assert max_raf(chain, ["a"]) == []                       # but no catalysts
    chain[1]["catalysts"] = ["c"]                            # catalyse only the second step
    assert max_raf(chain, ["a"]) == [], "r2 alone cannot make its own reactant b"
    chain[0]["catalysts"] = ["c"]                            # the end product catalyses both
    assert sorted(r["id"] for r in max_raf(chain, ["a"])) == ["r1", "r2"]


def test_maxraf_is_a_fixed_point_and_every_reaction_is_supported():
    reactions, food, _ = instance(7, 1.5, 11)
    raf = max_raf(reactions, food)
    assert raf, "this instance should contain an RAF"
    assert len(max_raf(raf, food)) == len(raf)               # running it again changes nothing
    assert is_raf(raf, food)
    cl = closure(food, raf)
    for r in raf:
        assert molecules_of(r) <= cl                          # F-generated
        assert any(c in cl for c in r["catalysts"])           # reflexively autocatalytic
    assert set(food) <= cl


def test_irreducible_raf_is_a_minimal_subraf():
    reactions, food, _ = instance(7, 1.5, 1)
    raf = max_raf(reactions, food)
    rng = np.random.default_rng(3)
    irr = irreducible_raf(reactions, food, rng)
    assert 0 < len(irr) < len(raf)
    assert {r["id"] for r in irr} <= {r["id"] for r in raf}
    assert is_raf(irr, food)
    for r in irr:                                             # no proper subset is an RAF
        assert max_raf([x for x in irr if x["id"] != r["id"]], food) == []


def test_different_runs_find_different_irreducible_rafs():
    """The deletion search visits the reactions in a random order, so a maxRAF
    with many irrRAFs gives different answers on different runs."""
    reactions, food, _ = instance(7, 1.5, 1)
    found = {tuple(sorted(r["id"] for r in irreducible_raf(reactions, food, np.random.default_rng(s))))
             for s in range(5)}
    assert len(found) > 1


# --- the published probability of finding an RAF ----------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("n, p_cat, published_P, published_size", [
    (8, 0.00041, 0.5, 375),          # Algorithms Mol. Biol. 10:15, "Performance" and "irrRAFs"
    (10, 0.0000792, 0.5, None),      # ibid., inhibition simulations
])
def test_published_catalysis_probabilities_give_a_coin_flip(n, p_cat, published_P, published_size):
    f = p_cat * REACTIONS[n]
    P, size = probability_of_raf(n, f, 40 if n == 8 else 30)
    assert abs(P - published_P) < 0.25
    if published_size is not None:
        assert 0.6 * published_size < size < 1.6 * published_size


@pytest.mark.slow
@pytest.mark.parametrize("n", [7, 8, 9])
def test_table_1_line_gives_probability_one_half(n):
    """Hordijk, Kauffman & Steel (2011), Table 1 case A: the level of catalysis at
    which half the instances contain an RAF is f(n) = 1.0970 + 0.0189 n."""
    f = 1.0970 + 0.0189 * n
    P, _ = probability_of_raf(n, f, 40)
    assert 0.2 <= P <= 0.8


@pytest.mark.slow
def test_only_one_to_two_reactions_per_molecule_are_needed_whatever_the_size():
    """The headline result (Entropy 12:1733): RAFs appear at 1-2 catalysed
    reactions per molecule, and that level does not grow with system size - here
    the molecule set grows eightfold from n = 7 to n = 10 while the transition
    stays inside [1, 2]."""
    low, high = {}, {}
    for n in (7, 8, 9, 10):
        low[n], _ = probability_of_raf(n, 1.0, 15, offset=101)
        high[n], _ = probability_of_raf(n, 2.0, 15, offset=202)
    assert all(P < 0.5 for P in low.values()), low         # f = 1 is below the transition
    assert all(P == 1.0 for P in high.values()), high      # f = 2 is above it, for every n
    assert low[10] <= low[7]                               # and it gets sharper with n


@pytest.mark.slow
def test_transition_at_n_ten_matches_the_minimal_networks_paper():
    """Steel, Hordijk & Smith (2013), Fig. 4: at n = 10 no RAFs at all below
    f = 1.20; they appear just above it and their size grows with f (about 1222
    reactions when they first appear, about 2000 by f = 1.45)."""
    assert probability_of_raf(10, 1.15, 12)[0] == 0.0
    P_low, size_low = probability_of_raf(10, 1.25, 20)
    P_high, size_high = probability_of_raf(10, 1.45, 20)
    assert P_low < P_high and P_high > 0.8
    assert 0.6 * 1222 < size_low < 1.6 * 1222
    assert size_high > size_low
    assert 0.6 * 2000 < size_high < 1.6 * 2000


@pytest.mark.slow
def test_maxraf_grows_with_catalysis_but_irreducible_rafs_do_not():
    """Steel, Hordijk & Smith (2013), Fig. 4: raising the level of catalysis makes
    the maximal RAF bigger while the size of an irreducible RAF stays flat."""
    means = {}
    for f in (1.4, 1.8):
        maxes, irrs = [], []
        for s in range(6):
            reactions, food, _ = instance(8, f, 248 + s)
            raf = max_raf(reactions, food)
            if not raf:
                continue
            irr = irreducible_raf(reactions, food, np.random.default_rng(s))
            assert is_raf(irr, food) and len(irr) <= len(raf)
            maxes.append(len(raf))
            irrs.append(len(irr))
        means[f] = (float(np.mean(maxes)), float(np.mean(irrs)))
    assert means[1.8][0] > 1.25 * means[1.4][0], means      # maxRAF grows with f
    assert means[1.8][1] < 1.6 * means[1.4][1], means       # irrRAF does not
    assert means[1.4][1] < 0.6 * means[1.4][0]              # and is a small part of it


@pytest.mark.slow
def test_irreducible_raf_size_at_n_ten():
    """Steel, Hordijk & Smith (2013), sec. 7: at n = 10 near the transition the
    maxRAFs hold about 1222 reactions and the irrRAFs inside them about 624."""
    reactions, food, _ = instance(10, 1.35, 310)
    raf = max_raf(reactions, food)
    irr = irreducible_raf(reactions, food, np.random.default_rng(0))
    assert is_raf(irr, food)
    assert 0.5 * 1222 < len(raf) < 2.0 * 1222
    assert 0.5 * 624 < len(irr) < 1.8 * 624


# --- the generated network --------------------------------------------------------
def test_default_network_is_the_polymer_model_with_its_raf_analysis():
    net = generate_network("raf", seed=1)
    analysis = net.extras["analysis"]
    assert len(net.species) == 2 ** 8 - 2                       # all binary strings up to length 7
    assert analysis["reactions_total"] == REACTIONS[7] == 1284
    assert net.extras["food"] == FOOD2 and net.extras["buffered"] == FOOD2
    assert analysis["catalysis_probability"] == pytest.approx(1.5 / 1284)
    assert abs(analysis["catalysis_events"] - 1.5 * len(net.species)) < 5 * np.sqrt(1.5 * len(net.species))
    assert analysis["raf_exists"] and analysis["max_raf_size"] > 0
    assert len(analysis["max_raf_reactions"]) == analysis["max_raf_size"]
    assert analysis["irreducible_rafs"][0]["size"] <= analysis["max_raf_size"]
    # every network reaction is a catalysed one: the catalyst is on both sides
    assert net.reactions and all(r.catalysts for r in net.reactions)
    assert len(net.reactions) == 2 * analysis["catalysis_events"]
    assert net.status == "complete" and all(r.rate is None for r in net.reactions)


def test_ligation_and_cleavage_are_mass_conserving():
    net = generate_network("raf", seed=2, n=5)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    assert net.extras["conservation"]
    for law in net.extras["conservation"]:
        assert not np.any(np.array([law["vector"].get(s, 0) for s in ids]) @ S)


def test_all_reactions_mode_adds_the_uncatalysed_pairs():
    catalysed = generate_network("raf", seed=1, n=6)
    every = generate_network("raf", seed=1, n=6, reactions="all")
    assert len(every.reactions) == len(catalysed.reactions) + 2 * REACTIONS[6]
    assert every.extras["analysis"] == catalysed.extras["analysis"]   # the analysis is the same


def test_analysis_of_a_user_supplied_system():
    """The same 4-reaction example, handed in as a reaction system."""
    system = {"food": ["00", "01", "10", "11"], "reactions": [
        {"id": "r1", "reactants": ["00", "01"], "products": ["0001"], "catalysts": ["1011"], "reversible": True},
        {"id": "r2", "reactants": ["10", "11"], "products": ["1011"], "catalysts": ["0001"], "reversible": True},
        {"id": "r3", "reactants": ["00", "11"], "products": ["0011"], "reversible": True},
        {"id": "r4", "reactants": ["01", "10"], "products": ["0110"], "reversible": True}]}
    net = generate_network("raf", system=system, reactions="all")
    analysis = net.extras["analysis"]
    assert analysis["max_raf_reactions"] == ["r1", "r2"]
    assert analysis["max_raf_molecules"] == ["00", "0001", "01", "10", "1011", "11"]
    assert analysis["irreducible_rafs"] == [{"size": 2, "reactions": ["r1", "r2"]}]
    assert "level_of_catalysis" not in analysis            # no polymer model behind it
    assert len(net.species) == 8 and len(net.reactions) == 12   # 4 pairs + 2 catalysed pairs
    assert "conservation" not in net.extras


def test_same_seed_same_analysis_and_bad_parameters_are_rejected():
    a = generate_network("raf", seed=5, n=6)
    b = generate_network("raf", seed=5, n=6)
    assert a.to_dict() == b.to_dict()
    assert generate_network("raf", seed=6, n=6).extras["analysis"] != a.extras["analysis"]
    with pytest.raises(ValueError, match="food set length"):
        generate_network("raf", n=3, t=4)
    with pytest.raises(ValueError, match="non-empty list"):
        generate_network("raf", system={"food": ["a"], "reactions": []})
    with pytest.raises(ValueError, match="invalid species ids"):
        generate_network("raf", system={"food": ["a"], "reactions": [{"reactants": ["a b"], "products": ["c"]}]})
