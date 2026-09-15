"""RBN (book 18.4.2; Drossel 2008) and RBN World (Faulconbridge thesis 2011; ALife XII 2010; Faulkner et al. 2018)."""

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries.rbn import (
    Node, World, _preorder, _satisfied, attractors, cycle_properties, in_degrees, link_trees,
    random_element, random_rbn, split_tree, successor_map,
)


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
    with pytest.raises(ValueError, match="K_distribution"):
        generate_network("rbn", model="rbn-world", K_distribution="poisson")


# --- RBN World --------------------------------------------------------------------
def test_bonding_properties_reproduce_table_2():
    """ALife XII table 2 / thesis table 7.1.2: the n = 4 example cycle."""
    T, F = 1, 0
    cycle = [[F, T, F, F], [F, F, F, F], [T, T, F, F], [T, T, T, F], [T, F, T, F], [F, T, T, F]]
    props = cycle_properties(cycle)
    assert {k: v for k, v in props.items() if k != "proportion"} == {
        "cycle-length": 6, "flashing": 3, "flashes": 8, "total": -4, "magnitude": 14}
    assert round(props["proportion"], 3) == 0.417


def test_brbn_atoms_have_two_bonding_sites():
    """Faulkner et al. 2018 sec. 5.1.2: each of b = 2 linking nodes replaces one input of a random node."""
    net = generate_network("rbn", seed=4, model="rbn-world", N=10, K=2)
    for el in net.extras["elements"]:
        refs = [(i, r) for i, ins in enumerate(el["inputs"]) for r in ins]
        for s in (0, 1):
            (node,) = [i for i, r in refs if r == f"site{s}"]
            assert el["bonding_site_targets"][s] == node
        assert all(len(ins) == 2 and len(set(ins)) == 2 for ins in el["inputs"])


def _leaf(ch):
    return Node(state={ord(ch): 0}, atom=ord(ch))


def _comp(*kids):
    return Node(kids, {a: 0 for k in kids for a in k.atoms})


def _shape(n):
    return chr(next(iter(n.atoms))) if not n.children else "(" + "".join(sorted(_shape(c) for c in n.children)) + ")"


def test_linking_zips_structure_trees():
    # Faulkner et al. 2018 fig. 12: (A(B(CDE))) + (((FG)H)J) -> ((A(B(CDE)FG)H)J)
    cde, F = _comp(_leaf("C"), _leaf("D"), _leaf("E")), _leaf("F")
    left = _comp(_leaf("A"), _comp(_leaf("B"), cde))
    right = _comp(_comp(_comp(F, _leaf("G")), _leaf("H")), _leaf("J"))
    assert [_shape(r) for r in link_trees([left, right], cde, F)] == [_shape(
        _comp(_comp(_comp(_leaf("B"), _comp(_leaf("C"), _leaf("D"), _leaf("E")), _leaf("F"), _leaf("G")), _leaf("A"), _leaf("H")), _leaf("J")))]
    # thesis 6.2.2: (A-B) + C bonding at B -> (A-B-C); 6.2.3: at (A-B) -> ((A-B)-C)
    B, C = _leaf("B"), _leaf("C")
    assert [_shape(r) for r in link_trees([_comp(_leaf("A"), B), C], B, C)] == ["(ABC)"]
    ab, C = _comp(_leaf("A"), _leaf("B")), _leaf("C")
    assert [_shape(r) for r in link_trees([ab, C], ab, C)] == ["((AB)C)"]
    # thesis 6.2.4: ((A-B)-(C-D)); 6.2.5: (A-(B-C-D)) + ((E-F)-G) -> ((A-(B-C-D)-E-F)-G)
    ab, cd = _comp(_leaf("A"), _leaf("B")), _comp(_leaf("C"), _leaf("D"))
    assert [_shape(r) for r in link_trees([ab, cd], ab, cd)] == ["((AB)(CD))"]
    bcd, E = _comp(_leaf("B"), _leaf("C"), _leaf("D")), _leaf("E")
    trees = [_comp(_leaf("A"), bcd), _comp(_comp(E, _leaf("F")), _leaf("G"))]
    assert [_shape(r) for r in link_trees(trees, bcd, E)] == ["(((BCD)AEF)G)"]


def test_breaking_bonds_splits_trees():
    """Thesis 6.2.8: (A-B-C) loses A-B, then B-C: A + (B-C), then A + B + C."""
    A, B, C = ord("A"), ord("B"), ord("C")
    bonds = {(A, 1): (B, 0), (B, 0): (A, 1), (B, 1): (C, 0), (C, 0): (B, 1)}
    abc = _comp(_leaf("A"), _leaf("B"), _leaf("C"))
    assert split_tree(abc, bonds) == [abc]
    del bonds[(A, 1)], bonds[(B, 0)]
    pieces = split_tree(abc, bonds)
    assert sorted(_shape(p) for p in pieces) == ["(BC)", "A"]
    del bonds[(B, 1)], bonds[(C, 0)]
    assert sorted(_shape(q) for p in pieces for q in split_tree(p, bonds)) == ["A", "B", "C"]


@pytest.mark.parametrize("rule", ["cycle-length-equal", "proportion-sum-one"])
def test_bonds_exist_only_while_the_criterion_holds(rule):
    """Thesis ch. 6: 'bonds only exist between bRBNs whose emergent bonding property meets the bonding criterion'."""
    rng, N, atoms = np.random.default_rng(4), 6, 24
    elements = [random_element(N, 2, 0.5, rng) for _ in range(3)]
    world = World(elements, N, rule)
    for a in range(atoms):
        world.elem[a] = a % 3
        leaf = Node(state={a: elements[a % 3].initial}, atom=a)
        world.roots.append(leaf)
        world.settle(leaf)
    for _ in range(300):
        i, j = rng.choice(len(world.roots), size=2, replace=False)
        free = [[(a, s) for a in sorted(r.atoms) for s in (0, 1) if (a, s) not in world.bonds]
                for r in (world.roots[i], world.roots[j])]
        world.react(*(f[int(rng.integers(len(f)))] for f in free))
    assert sorted(a for r in world.roots for a in r.atoms) == list(range(atoms))
    if rule == "cycle-length-equal":
        assert max(len(r.atoms) for r in world.roots) > 2
    for root in world.roots:
        inside = [k for k in world.bonds if k[0] in root.atoms]
        assert all(world.bonds[k][0] in root.atoms for k in inside)
        assert len(inside) // 2 == len(root.atoms) - 1, "bonds join distinct molecules, so each molecule is a tree"
        for T in _preorder(root):
            owner = {a: c for c in T.children for a in c.atoms}
            for (a, s), (b, t) in world.bonds.items():
                if a in owner and b in owner and owner[a] is not owner[b]:
                    assert _satisfied(rule, world.settle(owner[a]), world.settle(owner[b]))


def test_world_network_conserves_atoms_and_shows_synthesis_and_catalysis():
    """Thesis 6.2: synthesis A + B -> (A-B) and emergent catalysis A1 + B -> A2 + B; atoms are conserved."""
    kinds = set()
    for seed in range(6):
        net = generate_network("rbn", seed=seed, model="rbn-world")
        assert net.status == "observed"
        ids, R, P = net.matrices()
        S = (P - R).toarray()
        for law in net.extras["conservation"]:
            m = np.array([law["vector"][s] for s in ids])
            assert np.all(m @ S == 0)
        for r in net.reactions:
            assert sum(r.reactants.values()) == 2 and r.count >= 1
            n_out = sum(r.products.values())
            kinds.add("synthesis" if n_out == 1 else "catalysis" if r.catalysts else "other")
    assert {"synthesis", "catalysis"} <= kinds
