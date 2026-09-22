"""RBN World (Faulconbridge thesis 2011; ALife XII 2010; Faulkner et al. 2018): random Boolean
networks with bonding sites as atoms. The RBN machinery itself is tested in test_rbn.py."""

from collections import Counter

import numpy as np
import pytest

from chemart import evolve, generate_network
from chemart.chemistries.rbn_world import (
    Node, World, _preorder, _satisfied, cycle_properties, link_trees, random_element, split_tree,
)

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
    net = evolve("rbn-world", seed=4, N=10, K=2).network
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
        net = evolve("rbn-world", seed=seed).network
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


def test_frames_follow_the_soup_and_conserve_atoms():
    """A frame per elements * copies collisions; every frame's soup holds the same atoms."""
    traj = evolve("rbn-world", seed=2)
    net = traj.network
    assert traj.clock == "collisions"
    assert [f.t for f in traj.frames] == [100.0 * k for k in range(11)]
    assert traj.frames[0].state == net.initial_state and not traj.frames[0].fired
    assert traj.frames[-1].state == {s: float(n) for s, n in net.extras["analysis"]["final_population"].items()}
    for law in net.extras["conservation"]:
        totals = {sum(law["vector"][s] * n for s, n in f.state.items()) for f in traj.frames}
        assert totals == {20}, law["name"]
    fired = Counter()
    for f in traj.frames:
        for lhs, rhs, n in f.fired:
            fired[(tuple(sorted(lhs)), tuple(sorted(rhs)))] += n
    assert fired == Counter({(tuple(sorted(Counter(r.reactants).elements())),
                              tuple(sorted(Counter(r.products).elements()))): r.count for r in net.reactions})
    sizes = traj.series("largest_molecule_atoms")
    assert sizes[0] == 1 and max(sizes) <= net.extras["analysis"]["largest_molecule_atoms"]


def test_generate_network_runs_the_reactor_to_the_end():
    assert generate_network("rbn-world", seed=1).to_dict()["reactions"] == \
        evolve("rbn-world", seed=1).network.to_dict()["reactions"]


def test_invalid_atoms():
    with pytest.raises(ValueError, match="K must satisfy"):
        evolve("rbn-world", N=3, K=4)
    with pytest.raises(ValueError, match="input slots"):
        evolve("rbn-world", N=1, K=1)
    with pytest.raises(ValueError, match="unknown parameter 'K_distribution'"):
        evolve("rbn-world", K_distribution="poisson")
