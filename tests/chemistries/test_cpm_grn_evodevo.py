"""Cellular Potts + GRN evo-devo (book 18.6.1).

Reference values come from the equations of the papers themselves:

- Savill & Hogeweg (1997), J. theor. Biol. 184:229-235, eqns (1)-(2):
  H = sum J_cell,cell / 2 + sum J_cell,medium + lambda (v - V)^2, and
  p = 1 if dH <= -0.1 else exp(-(dH + 0.1)).
- Hogeweg (2000), J. theor. Biol. 203:317-333, eqn (1) (the same Hamiltonian),
  eqns (2a)-(2b) (lock-and-key adhesion), section 2.2 (the Boolean network and
  the maternal factors), section 3.1 (the Glazier & Graner engulfment
  condition) and Table 1 (the parameters and the developmental schedule).
- Hogeweg (2000), Artificial Life 6(1):85-101, section 2 (lambda = 0.5 chosen
  low so a squeezed cell can die) and section 3.2 (the morphogenetic
  mechanisms).

The Hamiltonians below are counted by hand on blocks small enough to check.
"""

import copy
import math
from functools import lru_cache
from types import SimpleNamespace

import numpy as np
import pytest

from chemart import catalog, generate_network
from chemart.chemistries import cpm_grn_evodevo as M

ID = "cpm-grn-evodevo"


@lru_cache(maxsize=1)
def _entry():
    return next(c for c in catalog.load() if c.id == ID)


def params(**over):
    """The catalog defaults as the generator sees them."""
    values = {p.name: copy.deepcopy(p.default) for p in _entry().params}
    values.update(over)
    return SimpleNamespace(**values)


DEFAULTS = params()
_CACHE: dict[int, object] = {}


def default_net(seed=1):
    if seed not in _CACHE:
        _CACHE[seed] = generate_network(ID, seed=seed)
    return _CACHE[seed]


def block(cpm, rows, cols):
    return [r * cpm.width + c for r in rows for c in cols]


# ===========================================================================
# the Cellular Potts Model: H, dH and the Metropolis rule
# ===========================================================================

def test_hamiltonian_counts_every_unlike_neighbour_pair_once():
    """H = sum of J over neighbour pairs with different sigma + lambda (v - V)^2.

    A 2x2 cell in the middle of a lattice has 20 Moore neighbour pairs with the
    medium: the four sites have 8 neighbours each, 3 of them inside the block.
    """
    cpm = M.CPM(7, 7, temperature=1.0, inelasticity=1.0, bond=lambda a, b: 3.0)
    cpm.place(1, block(cpm, (3, 4), (3, 4)), 4)
    pairs = sum(1 for x in range(cpm.size) for y in cpm.nbrs[x]
                if cpm.sigma[x] != cpm.sigma[y]) // 2
    assert pairs == 20
    assert cpm.hamiltonian() == 20 * 3.0            # v = V, so no volume term

    cpm.target[1] = 6                                # lambda (4 - 6)^2 = 4
    assert cpm.hamiltonian() == 20 * 3.0 + 4.0
    cpm.inelasticity = 0.5                           # Hogeweg's lambda
    assert cpm.hamiltonian() == 20 * 3.0 + 2.0


def test_hamiltonian_of_two_touching_cells():
    """Cell-cell and cell-medium bonds are both J, counted once per pair."""
    cpm = M.CPM(9, 9, temperature=1.0, inelasticity=1.0,
                bond=lambda a, b: 3.0 if M.MEDIUM in (a, b) else 2.0)
    cpm.place(1, block(cpm, (3, 4), (3, 4)), 4)
    cpm.place(2, block(cpm, (3, 4), (5, 6)), 4)
    # the two 2x2 cells form a 2x4 rectangle: 4 corner sites with 5 medium
    # neighbours and 4 edge sites with 3, so 32 cell-medium pairs; the four
    # sites either side of the seam give 4 cell-cell pairs
    assert cpm.hamiltonian() == 32 * 3.0 + 4 * 2.0 == 104.0


def test_delta_h_equals_a_full_recomputation_of_h():
    """The incremental dH is exactly H(after) - H(before) for every candidate copy."""
    cpm = M.CPM(9, 9, temperature=1.0, inelasticity=0.5,
                bond=lambda a, b: 3.0 if M.MEDIUM in (a, b) else 2.0)
    cpm.place(1, block(cpm, (3, 4), (3, 4)), 5)      # v != V, so the volume term bites
    cpm.place(2, block(cpm, (3, 4), (5, 6)), 3)
    checked = 0
    for site in cpm.sites:
        for neighbour in cpm.nbrs[site]:
            new, old = cpm.sigma[neighbour], cpm.sigma[site]
            if new == old:
                continue
            before = cpm.hamiltonian()
            predicted = cpm.delta_h(site, new)
            cpm.sigma[site] = new                     # apply
            if old != M.MEDIUM:
                cpm.volume[old] -= 1
            if new != M.MEDIUM:
                cpm.volume[new] += 1
            assert cpm.hamiltonian() - before == pytest.approx(predicted)
            cpm.sigma[site] = old                     # undo
            if old != M.MEDIUM:
                cpm.volume[old] += 1
            if new != M.MEDIUM:
                cpm.volume[new] -= 1
            checked += 1
    assert checked == 72                              # every candidate copy in the fixture
    assert cpm.delta_h(cpm.sites[0], cpm.sigma[cpm.sites[0]]) == 0.0


def test_acceptance_is_the_published_metropolis_rule():
    """Savill & Hogeweg eqn (2): p = 1 if dH < -0.1, else exp(-(dH + 0.1) / T)."""
    cpm = M.CPM(5, 5, temperature=3.0, inelasticity=1.0, bond=lambda a, b: 1.0)
    assert M.DELTA_H_OFFSET == 0.1
    assert cpm.acceptance(-5.0) == 1.0
    assert cpm.acceptance(-0.2) == 1.0
    # continuous at the switch: exp(0) = 1
    assert cpm.acceptance(-0.1) == pytest.approx(1.0)
    assert cpm.acceptance(0.0) == pytest.approx(math.exp(-0.1 / 3.0))
    assert cpm.acceptance(5.0) == pytest.approx(math.exp(-5.1 / 3.0))
    # hotter means more tolerant of an uphill move
    hot = M.CPM(5, 5, temperature=10.0, inelasticity=1.0, bond=lambda a, b: 1.0)
    assert hot.acceptance(5.0) > cpm.acceptance(5.0)


def test_division_cuts_across_the_longest_axis():
    """Hogeweg Table 1: 'perpendicular and halfway on longest axis of cell'."""
    cpm = M.CPM(14, 14, temperature=1.0, inelasticity=1.0, bond=lambda a, b: 1.0)
    cpm.place(1, block(cpm, (5, 6), range(3, 11)), 16)   # 2 rows x 8 columns
    assert cpm.divide(1, 2)
    left = [x % cpm.width for x in cpm.pixels(1)]
    right = [x % cpm.width for x in cpm.pixels(2)]
    assert cpm.volume[1] == cpm.volume[2] == 8           # halved
    assert max(left) < min(right)                        # the cut is across the long axis
    assert sorted({x // cpm.width for x in cpm.pixels(2)}) == [5, 6]   # both rows kept


# ===========================================================================
# the volume constraint
# ===========================================================================

def lone_cell(volume, j_medium, *, grid=40, steps=60, inelasticity=0.5,
              temperature=3.0, seed=0):
    """Final volume of a single cell facing the medium at surface energy j_medium."""
    rng = np.random.default_rng(seed)
    cpm = M.CPM(grid, grid, temperature=temperature, inelasticity=inelasticity,
                bond=lambda a, b: float(j_medium))
    side = max(1, round(math.sqrt(volume)))
    start = (grid - side) // 2
    cpm.place(1, block(cpm, range(start, start + side), range(start, start + side)), volume)
    for _ in range(steps):
        cpm.sweep(rng)
        if not cpm.volume.get(1):
            return 0
    return cpm.volume[1]


@pytest.mark.parametrize("j_medium", [0, 4, 8, 16, 24, 31])
def test_the_volume_constraint_holds_a_cell_near_its_target(j_medium):
    """lambda (v - V)^2 balances the surface energy, over the whole J range a GRN can express."""
    v = lone_cell(DEFAULTS.target_volume, j_medium)
    assert 0.7 * DEFAULTS.target_volume <= v <= 1.2 * DEFAULTS.target_volume


def test_a_squeezed_cell_dies_which_is_why_lambda_is_low():
    """Hogeweg: lambda = 0.5 is chosen low 'so that a cell can die because its volume goes to zero'."""
    # the same cell survives against a cheap medium and evaporates against a dear one
    assert lone_cell(20, 0) > 10
    assert lone_cell(20, 16) == 0
    # a stiffer cell resists: raising lambda saves it
    assert lone_cell(20, 16, inelasticity=8.0) > 10


# ===========================================================================
# lock-and-key adhesion (Hogeweg eqns 2a-2b)
# ===========================================================================

def test_lock_and_key_adhesion_matches_the_published_equations():
    a = M.Adhesion(DEFAULTS.adhesion_nodes)
    assert a.locks == (0, 1, 2, 3, 4) and a.keys == (5, 6, 7, 8, 9)
    assert a.medium_nodes == (3, 4, 5, 6, 7)          # straddles locks and keys
    assert a.maximum == 31.0                          # five bits

    nodes = DEFAULTS.grn_nodes
    off = (0,) * nodes
    on = (1,) * nodes
    # no key is complementary to its lock when both patterns are equal
    assert a.between(off, off) == 0.0
    assert a.between(on, on) == 0.0
    # every key of one is complementary to every lock of the other
    assert a.between(off, on) == 31.0 == a.between(on, off)
    assert a.medium(off) == 0.0 and a.medium(on) == 31.0

    # J is symmetric, and inside the range the five-bit place value allows
    rng = np.random.default_rng(0)
    for _ in range(200):
        x = tuple(int(b) for b in rng.integers(0, 2, nodes))
        y = tuple(int(b) for b in rng.integers(0, 2, nodes))
        assert a.between(x, y) == a.between(y, x)
        assert 0.0 <= a.between(x, y) <= 31.0
        assert 0.0 <= a.medium(x) <= 31.0
        assert (2 * a.between(x, y)) % 1 == 0          # a sum of two integers, halved


def test_a_cell_can_have_a_nonzero_adhesion_with_its_own_type():
    """Hogeweg: 'It is done in this way to allow all possible adhesion parameters
    also between identical cells.'"""
    a = M.Adhesion(10)
    nodes = DEFAULTS.grn_nodes
    # locks all 0, keys all 1: every key complements every lock, including its own
    pattern = tuple([0] * 5 + [1] * 5 + [0] * (nodes - 10))
    assert a.between(pattern, pattern) == 31.0
    # hand-computed single bits: only key 0 set, so only 2^0 matches, both ways
    one = tuple([0] * 5 + [1, 0, 0, 0, 0] + [0] * (nodes - 10))
    assert a._match(one, one) == 1
    assert a.between(one, one) == 1.0
    # J_im is the binary place value of the medium nodes 3..7
    medium_bits = tuple(1 if i in (3, 5) else 0 for i in range(nodes))
    assert a.medium(medium_bits) == (1 << 0) + (1 << 2)


# ===========================================================================
# the Boolean gene-regulation network
# ===========================================================================

def test_genome_is_two_inputs_and_one_of_the_sixteen_boolean_functions():
    rng = np.random.default_rng(3)
    nodes = DEFAULTS.grn_nodes
    genome = M.Genome.random(nodes, DEFAULTS.signal_nodes, rng)
    assert len(genome.inputs) == len(genome.functions) == nodes
    for (a, b), f in zip(genome.inputs, genome.functions):
        assert 0 <= f < M.FUNCTIONS == 16
        for i in (a, b):
            assert -nodes <= i <= nodes - 1          # -nodes..-1 and 0..nodes-1
    # inputs are drawn uniformly from -nodes..-1 and 1..nodes, so half are
    # environmental: "a 0.5 probability of an input from another gene"
    draws = [i for _ in range(200)
             for row in M.Genome.random(nodes, 2, rng).inputs for i in row]
    assert np.mean([i >= 0 for i in draws]) == pytest.approx(0.5, abs=0.03)


def test_a_node_applies_its_truth_table_to_its_two_inputs():
    nodes = 6
    #: node 0 reads nodes 1 and 2; node 1 reads environment entry 0 and node 1
    inputs = [[1, 2], [-1, 1], [0, 0], [0, 0], [0, 0], [0, 0]]
    for table in range(16):
        genome = M.Genome([list(r) for r in inputs], [table] * nodes, nodes, 2)
        for a in (0, 1):
            for b in (0, 1):
                state = (0, a, b, 0, 0, 0)
                assert genome.step(state, (0,) * nodes)[0] == (table >> (a * 2 + b)) & 1
    # AND (bit 3 only) and XOR (bits 1 and 2) behave as they should
    genome = M.Genome([list(r) for r in inputs], [0b1000] * nodes, nodes, 2)
    assert [genome.step((0, a, b, 0, 0, 0), (0,) * nodes)[0]
            for a in (0, 1) for b in (0, 1)] == [0, 0, 0, 1]
    genome = M.Genome([list(r) for r in inputs], [0b0110] * nodes, nodes, 2)
    assert [genome.step((0, a, b, 0, 0, 0), (0,) * nodes)[0]
            for a in (0, 1) for b in (0, 1)] == [0, 1, 1, 0]


def test_a_mutation_changes_one_input_or_one_function():
    rng = np.random.default_rng(1)
    genome = M.Genome.random(DEFAULTS.grn_nodes, 2, rng)
    changes = []
    for _ in range(60):
        child = genome.mutate(rng)
        differing_inputs = sum(a != b for a, b in zip(genome.inputs, child.inputs))
        differing_functions = sum(a != b for a, b in zip(genome.functions, child.functions))
        assert differing_inputs + differing_functions <= 1     # a point mutation
        changes.append((differing_inputs, differing_functions))
    assert any(i for i, _ in changes) and any(f for _, f in changes)


def test_signalling_reads_the_neighbours_and_other_negative_inputs_read_zero():
    """Hogeweg: inputs -1 and -2 read nodes 1 and 2 of the neighbours (OR-ed);
    'the other negative numbers provide an invariant input value of 0'."""
    dev = M.Development(M.Genome.random(DEFAULTS.grn_nodes, 2, np.random.default_rng(0)),
                        params(divisions=1, steps_per_stage=2, final_steps=2),
                        np.random.default_rng(0), growth=False)
    dev._cleavage(0)                                  # two touching cells
    dev.state = {c: tuple(1 for _ in range(dev.genome.nodes)) for c in dev.cpm.volume}
    env = dev._environment()
    assert len(env) == 2
    for vector in env.values():
        assert vector[:2] == (1, 1)                   # the two signalling nodes are OR-ed in
        assert set(vector[2:]) == {0}                 # every other entry is a constant 0


def test_the_maternal_factor_breaks_the_symmetry_of_the_first_cleavage():
    """Hogeweg Table 1: 'If first division, flip state of node 20 of daughter cell'."""
    p = params(divisions=2, steps_per_stage=2, final_steps=2)
    dev = M.Development(M.Genome.random(p.grn_nodes, p.signal_nodes, np.random.default_rng(0)),
                        p, np.random.default_rng(0), growth=False)
    assert dev.maternal == (p.grn_nodes - 5, p.grn_nodes - 4) == (19, 20)
    (zygote,) = dev.cpm.volume
    assert dev.state[zygote] == (0,) * p.grn_nodes    # "initiated in state 0"

    dev._cleavage(0)
    states = [dev.state[c] for c in sorted(dev.cpm.volume)]
    assert len(states) == 2
    differing = [i for i, (a, b) in enumerate(zip(*states)) if a != b]
    assert differing == [dev.maternal[0]]             # exactly the maternal node


# ===========================================================================
# the Glazier & Graner cell-sorting experiment
# ===========================================================================

def test_default_energies_satisfy_the_published_engulfment_condition():
    """Hogeweg sec. 3.1: 'celltype A engulfs celltype B if J_ab < J_mb and J_am < J_bm'."""
    J = DEFAULTS.sorting_energies
    assert J["ab"] < J["bm"] and J["am"] < J["bm"]
    # and heterotypic contacts must cost more than the average homotypic one,
    # which is what makes a mixed aggregate sort at all
    assert J["ab"] > (J["aa"] + J["bb"]) / 2
    # B is the more cohesive type, so B is the one that ends up engulfed
    assert J["bb"] < J["aa"]
    assert all(0 <= v <= 31 for v in J.values())      # inside the five-bit range


@pytest.mark.parametrize("seed", [0, 3, 7])
def test_cell_sorting_reproduces_the_glazier_graner_result(seed):
    """A mixed aggregate of two cell types sorts out: the heterotypic contact
    fraction falls, the energy falls, and the cohesive type is engulfed."""
    net = generate_network(ID, seed=seed, mode="cell-sorting")
    a = net.extras["analysis"]
    initial = a["heterotypic_contact_fraction"]["initial"]
    final = a["heterotypic_contact_fraction"]["final"]
    assert initial > 0.35                              # a well-mixed start
    assert final < 0.85 * initial                      # ... that sorts out
    assert a["energy"]["final"] < a["energy"]["initial"]
    assert a["engulfed_type"] == "B"                   # the cohesive type goes inside

    trace = a["per_step"]
    assert len(trace) == DEFAULTS.sorting_steps + 1
    assert trace[0]["heterotypic_contact_fraction"] == initial
    # the sorting is a relaxation, so the late energy is well below the early one
    assert min(s["energy"] for s in trace[-10:]) < min(s["energy"] for s in trace[:10])
    assert net.status == "observed"
    assert net.initial_state == {"A": 12.0, "B": 12.0}
    assert "interaction_law" in net.extras


# ===========================================================================
# the evo-devo network
# ===========================================================================

def test_network_is_the_observed_events_of_one_development():
    net = default_net()
    assert net.status == "observed"
    ids = {s.id for s in net.species}
    nodes = net.params["grn_nodes"]
    for s in net.species:
        assert set(s.structure) <= {"0", "1"} and len(s.structure) == nodes
    zygote = next(s for s in net.species if s.id == "T0")
    assert zygote.structure == "0" * nodes             # "initiated in state 0"
    assert net.initial_state == {"T0": 1.0}

    kinds = {"differentiation": 0, "division": 0, "death": 0}
    for r in net.reactions:
        assert r.count >= 1
        (reactant,), (n,) = r.reactants.keys(), r.reactants.values()
        assert reactant in ids and n == 1              # one cell in, always
        products = sum(r.products.values())
        if products == 0:
            kinds["death"] += r.count
        elif products == 1:
            assert not r.catalysts                     # Ta -> Tb really differentiates
            kinds["differentiation"] += r.count
        else:
            assert r.products == {reactant: 2}         # Ta -> 2 Ta, both daughters
            kinds["division"] += r.count
    a = net.extras["analysis"]
    assert kinds == {"differentiation": a["differentiations"],
                     "division": a["divisions"], "death": a["deaths"]}
    # the pre-scheduled cleavages alone take one zygote to 2^divisions cells
    assert a["divisions"] >= 2 ** net.params["divisions"] - 1


def test_the_population_balance_of_the_events_gives_the_surviving_cells():
    """Every cell that lived was born by a division and left by a death."""
    net = default_net()
    a = net.extras["analysis"]
    born = 1 + a["divisions"]            # the zygote, plus one new cell per division
    assert born - a["deaths"] == a["cells_alive"] == len(net.extras["compartments"]["cells"])
    assert sum(net.extras["final_state"].values()) == a["cells_alive"]


def test_extras_carry_the_lattice_the_cells_and_the_energies():
    net = default_net()
    space, cells = net.extras["space"], net.extras["compartments"]["cells"]
    grid = net.params["grid"]
    assert space["shape"] == [grid, grid] and len(space["sigma"]) == grid * grid
    assert space["dimensions"] == 2

    occupied = [s for s in space["sigma"] if s != 0]
    assert len(occupied) == net.extras["analysis"]["volume"]
    assert sum(c["volume"] for c in cells.values()) == len(occupied)
    assert set(cells) == {str(s) for s in set(occupied)}
    ids = {s.id for s in net.species}
    for cell in cells.values():
        assert cell["type"] in ids
        assert cell["volume"] > 0
        assert 0 <= cell["centroid"][0] < grid and 0 <= cell["centroid"][1] < grid
    # the outermost ring of sites is never a copy target, so it stays medium
    for i in range(grid):
        assert space["sigma"][i] == space["sigma"][(grid - 1) * grid + i] == 0
        assert space["sigma"][i * grid] == space["sigma"][i * grid + grid - 1] == 0

    energies = net.extras["energies"]
    assert energies["temperature"] == net.params["temperature"]
    assert energies["inelasticity"] == net.params["inelasticity"]
    assert energies["delta_h_offset"] == M.DELTA_H_OFFSET
    assert energies["receptors"]["locks"] == [0, 1, 2, 3, 4]
    assert energies["receptors"]["keys"] == [5, 6, 7, 8, 9]
    assert all(0 <= v <= 31 for v in energies["J"].values())
    assert all(0 <= v <= 31 for v in energies["J_medium"].values())
    # the J matrix really is the one the adhesion rule gives for those types
    adhesion = M.Adhesion(net.params["adhesion_nodes"])
    structure = {s.id: tuple(int(b) for b in s.structure) for s in net.species}
    for key, value in energies["J"].items():
        a, b = key.split("-")
        assert adhesion.between(structure[a], structure[b]) == value
    for key, value in energies["J_medium"].items():
        assert adhesion.medium(structure[key]) == value


def test_the_grn_in_extras_is_the_network_that_was_evolved():
    net = default_net()
    grn = net.extras["grn"]
    nodes = net.params["grn_nodes"]
    assert grn["nodes"] == nodes and grn["signal_nodes"] == net.params["signal_nodes"]
    assert len(grn["inputs"]) == len(grn["functions"]) == nodes
    assert all(0 <= f < 16 for f in grn["functions"])
    assert all(-nodes <= i <= nodes - 1 for row in grn["inputs"] for i in row)
    assert net.extras["maternal_nodes"] == [nodes - 5, nodes - 4]


def test_reproducible_with_a_seed():
    # the contract test already pins this at the catalog defaults; here it is
    # cheap parameters, and both experiments
    small = dict(population=2, generations=1, divisions=2, steps_per_stage=4,
                 final_steps=8, sorting_steps=6)
    for mode in ("evo-devo", "cell-sorting"):
        one = generate_network(ID, seed=2, mode=mode, **small)
        two = generate_network(ID, seed=2, mode=mode, **small)
        other = generate_network(ID, seed=3, mode=mode, **small)
        assert one.to_dict() == two.to_dict()
        assert one.to_dict() != other.to_dict()


def test_parameters_are_checked():
    with pytest.raises(ValueError, match="even"):
        generate_network(ID, adhesion_nodes=9)
    with pytest.raises(ValueError, match="maternal"):
        generate_network(ID, adhesion_nodes=20)
    with pytest.raises(ValueError, match="signal_nodes"):
        generate_network(ID, signal_nodes=12)
    with pytest.raises(ValueError, match="missing"):
        generate_network(ID, mode="cell-sorting", sorting_energies={"aa": 1})
    with pytest.raises(ValueError, match=">= 0"):
        generate_network(ID, mode="cell-sorting",
                         sorting_energies={"aa": 1, "ab": 1, "bb": 1, "am": 1, "bm": -1})
    with pytest.raises(ValueError, match="too small"):
        generate_network(ID, mode="cell-sorting", grid=12, sorting_cells=40)
    with pytest.raises(ValueError, match="3 sites a side"):
        M.CPM(2, 5, temperature=1.0, inelasticity=1.0, bond=lambda a, b: 1.0)


# ===========================================================================
# the published evolutionary result
# ===========================================================================

@pytest.mark.slow
def test_evolution_selects_for_cell_differentiation():
    """Hogeweg's fitness is cell differentiation alone, and it climbs: more of
    the population differentiates, and the best differentiation found grows.
    'About one-third of the evolutionary runs lead to extensive cell
    differentiation and morphogenesis' - so not every run must succeed."""
    first, last, differentiating = [], [], []
    for seed in (1, 2, 3):
        net = generate_network(ID, seed=seed, population=6, generations=4)
        history = net.extras["analysis"]["fitness"]["per_generation"]
        assert len(history) == 5
        first.append(history[0]["best_fitness"])
        last.append(history[-1]["best_fitness"])
        differentiating.append((history[0]["differentiating"], history[-1]["differentiating"]))
        # the best ever found is at least as good as any single generation
        assert net.extras["analysis"]["fitness"]["best_evolved"] >= max(
            g["best_fitness"] for g in history)
    assert sum(last) > sum(first)                      # differentiation is selected for
    assert sum(b for _, b in differentiating) > sum(a for a, _ in differentiating)
    assert sum(f > 0 for f in last) >= 2               # most runs do differentiate


@pytest.mark.slow
def test_morphogenesis_is_a_side_effect_of_differentiation():
    """The phenomena of Hogeweg (2000) sec. 3.2: cells move, divide, die and
    redifferentiate, all driven by surface energy alone, and several cell types
    arise from one zygote through the maternal signals and neighbour induction."""
    seen = {"types": 0, "deaths": 0, "divisions": 0, "differentiations": 0}
    for seed in (1, 2, 3, 4):
        net = generate_network(ID, seed=seed, population=6, generations=3,
                               final_steps=60)
        a = net.extras["analysis"]
        assert a["divisions"] >= 2 ** net.params["divisions"] - 1
        seen["types"] = max(seen["types"], a["types_seen"])
        for key in ("deaths", "divisions", "differentiations"):
            seen[key] += a[key]
        # growth and death are never encoded in the genome: they can only come
        # from the forces that differential adhesion generates
        assert set(a) >= {"heterotypic_contact_fraction", "surface", "volume"}
        assert 0.0 <= a["heterotypic_contact_fraction"] <= 1.0
    assert seen["types"] > 2                    # one zygote gives many expression patterns
    assert seen["divisions"] > 0
    assert seen["differentiations"] > 0
    assert seen["deaths"] > 0                   # squeezing really does kill cells
