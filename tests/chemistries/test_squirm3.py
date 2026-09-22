"""Squirm3 reproduces Hutton (2002) table 1 = book table 11.5, its reaction
counts and its self-replication, and the enzyme encoding of Hutton (2007)
(book eq. 11.15) with the worked examples of Hutton (2004).
"""

from collections import Counter

import pytest

from chemart import evolve, generate_network
from chemart.chemistries import squirm3 as sq

SEED = "e8-a1-b1-f1"


def atom_counts(species_id: str) -> Counter:
    return sq._mol_from_id(species_id).atom_counts


def atoms_in(stoichiometry: dict) -> Counter:
    """Atoms of each type on one side of a reaction."""
    total: Counter = Counter()
    for species_id, n in stoichiometry.items():
        for t, c in atom_counts(species_id).items():
            total[t] += c * n
    return total


# --- the rule set ---------------------------------------------------------------
def test_table_11_5_is_huttons_table_1():
    # book table 11.5 / Hutton (2002) table 1, verbatim
    assert [r.to_text() for r in sq.parse_rules(sq.REPLICATOR_RULES)] == [
        "R1: e8 + e0 -> e4e3",      # specific association
        "R2: x4y1 -> x2y5",         # general transformation
        "R3: x5 + x0 -> x7x6",      # general homo-association
        "R4: x3 + y6 -> x2y3",      # general hetero-association
        "R5: x7y3 -> x4y3",         # general transformation
        "R6: f4f3 -> f8 + f8",      # specific dissociation
        "R7: x2y8 -> x9y1",         # general transformation
        "R8: x9y9 -> x8 + y8",      # general dissociation
    ]


def test_the_rules_expand_to_188_explicit_reactions():
    # Hutton (2002) 2.3: "there would be 188 of them"; with six types and ten
    # states there are 60 x 60 x 2 = 7200 possible inputs, 7012 of them null.
    rules = sq.parse_rules(sq.REPLICATOR_RULES)
    assert len(sq.explicit_rules(rules)) == 188
    assert (6 * 10) ** 2 * 2 == 7200
    assert 7200 - 188 == 7012


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_reactions_needed_for_replication_scale_as_5n2_21n_24(n):
    # Hutton (2002) prints 5n^2 + 21n + 4 for n bases, but expanding R1-R8 over
    # the n + 2 types gives 2 + (n + 2) + 5(n + 2)^2 = 5n^2 + 21n + 24, which is
    # the 188 the same paragraph states for n = 4. The constant is 24, not 4.
    rules = sq.parse_rules(sq.REPLICATOR_RULES)
    assert len(sq.explicit_rules(rules, sq.ATOM_TYPES[:n + 2])) == 5 * n * n + 21 * n + 24
    assert 200 * (n + 2) ** 2 == (10 * (n + 2)) ** 2 * 2      # possible inputs


def test_type_variables_are_independent_but_repeats_are_not():
    # x and y may take the same type (36 each); x5 + x0 needs one type (6)
    one = sq.parse_rules("x5 + x0 -> x7x6")
    two = sq.parse_rules("x4y1 -> x2y5")
    assert len(sq.explicit_rules(one)) == 6
    assert len(sq.explicit_rules(two)) == 36


def test_membrane_set_is_the_41_rules_of_2007():
    rules = sq.membrane_rules()
    names = list(dict.fromkeys(r.name for r in rules))
    assert names == [f"R{i}" for i in range(1, 42)]
    text = [r.to_text() for r in rules]
    assert "R1: e1a37 -> e5a10" in text            # table 1, first rule
    assert "R34: a34a33 -> a37 + a37" in text      # table 1, last rule
    assert any(r.kind == "enzyme" and r.name == "R40" for r in rules)
    assert any(r.kind == "readout" and r.name == "R37" for r in rules)


def test_bad_rule_text():
    with pytest.raises(ValueError, match="type never changes"):
        sq.parse_rules("a1 + b0 -> a1c0")
    with pytest.raises(ValueError, match="changes nothing"):
        sq.parse_rules("a1 + b0 -> a1 + b0")
    with pytest.raises(ValueError, match="exactly two atoms"):
        sq.parse_rules("a1 -> a2")


# --- the enzyme encoding (book eq. 11.15 = Hutton 2007 R40) ---------------------
def test_encoding_matches_the_2004_worked_examples():
    # Hutton (2004): "the enzyme d5731 codes for e0+a0 -> e2a3", with S = 18
    # states, T = 6 types, a = 0 .. f = 5, bonds 0 unbonded and 1 bonded.
    assert sq.encode_enzyme(0, 0, 2, 3, x=4, y=0, b1=0, b2=1, n_states=18) == 5731
    assert sq.decode_enzyme(5731, n_states=18) == {
        "g": 0, "h": 0, "j": 2, "k": 3, "x": 4, "y": 0, "b1": 0, "b2": 1}
    # the enzyme that grows the membrane: d6740227 causes a8 + a0 -> a8a7
    assert sq.encode_enzyme(8, 0, 8, 7, x=0, y=0, b1=0, b2=1, n_states=18) == 6740227
    assert sq.decode_enzyme(6740227, n_states=18)["j"] == 8


def test_encoding_round_trips_over_the_whole_range():
    for n_states in (10, 18, 38):
        for spec in ({"g": 0, "h": 0, "j": 0, "k": 1, "x": 0, "y": 0, "b1": 0, "b2": 0},
                     {"g": n_states - 1, "h": n_states - 2, "j": 3, "k": 4,
                      "x": 5, "y": 2, "b1": 1, "b2": 0},
                     {"g": 7, "h": 1, "j": 2, "k": 9, "x": 1, "y": 5, "b1": 1, "b2": 1}):
            i = sq.encode_enzyme(**spec, n_states=n_states)
            assert i >= n_states
            assert sq.decode_enzyme(i, n_states=n_states) == spec


def test_encoding_is_the_books_formula():
    # i = 2(2(|T|(|T|(|S|(|S|(|S| g + h) + j) + k) + x) + y) + b1) + b2 + |S|
    S, T = 38, 6
    g, h, j, k, x, y, b1, b2 = 3, 5, 7, 11, 2, 4, 1, 0
    expected = 2 * (2 * (T * (T * (S * (S * (S * g + h) + j) + k) + x) + y) + b1) + b2 + S
    assert sq.encode_enzyme(g, h, j, k, x, y, b1, b2, n_states=S, n_types=T) == expected


def test_gene_reads_out_in_base_4():
    # Hutton (2007) R37 and fig. 13: bca -> 62; ebcafbdcaf -> enzymes 62 and 158
    assert sq.gene_to_enzyme("bca", n_states=38) == 62 == 0b011000 + 38
    assert sq.gene_to_enzyme("b", n_states=38) == 39
    assert sq.gene_to_enzyme("bc", n_states=38) == 44
    assert sq.gene_to_enzyme("bdca", n_states=38) == 158
    # the book's example "ebdcaf" is the gene bdca = 1320 base 4 = 120, plus |S|
    assert sq.gene_to_enzyme("bdca", n_states=38) == int("1320", 4) + 38
    with pytest.raises(ValueError, match="gene base"):
        sq.gene_to_enzyme("bze", n_states=38)


# --- molecules ------------------------------------------------------------------
def test_canonical_ids_are_invariant_under_the_molecules_symmetries():
    forward = sq.canonical([("e", 8), ("a", 1), ("b", 1), ("f", 1)], [(0, 1), (1, 2), (2, 3)])
    backward = sq.canonical([("f", 1), ("b", 1), ("a", 1), ("e", 8)], [(0, 1), (1, 2), (2, 3)])
    assert forward.id == backward.id == "e8-a1-b1-f1"
    ring = [(0, 1), (1, 2), (2, 3), (3, 0)]
    a = sq.canonical([("a", 36), ("a", 36), ("a", 37), ("a", 36)], ring)
    b = sq.canonical([("a", 37), ("a", 36), ("a", 36), ("a", 36)], ring)
    assert a.id == b.id and a.id.endswith("/ring")
    assert sq.canonical([("e", 0)], []).id == "e0"
    assert sq._mol_from_id(a.id).id == a.id
    assert sq._mol_from_id("e8-a1-b1-f1").id == "e8-a1-b1-f1"


# --- self-replication (the book's headline result) ------------------------------
@pytest.mark.slow
@pytest.mark.parametrize("gene", ["e8-a1-f1", "e8-a1-b1-f1", "e8-c1-a1-d1-f1"])
def test_the_replicator_copies_an_arbitrary_sequence(gene):
    # Hutton (2002) eq. 4: e8 {x1}* f1 + {x0}* -> 2 e8 {x1}* f1. Any string of
    # a1..d1 with e8 at one end and f1 at the other replicates in a soup of
    # atoms in state 0 (book fig. 11.17 uses e8-a1-b1-f1).
    net = evolve("squirm3", seed=3, seed_molecule=gene).network
    assert net.status == "observed"
    splits = [r for r in net.reactions if r.products.get(gene) == 2]
    assert splits, f"{gene} never split into two copies of itself"
    # the two copies come from one double strand, and nothing else is consumed
    assert all(sum(r.reactants.values()) == 1 for r in splits)
    assert all(atom_counts(next(iter(r.reactants))) == atoms_in({gene: 2})
               for r in splits)
    assert net.initial_state[gene] == 1.0


@pytest.mark.slow
def test_replication_consumes_the_soup_of_atoms_in_state_0():
    net = evolve("squirm3", seed=3).network
    used = [r for r in net.reactions if any(s.endswith("0") and "-" not in s
                                            for s in r.reactants)]
    assert used, "no free atom in state 0 was ever built into a molecule"
    assert {"R1", "R3", "R4"} <= set(net.extras["rule_counts"])


# --- conservation ----------------------------------------------------------------
def test_every_reaction_conserves_every_atom_type():
    for net in (evolve("squirm3", seed=1).network, generate_network("squirm3", seed=1)):
        for r in net.reactions:
            assert atoms_in(r.reactants) == atoms_in(r.products), r.to_text()


def test_the_conservation_vectors_are_exact():
    # S^T m = 0 for the atom count of each type: the positive control
    net = evolve("squirm3", seed=1).network
    ids, R, P = net.matrices()
    stoichiometry = (P - R).toarray()
    laws = {law["name"]: law for law in net.extras["conservation"]}
    assert set(laws) >= {"atoms_a", "atoms_e", "atoms_f", "atoms"}
    for law in laws.values():
        assert set(law["vector"]) == set(ids)
        m = [law["vector"][i] for i in ids]
        assert list(m @ stoichiometry) == [0] * len(net.reactions), law["name"]
    assert laws["atoms"]["vector"]["e8-a1-b1-f1"] == 4
    assert laws["atoms_a"]["vector"]["e8-a1-b1-f1"] == 1


# --- the closure ------------------------------------------------------------------
@pytest.mark.slow
def test_closure_from_the_seed_molecule_finds_the_replication_path():
    net = generate_network("squirm3", seed=1, max_species=40)
    ids = {s.id for s in net.species}
    assert SEED in ids and {"e0", "a0", "b0", "f0"} <= ids
    # R1 joins a free e0 to the seed, the first step of fig. 11.17
    assert any(r.reactants == {SEED: 1, "e0": 1} for r in net.reactions)
    assert net.status in ("complete", "truncated")
    for r in net.reactions:
        assert atoms_in(r.reactants) == atoms_in(r.products), r.to_text()


def test_the_closure_starts_from_a_molecule_not_a_cell():
    with pytest.raises(ValueError, match="not a cell"):
        generate_network("squirm3", rules="membrane", n_states=38,
                         seed_molecule="cell:e1-b1-c1-a1-f1")


# --- the membrane chemistry --------------------------------------------------------
@pytest.mark.slow
def test_the_2007_cell_starts_dividing():
    # the starting cell of fig. 2: a membrane loop of a36 with two a37 anchors,
    # the gene string bonded to both. R1 (e1a37 -> e5a10) is the first rule to
    # fire, then R6 and R2 attach a second e atom to the membrane.
    net = evolve("squirm3", seed=5, rules="membrane", n_states=38,
                 seed_molecule="cell:e1-b1-c1-a1-f1", space="lattice-moore",
                 width=16, height=16, food=90, steps=400).network
    cell = max(net.initial_state, key=lambda s: sum(atom_counts(s).values()))
    labels = Counter(sq._mol_from_id(cell).labels)
    assert labels[("a", 36)] == 18 and labels[("a", 37)] == 2     # the membrane loop
    assert labels[("e", 1)] == labels[("f", 1)] == 1              # the gene's two ends
    assert sum(labels.values()) == 25                             # 20 membrane + 5 gene
    # Hutton (2007): "Reactions R1, R6, R2 and R3 are the first to occur", joining
    # a second e atom to the membrane; R35 lets the membrane gain and lose atoms.
    fired = net.extras["rule_counts"]
    assert {"R1", "R6", "R2", "R3"} <= set(fired), f"the division never started: {fired}"
    assert "R35" in fired
    for r in net.reactions:
        assert atoms_in(r.reactants) == atoms_in(r.products), r.to_text()


# --- the world and its parameters ---------------------------------------------------
def test_flooding_dissolves_molecules_back_to_single_atoms():
    net = evolve("squirm3", seed=4, steps=1200, flood_period=300).network
    assert net.extras["floods"] == 3
    assert net.extras["dissolved_by_flood"] >= 1
    assert net.extras["space"]["model"] == "lattice-2d"
    assert net.extras["space"]["reaction_neighbourhood"] == "von-neumann"


def test_cosmic_rays_start_the_chemistry_without_a_seed_in_an_excited_state():
    # Hutton (2002) experiment 3: randomising states lets molecules form
    quiet = evolve("squirm3", seed=6, seed_molecule="a0", steps=400)
    assert not quiet.network.reactions and set(quiet.series("molecules")) == {0}
    struck = evolve("squirm3", seed=6, seed_molecule="a0", steps=400, cosmic_ray=0.01)
    assert struck.network.reactions and max(struck.series("molecules")) > 0


def test_same_seed_same_network():
    a = evolve("squirm3", seed=11)
    b = evolve("squirm3", seed=11)
    assert a.to_dict() == b.to_dict()
    assert evolve("squirm3", seed=12).network.to_dict() != a.network.to_dict()


def test_frames_follow_the_world():
    traj = evolve("squirm3", seed=1)
    assert traj.clock == "steps"
    # the initial state, after the first step, every 3000 // 200 = 15 steps, and the end
    assert traj.times()[:4] == [0.0, 1.0, 16.0, 31.0] and traj.times()[-1] == 3000.0
    net = traj.network
    assert traj.frames[0].state == net.initial_state and not traj.frames[0].fired
    assert traj.frames[-1].state == {k: float(v) for k, v in net.extras["final_state"].items()}
    assert sum(n for f in traj.frames for _, _, n in f.fired) == sum(r.count for r in net.reactions)
    assert traj.series("molecules")[0] == 1


def test_continuous_space_runs_the_same_chemistry():
    net = evolve("squirm3", seed=1, space="continuous", width=200,
                 height=200, food=60, steps=200).network
    assert net.extras["space"]["model"] == "continuous-2d"
    assert net.extras["space"]["reaction_radius"] == 15.0
    assert net.extras["rule_counts"].get("R1")


def test_bad_parameters():
    with pytest.raises(ValueError, match="n_states must be greater"):
        generate_network("squirm3", rules="membrane", n_states=10)
    with pytest.raises(ValueError, match="rule_text"):
        generate_network("squirm3", rules="custom")
    with pytest.raises(ValueError, match="not an atom"):
        generate_network("squirm3", seed_molecule="e8-zz")
    with pytest.raises(ValueError, match="lattice holds"):
        evolve("squirm3", width=5, height=5, food=400)
    with pytest.raises(ValueError, match="evolve face"):
        generate_network("squirm3", steps=10)
    with pytest.raises(ValueError, match="n_types"):
        generate_network("squirm3", n_types=3, seed_molecule="e8-a1-f1")
