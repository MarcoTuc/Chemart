"""Adleman's DNA Hamiltonian-path computation (Science 266:1021-1024, 1994).

Published facts reproduced here:

- fig. 1 and its legend: the seven-vertex graph has exactly one Hamiltonian path
  for v_in = 0, v_out = 6, namely 0->1, 1->2, 2->3, 3->4, 4->5, 5->6; removing the
  edge 2->3, or moving the designated vertices to v_in = 3, v_out = 5, leaves
  none (p. 1021);
- fig. 2: the code words O2, O3, O4, the edge oligonucleotides O2->3 and O3->4
  and the splint Obar3, with the halving rule that builds them;
- step 3 (p. 1022): the gel band at 140 bp for paths of seven vertices, and the
  120 bp contamination that the paper saw in fig. 3B;
- step 4: successive affinity purification on Obar1..Obar5 removes the paths that
  miss a vertex, among them the double-passage path 0->3, 3->2, 2->3, 3->4, 4->5,
  5->6 that survives the gel;
- the three graduated-PCR band patterns worked out on p. 1022;
- note 4: 50 pmol of each oligonucleotide, about 3e13 copies.
"""

from itertools import permutations

import pytest

import chemart
from chemart.chemistries import dna_hpp as m

ID = "dna-hpp"
# Published in the text of p. 1022 as the three paths whose bands superimpose in fig. 3B.
HAMILTONIAN = [0, 1, 2, 3, 4, 5, 6]
MISSING_VERTEX_2 = [0, 1, 3, 4, 5, 6]
DOUBLE_VERTEX_3 = [0, 3, 2, 3, 4, 5, 6]


@pytest.fixture(scope="module")
def net():
    return chemart.generate_network(ID, seed=1)


def sequences(net):
    return {s.id: s.structure for s in net.species}


def hamiltonian_paths(edges, v_in, v_out, n):
    """Brute force, to check the chemistry against the definition of the problem."""
    present = {tuple(e) for e in edges}
    middle = [v for v in range(n) if v not in (v_in, v_out)]
    return [[v_in, *order, v_out] for order in permutations(middle)
            if all(e in present for e in zip([v_in, *order], [*order, v_out]))]


# --- the graph of figure 1 ---------------------------------------------------
def test_figure_1_has_a_unique_hamiltonian_path(net):
    graph = net.extras["graph"]
    assert (graph["vertices"], graph["v_in"], graph["v_out"]) == (7, 0, 6)
    assert len(graph["edges"]) == 14
    assert graph["edges"] == m.ADLEMAN_GRAPH
    # fig. 1 legend
    assert hamiltonian_paths(graph["edges"], 0, 6, 7) == [HAMILTONIAN]
    # "there are no edges entering vertex 0"  (p. 1021)
    assert not [e for e in graph["edges"] if e[1] == 0]


@pytest.mark.parametrize("given, kept_by_pcr", [
    (dict(graph=[e for e in m.ADLEMAN_GRAPH if e != [2, 3]]), 3),   # "if the edge 2->3 were removed"
    (dict(v_in=3, v_out=5), 9),                                     # "if the designated vertices were changed"
])
def test_a_graph_without_a_hamiltonian_path_yields_none(given, kept_by_pcr):
    """p. 1021: both variants of the instance have no Hamiltonian path."""
    run = chemart.generate_network(ID, seed=1, **given)
    analysis = run.extras["analysis"]
    graph = run.extras["graph"]
    assert hamiltonian_paths(graph["edges"], graph["v_in"], graph["v_out"], 7) == []
    assert analysis["hamiltonian_paths"] == []
    assert analysis["solution_species"] == []
    # the earlier steps still find paths; the selection is what empties
    assert analysis["steps"][1]["kept"] == kept_by_pcr
    assert analysis["steps"][4]["kept"] == 0


# --- the encoding of figure 2 ------------------------------------------------
def test_figure_2_sequences_are_reproduced(net):
    seq = sequences(net)
    assert seq["O2"] == "TATCGGATCGGTATATCCGA"
    assert seq["O3"] == "GCTATTCGAGCTTAAAGCTA"
    assert seq["O4"] == "GGCTAGGTACCAGCATGCTT"
    # the two edge oligonucleotides fig. 2 prints
    assert seq["E2-3"] == "GTATATCCGAGCTATTCGAG" == seq["O2"][10:] + seq["O3"][:10]
    assert seq["E3-4"] == "CTTAAAGCTAGGCTAGGTAC" == seq["O3"][10:] + seq["O4"][:10]
    # "All oligonucleotides written 5' to 3' except Obar3": the splint is printed
    # base by base under O3, i.e. 3' -> 5'.
    assert m.complement(seq["O3"]) == "CGATAAGCTCGAATTTCGAT"
    assert seq["Obar3"] == m.revcomp(seq["O3"]) == "CGATAAGCTCGAATTTCGAT"[::-1]


def test_base_pairing_is_exact(net):
    seq = sequences(net)
    for i in range(7):
        code, splint = seq[f"O{i}"], seq[f"Obar{i}"]
        assert len(code) == len(splint) == 20
        assert set(code) <= set("ACGT")
        assert m.revcomp(splint) == code                      # involution
        assert all(m.COMPLEMENT[a] == b for a, b in zip(code, reversed(splint)))


def test_edge_oligos_join_the_correct_halves(net):
    """O_{i->j} = 3' 10-mer of O_i + 5' 10-mer of O_j, whole code words at the ends."""
    seq = sequences(net)
    for i, j in net.extras["graph"]["edges"]:
        left = seq[f"O{i}"] if i == 0 else seq[f"O{i}"][10:]
        right = seq[f"O{j}"] if j == 6 else seq[f"O{j}"][:10]
        assert seq[f"E{i}-{j}"] == left + right
        assert len(seq[f"E{i}-{j}"]) == 20 + 10 * ((i == 0) + (j == 6))
    # orientation is preserved: O_{2->3} is not O_{3->2}
    assert seq["E2-3"] != seq["E3-2"]


def test_a_splint_bridges_the_junction_it_ligates(net):
    """The splint of vertex j pairs with the 10-mer on each side of the junction."""
    seq = sequences(net)
    strand = seq["P0-1-2-3-4-5-6"]
    for slot, vertex in enumerate(HAMILTONIAN[1:-1], start=1):
        junction = strand[20 * slot: 20 * slot + 20]
        assert junction == seq[f"O{vertex}"]                  # the code word is restored
        assert m.revcomp(junction) == seq[f"Obar{vertex}"]     # and the splint covers it


# --- the ligation network (step 1) -------------------------------------------
def test_default_network_shape(net):
    assert (len(net.species), len(net.reactions)) == (797, 769)
    assert net.status == "truncated"          # a cyclic graph has infinitely many walks
    assert all(r.rate is None for r in net.reactions)
    assert net.extras["analysis"]["truncated_at_vertices"] == 7


def test_every_reaction_is_path_plus_splint_plus_edge(net):
    seq = sequences(net)
    for r in net.reactions:
        reactants = [s for s, n in r.reactants.items() for _ in range(n)]
        (product,) = [s for s, n in r.products.items() for _ in range(n)]
        assert len(reactants) == 3
        (splint,) = [s for s in reactants if s.startswith("Obar")]
        growing, edge = [s for s in reactants if s != splint]
        assert edge.startswith("E")
        junction = int(splint[4:])
        # the edge leaves the vertex the splint names, and the path ends there
        assert edge.startswith(f"E{junction}-")
        assert seq[growing].endswith(seq[f"O{junction}"][:10])
        assert seq[product] == seq[growing] + seq[edge]
        assert product.startswith("P")


def test_ligation_conserves_nucleotides(net):
    (law,) = net.extras["conservation"]
    assert law["name"] == "nucleotides"
    weight = dict(zip([s.id for s in net.species], law["vector"]))
    # a path duplex is its sense strand plus one splint per interior vertex
    assert weight["P0-1-2-3-4-5-6"] == 140 + 5 * 20
    for r in net.reactions:
        left = sum(weight[s] * n for s, n in r.reactants.items())
        assert left == sum(weight[s] * n for s, n in r.products.items()), r.to_text()[:120]


def test_initial_state_is_the_published_ligation_mix(net):
    """Note 4: 50 pmol of each edge oligonucleotide and of each splint except 0 and 6."""
    assert set(net.initial_state) == (
        {f"E{i}-{j}" for i, j in net.extras["graph"]["edges"]} | {f"Obar{i}" for i in range(1, 6)}
    )
    assert set(net.initial_state.values()) == {50.0}
    copies = net.extras["analysis"]["experiment"]["copies_per_oligonucleotide"]
    assert copies == pytest.approx(3e13, rel=0.05)     # "approximately 3 x 10^13 copies"


# --- the five steps ----------------------------------------------------------
def test_the_five_steps_filter_as_published(net):
    steps = net.extras["analysis"]["steps"]
    assert [s["step"] for s in steps] == [1, 2, 3, 4, 5]
    assert steps[1]["detail"] == "primers O0 and Obar6"
    assert steps[2]["band_bp"] == 140 == 20 * 7           # 20-mers, seven vertices
    assert steps[3]["order"] == [1, 2, 3, 4, 5]
    # each step keeps a subset of the previous one, and the last leaves one molecule
    kept = [steps[1]["kept"], steps[2]["kept"], steps[3]["kept"], steps[4]["kept"]]
    assert kept == sorted(kept, reverse=True) and kept[-1] == 1


def test_length_selection_keeps_the_140_mers(net):
    seq = sequences(net)
    analysis = net.extras["analysis"]
    survivors = analysis["gel_survivors"]
    # exactly the paths from 0 to 6 with seven vertex slots
    assert all(len(seq[m.molecule_id(w)]) == 140 for w in survivors)
    assert all(w[0] == 0 and w[-1] == 6 and len(w) == 7 for w in survivors)
    assert HAMILTONIAN in survivors
    # the paper's double-passage path is 140 bp too and survives the gel
    assert DOUBLE_VERTEX_3 in survivors
    # while the path that misses vertex 2 is 120 bp: the fig. 3B contamination
    missing = "".join(seq[f"E{i}-{j}"] for i, j in zip(MISSING_VERTEX_2, MISSING_VERTEX_2[1:]))
    assert len(missing) == 120
    assert MISSING_VERTEX_2 not in survivors


def test_affinity_purification_removes_paths_missing_a_vertex(net):
    analysis = net.extras["analysis"]
    step4 = analysis["steps"][3]
    # the first round (beads carrying Obar1) already removes 0->3, 3->2, 2->3, ...
    assert DOUBLE_VERTEX_3 in analysis["gel_survivors"]
    assert 1 not in DOUBLE_VERTEX_3
    assert step4["kept_after"] == [1, 1, 1, 1, 1]
    assert analysis["hamiltonian_paths"] == [HAMILTONIAN]


def test_final_readout_is_the_published_hamiltonian_path(net):
    analysis = net.extras["analysis"]
    assert analysis["solution_species"] == ["P0-1-2-3-4-5-6"]
    seq = sequences(net)
    assert len(seq["P0-1-2-3-4-5-6"]) == 140
    # the solution molecule is the ligation of the six edge oligonucleotides
    assert seq["P0-1-2-3-4-5-6"] == "".join(
        seq[f"E{i}-{j}"] for i, j in zip(HAMILTONIAN, HAMILTONIAN[1:])
    )


def test_graduated_pcr_reproduces_the_published_bands(net):
    """p. 1022 works out the bands of three paths; an empty lane is the paper's x."""
    assert m.graduated_pcr(HAMILTONIAN, 20, 7) == {
        "1": [40], "2": [60], "3": [80], "4": [100], "5": [120], "6": [140]}
    assert m.graduated_pcr(MISSING_VERTEX_2, 20, 7) == {
        "1": [40], "2": [], "3": [60], "4": [80], "5": [100], "6": [120]}
    # "80bp/40bp denotes that both a 40bp and an 80bp band will be produced in lane 3"
    assert m.graduated_pcr(DOUBLE_VERTEX_3, 20, 7) == {
        "1": [], "2": [60], "3": [40, 80], "4": [100], "5": [120], "6": [140]}
    assert net.extras["analysis"]["graduated_pcr"]["P0-1-2-3-4-5-6"]["6"] == [140]


# --- parameters --------------------------------------------------------------
def test_seed_reproducibility_and_what_the_seed_changes():
    again = chemart.generate_network(ID, seed=1)
    assert again.to_dict() == chemart.generate_network(ID, seed=1).to_dict()
    other = chemart.generate_network(ID, seed=2)
    one, two = sequences(again), sequences(other)
    # only the code words the paper does not publish are drawn from rng
    assert one["O0"] != two["O0"] and one["O1"] != two["O1"]
    assert one["O2"] == two["O2"] and one["O3"] == two["O3"] and one["O4"] == two["O4"]
    # the computation does not depend on which random code words were drawn
    assert other.extras["analysis"]["hamiltonian_paths"] == [HAMILTONIAN]


def test_an_acyclic_graph_gives_a_complete_closure():
    run = chemart.generate_network(
        ID, seed=1, graph=[[0, 1], [1, 2], [0, 2]], vertices=3, v_in=0, v_out=2,
        max_path_vertices=3, sequences={},
    )
    assert run.status == "complete"
    assert len(run.reactions) == 1                      # E0-1 + Obar1 + E1-2 -> P0-1-2
    assert run.extras["analysis"]["hamiltonian_paths"] == [[0, 1, 2]]
    assert len(sequences(run)["P0-1-2"]) == 60 == 20 * 3


def test_a_longer_closure_keeps_the_same_answer():
    """Longer walks are made and thrown away by the gel, as in the tube."""
    run = chemart.generate_network(ID, seed=1, max_path_vertices=9)
    assert len(run.species) > 3000 and run.status == "truncated"
    assert run.extras["analysis"]["hamiltonian_paths"] == [HAMILTONIAN]
    assert len(run.extras["analysis"]["gel_survivors"]) == 2


@pytest.mark.parametrize("given, message", [
    (dict(oligo_length=21), "oligo_length must be even"),
    (dict(v_in=6), "must be different"),
    (dict(graph=[[0, 9]]), "not in 0..6"),
    (dict(graph=[[0, 0]]), "self-loop"),
    (dict(sequences={"2": "ACGT"}), "sequences"),
    (dict(sequences={"x": "A" * 20}), "sequences keys"),
    (dict(sequences={"2": "ACGU" * 5}), "over ACGT"),
])
def test_rejects_inconsistent_parameters(given, message):
    with pytest.raises(ValueError, match=message):
        chemart.generate_network(ID, seed=1, **given)


def test_the_selection_finds_exactly_the_hamiltonian_paths_of_random_graphs():
    """The five steps are a decision procedure: check them against brute force."""
    import numpy as np

    rng = np.random.default_rng(4)
    found = 0
    for _ in range(25):
        n = int(rng.integers(4, 7))
        edges = [[i, j] for i in range(n) for j in range(n)
                 if i != j and j != 0 and i != n - 1 and rng.random() < 0.55]
        if not edges:
            continue
        run = chemart.generate_network(
            ID, seed=2, graph=edges, vertices=n, v_in=0, v_out=n - 1,
            max_path_vertices=n, sequences={},
        )
        expected = hamiltonian_paths(edges, 0, n - 1, n)
        assert sorted(run.extras["analysis"]["hamiltonian_paths"]) == sorted(expected)
        found += bool(expected)
    assert found, "the sample should contain graphs with a Hamiltonian path"
