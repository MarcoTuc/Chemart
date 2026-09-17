"""RNA-folding ribozyme AC.

Published facts reproduced here:

- book table 18.1 (Banzhaf & Yamamoto 2015, section 18.1.1): two 60 nt sequences
  with the structures the Vienna RNAfold server returned for them;
- ViennaRNA regression values (the folding is deterministic, so these are exact);
- the sequence -> structure -> function map of Flamm et al. (2010), J. Syst. Chem.
  1:4: the longest loop is the imaginary transition structure, its size is the
  loop cycle length, and a mono-cyclic ITS of size 2n alternates n broken and n
  formed bonds (Fujita 1986), so only even cycles catalyse;
- the neutral-network phenomenon of Ullrich & Flamm: many point mutants keep the
  structure, while a variety of alternative functions sits one mutation away.
"""

import numpy as np
import pytest
import RNA

import chemart
from chemart.chemistries import rna_folding_ac as m

# --- book table 18.1 ---------------------------------------------------------
BOOK_SEQ_1 = "UGUGGCCUUCUGUAGAGUGAACUUCACCACAUAUGCUGUCUCUGGCACGUGGAUGGUUUG"
BOOK_STRUCT_1 = "...((((.((......((((...))))..(((.((((......)))).))))).)))).."
BOOK_SEQ_2 = "UCCGUCAGAAGCGCCAUUCAGGAUCACGUUACCCCGAAAAAAAGGUACCAGGAGCUCUUC"
BOOK_STRUCT_2 = ".......((((.((..(((.((.......((((..........)))))).))))).))))"


@pytest.fixture(scope="module")
def net():
    return chemart.generate_network("rna-folding-ac", seed=1)


def catalyst_of(reaction, functions):
    """The species that sits on both sides and carries a catalytic function."""
    both = [s for s in reaction.catalysts if s in functions]
    assert len(both) == 1, reaction.to_text()[:120]
    return both[0]


def substrates_and_products(reaction, catalyst):
    """The reaction without one copy of its catalyst on each side."""
    left, right = dict(reaction.reactants), dict(reaction.products)
    for side in (left, right):
        side[catalyst] -= 1
        if side[catalyst] == 0:
            del side[catalyst]
    return left, right


# --- the book's table ---------------------------------------------------------
def test_book_table_18_1_second_sequence():
    # table 18.1, right-hand molecule of figure 18.1: reproduced exactly.
    structure, energy = m.fold(BOOK_SEQ_2)
    assert structure == BOOK_STRUCT_2
    assert energy == pytest.approx(-7.10, abs=0.005)


def test_book_table_18_1_first_sequence_is_a_suboptimal_structure():
    # The book's first structure is not the MFE under current Turner parameters
    # (identical under Turner 1999 and 2004); it evaluates 1.1 kcal/mol above it
    # and is in the 1.5 kcal/mol suboptimal ensemble. Recorded in `decisions`.
    structure, energy = m.fold(BOOK_SEQ_1)
    assert structure != BOOK_STRUCT_1
    assert energy == pytest.approx(-10.90, abs=0.005)
    assert RNA.energy_of_struct(BOOK_SEQ_1, BOOK_STRUCT_1) == pytest.approx(-9.80, abs=0.005)
    ensemble = {s.structure for s in RNA.subopt(BOOK_SEQ_1, 150)}
    assert BOOK_STRUCT_1 in ensemble


def test_vienna_folding_regression():
    assert m.fold("GGGGAAAACCCC") == ("((((....))))", pytest.approx(-5.40, abs=0.005))


# --- the loop decomposition (ITS level 1) ------------------------------------
def test_loop_cycle_sizes():
    # cycle size = unpaired bases + 2 per delimiting stem
    assert m.loops("((((....))))") == [(4, 2, ()), (4, 2, ()), (4, 2, ()), (6, 1, (5, 6, 7, 8))]
    # a bulge/interior loop is delimited by two stems
    interior = m.loops("((..((....))..))")
    assert (8, 2, (3, 4, 13, 14)) in interior
    assert (6, 1, (7, 8, 9, 10)) in interior
    # a multiloop is delimited by three
    assert (12, 3, (3, 4, 11, 12, 19, 20)) in m.loops("((..((..))..((..))..))")


def test_loop_decomposition_is_consistent_on_the_book_structure():
    for structure in (BOOK_STRUCT_1, BOOK_STRUCT_2):
        cycles = m.loops(structure)
        assert cycles, structure
        assert all(size == len(unpaired) + 2 * stems for size, stems, unpaired in cycles)
        # one loop per base pair that closes one
        assert len(cycles) == structure.count("(")


# --- the structure -> function map (Flamm et al. 2010) -----------------------
def test_hairpin_loop_gives_a_cleavase():
    f = m.function("GGGGAAAACCCC", "((((....))))", 4, 12, 3)
    assert f == {
        "its_size": 6, "stems": 1, "loop_sequence": "AAAA",
        "site": "UUU", "reaction": "cleavage",
    }


def test_interior_loop_gives_a_ligase():
    f = m.function("GGCCGGAAAACCGGCC", "((..((....))..))", 4, 12, 3)
    assert f["its_size"] == 8 and f["stems"] == 2
    assert f["reaction"] == "ligation"
    assert f["loop_sequence"] == "CCGG"          # positions 3, 4, 13, 14
    assert f["site"] == "CGG"                    # reverse complement of "CCG"


def test_only_even_mono_cyclic_transition_structures_catalyse():
    # odd cycle: n bonds cannot alternate around the ring (Fujita 1986)
    assert m.function("GGGAAACCC", "(((...)))", 4, 12, 3) is None
    # multiloop: a composite transition structure
    assert m.function("GGAACCAACCAAGGAACCAAGG", "((..((..))..((..))..))", 4, 12, 3) is None
    # outside the accepted ITS size window
    assert m.function("GGGGAAAACCCC", "((((....))))", 8, 12, 3) is None


def test_catalytic_criterion_holds_over_random_sequences():
    rng = np.random.default_rng(9)
    odd = 0
    for _ in range(200):
        seq = "".join(rng.choice(list("ACGU"), 30))
        structure, _ = m.fold(seq)
        cycles = m.loops(structure)
        f = m.function(seq, structure, 4, 12, 3)
        if not cycles:
            assert f is None
            continue
        size, stems, _ = max(cycles)
        if size % 2:
            odd += 1
        if f is not None:
            assert size % 2 == 0 and stems <= 2 and 4 <= size <= 12
            assert f["reaction"] == ("cleavage" if stems == 1 else "ligation")
    assert odd > 0, "the sample should contain odd-sized transition structures"


# --- the network --------------------------------------------------------------
def test_default_closure_shape(net):
    assert (len(net.species), len(net.reactions)) == (21, 8)
    assert net.status == "complete"
    analysis = net.extras["analysis"]
    assert analysis["catalytic_species"] == 4
    assert (analysis["cleavases"], analysis["ligases"]) == (3, 1)


def test_every_reaction_is_catalysed_and_splits_or_joins(net):
    functions = net.extras["functions"]
    assert net.reactions
    for reaction in net.reactions:
        catalyst = catalyst_of(reaction, functions)
        site = functions[catalyst]["site"]
        left, right = substrates_and_products(reaction, catalyst)
        if functions[catalyst]["reaction"] == "cleavage":
            (substrate,) = [s for s, n in left.items() for _ in range(n)] or [None]
            fragments = [s for s, n in right.items() for _ in range(n)]
            assert substrate is not None and len(fragments) == 2
            assert "".join(sorted(fragments, key=substrate.find)) == substrate
            assert site in substrate
        else:
            pieces = [s for s, n in left.items() for _ in range(n)]
            (product,) = [s for s, n in right.items() for _ in range(n)]
            assert len(pieces) == 2
            assert "".join(pieces) == product or "".join(reversed(pieces)) == product


def test_nucleotides_are_conserved(net):
    (law,) = net.extras["conservation"]
    assert law["name"] == "nucleotides"
    weight = dict(zip([s.id for s in net.species], law["vector"]))
    assert all(weight[s.id] == len(s.id) for s in net.species)
    for reaction in net.reactions:
        left = sum(weight[s] * n for s, n in reaction.reactants.items())
        right = sum(weight[s] * n for s, n in reaction.products.items())
        assert left == right, reaction.to_text()[:120]


def test_species_carry_their_fold_and_energy(net):
    for species in net.species:
        sequence, structure = species.structure.split(" ")
        assert sequence == species.id
        assert (structure, net.extras["energies"][species.id]) == (
            m.fold(species.id)[0], pytest.approx(m.fold(species.id)[1], abs=0.005)
        )


def test_seed_is_reproducible_and_matters():
    again = chemart.generate_network("rna-folding-ac", seed=1)
    assert again.to_dict() == chemart.generate_network("rna-folding-ac", seed=1).to_dict()
    other = chemart.generate_network("rna-folding-ac", seed=2)
    assert [s.id for s in other.species] != [s.id for s in again.species]


def test_well_stirred_run_reports_the_reactions_that_fired():
    run = chemart.generate_network(
        "rna-folding-ac", seed=1, mode="well-stirred", pool=12, steps=400
    )
    assert run.status == "observed"
    assert run.reactions
    assert all(r.count >= 1 for r in run.reactions)
    functions = run.extras["functions"]
    weight = dict(zip([s.id for s in run.species], run.extras["conservation"][0]["vector"]))
    for reaction in run.reactions:
        # a collision always draws three molecules, so a cleavage recorded here
        # carries the third one through unchanged alongside the catalyst
        assert any(s in functions for s in reaction.catalysts), reaction.to_text()[:120]
        left = sum(weight[s] * n for s, n in reaction.reactants.items())
        assert left == sum(weight[s] * n for s, n in reaction.products.items())
    assert set(run.initial_state) <= {s.id for s in run.species}


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError, match="its_min"):
        chemart.generate_network("rna-folding-ac", seed=1, its_min=10, its_max=6)
    with pytest.raises(ValueError, match="pool"):
        chemart.generate_network("rna-folding-ac", seed=1, mode="well-stirred", pool=2)
    with pytest.raises(ValueError, match="min_recognition"):
        chemart.generate_network("rna-folding-ac", seed=1, seq_length=10, min_recognition=11)


# --- the neutral network (Ullrich & Flamm; book 18.1.1) ----------------------
def test_neutral_network_is_large_and_neighbours_carry_new_functions():
    """Book 18.1.1 / ref [874]: many neutral variants, and a variety of
    alternative functions within a short distance in sequence space."""
    rng = np.random.default_rng(5)
    neutral = total = 0
    functions_nearby = []
    for _ in range(20):
        seq = "".join(rng.choice(list("ACGU"), 30))
        reference, _ = m.fold(seq)
        variants = set()
        for i in range(30):
            for base in "ACGU":
                if base == seq[i]:
                    continue
                mutant = seq[:i] + base + seq[i + 1:]
                structure, _ = m.fold(mutant)
                total += 1
                neutral += structure == reference
                f = m.function(mutant, structure, 4, 12, 3)
                variants.add((f["reaction"], f["its_size"]) if f else None)
        functions_nearby.append(len(variants))
    assert neutral / total == pytest.approx(0.40, abs=0.10)
    assert min(functions_nearby) >= 3


@pytest.mark.slow
def test_larger_closure_truncates_and_stays_conservative():
    run = chemart.generate_network("rna-folding-ac", seed=1, pool=12, max_species=100)
    assert run.status == "truncated"
    assert len(run.species) == 100
    weight = {s.id: len(s.id) for s in run.species}
    for reaction in run.reactions:
        left = sum(weight[s] * n for s, n in reaction.reactants.items())
        assert left == sum(weight[s] * n for s, n in reaction.products.items())


@pytest.mark.slow
def test_paper_scale_sequences_fold_and_map():
    """Flamm et al. (2010) use tRNA-size genes of about 100 nt."""
    rng = np.random.default_rng(3)
    catalytic = 0
    for _ in range(40):
        seq = "".join(rng.choice(list("ACGU"), 100))
        structure, energy = m.fold(seq)
        assert len(structure) == 100 and energy < 0
        catalytic += m.function(seq, structure, 4, 12, 3) is not None
    assert catalytic > 0
