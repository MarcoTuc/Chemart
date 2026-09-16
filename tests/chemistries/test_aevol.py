"""Aevol (book 18.1.2; Knibbe et al. 2007; Parsons, Knibbe & Beslon 2011, book [654]).

The model is a port of the aevol source (version 9.4.0, commit bea4a25,
https://gitlab.inria.fr/aevol/aevol), BASE_2 flavour. The reference values
below are those of aevol's own unit tests: the hand-built chromosome of
test/gtest/DnaTest.cpp and test/gtest/EukTest_withrnas.cpp (promoters,
transcription levels, RNA lengths, gene positions and sizes) and the triangle
and area values of test/gtest/DiscreteDoubleFuzzyTest.cpp.
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import aevol as A

ID = "aevol"

# --- the hand-built chromosome of aevol's DnaTest.cpp ------------------------
AS = ["0011", "11101", "110011", "11000", "000101"]          # arbitrary sequences
#: Shine-Dalgarno, 4-base spacer, START, M0, H0, H0, W0, STOP
GENE = A.SHINE_DAL_SEQ + "0011" + "000" + "100" + "110" + "110" + "010" + "001"
TERM = "01000001101"                                          # stem-loop terminator
#: 2 and 1 mismatches from the consensus => transcription levels 0.6 and 0.8
PROM = ["0101010001110110010110", "0101011001110010010010"]
FIXTURE = AS[0] + PROM[0] + AS[1] + GENE + AS[2] + TERM + AS[3] + PROM[1] + AS[4]

#: the same gene and promoter on the lagging strand (EukTest_withrnas.cpp)
LAG_GENE = "0111011001001101110011" + "001001"
LAG_PROM = "1001011001000111010101"


def decode(genome, w_max=0.1):
    g = A.bits(genome)
    rnas = A.transcribe(g)
    return g, rnas, A.translate(g, rnas, w_max)


# --- transcription ----------------------------------------------------------
def test_promoters_and_transcription_levels_of_the_aevol_fixture():
    assert len(FIXTURE) == 109 and len(GENE) == 28 and len(TERM) == A.TERM_SIZE
    g = A.bits(FIXTURE)
    lead, lag = A.promoter_distances(g)
    # DnaTest.cpp: promoters at 4 (2 mismatches) and 81 (1 mismatch), none lagging
    assert [(int(p), int(lead[p])) for p in np.flatnonzero(lead <= A.PROM_MAX_DIFF)] == [(4, 2), (81, 1)]
    assert np.flatnonzero(lag <= A.PROM_MAX_DIFF).size == 0
    assert A.PROM_SEQ == "0101011001110010010110" and A.PROM_MAX_DIFF == 4
    # basal level 1 - d/(PROM_MAX_DIFF + 1): 1.0, 0.8, 0.6, 0.4, 0.2
    _, rnas, _ = decode(FIXTURE)
    assert [(r.pos, r.strand, r.basal) for r in rnas] == [(4, A.LEADING, 0.6), (81, A.LEADING, 0.8)]


def test_terminator_is_a_stem_loop_of_four_complementary_pairs():
    g = A.bits(TERM + "0" * 40)
    lead, lag = A.terminators(g)
    assert lead[0]
    for t in range(A.TERM_STEM_SIZE):                       # the four pairs really are complementary
        assert g[t] != g[A.TERM_SIZE - 1 - t]
    # a lagging terminator at q is a leading one at q - (TERM_SIZE - 1)
    assert lag[A.TERM_SIZE - 1]
    assert not A.terminators(A.bits("0" * 40))[0].any()      # a run of zeros pairs with itself


def test_transcripts_run_from_the_promoter_to_the_first_terminator():
    _, rnas, _ = decode(FIXTURE)
    # the terminator sits at 65; transcription starts at 4 + 22 = 26 and at 81 + 22 = 103 (wrapping)
    assert [r.length for r in rnas] == [50, 82]
    assert [r.start for r in rnas] == [26, 103]
    assert rnas[0].length == (65 - 26) + A.TERM_SIZE

    # EukTest_withrnas.cpp indiv1 and indiv2 (their terminator is found before the end,
    # so the linear chromosome of that test and our circular one agree)
    _, rnas, prots = decode(PROM[0] + AS[1] + GENE + AS[2] + TERM)
    assert len(rnas) == 1 and (rnas[0].pos, rnas[0].basal, rnas[0].length) == (0, 0.6, 50)
    assert [(p.first_aa, len(p.codons), p.e) for p in prots] == [(27 + 13, 4, 0.6)]
    _, rnas, prots = decode(PROM[0] + AS[2] + TERM)
    assert len(rnas) == 1 and rnas[0].length == 17 and prots == []


def test_lagging_strand_is_read_backwards():
    # EukTest_withrnas.cpp indiv3: promoter at 71, RNA of 47 bases, gene at 43 - 13
    genome = TERM + AS[1] + LAG_GENE + AS[4] + LAG_PROM
    _, rnas, prots = decode(genome)
    assert len(genome) == 72
    assert [(r.pos, r.strand, r.basal, r.length) for r in rnas] == [(71, A.LAGGING, 0.6, 47)]
    assert [(p.first_aa, p.strand, len(p.codons), p.e) for p in prots] == [(43 - 13, A.LAGGING, 4, 0.6)]


# --- translation and folding ------------------------------------------------
def test_gene_starts_after_the_rbs_and_two_transcripts_pool_their_levels():
    _, rnas, prots = decode(FIXTURE)
    # DnaTest.cpp: two proteins, both at position_first_aa 44, the second a duplicate
    # of the first; here the duplicate is merged and its transcription level added
    assert len(prots) == 1
    prot = prots[0]
    assert prot.first_aa == 31 + A.RBS_SIZE == 44 and A.RBS_SIZE == 13
    assert len(prot.codons) == 4
    assert prot.codons == [A.CODON_M0, A.CODON_H0, A.CODON_H0, A.CODON_W0]
    assert prot.e == pytest.approx(0.6 + 0.8)                # both transcripts carry the gene


def test_codons_fold_into_the_triangle_m_w_h():
    w_max = 0.033333333
    # M0 H0 H0 W0: every digit 0, so m = 0, h = -1 (rescaled from 0 in [-1, 1]) and w = 0;
    # w = 0 makes the protein non-functional
    assert A.fold([A.CODON_M0, A.CODON_H0, A.CODON_H0, A.CODON_W0], w_max) == (0.0, 0.0, -1.0, False)
    # all ones: m = 1, w = w_max, h = +1
    assert A.fold([A.CODON_M1, A.CODON_W1, A.CODON_H1], w_max) == (1.0, w_max, 1.0, True)
    # Gray code: M1 M0 reads 11 -> binary 10 = 2, over 2 codons -> 2/3
    m, _, _, _ = A.fold([A.CODON_M1, A.CODON_M0, A.CODON_W1, A.CODON_H1], w_max)
    assert m == pytest.approx(1.0)                            # gray 11 -> 11 -> 3/3
    m, _, _, _ = A.fold([A.CODON_M1, A.CODON_M1, A.CODON_W1, A.CODON_H1], w_max)
    assert m == pytest.approx(2 / 3)                          # gray 11 -> 10 -> 2/3
    # the START codon codes for the same amino acid as H0
    assert A.fold([A.CODON_START, A.CODON_M1, A.CODON_W1], w_max) == (1.0, w_max, -1.0, True)
    # a missing kind takes its default: m = 0.5, w = 0, h = 0.5 -> h rescaled to 0
    assert A.fold([A.CODON_W1], w_max) == (0.5, w_max, 0.0, False)
    assert A.fold([], w_max) == (0.5, 0.0, 0.0, False)


def test_a_gene_needs_a_stop_codon_inside_its_transcript():
    # the fixture gene without its STOP codon yields no protein
    truncated = AS[0] + PROM[0] + AS[1] + GENE[:-3] + "110" + AS[2] + TERM + AS[3] + PROM[1] + AS[4]
    _, _, prots = decode(truncated)
    assert prots == []


# --- the fuzzy sets (DiscreteDoubleFuzzyTest.cpp) ---------------------------
def test_add_triangle_matches_the_aevol_unit_test():
    # 5 points => x = 0, .25, .5, .75, 1
    points = np.zeros(5)
    A.add_triangle(points, 0.5, 0.25, 0.1)
    assert points.tolist() == [0, 0, 0.1, 0, 0]

    m, w, h = 0.5, 0.3, 0.1
    points = np.zeros(5)
    A.add_triangle(points, m, w, h)
    assert points == pytest.approx([0, h / w * (0.25 - (m - w)), h, h / w * ((m + w) - 0.75), 0])

    m, w, h = 0.6, 0.3, 0.25
    points = np.zeros(5)
    A.add_triangle(points, m, w, h)
    assert points == pytest.approx([0, 0, h / w * (0.5 - (m - w)), h / w * ((m + w) - 0.75), 0])


def test_geometric_area_of_a_flat_fuzzy_set():
    # DiscreteDoubleFuzzyTest Flat: 300 points at 0.5 have area 0.5
    assert A.geometric_area(np.full(300, 0.5)) == pytest.approx(0.5)
    assert A.geometric_area(np.zeros(300)) == 0.0
    # a triangle of half-width w and height h has area w * h
    points = np.zeros(10001)
    A.add_triangle(points, 0.5, 0.2, 0.7)
    assert A.geometric_area(points) == pytest.approx(0.2 * 0.7, rel=1e-3)


def test_phenotype_sums_the_protein_triangles_and_clips_them():
    class P:
        def __init__(self, m, w, h, e, functional=True):
            self.m, self.w, self.h, self.e, self.functional = m, w, h, e, functional

    sampling = 301
    activator, inhibitor = P(0.5, 0.1, 0.8, 1.0), P(0.5, 0.05, -0.5, 1.0)
    expected = np.zeros(sampling)
    A.add_triangle(expected, 0.5, 0.1, 0.8)
    minus = np.zeros(sampling)
    A.add_triangle(minus, 0.5, 0.05, -0.5)
    assert A.phenotype([activator, inhibitor], sampling) == pytest.approx(np.maximum(expected + minus, 0.0))

    # heights are scaled by the transcription level, and clipped at 1
    tall = A.phenotype([P(0.5, 0.1, 1.0, 1.0), P(0.5, 0.1, 1.0, 1.0)], sampling)
    assert tall.max() == pytest.approx(1.0)
    # an inhibitor alone leaves a flat zero phenotype (the sum is clipped below at 0)
    assert A.phenotype([inhibitor], sampling).max() == 0.0
    # a non-functional protein contributes nothing
    assert A.phenotype([P(0.5, 0.1, 0.8, 1.0, functional=False)], sampling).max() == 0.0


def test_target_is_the_sum_of_gaussians_clipped_to_the_unit_square():
    target = A.target_function(A.DEFAULT_TARGET, 300)
    assert target.size == 300
    assert target.min() >= 0.0 and target.max() <= 1.0
    x = np.arange(300) / 299
    raw = sum(h * np.exp(-(x - m) ** 2 / (2 * w * w)) for h, m, w in A.DEFAULT_TARGET)
    assert target == pytest.approx(np.clip(raw, 0.0, 1.0))
    # the -1.4 gaussian at 0.5 carves a hole out of the 1.2 gaussian at 0.52
    assert target[int(np.argmin(np.abs(x - 0.5)))] == 0.0
    assert A.geometric_area(target) == pytest.approx(0.152482671713, rel=1e-9)


def test_a_genome_without_a_gene_has_a_flat_phenotype_and_zero_fitness():
    world = A.World(np.random.default_rng(0), width=1, height=1)
    indiv = world.evaluate(np.zeros(1000, dtype=np.uint8))
    assert indiv.proteins == []
    assert indiv.phenotype.max() == 0.0
    assert indiv.metabolic_error == pytest.approx(world.target_area)
    assert indiv.fitness == 0.0                               # a flat phenotype is sterile


def test_fitness_is_an_exponential_of_the_metabolic_error():
    world = A.World(np.random.default_rng(0), width=1, height=1, selection_pressure=10.0)
    indiv = world.random_individual(5000)
    assert indiv.metabolic_error < world.target_area          # make_random beats a flat phenotype
    assert indiv.fitness == pytest.approx(np.exp(-10.0 * indiv.metabolic_error))
    assert indiv.phenotype.size == world.target.size == 300


# --- mutation operators -----------------------------------------------------
def genome_of(n, seed=0):
    return np.random.default_rng(seed).integers(0, 2, n, dtype=np.uint8)


def test_local_mutations():
    g = genome_of(40)
    switched = A.do_switch(g.copy(), 7)
    assert switched[7] != g[7] and (np.delete(switched, 7) == np.delete(g, 7)).all()

    seq = np.array([1, 0, 1], dtype=np.uint8)
    inserted = A.do_small_insertion(g, 10, seq)
    assert inserted.size == 43 and (inserted[10:13] == seq).all()
    assert (np.concatenate([inserted[:10], inserted[13:]]) == g).all()

    assert (A.do_small_deletion(g, 10, 4) == np.concatenate([g[:10], g[14:]])).all()
    wrapped = A.do_small_deletion(g, 38, 5)                   # spans the origin
    assert wrapped.size == 35 and (wrapped == g[3:38]).all()


def test_rearrangements_preserve_or_move_whole_segments():
    g = genome_of(40, seed=1)
    dup = A.do_duplication(g, 5, 15, 30)
    assert dup.size == 50 and (dup[30:40] == g[5:15]).all()
    assert (np.concatenate([dup[:30], dup[40:]]) == g).all()
    wrap = A.do_duplication(g, 35, 5, 0)                      # the segment spans the origin
    assert wrap.size == 50 and (wrap[:10] == np.concatenate([g[35:], g[:5]])).all()

    assert (A.do_deletion(g, 5, 15) == np.concatenate([g[:5], g[15:]])).all()
    assert (A.do_deletion(g, 35, 5) == g[5:35]).all()

    inv = A.do_inversion(g, 5, 15)
    assert inv.size == 40 and (inv[5:15] == 1 - g[5:15][::-1]).all()
    assert (A.do_inversion(inv, 5, 15) == g).all()            # an inversion is an involution

    for invert in (False, True):
        trans = A.do_translocation(g, 5, 20, 12, 30, invert)
        assert trans.size == 40
        assert Counter(trans.tolist()) == Counter(g.tolist()) or invert
        assert (trans[:5] == g[:5]).all() and (trans[30:] == g[30:]).all()


def test_mutations_keep_the_genome_within_its_limits():
    rng = np.random.default_rng(3)
    rates = dict.fromkeys(A.MUTATION_TYPES, 2e-3)
    mutator = A.Mutator(rng, rates, 6, (50, 400))
    g = genome_of(200, seed=2)
    drawn = 0
    for _ in range(300):
        events = mutator.draw(g.size)
        drawn += len(events)
        g = A.apply_mutations(g, events)
        assert 50 <= g.size <= 400
        assert g.dtype == np.uint8 and set(np.unique(g).tolist()) <= {0, 1}
    assert drawn > 100                                        # the operators really did fire
    assert set(A.MUTATION_TYPES) == {
        "switch", "small_insertion", "small_deletion",
        "duplication", "deletion", "translocation", "inversion"}


def test_mutation_counts_are_binomial_in_the_genome_length():
    rng = np.random.default_rng(5)
    rates = {k: 0.0 for k in A.MUTATION_TYPES}
    rates["switch"] = 1e-3
    mutator = A.Mutator(rng, rates, 6, (1, 10 ** 6))
    counts = [len(mutator.draw(1000)) for _ in range(400)]
    assert np.mean(counts) == pytest.approx(1.0, rel=0.25)    # 1000 x 1e-3 switches
    assert A.Mutator(rng, {k: 0.0 for k in A.MUTATION_TYPES}, 6, (1, 10 ** 6)).draw(1000) == []


# --- the observed network ---------------------------------------------------
SMALL = dict(grid_width=3, grid_height=3, generations=8, max_genome_length=10000)


def genotypes(net):
    return [s for s in net.species if s.id.startswith("G")]


def is_expression(reaction):
    return any(s.startswith("P") for s in reaction.products)


def test_network_is_a_population_run_of_genotypes_and_their_proteomes():
    net = generate_network(ID, seed=3, **SMALL)
    assert net.status == "observed"
    # every genotype carries its genome, and its id is the hash of that genome
    for s in genotypes(net):
        assert A.genotype_id(A.bits(s.structure)) == s.id
    assert len(genotypes(net)) == net.extras["analysis"]["genotypes_seen"]
    # every protein carries its coding bits and the triple they fold into
    proteins = [s for s in net.species if s.id.startswith("P")]
    assert proteins
    for s in proteins:
        fields = dict(item.split("=") for item in s.structure.split())
        codons = [int(fields["code"][i:i + 3], 2) for i in range(0, len(fields["code"]), 3)]
        m, w, h, functional = A.fold(codons, net.params["w_max"])
        exact = net.extras["proteins"][s.id]
        assert (m, w, h) == pytest.approx((exact["m"], exact["w"], exact["h"]))
        assert functional == exact["functional"] == bool(int(fields["functional"]))
        # the structure prints the triple to six significant digits
        assert (m, w, h) == pytest.approx(
            (float(fields["m"]), float(fields["w"]), float(fields["h"])), rel=1e-5)
        assert len(codons) == exact["codons"] and s.id.startswith(f"P{len(codons)}-")

    # expression: a genotype catalyses the production of its own proteome
    expressed = 0
    for r in net.reactions:
        if is_expression(r):
            gid, = r.reactants
            assert r.reactants == {gid: 1} and r.products[gid] == 1
            assert {s: n for s, n in r.products.items() if s != gid} == net.extras["proteome"][gid]
            assert r.count >= 1                               # once per individual born with it
            expressed += r.count
    # every individual ever created expresses its genome, unless it has no gene at all
    assert 0 < expressed <= 9 * (SMALL["generations"] + 1)


def test_replication_events_conserve_the_population():
    net = generate_network(ID, seed=3, **SMALL)
    population = Counter({s: int(n) for s, n in net.initial_state.items()})
    assert sum(population.values()) == 9
    kinds = Counter()
    for r in net.reactions:
        if is_expression(r):
            continue
        # P + V -> P + O (the offspring replaces a neighbour) or P -> O (in place)
        kinds[(sum(r.reactants.values()), sum(r.products.values()))] += r.count
        for _ in range(r.count):
            for s, n in r.reactants.items():
                population[s] -= n
            for s, n in r.products.items():
                population[s] += n
    assert set(kinds) <= {(1, 1), (2, 2)}
    assert sum(kinds.values()) == 9 * SMALL["generations"]     # one replication per cell and generation
    population += Counter()                                    # drop the zeros
    assert population == Counter({s: int(n) for s, n in net.extras["final_state"].items()})
    assert sum(population.values()) == 9
    assert Counter(net.extras["space"]["final_grid"]) == population


def test_extras_carry_the_phenotype_the_target_and_the_fitness():
    net = generate_network(ID, seed=3, **SMALL)
    a = net.extras["analysis"]
    assert a["phenotypic_target"]["gaussians"] == A.DEFAULT_TARGET
    target = A.target_function(A.DEFAULT_TARGET, net.params["env_sampling"])
    assert a["phenotypic_target"]["points"] == pytest.approx(target)
    assert a["phenotypic_target"]["area"] == pytest.approx(A.geometric_area(target))
    best = a["best"]
    phenotype = np.array(best["phenotype"])
    assert phenotype.size == net.params["env_sampling"]
    assert best["metabolic_error"] == pytest.approx(A.geometric_area(phenotype - target))
    assert best["fitness"] == pytest.approx(np.exp(-net.params["selection_pressure"] * best["metabolic_error"]))
    assert 0.0 < best["metabolic_error"] < a["phenotypic_target"]["area"]
    assert best["functional_proteins"] >= 1
    assert sum(a["mutations"].values()) > 0
    assert len(a["per_generation"]) == SMALL["generations"] + 1


def test_reproducible_with_a_seed():
    assert generate_network(ID, seed=5, **SMALL).to_dict() == generate_network(ID, seed=5, **SMALL).to_dict()
    assert generate_network(ID, seed=5, **SMALL).to_dict() != generate_network(ID, seed=6, **SMALL).to_dict()


def test_parameters_are_checked():
    with pytest.raises(ValueError, match="gaussian"):
        generate_network(ID, target=[[1.0, 0.5]], **SMALL)
    with pytest.raises(ValueError, match="width"):
        generate_network(ID, target=[[1.0, 0.5, 0.0]], **SMALL)
    with pytest.raises(ValueError, match="odd"):
        generate_network(ID, selection_patch_size=2, **SMALL)
    with pytest.raises(ValueError, match="genome_length"):
        generate_network(ID, genome_length=50, min_genome_length=100, **SMALL)


@pytest.mark.slow
def test_fitness_improves_and_genes_accumulate_over_a_run():
    # the phenomenon of the platform: selection on the distance between phenotype
    # and target drives the metabolic error down while mutations add and duplicate genes
    for seed in (1, 3):
        net = generate_network(ID, seed=seed, grid_width=4, grid_height=4, generations=60,
                               max_genome_length=10000)
        history = net.extras["analysis"]["per_generation"]
        first, last = history[0], history[-1]
        assert first["mean_functional_proteins"] == 1.0       # one clone in every cell
        assert last["best_metabolic_error"] < 0.97 * first["best_metabolic_error"]
        assert last["mean_metabolic_error"] < first["mean_metabolic_error"]
        assert last["best_fitness"] > first["best_fitness"]
        assert last["mean_functional_proteins"] > first["mean_functional_proteins"]
        assert last["mean_genome_length"] > first["mean_genome_length"]   # duplications win
