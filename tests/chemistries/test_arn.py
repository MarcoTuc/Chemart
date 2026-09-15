"""ARN (book 11.3.4; Banzhaf 2003; Kuo & Banzhaf 2004; Kuo, Leier & Banzhaf 2004, 2006)."""

import numpy as np
import pytest
from odes import integrate, rhs

from chemart import generate_network
from chemart.chemistries.arn import bits, duplication_divergence, genes, majority, promoters

PROMOTER = bits("01010101")


def structure(sp):
    return dict(item.split("=") for item in sp.structure.split())


# --- genome decoding (Kuo & Banzhaf 2004; Kuo, Leier & Banzhaf 2004 Fig. 1) --
def test_promoter_rules_and_gene_layout():
    # "0101010101" (and any longer periodic run) is one promoter at its first bit
    assert promoters(bits("1101010101011"), PROMOTER).tolist() == [2]
    assert promoters(bits("0101010101010101"), PROMOTER).tolist() == [0]
    assert promoters(bits("01010101" + "11" + "01010101"), PROMOTER).tolist() == [0, 10]

    # [inhibitor 32][enhancer 32][promoter 8][5 x 32 coding bits]
    inhibitor, enhancer = bits("1" * 32), bits("0011" * 8)
    coding = bits("1110" * 40)
    genome = np.concatenate([inhibitor, enhancer, PROMOTER, coding, bits("111")])
    (pos, inh, enh, cod), = genes(genome, PROMOTER)
    assert pos == 64
    assert (inh == inhibitor).all() and (enh == enhancer).all() and (cod == coding).all()
    # a promoter without room for both sites upstream or 160 coding bits downstream is no gene
    assert genes(np.concatenate([PROMOTER, coding]), PROMOTER) == []
    assert genes(genome[:64 + 8 + 159], PROMOTER) == []


def test_majority_voting_translation():
    # book Fig. 11.16: protein bit k = majority of bit k over the 5 segments
    segments = np.array([[1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]], dtype=np.uint8)
    assert majority(segments.ravel(), 4, np.random.default_rng(0)).tolist() == [1, 1, 0, 0]


# --- eqs. 11.12-11.14 -------------------------------------------------------
@pytest.mark.parametrize("seed", [1, 5])
def test_reactions_reproduce_equations_11_12_to_11_14(seed):
    beta, delta = 0.7, 2.5
    net = generate_network("arn", seed=seed, beta=beta, delta=delta)
    N = len(net.species)
    assert N >= 3
    s = [structure(sp) for sp in net.species]
    prot = np.array([bits(x["protein"]) for x in s])
    enh = np.array([bits(x["enhancer"]) for x in s])
    inh = np.array([bits(x["inhibitor"]) for x in s])
    # d_j: number of complementary bits between protein j and the site of gene i
    Ea = np.exp(beta * ((enh[:, None, :] != prot[None, :, :]).sum(2) - 32))
    Eh = np.exp(beta * ((inh[:, None, :] != prot[None, :, :]).sum(2) - 32))

    c = np.random.default_rng(seed).dirichlet(np.ones(N))
    a, h = Ea @ c / N, Eh @ c / N
    book = delta * (a - h) * c / c.sum()                               # eq. 11.12

    ids, f = rhs(net)
    assert ids == [sp.id for sp in net.species]
    assert f(0.0, c) == pytest.approx(book - c * book.sum(), rel=1e-9, abs=1e-18)  # on the simplex
    net.outflow = None
    assert rhs(net)[1](0.0, c) == pytest.approx(book, rel=1e-9, abs=1e-18)       # production terms alone
    assert net.initial_state == {sp.id: 1.0 / N for sp in net.species}


# --- genome statistics (Kuo & Banzhaf 2004 Figs. 3-5; Kuo et al. 2006 Figs. 5-7)
def test_random_genomes_of_131072_bits_carry_340_to_440_genes():
    counts = [len(genes(np.random.default_rng(s).integers(0, 2, 131072, dtype=np.uint8), PROMOTER))
              for s in range(30)]
    assert 340 <= min(counts) and max(counts) <= 440
    net = generate_network("arn", seed=0, genome_length=131072, link_threshold=33)
    assert 340 <= net.extras["analysis"]["genes"] <= 440


def test_duplication_divergence_gene_numbers_broad_at_1_percent_random_like_at_5():
    def counts(rate):
        return np.array([len(genes(duplication_divergence(np.random.default_rng(s), 12, rate), PROMOTER))
                         for s in range(40)])

    one, five = counts(0.01), counts(0.05)
    assert duplication_divergence(np.random.default_rng(0), 12, 0.01).size == 131072
    # Fig. 3: 0-3000 genes, power-law-like; Fig. 4: 200-700 genes, close to random genomes
    assert one.max() <= 3000 and one.min() < 200 and one.max() > 700
    assert 200 <= five.min() and five.max() <= 700
    assert one.std() / one.mean() > 3 * five.std() / five.mean()


# --- static topology (book 11.3.4 "Exploring ARNs"; Kuo et al. 2006 Figs. 2-4)
def test_link_threshold_sharp_transition_and_nested_graphs():
    def links(th):
        net = generate_network("arn", seed=3, genome_length=16384, link_threshold=th)
        # products tell regulator from target: P_j + P_i -> P_j + 2 P_i or -> P_j
        return net, {(tuple(sorted(r.reactants.items())), tuple(sorted(r.products.items())), r.rate["site"])
                     for r in net.reactions}

    full, all_links = links(0)
    N = len(full.species)
    assert len(all_links) == len(full.reactions) == 2 * N * N          # fully connected at threshold 0
    fractions = [len(links(th)[1]) / (2 * N * N) for th in range(0, 33, 4)]
    assert fractions[0] == 1.0 and fractions[-1] == 0.0
    assert all(x >= y for x, y in zip(fractions, fractions[1:]))
    assert fractions[2] > 0.3 and fractions[6] < 0.01                  # 1 -> 0 between 8 and 24 bits
    assert links(22)[1] <= links(21)[1]                                # Fig. 3 is a subgraph of Fig. 2
    assert all(r.rate["match"] >= 21 for r in links(21)[0].reactions)


# --- dynamics (Banzhaf 2003 section 3, Figs. 3-4) ---------------------------
def test_equal_start_settles_with_a_few_dominant_proteins():
    for seed in range(4):
        net = generate_network("arn", seed=seed, genome_length=2048)
        N = len(net.species)
        xs, t = integrate(net, 2e7, t_eval=[0.0, 1.8e7, 2e7])
        X = np.array([xs[sp.id] for sp in net.species])
        assert np.allclose(X.sum(axis=0), 1.0, atol=1e-6)               # stays on the simplex
        assert X[:, 2].max() > max(0.4, 2.0 / N)                        # a protein dominates
        assert np.abs(X[:, 2] - X[:, 1]).max() < 1e-2                   # point attractor reached
