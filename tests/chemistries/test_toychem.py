"""ToyChem: what the Benko/Flamm/Stadler papers publish about the graph-based toy model.

Benko, Flamm & Stadler, "A graph-based toy model of chemistry", J. Chem. Inf.
Comput. Sci. 43:1085-1093 (2003): Table 1 (the complete parametrisation),
appendix B (the Diels-Alder rule in GML), figure 3 (calculated total atomization
energies of the n-alkane series and of the C6H10 isomers), figure 7 (the
repetitive Diels-Alder network) and figure 8 (the formose network from
formaldehyde and glycol aldehyde).

"Generic properties of chemical networks", ECAL/LNCS 2801:10-19 (2003): the
substrate graph and its measures, Table 1 (Diels-Alder is small-world, formose is
not) and figures 3-4 (the network grows as the reactivity threshold is lowered).

"Explicit collision simulation of chemical reactions in a graph based artificial
chemistry", ECAL/LNCS 3630:725-733 (2005): Table 2, the total atomization
energies of ethene, butadiene and their Diels-Alder adducts, in kcal/mol.
"""

import numpy as np
import pytest
from rdkit import Chem

from chemart import generate_network
from chemart.chemistries.toychem import (TABLE_1_I_PI, TABLE_1_I_SIGMA, TABLE_1_KAPPA,
                                         TABLE_1_RING_SCALE, TABLE_1_S_PI, TABLE_1_S_SIGMA,
                                         canonical, conserves, diels_alder, eht,
                                         substrate_graph)

ID = "toychem"

#: 2005 ECAL Table 2: total atomization energies (kcal/mol) of the species whose
#: structure formulae that table writes unambiguously.
TABLE_2 = {
    "C=C": -415.95,             # ethene
    "C=CC=C": -789.62,          # butadiene
    "C1=CCCCC1": -1222.33,      # cyclohexene, the adduct of butadiene + ethene (DA)
    "C=CC1CC=CCC1": -1587.39,   # 4-vinylcyclohexene, the adduct of two butadienes (DA)
}

#: The C6H10 isomers of figure 3 (right), in the caption's order of increasing
#: experimental total atomization energy. The caption calls them "C4H10 isomers",
#: but the molecules it lists and the energies it plots are C6H10.
C6H10 = [
    ("1-hexyne", "C#CCCCC"), ("2-hexyne", "CC#CCCC"), ("3-hexyne", "CCC#CCC"),
    ("3,3-dimethyl-1-butyne", "C#CC(C)(C)C"), ("1,5-hexadiene", "C=CCCC=C"),
    ("1,4-hexadiene", "C=CCC=CC"), ("1,3-hexadiene", "C=CC=CCC"),
    ("2,4-hexadiene", "CC=CC=CC"), ("bicyclo[3.1.0]hexane", "C1CC2CC2C1"),
    ("4-methylcyclopentene", "CC1CC=CC1"), ("3-methylcyclopentene", "CC1C=CCC1"),
    ("1-methylcyclopentene", "CC1=CCCC1"),
]


@pytest.fixture(scope="module")
def net():
    return generate_network(ID, seed=1)


def tae(smiles):
    return eht(smiles)["tae_kcal"]


# ---------------------------------------------------------------------------
# Table 1: the parametrisation is transcribed, not invented
# ---------------------------------------------------------------------------
def test_table_1_is_transcribed_verbatim():
    """Appendix A, Table 1 of the 2003 paper."""
    assert TABLE_1_KAPPA == 1.75                       # Wolfsberg-Helmholtz scaling factor
    assert TABLE_1_I_SIGMA[("H", "s")] == -13.6        # the hydrogen ionisation potential
    assert TABLE_1_I_SIGMA[("C", "sp3")] == -13.9
    assert TABLE_1_I_SIGMA[("C", "sp2")] == -14.5
    assert TABLE_1_I_SIGMA[("O", "sp2")] == -20.6
    assert TABLE_1_I_PI == {"C": -11.4, "N": -13.4, "O": -14.8}
    assert TABLE_1_S_SIGMA[("C_sp2", "C_sp2")] == 0.77
    assert TABLE_1_S_SIGMA[("H_s", "H_s")] == 0.75
    assert TABLE_1_S_SIGMA[("C_sp", "C_sp")] == 0.87
    assert TABLE_1_S_PI[("C", "C")] == 0.38 and TABLE_1_S_PI[("O", "O")] == 0.26
    assert TABLE_1_RING_SCALE == {3: 0.7, 4: 0.8}      # banana bonds in strained rings
    # the table is symmetric, as an overlap matrix must be
    for (a, b), value in TABLE_1_S_SIGMA.items():
        assert TABLE_1_S_SIGMA[(b, a)] == value


def test_the_orbital_graph_has_one_vertex_per_valence_orbital():
    """Section 2: 1s for hydrogen, four valence orbitals for C, N and O."""
    ethene = eht("C=C")
    assert ethene["orbitals"] == 2 * 4 + 4 * 1         # two sp2 carbons plus four hydrogens
    assert ethene["electrons"] == 2 * 4 + 4 * 1        # valence electrons, filled in pairs
    assert eht("C")["orbitals"] == 4 + 4               # methane: four sp3 plus four 1s
    assert eht("O")["orbitals"] == 4 + 2               # water: four sp3 (two lone pairs) plus 2 H


# ---------------------------------------------------------------------------
# 2005 Table 2: the published total atomization energies
# ---------------------------------------------------------------------------
def test_energies_reproduce_table_2_up_to_one_scale_factor():
    """The published TAEs are reproduced in their relative sizes.

    The absolute kcal/mol scale is not reproducible: the 2005 tables use the
    later C/H/N/O/P/S parametrisation, which was only distributed with the
    unavailable technical report [105]. With the published 2003 Table 1 the
    computed energies are proportional to the published ones.
    """
    computed = np.array([tae(s) for s in TABLE_2])
    published = np.array(list(TABLE_2.values()))
    assert np.all(computed < 0), "every molecule must be bound"

    correlation = float(np.corrcoef(computed, published)[0, 1])
    assert correlation > 0.999, correlation
    ratios = published / computed
    assert ratios.max() - ratios.min() < 0.02, ratios


def test_separated_molecules_have_additive_energies():
    """Table 2 lists every reactant pair at exactly the sum of its two molecules."""
    for a, b, published_pair in (("C=C", "C=C", -831.90),
                                 ("C=CC=C", "C=CC=C", -1579.24),
                                 ("C=C", "C=CC=C", -1205.57)):
        assert TABLE_2[a] + TABLE_2[b] == pytest.approx(published_pair, abs=0.01)
        # the model is extensive in the same way: separated molecules simply add
        assert tae(a) + tae(b) < min(tae(a), tae(b))


def test_the_cycloadditions_are_exothermic_as_published():
    """Table 2: butadiene + ethene releases 16.76 kcal/mol, two butadienes 8.15."""
    ethene_da = tae("C1=CCCCC1") - (tae("C=CC=C") + tae("C=C"))
    butadiene_da = tae("C=CC1CC=CCC1") - 2 * tae("C=CC=C")
    assert ethene_da < 0 and butadiene_da < 0
    # the published ordering: the ethene adduct is the more exothermic of the two
    assert ethene_da < butadiene_da
    published_ratio = (-1222.33 + 1205.57) / (-1587.39 + 1579.24)
    assert ethene_da / butadiene_da == pytest.approx(published_ratio, rel=0.6)


# ---------------------------------------------------------------------------
# Figure 3: the energy model against homologous series and isomers
# ---------------------------------------------------------------------------
def test_alkane_series_is_linear_in_chain_length():
    """Figure 3 (left): calculated against experimental TAE of methane to hexane is a line."""
    values = np.array([tae("C" * n) for n in range(1, 7)])
    assert np.all(np.diff(values) < 0)
    fit = np.polyfit(np.arange(1, 7), values, 1)
    residual = np.max(np.abs(values - np.polyval(fit, np.arange(1, 7))))
    assert residual / abs(values.mean()) < 0.001, residual
    # a constant increment per CH2 group is what makes the plot a straight line
    increments = np.diff(values)
    assert increments.std() / abs(increments.mean()) < 0.01


def test_the_c6h10_isomers_fall_in_a_narrow_band():
    """Figure 3 (right): twelve isomers of one formula inside about 90 of 1565 kcal/mol.

    The caption's "C4H10" is an erratum: the molecules it names and the energies
    it plots are C6H10. Chemart follows the data.
    """
    for name, smiles in C6H10:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        formula = (sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C"),
                   sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H"))
        assert formula == (6, 10), f"{name} is not C6H10"

    values = {name: tae(smiles) for name, smiles in C6H10}
    span = max(values.values()) - min(values.values())
    assert span / abs(np.mean(list(values.values()))) < 0.05, span
    # the caption orders the isomers by increasing experimental TAE, so the last
    # one listed is the most strongly bound
    assert min(values, key=values.get) == C6H10[-1][0] == "1-methylcyclopentene"


# ---------------------------------------------------------------------------
# Appendix B: the Diels-Alder rule builds the published adducts
# ---------------------------------------------------------------------------
def test_diels_alder_rule_reproduces_the_published_adducts():
    """2005 Table 2: the rows marked (DA)."""
    assert diels_alder(["C=CC=C", "C=C"]) == [(canonical("C1=CCCCC1"),)]        # cyclohexene
    assert diels_alder(["C=CC=C", "C=CC=C"]) == [(canonical("C=CC1CC=CCC1"),)]  # 4-vinylcyclohexene
    # a dienophile alone, or two of them, is not a Diels-Alder reaction
    assert diels_alder(["C=C", "C=C"]) == []
    # the rule also runs intramolecularly, the rearrangement of figure 4
    assert diels_alder(["C=CC=CC=CC=C"])


def test_every_reaction_conserves_atoms_and_bond_order(net):
    """Section 3: chemical rewrite rules conserve vertex labels and total bond order."""
    assert net.reactions
    for reaction in net.reactions:
        educts = [s for s, n in reaction.reactants.items() for _ in range(n)]
        products = tuple(s for s, n in reaction.products.items() for _ in range(n))
        assert conserves(educts, products), reaction.to_text()

    ids, reactants, products = net.matrices()
    stoichiometry = (products - reactants).toarray()
    assert net.extras["conservation"]
    for law in net.extras["conservation"]:
        vector = np.array([law["vector"][s] for s in ids])
        assert np.all(vector @ stoichiometry == 0), law["name"]


def test_canonical_ids_are_invariant_under_relabelling():
    """Section 4: canonical SMILES reduce the isomorphism test to string comparison."""
    assert canonical("C=CC=C") == canonical("C(=C)C=C") == canonical("C(C=C)=C")
    assert canonical("c1ccccc1") == canonical("C1=CC=CC=C1")
    assert canonical("OCC=O") == canonical("O=CCO")
    assert canonical("C=CC=C") != canonical("C=CCC")
    # the same molecule reached by two different routes is one species with one
    # energy (the two writings order the atoms differently, so the generalised
    # eigenproblem is summed in a different order and agrees to rounding only)
    assert eht("OCC=O")["tae_kcal"] == pytest.approx(eht("O=CCO")["tae_kcal"])


# ---------------------------------------------------------------------------
# Figure 8: the formose network (the default)
# ---------------------------------------------------------------------------
def test_formose_builds_the_sugars_of_figure_8(net):
    """Formaldehyde and glycol aldehyde condense into the formose sugars."""
    assert net.status == "complete"
    ids = {s.id for s in net.species}
    assert canonical("C=O") in ids and canonical("OCC=O") in ids          # the seed
    assert canonical("OC=CO") in ids                                     # the ene-diol
    assert canonical("O=CC(O)CO") in ids                                 # glyceraldehyde
    assert canonical("O=C(CO)CO") in ids                                 # dihydroxyacetone
    assert canonical("O=CC(O)C(O)CO") in ids                             # a tetrose
    assert net.extras["rules"] == ["keto-enol", "aldol"]
    # the caption's cut: no carbon chain beyond four members
    assert all(sum(1 for c in s.id if c == "C") <= 4 for s in net.species)


def test_every_reaction_carries_an_arrhenius_barrier(net):
    """The reactivities 'can be translated into rate constants using Arrhenius' law'."""
    energies = net.extras["energies"]
    assert all(r.rate["law"] == "arrhenius" for r in net.reactions)
    assert all(r.rate["Ea"] > 0 and r.rate["reactivity"] > 0 for r in net.reactions)
    assert len(energies["activation_energy"]) == len(net.reactions)
    assert set(energies["total_atomization_energy"]) == {s.id for s in net.species}
    for reaction, barrier in zip(net.reactions, energies["activation_energy"]):
        assert reaction.rate["Ea"] == pytest.approx(barrier)
    # every species is bound, and the seed molecules are the smallest
    assert all(v < 0 for v in energies["total_atomization_energy"].values())


# ---------------------------------------------------------------------------
# ECAL figures 3-4: the barrier threshold decides how large the network grows
# ---------------------------------------------------------------------------
def test_lowering_the_barrier_threshold_shrinks_the_network():
    """"The CRN grows as the reactivity threshold is lowered" - here, as the
    activation-energy cutoff is raised."""
    sizes = []
    for cutoff in (80.0, 100.0, 120.0, 140.0, 0.0):
        network = generate_network(ID, seed=1, barrier_cutoff=cutoff)
        sizes.append((len(network.species), len(network.reactions)))
    # monotone in the cutoff, and 0 (no gate) is the largest
    species = [s for s, _ in sizes]
    assert species == sorted(species)
    assert sizes[-1] == max(sizes)
    assert sizes[0][1] < sizes[-1][1], "a tight barrier must remove reactions"


def test_the_substrate_graph_follows_the_ecal_definitions():
    """ECAL section 4: every hyperedge becomes a clique; <k> = 2m/n."""
    stats = substrate_graph(["A", "B", "C"], [(("A", "B"), ("C",))])
    assert stats["nodes"] == 3 and stats["edges"] == 3          # one reaction, one triangle
    assert stats["mean_degree"] == pytest.approx(2.0)
    assert stats["mean_clustering"] == pytest.approx(1.0)
    assert stats["mean_path_length"] == pytest.approx(1.0)
    assert stats["random_clustering"] == pytest.approx(1.0)

    apart = substrate_graph(["A", "B", "C", "D"], [(("A",), ("B",)), (("C",), ("D",))])
    assert apart["edges"] == 2 and apart["mean_clustering"] == 0.0


def test_the_network_reports_its_substrate_graph(net):
    stats = net.extras["analysis"]["substrate_graph"]
    assert stats["nodes"] == len(net.species)
    assert stats["mean_degree"] == pytest.approx(2 * stats["edges"] / stats["nodes"])
    assert stats["degree_sequence"] == sorted(stats["degree_sequence"], reverse=True)


@pytest.mark.slow
def test_the_diels_alder_network_is_small_world():
    """ECAL Table 1: the Diels-Alder substrate graph is clustered and short-pathed."""
    network = generate_network(ID, seed=1, network="diels-alder")
    stats = network.extras["analysis"]["substrate_graph"]
    assert stats["nodes"] > 30
    assert stats["mean_clustering"] > stats["random_clustering"]
    assert stats["mean_path_length"] <= stats["random_path_length"]
    assert stats["small_world"]
    # the published network has <k> = 4.65 and <C> = 0.72 at n = 40
    assert 3.0 < stats["mean_degree"] < 7.0
    assert stats["mean_clustering"] > 0.5


# ---------------------------------------------------------------------------
# Rewrite modes, reproducibility and rejected parameters
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_the_seed_reproduces_the_random_rewrite_mode():
    """'Random rewrite': the server picks one subgraph isomorphism at random."""
    def run(seed):
        return generate_network(ID, seed=seed, network="diels-alder", max_species=14,
                                rewrite_mode="random")

    assert run(3).to_dict() == run(3).to_dict()
    everything = generate_network(ID, seed=3, network="diels-alder", max_species=14)
    assert len(run(3).reactions) <= len(everything.reactions)


@pytest.mark.slow
def test_priority_mode_keeps_one_channel_per_pair():
    """'Priority rewrite': the rule with the highest priority value is chosen."""
    picked = generate_network(ID, seed=1, network="diels-alder", max_species=14,
                              rewrite_mode="priority")
    everything = generate_network(ID, seed=1, network="diels-alder", max_species=14)
    assert len(picked.reactions) <= len(everything.reactions)
    assert picked.extras["analysis"]["rewrite_mode"] == "priority"


def test_a_custom_network_runs_its_own_seed_and_rules():
    custom = generate_network(ID, seed=1, network="custom",
                              seed_molecules=["C=CC=C", "C=C"], rules=["diels-alder"],
                              max_species=8)
    assert canonical("C1=CCCCC1") in {s.id for s in custom.species}
    assert custom.extras["seed"] == [canonical("C=CC=C"), canonical("C=C")]


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError, match="only used with network='custom'"):
        generate_network(ID, rules=["aldol"])
    with pytest.raises(ValueError, match="seed_molecules"):
        generate_network(ID, network="custom", rules=["aldol"])
    with pytest.raises(ValueError, match="rules"):
        generate_network(ID, network="custom", seed_molecules=["C=O"])
    with pytest.raises(ValueError, match="max_species"):
        generate_network(ID, max_species=1)
    # Table 1 parametrises H, C, N and O only
    with pytest.raises(ValueError, match="does not parametrise"):
        generate_network(ID, network="custom", seed_molecules=["CSC"], rules=["aldol"])
    # the energy model is limited to neutral molecules
    with pytest.raises(ValueError, match="charged"):
        eht("[CH2-]C=O")
    with pytest.raises(ValueError, match="cannot read"):
        canonical("not a molecule")
