"""Synthon reproduces the published parts of the Lenaerts-Bersini framework.

Sources of the expectations:

- Banzhaf & Yamamoto (2015), section 18.3.2: molecules as graphs of atoms with
  explicit charges and covalent and ionic bonds; reactions as isomerisations
  between ensembles of molecules.
- Lenaerts & Bersini (2005), "On the generation and analysis of complex reaction
  networks in interstellar chemistry" (the authors' own description of the same
  Synthon framework): table 1 (the eleven reaction classes and their rate
  constants, from Duley & Williams 1984 p. 143), table 2 (generated reaction
  objects), figure 3 (the reaction graphs of radiative association, charge
  transfer and dissociative recombination), figure 4 and its caption (the H/O
  network and the molecule translation), figures 5 and 6 (the MCNG run from
  n(H) = 1000 cm^-3 and n(O) = 0.44 cm^-3, in which H and H2 dominate).

The Artificial Life 15(1) paper itself could not be obtained, so none of its
numbers are asserted here (see the catalog entry's decisions).
"""

from collections import Counter

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import synthon as S

ID = "synthon"

#: figure 4 caption: the paper's own notation -> the chemical formula it means.
FIGURE_4_CAPTION = {
    "H": "H^+", "*H": "H", "HH": "H2", "H.*H": "H2^+", "H.HH": "H3^+",
    "*O*": "O", "HO": "HO^+", "*OH": "HO", "HOH": "H2O", "*OO*": "O2",
    "*O(H)H": "H2O^+", "HO(H)H": "H3O^+", "*": "e-",
}


def formulas(net) -> set[str]:
    return set(net.extras["formulas"].values())


def reactions_by_formula(net) -> set[tuple[tuple[str, ...], tuple[str, ...]]]:
    out = set()
    f = net.extras["formulas"]
    for r in net.reactions:
        left = tuple(sorted(f[s] for s, n in r.reactants.items() for _ in range(n)))
        right = tuple(sorted(f[s] for s, n in r.products.items() for _ in range(n)))
        out.add((left, right))
    return out


# --- the published reaction classes -------------------------------------------------
def test_reaction_classes_are_table_1():
    """The eleven classes G1-G11 with the templates and rates of table 1."""
    net = generate_network(ID)
    printed = [(c["class"], c["name"], c["template"], c["k"], c["units"])
               for c in net.extras["reaction_classes"]]
    assert printed == [
        ("G1", "cosmic ray ionization", "A + c.r. -> A+ + e-", 1e-17, "s^-1"),
        ("G2", "cosmic ray ionization", "AB + c.r. -> AB+ + e-", 1e-17, "s^-1"),
        ("G3", "cosmic ray ionization", "AB + c.r. -> A + B+ + e-", 1e-19, "s^-1"),
        ("G4", "ion-molecule exchange", "A+ + BC -> AB+ + C", 1e-9, "cm^3 s^-1"),
        ("G5", "charge rearrangement", "A+ + B -> A + B+", 1e-9, "cm^3 s^-1"),
        ("G6", "dissociative recombination", "AB+ + e- -> A + B", 1e-6, "cm^3 s^-1"),
        ("G7", "dissociative recombination", "A+ + e- -> A", 1e-11, "cm^3 s^-1"),
        ("G8", "neutral reactions", "A + BC -> AB + C", 1e-11, "cm^3 s^-1"),
        ("G9", "photodissociation", "AB + hv -> A + B", 1e-11, "s^-1"),
        ("G10", "photodissociation", "H2 + hv -> 2 H", 1e-14, "s^-1"),
        ("G11", "grain surface reaction", "H + H :g -> H2 + g", 1e-17, "cm^3 s^-1"),
    ]
    # every class contributes at least one reaction object to the H/O network
    used = {c for classes in net.extras["reaction_classes_used"] for c in classes}
    assert used == {c[0] for c in printed}
    # the constants reach the reactions themselves
    published = {c["class"]: (c["k"], c["units"]) for c in net.extras["reaction_classes"]}
    for r in net.reactions:
        classes = r.rate["classes"].split(",")
        assert r.rate["k"] == pytest.approx(sum(published[c][0] for c in classes))
        assert r.rate["units"] == published[classes[0]][1]
    single = {r.rate["classes"]: r.rate["k"] for r in net.reactions if "," not in r.rate["classes"]}
    assert single["G1"] == 1e-17 and single["G11"] == 1e-17 and single["G10"] == 1e-14


# --- figure 4: the H/O network -------------------------------------------------------
def test_figure_4_molecules_are_generated_from_H_and_O():
    """From the atoms H and O alone, the thirteen molecules of figure 4 appear."""
    net = generate_network(ID, seed=1)
    assert net.status == "complete"
    assert set(FIGURE_4_CAPTION.values()) <= formulas(net)
    # the observational constraint of the figure holds: nothing else is generated
    assert formulas(net) == set(S.FIGURE_4)
    assert net.initial_state == {"H*": 1000.0, "O::**": 0.44}      # Duley & Williams p. 143


def test_species_ids_are_the_papers_own_notation():
    """The caption's notation: HH = H2, H.*H = H2+, H.HH = H3+, HOH = H2O, HO = OH+."""
    net = generate_network(ID, seed=1)
    ids = {s.id for s in net.species}
    assert {"H*", "H", "HH", "H.H*", "H.HH", "HO::H", "HO::", "HO::*", "O::**", "e-"} <= ids
    formula = net.extras["formulas"]
    assert formula["H"] == "H^+" and formula["H*"] == "H"          # *H = H, H = H+
    assert formula["H.H*"] == "H2^+" and formula["H.HH"] == "H3^+"  # the ionic representation
    assert formula["HO::"] == "HO^+" and formula["HO::*"] == "HO"
    assert formula["HO::H"] == "H2O" and formula["O::*O::*"] == "O2"
    # id and structure are the same canonical code
    assert all(s.structure == s.id for s in net.species)


def test_table_2_reaction_objects():
    """Table 2: object 1 (H + c.r. -> H+ + e-) and object 32 (H2+ + H2 -> H3+ + H)."""
    net = generate_network(ID, seed=1)
    found = reactions_by_formula(net)
    assert (("H",), ("H^+", "e-")) in found
    assert (("H2", "H2^+"), ("H", "H3^+")) in found
    # object 8, H+ + H -> H + H+, is a reaction object of class G5 but the identity
    # as a stoichiometric reaction, so it is not part of the network.
    assert (("H", "H^+"), ("H", "H^+")) not in found
    proton, hydrogen = S.Mol(((("H"), 0, 0),)), S.atom("H")
    limits = S.Limits(4, 10, 7, 1, None)
    objects = [tuple(sorted(m.formula for m in products))
               for products, cls in S._outcomes((proton, hydrogen), {"G5"}, limits)]
    assert ("H", "H^+") in objects


# --- isomerisation: the defining property (book 18.3.2) -------------------------------
def test_every_reaction_is_an_isomerisation():
    """Educt and product ensembles are isomers: same atoms, same electrons, same charge."""
    net = generate_network(ID, seed=1)
    by_id = {s.id: s for s in net.species}
    atoms = net.extras["conservation"]
    electrons, charges = net.extras["electrons"], net.extras["charges"]
    assert len(by_id) == len(net.species)
    for r in net.reactions:
        for law in atoms:
            left = sum(law["vector"][s] * n for s, n in r.reactants.items())
            right = sum(law["vector"][s] * n for s, n in r.products.items())
            assert left == right, (r.to_text(), law["name"])
        assert (sum(electrons[s] * n for s, n in r.reactants.items())
                == sum(electrons[s] * n for s, n in r.products.items()))
        assert (sum(charges[s] * n for s, n in r.reactants.items())
                == sum(charges[s] * n for s, n in r.products.items()))


def test_conservation_vectors_annihilate_the_stoichiometric_matrix():
    net = generate_network(ID, seed=1)
    ids, R, P = net.matrices()
    stoichiometry = (P - R).toarray()
    names = {law["name"] for law in net.extras["conservation"]}
    assert names == {"atoms of H", "atoms of O", "electrons", "charge"}
    for law in net.extras["conservation"]:
        vector = np.array([law["vector"][i] for i in ids])
        assert np.all(vector @ stoichiometry == 0), law["name"]


# --- canonical identity (the role of CANGEN) ------------------------------------------
def test_canonical_code_does_not_depend_on_how_the_molecule_is_built():
    water = S.molecule([("O", 2, 0), ("H", 0, 0), ("H", 0, 0)], {(0, 1): 1, (0, 2): 1}, set())
    same = S.molecule([("H", 0, 0), ("O", 2, 0), ("H", 0, 0)], {(1, 2): 1, (0, 1): 1}, set())
    assert water == same and water.id == same.id == "HO::H"
    assert water.formula == "H2O" and water.charge == 0 and water.electrons == 8


def test_electronic_isomers_are_different_species():
    """Same atoms and electrons, a different lone-pair/radical split: another synthon."""
    water = S.molecule([("O", 2, 0), ("H", 0, 0), ("H", 0, 0)], {(0, 1): 1, (0, 2): 1}, set())
    biradical = S.molecule([("O", 1, 2), ("H", 0, 0), ("H", 0, 0)], {(0, 1): 1, (0, 2): 1}, set())
    assert water.formula == biradical.formula == "H2O"
    assert water.electrons == biradical.electrons
    assert water.id != biradical.id and water != biradical


def test_charges_follow_the_dugundji_ugi_bookkeeping():
    # q = v - 2*lone_pairs - radicals - sum of bond orders
    assert S.atom("H").charge == 0 and S.atom("O").charge == 0
    assert S.Mol((("H", 0, 0),)).charge == 1                        # the proton
    assert S.Mol((("O", 2, 0), ("H", 0, 0)), ((0, 1, 1),)).charge == 1      # OH+
    assert S.Mol((("O", 2, 1), ("H", 0, 0)), ((0, 1, 1),)).charge == 0      # OH
    assert S.ELECTRON.charge == -1 and S.ELECTRON.electrons == 1


# --- the generators --------------------------------------------------------------------
def test_closure_is_complete_and_reproducible():
    a = generate_network(ID, seed=5)
    b = generate_network(ID, seed=5)
    assert a.to_dict() == b.to_dict()
    assert a.status == "complete"
    # N = 6, the value of the published experiment, gives the same network as the default
    assert generate_network(ID, max_atoms=6).to_dict()["species"] == a.to_dict()["species"]


def test_species_budget_truncates_the_closure():
    small = generate_network(ID, max_species=8, seed=1)
    assert small.status == "truncated" and len(small.species) <= 8


def test_kinetics_prunes_the_network_and_H_and_H2_dominate():
    """MCNG (figures 5 and 6): few reactions fire, and H and H2 are the dominant species."""
    net = generate_network(ID, method="kinetic", steps=2000, seed=3)
    assert net.status == "observed"
    assert all(r.count is not None and r.count > 0 for r in net.reactions)
    assert sum(r.count for r in net.reactions) == 2000
    # the observed network is much smaller than the closure it is drawn from
    assert len(net.reactions) < len(generate_network(ID, seed=3).reactions) / 4
    counts: Counter = Counter()
    for species, n in net.extras["final_state"].items():
        counts[net.extras["formulas"][species]] += n
    assert [f for f, _ in counts.most_common(2)] == ["H", "H2"]
    assert net.extras["elapsed_time"] > 0


def test_figure_3_classes_have_no_published_rates():
    net = generate_network(ID, templates="fig3", seed=1)
    assert [c["class"] for c in net.extras["reaction_classes"]] == ["ra", "ct", "dr"]
    assert all(c["k"] is None for c in net.extras["reaction_classes"])
    assert all(r.rate is None for r in net.reactions)
    # radiative association A + B -> AB builds H2, OH, H2O and O2 from the atoms
    assert {"H2", "HO", "H2O", "O2"} <= formulas(net)
    with pytest.raises(ValueError, match="publishes none"):
        generate_network(ID, templates="fig3", method="kinetic")


@pytest.mark.slow
def test_without_the_observational_constraint_the_closure_explodes():
    """The drawback the paper reports: the formal generator grows combinatorially."""
    free = generate_network(ID, allowed=[], max_atoms=3, max_species=400, seed=1)
    constrained = generate_network(ID, max_atoms=3, seed=1)
    assert len(free.species) > 4 * len(constrained.species)
    assert len(free.reactions) > 10 * len(constrained.reactions)
    # species the observational constraint of figure 4 hides
    assert {"O3", "HO2", "H^-"} <= formulas(free)


# --- parameters -------------------------------------------------------------------------
def test_bad_parameters():
    with pytest.raises(ValueError, match="unknown parameter"):
        generate_network(ID, not_a_parameter=1)
    with pytest.raises(ValueError, match="max_atoms"):
        generate_network(ID, max_atoms=0)
    with pytest.raises(ValueError, match="unknown element"):
        generate_network(ID, initial=["Xx"])
    with pytest.raises(ValueError, match="not in allowed"):
        generate_network(ID, initial=["H", "N"])
    with pytest.raises(ValueError, match="densities"):
        generate_network(ID, densities=[1.0])
