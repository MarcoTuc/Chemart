"""String Metabolic Network (book 18.3.1, Ono, Fujiwara & Yuta 2005): reaction classes, conservation, mutation, flow."""

from collections import Counter

import numpy as np
import pytest
from chemart.simulate import integrate

from chemart import generate_network
from chemart.chemistries.smn import letter_mass


def stoich(side):
    return dict(Counter(side))


# --- eqs. 18.17-18.20: the book's example reactions ---------------------------------
def test_book_examples_eq_18_19_and_18_20():
    net = generate_network("smn", genome=["adbg + ef", "facb|aha + ee|fg"])
    ligation, recombination = net.extras["genome"]
    assert ligation["reaction"] == "adbg + ef <-> adbgef"                       # eq. 18.19
    assert recombination["reaction"] == "facbaha + eefg <-> facbfg + eeaha"      # eq. 18.20
    texts = set(net.to_text().splitlines())
    assert texts == {
        "adbg + ef -> adbgef", "adbgef -> adbg + ef",
        "facbaha + eefg -> facbfg + eeaha", "facbfg + eeaha -> facbaha + eefg",
    }
    assert net.outflow == "constant-total"
    assert net.initial_state == {"adbg": 1.0, "ef": 1.0, "facbaha": 1.0, "eefg": 1.0}
    assert all(r.rate is None for r in net.reactions)


# --- eqs. 18.17-18.18: every enzyme is a reversible ligation or recombination -------
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_every_reaction_is_a_ligation_or_recombination(seed):
    net = generate_network("smn", seed=seed, n_enzymes=30, mutations=40)
    species = {s.id for s in net.species}
    assert len(net.extras["genome"]) == 70
    for e in net.extras["genome"]:
        fwd, bwd = (net.reactions[i] for i in e["reactions"])
        assert (fwd.reactants, fwd.products) == (bwd.products, bwd.reactants)   # reversible
        left, right = (part.strip() for part in e["educts"].split(" + "))
        if e["kind"] == "ligation":
            assert "|" not in e["educts"]
            assert fwd.reactants == stoich([left, right]) and fwd.products == stoich([left + right])
        else:
            a, b = left.split("|")
            c, d = right.split("|")
            assert fwd.reactants == stoich([a + b, c + d])
            assert fwd.products == stoich([a + d, c + b])                          # AB + CD -> AD + CB
        assert set(fwd.reactants) | set(fwd.products) <= species
    assert {e["kind"] for e in net.extras["genome"]} == {"ligation", "recombination"}
    assert max(len(s.id) for s in net.species) <= 12
    assert set("".join(s.id for s in net.species)) <= set("abcdefgh")


# --- letter counts are conserved -----------------------------------------------------
def test_letter_counts_are_conserved():
    net = generate_network("smn", seed=4, n_enzymes=40, mutations=60)
    ids, R, P = net.matrices()
    S = (P - R).toarray()
    laws = net.extras["conservation"]
    assert len(laws) == 8
    for law in laws:
        vector = np.array([law["vector"].get(s, 0) for s in ids])
        assert not np.any(vector @ S)
    lengths = np.array([len(s) for s in ids])
    assert not np.any(lengths @ S)                                                 # total letter mass too


# --- book 18.3.1: the mutation operator -----------------------------------------------
def test_mutation_is_duplication_then_educt_replacement():
    net = generate_network("smn", seed=5, n_enzymes=8, mutations=50)
    genome = net.extras["genome"]
    by_id = {e["enzyme"]: e for e in genome}
    initial = net.extras["initial_metabolites"]
    compounds = list(initial)
    for e in genome:
        left, right = (part.strip() for part in e["educts"].split(" + "))
        educts = [left.replace("|", ""), right.replace("|", "")]
        if e["parent"] is None:
            assert all(x in compounds for x in educts)                             # acts on known compounds
        else:
            parent = by_id[e["parent"]]
            assert int(e["parent"][1:]) < int(e["enzyme"][1:])
            assert e["kind"] == parent["kind"]
            p_left, p_right = (part.strip() for part in parent["educts"].split(" + "))
            p_educts = [p_left.replace("|", ""), p_right.replace("|", "")]
            changed = [k for k in range(2) if educts[k] != p_educts[k]]
            assert len(changed) == 1                                               # one educt replaced
            k = changed[0]
            assert educts[k] in compounds                                          # by a compound already present
            if e["kind"] == "recombination":
                other = 1 - k
                assert [left, right][other] == [p_left, p_right][other]            # other educt keeps its cut
        for s in (*net.reactions[e["reactions"][0]].reactants, *net.reactions[e["reactions"][0]].products):
            if s not in compounds:
                compounds.append(s)
    assert [s.id for s in net.species] == compounds
    assert sum(e["parent"] is not None for e in genome) == 50


def test_no_mutations_means_initial_genome_only():
    net = generate_network("smn", seed=3)
    assert len(net.extras["genome"]) == 12 and all(e["parent"] is None for e in net.extras["genome"])
    assert len(net.extras["initial_metabolites"]) == 6
    assert set(net.initial_state) == set(net.extras["initial_metabolites"])


# --- book 18.3.1: constant amount of substance, so letter mass measures long compounds --
def test_regulated_flow_keeps_amount_constant_while_mass_grows():
    net = generate_network("smn", genome=["a + b", "ab + c", "a|bc + ab|c"], initial_metabolites=["a", "b", "c"])
    for r in net.reactions:                                          # test-only kinetics: the book gives none
        r.rate = {"law": "mass-action", "k": 1.0}
    x, _ = integrate(net, 20.0)
    total = sum(v for v in x.values())
    assert np.allclose(total, 3.0, rtol=1e-6)                         # amount of substance held constant
    mass = [letter_mass({s: v[i] for s, v in x.items()}) for i in (0, -1)]
    assert mass[0] == pytest.approx(3.0) and mass[1] > mass[0] + 0.5  # ligation raises the letter mass


# --- parameters -----------------------------------------------------------------------
def test_alphabet_size_and_max_length_bound_the_compounds():
    net = generate_network("smn", seed=8, alphabet_size=2, n_enzymes=50, mutations=100, max_length=6)
    assert set("".join(s.id for s in net.species)) <= {"a", "b"}
    assert max(len(s.id) for s in net.species) <= 6


def test_explicit_genome_is_validated():
    with pytest.raises(ValueError, match="both educts"):
        generate_network("smn", genome=["ab|c + de"])
    with pytest.raises(ValueError, match="alphabet"):
        generate_network("smn", genome=["ab + cz"])
    with pytest.raises(ValueError, match="elastic"):
        generate_network("smn", genome=["a|b + c|b"])
    with pytest.raises(ValueError, match="max_length"):
        generate_network("smn", genome=["abcdefgh + abcdefgh"])
    with pytest.raises(ValueError, match="n_metabolites"):
        generate_network("smn", alphabet_size=2, max_metabolite_length=1, n_metabolites=3)
