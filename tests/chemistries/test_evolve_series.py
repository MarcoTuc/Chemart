"""The EVOLVE series of virtual ecosystems (book 8.2.3). Catalog id: evolve-series.

None of the six EVOLVE-series publications is openly available (all closed
access, no repository copy, the Internet Archive offline), so what is checked
here is (a) the book's paragraph, sentence by sentence - the weighted critical
section, the 2D world, energy and matter conservation, and the scavenger niche
that emerges, goes extinct and reappears - and (b) the two findings the
published abstracts state: EVOLVE II's evolvable mutational sensitivity
("the magnitude of phenotypic change resulting from mutation is itself a
property of the gene", and resistance to phenotypic change evolving in slowly
varying environments) and the coexistence of lineages with different survival
strategies.

Where the model does *not* do what one might hope, the test says so: the enzyme
match does not improve over a default run, so what is asserted is the weaker
true statement - selection on the match is real (a badly matched founder leaves
no descendants) while drift erodes it from near-optimal founders.
"""

from types import SimpleNamespace

import pytest

from chemart import describe_chemistry, generate_network
from chemart.chemistries import evolve_series as E

ID = "evolve-series"
DEFAULTS = {k: v.get("default")
            for k, v in describe_chemistry(ID)["params"]["properties"].items()}


def chemistry(**over) -> E.Chemistry:
    return E.Chemistry(SimpleNamespace(**{**DEFAULTS, **over}))


@pytest.fixture(scope="module")
def net():
    return generate_network(ID, seed=1)


def matter_vector(network) -> dict[str, int]:
    laws = {law["name"]: law["vector"] for law in network.extras["conservation"]}
    return laws["matter"]


def total(vector, state) -> int:
    return sum(vector[s] * int(n) for s, n in state.items())


# ---------------------------------------------------------------------------
# "Enzymes and other proteins are represented as strings with critical sections
#  that determine their shape. Their function is then fetched from a table by a
#  matching mechanism that weights the critical section more heavily than the
#  rest of the string." (book 8.2.3)
# ---------------------------------------------------------------------------
def test_the_function_table_is_deterministic():
    chem = chemistry()
    # ref_i[k] = (i + k (i + 1)) mod alphabet, for alphabet 4 and length 8
    assert E.reference(0, 8, 4) == (0, 1, 2, 3, 0, 1, 2, 3)
    assert E.reference(2, 8, 4) == (2, 1, 0, 3, 2, 1, 0, 3)
    assert chem.spell("light") == "abcdabcd"
    assert chem.spell("fix") == "bdbdbdbd"
    assert chem.spell("respire") == "cbadcbad"
    assert chem.spell("replicate") == "dddddddd"
    # the critical section is the centred window of the protein string
    assert E.critical_slice(8, 3) == slice(2, 5)
    assert (chem.critical.start, chem.critical.stop) == (2, 5)
    assert chem.max_score == pytest.approx(3.0 * 3 + 5)


def test_the_critical_section_is_weighted_more_heavily_than_the_rest():
    chem = chemistry()
    ref = chem.table["light"]
    protein = tuple(ord(c) - ord("a") for c in "abcdabcd")
    assert chem.score(protein, ref) == pytest.approx(1.0)      # a perfect key

    # matches the reference in the critical section only (3 of 8 symbols)
    inside = tuple(ord(c) - ord("a") for c in "bccdacda")
    # matches everywhere *but* the critical section (5 of 8 symbols)
    outside = tuple(ord(c) - ord("a") for c in "abaabbcd")
    assert chem.score(inside, ref) == pytest.approx(9 / 14)
    assert chem.score(outside, ref) == pytest.approx(5 / 14)
    # the book's point: the critical section decides the shape, so the protein
    # matching 3 critical symbols beats the one matching 5 peripheral ones
    assert chem.score(inside, ref) > chem.score(outside, ref)

    # and it is the weighting that does it: with an unweighted match the
    # ordering reverses, because then only the number of matches counts
    plain = chemistry(critical_weight=1.0)
    assert plain.score(inside, ref) == pytest.approx(3 / 8)
    assert plain.score(outside, ref) == pytest.approx(5 / 8)
    assert plain.score(inside, ref) < plain.score(outside, ref)


def test_a_protein_below_the_threshold_has_no_function():
    chem = chemistry()
    assert chem.function(tuple(ord(c) - ord("a") for c in "abcdabcd")) == ("light", 1.0)
    # best score 4/14 = 0.286, under the 0.55 threshold: no function at all
    name, quality = chem.function(tuple(ord(c) - ord("a") for c in "aaaaaaaa"))
    assert name is None and quality == pytest.approx(4 / 14)


def test_the_sensitivity_locus_is_read_from_the_gene():
    # EVOLVE II (Conrad & Strizich 1985, abstract): "the magnitude of phenotypic
    # change resulting from mutation is itself a property of the gene".
    chem = chemistry()
    assert chem.gene_length == 9 and chem.genome_length == 54
    assert chem.sensitivity("a" + "abcdabcd") == 1
    assert chem.sensitivity("c" + "abcdabcd") == 3


# ---------------------------------------------------------------------------
# "virtual organisms are also subject to energy and matter conservation"
# ---------------------------------------------------------------------------
def test_matter_is_conserved_by_every_event(net):
    vector = matter_vector(net)
    assert vector[E.MINERAL] == vector[E.ORGANIC] == 1
    assert vector[net.species[0].id] == net.params["body"]
    for r in net.reactions:
        left = sum(vector[s] * n for s, n in r.reactants.items())
        right = sum(vector[s] * n for s, n in r.products.items())
        assert left == right, f"matter is not conserved by {r.to_text()}"


def test_the_world_inventory_never_changes(net):
    vector = matter_vector(net)
    p = net.params
    # every cell starts with initial_mineral and initial_organic units, and each
    # founder holds one body's worth
    expected = (p["width"] * p["height"] * (p["initial_mineral"] + p["initial_organic"])
                + p["n_founders"] * p["body"])
    start = total(vector, net.initial_state)
    end = total(vector, net.extras["final_state"])
    pools = net.extras["analysis"]["history"]
    assert start == expected == net.extras["analysis"]["matter_total"]
    # the survivors' bodies plus what is left in the two pools is the same total
    assert end + pools["mineral"][-1] + pools["organic"][-1] == expected
    assert start == end + pools["mineral"][-1] + pools["organic"][-1]


def test_energy_is_accounted_for_to_the_last_unit(net):
    ledger = net.extras["energies"]["ledger"]
    # light is the only input, and what is not harvested is lost
    assert ledger["light_incident"] == ledger["light_harvested"] + ledger["light_lost"]
    # what the organisms took in is what they spent, banked or still hold
    assert (ledger["founders_endowment"] + ledger["light_harvested"] + ledger["respired"]
            == ledger["fix_spent"] + ledger["maintenance_spent"]
            + ledger["repro_dissipated"] + ledger["death_deposited"]
            + ledger["organisms_energy"])
    # and the energy banked in organic matter is what was fixed into it plus what
    # the dead left in it, less what was respired out of it and dissipated
    assert (ledger["fix_spent"] + ledger["death_deposited"]
            == ledger["respired"] + ledger["dissipated"] + ledger["organic_energy"])
    assert ledger["organisms_balanced"] is True
    assert ledger["organic_balanced"] is True
    assert ledger["light_harvested"] > 0 and ledger["respired"] > 0


def test_no_light_no_ecosystem():
    # energy enters only as light; with the sun off the founders burn the energy
    # they were endowed with, never reproduce, and starve
    dark = generate_network(ID, seed=1, light_input=0)
    counts = dark.extras["analysis"]["event_counts"]
    ledger = dark.extras["energies"]["ledger"]
    assert ledger["light_incident"] == 0 and ledger["light_harvested"] == 0
    assert "reproduction" not in counts and "conjugation" not in counts
    assert dark.extras["analysis"]["survivors"] == 0
    assert max(dark.extras["analysis"]["history"]["organisms"]) == dark.params["n_founders"]
    # matter still balances exactly in a dead world
    vector = matter_vector(dark)
    for r in dark.reactions:
        assert (sum(vector[s] * n for s, n in r.reactants.items())
                == sum(vector[s] * n for s, n in r.products.items()))


# ---------------------------------------------------------------------------
# The trophic structure and the events it produces
# ---------------------------------------------------------------------------
def test_every_event_is_catalysed_by_the_organism_that_performs_it(net):
    body = net.params["body"]
    kinds = {}
    for r, event in zip(net.reactions, net.extras["events"]):
        kinds.setdefault(event["kind"], []).append(r)
    assert {"fix", "reproduction", "death"} <= set(kinds)

    for r in kinds["fix"]:                      # G + mineral -> G + organic
        assert r.reactants[E.MINERAL] == 1 and r.products[E.ORGANIC] == 1
        assert r.catalysts and E.ORGANIC not in r.reactants
    for r in kinds.get("respire", []):          # G + organic -> G + mineral
        assert r.reactants[E.ORGANIC] == 1 and r.products[E.MINERAL] == 1
        assert r.catalysts
    for r in kinds["reproduction"]:             # G + b organic -> G + O
        assert r.reactants[E.ORGANIC] == body
        assert sum(n for s, n in r.products.items() if s not in E.POOLS) == 2
        assert r.catalysts
    for r in kinds["death"]:                    # G -> b organic
        assert r.products == {E.ORGANIC: body} and len(r.reactants) == 1
    for r in kinds.get("conjugation", []):      # G + H + b organic -> G + H + O
        # two parents in, two parents plus the offspring out; the donor may be a
        # clone of the parent, in which case both are the same species
        assert sum(n for s, n in r.reactants.items() if s not in E.POOLS) == 2
        assert sum(n for s, n in r.products.items() if s not in E.POOLS) == 3
        assert r.reactants[E.ORGANIC] == body and r.catalysts


def test_the_founders_are_pure_autotrophs(net):
    # "Starting with a population of only autotrophs..." - scavenging is the one
    # table entry the founders lack, so it can only arise by mutation
    phenotypes = net.extras["phenotypes"]
    founders = [s for s in net.initial_state if s not in E.POOLS]
    assert len(founders) >= 1
    for gid in founders:
        assert phenotypes[gid]["guild"] == "autotroph"
        assert phenotypes[gid]["efficiency"]["fix"] > 0.0
        assert phenotypes[gid]["efficiency"]["respire"] == 0.0
    assert net.extras["analysis"]["history"]["scavengers"][0] == 0


def test_the_scavenger_niche_emerges(net):
    # "...scavenger populations emerge, go extinct, and reappear, showing the
    # emergence of an ecological niche aimed at decomposing the metabolic
    # remains of autotrophs."
    analysis = net.extras["analysis"]
    history = analysis["history"]
    assert max(history["scavengers"]) > 0, "no scavenger ever appeared"
    assert analysis["event_counts"]["respire"] > 0
    assert analysis["scavenger_episodes"], "the guild never established itself"
    # the only route from organic matter back to mineral matter is scavenging,
    # and the only routes into organic matter are fixing and death
    trophic = net.extras["chemistry"]["trophic"]
    assert set(trophic) == {"fix", "respire", "reproduce", "death"}
    # lineages with different survival strategies coexist (EVOLVE II abstract)
    guilds = analysis["final_guilds"]
    assert guilds.get("autotroph", 0) > 0
    assert guilds.get("mixotroph", 0) + guilds.get("scavenger", 0) > 0


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_resistance_to_phenotypic_change_evolves(seed):
    # EVOLVE II (Conrad & Strizich 1985, abstract): "Organisms developed a
    # resistance to phenotypic change in response to mutation in slowly varying
    # environments."  The default environment is constant (light_period 0).
    analysis = generate_network(ID, seed=seed).extras["analysis"]
    assert analysis["founder_sensitivity"] == pytest.approx(DEFAULTS["initial_sensitivity"])
    assert analysis["final_sensitivity"] < analysis["founder_sensitivity"]


# ---------------------------------------------------------------------------
# Selection on the lock-and-key match
# ---------------------------------------------------------------------------
def genomes(chem):
    """A perfectly matched autotroph, and one whose critical sections are dead."""
    locus = chr(ord("a") + DEFAULTS["initial_sensitivity"] - 1)
    genes = [locus + chem.spell(n) for n in E.FUNCTIONS if n != "respire"]
    genes += [locus + "aaaaaaaa"] * (DEFAULTS["n_genes"] - len(genes))
    good = "".join(genes)
    symbols = list(good)
    for gene in (0, 1):                       # wreck the light and fix keys
        for k in range(chem.critical.start, chem.critical.stop):
            i = gene * chem.gene_length + 1 + k
            symbols[i] = chr(ord("a") + (ord(symbols[i]) - ord("a") + 1) % chem.alphabet)
    return good, "".join(symbols)


def test_a_dead_key_is_a_dead_enzyme():
    chem = chemistry()
    good, bad = genomes(chem)
    assert chem.function(tuple(ord(c) - ord("a") for c in "aaaaaaaa"))[0] is None
    assert chem.phenotype(good)["efficiency"] == {
        "light": 1.0, "fix": 1.0, "respire": 0.0, "replicate": 1.0}
    dead = chem.phenotype(bad)["efficiency"]
    assert dead["light"] == 0.0 and dead["fix"] == 0.0
    assert dead["replicate"] == 1.0, "only the light and fix genes were touched"


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_selection_acts_on_the_match(seed):
    # Seeded half and half, the badly matched founder leaves no descendants and
    # the population's mean light+fix efficiency climbs well above the 1.0 it
    # started from (the average of the perfect 2.0 and the dead 0.0).
    chem = chemistry()
    good, bad = genomes(chem)
    result = generate_network(ID, seed=seed, founder_genomes=[good, bad])
    history = result.extras["analysis"]["history"]
    start = history["mean_light_efficiency"][0] + history["mean_fix_efficiency"][0]
    end = history["mean_light_efficiency"][-1] + history["mean_fix_efficiency"][-1]
    assert start == pytest.approx(1.0)
    assert end > 1.4
    assert result.extras["final_state"].get(E.genotype_id(bad), 0) == 0


def test_the_match_is_not_improved_by_a_default_run(net):
    # Honest negative result: from the near-optimal default founders, mutation at
    # the default rate erodes the enzyme match faster than selection repairs it.
    # Nothing in the sources fixes a mutation rate, so this is reported, not tuned
    # away; at mutation_rate 0.02 the match is held near its initial value instead.
    analysis = net.extras["analysis"]
    assert analysis["founder_fix_efficiency"] > analysis["final_fix_efficiency"]
    slow_mutation = generate_network(ID, seed=1, mutation_rate=0.02).extras["analysis"]
    assert (slow_mutation["final_fix_efficiency"] - slow_mutation["founder_fix_efficiency"]
            > analysis["final_fix_efficiency"] - analysis["founder_fix_efficiency"])


# ---------------------------------------------------------------------------
# Reproducibility and rejected parameters
# ---------------------------------------------------------------------------
def test_the_seed_reproduces_the_run():
    a = generate_network(ID, seed=11)
    assert a.to_dict() == generate_network(ID, seed=11).to_dict()
    assert a.to_dict() != generate_network(ID, seed=12).to_dict()


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError, match="critical_length"):
        generate_network(ID, critical_length=9, protein_length=8)
    with pytest.raises(ValueError, match="critical_weight"):
        generate_network(ID, critical_weight=0.5)
    with pytest.raises(ValueError, match="n_genes"):
        generate_network(ID, n_genes=3)
    with pytest.raises(ValueError, match="initial_sensitivity"):
        generate_network(ID, initial_sensitivity=5, alphabet=4)
    with pytest.raises(ValueError, match="n_founders"):
        generate_network(ID, width=2, height=2, n_founders=5)
    with pytest.raises(ValueError, match="founder genome"):
        generate_network(ID, founder_genomes=["abc"])
    with pytest.raises(ValueError, match="alphabet"):
        generate_network(ID, founder_genomes=["z" * 54])


# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_the_niche_comes_and_goes_over_a_long_run():
    # The book's full claim needs more than the default 120 steps: over 600 steps
    # the scavengers appear, die out and come back, so their presence breaks into
    # several separate episodes.
    long = generate_network(ID, seed=4, steps=600)
    analysis = long.extras["analysis"]
    episodes = analysis["scavenger_episodes"]
    assert len(episodes) >= 2, f"the guild never came back: {episodes}"
    scavengers = analysis["history"]["scavengers"]
    for before, after in zip(episodes, episodes[1:]):
        assert after["from"] > before["to"] + 1
        assert scavengers[before["to"] + 1] == 0    # it really did go extinct
    assert max(scavengers) > 0

    # matter and energy still balance exactly over the longer run
    vector = matter_vector(long)
    for r in long.reactions:
        assert (sum(vector[s] * n for s, n in r.reactants.items())
                == sum(vector[s] * n for s, n in r.products.items()))
    ledger = long.extras["energies"]["ledger"]
    assert ledger["organisms_balanced"] is True and ledger["organic_balanced"] is True
