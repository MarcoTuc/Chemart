"""Dorin & Korb's virtual ecosystem chemistry (book 8.2.3, ref [240]).

Everything checked here comes from Dorin & Korb, "Building Virtual Ecosystems
from Artificial Chemistry", Monash TR 2007/212 (the full version of ECAL 2007,
LNCS 4648:103-112): table 1 (the bond table), the appendix (bonding rules,
electron shells, known atoms and molecules), sections 2.1-2.2 (catalysts,
energy neighbourhoods, sunlight) and sections 3.1.1-3.1.3 with figures 5-6
(photosynthesis, respiration, biosynthesis, the autotroph, the decomposer).

The paper publishes no numbers and states (section 4) that no run with all
organism types had been performed, so what is reproduced is the table itself,
the exact atom and energy bookkeeping the model claims, what each of the four
catalysts is allowed to do, and the trophic structure.
"""

from collections import Counter

import pytest

from chemart import evolve, generate_network
from chemart.chemistries import dorin_korb_ecosystem as D

ID = "dorin-korb-ecosystem"

#: Table 1 of the paper: make / break / catalysed make / catalysed break / bond energy.
TABLE_1 = {
    "A-B": ("low",      "low",      {"K": "high"},   {"EAB": "high"},             (-1, "high")),
    "C-C": ("moderate", "low",      None,            {"ECC": "high"},             (-1, "low")),
    "A-O": ("high",     "low",      None,            {"K": "high", "EO": "high"}, (+1, "low")),
    "B-O": ("high",     "low",      None,            {"K": "high", "EO": "high"}, (+1, "low")),
    "C-O": ("low",      "moderate", {"ECC": "high"}, None,                        (+1, "low")),
}

#: What each of the four catalysts of section 2.1 may do, read off table 1.
SPEC = {
    "K": {("make", "A-B"), ("break", "A-O"), ("break", "B-O")},
    "EAB": {("break", "A-B")},
    "ECC": {("break", "C-C"), ("make", "C-O")},
    "EO": {("break", "A-O"), ("break", "B-O")},
}


@pytest.fixture(scope="module")
def traj():
    return evolve(ID, seed=1)


@pytest.fixture(scope="module")
def net(traj):
    return traj.network


def catalysis(network):
    """Observed (action, bond) pairs per catalyst."""
    out = {}
    for key in network.extras["analysis"]["event_counts"]:
        action, bond, catalyst = key.split(":")
        if catalyst != "none":
            out.setdefault(catalyst, set()).add((action, bond))
    return out


def atoms_of(network):
    """Atom content per species, from the conservation laws."""
    out = {}
    for law in network.extras["conservation"]:
        for species, n in law["vector"].items():
            out.setdefault(species, {})[law["name"].removeprefix("atom ")] = n
    return out


# ---------------------------------------------------------------------------
# The chemistry as published
# ---------------------------------------------------------------------------
def test_valences_come_from_the_appendix_electron_shells():
    # Appendix "Known atoms": A:1, B:2:3, O:2:1, C:2:4:4, K:2:2, enzymes 2:2, with
    # rule 2 (a bond needs a free outer-shell electron *and* a free slot) and rule 3
    # (shells hold 2, 4, 8 electrons from inner to outer).
    assert D.SHELLS["A"] == (1,) and D.SHELLS["B"] == (2, 3) and D.SHELLS["O"] == (2, 1)
    assert D.SHELLS["C"] == (2, 4, 4) and D.SHELLS["K"] == (2, 2)
    assert D.CAPACITY == (2, 4, 8)
    # "A-B (A:2=B:2:4) - one only"; "A-O ... one only since O has run out of electrons"
    assert D.valence("A") == D.valence("B") == D.valence("O") == 1
    # "C-C-C-C-C (Internal C is C:2:4:8)": an internal C carries four bonds
    assert D.valence("C") == 4
    # the enzymes are "all two-bond maximum atoms", and K bonds twice in K-C-K
    assert D.valence("K") == D.valence("EAB") == D.valence("ECC") == D.valence("EO") == 2
    assert D.VALENCE == {"A": 1, "B": 1, "O": 1, "C": 4, "K": 2, "EAB": 2, "ECC": 2, "EO": 2}


def test_reaction_table_is_table_1(net):
    assert D.TABLE == TABLE_1
    # Table 1 caption: the bond energy is released when the bond is made and must be
    # supplied to break it, so the negative bonds are the ones that store energy.
    # A-B is "- high" (sugar), C-C is "- low" (biomass), the rest are "+ low".
    assert net.extras["energies"]["bond_energy"] == {
        "A-B": -8, "C-C": -1, "A-O": 1, "B-O": 1, "C-O": 1,
    }
    assert net.params["energy_high"] == 8 and net.params["energy_low"] == 1
    reported = net.extras["chemistry"]["table"]
    assert reported["A-B"]["energy"] == "- high" and reported["C-O"]["energy"] == "+ low"
    assert set(net.extras["chemistry"]["catalysts"]) == set(SPEC)


# ---------------------------------------------------------------------------
# The bookkeeping the model is built on
# ---------------------------------------------------------------------------
def test_atoms_are_never_created_or_destroyed(net):
    laws = {law["name"]: law["vector"] for law in net.extras["conservation"]}
    assert set(laws) == {f"atom {a}" for a in D.ATOM_ORDER}
    for name, vector in laws.items():
        for r in net.reactions:
            left = sum(vector.get(s, 0) * n for s, n in r.reactants.items())
            right = sum(vector.get(s, 0) * n for s, n in r.products.items())
            assert left == right, f"{name} is not conserved by {r.to_text()}"

    # The reactor is closed, so the inventory at the end is the one seeded at the
    # start: the free atoms, the free catalysts, the figure 5 body (16 C walls, one
    # K, one EAB, an A-O and a B-O in the vacuole) and the figure 6 body (4 C, an O
    # spacer, one ECC).
    expected = (Counter(net.params["atoms"]) + Counter(net.params["catalysts"])
                + Counter({"C": 16, "K": 1, "EAB": 1, "A": 1, "B": 1, "O": 2})
                + Counter({"C": 4, "O": 1, "ECC": 1}))
    for atom, vector in ((n.removeprefix("atom "), v) for n, v in laws.items()):
        start = sum(vector.get(s, 0) * n for s, n in net.initial_state.items())
        end = sum(vector.get(s, 0) * n for s, n in net.extras["final_state"].items())
        assert start == end == expected[atom]


def test_energy_balances_exactly_across_every_reaction(net):
    energy = net.extras["energies"]["bond_energy"]
    per_species = net.extras["energies"]["species_bonds"]

    def bonds(side):
        total = Counter()
        for species, n in side.items():
            for bond, m in per_species[species].items():
                total[bond] += n * m
        return total

    for r, event in zip(net.reactions, net.extras["events"]):
        made = bonds(r.products) - bonds(r.reactants)
        broken = bonds(r.reactants) - bonds(r.products)
        # section 2.2: one bond is made or broken per reaction, catalysts survive
        assert sum(made.values()) + sum(broken.values()) == 1
        assert dict(made or broken) == {event["bond"]: 1}
        released = (sum(energy[b] * n for b, n in made.items())
                    - sum(energy[b] * n for b, n in broken.items()))
        assert released == event["energy_released"]

    # "Energy released from a chemical bond must be used in that time step or it is
    # released in non-recoverable form": every unit is spent on a bond or dissipated.
    ledger = net.extras["energies"]["ledger"]
    assert ledger["released"] + ledger["light_spent"] == ledger["consumed"] + ledger["dissipated"]
    assert ledger["light_incident"] == ledger["light_spent"] + ledger["light_lost"]
    assert ledger["balanced"] is True
    assert ledger["released"] > 0 and ledger["consumed"] > 0 and ledger["light_spent"] > 0


# ---------------------------------------------------------------------------
# The four catalysts
# ---------------------------------------------------------------------------
def test_each_catalyst_does_exactly_its_table_1_reactions(net):
    # section 2.1: "four types of catalyst are required" - chlorophyll makes A-B and
    # breaks A-O and B-O, one enzyme breaks A-B, one decomposes C-C, one the
    # inorganic bonds. Table 1 adds that enzCC also catalyses making C-O.
    from_table = {}
    for bond, row in D.TABLE.items():
        for action, cats in (("make", row[D.MAKE_CAT]), ("break", row[D.BREAK_CAT])):
            for catalyst in cats or {}:
                from_table.setdefault(catalyst, set()).add((action, bond))
    assert from_table == SPEC
    assert catalysis(net) == SPEC
    # a catalyst is not consumed: it stands on both sides of the reaction it enables
    assert any(r.catalysts for r in net.reactions)


@pytest.mark.parametrize("seed", [2, 3])
def test_catalysts_never_step_outside_the_table(seed):
    for catalyst, actions in catalysis(generate_network(ID, seed=seed)).items():
        assert actions <= SPEC[catalyst]


# ---------------------------------------------------------------------------
# The three reactions the paper writes out
# ---------------------------------------------------------------------------
def test_photosynthesis_respiration_and_biosynthesis(net):
    energy = net.extras["energies"]["bond_energy"]
    counts = net.extras["analysis"]["event_counts"]
    high, low = net.params["energy_high"], net.params["energy_low"]

    def released(steps):
        return sum(energy[b] if a == "make" else -energy[b] for a, b in steps)

    # 3.1.1  AO + BO --(chlorophyll & sunlight)--> AB + 2 O
    photosynthesis = [("break", "A-O"), ("break", "B-O"), ("make", "A-B")]
    assert all(("K" in (D.TABLE[b][D.MAKE_CAT if a == "make" else D.BREAK_CAT] or {}))
               for a, b in photosynthesis)
    assert released(photosynthesis) == -(high + 2 * low)      # paid for by sunlight
    assert all(counts.get(f"{a}:{b}:K", 0) > 0 for a, b in photosynthesis)
    assert Counter({"A": 1, "B": 1, "K": 1}) == Counter(
        next(r for r in net.reactions if r.products.get("AB")).reactants)

    # 3.1.2  O + AB --(enzyme)--> A + BO + energy
    respiration = [("break", "A-B"), ("make", "B-O")]
    assert released(respiration) == high + low > 0            # "releases energy"
    assert counts.get("break:A-B:EAB", 0) > 0 and counts.get("make:B-O:none", 0) > 0

    # 3.1.3  C + C --(energy)--> C2
    assert released([("make", "C-C")]) == -low                # "consuming energy"
    assert counts.get("make:C-C:none", 0) > 0


# ---------------------------------------------------------------------------
# The trophic structure
# ---------------------------------------------------------------------------
def test_no_sugar_is_made_without_sunlight():
    # "the only way energy can be stored is in complex molecules", and only sunlight
    # can pay the high A-B bond energy, so in the dark no sugar bond ever exists
    run = evolve(ID, seed=1, light_amplitude=0.0)
    dark = run.network
    analysis = dark.extras["analysis"]
    assert analysis["trophic"]["sugar_made"] == 0
    assert max(run.series("sugar_bonds")) == 0
    assert dark.extras["energies"]["ledger"]["light_incident"] == 0
    assert generate_network(ID, seed=1).extras["analysis"]["trophic"]["sugar_made"] > 0


def test_trophic_structure_emerges(traj, net):
    # the paper's claim: organisms "naturally fall into trophic levels, generate
    # energy from chemical bonds and transform material elements in the process"
    trophic = net.extras["analysis"]["trophic"]
    assert trophic["sugar_made"] >= 5            # producers fix sunlight in A-B bonds
    assert trophic["sugar_respired"] >= 2        # the sugar enzyme releases it again
    assert trophic["biomass_built"] >= 5         # that energy pays for C-C biosynthesis
    assert trophic["biomass_decomposed"] >= 1    # the organic decomposer breaks C-C
    assert trophic["inorganic_split"] >= 20      # A-O and B-O go back to free atoms

    # decomposition splits biomass, returning atoms to the pool: the decomposer
    # stands on both sides (it is a catalyst), the C-C molecule it attacks carries at
    # least two carbons, and it leaves the reaction as more, smaller molecules
    content = atoms_of(net)
    breaks = [r for r, e in zip(net.reactions, net.extras["events"])
              if e["bond"] == "C-C" and e["action"] == "break" and "ECC" in e["catalysts"]]
    assert breaks
    for r in breaks:
        assert sum(r.products.values()) > sum(r.reactants.values())
        assert max(content[s].get("C", 0) for s in r.reactants) >= 2
    biomass = traj.series("biomass_bonds")
    assert any(b < a for a, b in zip(biomass, biomass[1:]))


def test_without_the_organic_decomposer_nothing_breaks_C_C_bonds():
    # remove every ECC atom: the free ones and the one the figure 6 body carries
    none = generate_network(ID, seed=1, structures="photoautotroph",
                            catalysts={"K": 4, "EAB": 4, "ECC": 0, "EO": 4})
    analysis = none.extras["analysis"]
    assert analysis["trophic"]["biomass_decomposed"] == 0
    assert not any(k.endswith(":ECC") for k in analysis["event_counts"])
    assert "ECC" not in catalysis(none)


def test_figure_5_autotroph_makes_and_respires_sugar_on_its_own():
    # nothing on the grid but the figure 5 body: a C box with chlorophyll on one
    # inner wall, the sugar enzyme on the opposite one, and A-O and B-O trapped in
    # the vacuole. It photosynthesises and then respires what it made.
    run = evolve(ID, seed=5, atoms={}, catalysts={}, structures="photoautotroph")
    alone = run.network
    counts = alone.extras["analysis"]["event_counts"]
    assert counts.get("break:A-O:K", 0) > 0 and counts.get("break:B-O:K", 0) > 0
    assert counts.get("make:A-B:K", 0) > 0 and counts.get("break:A-B:EAB", 0) > 0
    assert max(run.series("sugar_bonds")) > 0

    # the body itself never falls apart: its bonds are anchors that do not react
    content = atoms_of(alone)
    bodies = [s for s, a in content.items()
              if a.get("K") and a.get("EAB") and s in alone.extras["final_state"]]
    assert len(bodies) == 1 and content[bodies[0]]["C"] >= 16


def test_a_frame_per_step(traj, net):
    # frame 0 is the seeded grid; each later frame is one movement and one reaction phase
    assert traj.clock == "steps" and traj.times() == [float(t) for t in range(net.params["steps"] + 1)]
    assert traj.frames[0].fired == [] and traj.frames[0].state == net.initial_state
    assert traj.frames[-1].state == {s: float(n) for s, n in net.extras["final_state"].items()}
    # the reactions of each step turn one frame's molecules into the next one's
    for before, after in zip(traj.frames, traj.frames[1:]):
        state = Counter(before.state)
        for lhs, rhs, n in after.fired:
            state.subtract({s: n * k for s, k in Counter(lhs).items()})
            state.update({s: n * k for s, k in Counter(rhs).items()})
        assert +state == Counter(after.state)
    assert set(traj.frames[0].observables) == {"sugar_bonds", "biomass_bonds", "inorganic_bonds", "free_atoms"}


# ---------------------------------------------------------------------------
# Reproducibility and rejected parameters
# ---------------------------------------------------------------------------
def test_the_seed_reproduces_the_run():
    a = generate_network(ID, seed=11)
    assert a.to_dict() == generate_network(ID, seed=11).to_dict()
    assert a.to_dict() != generate_network(ID, seed=12).to_dict()


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError, match="energy_high"):
        generate_network(ID, energy_high=2, energy_low=2)
    with pytest.raises(ValueError, match="p_low"):
        generate_network(ID, p_low=0.5, p_moderate=0.2)
    with pytest.raises(ValueError, match="unknown atom"):
        generate_network(ID, atoms={"Z": 1})
    with pytest.raises(ValueError, match="do not fit"):
        generate_network(ID, width=6, height=6, structures="none", atoms={"A": 400})
    with pytest.raises(ValueError, match="no free room"):
        generate_network(ID, width=6, height=6)


@pytest.mark.slow
def test_the_books_stay_exact_on_a_larger_run():
    big = generate_network(ID, seed=4, width=32, height=32, steps=400,
                           atoms={"A": 90, "B": 90, "O": 110, "C": 80},
                           catalysts={"K": 10, "EAB": 10, "ECC": 10, "EO": 10})
    laws = [law["vector"] for law in big.extras["conservation"]]
    for vector in laws:
        for r in big.reactions:
            assert (sum(vector.get(s, 0) * n for s, n in r.reactants.items())
                    == sum(vector.get(s, 0) * n for s, n in r.products.items()))
    ledger = big.extras["energies"]["ledger"]
    assert ledger["released"] + ledger["light_spent"] == ledger["consumed"] + ledger["dissipated"]
    assert catalysis(big) == SPEC
    trophic = big.extras["analysis"]["trophic"]
    assert trophic["sugar_made"] > 20 and trophic["sugar_respired"] > 10
    assert trophic["biomass_built"] > 10 and trophic["biomass_decomposed"] > 5
