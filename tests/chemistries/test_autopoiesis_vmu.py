"""Autopoiesis of Varela, Maturana & Uribe (1974), in McMullin's SCL reconstruction.

Book 6.1.5 and figure 6.3 for the reactions and the phenomena; the rule set and
every published parameter value are those of the original SCL system as
documented by Von Kamp (2002); chain-based bond inhibition is the rule the 1974
paper omitted, restored by McMullin & Varela (1997). See the catalog sources.
"""

import numpy as np
import pytest

from chemart import generate_network
from chemart.chemistries import autopoiesis_vmu as A

ID = "autopoiesis-vmu"
NO_DECAY = {"disintegration_probability": 0.0}


def ruptured_cell(seed, inhibition=True, warmup=20, size=13):
    """A formed cell whose membrane has just lost one of its constituents."""
    world = A.World(np.random.default_rng(seed), size, size,
                    disintegration_probability=0.0, bond_inhibition=inhibition,
                    record=False)
    centre = world.seed_catalysts(1)[0]
    ring = world.seed_ring(centre)
    for _ in range(warmup):
        world.step()
    assert world.enclosed()[1], "the cell must be formed before its membrane decays"
    world.disintegrate(ring[0])                     # [O] -> O + O, bonds dropped
    return world


def encloses_membrane(world):
    """Is a catalyst enclosed by a closed chain of six or more links (a membrane)?"""
    cycles, interiors = world.enclosed()
    return bool(interiors) and any(len(c) >= 6 for c in cycles)


def repair_time(world, steps=60):
    for t in range(steps):
        world.step()
        if encloses_membrane(world):
            return t + 1
    return None


# --- the reactions the chemistry defines -----------------------------------------
def test_defined_network_is_the_books_three_reactions():
    # book 6.1.5 and figure 6.3: * + 2 O -> * + [O], links bond to links, [O] decays
    net = generate_network(ID, mode="reactions")
    assert net.status == "complete"
    text = {r.to_text() for r in net.reactions}
    assert "C + 2 S -> C + L0" in text                     # production, the catalyst is a catalyst
    assert {"2 L0 -> 2 L1", "L0 + L1 -> L1 + L2", "2 L1 -> 2 L2"} <= text   # bonding
    assert {"L0 -> 2 S", "L1 -> 2 S", "L2 -> 2 S"} <= text                  # disintegration
    assert {"L2 + S -> L2S", "L2S -> L2 + S"} <= text      # what makes a membrane permeable
    (production,) = [r for r in net.reactions if r.reactants == {"C": 1, "S": 2}]
    assert production.catalysts == {"C": 1}
    # no rate law is published for the model
    assert all(r.rate is None for r in net.reactions)


def test_observed_run_fires_those_reactions_with_counts():
    net = generate_network(ID, seed=1)
    assert net.status == "observed"
    events = net.extras["analysis"]["events"]
    (production,) = [r for r in net.reactions if r.reactants == {"C": 1, "S": 2}]
    assert production.products == {"C": 1, "L0": 1}
    assert production.count == events["production"] > 0
    assert all(r.count and r.count > 0 for r in net.reactions)
    assert events["bond"] > 0 and events["disintegration"] > 0
    # a link disintegrates into substrate, and links pass substrate in and out
    assert any(set(r.products) == {"S"} for r in net.reactions)
    assert events["absorption"] > 0 and events["emission"] > 0
    space = net.extras["space"]
    assert space["shape"] == [30, 30] and len(space["grid"]) == 30
    assert all(len(row) == 30 for row in space["grid"])
    assert net.initial_state["C"] == 1.0 and net.initial_state["S"] == 899.0


# --- the published phenomena -------------------------------------------------------
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_membrane_is_permeable_to_substrate_but_not_to_links_or_catalysts(seed):
    """Book 6.1.5: closed membranes are 'nevertheless permeable to substrate molecules'.

    Substrate crosses by being absorbed and re-emitted on the other side; bonded
    links are immobile, so nothing else crosses.
    """
    net = generate_network(ID, seed=seed, initial="cell", steps=120, **NO_DECAY)
    membrane = net.extras["analysis"]["membrane"]
    assert membrane["membranes"] == 1 and membrane["chain_lengths"] == [12]
    assert membrane["interior_sizes"] == [9]
    assert membrane["steps_enclosed"] == 121            # closed at every step
    crossings = net.extras["analysis"]["permeability"]
    assert crossings["substrate_crossings"] > 0
    assert crossings["link_crossings"] == 0
    assert crossings["catalyst_crossings"] == 0


def test_boundary_is_repaired_after_one_of_its_links_decays():
    """The headline result of Varela et al. (1974): the membrane regenerates.

    A free link produced inside the cell moves into the gap and bonds to the two
    chain ends, closing the membrane again.
    """
    world = ruptured_cell(seed=0)
    assert not encloses_membrane(world), "the decayed link must open the membrane"
    assert repair_time(world) is not None

    # and over a whole run with links decaying at the published rate, the cell is
    # repeatedly ruptured and repaired before it finally comes apart
    net = generate_network(ID, seed=0, initial="cell", steps=120)
    membrane = net.extras["analysis"]["membrane"]
    assert membrane["ruptures"] > 0 and membrane["repairs"] > 0


def test_a_closed_boundary_forms_spontaneously_around_the_catalyst():
    """Book figure 6.3: from one catalyst in substrate, links bond into a boundary."""
    world = A.World(np.random.default_rng(3), 14, 14, record=False)
    world.seed_catalysts(1)
    first, longest = None, 0
    for t in range(2000):
        world.step()
        cycles, interiors = world.enclosed()
        if cycles:
            longest = max(longest, max(len(c) for c in cycles))
        if first is None and encloses_membrane(world):
            first = t
    assert first is not None, "no membrane ever closed around the catalyst"
    assert longest >= 6


# --- the rule the 1974 publication omitted ------------------------------------------
def test_chain_based_bond_inhibition_keeps_the_links_inside_a_cell_free():
    """McMullin & Varela (1997): the interaction missing from the published model.

    Without it the free links inside a cell bond to each other and stop moving,
    so none is available to repair the membrane.
    """
    free, bonded = [], []
    for seed in range(6):
        with_rule = generate_network(ID, seed=seed, initial="cell", steps=80, **NO_DECAY)
        without = generate_network(ID, seed=seed, initial="cell", steps=80,
                                   bond_inhibition=False, **NO_DECAY)
        free.append((with_rule.extras["analysis"]["per_step"]["free_links"][-1],
                     with_rule.extras["analysis"]["events"].get("bond", 0)))
        bonded.append((without.extras["analysis"]["per_step"]["free_links"][-1],
                       without.extras["analysis"]["events"].get("bond", 0)))
    # with the rule the interior links stay free and hardly ever bond
    assert all(links >= 5 and bonds <= 2 for links, bonds in free), free
    # without it they bond into clusters and almost no free link is left
    assert all(links <= 3 and bonds >= 4 for links, bonds in bonded), bonded


def test_without_chain_based_bond_inhibition_the_membrane_is_not_repaired():
    seeds = range(10)
    with_rule = [repair_time(ruptured_cell(s)) for s in seeds]
    without = [repair_time(ruptured_cell(s, inhibition=False)) for s in seeds]
    assert sum(t is not None for t in with_rule) >= 6, with_rule
    assert sum(t is not None for t in without) <= 3, without


# --- mechanics ------------------------------------------------------------------------
def test_a_closed_chain_encloses_its_interior_but_a_cluster_encloses_nothing():
    world = A.World(np.random.default_rng(0), 15, 15)
    centre = world.seed_catalysts(1)[0]
    world.seed_ring(centre)                       # Von Kamp (2002), figure 2
    cycles, interiors = world.enclosed()
    assert [len(c) for c in cycles] == [12]
    assert [len(r) for r in interiors.values()] == [9] and centre in interiors[centre]
    # three mutually bonded links are a cluster: a closed chain that holds nothing
    other = A.World(np.random.default_rng(0), 15, 15)
    other.seed_catalysts(1)
    sites = [other.index(0, 0), other.index(1, 0), other.index(0, 1)]
    for s in sites:
        other.put(s, A.LINK)
    other.bond(sites[0], sites[1])
    other.bond(sites[1], sites[2])
    other.bond(sites[2], sites[0])
    cycles, interiors = other.enclosed()
    assert [len(c) for c in cycles] == [3] and interiors == {}


def test_reproducible_with_a_seed():
    kw = dict(width=12, height=12, steps=40)
    assert generate_network(ID, seed=5, **kw).to_dict() == generate_network(ID, seed=5, **kw).to_dict()
    assert generate_network(ID, seed=5, **kw).to_dict() != generate_network(ID, seed=6, **kw).to_dict()


def test_parameters_are_checked():
    with pytest.raises(ValueError, match="unknown particle type"):
        generate_network(ID, mobility={"proton": 0.5})
    with pytest.raises(ValueError, match=r"mobility\['link'\]"):
        generate_network(ID, mobility={"link": 2.0})
    with pytest.raises(ValueError, match="exceeds"):
        generate_network(ID, width=5, height=5, n_catalysts=30, steps=0)
