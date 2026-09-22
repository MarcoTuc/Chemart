"""The contract every implemented chemistry must satisfy.

Parametrized over the catalog: a chemistry enters these tests as soon as its
module `chemart/chemistries/<id>.py` exists.
"""

import os

import pytest

import chemart
import chemart.api
from chemart import catalog, contract
from chemart.kinetics import RATE_LAWS

ENTRIES = {c.id: c for c in catalog.load()}
ALL_IMPLEMENTED = sorted(i for i, c in ENTRIES.items() if c.implemented)
IMPLEMENTED = ALL_IMPLEMENTED
# CHEMART_ONLY=id1,id2 restricts the *parametrized* contract tests to those
# chemistries (used by agents implementing one chemistry while others are in
# progress). The coverage gate below always reads the unfiltered list, so a
# single-chemistry run never makes it look as though entries went missing.
if os.environ.get("CHEMART_ONLY"):
    IMPLEMENTED = [i for i in os.environ["CHEMART_ONLY"].split(",") if i in ENTRIES]


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_contract(cid):
    """Defaults run in seconds, the network is plain JSON, same seed gives the
    same network, and the catalog claims everything the network contains."""
    entry = ENTRIES[cid]
    evolve = chemart.api.evolver_for(entry) if "evolve" in chemart.api.faces(entry) else None
    assert contract.problems(entry, chemart.api.generator_for(entry), evolve=evolve) == []


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_network_record_is_well_formed(cid):
    net = chemart.generate_network(cid, seed=1)
    assert net.chemistry == cid
    ids, R, P = net.matrices()
    assert R.shape == P.shape == (len(net.species), len(net.reactions))
    assert all(r.rate is None or r.rate["law"] in RATE_LAWS for r in net.reactions)


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_rejects_bad_parameters(cid):
    with pytest.raises(ValueError, match="unknown parameter"):
        chemart.generate_network(cid, not_a_parameter=1)
    for p in ENTRIES[cid].params:
        if p.type in ("int", "float") and p.min is not None:
            with pytest.raises(ValueError, match=p.name):
                chemart.generate_network(cid, **{p.name: p.min - 1})
            break


def test_every_catalog_entry_is_implemented():
    missing = sorted(set(ENTRIES) - set(ALL_IMPLEMENTED))
    assert not missing, f"catalogued but not implemented: {missing}"
