"""The contract every implemented chemistry must satisfy.

Parametrized over the catalog: a chemistry enters these tests as soon as its
module `chemart/chemistries/<id>.py` exists.
"""

import json
import os
import time

import pytest

import chemart
from chemart import catalog
from chemart.kinetics import RATE_LAWS
from chemart.network import Network

ENTRIES = {c.id: c for c in catalog.load()}
IMPLEMENTED = sorted(i for i, c in ENTRIES.items() if c.implemented)
# CHEMART_ONLY=id1,id2 restricts the contract to those chemistries (used by
# agents implementing one chemistry while others are in progress).
if os.environ.get("CHEMART_ONLY"):
    IMPLEMENTED = [i for i in os.environ["CHEMART_ONLY"].split(",") if i in ENTRIES]


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_zero_argument_generation(cid):
    start = time.perf_counter()
    net = chemart.generate_network(cid, seed=1)
    assert time.perf_counter() - start < 5, "default parameters must run in seconds"
    assert net.chemistry == cid

    restored = Network.from_dict(json.loads(json.dumps(net.to_dict())))
    assert restored == net, "network must be plain JSON data and round-trip exactly"

    ids, R, P = net.matrices()
    assert R.shape == P.shape == (len(net.species), len(net.reactions))
    assert all(r.rate is None or r.rate["law"] in RATE_LAWS for r in net.reactions)


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_same_seed_same_network(cid):
    a = chemart.generate_network(cid, seed=7)
    b = chemart.generate_network(cid, seed=7)
    assert a.to_dict() == b.to_dict()


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_catalog_claims_cover_content(cid):
    net = chemart.generate_network(cid, seed=1)
    unclaimed = set(net.provides) - set(ENTRIES[cid].provides)
    assert not unclaimed, f"network contains {unclaimed} but the catalog does not claim it"


@pytest.mark.parametrize("cid", IMPLEMENTED)
def test_rejects_bad_parameters(cid):
    with pytest.raises(ValueError, match="unknown parameter"):
        chemart.generate_network(cid, not_a_parameter=1)
    for p in ENTRIES[cid].params:
        if p.type in ("int", "float") and p.min is not None:
            with pytest.raises(ValueError, match=p.name):
                chemart.generate_network(cid, **{p.name: p.min - 1})
            break


@pytest.mark.xfail(reason="implementation waves W1-W7 in progress", strict=False)
def test_every_catalog_entry_is_implemented():
    assert set(IMPLEMENTED) == set(ENTRIES)
