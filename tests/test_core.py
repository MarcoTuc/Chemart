"""Core machinery: the network record, closure, soup, parameters, tools, CLI."""

import json

import numpy as np
import pytest

import chemart
from chemart import api, catalog
from chemart.cli import main as cli_main
from chemart.expand import expand
from chemart.kinetics import k_to_c, rate_problem
from chemart.network import Network, Reaction, Species
from chemart.soup import soup


def small_network(**kwargs):
    return Network(
        species=[Species("A"), Species("B"), Species("C")],
        reactions=[
            Reaction.of(["A", "B"], ["C"], rate={"law": "mass-action", "k": 2.0}),
            Reaction.of(["A", "C"], ["A", "B"]),
        ],
        **kwargs,
    )


# --- network -----------------------------------------------------------------
def test_network_round_trip_and_text():
    net = small_network(inflow={"A": 1.0})
    assert Network.from_dict(json.loads(json.dumps(net.to_dict()))) == net
    assert net.to_text().splitlines() == ["A + B -> C  [mass-action k=2.0]", "A + C -> A + B"]


def test_network_provides_is_derived_from_content():
    net = small_network(inflow={"A": 1.0}, extras={"space": {"lattice": [10, 10]}})
    assert net.provides == [
        "catalysts", "flow", "rate-constants", "space", "stoichiometry", "topology",
    ]
    assert Network([Species("X")], []).provides == ["topology"]


def test_network_matrices():
    ids, R, P = small_network().matrices()
    assert ids == ["A", "B", "C"]
    assert R.toarray().tolist() == [[1, 1], [1, 0], [0, 1]]
    assert (P - R).toarray().tolist() == [[-1, 0], [-1, 1], [1, -1]]


@pytest.mark.parametrize(
    "bad, message",
    [
        (dict(reactions=[Reaction.of(["A"], ["Z"])]), "unknown species 'Z'"),
        (dict(reactions=[Reaction({"A": 0}, {"B": 1})]), "positive int"),
        (dict(reactions=[Reaction.of(["A"], ["B"], rate={"law": "magic"})]), "rate law"),
        (dict(reactions=[], status="finished"), "status"),
    ],
)
def test_network_rejects_inconsistent_content(bad, message):
    with pytest.raises(ValueError, match=message):
        Network(species=[Species("A"), Species("B")], **bad)


def test_rate_vocabulary():
    assert rate_problem({"law": "hill", "vmax": 1, "K": 2, "n": 4, "mode": "repression"}) is None
    assert "missing" in rate_problem({"law": "michaelis-menten", "vmax": 1})


def test_k_to_c():
    # bimolecular heterodimer: c = k / (N_A V); homodimer 2A doubles it (2! = 2)
    assert k_to_c(1.0, {"A": 1, "B": 1}, volume=1.0, avogadro=10.0) == pytest.approx(0.1)
    assert k_to_c(1.0, {"A": 2}, volume=1.0, avogadro=10.0) == pytest.approx(0.2)
    assert k_to_c(3.0, {"A": 1}, volume=5.0) == 3.0


# --- expand ------------------------------------------------------------------
def divide(a, b):
    """Prime-number chemistry: a + b -> a + b/a when a divides b."""
    if a < b and b % a == 0:
        return (a, b // a)
    return None


def test_expand_reaches_closure():
    species, reactions, status = expand(divide, [12, 2, 3])
    assert status == "complete"
    assert set(species) == {12, 2, 3, 6, 4}
    assert ((2, 12), (2, 6)) in reactions and ((3, 12), (3, 4)) in reactions
    assert ((2, 4), (2, 2)) in reactions  # multiset (2,4) -> (2,2)


def test_expand_truncates_on_budget():
    species, _, status = expand(divide, [64, 2], max_species=4)
    assert status == "truncated" and len(species) <= 4


def test_expand_tries_each_combination_once():
    calls = []

    def record(*mols):
        calls.append(mols)
        return None

    expand(record, ["a", "b", "c"], arity=2, ordered=True)
    assert len(calls) == len(set(calls)) == 9
    calls.clear()
    expand(record, ["a", "b", "c"], arity=(1, 2), ordered=False)
    assert len(calls) == len(set(calls)) == 3 + 6


# --- soup --------------------------------------------------------------------
def test_soup_is_reproducible_and_dilutes():
    def run(seed):
        return soup(divide, list(range(2, 60)), 2000, np.random.default_rng(seed), dilution="constant")

    (r1, p1), (r2, _) = run(3), run(3)
    assert r1 == r2 and len(p1) == 58
    assert all(count > 0 for _, _, count in r1)


# --- parameters --------------------------------------------------------------
def entry_with(*params):
    return catalog.Chemistry(
        id="toy", name="Toy", family="core", kind="generator", constructive=False,
        book="-", S={"repr": "x"}, R={"scheme": "x"}, A={"reactor": "ode"},
        params=list(params),
    )


N = catalog.Param(name="N", role="structural", type="int", default=4, min=1, max=10, meaning="size")
MODE = catalog.Param(name="mode", role="structural", type="enum", default="a", choices=["a", "b"], meaning="mode")
RATE = catalog.Param(name="k", role="kinetic", type="float", default=1, meaning="rate")


def test_resolve_params_fills_defaults_and_widens_floats():
    assert api.resolve_params(entry_with(N, MODE, RATE), {"N": 7}) == {"N": 7, "mode": "a", "k": 1.0}


@pytest.mark.parametrize(
    "given, message",
    [
        ({"M": 3}, "unknown parameter 'M'"),
        ({"N": 11}, "N=11 is invalid: must be <= 10"),
        ({"N": 2.5}, "expected an integer"),
        ({"N": True}, "expected an integer"),
        ({"mode": "c"}, r"one of \['a', 'b'\]"),
    ],
)
def test_resolve_params_errors_are_actionable(given, message):
    with pytest.raises(ValueError, match=message):
        api.resolve_params(entry_with(N, MODE, RATE), given)


def test_param_json_schema():
    assert N.json_schema() == {"type": "integer", "default": 4, "minimum": 1, "maximum": 10, "description": "size"}
    assert MODE.json_schema()["enum"] == ["a", "b"]


# --- public interface --------------------------------------------------------
def test_list_and_describe_cover_the_catalog():
    listing = chemart.list_chemistries()
    assert len(listing) == len(catalog.load())
    described = chemart.describe_chemistry("matrix-chemistry")
    assert described["params"]["type"] == "object" and "N" in described["params"]["properties"]
    json.dumps(described)


def test_unknown_chemistry_suggests_close_ids():
    with pytest.raises(ValueError, match="matrix-chemistry"):
        chemart.describe_chemistry("matrix-chemstry")


def test_tool_definitions_match_functions():
    names = [t["name"] for t in chemart.tool_definitions()]
    assert names == ["list_chemistries", "describe_chemistry", "generate_network"]
    assert chemart.call_tool("list_chemistries") == chemart.list_chemistries()
    assert chemart.call_tool("describe_chemistry", {"chemistry": "gamma"}) == chemart.describe_chemistry("gamma")


def test_cli_list_and_errors(capsys):
    assert cli_main(["list"]) == 0
    assert len(json.loads(capsys.readouterr().out)) == len(catalog.load())
    assert cli_main(["describe", "no-such-thing"]) == 2
    assert "unknown chemistry" in capsys.readouterr().err
