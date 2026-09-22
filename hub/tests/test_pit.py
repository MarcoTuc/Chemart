"""The simulation pit: local only, guarded, and streaming runs."""

import json

import pytest
from fastapi.testclient import TestClient

from chemart_hub.pit import create_pit_app

POST = {"X-Chemart-Pit": "1"}


@pytest.fixture(scope="module")
def pit():
    return TestClient(create_pit_app(), base_url="http://127.0.0.1:8765")


def lines(response):
    return [json.loads(line) for line in response.text.splitlines() if line.strip()]


def test_serves_the_page_with_a_strict_policy(pit):
    r = pit.get("/")
    assert r.status_code == 200 and "Chemart pit" in r.text
    assert "script-src 'self'" in r.headers["content-security-policy"]
    assert pit.get("/static/pit/pit.js").status_code == 200


def test_refuses_other_hosts_and_unmarked_posts(pit):
    assert pit.get("/", headers={"Host": "evil.example"}).status_code == 403     # DNS rebinding
    r = pit.post("/api/run", json={"chemistry": "brusselator"})
    assert r.status_code == 403 and "X-Chemart-Pit" in r.json()["error"].replace("x-chemart-pit", "X-Chemart-Pit")


def test_lists_chemistries_with_their_faces(pit):
    rows = {r["id"]: r for r in pit.get("/api/chemistries").json()}
    assert rows["brusselator"]["type"] == "given" and rows["brusselator"]["faces"] == ["generate"]
    assert rows["alchemy"]["faces"] == ["generate", "evolve"]
    info = pit.get("/api/chemistry/alchemy").json()
    assert "collisions" in info["evolve_params"]["properties"] and info["measures"]


def test_streams_an_ode_run_and_keeps_it(pit):
    r = pit.post("/api/run", headers=POST, json={"chemistry": "brusselator", "method": "ode",
                                                 "t_end": 10, "points": 11, "seed": 1})
    msgs = lines(r)
    assert [m["type"] for m in msgs][:2] == ["network", "frame"] and msgs[-1]["type"] == "done"
    assert sum(m["type"] == "frame" for m in msgs) == 11
    assert msgs[0]["measures"]["values"]["deficiency"] == 0
    traj = pit.get(f"/api/runs/{msgs[-1]['run']}.json").json()
    assert traj["method"] == "ode" and len(traj["frames"]) == 11


def test_streams_an_evolving_gas(pit):
    r = pit.post("/api/run", headers=POST, json={"chemistry": "alchemy", "method": "evolve", "seed": 1,
                                                 "measures": ["richness", "n_reactions"], "window": 2})
    msgs = lines(r)
    frames = [m for m in msgs if m["type"] == "frame"]
    assert frames[0]["measures"]["richness"] == 100 and frames[-1]["measures"]["richness"] == 32
    assert "n_reactions" in frames[1]["measures"]
    assert msgs[-2]["type"] == "network" and msgs[-1]["type"] == "done"


def test_rates_from_a_distribution_and_errors_in_the_stream(pit):
    body = {"chemistry": "kauffman-autocatalytic-sets", "method": "ode", "t_end": 1, "points": 3,
            "seed": 0, "x0": 1.0}
    assert "no rate constants" in lines(pit.post("/api/run", headers=POST, json=body))[-1]["message"]
    body["rates"] = {"dist": "lognormal", "mean": 0, "sigma": 1}
    assert lines(pit.post("/api/run", headers=POST, json=body))[-1]["type"] == "done"


def test_measure_and_sweep(pit):
    r = pit.post("/api/measure", headers=POST, json={"chemistry": "michaelis-menten", "names": ["deficiency"]})
    assert r.json()["values"] == {"deficiency": 0}
    msgs = lines(pit.post("/api/sweep", headers=POST, json={
        "chemistry": "random-catalytic-networks", "param": "n", "values": [5, 10], "seeds": [0, 1],
        "names": ["n_species"]}))
    rows = [m["row"] for m in msgs if m["type"] == "row"]
    assert [(r["n"], r["n_species"]) for r in rows] == [(5, 5), (5, 5), (10, 10), (10, 10)]
