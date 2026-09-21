"""Fixtures: an app on a temporary data folder, users with tokens, and the
files of a small generator repo and a network repo."""

from __future__ import annotations

import hashlib
import json
import textwrap
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

import chemart
from chemart_hub import auth, db, service
from chemart_hub.app import create_app
from chemart_hub.config import Settings

TINY_YAML = textwrap.dedent("""\
    hub:
      repo_type: generator
      license: MIT
      tags: [toy, autocatalysis]
    chemistries:
    - id: tiny-chem
      name: Tiny autocatalytic chemistry
      intuition: >
        A replicator X copies itself from food F and decays; the smallest
        network with autocatalysis in it.
      family: core
      kind: generator
      constructive: false
      S: {repr: "F (food), X (replicator)", definition: explicit}
      R: {definition: explicit, arity: [1, 2], scheme: "F + X -> 2X ; X -> 0"}
      A: {reactor: [ode], dilution: none}
      params:
        - {name: k, type: float, default: 1.0, min: 0, role: kinetic, meaning: "replication rate"}
        - {name: d, type: float, default: 0.1, min: 0, role: kinetic, meaning: "decay rate"}
      provides: [topology, stoichiometry, catalysts, rate-constants]
      fidelity: original
    """)

TINY_GENERATOR = textwrap.dedent('''\
    """Tiny autocatalytic chemistry."""
    from chemart.helpers.explicit import network


    def generate(p, rng):
        return network([("F + X -> 2 X", p.k), ("X -> ", p.d)], species=["F", "X"])
    ''')


def tiny_generate(p, rng):
    from chemart.helpers.explicit import network

    return network([("F + X -> 2 X", p.k), ("X -> ", p.d)], species=["F", "X"])


def tiny_files(**overrides: bytes | None) -> dict[str, bytes]:
    """A valid generator repo; overrides replace (or with None, drop) files."""
    from chemart import api
    from chemart.hub import _format

    _, entry = _format.parse_chemart_yaml(TINY_YAML.encode())
    preview = api.run_generator(entry, tiny_generate, 0, {})
    preview.chemistry = "alice/tiny-chem"
    files = {
        "chemart.yaml": TINY_YAML.encode(),
        "generator.py": TINY_GENERATOR.encode(),
        "preview.json": json.dumps(preview.to_dict()).encode(),
        "README.md": b"# Tiny chem\n\nA replicator and its food.\n",
    }
    for path, data in overrides.items():
        if data is None:
            files.pop(path, None)
        else:
            files[path] = data
    return files


def network_files() -> dict[str, bytes]:
    net = chemart.generate_network("brusselator", seed=1)
    return {
        "network.json": json.dumps(net.to_dict()).encode(),
        "chemart.yaml": b"hub:\n  repo_type: network\n  title: A Brusselator snapshot\n  tags: [oscillator]\n",
        "README.md": b"The Brusselator at default parameters.\n",
    }


class Hub:
    """A TestClient plus helpers for pushing files as a given user."""

    def __init__(self, client: TestClient, settings: Settings):
        self.client = client
        self.settings = settings
        self.tokens: dict[str, str] = {}

    def user(self, name: str, *, admin: bool = False) -> str:
        c = db.connect(self.settings.db_path)
        try:
            row = service.create_user(c, name, "correct horse battery", admin=admin)
            self.tokens[name] = auth.create_token(c, row["id"], "test", "write")
        finally:
            c.close()
        return self.tokens[name]

    def headers(self, user: str | None) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.tokens[user]}"} if user else {}

    def create(self, user: str, name: str, repo_type: str, namespace: str | None = None):
        body = {"name": name, "repo_type": repo_type}
        if namespace:
            body["namespace"] = namespace
        return self.client.post("/api/repos", json=body, headers=self.headers(user))

    def push(self, user: str, repo: str, files: dict[str, bytes], *, parent="head",
             replace: bool = True, message: str = "push"):
        ops = []
        for path, data in files.items():
            sha = hashlib.sha256(data).hexdigest()
            r = self.client.put(f"/api/repos/{repo}/blobs/{sha}", content=data, headers=self.headers(user))
            if r.status_code >= 400:
                return r
            ops.append({"op": "add", "path": path, "sha256": sha})
        if parent == "head":
            parent = self.client.get(f"/api/repos/{repo}").json()["head"]
        return self.client.post(
            f"/api/repos/{repo}/commit/main",
            json={"message": message, "parent_commit": parent, "operations": ops, "replace": replace},
            headers=self.headers(user),
        )


@pytest.fixture
def settings(tmp_path: Path) -> Settings:
    return Settings(tmp_path / "data", public_url="http://hub.test")


@pytest.fixture
def hub(settings: Settings) -> Hub:
    with TestClient(create_app(settings)) as client:
        yield Hub(client, settings)
