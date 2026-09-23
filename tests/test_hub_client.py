"""The hub client without a server: ids, cache integrity, offline use, the gate.

`chemart.hub._http.transport` is replaced by a fake hub that answers from a
dict, so these tests need neither the network nor the chemart-hub package.
"""

from __future__ import annotations

import hashlib
import json

import pytest

import chemart
from chemart import hub
from chemart.hub import _http
from chemart.hub._ids import check_path, parse_repo_id

COMMIT = "c" * 64


class FakeHub:
    """Serves one generator repo and one network repo at COMMIT."""

    def __init__(self):
        net = chemart.generate_network("brusselator", seed=1)
        self.files = {
            "bob/snap": {"network.json": json.dumps(net.to_dict()).encode()},
            "alice/tiny": {
                "chemart.yaml": _tiny_yaml(),
                "generator.py": b"def generate(p, rng):\n    raise AssertionError('must not run')\n",
            },
        }
        self.types = {"bob/snap": "network", "alice/tiny": "generator"}
        self.tamper: set[str] = set()
        self.requests: list[str] = []
        self.network = net

    def __call__(self, method, url, headers, body, timeout):
        path = url.split("//", 1)[1].split("/", 1)[1]
        self.requests.append(f"{method} /{path}")
        parts = path.split("/")
        if parts[:2] == ["api", "repos"] and len(parts) == 6 and parts[4] == "revision":
            repo = f"{parts[2]}/{parts[3]}"
            if repo not in self.files:
                return _http.Response(404, b'{"error": "repo not found"}')
            files = [{"path": p, "sha256": hashlib.sha256(d).hexdigest(), "size": len(d)}
                     for p, d in self.files[repo].items()]
            answer = {"commit": COMMIT, "repo_type": self.types[repo], "files": files}
            if "hostile-path" in self.tamper:
                answer["files"].append({"path": "../../.bashrc", "sha256": "0" * 64, "size": 1})
            return _http.Response(200, json.dumps(answer).encode())
        if len(parts) >= 5 and parts[2] == "resolve":
            data = self.files[f"{parts[0]}/{parts[1]}"][parts[4]]
            return _http.Response(200, data + (b" " if "content" in self.tamper else b""))
        return _http.Response(404, b'{"error": "no such route"}')


def _tiny_yaml() -> bytes:
    import yaml

    from chemart.catalog import CATALOG_DIR

    entry = yaml.safe_load((CATALOG_DIR / "brusselator.yaml").read_text())["chemistries"][0]
    entry.update(id="tiny", fidelity="original")
    return yaml.safe_dump({"hub": {"repo_type": "generator"}, "chemistries": [entry]}).encode()


@pytest.fixture
def fake(monkeypatch, tmp_path):
    server = FakeHub()
    monkeypatch.setattr(_http, "transport", server)
    monkeypatch.setenv("CHEMART_HUB_URL", "http://hub.test")
    monkeypatch.setenv("CHEMART_HOME", str(tmp_path))
    monkeypatch.delenv("CHEMART_HUB_TOKEN", raising=False)
    monkeypatch.delenv("CHEMART_HUB_OFFLINE", raising=False)
    return server


@pytest.mark.parametrize("text, expected", [
    ("alice/my-chem", ("alice", "my-chem", "main")),
    ("alice/my-chem@3f2a9c1", ("alice", "my-chem", "3f2a9c1")),
])
def test_repo_ids(text, expected):
    repo, rev = parse_repo_id(text)
    assert (repo.namespace, repo.name, rev) == expected


@pytest.mark.parametrize("bad", ["alice", "alice/My_Chem", "a/b/c", "alice/x@HEAD", "alice/x@abc"])
def test_bad_repo_ids(bad):
    with pytest.raises(ValueError):
        parse_repo_id(bad)


def test_conflicting_revisions():
    with pytest.raises(ValueError, match="also given"):
        parse_repo_id("alice/x@1234567", "7654321")


@pytest.mark.parametrize("bad", ["../x.py", ".env", "a/b/c.py", "x.html", "dir/.hidden.py", "", "x/../y.py"])
def test_hostile_paths(bad):
    with pytest.raises(ValueError):
        check_path(bad)


def test_load_network_and_cache(fake):
    assert chemart.load_network("bob/snap") == fake.network
    n = len(fake.requests)
    assert chemart.load_network(f"bob/snap@{COMMIT}") == fake.network
    assert len(fake.requests) == n, "a cached full commit id needs no network"


def test_offline_mode(fake, monkeypatch):
    monkeypatch.setenv("CHEMART_HUB_OFFLINE", "1")
    with pytest.raises(hub.HubConnectionError, match="not in the cache"):
        chemart.load_network("bob/snap")
    monkeypatch.delenv("CHEMART_HUB_OFFLINE")
    chemart.load_network("bob/snap")
    monkeypatch.setenv("CHEMART_HUB_OFFLINE", "1")
    fake.requests.clear()
    assert chemart.load_network("bob/snap") == fake.network
    assert chemart.load_network("bob/snap@ccccccc") == fake.network
    assert fake.requests == []


def test_integrity_check(fake):
    fake.tamper.add("content")
    with pytest.raises(hub.HubError, match="integrity check"):
        chemart.load_network("bob/snap")
    fake.tamper.clear()
    assert chemart.load_network("bob/snap") == fake.network      # nothing bad was cached


def test_hostile_server_paths_are_refused(fake):
    fake.tamper.add("hostile-path")
    with pytest.raises(hub.HubError, match="malformed answer"):
        chemart.load_network("bob/snap")


def test_not_found(fake):
    with pytest.raises(hub.RepoNotFoundError):
        chemart.load_network("bob/nothing")


def test_trust_gate_runs_nothing(fake):
    info = chemart.describe_chemistry("alice/tiny")
    assert info["hub"]["commit"] == COMMIT
    with pytest.raises(ValueError, match="trust_remote_code=True") as err:
        chemart.generate_network("alice/tiny")
    assert "http://hub.test/alice/tiny/blob/" in str(err.value)
    assert not any("generator.py" in r for r in fake.requests), "code was downloaded before trust"


def test_tool_calls_never_trust_remote_code(fake):
    for tool in ("generate_network", "simulate_network", "measure_network"):
        with pytest.raises(ValueError, match="trust_remote_code=True"):
            chemart.call_tool(tool, {"chemistry": "alice/tiny"})


def test_saved_tokens_are_per_hub_and_private(fake, tmp_path):
    from chemart.hub import _config

    _config.save_token("chm_one", "http://hub.test")
    _config.save_token("chm_two", "https://other.example")
    assert _config.token() == "chm_one"
    assert (tmp_path / "tokens.json").stat().st_mode & 0o777 == 0o600
    _config.save_token(None, "http://hub.test")
    assert _config.token() is None and _config.token("https://other.example") == "chm_two"


def test_version_specs():
    from chemart.hub._format import version_satisfies

    assert version_satisfies("0.2.1", ">=0.2, <1")
    assert not version_satisfies("1.0", ">=0.2, <1")
    assert version_satisfies("0.1.0", "==0.1")
