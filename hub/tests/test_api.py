"""The JSON API: accounts, repos, commits, revisions, validation, search."""

from __future__ import annotations

import hashlib
import json

import pytest

from conftest import network_files, tiny_files


def test_whoami_needs_a_token(hub):
    hub.user("alice")
    assert hub.client.get("/api/whoami").status_code == 401
    r = hub.client.get("/api/whoami", headers=hub.headers("alice"))
    assert r.json()["name"] == "alice" and r.json()["scope"] == "write"
    bad = hub.client.get("/api/whoami", headers={"Authorization": "Bearer chm_nope"})
    assert bad.status_code == 401 and "invalid" in bad.json()["error"]


def test_generator_repo_round_trip(hub):
    hub.user("alice")
    assert hub.create("alice", "tiny-chem", "generator").status_code == 201
    r = hub.push("alice", "alice/tiny-chem", tiny_files())
    assert r.status_code == 200, r.json()
    commit = r.json()["commit"]
    assert len(commit) == 64 and r.json()["url"] == "http://hub.test/alice/tiny-chem"

    info = hub.client.get("/api/repos/alice/tiny-chem").json()
    assert info["head"] == commit
    assert info["title"] == "Tiny autocatalytic chemistry"
    assert info["has_code"] is True and info["fidelity"] == "original"
    assert info["tags"] == ["toy", "autocatalysis"]
    assert "catalysts" in info["provides"]

    rev = hub.client.get(f"/api/repos/alice/tiny-chem/revision/{commit[:7]}").json()
    assert rev["commit"] == commit and rev["author"] == "alice"
    assert {f["path"] for f in rev["files"]} == {"chemart.yaml", "generator.py", "preview.json", "README.md"}

    raw = hub.client.get(f"/alice/tiny-chem/resolve/{commit}/generator.py")
    assert raw.status_code == 200 and raw.content == tiny_files()["generator.py"]
    assert raw.headers["x-content-type-options"] == "nosniff"
    assert raw.headers["content-type"].startswith("text/plain")
    assert "immutable" in raw.headers["cache-control"]
    assert raw.headers["x-chemart-commit"] == commit
    assert hub.client.get("/alice/tiny-chem/resolve/main/generator.py").headers["cache-control"] == "no-cache"


def test_network_repo_and_downloads(hub):
    hub.user("bob")
    hub.create("bob", "bruss-snapshot", "network")
    assert hub.push("bob", "bob/bruss-snapshot", network_files()).status_code == 200
    info = hub.client.get("/api/repos/bob/bruss-snapshot").json()
    assert info["title"] == "A Brusselator snapshot" and info["has_code"] is False
    assert info["downloads"] == 0
    hub.client.get("/bob/bruss-snapshot/resolve/main/network.json")
    hub.client.get("/bob/bruss-snapshot/resolve/main/README.md")  # not the primary file
    assert hub.client.get("/api/repos/bob/bruss-snapshot").json()["downloads"] == 1


def test_revisions_are_immutable_and_main_moves(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    first = hub.push("alice", "alice/tiny-chem", tiny_files()).json()["commit"]
    second = hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": b"# v2\n"})).json()["commit"]
    assert first != second
    assert hub.client.get("/alice/tiny-chem/resolve/main/README.md").content == b"# v2\n"
    assert hub.client.get(f"/alice/tiny-chem/resolve/{first}/README.md").content.startswith(b"# Tiny chem")
    history = hub.client.get("/api/repos/alice/tiny-chem/commits").json()["commits"]
    assert [c["commit"] for c in history] == [second, first]
    assert history[0]["parent"] == first


def test_identical_push_is_a_no_op(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    first = hub.push("alice", "alice/tiny-chem", tiny_files()).json()
    again = hub.push("alice", "alice/tiny-chem", tiny_files()).json()
    assert again == {**first, "unchanged": True}


def test_stale_parent_is_a_conflict(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    first = hub.push("alice", "alice/tiny-chem", tiny_files()).json()["commit"]
    hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": b"v2"}))
    r = hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": b"v3"}), parent=first)
    assert r.status_code == 409 and "has moved" in r.json()["error"]


def test_partial_commits_and_deletes(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files())
    r = hub.push("alice", "alice/tiny-chem", {"LICENSE": b"MIT"}, replace=False)
    assert r.status_code == 200
    head = hub.client.get("/api/repos/alice/tiny-chem").json()["head"]
    files = {f["path"] for f in hub.client.get(f"/api/repos/alice/tiny-chem/revision/{head}").json()["files"]}
    assert "LICENSE" in files and "generator.py" in files
    r = hub.client.post("/api/repos/alice/tiny-chem/commit/main", headers=hub.headers("alice"),
                        json={"parent_commit": head, "operations": [{"op": "delete", "path": "generator.py"}]})
    assert r.status_code == 422 and any("generator.py" in p for p in r.json()["problems"])


def test_only_owners_and_org_members_push(hub):
    hub.user("alice")
    hub.user("mallory")
    hub.create("alice", "tiny-chem", "generator")
    r = hub.push("mallory", "alice/tiny-chem", tiny_files())
    assert r.status_code == 403
    assert hub.create("mallory", "x", "network", namespace="alice").status_code == 403

    from chemart_hub import db, service
    c = db.connect(hub.settings.db_path)
    service.create_org(c, "lab", "alice")
    c.close()
    assert hub.create("alice", "tiny-chem", "generator", namespace="lab").status_code == 201
    assert hub.push("alice", "lab/tiny-chem", tiny_files()).status_code == 200
    assert hub.push("mallory", "lab/tiny-chem", tiny_files(**{"README.md": b"x"})).status_code == 403


def test_read_only_tokens_cannot_write(hub):
    from chemart_hub import auth, db, service
    c = db.connect(hub.settings.db_path)
    row = service.create_user(c, "reader", "correct horse battery")
    hub.tokens["reader"] = auth.create_token(c, row["id"], "ro", "read")
    c.close()
    r = hub.create("reader", "x", "network")
    assert r.status_code == 403 and "read-only" in r.json()["error"]


@pytest.mark.parametrize(
    "overrides, fragment",
    [
        ({"chemart.yaml": b"hub: [unclosed"}, "not valid YAML"),
        ({"chemart.yaml": b"a: &x [1]\nb: *x\n"}, "anchors and aliases"),
        ({"generator.py": b"def generate(p, rng:\n"}, "not valid Python"),
        ({"generator.py": b"def make(p, rng):\n    pass\n"}, "must define a top-level function generate"),
        ({"generator.py": None}, "needs generator.py"),
        ({"preview.json": None}, "needs preview.json"),
        ({"preview.json": b'{"species": NaN}'}, "not valid JSON"),
        ({"preview.json": b'{"species": [], "reactions": [{"reactants": ["X"]}]}'}, "preview.json"),
        ({"chemart.yaml": tiny_files()["chemart.yaml"].replace(b"id: tiny-chem", b"id: other-chem")},
         "must equal the repo name"),
        ({"chemart.yaml": tiny_files()["chemart.yaml"].replace(b"fidelity: original", b"fidelity: vibes")},
         "fidelity must be one of"),
        ({"chemart.yaml": tiny_files()["chemart.yaml"].replace(b"name: d,", b"name: seed,")},
         "'seed' is reserved"),
        ({"chemart.yaml": tiny_files()["chemart.yaml"].replace(b"  license: MIT", b"  builtin: tiny-chem")},
         "reserved for the official"),
        ({"notes.html": b"<script>"}, "file type not allowed"),
        ({"big.md": b"x" * (1024 * 1024 + 1)}, "too large"),
    ],
)
def test_invalid_generator_repos_are_refused(hub, overrides, fragment):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    r = hub.push("alice", "alice/tiny-chem", tiny_files(**overrides))
    assert r.status_code == 422, r.json()
    assert any(fragment in p for p in r.json().get("problems", [r.json()["error"]])), r.json()
    assert hub.client.get("/api/repos/alice/tiny-chem").json()["head"] is None


def test_preview_must_use_default_parameters(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    preview = json.loads(tiny_files()["preview.json"])
    preview["params"]["k"] = 2.0
    r = hub.push("alice", "alice/tiny-chem", tiny_files(**{"preview.json": json.dumps(preview).encode()}))
    assert r.status_code == 422 and any("default parameters" in p for p in r.json()["problems"])


def test_network_repos_carry_no_code(hub):
    hub.user("bob")
    hub.create("bob", "snap", "network")
    r = hub.push("bob", "bob/snap", {**network_files(), "run.py": b"print(1)\n"})
    assert r.status_code == 422 and any("carries no code" in p for p in r.json()["problems"])


def test_blob_hash_is_checked(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    wrong = hashlib.sha256(b"something else").hexdigest()
    r = hub.client.put(f"/api/repos/alice/tiny-chem/blobs/{wrong}", content=b"data", headers=hub.headers("alice"))
    assert r.status_code == 422 and "hashes to" in r.json()["error"]


def test_unknown_blob_in_commit(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    r = hub.client.post("/api/repos/alice/tiny-chem/commit/main", headers=hub.headers("alice"),
                        json={"parent_commit": None,
                              "operations": [{"op": "add", "path": "README.md", "sha256": "0" * 64}]})
    assert r.status_code == 422 and "was not uploaded" in r.json()["problems"][0]


def test_revision_errors(hub):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    assert hub.client.get("/api/repos/alice/tiny-chem/revision/main").status_code == 404  # empty repo
    hub.push("alice", "alice/tiny-chem", tiny_files())
    assert hub.client.get("/api/repos/alice/tiny-chem/revision/abc").status_code == 400
    assert hub.client.get("/api/repos/alice/tiny-chem/revision/0000000").status_code == 404
    assert hub.client.get("/api/repos/nobody/nothing").status_code == 404


def test_search_filters_and_likes(hub):
    hub.user("alice")
    hub.user("bob")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files())
    hub.create("bob", "bruss-snapshot", "network")
    hub.push("bob", "bob/bruss-snapshot", network_files())
    hub.create("bob", "empty", "network")

    def ids(**params):
        return [r["id"] for r in hub.client.get("/api/repos", params=params).json()["repos"]]

    assert set(ids()) == {"alice/tiny-chem", "bob/bruss-snapshot"}      # empty repos hidden
    assert ids(repo_type="network") == ["bob/bruss-snapshot"]
    assert ids(search="replicator") == ["alice/tiny-chem"]
    assert ids(search="brusselator") == ["bob/bruss-snapshot"]
    assert ids(tag="oscillator") == ["bob/bruss-snapshot"]
    assert ids(provides="catalysts,rate-constants", repo_type="generator") == ["alice/tiny-chem"]
    assert ids(author="alice") == ["alice/tiny-chem"]
    assert ids(has_code="true") == ["alice/tiny-chem"]

    r = hub.client.post("/api/repos/bob/bruss-snapshot/like", headers=hub.headers("alice"))
    assert r.json() == {"likes": 1, "liked": True}
    hub.client.post("/api/repos/bob/bruss-snapshot/like", headers=hub.headers("alice"))  # idempotent
    assert ids(sort="likes")[0] == "bob/bruss-snapshot"
    assert hub.client.get("/api/repos/bob/bruss-snapshot").json()["likes"] == 1
    assert hub.client.get("/api/repos", params={"sort": "sideways"}).status_code == 422


def test_names_are_checked(hub):
    hub.user("alice")
    assert hub.create("alice", "Bad_Name", "network").status_code == 422
    assert hub.create("alice", "ok-name", "dataset").status_code == 422
    assert hub.create("alice", "ok-name", "network").status_code == 201
    assert hub.create("alice", "ok-name", "network").status_code == 409


def test_delete_repo(hub):
    hub.user("alice")
    hub.user("mallory")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files())
    assert hub.client.delete("/api/repos/alice/tiny-chem", headers=hub.headers("mallory")).status_code == 403
    assert hub.client.delete("/api/repos/alice/tiny-chem", headers=hub.headers("alice")).status_code == 200
    assert hub.client.get("/api/repos/alice/tiny-chem").status_code == 404
    assert hub.client.get("/api/repos", params={"search": "replicator"}).json()["total"] == 0
