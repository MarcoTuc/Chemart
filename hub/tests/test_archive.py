"""The archive: deleting hides and keeps; only a purge removes, from disk too."""

from __future__ import annotations

import hashlib
import os
import sqlite3
import time

import pytest

from chemart_hub import db, service
from chemart_hub.storage import BlobStore

from conftest import network_files, tiny_files
from test_admin import HTML, _csrf, login, post

SECRET = b"# Tiny chem\n\nA sentence that exists in exactly one repo.\n"
SHARED = b"# Shared\n\nThe same README in two repos.\n"


@pytest.fixture
def site(hub):
    hub.user("root", admin=True)
    hub.user("alice")
    hub.user("bob")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": SECRET}))
    hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": SECRET, "NOTES.md": SHARED}))
    hub.create("bob", "bruss-snapshot", "network")
    hub.push("bob", "bob/bruss-snapshot", {**network_files(), "README.md": SHARED})
    return hub


def blob_exists(hub, data: bytes) -> bool:
    return BlobStore(hub.settings.blobs_dir).has(hashlib.sha256(data).hexdigest())


def in_database_files(hub, needle: str) -> bool:
    raw = b"".join(p.read_bytes() for p in hub.settings.data_dir.glob("hub.sqlite3*"))
    return needle.encode() in raw


def archived_ids(hub, table: str = "repos") -> list[int]:
    c = db.connect(hub.settings.db_path)
    try:
        return [r[0] for r in c.execute(f"SELECT id FROM {table} WHERE archived_at IS NOT NULL")]
    finally:
        c.close()


def test_deleting_archives_and_restore_brings_it_back(site):
    site.client.post("/api/repos/alice/tiny-chem/like", headers=site.headers("bob"))
    head = site.client.get("/api/repos/alice/tiny-chem").json()["head"]
    r = site.client.delete("/api/repos/alice/tiny-chem", headers=site.headers("alice"))
    assert r.json() == {"deleted": "alice/tiny-chem", "archived": True}

    # gone from everywhere...
    assert site.client.get("/api/repos/alice/tiny-chem").status_code == 404
    assert site.client.get(f"/alice/tiny-chem/resolve/{head}/README.md").status_code == 404
    assert site.client.get("/alice/tiny-chem", headers=HTML).status_code == 404
    assert site.client.get("/api/repos", params={"search": "replicator"}).json()["total"] == 0
    assert "tiny-chem" not in site.client.get("/alice", headers=HTML).text
    # ...but nothing is erased
    assert blob_exists(site, SECRET) and len(archived_ids(site)) == 1

    csrf = login(site, "root")
    page = site.client.get("/admin/archive", headers=HTML).text
    assert "alice/tiny-chem" in page and "by alice" in page
    rid = archived_ids(site)[0]
    assert post(site, f"/admin/archive/repos/{rid}/restore", csrf).status_code == 303
    info = site.client.get("/api/repos/alice/tiny-chem").json()
    assert info["head"] == head and info["likes"] == 1                  # exactly as it was
    assert site.client.get(f"/alice/tiny-chem/resolve/{head}/README.md").content == SECRET


def test_names_are_free_again_but_restore_does_not_clobber(site):
    site.client.delete("/api/repos/alice/tiny-chem", headers=site.headers("alice"))
    assert site.create("alice", "tiny-chem", "generator").status_code == 201     # owner reuses the name
    csrf = login(site, "root")
    assert "name in use again" in site.client.get("/admin/archive", headers=HTML).text
    r = post(site, f"/admin/archive/repos/{archived_ids(site)[0]}/restore", csrf)
    assert r.status_code == 409 and "exists" in r.text


def test_purge_removes_rows_files_and_traces(site):
    head = site.client.get("/api/repos/alice/tiny-chem").json()["head"]
    site.client.delete("/api/repos/alice/tiny-chem", headers=site.headers("alice"))
    rid = archived_ids(site)[0]
    csrf = login(site, "root")

    confirm = site.client.get("/admin/archive/purge", params={"kind": "repo", "id": rid}, headers=HTML)
    assert confirm.status_code == 200 and "cannot be undone" in confirm.text
    assert post(site, "/admin/archive/purge", csrf, kind="repo", id=str(rid), confirm="nope").status_code == 422
    assert blob_exists(site, SECRET)                                     # a wrong confirmation purges nothing
    assert in_database_files(site, head)                                 # (so the check below means something)

    done = post(site, "/admin/archive/purge", csrf, kind="repo", id=str(rid), confirm="alice/tiny-chem")
    assert done.status_code == 303 and "freed=" in done.headers["location"]
    assert not blob_exists(site, SECRET)                                 # its own file: gone from disk
    assert blob_exists(site, SHARED)                                     # still used by bob's repo
    assert not in_database_files(site, head)                             # no trace left in the database
    assert archived_ids(site) == []
    c = db.connect(site.settings.db_path)
    assert c.execute("SELECT COUNT(*) FROM commits WHERE id = ?", (head,)).fetchone()[0] == 0
    assert [r["action"] for r in c.execute("SELECT action FROM admin_log")][-1] == "purge repo"
    c.close()
    assert site.client.get("/bob/bruss-snapshot/resolve/main/README.md").content == SHARED


def test_deleted_accounts_are_archived_with_their_repos(site):
    csrf = login(site, "root")
    assert post(site, "/admin/users/alice/delete", csrf, confirm="alice").status_code == 303
    assert site.client.get("/api/whoami", headers=site.headers("alice")).status_code == 401
    assert site.client.get("/alice", headers=HTML).status_code == 404
    assert site.client.get("/api/repos/alice/tiny-chem").status_code == 404
    site.client.cookies.clear()
    page = site.client.get("/login", headers=HTML)
    bad = site.client.post("/login", data={"username": "alice", "password": "correct horse battery",
                                           "csrf": _csrf(page.text)}, headers=HTML)
    assert bad.status_code == 400                                        # as if it did not exist
    signup = site.client.get("/signup", headers=HTML)
    taken = site.client.post("/signup", data={"username": "alice", "password": "correct horse battery",
                                              "csrf": _csrf(signup.text)}, headers=HTML)
    assert taken.status_code == 400 and "taken" in taken.text           # the namespace stays reserved

    csrf = login(site, "root")
    page = site.client.get("/admin/archive", headers=HTML).text
    assert "restore alice to restore it" in page
    aid = archived_ids(site, "accounts")[0]
    assert post(site, f"/admin/archive/accounts/{aid}/restore", csrf).status_code == 303
    assert site.client.get("/api/whoami", headers=site.headers("alice")).json()["name"] == "alice"
    assert site.client.get("/api/repos/alice/tiny-chem").json()["head"]

    post(site, "/admin/users/alice/delete", csrf, confirm="alice")
    aid = archived_ids(site, "accounts")[0]
    done = post(site, "/admin/archive/purge", csrf, kind="account", id=str(aid), confirm="alice")
    assert done.status_code == 303
    assert not blob_exists(site, SECRET) and blob_exists(site, SHARED)
    c = db.connect(site.settings.db_path)
    assert service.account(c, "alice", include_archived=True) is None     # the name is free now
    c.close()


def test_purge_everything(site):
    site.client.delete("/api/repos/alice/tiny-chem", headers=site.headers("alice"))
    site.client.delete("/api/repos/bob/bruss-snapshot", headers=site.headers("bob"))
    csrf = login(site, "root")
    assert "Purge the whole archive" in site.client.get("/admin/archive", headers=HTML).text
    assert post(site, "/admin/archive/purge", csrf, kind="all", confirm="purge everything").status_code == 303
    assert not blob_exists(site, SECRET) and not blob_exists(site, SHARED)
    assert archived_ids(site) == []
    c = db.connect(site.settings.db_path)
    assert c.execute("SELECT COUNT(*) FROM commits").fetchone()[0] == 0
    c.close()


def test_orphans_wait_for_the_grace_period(site):
    store = BlobStore(site.settings.blobs_dir)
    old, fresh = store.put_bytes(b"orphan from a failed push"), store.put_bytes(b"upload still in flight")
    hour_ago = time.time() - 7200
    os.utime(store.path(old), (hour_ago, hour_ago))
    csrf = login(site, "root")
    assert post(site, "/admin/archive/gc", csrf).status_code == 303
    assert not store.has(old) and store.has(fresh)
    assert blob_exists(site, SECRET)                                     # live repos untouched


def test_only_superadmins_touch_the_archive(site):
    site.client.delete("/api/repos/bob/bruss-snapshot", headers=site.headers("bob"))
    csrf = login(site, "alice")
    assert site.client.get("/admin/archive", headers=HTML).status_code == 403
    rid = archived_ids(site)[0]
    assert post(site, f"/admin/archive/repos/{rid}/restore", csrf).status_code == 403
    assert post(site, "/admin/archive/purge", csrf, kind="all", confirm="purge everything").status_code == 403
    assert blob_exists(site, SHARED)


def test_upgrade_from_a_v2_database(tmp_path):
    """An existing hub keeps its data when the archive schema arrives."""
    from chemart_hub.config import Settings

    path = tmp_path / "data" / "hub.sqlite3"
    path.parent.mkdir(parents=True)
    conn = sqlite3.connect(path)
    conn.executescript(db._SCHEMA_V1)
    conn.executescript(db._SCHEMA_V2)
    conn.execute("PRAGMA user_version = 2")
    conn.execute("INSERT INTO accounts (id, name, created_at) VALUES (1, 'old', 'x')")
    conn.execute("INSERT INTO repos (id, owner_id, name, repo_type, head, created_at, updated_at, featured) "
                 "VALUES (7, 1, 'kept', 'network', 'c0ffee', 'x', 'x', 1)")
    conn.execute("INSERT INTO commits (id, repo_id, message, created_at, manifest) "
                 "VALUES ('c0ffee', 7, 'm', 'x', '{\"files\": []}')")
    conn.commit()
    conn.close()

    db.migrate(Settings(tmp_path / "data").db_path)
    c = db.connect(path)
    assert c.execute("PRAGMA user_version").fetchone()[0] == db.SCHEMA_VERSION
    row = c.execute("SELECT * FROM repos WHERE id = 7").fetchone()
    assert row["name"] == "kept" and row["featured"] == 1 and row["archived_at"] is None
    assert c.execute("SELECT repo_id FROM commits WHERE id = 'c0ffee'").fetchone()[0] == 7
    c.execute("UPDATE repos SET archived_at = 'x' WHERE id = 7")
    c.execute("INSERT INTO repos (owner_id, name, repo_type, created_at, updated_at) VALUES (1, 'kept', 'network', 'y', 'y')")
    with pytest.raises(sqlite3.IntegrityError):                          # only one *live* repo per name
        c.execute("INSERT INTO repos (owner_id, name, repo_type, created_at, updated_at) VALUES (1, 'kept', 'network', 'z', 'z')")
    c.close()


def test_upgrade_from_a_v4_database_renames_network_to_type(tmp_path):
    """v5: the old given/generated facet is cleared, to be refilled as a type."""
    from chemart_hub.config import Settings

    path = tmp_path / "data" / "hub.sqlite3"
    path.parent.mkdir(parents=True)
    conn = sqlite3.connect(path)
    for script in (db._SCHEMA_V1, db._SCHEMA_V2, db._SCHEMA_V3, db._SCHEMA_V4):
        conn.executescript(script)
    conn.execute("PRAGMA user_version = 4")
    conn.execute("INSERT INTO accounts (id, name, created_at) VALUES (1, 'old', 'x')")
    conn.execute("INSERT INTO repos (id, owner_id, name, repo_type, created_at, updated_at, network) "
                 "VALUES (7, 1, 'kept', 'generator', 'x', 'x', 'generated')")
    conn.commit()
    conn.close()

    db.migrate(Settings(tmp_path / "data").db_path)
    c = db.connect(path)
    assert c.execute("PRAGMA user_version").fetchone()[0] == db.SCHEMA_VERSION
    row = c.execute("SELECT * FROM repos WHERE id = 7").fetchone()
    assert row["name"] == "kept" and row["chem_type"] is None and "network" not in row.keys()
    c.close()


def test_cli_gc(site, capsys):
    from chemart_hub.cli import main

    assert main(["--data-dir", str(site.settings.data_dir), "gc"]) == 0
    assert "orphaned files" in capsys.readouterr().out
    assert blob_exists(site, SECRET)
