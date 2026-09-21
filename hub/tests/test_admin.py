"""Curation: the superadmin panel and repo editing from the browser."""

from __future__ import annotations

import re

import pytest

from chemart_hub import auth, db, service

from conftest import network_files, tiny_files

HTML = {"Accept": "text/html"}


def _csrf(html: str) -> str:
    return re.search(r'name="csrf" value="([^"]*)"', html).group(1)


def login(hub, name: str) -> str:
    """Log the test client in as `name`; return the session's CSRF token."""
    c = hub.client
    c.cookies.clear()
    page = c.get("/login", headers=HTML)
    r = c.post("/login", data={"username": name, "password": "correct horse battery", "csrf": _csrf(page.text)},
               headers=HTML, follow_redirects=False)
    assert r.status_code == 303, r.text
    return _csrf(c.get("/settings/account", headers=HTML).text)


@pytest.fixture
def site(hub):
    hub.user("root", admin=True)
    hub.user("alice")
    hub.user("bob")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files())
    hub.create("bob", "bruss-snapshot", "network")
    hub.push("bob", "bob/bruss-snapshot", network_files())
    return hub


def post(hub, url, csrf, **data):
    return hub.client.post(url, data={"csrf": csrf, **data}, headers=HTML, follow_redirects=False)


def test_only_superadmins_reach_the_panel(site):
    assert site.client.get("/admin", headers=HTML, follow_redirects=False).status_code == 303   # to login
    login(site, "alice")
    assert site.client.get("/admin", headers=HTML).status_code == 403
    assert site.client.get("/admin/users", headers=HTML).status_code == 403
    csrf = _csrf(site.client.get("/settings/account", headers=HTML).text)
    assert post(site, "/admin/repos/bob/bruss-snapshot/flags", csrf, hidden="1").status_code == 403
    assert "Admin</a>" not in site.client.get("/", headers=HTML).text

    login(site, "root")
    for url in ("/admin", "/admin/repos", "/admin/users", "/admin/site", "/admin/log"):
        assert site.client.get(url, headers=HTML).status_code == 200, url
    assert 'href="/admin"' in site.client.get("/", headers=HTML).text


def test_feature_hide_and_delete_any_repo(site):
    csrf = login(site, "root")
    assert post(site, "/admin/repos/alice/tiny-chem/flags", csrf, featured="1").status_code == 303
    home = site.client.get("/", headers=HTML).text
    assert "Curated picks" in home and "alice/tiny-chem" in home

    post(site, "/admin/repos/bob/bruss-snapshot/flags", csrf, hidden="1")
    site.client.cookies.clear()
    assert "bruss-snapshot" not in site.client.get("/browse", headers=HTML).text
    assert "bruss-snapshot" not in site.client.get("/api/repos").text
    assert site.client.get("/bob/bruss-snapshot", headers=HTML).status_code == 200    # still loads by id
    assert site.client.get("/bob/bruss-snapshot/resolve/main/network.json").status_code == 200

    csrf = login(site, "root")
    assert "bruss-snapshot" in site.client.get("/browse", headers=HTML).text       # admins still see it
    bad = post(site, "/alice/tiny-chem/edit/delete", csrf, confirm="wrong")
    assert bad.status_code == 422
    assert post(site, "/alice/tiny-chem/edit/delete", csrf, confirm="alice/tiny-chem").status_code == 303
    assert site.client.get("/api/repos/alice/tiny-chem").status_code == 404
    actions = [e["action"] for e in service_log(site)]
    assert {"feature", "hide", "archive repo"} <= set(actions)


def service_log(hub):
    c = db.connect(hub.settings.db_path)
    try:
        return [dict(r) for r in c.execute("SELECT * FROM admin_log ORDER BY id")]
    finally:
        c.close()


def test_edit_files_in_the_browser(site):
    csrf = login(site, "root")
    editor = site.client.get("/alice/tiny-chem/edit/file", params={"path": "README.md"}, headers=HTML)
    assert editor.status_code == 200 and "Tiny chem" in editor.text
    r = post(site, "/alice/tiny-chem/edit/file", csrf, path="README.md", original="README.md",
             content="# Curated\r\n\r\nRewritten by the curators.", message="Tidy the card")
    assert r.status_code == 303
    assert site.client.get("/alice/tiny-chem/resolve/main/README.md").content == b"# Curated\n\nRewritten by the curators.\n"
    history = site.client.get("/api/repos/alice/tiny-chem/commits").json()["commits"]
    assert history[0]["message"] == "Tidy the card" and history[0]["author"] == "root"

    # an edit that breaks the repo is refused, with reasons, and the text is kept
    broken = post(site, "/alice/tiny-chem/edit/file", csrf, path="chemart.yaml", original="chemart.yaml",
                  content="hub: [broken", message="")
    assert broken.status_code == 400 and "not valid YAML" in broken.text and "hub: [broken" in broken.text
    assert len(site.client.get("/api/repos/alice/tiny-chem/commits").json()["commits"]) == 2

    # add, rename and delete files
    post(site, "/alice/tiny-chem/edit/file", csrf, path="NOTES.md", original="", content="notes")
    post(site, "/alice/tiny-chem/edit/file", csrf, path="notes.txt", original="NOTES.md", content="notes")
    files = {f["path"] for f in site.client.get("/api/repos/alice/tiny-chem/revision/main").json()["files"]}
    assert "notes.txt" in files and "NOTES.md" not in files
    post(site, "/alice/tiny-chem/edit/file/delete", csrf, path="notes.txt")
    files = {f["path"] for f in site.client.get("/api/repos/alice/tiny-chem/revision/main").json()["files"]}
    assert "notes.txt" not in files
    refused = post(site, "/alice/tiny-chem/edit/file/delete", csrf, path="generator.py")
    assert refused.status_code == 400 and "needs generator.py" in refused.text


def test_edit_card_and_restore(site):
    csrf = login(site, "root")
    first = site.client.get("/api/repos/bob/bruss-snapshot").json()["head"]
    r = post(site, "/bob/bruss-snapshot/edit/metadata", csrf, title="The Brusselator, curated",
             description="A hand-picked snapshot.", tags="oscillator, classic", license="CC0-1.0")
    assert r.status_code == 303
    info = site.client.get("/api/repos/bob/bruss-snapshot").json()
    assert info["title"] == "The Brusselator, curated" and info["summary"] == "A hand-picked snapshot."
    assert info["tags"] == ["oscillator", "classic"] and info["license"] == "CC0-1.0"

    assert post(site, "/bob/bruss-snapshot/edit/restore", csrf, commit=first).status_code == 303
    info = site.client.get("/api/repos/bob/bruss-snapshot").json()
    assert info["title"] == "A Brusselator snapshot" and info["head"] != first     # a new commit
    commits = site.client.get("/bob/bruss-snapshot/commits", headers=HTML).text
    assert "Restore" in commits and "current" in commits


def test_owners_edit_their_own_but_do_not_curate(site):
    csrf = login(site, "alice")
    assert site.client.get("/alice/tiny-chem/edit", headers=HTML).status_code == 200
    assert "Curation" not in site.client.get("/alice/tiny-chem/edit", headers=HTML).text
    assert site.client.get("/bob/bruss-snapshot/edit", headers=HTML).status_code == 403
    r = post(site, "/bob/bruss-snapshot/edit/file", csrf, path="README.md", original="README.md", content="mine now")
    assert r.status_code == 403
    assert post(site, "/alice/tiny-chem/edit/metadata", csrf, title="Mine", description="", tags="",
                license="").status_code == 303
    assert service_log(site) == []                     # owners editing their own repos leave no admin trace


def test_manage_accounts(site):
    csrf = login(site, "root")
    # suspend: login refused, tokens dead, repos kept
    assert post(site, "/admin/users/alice/suspend", csrf, flag="1").status_code == 303
    assert site.client.get("/api/whoami", headers=site.headers("alice")).status_code == 401
    assert site.client.get("/api/repos/alice/tiny-chem").status_code == 200
    post(site, "/admin/users/alice/suspend", csrf, flag="0")
    assert site.client.get("/api/whoami", headers=site.headers("alice")).status_code == 200

    # reset a password, promote, demote
    post(site, "/admin/users/bob/password", csrf, password="a brand new password")
    c = db.connect(site.settings.db_path)
    assert auth.verify_password("a brand new password", service.account(c, "bob")["password_hash"])
    c.close()
    post(site, "/admin/users/alice/superadmin", csrf, flag="1")
    assert site.client.get("/api/whoami", headers=site.headers("alice")).json()["is_admin"] is True
    post(site, "/admin/users/alice/superadmin", csrf, flag="0")

    # the last superadmin cannot be removed, suspended or deleted
    assert post(site, "/admin/users/root/superadmin", csrf, flag="0").status_code == 409
    assert post(site, "/admin/users/root/suspend", csrf, flag="1").status_code == 409

    # delete an account with its repos
    assert post(site, "/admin/users/bob/delete", csrf, confirm="nope").status_code == 422
    assert post(site, "/admin/users/bob/delete", csrf, confirm="bob").status_code == 303
    assert site.client.get("/api/repos/bob/bruss-snapshot").status_code == 404
    assert site.client.get("/bob", headers=HTML).status_code == 404
    assert site.client.get("/api/repos", params={"search": "brusselator"}).json()["total"] == 0


def test_profiles_orgs_and_transfer(site):
    csrf = login(site, "root")
    post(site, "/settings/profile/alice", csrf, fullname="Alice Liddell", bio="Curious about **autocatalysis**.",
         email="alice@example.org")
    page = site.client.get("/alice", headers=HTML).text
    assert "Alice Liddell" in page and "<strong>autocatalysis</strong>" in page

    assert post(site, "/admin/orgs", csrf, name="lab", owner="alice").status_code == 303
    post(site, "/settings/profile/lab/members", csrf, user="bob", action="add")
    assert "bob" in site.client.get("/lab", headers=HTML).text
    assert post(site, "/admin/repos/alice/tiny-chem/transfer", csrf, owner="lab").status_code == 303
    assert site.client.get("/api/repos/lab/tiny-chem").json()["head"]
    assert site.client.get("/api/repos/alice/tiny-chem").status_code == 404


def test_site_text_is_editable_and_sanitised(site):
    csrf = login(site, "root")
    r = post(site, "/admin/site", csrf, home_title="Welcome to the chemical mart",
             announcement="**New:** the [autumn catalog](https://example.org) <script>alert(1)</script>",
             about="# Who we are\n\nCurators of artificial chemistries.", footer="")
    assert r.status_code == 303
    site.client.cookies.clear()
    home = site.client.get("/", headers=HTML).text
    assert "Welcome to the chemical mart" in home
    assert "<strong>New:</strong>" in home and "<script>alert" not in home
    assert "Who we are" in site.client.get("/about", headers=HTML).text
    assert "Chemart Hub · the one stop shop" in home                     # empty footer = default

    csrf = login(site, "root")
    post(site, "/admin/site", csrf, announcement="")
    assert "class=\"announcement\"" not in site.client.get("/", headers=HTML).text
    assert "edit site text" in [e["action"] for e in service_log(site)]


def test_change_own_password(site):
    csrf = login(site, "alice")
    bad = post(site, "/settings/account/password", csrf, current="wrong", password="new password 1",
               confirm="new password 1")
    assert bad.status_code == 400 and "current password is wrong" in bad.text
    ok = post(site, "/settings/account/password", csrf, current="correct horse battery",
              password="new password 1", confirm="new password 1")
    assert ok.status_code == 303
    c = db.connect(site.settings.db_path)
    assert auth.verify_password("new password 1", service.account(c, "alice")["password_hash"])
    c.close()


def test_cli_set_admin_and_password(site, capsys):
    from chemart_hub.cli import main

    data = str(site.settings.data_dir)
    assert main(["--data-dir", data, "set-admin", "alice"]) == 0
    assert main(["--data-dir", data, "set-password", "alice", "--password", "another long one"]) == 0
    c = db.connect(site.settings.db_path)
    row = service.account(c, "alice")
    assert row["is_admin"] == 1 and auth.verify_password("another long one", row["password_hash"])
    c.close()
    assert main(["--data-dir", data, "set-admin", "nobody"]) == 2
