"""The web site: pages render, forms need CSRF, user content is sanitised."""

from __future__ import annotations

import re

import pytest

from conftest import network_files, tiny_files

HTML = {"Accept": "text/html"}


def _csrf(page_html: str) -> str:
    return re.search(r'name="csrf" value="([^"]*)"', page_html).group(1)


def _login(hub, name="alice", password="correct horse battery"):
    page = hub.client.get("/login", headers=HTML)
    r = hub.client.post("/login", data={"username": name, "password": password, "csrf": _csrf(page.text), "next": "/"},
                        headers=HTML, follow_redirects=False)
    assert r.status_code == 303, r.text
    return r


@pytest.fixture
def stocked(hub):
    hub.user("alice")
    hub.user("bob")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files())
    hub.create("bob", "bruss-snapshot", "network")
    hub.push("bob", "bob/bruss-snapshot", network_files())
    return hub


def test_pages_render(stocked):
    c = stocked.client
    home = c.get("/", headers=HTML)
    assert home.status_code == 200
    assert "░░╔════╝" in home.text and "1</b> chemistry in stock" in home.text
    assert "alice/" in home.text and "tiny-chem" in home.text

    repo = c.get("/alice/tiny-chem", headers=HTML)
    assert repo.status_code == 200
    for fragment in ("Tiny autocatalytic chemistry", "How it works", "replication rate",
                     "trust_remote_code=True", "runs code", "Preview at default parameters", "F + X -&gt; 2 X"):
        assert fragment in repo.text, fragment

    net = c.get("/bob/bruss-snapshot", headers=HTML)
    assert "chemart.load_network(\n    &#34;bob/bruss-snapshot&#34;)" in net.text
    assert "runs code" not in net.text

    assert c.get("/alice/tiny-chem/tree/main", headers=HTML).status_code == 200
    blob = c.get("/alice/tiny-chem/blob/main/generator.py", headers=HTML)
    assert blob.status_code == 200 and "runs on the machine" in blob.text
    assert c.get("/alice/tiny-chem/commits", headers=HTML).status_code == 200
    assert c.get("/alice", headers=HTML).status_code == 200
    browse = c.get("/browse", params={"provides": "initial-state"}, headers=HTML)
    assert browse.status_code == 200 and "bruss-snapshot" in browse.text and "tiny-chem" not in browse.text
    assert c.get("/browse", params={"q": "brusselator"}, headers=HTML).text.count("bruss-snapshot") >= 1
    assert c.get("/new", headers=HTML).status_code == 200
    assert c.get("/static/pygments.css").headers["content-type"].startswith("text/css")


def test_missing_pages_are_html_404(stocked):
    r = stocked.client.get("/alice/nothing-here", headers=HTML)
    assert r.status_code == 404 and "Not on the shelf" in r.text
    assert stocked.client.get("/nobody", headers=HTML).status_code == 404
    assert stocked.client.get("/api/repos/nobody/x").json()["error"]


def test_security_headers(stocked):
    r = stocked.client.get("/", headers=HTML)
    assert "script-src 'self'" in r.headers["content-security-policy"]
    assert r.headers["x-content-type-options"] == "nosniff"
    raw = stocked.client.get("/alice/tiny-chem/resolve/main/README.md")
    assert "sandbox" in raw.headers["content-security-policy"]


def test_login_session_tokens_and_logout(stocked):
    c = stocked.client
    bad = c.post("/login", data={"username": "alice", "password": "nope", "csrf": _csrf(c.get("/login", headers=HTML).text)},
                 headers=HTML)
    assert bad.status_code == 400 and "Wrong username or password" in bad.text
    _login(stocked)
    page = c.get("/settings/tokens", headers=HTML)
    assert page.status_code == 200 and "API tokens" in page.text
    made = c.post("/settings/tokens", data={"name": "laptop", "scope": "write", "csrf": _csrf(page.text)}, headers=HTML)
    token = re.search(r'value="(chm_[^"]+)"', made.text).group(1)
    assert c.get("/api/whoami", headers={"Authorization": f"Bearer {token}"}).json()["name"] == "alice"

    # the session can like through the web form
    repo = c.get("/bob/bruss-snapshot", headers=HTML)
    c.post("/bob/bruss-snapshot/like", data={"csrf": _csrf(repo.text)}, headers=HTML)
    assert c.get("/api/repos/bob/bruss-snapshot").json()["likes"] == 1

    c.post("/logout", data={"csrf": _csrf(c.get("/", headers=HTML).text)}, headers=HTML)
    assert c.get("/settings/tokens", headers=HTML, follow_redirects=False).status_code == 303


def test_forms_need_csrf(stocked):
    c = stocked.client
    r = c.post("/login", data={"username": "alice", "password": "correct horse battery", "csrf": "forged"}, headers=HTML)
    assert r.status_code == 403
    _login(stocked)
    r = c.post("/settings/tokens", data={"name": "x", "scope": "write", "csrf": "forged"}, headers=HTML)
    assert r.status_code == 403
    # a session alone cannot write through the JSON API without the CSRF header
    r = c.post("/api/repos", json={"name": "sneaky", "repo_type": "network"})
    assert r.status_code == 403 and "CSRF" in r.json()["error"]


def test_signup_and_new_repo(hub):
    c = hub.client
    page = c.get("/signup", headers=HTML)
    r = c.post("/signup", data={"username": "dora", "password": "correct horse battery", "csrf": _csrf(page.text)},
               headers=HTML, follow_redirects=False)
    assert r.status_code == 303 and r.headers["location"].startswith("/settings/tokens")
    taken = c.post("/signup", data={"username": "browse", "password": "correct horse battery",
                                    "csrf": _csrf(page.text)}, headers=HTML)
    assert taken.status_code == 400 and "reserved" in taken.text
    new = c.get("/new", headers=HTML)
    r = c.post("/new", data={"namespace": "dora", "name": "my-net", "repo_type": "network", "csrf": _csrf(new.text)},
               headers=HTML, follow_redirects=False)
    assert r.headers["location"] == "/dora/my-net"
    assert "This repo is empty" in c.get("/dora/my-net", headers=HTML).text


def test_login_is_throttled(hub):
    from chemart_hub.routes import web

    hub.user("alice")
    web._throttle.failures.clear()
    csrf = _csrf(hub.client.get("/login", headers=HTML).text)
    for _ in range(10):
        hub.client.post("/login", data={"username": "alice", "password": "wrong", "csrf": csrf}, headers=HTML)
    r = hub.client.post("/login", data={"username": "alice", "password": "correct horse battery", "csrf": csrf},
                        headers=HTML)
    assert r.status_code == 429
    web._throttle.failures.clear()


@pytest.mark.parametrize("payload", [
    b"<script>alert(1)</script>",
    b"[click](javascript:alert(1))",
    b'<img src=x onerror="alert(1)">',
    b'<a href="javascript:alert(1)">x</a>',
    b"<iframe src=https://evil.example></iframe>",
])
def test_readme_is_sanitised(hub, payload):
    hub.user("alice")
    hub.create("alice", "tiny-chem", "generator")
    hub.push("alice", "alice/tiny-chem", tiny_files(**{"README.md": b"# Hi\n\n" + payload + b"\n"}))
    for url in ("/alice/tiny-chem", "/alice/tiny-chem/blob/main/README.md"):
        html = hub.client.get(url, headers=HTML).text
        body = html.split("<main>", 1)[1]
        assert "<script>alert" not in body and "javascript:" not in body
        assert "onerror" not in body and "<iframe" not in body


def test_open_redirects_are_refused(stocked):
    c = stocked.client
    page = c.get("/login", headers=HTML)
    r = c.post("/login", data={"username": "alice", "password": "correct horse battery", "csrf": _csrf(page.text),
                               "next": "//evil.example/"}, headers=HTML, follow_redirects=False)
    assert r.headers["location"] == "/"
