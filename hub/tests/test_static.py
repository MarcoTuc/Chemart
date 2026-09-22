"""The static hub: a registry on GitHub, built into a site for GitHub Pages.

A registry is a git repository of repos/<ns>/<name>/ folders. These tests
build one with two official chemistries, a user's chemistry and a network
with two revisions, then check the built site, read it with the real client
over HTTP, and exercise the pull-request rules.
"""

from __future__ import annotations

import functools
import http.server
import json
import os
import subprocess
import threading
from pathlib import Path

import pytest

import chemart
from chemart import hub
from chemart.hub import _pr, _static
from chemart.hub._ids import RepoId
from chemart_hub import static

SITE = "http://127.0.0.1:{port}/hub"


def git(root: Path, *args: str, date: str | None = None, author: str = "Ada") -> None:
    env = {**os.environ, "GIT_AUTHOR_NAME": author, "GIT_AUTHOR_EMAIL": f"{author}@example.com",
           "GIT_COMMITTER_NAME": author, "GIT_COMMITTER_EMAIL": f"{author}@example.com"}
    if date:
        env.update(GIT_AUTHOR_DATE=date, GIT_COMMITTER_DATE=date)
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True, env=env)


def network_folder(root: Path, repo: str, seed: int) -> None:
    net = chemart.generate_network("brusselator", seed=seed)
    ns, name = repo.split("/")
    _pr.stage(RepoId(ns, name), hub.network_files(net, repo, title="A Brusselator snapshot",
                                                   tags=["oscillator"]), root)


@pytest.fixture(scope="module")
def registry(tmp_path_factory) -> Path:
    root = tmp_path_factory.mktemp("registry")
    git(root, "init", "-q", "-b", "main")
    (root / "namespaces.yaml").write_text("orgs:\n  chemart:\n    members: [maintainer]\n"
                                          "  lab:\n    members: [ada]\n")
    (root / "site.yaml").write_text("featured: [ada/tiny-chem]\ntexts:\n  home_title: A test hub\n")
    static.sync_official(root, only=["brusselator", "gamma", "tierra"])   # tierra is archived
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "Stock the official shelf", date="2026-09-20T10:00:00+00:00")

    work = tmp_path_factory.mktemp("work") / "tiny-chem"
    hub.new(work, "tiny-chem", "Tiny chemistry")
    _pr.stage(RepoId("ada", "tiny-chem"), hub.generator_files(work, "ada/tiny-chem"), root)
    network_folder(root, "ada/snapshot", seed=1)
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "Share tiny-chem and a snapshot", date="2026-09-21T09:00:00+00:00")

    network_folder(root, "ada/snapshot", seed=2)
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "Snapshot at seed 2", date="2026-09-21T11:30:00+00:00")
    return root


class _Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


@pytest.fixture(scope="module")
def site(registry, tmp_path_factory):
    """The built site, served at /hub/ like a GitHub project page."""
    serve = tmp_path_factory.mktemp("serve")
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(_Quiet, directory=str(serve)))
    url = SITE.format(port=server.server_address[1])
    report = static.build(registry, serve / "hub", site_url=url, registry_repo="acme/chemart-hub")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield {"out": serve / "hub", "url": url, "report": report}
    server.shutdown()


# --------------------------------------------------------------------------
# The build
# --------------------------------------------------------------------------

def test_the_build_holds_every_page_file_and_answer(site):
    out, report = site["out"], site["report"]
    assert report.problems == []
    assert sorted(report.repos) == ["ada/snapshot", "ada/tiny-chem", "chemart/brusselator", "chemart/gamma"]
    for page in ["index.html", "browse/index.html", "about/index.html", "new/index.html", "404.html",
                 "ada/index.html", "ada/tiny-chem/index.html", "ada/snapshot/commits/index.html",
                 "ada/snapshot/tree/main/index.html", "static/chemart.css", "static/pygments.css",
                 "static/static-browse.js", "hub.json", "api/v1/index.json", ".nojekyll"]:
        assert (out / page).is_file(), page
    meta = json.loads((out / "hub.json").read_text())
    assert meta["format"] == "chemart-static-hub" and meta["registry"] == "acme/chemart-hub"
    commits = json.loads((out / "api/v1/repos/ada/snapshot.json").read_text())["commits"]
    assert [c["message"] for c in commits] == ["Snapshot at seed 2", "Share tiny-chem and a snapshot"]
    assert commits[0]["parent"] == commits[1]["commit"] and commits[0]["created_at"] == "2026-09-21T11:30:00Z"
    for c in commits:
        assert (out / f"ada/snapshot/resolve/{c['commit']}/network.json").is_file()


def test_pages_are_static_and_link_under_the_project_path(site):
    home = (site["out"] / "index.html").read_text()
    assert "A test hub" in home and 'href="/hub/browse/' in home and 'href="/hub/static/chemart.css"' in home
    assert "Log in" not in home and "/login" not in home and "opens a pull request" in home
    repo = (site["out"] / "ada/tiny-chem/index.html").read_text()
    assert "https://github.com/acme/chemart-hub/tree/main/repos/ada/tiny-chem" in repo
    assert "/like" not in repo
    assert 'href="/hub/ada/"' in repo                  # folders get their trailing slash


def test_the_stock_line_splits_the_chemistries(site):
    """Under the search bar: three chemistries in stock, of which gamma and
    ada/tiny-chem are generators and the Brusselator a given network; the
    network repo ada/snapshot is a shared network."""
    home = (site["out"] / "index.html").read_text()
    assert "<b>3</b> chemistries in stock" in home
    assert "<b>2</b> generators" in home
    assert "<b>1</b> given network<" in home
    assert "<b>1</b> shared network<" in home
    index = {r["id"]: r for r in json.loads((site["out"] / "api/v1/index.json").read_text())["repos"]}
    assert index["chemart/brusselator"]["type"] == "given"
    assert index["chemart/gamma"]["type"] == "generator"
    assert "<b>0</b> gases" in home


def test_builds_are_reproducible(registry, site, tmp_path):
    again = tmp_path / "again"
    static.build(registry, again, site_url=site["url"], registry_repo="acme/chemart-hub")
    for path in (site["out"] / "api").rglob("*.json"):
        assert (again / path.relative_to(site["out"])).read_bytes() == path.read_bytes(), path


def test_the_live_hub_is_untouched_by_a_build(site):
    from chemart_hub.routes import web

    assert web.templates.env.globals["STATIC"] is False
    assert web.PAGE == 24


def test_archived_chemistries_stay_off_the_official_shelf(registry):
    assert (registry / "repos/chemart/brusselator").is_dir()
    assert not (registry / "repos/chemart/tierra").exists()


# --------------------------------------------------------------------------
# The client, reading the site over HTTP
# --------------------------------------------------------------------------

@pytest.fixture
def client_env(site, tmp_path, monkeypatch):
    monkeypatch.setenv("CHEMART_HUB_URL", site["url"])
    monkeypatch.setenv("CHEMART_HOME", str(tmp_path / "home"))
    _static._known.clear()
    yield site
    _static._known.clear()


def test_the_client_reads_a_static_hub(client_env):
    assert _static.is_static()
    assert [r["id"] for r in hub.search("brusselator")] == ["ada/snapshot", "chemart/brusselator"]
    assert [r["id"] for r in hub.search(author="ada", repo_type="network")] == ["ada/snapshot"]
    assert hub.repo_info("ada/snapshot")["title"] == "A Brusselator snapshot"

    assert chemart.load_network("ada/snapshot").seed == 2
    first = json.loads((client_env["out"] / "api/v1/repos/ada/snapshot.json").read_text())["commits"][-1]
    assert chemart.load_network(f"ada/snapshot@{first['commit'][:8]}").seed == 1

    official = chemart.generate_network("chemart/brusselator", seed=1)
    assert official.to_dict()["reactions"] == chemart.generate_network("brusselator", seed=1).to_dict()["reactions"]

    with pytest.raises(ValueError, match="trust_remote_code"):
        chemart.generate_network("ada/tiny-chem", seed=0)
    with pytest.warns(UserWarning):
        net = chemart.generate_network("ada/tiny-chem", seed=0, trust_remote_code=True)
    assert net.chemistry.startswith("ada/tiny-chem@")


def test_a_missing_repo_is_a_clean_error(client_env):
    with pytest.raises(hub.RepoNotFoundError):
        chemart.load_network("ada/no-such-thing")


def test_login_explains_there_is_none(client_env):
    with pytest.raises(hub.HubAuthError, match="static hub"):
        hub.login("chm_whatever")


def test_push_without_the_github_cli_writes_the_folder_and_the_steps(client_env, tmp_path, monkeypatch):
    monkeypatch.setattr(_pr, "gh_ready", lambda: False)
    monkeypatch.chdir(tmp_path)
    net = chemart.generate_network("brusselator", seed=3)
    result = hub.push_network(net, "bob/another-snapshot", title="Another")
    folder = Path(result["folder"])
    assert folder == tmp_path / "chemart-pull-request/repos/bob/another-snapshot"
    assert {p.name for p in folder.iterdir()} == {"network.json", "chemart.yaml", "README.md"}
    assert "acme/chemart-hub" in " ".join(result["steps"])


# --------------------------------------------------------------------------
# Pull requests
# --------------------------------------------------------------------------

@pytest.fixture
def clone(registry, tmp_path):
    root = tmp_path / "clone"
    subprocess.run(["git", "clone", "-q", str(registry), str(root)], check=True, capture_output=True)
    git(root, "checkout", "-q", "-b", "pr")
    return root


def _commit(root: Path, author: str) -> None:
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "change", author=author)


def test_a_contributor_may_change_their_own_namespace(clone):
    network_folder(clone, "bob/snap", seed=4)
    _commit(clone, "bob")
    assert static.validate_pr(clone, "main", "HEAD", "bob") == []


def test_but_not_someone_elses(clone):
    network_folder(clone, "ada/snapshot", seed=5)
    _commit(clone, "mallory")
    problems = static.validate_pr(clone, "main", "HEAD", "mallory")
    assert len(problems) == 1 and "ada/snapshot" in problems[0] and "mallory" in problems[0]


def test_organisation_members_may_change_its_repos(clone):
    network_folder(clone, "lab/snap", seed=6)
    _commit(clone, "ada")
    assert static.validate_pr(clone, "main", "HEAD", "ada") == []
    assert static.validate_pr(clone, "main", "HEAD", "bob") != []


def test_only_maintainers_change_files_outside_repos(clone):
    (clone / "site.yaml").write_text("featured: [mallory/spam]\n")
    _commit(clone, "mallory")
    assert any("outside" in p for p in static.validate_pr(clone, "main", "HEAD", "mallory"))
    assert static.validate_pr(clone, "main", "HEAD", "maintainer") == []


def test_an_invalid_repo_is_refused(clone):
    folder = clone / "repos/bob/broken"
    folder.mkdir(parents=True)
    (folder / "generator.py").write_text("def generate(p, rng):\n    return None\n")
    _commit(clone, "bob")
    assert static.validate_pr(clone, "main", "HEAD", "bob", run_contract=False) != []


def test_a_chemistry_must_pass_its_contract(clone, tmp_path):
    work = tmp_path / "chem"
    hub.new(work, "slow-chem", "Slow")
    files = hub.generator_files(work, "bob/slow-chem")
    files["generator.py"] = files["generator.py"].replace(b"def generate(p, rng):",
                                                          b"def generate(p, rng):\n    p = None")
    _pr.stage(RepoId("bob", "slow-chem"), files, clone)
    _commit(clone, "bob")
    assert static.validate_pr(clone, "main", "HEAD", "bob") != []


# --------------------------------------------------------------------------
# Starting a registry
# --------------------------------------------------------------------------

def test_init_registry_writes_the_template(tmp_path):
    written = static.init_registry(tmp_path, owner="alice", site_url="https://alice.github.io/chemart-hub/",
                                   official=False)
    names = {p.relative_to(tmp_path).as_posix() for p in written}
    assert {"README.md", "CONTRIBUTING.md", "namespaces.yaml", "site.yaml", ".gitignore",
            ".github/workflows/validate.yml", ".github/workflows/publish.yml",
            ".github/workflows/sync-official.yml"} <= names
    assert "members: [alice]" in (tmp_path / "namespaces.yaml").read_text()
    assert "https://alice.github.io/chemart-hub" in (tmp_path / "README.md").read_text()
    publish = (tmp_path / ".github/workflows/publish.yml").read_text()
    assert "{{" not in publish.replace("${{", "")
    with pytest.raises(FileExistsError):
        static.init_registry(tmp_path, owner="alice", site_url="https://x", official=False)
