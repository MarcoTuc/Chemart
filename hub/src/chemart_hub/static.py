"""The Chemart Hub as a static site, built from a registry on GitHub.

A registry is a git repository of hub repos:

    repos/<namespace>/<name>/   one hub repo per folder: the same files a hub
                                repo holds (chemart.yaml, generator.py, preview.json,
                                README.md, or network.json)
    namespaces.yaml             optional: organisations (members, full name, bio)
                                and user profiles
    site.yaml                   optional: featured and hidden repos, site texts

People contribute by pull request. Three commands serve the registry:

    chemart-hub validate-pr     the check that runs on every pull request
    chemart-hub sync-official   refresh repos/chemart/ from the installed catalog
    chemart-hub build-static    turn the registry into a site for GitHub Pages

The site's pages are the live hub's own. The registry is loaded into a
throwaway hub (SQLite and a blob store), then every public page is requested
from the real app as an anonymous visitor and saved. Each git commit that
changed a repo's folder becomes a hub commit carrying the git commit's date
and message, so revision ids are the same on every build and a pinned
revision keeps working. The site holds:

    index.html, browse/, about/, new/, <ns>/, <ns>/<name>/{commits,tree,blob}/...
    <ns>/<name>/resolve/<commit>/<path>   the raw files, as the live hub serves them
    static/                     the hub's CSS and scripts
    hub.json                    marks the site as a static hub, for the chemart client
    api/v1/index.json           every repo, for search
    api/v1/repos/<ns>/<name>.json   one repo, its commits and their file lists
    404.html, .nojekyll
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlsplit

import yaml

from chemart.hub import _format
from chemart.hub._ids import RepoId, check_path

STATIC_FORMAT = "chemart-static-hub"
STATIC_VERSION = 1
OFFICIAL = "chemart"
REPOS = "repos"


# --------------------------------------------------------------------------
# Reading the registry from git
# --------------------------------------------------------------------------

def _git(registry: Path, *args: str, stdin: bytes | None = None) -> bytes:
    done = subprocess.run(["git", "-C", str(registry), *args], input=stdin,
                          capture_output=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {done.stderr.decode(errors='replace').strip()}")
    return done.stdout


def _utc(stamp: str) -> str:
    """A git ISO date, as the hub stores dates: UTC, second precision."""
    return datetime.fromisoformat(stamp).astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def repo_ids(registry: Path, rev: str = "HEAD") -> list[RepoId]:
    """The repos present at `rev`: every repos/<ns>/<name>/ folder with a file in it."""
    names = _git(registry, "ls-tree", "-r", "--name-only", rev, "--", f"{REPOS}/").decode().splitlines()
    found = {tuple(n.split("/")[1:3]) for n in names if n.count("/") >= 3}
    return [RepoId(ns, name) for ns, name in sorted(found)]


def _snapshot(registry: Path, rev: str, repo: RepoId, cat) -> dict[str, bytes]:
    """The files of one repo folder at `rev` (paths relative to the folder)."""
    prefix = f"{REPOS}/{repo.namespace}/{repo.name}/"
    listing = _git(registry, "ls-tree", "-r", "-z", rev, "--", prefix).split(b"\0")
    files: dict[str, bytes] = {}
    for item in filter(None, listing):
        meta, _, path = item.decode().partition("\t")
        _, kind, obj = meta.split()
        if kind != "blob":
            continue
        files[path[len(prefix):]] = cat(obj)
    return files


class _Cat:
    """`git cat-file --batch`, one process for many blobs."""

    def __init__(self, registry: Path):
        self.proc = subprocess.Popen(["git", "-C", str(registry), "cat-file", "--batch"],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    def __call__(self, obj: str) -> bytes:
        self.proc.stdin.write(obj.encode() + b"\n")
        self.proc.stdin.flush()
        header = self.proc.stdout.readline().split()
        size = int(header[2])
        data = self.proc.stdout.read(size)
        self.proc.stdout.read(1)                   # the newline after each object
        return data

    def close(self) -> None:
        self.proc.stdin.close()
        self.proc.wait()


@dataclass
class Version:
    """One state of a repo folder, from the git commit that produced it."""
    git: str
    created_at: str
    message: str
    files: dict[str, bytes]


def history(registry: Path, repo: RepoId, rev: str = "HEAD") -> list[Version]:
    """The states of a repo folder, oldest first; empty states (deletions) are skipped."""
    path = f"{REPOS}/{repo.namespace}/{repo.name}"
    log = _git(registry, "log", "--reverse", "--format=%H%x1f%cI%x1f%s", rev, "--", path).decode()
    cat = _Cat(registry)
    try:
        out = []
        for line in filter(None, log.splitlines()):
            sha, stamp, subject = line.split("\x1f", 2)
            files = _snapshot(registry, sha, repo, cat)
            if files:
                out.append(Version(sha, _utc(stamp), subject.strip() or "Update", files))
        return out
    finally:
        cat.close()


def _repo_type(files: dict[str, bytes]) -> str:
    hub: dict[str, Any] = {}
    if "chemart.yaml" in files:
        try:
            hub, _ = _format.parse_chemart_yaml(files["chemart.yaml"])
        except _format.FormatError:
            pass
    return hub.get("repo_type") or ("network" if "network.json" in files else "generator")


def _yaml(registry: Path, name: str) -> dict[str, Any]:
    path = registry / name
    if not path.is_file():
        return {}
    data = yaml.safe_load(path.read_text()) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{name} must be a mapping")
    return data


# --------------------------------------------------------------------------
# Loading the registry into a throwaway hub
# --------------------------------------------------------------------------

@dataclass
class Report:
    repos: list[str] = field(default_factory=list)
    commits: int = 0
    pages: int = 0
    files: int = 0
    problems: list[str] = field(default_factory=list)


def _principal(row):
    from chemart_hub import auth

    return auth.Principal(user=row, scope="write", via="token")


def load(registry: Path, settings, report: Report) -> None:
    """Fill a fresh hub database and blob store from the registry's git history."""
    import secrets

    from chemart_hub import db, service
    from chemart_hub.errors import ApiError
    from chemart_hub.storage import BlobStore

    c = db.connect(settings.db_path)
    store = BlobStore(settings.blobs_dir)
    people = _yaml(registry, "namespaces.yaml")
    orgs = people.get("orgs") or {}
    users = people.get("users") or {}
    try:
        # Organisations need an owner; a builder account owns them and then leaves.
        builder = service.create_user(c, "registry-builder", secrets.token_urlsafe(24))
        for repo in repo_ids(registry):
            ns = repo.namespace
            if service.account(c, ns) is not None:
                continue
            try:
                if ns == settings.official_namespace or ns in orgs:
                    service.create_org(c, ns, builder["name"], allow_reserved=ns == settings.official_namespace)
                else:
                    service.create_user(c, ns, secrets.token_urlsafe(24))
            except ApiError as err:
                report.problems.append(f"namespace {ns}: {err.message}")
                continue
            profile = orgs.get(ns) or users.get(ns) or {}
            with c:
                c.execute("UPDATE accounts SET fullname = ?, bio = ? WHERE name = ?",
                          (profile.get("fullname"), profile.get("bio"), ns))
        with c:
            c.execute("DELETE FROM org_members WHERE user_id = ?", (builder["id"],))
            c.execute("UPDATE accounts SET archived_at = ? WHERE id = ?", (db.now(), builder["id"]))

        for repo in repo_ids(registry):
            owner = service.account(c, repo.namespace)
            if owner is None:
                continue
            versions = history(registry, repo)
            if not versions:
                continue
            principal = _principal(owner)
            repo_type = _repo_type(versions[-1].files)
            try:
                row = service.create_repo(c, principal, repo.name, repo_type, repo.namespace)
            except ApiError as err:
                report.problems.append(f"{repo}: {err.message}")
                continue
            with c:                                   # born with its first commit, not with the build
                c.execute("UPDATE repos SET created_at = ? WHERE id = ?", (versions[0].created_at, row["id"]))
            for v in versions:
                operations = []
                for path, data in sorted(v.files.items()):
                    try:
                        check_path(path)
                    except ValueError:
                        continue
                    operations.append({"op": "add", "path": path, "sha256": store.put_bytes(data)})
                try:
                    result = service.commit(c, store, settings, principal, row, message=v.message,
                                            parent=row["head"], operations=operations, replace=True,
                                            created_at=v.created_at)
                except ApiError as err:
                    detail = "".join(f"; {p}" for p in err.problems)
                    report.problems.append(f"{repo} at {v.git[:10]}: {err.message}{detail}")
                    continue
                if not result["unchanged"]:
                    report.commits += 1
                row = service.get_repo(c, repo.namespace, repo.name)
            if row["head"]:
                report.repos.append(str(repo))

        site = _yaml(registry, "site.yaml")
        with c:
            for key, flag in (("featured", "featured"), ("hidden", "hidden")):
                for rid in site.get(key) or []:
                    ns, _, name = str(rid).partition("/")
                    c.execute(f"UPDATE repos SET {flag} = 1 WHERE name = ? AND owner_id = "
                              "(SELECT id FROM accounts WHERE name = ?)", (name, ns))
            from chemart_hub import curation

            for key, value in (site.get("texts") or {}).items():
                if key in curation.SITE_TEXT:
                    c.execute("INSERT OR REPLACE INTO site_text (key, value, updated_at, updated_by) "
                              "VALUES (?, ?, ?, ?)", (key, str(value), db.now(), "registry"))
    finally:
        c.close()


# --------------------------------------------------------------------------
# Rendering the pages
# --------------------------------------------------------------------------

_ATTR = re.compile(r'(\s(?:href|src|action))="(/[^"]*)"')


def _rewrite(html: str, base: str, dirs: set[str]) -> str:
    """Point site-absolute links at the project path, with the trailing slash
    GitHub Pages wants for folders (otherwise it answers with a redirect)."""
    def fix(m: re.Match) -> str:
        url = m.group(2)
        if url.startswith("//"):
            return m.group(0)
        path, rest = re.match(r"([^?#]*)(.*)", url, re.S).groups()
        if path in dirs and path != "/":
            path += "/"
        return f'{m.group(1)}="{base}{path}{rest}"'
    return _ATTR.sub(fix, html)


def _date_only(stamp: str | None) -> str:
    """The static site is built once, so relative times would go stale."""
    return f"on {stamp[:10]}" if stamp else ""


def _page_path(url: str) -> str:
    return url.split("?", 1)[0]


def build(registry: str | Path, out: str | Path, *, site_url: str, registry_repo: str,
          branch: str = "main", clean: bool = True) -> Report:
    """Render the registry at HEAD as a static hub in `out`; return what was built."""
    from fastapi.testclient import TestClient

    from chemart_hub import db, service
    from chemart_hub.app import create_app
    from chemart_hub.config import Settings
    from chemart_hub.routes import web
    from chemart_hub.storage import BlobStore

    registry, out = Path(registry).resolve(), Path(out).resolve()
    site_url = site_url.rstrip("/")
    base = urlsplit(site_url).path.rstrip("/")
    registry_url = f"https://github.com/{registry_repo}"
    report = Report()

    with tempfile.TemporaryDirectory(prefix="chemart-static-") as tmp:
        settings = Settings(data_dir=Path(tmp), public_url=site_url, allow_signup=False)
        app = create_app(settings)
        load(registry, settings, report)

        env = web.templates.env
        saved = (dict(env.globals), env.filters["ago"], web.PAGE)
        env.globals.update(STATIC=True, REGISTRY_URL=registry_url, REGISTRY_BRANCH=branch)
        env.filters["ago"] = _date_only
        web.PAGE = 1_000_000                         # one browse page with every repo
        try:
            if clean and out.exists():
                shutil.rmtree(out)
            out.mkdir(parents=True, exist_ok=True)
            client = TestClient(app, headers={"accept": "text/html"})
            c = db.connect(settings.db_path)
            store = BlobStore(settings.blobs_dir)
            try:
                _render(client, c, store, settings, out, base, report)
                _api(c, settings, out, registry_repo, branch)
            finally:
                c.close()
        finally:
            env.globals.clear()
            env.globals.update(saved[0])
            env.filters["ago"], web.PAGE = saved[1], saved[2]
    return report


def _render(client, c, store, settings, out: Path, base: str, report: Report) -> None:
    from chemart_hub import service

    pages = ["/", "/browse?sort=updated", "/about", "/new"]
    raw: list[tuple[str, str]] = []                  # (site path, sha256)
    namespaces: set[str] = set()
    rows = c.execute("SELECT r.*, a.name AS namespace FROM repos r JOIN accounts a ON a.id = r.owner_id "
                     "WHERE r.head IS NOT NULL AND r.archived_at IS NULL").fetchall()
    for row in rows:
        rid = f"/{row['namespace']}/{row['name']}"
        namespaces.add(row["namespace"])
        pages += [rid, f"{rid}/commits", f"{rid}/tree/main"]
        cid = row["head"]
        while cid:
            commit = c.execute("SELECT * FROM commits WHERE id = ?", (cid,)).fetchone()
            pages.append(f"{rid}/tree/{cid}")
            for path, meta in service.manifest(commit).items():
                pages.append(f"{rid}/blob/{cid}/{quote(path)}")
                raw.append((f"{rid}/resolve/{cid}/{path}", meta["sha256"]))
            cid = commit["parent"]
    pages += [f"/{ns}" for ns in sorted(namespaces)]

    dirs = {_page_path(p) for p in pages}
    for url in pages:
        response = client.get(url)
        if response.status_code != 200:
            report.problems.append(f"page {url}: HTTP {response.status_code}")
            continue
        path = _page_path(url)
        target = out / path.strip("/") / "index.html" if path != "/" else out / "index.html"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(_rewrite(response.text, base, dirs))
        report.pages += 1

    missing = client.get("/-/-")
    (out / "404.html").write_text(_rewrite(missing.text, base, dirs))

    for site_path, sha in raw:
        target = out / site_path.strip("/")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(store.get(sha))
        report.files += 1

    here = Path(__file__).resolve().parent
    shutil.copytree(here / "static", out / "static", dirs_exist_ok=True)
    (out / "static" / "pygments.css").write_text(client.get("/static/pygments.css").text)
    (out / ".nojekyll").write_text("")


def _api(c, settings, out: Path, registry_repo: str, branch: str) -> None:
    """The JSON the chemart client reads from a static hub."""
    from chemart_hub import db, service

    api = out / "api" / "v1"
    rows = c.execute("SELECT r.*, a.name AS namespace FROM repos r JOIN accounts a ON a.id = r.owner_id "
                     "WHERE r.head IS NOT NULL AND r.archived_at IS NULL ORDER BY a.name, r.name").fetchall()
    index = []
    for row in rows:
        info = service.repo_json(c, row, settings)
        fts = c.execute("SELECT repo, title, summary, tags, body FROM repos_fts WHERE rowid = ?",
                        (row["id"],)).fetchone()
        info["text"] = " ".join(v for v in (fts or []) if v)[:4000].lower()
        if not info["hidden"]:
            index.append(info)

        commits, cid = [], row["head"]
        while cid:
            commit = c.execute("SELECT * FROM commits WHERE id = ?", (cid,)).fetchone()
            commits.append({**service.commit_json(c, commit), "repo_type": row["repo_type"]})
            cid = commit["parent"]
        target = api / "repos" / row["namespace"] / f"{row['name']}.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        public = {k: v for k, v in info.items() if k != "text"}
        target.write_text(json.dumps({"repo": public, "commits": commits}, ensure_ascii=False))

    (api / "index.json").write_text(json.dumps({"repos": index}, ensure_ascii=False))
    (out / "hub.json").write_text(json.dumps({
        "format": STATIC_FORMAT, "version": STATIC_VERSION, "site": settings.public_url,
        "registry": registry_repo, "branch": branch, "built_at": db.now(),
    }, indent=2))


# --------------------------------------------------------------------------
# The official shelf
# --------------------------------------------------------------------------

def sync_official(registry: str | Path, only: list[str] | None = None) -> list[str]:
    """Write repos/chemart/<id>/ for every catalog chemistry (archived ones
    excluded) and remove folders of chemistries that left the catalog.
    `only` restricts the shelf to those ids (still excluding archived ones).
    Returns the ids whose files changed."""
    from chemart import catalog

    from chemart_hub.seed import official_files

    root = Path(registry) / REPOS / OFFICIAL
    wanted = {c.id for c in catalog.active(catalog.load()) if c.implemented}
    if only is not None:
        wanted &= set(only)
    changed = []
    for cid in sorted(wanted):
        folder = root / cid
        files = official_files(cid)
        current = {p.name: p.read_bytes() for p in folder.glob("*") if p.is_file()} if folder.is_dir() else {}
        if current == files:
            continue
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True)
        for name, data in files.items():
            (folder / name).write_bytes(data)
        changed.append(cid)
    for folder in sorted(root.glob("*")) if root.is_dir() else []:
        if folder.is_dir() and folder.name not in wanted:
            shutil.rmtree(folder)
            changed.append(folder.name)
    return changed


# --------------------------------------------------------------------------
# Starting a registry
# --------------------------------------------------------------------------

TEMPLATE = Path(__file__).resolve().parent / "registry_template"


def init_registry(folder: str | Path, *, owner: str, site_url: str, chemart_repo: str = "MarcoTuc/Chemart",
                  chemart_ref: str = "chemart-library", branch: str = "main",
                  official: bool = True) -> list[Path]:
    """Write a new registry into `folder`: README, contributing guide,
    namespaces.yaml (with `owner` as the first maintainer), site.yaml, the
    validate / publish / sync-official workflows and, with `official`, the
    official shelf. Refuses to overwrite existing files."""
    folder = Path(folder)
    values = {"{{OWNER}}": owner, "{{SITE_URL}}": site_url.rstrip("/"), "{{CHEMART_REPO}}": chemart_repo,
              "{{CHEMART_REF}}": chemart_ref, "{{BRANCH}}": branch}
    sources = [p for p in sorted(TEMPLATE.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
    targets = [folder / p.relative_to(TEMPLATE) for p in sources]
    clashes = [str(t) for t in targets if t.exists()]
    if clashes:
        raise FileExistsError(f"these files already exist: {', '.join(clashes)}")
    for source, target in zip(sources, targets):
        text = source.read_text()
        for key, value in values.items():
            text = text.replace(key, value)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
    if official:
        sync_official(folder)
    return targets


# --------------------------------------------------------------------------
# Checking a pull request
# --------------------------------------------------------------------------

def validate_pr(registry: str | Path, base: str, head: str, author: str, *,
                run_contract: bool = True) -> list[str]:
    """Problems with the changes from `base` to `head` made by GitHub user `author`.

    A contributor may change only repos/<their name>/ or the folder of an
    organisation that lists them as a member in namespaces.yaml, and every
    changed repo must be valid. Anything outside repos/ is for maintainers,
    the members of the official organisation. Generator repos are imported
    and run through `chemart check` when `run_contract` is true: the pull
    request's own code runs, so do this only on a throwaway CI runner.
    """
    from chemart.hub import check as contract_check
    from chemart_hub.service import builtin_ids

    registry = Path(registry)
    author = author.lower()
    people = _yaml(registry, "namespaces.yaml")
    orgs = {k.lower(): v or {} for k, v in (people.get("orgs") or {}).items()}

    def members(ns: str) -> set[str]:
        return {m.lower() for m in (orgs.get(ns, {}).get("members") or [])}

    maintainer = author in members(OFFICIAL)
    changed = _git(registry, "diff", "--name-only", f"{base}...{head}").decode().splitlines()
    problems: list[str] = []
    touched: set[RepoId] = set()
    for path in changed:
        parts = path.split("/")
        if parts[0] != REPOS or len(parts) < 4:
            if not maintainer:
                problems.append(f"{path}: only maintainers may change files outside {REPOS}/<namespace>/<name>/")
            continue
        touched.add(RepoId(parts[1], parts[2]))

    for repo in sorted(touched, key=str):
        ns = repo.namespace.lower()
        allowed = maintainer or ns == author or author in members(ns)
        if not allowed:
            problems.append(f"{repo}: you ({author}) may only change repos under {REPOS}/{author}/"
                            + (f" or an organisation that lists you in namespaces.yaml" if orgs else ""))
            continue
        cat = _Cat(registry)
        try:
            files = _snapshot(registry, head, repo, cat)
        finally:
            cat.close()
        if not files:
            continue                                  # the repo was deleted
        repo_type = _repo_type(files)
        try:
            card = _format.inspect(files, repo, repo_type, official_namespace=OFFICIAL,
                                   builtin_ids=builtin_ids())
        except _format.FormatError as err:
            problems += [f"{repo}: {p}" for p in err.problems]
            continue
        if run_contract and repo_type == "generator" and not card.builtin:
            with tempfile.TemporaryDirectory() as tmp:
                for path, data in files.items():
                    target = Path(tmp) / path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(data)
                problems += [f"{repo}: {p}" for p in contract_check(tmp)]
    return problems
