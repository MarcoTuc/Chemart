"""What the hub does, independent of HTTP: accounts, repos, commits, search.

Routes (API and web) are thin wrappers over these functions. Every commit
runs `chemart.hub._format.inspect` over the *whole* new snapshot, the same
check the client runs before uploading; nothing uploaded is ever executed.
"""

from __future__ import annotations

import json
import re
import sqlite3
from datetime import UTC, datetime, timedelta
from functools import lru_cache
from typing import Any

from chemart import catalog
from chemart.hub import _format
from chemart.hub._ids import (
    RESERVED_NAMESPACES, RepoId, check_name, check_namespace, check_path,
)

from chemart_hub import auth, db
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError, conflict, forbidden, invalid, not_found
from chemart_hub.storage import SHA256_RE, BlobStore, commit_id

SORTS = ("trending", "likes", "downloads", "updated", "created", "name")
PAGE_MAX = 100


@lru_cache(maxsize=1)
def builtin_ids() -> frozenset[str]:
    return frozenset(c.id for c in catalog.load())


# --------------------------------------------------------------------------
# Accounts
# --------------------------------------------------------------------------

def account(conn: sqlite3.Connection, name: str, *, include_archived: bool = False) -> sqlite3.Row | None:
    """A user or organisation by name. Archived accounts are invisible unless
    asked for; their names stay reserved until they are purged."""
    row = conn.execute("SELECT * FROM accounts WHERE name = ?", (name.lower(),)).fetchone()
    if row is not None and row["archived_at"] and not include_archived:
        return None
    return row


def create_user(conn: sqlite3.Connection, name: str, password: str, email: str | None = None,
                *, admin: bool = False, allow_reserved: bool = False) -> sqlite3.Row:
    name = name.strip().lower()
    try:
        check_namespace(name)
    except ValueError as err:
        raise invalid(str(err)) from None
    if name in RESERVED_NAMESPACES and not allow_reserved:
        raise invalid(f"the name {name!r} is reserved")
    auth.check_password_strength(password)
    if account(conn, name, include_archived=True) is not None:
        raise conflict(f"the name {name!r} is taken")
    with conn:
        conn.execute(
            "INSERT INTO accounts (name, email, password_hash, is_admin, created_at) VALUES (?, ?, ?, ?, ?)",
            (name, email, auth.hash_password(password), int(admin), db.now()),
        )
    return account(conn, name)


def create_org(conn: sqlite3.Connection, name: str, owner: str, *, allow_reserved: bool = False) -> sqlite3.Row:
    name = name.strip().lower()
    try:
        check_namespace(name)
    except ValueError as err:
        raise invalid(str(err)) from None
    if name in RESERVED_NAMESPACES and not allow_reserved:
        raise invalid(f"the name {name!r} is reserved")
    user = account(conn, owner)
    if user is None or user["is_org"]:
        raise not_found(f"user {owner!r}")
    if account(conn, name, include_archived=True) is not None:
        raise conflict(f"the name {name!r} is taken")
    with conn:
        cur = conn.execute("INSERT INTO accounts (name, is_org, created_at) VALUES (?, 1, ?)", (name, db.now()))
        conn.execute("INSERT INTO org_members (org_id, user_id, role) VALUES (?, ?, 'admin')",
                     (cur.lastrowid, user["id"]))
    return account(conn, name)


def add_member(conn: sqlite3.Connection, org: str, user: str) -> None:
    o, u = account(conn, org), account(conn, user)
    if o is None or not o["is_org"]:
        raise not_found(f"organisation {org!r}")
    if u is None or u["is_org"]:
        raise not_found(f"user {user!r}")
    with conn:
        conn.execute("INSERT OR IGNORE INTO org_members (org_id, user_id) VALUES (?, ?)", (o["id"], u["id"]))


# --------------------------------------------------------------------------
# Repos
# --------------------------------------------------------------------------

_REPO_SELECT = "SELECT r.*, a.name AS namespace FROM repos r JOIN accounts a ON a.id = r.owner_id"
#: Archived repos (deleted, waiting in the archive) are invisible everywhere
#: except the archive itself.
LIVE = "r.archived_at IS NULL"


def get_repo(conn: sqlite3.Connection, namespace: str, name: str) -> sqlite3.Row:
    row = conn.execute(f"{_REPO_SELECT} WHERE a.name = ? AND r.name = ? AND {LIVE}",
                       (namespace.lower(), name.lower())).fetchone()
    if row is None:
        raise not_found(f"repo {namespace}/{name}")
    return row


def repo_json(conn: sqlite3.Connection, repo: sqlite3.Row, settings: Settings) -> dict[str, Any]:
    tags = [r["tag"] for r in conn.execute("SELECT tag FROM repo_tags WHERE repo_id = ? ORDER BY rowid", (repo["id"],))]
    provides = [r["tag"] for r in conn.execute("SELECT tag FROM repo_provides WHERE repo_id = ? ORDER BY tag", (repo["id"],))]
    rid = f"{repo['namespace']}/{repo['name']}"
    return {
        "id": rid,
        "namespace": repo["namespace"],
        "name": repo["name"],
        "repo_type": repo["repo_type"],
        "head": repo["head"],
        "created_at": repo["created_at"],
        "updated_at": repo["updated_at"],
        "likes": repo["likes"],
        "downloads": repo["downloads"],
        "title": repo["title"],
        "summary": repo["summary"],
        "family": repo["family"],
        "kind": repo["kind"],
        "constructive": None if repo["constructive"] is None else bool(repo["constructive"]),
        "fidelity": repo["fidelity"],
        "license": repo["license"],
        "has_code": bool(repo["has_code"]),
        "builtin": repo["builtin"],
        "featured": bool(repo["featured"]),
        "hidden": bool(repo["hidden"]),
        "provides": provides,
        "tags": tags,
        "url": f"{settings.public_url}/{rid}",
    }


def create_repo(conn: sqlite3.Connection, principal: auth.Principal, name: str, repo_type: str,
                namespace: str | None = None, *, exist_ok: bool = False) -> sqlite3.Row:
    namespace = (namespace or principal.name).lower()
    try:
        check_name(name)
    except ValueError as err:
        raise invalid(str(err)) from None
    if repo_type not in _format.REPO_TYPES:
        raise invalid(f"repo_type must be one of {list(_format.REPO_TYPES)}")
    owner = account(conn, namespace)
    if owner is None:
        raise not_found(f"namespace {namespace!r}")
    if not auth.can_write(conn, principal, owner):
        raise forbidden(f"you cannot create repos under {namespace!r}")
    existing = conn.execute("SELECT * FROM repos WHERE owner_id = ? AND name = ? AND archived_at IS NULL",
                            (owner["id"], name)).fetchone()
    if existing is not None:
        if exist_ok and existing["repo_type"] == repo_type:
            return get_repo(conn, namespace, name)
        if exist_ok:
            raise conflict(f"{namespace}/{name} exists and is a {existing['repo_type']} repo")
        raise conflict(f"{namespace}/{name} already exists")
    stamp = db.now()
    with conn:
        conn.execute(
            "INSERT INTO repos (owner_id, name, repo_type, created_at, updated_at, title) VALUES (?, ?, ?, ?, ?, ?)",
            (owner["id"], name, repo_type, stamp, stamp, name),
        )
    return get_repo(conn, namespace, name)


def delete_repo(conn: sqlite3.Connection, principal: auth.Principal, repo: sqlite3.Row) -> None:
    """Deleting moves a repo to the archive: gone from the site and the API,
    but intact until a superadmin restores or purges it."""
    owner = conn.execute("SELECT * FROM accounts WHERE id = ?", (repo["owner_id"],)).fetchone()
    if not auth.can_write(conn, principal, owner):
        raise forbidden("you cannot delete this repo")
    with conn:
        conn.execute("UPDATE repos SET archived_at = ?, archived_by = ?, archived_reason = 'deleted' "
                     "WHERE id = ? AND archived_at IS NULL", (db.now(), principal.name, repo["id"]))


def require_write_access(conn: sqlite3.Connection, principal: auth.Principal, repo: sqlite3.Row) -> None:
    owner = conn.execute("SELECT * FROM accounts WHERE id = ?", (repo["owner_id"],)).fetchone()
    if not auth.can_write(conn, principal, owner):
        raise forbidden(f"you cannot push to {repo['namespace']}/{repo['name']}")


# --------------------------------------------------------------------------
# Revisions
# --------------------------------------------------------------------------

def resolve_revision(conn: sqlite3.Connection, repo: sqlite3.Row, revision: str) -> sqlite3.Row:
    """'main', a full commit id, or a unique prefix of at least 7 hex digits."""
    rid = f"{repo['namespace']}/{repo['name']}"
    if revision == "main":
        if repo["head"] is None:
            raise not_found(f"{rid} has no commits yet; revision 'main'")
        revision = repo["head"]
    if not re.fullmatch(r"[0-9a-f]{7,64}", revision or ""):
        raise ApiError(400, f"invalid revision {revision!r}: use 'main' or a commit id (>= 7 hex digits)")
    rows = conn.execute(
        "SELECT * FROM commits WHERE repo_id = ? AND id LIKE ? LIMIT 2", (repo["id"], revision + "%")
    ).fetchall()
    if not rows:
        raise not_found(f"revision {revision!r} of {rid}")
    if len(rows) > 1:
        raise ApiError(400, f"revision {revision!r} is ambiguous in {rid}; give more digits")
    return rows[0]


def manifest(commit: sqlite3.Row) -> dict[str, dict[str, Any]]:
    """path -> {sha256, size}"""
    return {p: {"sha256": sha, "size": size} for p, sha, size in json.loads(commit["manifest"])["files"]}


def commit_json(conn: sqlite3.Connection, commit: sqlite3.Row) -> dict[str, Any]:
    author = conn.execute("SELECT name FROM accounts WHERE id = ?", (commit["author_id"],)).fetchone()
    files = manifest(commit)
    return {
        "commit": commit["id"],
        "parent": commit["parent"],
        "message": commit["message"],
        "author": author["name"] if author else None,
        "created_at": commit["created_at"],
        "files": [{"path": p, **meta} for p, meta in sorted(files.items())],
    }


def history(conn: sqlite3.Connection, repo: sqlite3.Row, limit: int = 100) -> list[dict[str, Any]]:
    """Commits reachable from head, newest first."""
    out: list[dict[str, Any]] = []
    cid = repo["head"]
    while cid and len(out) < limit:
        row = conn.execute("SELECT * FROM commits WHERE id = ?", (cid,)).fetchone()
        if row is None:
            break
        entry = commit_json(conn, row)
        entry.pop("files")
        out.append(entry)
        cid = row["parent"]
    return out


def read_snapshot(store: BlobStore, files: dict[str, dict[str, Any]]) -> dict[str, bytes]:
    return {path: store.get(meta["sha256"]) for path, meta in files.items()}


def card_for(conn: sqlite3.Connection, store: BlobStore, repo: sqlite3.Row,
             commit: sqlite3.Row, settings: Settings) -> _format.RepoCard:
    return _format.inspect(
        read_snapshot(store, manifest(commit)),
        RepoId(repo["namespace"], repo["name"]),
        repo["repo_type"],
        official_namespace=settings.official_namespace,
    )


def commit(
    conn: sqlite3.Connection,
    store: BlobStore,
    settings: Settings,
    principal: auth.Principal,
    repo: sqlite3.Row,
    *,
    message: str,
    parent: str | None,
    operations: list[dict[str, Any]],
    replace: bool = False,
) -> dict[str, Any]:
    """Apply `operations` on top of `parent` (which must be the current head),
    validate the whole resulting snapshot, and move `main` to the new commit."""
    require_write_access(conn, principal, repo)
    rid = RepoId(repo["namespace"], repo["name"])
    if parent != repo["head"]:
        raise conflict(
            f"{rid} has moved: its head is {repo['head']}, not {parent}; fetch it and retry"
        )
    message = (message or "").strip()[:2000] or "Update"
    files: dict[str, dict[str, Any]] = {}
    if repo["head"] is not None and not replace:
        files = manifest(conn.execute("SELECT * FROM commits WHERE id = ?", (repo["head"],)).fetchone())

    problems: list[str] = []
    for op in operations:
        if not isinstance(op, dict):
            problems.append("each operation must be an object")
            continue
        kind, path = op.get("op"), op.get("path")
        try:
            check_path(path)
        except ValueError as err:
            problems.append(str(err))
            continue
        if kind == "add":
            sha = op.get("sha256")
            if not isinstance(sha, str) or not SHA256_RE.match(sha):
                problems.append(f"{path}: sha256 must be 64 lowercase hex digits")
            elif not store.has(sha):
                problems.append(f"{path}: blob {sha} was not uploaded")
            else:
                files[path] = {"sha256": sha, "size": store.size(sha)}
        elif kind == "delete":
            if files.pop(path, None) is None:
                problems.append(f"cannot delete {path}: no such file")
        else:
            problems.append(f"unknown operation {kind!r}; use 'add' or 'delete'")
    if problems:
        raise invalid("the commit is malformed", problems)

    if repo["head"] is not None:
        head_files = manifest(conn.execute("SELECT * FROM commits WHERE id = ?", (repo["head"],)).fetchone())
        if head_files == files:
            return {"commit": repo["head"], "unchanged": True}

    try:
        card = _format.inspect(
            read_snapshot(store, files), rid, repo["repo_type"],
            official_namespace=settings.official_namespace, builtin_ids=builtin_ids(),
        )
    except _format.FormatError as err:
        raise invalid(f"{rid}: the files do not form a valid {repo['repo_type']} repo", err.problems) from None

    stamp = db.now()
    triples = [[p, m["sha256"], m["size"]] for p, m in sorted(files.items())]
    cid = commit_id(repo=str(rid), parent=parent, files=triples, message=message,
                    author=principal.name, created_at=stamp)
    entry = card.entry
    with conn:
        conn.execute(
            "INSERT INTO commits (id, repo_id, parent, author_id, message, created_at, manifest) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (cid, repo["id"], parent, principal.id, message, stamp, json.dumps({"files": triples})),
        )
        moved = conn.execute(
            "UPDATE repos SET head = ?, updated_at = ?, title = ?, summary = ?, family = ?, kind = ?, "
            "constructive = ?, fidelity = ?, license = ?, has_code = ?, builtin = ? "
            "WHERE id = ? AND head IS ?",
            (
                cid, stamp, card.title or repo["name"], card.summary,
                (entry.family if entry else card.hub.get("family")),
                entry.kind if entry else None,
                None if entry is None else int(entry.constructive),
                entry.fidelity if entry else None,
                card.hub.get("license"), int(card.has_code), card.builtin,
                repo["id"], parent,
            ),
        ).rowcount
        if moved != 1:
            raise conflict(f"{rid} moved while committing; fetch it and retry")
        conn.execute("DELETE FROM repo_provides WHERE repo_id = ?", (repo["id"],))
        conn.execute("DELETE FROM repo_tags WHERE repo_id = ?", (repo["id"],))
        conn.executemany("INSERT INTO repo_provides (repo_id, tag) VALUES (?, ?)",
                         [(repo["id"], t) for t in sorted(set(card.provides))])
        conn.executemany("INSERT INTO repo_tags (repo_id, tag) VALUES (?, ?)",
                         [(repo["id"], t) for t in dict.fromkeys(card.hub.get("tags", []))])
        conn.execute("DELETE FROM repos_fts WHERE rowid = ?", (repo["id"],))
        body = " ".join(filter(None, [
            card.readme or "",
            (entry.intuition or "") if entry else "",
            " ".join(entry.aliases) if entry else "",
            (entry.origin or "") if entry else "",
        ]))
        conn.execute(
            "INSERT INTO repos_fts (rowid, repo, title, summary, tags, body) VALUES (?, ?, ?, ?, ?, ?)",
            (repo["id"], f"{rid.namespace} {rid.name.replace('-', ' ')}", card.title, card.summary,
             " ".join(card.hub.get("tags", [])), body),
        )
    return {"commit": cid, "unchanged": False}


# --------------------------------------------------------------------------
# Search, likes, downloads
# --------------------------------------------------------------------------

def _fts_query(text: str) -> str | None:
    words = re.findall(r"\w+", text.lower())
    return " ".join(f'"{w}"*' for w in words[:12]) or None


def _week_ago() -> tuple[str, str]:
    t = datetime.now(UTC) - timedelta(days=7)
    return t.strftime("%Y-%m-%dT%H:%M:%SZ"), t.strftime("%Y-%m-%d")


def search(
    conn: sqlite3.Connection,
    *,
    q: str | None = None,
    repo_type: str | None = None,
    family: str | None = None,
    kind: str | None = None,
    provides: list[str] | None = None,
    tag: str | None = None,
    author: str | None = None,
    fidelity: str | None = None,
    constructive: bool | None = None,
    has_code: bool | None = None,
    featured: bool | None = None,
    sort: str = "trending",
    limit: int = 30,
    offset: int = 0,
    include_empty: bool = False,
    include_hidden: bool = False,
) -> tuple[list[sqlite3.Row], int]:
    """Repos matching every given filter; returns (page, total).

    Hidden repos (unlisted by an admin) are left out unless `include_hidden`;
    they still load by id.
    """
    if sort not in SORTS:
        raise invalid(f"sort must be one of {list(SORTS)}")
    where, args = [LIVE], []
    if not include_empty:
        where.append("r.head IS NOT NULL")
    if not include_hidden:
        where.append("r.hidden = 0")
    if featured is not None:
        where.append("r.featured = ?")
        args.append(int(featured))
    for column, value in (("r.repo_type", repo_type), ("r.family", family), ("r.kind", kind),
                          ("r.fidelity", fidelity)):
        if value:
            where.append(f"{column} = ?")
            args.append(value)
    if author:
        where.append("a.name = ?")
        args.append(author.lower())
    if constructive is not None:
        where.append("r.constructive = ?")
        args.append(int(constructive))
    if has_code is not None:
        where.append("r.has_code = ?")
        args.append(int(has_code))
    for p in provides or []:
        where.append("EXISTS (SELECT 1 FROM repo_provides x WHERE x.repo_id = r.id AND x.tag = ?)")
        args.append(p)
    if tag:
        where.append("EXISTS (SELECT 1 FROM repo_tags x WHERE x.repo_id = r.id AND x.tag = ?)")
        args.append(tag)
    match = _fts_query(q) if q else None
    if q and match is None:
        return [], 0
    if match:
        where.append("r.id IN (SELECT rowid FROM repos_fts WHERE repos_fts MATCH ?)")
        args.append(match)
    clause = " WHERE " + " AND ".join(where)

    since, since_day = _week_ago()
    trending = (
        "(5 * (SELECT COUNT(*) FROM likes l WHERE l.repo_id = r.id AND l.created_at >= ?)"
        " + (SELECT COALESCE(SUM(d.count), 0) FROM downloads_daily d WHERE d.repo_id = r.id AND d.day >= ?))"
    )
    order = {
        "trending": f"{trending} DESC, r.likes DESC, r.updated_at DESC",
        "likes": "r.likes DESC, r.updated_at DESC",
        "downloads": "r.downloads DESC, r.updated_at DESC",
        "updated": "r.updated_at DESC",
        "created": "r.created_at DESC",
        "name": "a.name, r.name",
    }[sort]
    order_args = [since, since_day] if sort == "trending" else []
    total = conn.execute(f"SELECT COUNT(*) FROM repos r JOIN accounts a ON a.id = r.owner_id{clause}", args).fetchone()[0]
    limit = max(1, min(int(limit), PAGE_MAX))
    rows = conn.execute(
        f"{_REPO_SELECT}{clause} ORDER BY {order} LIMIT ? OFFSET ?",
        [*args, *order_args, limit, max(0, int(offset))],
    ).fetchall()
    return rows, total


def facet_counts(conn: sqlite3.Connection) -> dict[str, list[tuple[str, int]]]:
    """Counts for the browse sidebar, over non-empty repos."""
    out: dict[str, list[tuple[str, int]]] = {}
    for column in ("repo_type", "family", "fidelity", "license"):
        out[column] = [
            (r[0], r[1]) for r in conn.execute(
                f"SELECT {column}, COUNT(*) FROM repos WHERE head IS NOT NULL AND hidden = 0 "
                f"AND archived_at IS NULL AND {column} IS NOT NULL "
                f"GROUP BY {column} ORDER BY COUNT(*) DESC, {column}"
            )
        ]
    for table, key in (("repo_provides", "provides"), ("repo_tags", "tag")):
        out[key] = [
            (r[0], r[1]) for r in conn.execute(
                f"SELECT x.tag, COUNT(*) FROM {table} x JOIN repos r ON r.id = x.repo_id "
                f"WHERE r.head IS NOT NULL AND r.hidden = 0 AND r.archived_at IS NULL "
                f"GROUP BY x.tag ORDER BY COUNT(*) DESC, x.tag LIMIT 40"
            )
        ]
    return out


def set_like(conn: sqlite3.Connection, principal: auth.Principal, repo: sqlite3.Row, liked: bool) -> int:
    with conn:
        if liked:
            conn.execute("INSERT OR IGNORE INTO likes (user_id, repo_id, created_at) VALUES (?, ?, ?)",
                         (principal.id, repo["id"], db.now()))
        else:
            conn.execute("DELETE FROM likes WHERE user_id = ? AND repo_id = ?", (principal.id, repo["id"]))
        count = conn.execute("SELECT COUNT(*) FROM likes WHERE repo_id = ?", (repo["id"],)).fetchone()[0]
        conn.execute("UPDATE repos SET likes = ? WHERE id = ?", (count, repo["id"]))
    return count


def has_liked(conn: sqlite3.Connection, principal: auth.Principal | None, repo: sqlite3.Row) -> bool:
    if principal is None:
        return False
    return conn.execute("SELECT 1 FROM likes WHERE user_id = ? AND repo_id = ?",
                        (principal.id, repo["id"])).fetchone() is not None


def record_download(conn: sqlite3.Connection, repo: sqlite3.Row) -> None:
    day = datetime.now(UTC).strftime("%Y-%m-%d")
    with conn:
        conn.execute(
            "INSERT INTO downloads_daily (repo_id, day, count) VALUES (?, ?, 1) "
            "ON CONFLICT (repo_id, day) DO UPDATE SET count = count + 1",
            (repo["id"], day),
        )
        conn.execute("UPDATE repos SET downloads = downloads + 1 WHERE id = ?", (repo["id"],))
