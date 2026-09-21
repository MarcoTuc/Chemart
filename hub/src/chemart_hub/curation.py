"""Curating the hub by hand: what superadmins (and repo owners) can change
from the web interface.

Everything that changes a repo's *content* goes through `service.commit`, so
edits made in the browser are ordinary commits: validated by the same rules as
a push, recorded in the history, and undoable with `restore`. Everything a
superadmin does to other people's things is written to the admin log.

Nothing is destroyed by deleting. Deleting a repo or an account moves it to
the **archive**: it vanishes from the site and the API, but every row and
file is kept, and a superadmin can restore it exactly as it was. Only a
**purge** (from the archive) removes it: the database rows go, the stored
files that no remaining commit uses are deleted from disk, and the database
is compacted so no trace of the purged rows survives in its free pages.
"""

from __future__ import annotations

import json
import sqlite3
import time
from typing import Any

import yaml

from chemart.hub import _format
from chemart.hub._ids import check_path

from chemart_hub import auth, db, service
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError, conflict, forbidden, invalid, not_found
from chemart_hub.storage import BlobStore

# --------------------------------------------------------------------------
# Site text
# --------------------------------------------------------------------------

#: key -> (label, default, kind, help). kind: "line" (one line of plain text),
#: "text" (a paragraph of plain text) or "markdown" (rendered, then sanitised).
SITE_TEXT: dict[str, tuple[str, str, str, str]] = {
    "announcement": (
        "Announcement banner", "", "markdown",
        "Shown at the top of every page while it is not empty.",
    ),
    "home_title": (
        "Home page headline", "The one stop shop for artificial chemistry & artificial life", "line", "",
    ),
    "home_lede": (
        "Home page subtitle",
        "Browse chemistries and reaction networks that people share, and load any of them into "
        "your code in one line.", "text", "",
    ),
    "featured_title": ("Featured shelf heading", "Curated picks", "line", "Shown above featured repos."),
    "official_title": ("Official shelf heading", "The book chemistries", "line", ""),
    "official_note": (
        "Official shelf note", "official repos from Banzhaf & Yamamoto", "line",
        "Shown after the count, e.g. “98 official repos from Banzhaf & Yamamoto”.",
    ),
    "share_intro": (
        "Share page introduction",
        "There are two kinds of repo. A **chemistry** holds a generator: a catalog entry plus the "
        "Python that builds its reaction network, so people can run it with their own parameters. "
        "A **network** holds one reaction network as plain data, with no code.",
        "markdown", "",
    ),
    "about": (
        "About page",
        "# About Chemart Hub\n\nChemart Hub is the shared shelf of Chemart, the one stop shop "
        "for artificial chemistry and artificial life. Every repo here loads in one line of "
        "Python: `chemart.generate_network(\"namespace/name\")` for a chemistry, "
        "`chemart.load_network(\"namespace/name\")` for a network.\n",
        "markdown", "The page at /about.",
    ),
    "footer": (
        "Footer", "Chemart Hub · the one stop shop for artificial chemistry & artificial life", "line", "",
    ),
}


def site_texts(conn: sqlite3.Connection) -> dict[str, str]:
    """Every site text, with defaults for keys nobody has edited."""
    texts = {key: spec[1] for key, spec in SITE_TEXT.items()}
    texts.update({r["key"]: r["value"] for r in conn.execute("SELECT key, value FROM site_text")
                  if r["key"] in SITE_TEXT})
    return texts


def set_site_text(conn: sqlite3.Connection, actor: auth.Principal, key: str, value: str) -> bool:
    """Store one text; an empty value or the default restores the default. Returns True if it changed."""
    require_admin(actor)
    if key not in SITE_TEXT:
        raise invalid(f"unknown site text {key!r}")
    value = value.replace("\r\n", "\n").strip()
    if len(value) > 20000:
        raise invalid("site texts are limited to 20000 characters")
    before = site_texts(conn)[key]
    default = SITE_TEXT[key][1]
    with conn:
        if value == default or (not value and key != "announcement"):
            conn.execute("DELETE FROM site_text WHERE key = ?", (key,))
        else:
            conn.execute(
                "INSERT INTO site_text (key, value, updated_at, updated_by) VALUES (?, ?, ?, ?) "
                "ON CONFLICT (key) DO UPDATE SET value = excluded.value, updated_at = excluded.updated_at, "
                "updated_by = excluded.updated_by",
                (key, value, db.now(), actor.name),
            )
    changed = site_texts(conn)[key] != before
    if changed:
        log(conn, actor, "edit site text", key)
    return changed


# --------------------------------------------------------------------------
# The admin log
# --------------------------------------------------------------------------

def log(conn: sqlite3.Connection, actor: auth.Principal, action: str, target: str, detail: str = "") -> None:
    with conn:
        conn.execute("INSERT INTO admin_log (at, actor, action, target, detail) VALUES (?, ?, ?, ?, ?)",
                     (db.now(), actor.name, action, target, detail[:500]))


def recent_log(conn: sqlite3.Connection, limit: int = 200) -> list[sqlite3.Row]:
    return conn.execute("SELECT * FROM admin_log ORDER BY id DESC LIMIT ?", (limit,)).fetchall()


def require_admin(principal: auth.Principal | None) -> auth.Principal:
    if principal is None:
        raise ApiError(401, "log in as a superadmin")
    if not principal.is_admin:
        raise forbidden("only superadmins can do this")
    return principal


def _logged_if_not_own(conn, actor: auth.Principal, repo: sqlite3.Row, action: str, detail: str = "") -> None:
    """Admins acting on someone else's repo leave a trace; owners editing their own do not."""
    if actor.is_admin and repo["owner_id"] != actor.id:
        log(conn, actor, action, f"{repo['namespace']}/{repo['name']}", detail)


# --------------------------------------------------------------------------
# Accounts
# --------------------------------------------------------------------------

def _admin_count(conn: sqlite3.Connection) -> int:
    return conn.execute("SELECT COUNT(*) FROM accounts WHERE is_admin = 1 AND disabled = 0 AND is_org = 0 "
                        "AND archived_at IS NULL").fetchone()[0]


def _user(conn: sqlite3.Connection, name: str) -> sqlite3.Row:
    row = service.account(conn, name)
    if row is None:
        raise not_found(f"account {name!r}")
    return row


def set_admin(conn: sqlite3.Connection, actor: auth.Principal, name: str, flag: bool) -> None:
    require_admin(actor)
    user = _user(conn, name)
    if user["is_org"]:
        raise invalid("organisations cannot be superadmins")
    if not flag and user["is_admin"] and _admin_count(conn) <= 1:
        raise conflict("this is the last superadmin; promote someone else first")
    with conn:
        conn.execute("UPDATE accounts SET is_admin = ? WHERE id = ?", (int(flag), user["id"]))
    log(conn, actor, "grant superadmin" if flag else "revoke superadmin", user["name"])


def set_disabled(conn: sqlite3.Connection, actor: auth.Principal, name: str, flag: bool) -> None:
    """Suspend (or restore) an account: no log-in, no token works, repos stay."""
    require_admin(actor)
    user = _user(conn, name)
    if user["is_org"]:
        raise invalid("organisations do not log in; suspend their members instead")
    if flag and user["id"] == actor.id:
        raise conflict("you cannot suspend yourself")
    if flag and user["is_admin"] and _admin_count(conn) <= 1:
        raise conflict("this is the last superadmin")
    with conn:
        conn.execute("UPDATE accounts SET disabled = ? WHERE id = ?", (int(flag), user["id"]))
        if flag:
            conn.execute("DELETE FROM sessions WHERE user_id = ?", (user["id"],))
    log(conn, actor, "suspend" if flag else "unsuspend", user["name"])


def set_password(conn: sqlite3.Connection, actor: auth.Principal, name: str, password: str,
                 *, current: str | None = None) -> None:
    """Superadmins reset anyone's password; users change their own with the current one."""
    user = _user(conn, name)
    if user["is_org"]:
        raise invalid("organisations have no password")
    if actor.id == user["id"] and not actor.is_admin:
        if not auth.verify_password(current or "", user["password_hash"]):
            raise forbidden("the current password is wrong")
    elif actor.id != user["id"]:
        require_admin(actor)
    auth.check_password_strength(password)
    with conn:
        conn.execute("UPDATE accounts SET password_hash = ? WHERE id = ?", (auth.hash_password(password), user["id"]))
        if actor.id != user["id"]:
            conn.execute("DELETE FROM sessions WHERE user_id = ?", (user["id"],))
    if actor.id != user["id"]:
        log(conn, actor, "reset password", user["name"])


def update_profile(conn: sqlite3.Connection, actor: auth.Principal, name: str, *,
                   fullname: str, bio: str, email: str | None = None) -> None:
    account = _user(conn, name)
    if actor.id != account["id"] and not auth.can_write(conn, actor, account):
        raise forbidden("you cannot edit this profile")
    fullname, bio = fullname.strip()[:100], bio.replace("\r\n", "\n").strip()[:2000]
    with conn:
        conn.execute("UPDATE accounts SET fullname = ?, bio = ? WHERE id = ?",
                     (fullname or None, bio or None, account["id"]))
        if email is not None and not account["is_org"]:
            conn.execute("UPDATE accounts SET email = ? WHERE id = ?", (email.strip()[:200] or None, account["id"]))
    if actor.is_admin and actor.id != account["id"]:
        log(conn, actor, "edit profile", account["name"])


def delete_account(conn: sqlite3.Connection, actor: auth.Principal, name: str) -> int:
    """Archive a user or organisation and every repo it owns. Returns the repo count.

    The account disappears (no log-in, no tokens, no profile) but its name
    stays reserved, so nobody can take over its namespace, until it is purged.
    """
    require_admin(actor)
    account = _user(conn, name)
    if account["id"] == actor.id:
        raise conflict("you cannot delete your own account from here")
    if account["is_admin"] and _admin_count(conn) <= 1:
        raise conflict("this is the last superadmin")
    stamp = db.now()
    with conn:
        n = conn.execute("UPDATE repos SET archived_at = ?, archived_by = ?, archived_reason = 'account' "
                         "WHERE owner_id = ? AND archived_at IS NULL", (stamp, actor.name, account["id"])).rowcount
        conn.execute("UPDATE accounts SET archived_at = ?, archived_by = ? WHERE id = ?",
                     (stamp, actor.name, account["id"]))
        conn.execute("DELETE FROM sessions WHERE user_id = ?", (account["id"],))
    log(conn, actor, "archive account", account["name"], f"with {n} repos")
    return n


def set_membership(conn: sqlite3.Connection, actor: auth.Principal, org: str, user: str, member: bool) -> None:
    org_row, user_row = _user(conn, org), _user(conn, user)
    if not org_row["is_org"] or user_row["is_org"]:
        raise invalid("members join organisations, and members are users")
    if not auth.can_write(conn, actor, org_row):
        raise forbidden("you cannot manage this organisation")
    with conn:
        if member:
            conn.execute("INSERT OR IGNORE INTO org_members (org_id, user_id) VALUES (?, ?)", (org_row["id"], user_row["id"]))
        else:
            conn.execute("DELETE FROM org_members WHERE org_id = ? AND user_id = ?", (org_row["id"], user_row["id"]))
    if actor.is_admin:
        log(conn, actor, "add member" if member else "remove member", org_row["name"], user_row["name"])


def create_org(conn: sqlite3.Connection, actor: auth.Principal, name: str, owner: str) -> sqlite3.Row:
    require_admin(actor)
    row = service.create_org(conn, name, owner, allow_reserved=True)
    log(conn, actor, "create organisation", row["name"], f"owner {owner}")
    return row


# --------------------------------------------------------------------------
# Repos: curation
# --------------------------------------------------------------------------

def set_flags(conn: sqlite3.Connection, actor: auth.Principal, repo: sqlite3.Row, *,
              featured: bool | None = None, hidden: bool | None = None) -> None:
    require_admin(actor)
    rid = f"{repo['namespace']}/{repo['name']}"
    with conn:
        if featured is not None:
            conn.execute("UPDATE repos SET featured = ? WHERE id = ?", (int(featured), repo["id"]))
        if hidden is not None:
            conn.execute("UPDATE repos SET hidden = ? WHERE id = ?", (int(hidden), repo["id"]))
    if featured is not None:
        log(conn, actor, "feature" if featured else "unfeature", rid)
    if hidden is not None:
        log(conn, actor, "hide" if hidden else "unhide", rid)


def delete_repo(conn: sqlite3.Connection, actor: auth.Principal, repo: sqlite3.Row) -> None:
    """Move a repo to the archive (see the module docstring)."""
    service.delete_repo(conn, actor, repo)
    _logged_if_not_own(conn, actor, repo, "archive repo")


def transfer(conn: sqlite3.Connection, actor: auth.Principal, repo: sqlite3.Row, new_owner: str) -> sqlite3.Row:
    """Move a repo to another user or organisation, keeping its name and history."""
    require_admin(actor)
    target = _user(conn, new_owner)
    if conn.execute("SELECT 1 FROM repos WHERE owner_id = ? AND name = ? AND archived_at IS NULL",
                    (target["id"], repo["name"])).fetchone():
        raise conflict(f"{target['name']}/{repo['name']} already exists")
    with conn:
        conn.execute("UPDATE repos SET owner_id = ? WHERE id = ?", (target["id"], repo["id"]))
        conn.execute("UPDATE repos_fts SET repo = ? WHERE rowid = ?",
                     (f"{target['name']} {repo['name'].replace('-', ' ')}", repo["id"]))
    log(conn, actor, "transfer", f"{repo['namespace']}/{repo['name']}", f"to {target['name']}")
    return service.get_repo(conn, target["name"], repo["name"])


# --------------------------------------------------------------------------
# Repos: editing content in the browser
# --------------------------------------------------------------------------

def head_files(conn: sqlite3.Connection, store: BlobStore, repo: sqlite3.Row) -> dict[str, bytes]:
    if repo["head"] is None:
        return {}
    commit = conn.execute("SELECT * FROM commits WHERE id = ?", (repo["head"],)).fetchone()
    return service.read_snapshot(store, service.manifest(commit))


def edit_files(conn: sqlite3.Connection, store: BlobStore, settings: Settings, actor: auth.Principal,
               repo: sqlite3.Row, changes: dict[str, bytes | None], message: str) -> dict[str, Any]:
    """Commit `changes` (path -> new bytes, or None to delete) on top of head."""
    operations = []
    for path, data in changes.items():
        try:
            check_path(path)
        except ValueError as err:
            raise invalid(str(err)) from None
        if data is None:
            operations.append({"op": "delete", "path": path})
        else:
            operations.append({"op": "add", "path": path, "sha256": store.put_bytes(data)})
    result = service.commit(conn, store, settings, actor, repo, message=message,
                            parent=repo["head"], operations=operations)
    if not result.get("unchanged"):
        _logged_if_not_own(conn, actor, repo, "edit files", f"{', '.join(changes)}: {message}")
    return result


def text_from_form(text: str) -> bytes:
    """Browsers send CRLF from a textarea; files on the hub use LF."""
    text = text.replace("\r\n", "\n")
    if text and not text.endswith("\n"):
        text += "\n"
    return text.encode("utf-8")


class _NoAliases(yaml.SafeDumper):
    def ignore_aliases(self, data):  # the hub refuses YAML anchors and aliases
        return True


def update_metadata(conn: sqlite3.Connection, store: BlobStore, settings: Settings, actor: auth.Principal,
                    repo: sqlite3.Row, *, title: str, description: str, tags: str, license: str) -> dict[str, Any]:
    """Rewrite the `hub:` block of chemart.yaml (creating the file for a network repo)."""
    files = head_files(conn, store, repo)
    doc: dict[str, Any] = {}
    if "chemart.yaml" in files:
        loaded = _format.load_yaml(files["chemart.yaml"])
        doc = loaded if isinstance(loaded, dict) else {}
    hub = dict(doc.get("hub") or {})
    hub["repo_type"] = repo["repo_type"]
    values = {
        "title": " ".join(title.split()),
        "description": " ".join(description.split()),
        "license": license.strip(),
        "tags": [t.strip().lower() for t in tags.replace(",", " ").split() if t.strip()],
    }
    for key, value in values.items():
        if value:
            hub[key] = value
        else:
            hub.pop(key, None)
    new_doc = {"hub": hub, **{k: v for k, v in doc.items() if k != "hub"}}
    data = yaml.dump(new_doc, Dumper=_NoAliases, sort_keys=False, allow_unicode=True, width=100).encode()
    return edit_files(conn, store, settings, actor, repo, {"chemart.yaml": data}, "Edit title, description and tags")


def restore(conn: sqlite3.Connection, store: BlobStore, settings: Settings, actor: auth.Principal,
            repo: sqlite3.Row, commit_id: str) -> dict[str, Any]:
    """Make an old commit's files the new head (a new commit; nothing is erased)."""
    old = service.resolve_revision(conn, repo, commit_id)
    files = service.manifest(old)
    operations = [{"op": "add", "path": p, "sha256": m["sha256"]} for p, m in files.items()]
    result = service.commit(conn, store, settings, actor, repo, message=f"Restore {old['id'][:12]}",
                            parent=repo["head"], operations=operations, replace=True)
    if not result.get("unchanged"):
        _logged_if_not_own(conn, actor, repo, "restore", old["id"][:12])
    return result


# --------------------------------------------------------------------------
# The archive: restore, purge, and reclaiming disk space
# --------------------------------------------------------------------------

def _archived_repo(conn: sqlite3.Connection, repo_id: int) -> sqlite3.Row:
    row = conn.execute("SELECT r.*, a.name AS namespace FROM repos r JOIN accounts a ON a.id = r.owner_id "
                       "WHERE r.id = ? AND r.archived_at IS NOT NULL", (repo_id,)).fetchone()
    if row is None:
        raise not_found(f"archived repo #{repo_id}")
    return row


def _archived_account(conn: sqlite3.Connection, account_id: int) -> sqlite3.Row:
    row = conn.execute("SELECT * FROM accounts WHERE id = ? AND archived_at IS NOT NULL", (account_id,)).fetchone()
    if row is None:
        raise not_found(f"archived account #{account_id}")
    return row


def _manifests(conn: sqlite3.Connection, where: str = "1", args: tuple = ()) -> dict[int, set[str]]:
    """repo id -> the sha256 of every file in any of its commits."""
    out: dict[int, set[str]] = {}
    for row in conn.execute(f"SELECT c.repo_id, c.manifest FROM commits c JOIN repos r ON r.id = c.repo_id "
                            f"WHERE {where}", args):
        out.setdefault(row["repo_id"], set()).update(f[1] for f in json.loads(row["manifest"])["files"])
    return out


def archive(conn: sqlite3.Connection, store: BlobStore) -> dict[str, Any]:
    """What the archive holds, and how much disk a purge of each item frees."""
    live: set[str] = set().union(*_manifests(conn, "r.archived_at IS NULL").values())
    archived = _manifests(conn, "r.archived_at IS NOT NULL")

    def freed(shas: set[str]) -> int:
        return sum(store.size(s) for s in shas - live if store.has(s))

    repos = []
    for row in conn.execute(
            "SELECT r.*, a.name AS namespace, a.archived_at AS owner_archived, "
            "(SELECT COUNT(*) FROM commits c WHERE c.repo_id = r.id) AS n_commits "
            "FROM repos r JOIN accounts a ON a.id = r.owner_id WHERE r.archived_at IS NOT NULL "
            "ORDER BY r.archived_at DESC, r.id DESC"):
        taken = conn.execute("SELECT 1 FROM repos WHERE owner_id = ? AND name = ? AND archived_at IS NULL",
                             (row["owner_id"], row["name"])).fetchone() is not None
        repos.append({**dict(row), "id_text": f"{row['namespace']}/{row['name']}", "freed": freed(archived.get(row["id"], set())),
                      "name_taken": taken})
    accounts = []
    for row in conn.execute(
            "SELECT a.*, (SELECT COUNT(*) FROM repos r WHERE r.owner_id = a.id) AS n_repos "
            "FROM accounts a WHERE a.archived_at IS NOT NULL ORDER BY a.archived_at DESC"):
        shas: set[str] = set().union(*(archived.get(r["id"], set()) for r in conn.execute(
            "SELECT id FROM repos WHERE owner_id = ?", (row["id"],))))
        accounts.append({**dict(row), "freed": freed(shas)})
    blobs = list(store.blobs())
    everything: set[str] = set().union(*archived.values())
    referenced = live | everything
    orphans = [(sha, path) for sha, path in blobs if sha not in referenced]
    return {
        "repos": repos,
        "accounts": accounts,
        "disk": sum(path.stat().st_size for _, path in blobs),
        "blob_count": len(blobs),
        "freed_by_all": freed(everything) + sum(path.stat().st_size for _, path in orphans),
        "orphans": len(orphans) + len(store.partial_uploads()),
    }


def restore_repo(conn: sqlite3.Connection, actor: auth.Principal, repo_id: int) -> sqlite3.Row:
    require_admin(actor)
    repo = _archived_repo(conn, repo_id)
    rid = f"{repo['namespace']}/{repo['name']}"
    owner = conn.execute("SELECT * FROM accounts WHERE id = ?", (repo["owner_id"],)).fetchone()
    if owner["archived_at"]:
        raise conflict(f"{owner['name']} itself is archived; restore the account first")
    if conn.execute("SELECT 1 FROM repos WHERE owner_id = ? AND name = ? AND archived_at IS NULL",
                    (repo["owner_id"], repo["name"])).fetchone():
        raise conflict(f"a new {rid} exists; delete (or transfer) it before restoring this one")
    with conn:
        conn.execute("UPDATE repos SET archived_at = NULL, archived_by = NULL, archived_reason = NULL "
                     "WHERE id = ?", (repo["id"],))
    log(conn, actor, "restore repo", rid)
    return service.get_repo(conn, repo["namespace"], repo["name"])


def restore_account(conn: sqlite3.Connection, actor: auth.Principal, account_id: int) -> list[str]:
    """Bring back an account and the repos archived along with it; returns
    the repos that could not come back because their name was reused."""
    require_admin(actor)
    account = _archived_account(conn, account_id)
    skipped = []
    with conn:
        conn.execute("UPDATE accounts SET archived_at = NULL, archived_by = NULL WHERE id = ?", (account["id"],))
        for repo in conn.execute("SELECT * FROM repos WHERE owner_id = ? AND archived_reason = 'account' "
                                 "AND archived_at IS NOT NULL", (account["id"],)).fetchall():
            if conn.execute("SELECT 1 FROM repos WHERE owner_id = ? AND name = ? AND archived_at IS NULL",
                            (repo["owner_id"], repo["name"])).fetchone():
                skipped.append(f"{account['name']}/{repo['name']}")
                continue
            conn.execute("UPDATE repos SET archived_at = NULL, archived_by = NULL, archived_reason = NULL "
                         "WHERE id = ?", (repo["id"],))
    log(conn, actor, "restore account", account["name"], f"skipped {', '.join(skipped)}" if skipped else "")
    return skipped


def _secure(conn: sqlite3.Connection) -> None:
    """Zero the pages of deleted rows as they are freed (VACUUM then drops them)."""
    conn.execute("PRAGMA secure_delete = ON")


def _recount_likes(conn: sqlite3.Connection) -> None:
    conn.execute("UPDATE repos SET likes = (SELECT COUNT(*) FROM likes l WHERE l.repo_id = repos.id)")


def _files_of(conn: sqlite3.Connection, repo_ids: list[int]) -> set[str]:
    if not repo_ids:
        return set()
    marks = ",".join("?" * len(repo_ids))
    return set().union(*_manifests(conn, f"r.id IN ({marks})", tuple(repo_ids)).values())


def purge_repo(conn: sqlite3.Connection, store: BlobStore, actor: auth.Principal, repo_id: int) -> dict[str, int]:
    require_admin(actor)
    repo = _archived_repo(conn, repo_id)
    started, doomed = time.time(), _files_of(conn, [repo["id"]])
    _secure(conn)
    with conn:
        conn.execute("DELETE FROM repos_fts WHERE rowid = ?", (repo["id"],))
        conn.execute("DELETE FROM repos WHERE id = ?", (repo["id"],))
    log(conn, actor, "purge repo", f"{repo['namespace']}/{repo['name']}")
    return collect_garbage(conn, store, doomed=doomed, since=started)


def purge_account(conn: sqlite3.Connection, store: BlobStore, actor: auth.Principal, account_id: int) -> dict[str, int]:
    require_admin(actor)
    account = _archived_account(conn, account_id)
    repo_ids = [r["id"] for r in conn.execute("SELECT id FROM repos WHERE owner_id = ?", (account["id"],))]
    started, doomed = time.time(), _files_of(conn, repo_ids)
    _secure(conn)
    with conn:
        conn.executemany("DELETE FROM repos_fts WHERE rowid = ?", [(i,) for i in repo_ids])
        conn.execute("DELETE FROM accounts WHERE id = ?", (account["id"],))   # cascades to its repos
        _recount_likes(conn)                                                 # its likes went too
    log(conn, actor, "purge account", account["name"], f"with {len(repo_ids)} repos")
    return collect_garbage(conn, store, doomed=doomed, since=started)


def purge_all(conn: sqlite3.Connection, store: BlobStore, actor: auth.Principal) -> dict[str, int]:
    require_admin(actor)
    accounts = [r["id"] for r in conn.execute("SELECT id FROM accounts WHERE archived_at IS NOT NULL")]
    # every repo that goes: archived ones, and all repos of archived accounts
    going = [r["id"] for r in conn.execute(
        "SELECT r.id FROM repos r JOIN accounts a ON a.id = r.owner_id "
        "WHERE r.archived_at IS NOT NULL OR a.archived_at IS NOT NULL")]
    repos = [r["id"] for r in conn.execute(
        "SELECT r.id FROM repos r JOIN accounts a ON a.id = r.owner_id "
        "WHERE r.archived_at IS NOT NULL AND a.archived_at IS NULL")]
    started, doomed = time.time(), _files_of(conn, going)
    _secure(conn)
    with conn:
        conn.executemany("DELETE FROM repos_fts WHERE rowid = ?", [(i,) for i in going])
        conn.executemany("DELETE FROM repos WHERE id = ?", [(i,) for i in repos])
        conn.executemany("DELETE FROM accounts WHERE id = ?", [(i,) for i in accounts])
        _recount_likes(conn)
    log(conn, actor, "purge archive", "everything", f"{len(going)} repos, {len(accounts)} accounts")
    return collect_garbage(conn, store, doomed=doomed, since=started)


#: A blob uploaded this recently may belong to a push whose commit has not
#: arrived yet, so garbage collection leaves it alone.
UPLOAD_GRACE_SECONDS = 3600


def collect_garbage(conn: sqlite3.Connection, store: BlobStore, *, doomed: set[str] = frozenset(),
                    since: float | None = None, grace: float = UPLOAD_GRACE_SECONDS) -> dict[str, int]:
    """Delete stored files no commit uses any more, then compact the database.

    Blobs are shared between repos, so a blob is removed only when *no*
    remaining commit (live or archived) references it. An unreferenced blob
    is either the file of something just purged (`doomed`: removed, unless it
    was uploaded again after the purge began at `since`) or an orphan of an
    interrupted upload (removed once it is older than `grace`, so a push
    whose commit is still on its way is not broken). Afterwards the SQLite
    file is vacuumed and its write-ahead log truncated, so the purged rows do
    not survive in free pages either.
    """
    referenced: set[str] = set().union(*_manifests(conn).values())
    cutoff = time.time() - grace
    since = time.time() if since is None else since
    removed = freed = 0
    for sha, path in list(store.blobs()):
        if sha in referenced:
            continue
        mtime = path.stat().st_mtime
        if (sha in doomed and mtime < since) or mtime < cutoff:
            freed += store.remove(sha)
            removed += 1
    for tmp in store.partial_uploads():
        if tmp.stat().st_mtime < cutoff:
            freed += tmp.stat().st_size
            tmp.unlink(missing_ok=True)
            removed += 1
    conn.commit()
    conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    conn.execute("VACUUM")
    conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    return {"files_removed": removed, "bytes_freed": freed}
