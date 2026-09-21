"""SQLite schema and connections.

One file, WAL mode, a connection per request. Users and organisations share
the `accounts` table because both own repos under their name. A repo's `head`
column is its only ref (`main`); it moves by compare-and-swap. The facet
columns (family, kind, ...) are copied from the head commit's card so that
browsing never has to parse files.
"""

from __future__ import annotations

import sqlite3
from datetime import UTC, datetime
from pathlib import Path

SCHEMA_VERSION = 3

_SCHEMA_V1 = """
CREATE TABLE accounts (
    id            INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE,
    is_org        INTEGER NOT NULL DEFAULT 0,
    email         TEXT,
    password_hash TEXT,
    is_admin      INTEGER NOT NULL DEFAULT 0,
    fullname      TEXT,
    created_at    TEXT NOT NULL
);
CREATE TABLE org_members (
    org_id  INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    user_id INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    role    TEXT NOT NULL DEFAULT 'member',
    PRIMARY KEY (org_id, user_id)
);
CREATE TABLE sessions (
    id_hash    TEXT PRIMARY KEY,
    user_id    INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    csrf       TEXT NOT NULL,
    expires_at TEXT NOT NULL
);
CREATE TABLE tokens (
    id           INTEGER PRIMARY KEY,
    user_id      INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    name         TEXT NOT NULL,
    token_hash   TEXT NOT NULL UNIQUE,
    scope        TEXT NOT NULL CHECK (scope IN ('read', 'write')),
    created_at   TEXT NOT NULL,
    last_used_at TEXT,
    revoked      INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE repos (
    id           INTEGER PRIMARY KEY,
    owner_id     INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    name         TEXT NOT NULL,
    repo_type    TEXT NOT NULL CHECK (repo_type IN ('generator', 'network')),
    head         TEXT,
    created_at   TEXT NOT NULL,
    updated_at   TEXT NOT NULL,
    likes        INTEGER NOT NULL DEFAULT 0,
    downloads    INTEGER NOT NULL DEFAULT 0,
    title        TEXT,
    summary      TEXT,
    family       TEXT,
    kind         TEXT,
    constructive INTEGER,
    fidelity     TEXT,
    license      TEXT,
    has_code     INTEGER NOT NULL DEFAULT 0,
    builtin      TEXT,
    UNIQUE (owner_id, name)
);
CREATE TABLE repo_provides (repo_id INTEGER NOT NULL REFERENCES repos(id) ON DELETE CASCADE, tag TEXT NOT NULL);
CREATE TABLE repo_tags     (repo_id INTEGER NOT NULL REFERENCES repos(id) ON DELETE CASCADE, tag TEXT NOT NULL);
CREATE INDEX repo_provides_tag ON repo_provides(tag, repo_id);
CREATE INDEX repo_tags_tag ON repo_tags(tag, repo_id);
CREATE TABLE commits (
    id         TEXT PRIMARY KEY,
    repo_id    INTEGER NOT NULL REFERENCES repos(id) ON DELETE CASCADE,
    parent     TEXT,
    author_id  INTEGER REFERENCES accounts(id) ON DELETE SET NULL,
    message    TEXT NOT NULL,
    created_at TEXT NOT NULL,
    manifest   TEXT NOT NULL
);
CREATE INDEX commits_repo ON commits(repo_id, created_at);
CREATE TABLE likes (
    user_id    INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    repo_id    INTEGER NOT NULL REFERENCES repos(id) ON DELETE CASCADE,
    created_at TEXT NOT NULL,
    PRIMARY KEY (user_id, repo_id)
);
CREATE TABLE downloads_daily (
    repo_id INTEGER NOT NULL REFERENCES repos(id) ON DELETE CASCADE,
    day     TEXT NOT NULL,
    count   INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (repo_id, day)
);
CREATE VIRTUAL TABLE repos_fts USING fts5(repo, title, summary, tags, body);
"""

#: v2: curation. Featured and hidden repos, suspended accounts, profile bios,
#: editable site text, and a log of every admin action.
_SCHEMA_V2 = """
ALTER TABLE repos ADD COLUMN featured INTEGER NOT NULL DEFAULT 0;
ALTER TABLE repos ADD COLUMN hidden INTEGER NOT NULL DEFAULT 0;
ALTER TABLE accounts ADD COLUMN disabled INTEGER NOT NULL DEFAULT 0;
ALTER TABLE accounts ADD COLUMN bio TEXT;
CREATE TABLE site_text (
    key        TEXT PRIMARY KEY,
    value      TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    updated_by TEXT
);
CREATE TABLE admin_log (
    id     INTEGER PRIMARY KEY,
    at     TEXT NOT NULL,
    actor  TEXT NOT NULL,
    action TEXT NOT NULL,
    target TEXT NOT NULL,
    detail TEXT
);
"""

#: v3: the archive. Deleting a repo or an account archives it; only a purge
#: removes it. A repo's name is unique among *live* repos only, so its owner
#: can reuse the name while the old repo waits in the archive. SQLite cannot
#: drop a table constraint, so `repos` is rebuilt (the documented procedure:
#: foreign keys off, copy, drop, rename; ids are kept, so every reference
#: stays valid).
_SCHEMA_V3 = """
PRAGMA foreign_keys = OFF;
BEGIN;
CREATE TABLE repos_v3 (
    id              INTEGER PRIMARY KEY,
    owner_id        INTEGER NOT NULL REFERENCES accounts(id) ON DELETE CASCADE,
    name            TEXT NOT NULL,
    repo_type       TEXT NOT NULL CHECK (repo_type IN ('generator', 'network')),
    head            TEXT,
    created_at      TEXT NOT NULL,
    updated_at      TEXT NOT NULL,
    likes           INTEGER NOT NULL DEFAULT 0,
    downloads       INTEGER NOT NULL DEFAULT 0,
    title           TEXT,
    summary         TEXT,
    family          TEXT,
    kind            TEXT,
    constructive    INTEGER,
    fidelity        TEXT,
    license         TEXT,
    has_code        INTEGER NOT NULL DEFAULT 0,
    builtin         TEXT,
    featured        INTEGER NOT NULL DEFAULT 0,
    hidden          INTEGER NOT NULL DEFAULT 0,
    archived_at     TEXT,
    archived_by     TEXT,
    archived_reason TEXT
);
INSERT INTO repos_v3 (id, owner_id, name, repo_type, head, created_at, updated_at, likes, downloads,
                      title, summary, family, kind, constructive, fidelity, license, has_code,
                      builtin, featured, hidden)
    SELECT id, owner_id, name, repo_type, head, created_at, updated_at, likes, downloads,
           title, summary, family, kind, constructive, fidelity, license, has_code,
           builtin, featured, hidden FROM repos;
DROP TABLE repos;
ALTER TABLE repos_v3 RENAME TO repos;
CREATE UNIQUE INDEX repos_live_name ON repos(owner_id, name) WHERE archived_at IS NULL;
CREATE INDEX repos_archived ON repos(archived_at);
ALTER TABLE accounts ADD COLUMN archived_at TEXT;
ALTER TABLE accounts ADD COLUMN archived_by TEXT;
COMMIT;
PRAGMA foreign_keys = ON;
"""


def now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path, check_same_thread=False, timeout=10)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA busy_timeout = 10000")
    return conn


def migrate(path: Path) -> None:
    """Create or upgrade the database at `path`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect(path)
    try:
        conn.execute("PRAGMA journal_mode = WAL")
        version = conn.execute("PRAGMA user_version").fetchone()[0]
        if version > SCHEMA_VERSION:
            raise RuntimeError(f"database schema v{version} is newer than this server (v{SCHEMA_VERSION})")
        if version < 1:
            conn.executescript(_SCHEMA_V1)
            conn.execute("PRAGMA user_version = 1")
        if version < 2:
            conn.executescript(_SCHEMA_V2)
            conn.execute("PRAGMA user_version = 2")
        if version < 3:
            conn.executescript(_SCHEMA_V3)
            problems = conn.execute("PRAGMA foreign_key_check").fetchall()
            if problems:
                raise RuntimeError(f"migration to v3 left dangling references: {problems[:5]}")
            conn.execute("PRAGMA user_version = 3")
        conn.commit()
    finally:
        conn.close()
