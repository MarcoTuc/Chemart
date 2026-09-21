"""Passwords, API tokens, browser sessions, and who may write where.

Passwords are hashed with scrypt (standard library). API tokens look like
``chm_<random>`` and only their sha256 is stored, so a database leak leaks no
usable token. The API accepts a token as ``Authorization: Bearer``; the web
site uses a session cookie, and any state-changing request made with a
session must carry its CSRF token.
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import secrets
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta

from fastapi import Request

from chemart_hub import db
from chemart_hub.errors import ApiError, forbidden, unauthorized

SESSION_COOKIE = "chemart_session"
CSRF_HEADER = "x-csrf-token"
MIN_PASSWORD = 8
_SCRYPT = {"n": 2**14, "r": 8, "p": 1}


# --- passwords ------------------------------------------------------------

def hash_password(password: str) -> str:
    salt = secrets.token_bytes(16)
    key = hashlib.scrypt(password.encode(), salt=salt, dklen=32, **_SCRYPT)
    b64 = lambda b: base64.b64encode(b).decode()  # noqa: E731
    return f"scrypt${_SCRYPT['n']}${_SCRYPT['r']}${_SCRYPT['p']}${b64(salt)}${b64(key)}"


def verify_password(password: str, stored: str | None) -> bool:
    if not stored:
        return False
    try:
        scheme, n, r, p, salt, key = stored.split("$")
        if scheme != "scrypt":
            return False
        want = base64.b64decode(key)
        got = hashlib.scrypt(password.encode(), salt=base64.b64decode(salt),
                             n=int(n), r=int(r), p=int(p), dklen=len(want))
    except (ValueError, TypeError):
        return False
    return hmac.compare_digest(got, want)


def check_password_strength(password: str) -> None:
    if len(password) < MIN_PASSWORD:
        raise ApiError(422, f"password must have at least {MIN_PASSWORD} characters")


# --- tokens and sessions ----------------------------------------------------

def digest(secret: str) -> str:
    return hashlib.sha256(secret.encode()).hexdigest()


def create_token(conn: sqlite3.Connection, user_id: int, name: str, scope: str = "write") -> str:
    """Create an API token and return it; this is the only time it is visible."""
    if scope not in ("read", "write"):
        raise ApiError(422, "scope must be 'read' or 'write'")
    secret = "chm_" + secrets.token_urlsafe(32)
    with conn:
        conn.execute(
            "INSERT INTO tokens (user_id, name, token_hash, scope, created_at) VALUES (?, ?, ?, ?, ?)",
            (user_id, name[:64] or "token", digest(secret), scope, db.now()),
        )
    return secret


def create_session(conn: sqlite3.Connection, user_id: int, days: int) -> str:
    secret = secrets.token_urlsafe(32)
    expires = (datetime.now(UTC) + timedelta(days=days)).strftime("%Y-%m-%dT%H:%M:%SZ")
    with conn:
        conn.execute("DELETE FROM sessions WHERE expires_at < ?", (db.now(),))
        conn.execute(
            "INSERT INTO sessions (id_hash, user_id, csrf, expires_at) VALUES (?, ?, ?, ?)",
            (digest(secret), user_id, secrets.token_urlsafe(24), expires),
        )
    return secret


def end_session(conn: sqlite3.Connection, secret: str | None) -> None:
    if secret:
        with conn:
            conn.execute("DELETE FROM sessions WHERE id_hash = ?", (digest(secret),))


# --- who is asking ----------------------------------------------------------

@dataclass
class Principal:
    user: sqlite3.Row
    scope: str            # "read" | "write"
    via: str              # "token" | "session"
    csrf: str | None = None

    @property
    def name(self) -> str:
        return self.user["name"]

    @property
    def id(self) -> int:
        return self.user["id"]

    @property
    def is_admin(self) -> bool:
        """A superadmin: may change anything on the site."""
        return bool(self.user["is_admin"])


def identify(conn: sqlite3.Connection, request: Request) -> Principal | None:
    header = request.headers.get("authorization", "")
    if header:
        kind, _, secret = header.partition(" ")
        if kind.lower() != "bearer" or not secret.strip():
            raise unauthorized("malformed Authorization header; expected 'Bearer chm_...'")
        row = conn.execute(
            "SELECT t.id AS token_id, t.scope, a.* FROM tokens t JOIN accounts a ON a.id = t.user_id "
            "WHERE t.token_hash = ? AND t.revoked = 0 AND a.disabled = 0 AND a.archived_at IS NULL",
            (digest(secret.strip()),),
        ).fetchone()
        if row is None:
            raise unauthorized("invalid or revoked API token (or the account is suspended or deleted)")
        with conn:
            conn.execute("UPDATE tokens SET last_used_at = ? WHERE id = ?", (db.now(), row["token_id"]))
        return Principal(row, row["scope"], "token")
    secret = request.cookies.get(SESSION_COOKIE)
    if secret:
        row = conn.execute(
            "SELECT s.csrf, a.* FROM sessions s JOIN accounts a ON a.id = s.user_id "
            "WHERE s.id_hash = ? AND s.expires_at > ? AND a.disabled = 0 AND a.archived_at IS NULL",
            (digest(secret), db.now()),
        ).fetchone()
        if row is not None:
            return Principal(row, "write", "session", row["csrf"])
    return None


def require_writer(principal: Principal | None, request: Request) -> Principal:
    """The caller may change things: a write token, or a session plus its CSRF token."""
    if principal is None:
        raise unauthorized()
    if principal.scope != "write":
        raise forbidden("this token is read-only; create a write token")
    if principal.via == "session" and request.method not in ("GET", "HEAD", "OPTIONS"):
        sent = request.headers.get(CSRF_HEADER, "")
        if not (principal.csrf and hmac.compare_digest(sent, principal.csrf)):
            raise forbidden("missing or wrong CSRF token")
    return principal


def can_write(conn: sqlite3.Connection, principal: Principal | None, owner: sqlite3.Row) -> bool:
    """May `principal` push to repos owned by `owner` (a user or an org)?"""
    if principal is None:
        return False
    if owner["archived_at"]:
        return bool(principal.user["is_admin"])
    if principal.user["is_admin"] or principal.id == owner["id"]:
        return True
    return conn.execute(
        "SELECT 1 FROM org_members WHERE org_id = ? AND user_id = ?", (owner["id"], principal.id)
    ).fetchone() is not None


def namespaces_of(conn: sqlite3.Connection, user_id: int) -> list[str]:
    rows = conn.execute(
        "SELECT a.name FROM org_members m JOIN accounts a ON a.id = m.org_id "
        "WHERE m.user_id = ? AND a.archived_at IS NULL ORDER BY a.name",
        (user_id,),
    ).fetchall()
    return [r["name"] for r in rows]
