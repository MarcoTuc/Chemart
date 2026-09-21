"""The JSON API under /api, which `chemart.hub` (the client) talks to.

Uploads are two-step: PUT each file's bytes to /blobs/<sha256> (idempotent,
hash-checked), then POST a commit that names them. The commit carries the
head it was based on, and fails with 409 if someone pushed in between.
"""

from __future__ import annotations

import sqlite3
from typing import Any

from fastapi import APIRouter, Body, Depends, Query, Request

from chemart.hub import _config

from chemart_hub import auth, deps, service
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError, forbidden, invalid, unauthorized
from chemart_hub.storage import SHA256_RE, BlobMismatch, BlobStore, BlobTooLarge

router = APIRouter(prefix="/api")


def _bool(value: str | None) -> bool | None:
    if value is None or value == "":
        return None
    return value.lower() in {"1", "true", "yes"}


@router.get("/whoami")
def whoami(c: sqlite3.Connection = Depends(deps.conn),
           p: auth.Principal | None = Depends(deps.principal)) -> dict[str, Any]:
    if p is None:
        raise unauthorized()
    return {"name": p.name, "orgs": auth.namespaces_of(c, p.id), "scope": p.scope,
            "is_admin": bool(p.user["is_admin"])}


@router.get("/repos")
def list_repos(
    search: str | None = None,
    repo_type: str | None = None,
    family: str | None = None,
    kind: str | None = None,
    provides: str | None = Query(None, description="comma-separated; all must match"),
    tag: str | None = None,
    author: str | None = None,
    fidelity: str | None = None,
    constructive: str | None = None,
    has_code: str | None = None,
    sort: str = "trending",
    limit: int = 30,
    offset: int = 0,
    c: sqlite3.Connection = Depends(deps.conn),
    settings: Settings = Depends(deps.settings),
) -> dict[str, Any]:
    rows, total = service.search(
        c, q=search, repo_type=repo_type, family=family, kind=kind,
        provides=[p for p in (provides or "").split(",") if p], tag=tag, author=author,
        fidelity=fidelity, constructive=_bool(constructive), has_code=_bool(has_code),
        sort=sort, limit=limit, offset=offset,
    )
    return {"total": total, "repos": [service.repo_json(c, r, settings) for r in rows]}


@router.post("/repos", status_code=201)
def create_repo(
    body: dict[str, Any] = Body(...),
    c: sqlite3.Connection = Depends(deps.conn),
    p: auth.Principal = Depends(deps.writer),
    settings: Settings = Depends(deps.settings),
) -> dict[str, Any]:
    repo = service.create_repo(
        c, p, str(body.get("name", "")), str(body.get("repo_type", "")),
        body.get("namespace"), exist_ok=bool(body.get("exist_ok")),
    )
    return service.repo_json(c, repo, settings)


@router.get("/repos/{namespace}/{name}")
def repo_info(namespace: str, name: str, c: sqlite3.Connection = Depends(deps.conn),
              settings: Settings = Depends(deps.settings)) -> dict[str, Any]:
    return service.repo_json(c, service.get_repo(c, namespace, name), settings)


@router.delete("/repos/{namespace}/{name}")
def delete_repo(namespace: str, name: str, c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal = Depends(deps.writer)) -> dict[str, Any]:
    service.delete_repo(c, p, service.get_repo(c, namespace, name))
    return {"deleted": f"{namespace}/{name}", "archived": True}


@router.get("/repos/{namespace}/{name}/revision/{revision}")
def revision(namespace: str, name: str, revision: str,
             c: sqlite3.Connection = Depends(deps.conn)) -> dict[str, Any]:
    repo = service.get_repo(c, namespace, name)
    info = service.commit_json(c, service.resolve_revision(c, repo, revision))
    info["repo_type"] = repo["repo_type"]
    return info


@router.get("/repos/{namespace}/{name}/commits")
def commits(namespace: str, name: str, limit: int = 100,
            c: sqlite3.Connection = Depends(deps.conn)) -> dict[str, Any]:
    repo = service.get_repo(c, namespace, name)
    return {"commits": service.history(c, repo, max(1, min(limit, 1000)))}


@router.put("/repos/{namespace}/{name}/blobs/{sha256}")
async def put_blob(
    namespace: str, name: str, sha256: str, request: Request,
    c: sqlite3.Connection = Depends(deps.conn),
    p: auth.Principal = Depends(deps.writer),
    store: BlobStore = Depends(deps.store),
) -> dict[str, Any]:
    repo = service.get_repo(c, namespace, name)
    service.require_write_access(c, p, repo)
    if not SHA256_RE.match(sha256):
        raise ApiError(400, "the blob name must be the sha256 of its content (64 lowercase hex digits)")
    limit = _config.MAX_NETWORK
    declared = request.headers.get("content-length")
    if declared and declared.isdigit() and int(declared) > limit:
        raise ApiError(413, f"file exceeds the {limit}-byte limit")
    if store.has(sha256):
        store.touch(sha256)
        return {"sha256": sha256, "size": store.size(sha256), "existed": True}
    try:
        size = await store.put_stream(sha256, request.stream(), limit)
    except BlobTooLarge as err:
        raise ApiError(413, str(err)) from None
    except BlobMismatch as err:
        raise invalid(str(err)) from None
    return {"sha256": sha256, "size": size, "existed": False}


@router.post("/repos/{namespace}/{name}/commit/main")
def commit(
    namespace: str, name: str,
    body: dict[str, Any] = Body(...),
    c: sqlite3.Connection = Depends(deps.conn),
    p: auth.Principal = Depends(deps.writer),
    store: BlobStore = Depends(deps.store),
    settings: Settings = Depends(deps.settings),
) -> dict[str, Any]:
    repo = service.get_repo(c, namespace, name)
    operations = body.get("operations")
    if not isinstance(operations, list):
        raise invalid("'operations' must be a list of {op, path, sha256} objects")
    result = service.commit(
        c, store, settings, p, repo,
        message=str(body.get("message") or ""),
        parent=body.get("parent_commit"),
        operations=operations,
        replace=bool(body.get("replace")),
    )
    result["url"] = f"{settings.public_url}/{namespace}/{name}"
    return result


@router.post("/repos/{namespace}/{name}/like")
def like(namespace: str, name: str, c: sqlite3.Connection = Depends(deps.conn),
         p: auth.Principal = Depends(deps.writer)) -> dict[str, Any]:
    return {"likes": service.set_like(c, p, service.get_repo(c, namespace, name), True), "liked": True}


@router.delete("/repos/{namespace}/{name}/like")
def unlike(namespace: str, name: str, c: sqlite3.Connection = Depends(deps.conn),
           p: auth.Principal = Depends(deps.writer)) -> dict[str, Any]:
    return {"likes": service.set_like(c, p, service.get_repo(c, namespace, name), False), "liked": False}


@router.post("/tokens", status_code=201)
def create_token(
    body: dict[str, Any] = Body(...),
    c: sqlite3.Connection = Depends(deps.conn),
    p: auth.Principal = Depends(deps.writer),
) -> dict[str, Any]:
    if p.via != "session":
        raise forbidden("tokens are created from a logged-in browser session (Settings → Tokens)")
    scope = str(body.get("scope", "write"))
    return {"token": auth.create_token(c, p.id, str(body.get("name", "token")), scope), "scope": scope}
