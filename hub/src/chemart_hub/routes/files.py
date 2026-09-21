"""Raw file downloads: /{namespace}/{name}/resolve/{revision}/{path}.

Files at a full commit id never change, so they are served as immutable.
Code and text are always text/plain with nosniff, so a browser never runs
or renders an uploaded file.
"""

from __future__ import annotations

import sqlite3

from fastapi import APIRouter, Depends
from fastapi.responses import Response

from chemart.hub._ids import check_path

from chemart_hub import deps, service
from chemart_hub.errors import not_found
from chemart_hub.storage import BlobStore

router = APIRouter()

#: Fetching this file is what "loading the repo" means; it counts as a download.
PRIMARY = {"generator": "chemart.yaml", "network": "network.json"}


@router.get("/{namespace}/{name}/resolve/{revision}/{path:path}")
def resolve(namespace: str, name: str, revision: str, path: str,
            c: sqlite3.Connection = Depends(deps.conn),
            store: BlobStore = Depends(deps.store)) -> Response:
    repo = service.get_repo(c, namespace, name)
    commit = service.resolve_revision(c, repo, revision)
    try:
        check_path(path)
    except ValueError:
        raise not_found(f"file {path!r}") from None
    meta = service.manifest(commit).get(path)
    if meta is None:
        raise not_found(f"file {path!r} at {commit['id'][:12]}")
    if path == PRIMARY[repo["repo_type"]]:
        service.record_download(c, repo)
    media = "application/json" if path.endswith(".json") else "text/plain; charset=utf-8"
    pinned = revision == commit["id"]
    return Response(
        store.get(meta["sha256"]),
        media_type=media,
        headers={
            "ETag": f'"{meta["sha256"]}"',
            "X-Chemart-Commit": commit["id"],
            "X-Content-Type-Options": "nosniff",
            "Content-Security-Policy": "default-src 'none'; sandbox",
            "Cache-Control": "public, max-age=31536000, immutable" if pinned else "no-cache",
        },
    )
