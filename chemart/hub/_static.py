"""Reading a static hub: a site of plain files, such as GitHub Pages, built from
a registry repository by ``chemart-hub build-static``.

A static hub answers the client's read questions from files:

    hub.json                              {"format": "chemart-static-hub", "version": 1,
                                           "registry": "owner/repo", "branch": "main", ...}
    api/v1/index.json                     every listed repo, for search
    api/v1/repos/<ns>/<name>.json         {"repo": {...}, "commits": [newest first, with files]}
    <ns>/<name>/resolve/<commit>/<path>   raw files, at the same URL as on a live hub

It takes no writes: sharing opens a pull request on the registry instead
(`chemart.hub._pr`). Whether the configured hub is static is found out once
per hub URL, by asking for hub.json.
"""

from __future__ import annotations

import re
from typing import Any

from chemart.hub import _config, _http
from chemart.hub._ids import COMMIT_RE, RepoId

FORMAT = "chemart-static-hub"
SUPPORTED_VERSION = 1

_known: dict[str, dict[str, Any] | None] = {}
_HEX = re.compile(r"[0-9a-f]{7,64}")


def meta() -> dict[str, Any] | None:
    """The configured hub's hub.json if it is a static hub, else None.

    A hub that cannot be reached raises HubConnectionError and is asked again
    next time; a live hub (no hub.json) is remembered as such.
    """
    url = _config.hub_url()
    if url not in _known:
        try:
            data = _http.request("GET", "/hub.json", auth=False).json()
        except _http.RepoNotFoundError:
            data = None
        except ValueError:                       # a 200 that is not JSON: not ours
            data = None
        if isinstance(data, dict) and data.get("format") == FORMAT:
            if int(data.get("version", 0)) > SUPPORTED_VERSION:
                raise _http.HubError(f"the hub at {url} uses static format {data['version']}, "
                                     "newer than this chemart understands; upgrade chemart")
        else:
            data = None
        _known[url] = data
    return _known[url]


def is_static() -> bool:
    return meta() is not None


def _repo_file(repo: RepoId) -> dict[str, Any]:
    try:
        return _http.request("GET", f"/api/v1/repos/{repo.namespace}/{repo.name}.json", auth=False).json()
    except _http.RepoNotFoundError:
        raise _http.RepoNotFoundError(f"repo {repo} not found on {_config.hub_url()}", 404) from None


def revision(repo: RepoId, revision: str) -> dict[str, Any]:
    """What a live hub's /revision/ endpoint returns: {commit, repo_type, files: [...]}.

    `revision` is 'main', a full commit id, or a unique prefix of at least 7 hex digits.
    """
    commits = _repo_file(repo).get("commits") or []
    if not commits:
        raise _http.RepoNotFoundError(f"{repo} has no commits yet", 404)
    if revision == "main":
        found = commits[:1]
    elif COMMIT_RE.fullmatch(revision):
        found = [c for c in commits if c["commit"] == revision]
    elif _HEX.fullmatch(revision):
        found = [c for c in commits if c["commit"].startswith(revision)]
        if len(found) > 1:
            raise _http.HubError(f"revision {revision!r} of {repo} is ambiguous; give more digits", 400)
    else:
        raise _http.HubError(f"bad revision {revision!r}: use 'main', a commit id or a prefix of 7+ hex digits", 400)
    if not found:
        raise _http.RepoNotFoundError(f"{repo} has no revision {revision!r}", 404)
    c = found[0]
    return {"commit": c["commit"], "repo_type": c.get("repo_type"), "files": c["files"]}


def repo_info(repo: RepoId) -> dict[str, Any]:
    return _repo_file(repo)["repo"]


_SORT_KEYS = {"updated": "updated_at", "created": "created_at"}


def search(query: str | None = None, *, repo_type: str | None = None, family: str | None = None,
           kind: str | None = None, provides: list[str] | None = None, tag: str | None = None,
           author: str | None = None, sort: str = "updated", limit: int = 30) -> list[dict[str, Any]]:
    """Search the hub's index the way the live hub does, in memory.

    Every word of `query` must appear in the repo's id, title, summary, tags
    or text. There are no likes or downloads on a static hub, so the
    'trending', 'likes' and 'downloads' sorts fall back to 'updated'.
    """
    repos = _http.request("GET", "/api/v1/index.json", auth=False).json().get("repos") or []
    terms = (query or "").lower().split()
    out = []
    for r in repos:
        if repo_type and r.get("repo_type") != repo_type:
            continue
        if family and r.get("family") != family:
            continue
        if kind and r.get("kind") != kind:
            continue
        if tag and tag not in (r.get("tags") or []):
            continue
        if author and r.get("namespace") != author.lower():
            continue
        if provides and not set(provides) <= set(r.get("provides") or []):
            continue
        if terms:
            hay = " ".join([r.get("id") or "", r.get("title") or "", r.get("summary") or "",
                            " ".join(r.get("tags") or []), r.get("text") or ""]).lower()
            if not all(t in hay for t in terms):
                continue
        out.append({k: v for k, v in r.items() if k != "text"})
    if sort == "name":
        out.sort(key=lambda r: (r["name"], r["namespace"]))
    else:
        key = _SORT_KEYS.get(sort, "updated_at")
        out.sort(key=lambda r: r.get(key) or "", reverse=True)
    return out[:limit]
