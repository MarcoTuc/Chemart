"""The local cache of hub repos, and how a revision becomes files on disk.

    $CHEMART_HOME/hub/<hub host>/<ns>--<name>/
        refs/main                     the commit 'main' last resolved to
        manifests/<commit>.json       that commit's file list (path, sha256, size)
        snapshots/<commit>/<path>     the files, downloaded on demand

Every file is checked against its sha256 before it is moved into place, so a
file that exists in a snapshot is a verified copy. A full commit id that is
already cached never touches the network; 'main' and short ids do, unless the
hub is unreachable or CHEMART_HUB_OFFLINE=1, in which case the cache answers.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
import warnings
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote, urlsplit

from chemart.hub import _config, _http
from chemart.hub._ids import COMMIT_RE, RepoId, check_path


def _hub_folder() -> str:
    """One cache per hub, so alice/x on two hubs never mix."""
    netloc = urlsplit(_config.hub_url()).netloc or "hub"
    return re.sub(r"[^A-Za-z0-9._-]", "_", netloc)


def repo_dir(repo: RepoId) -> Path:
    return _config.home() / "hub" / _hub_folder() / repo.folder


def snapshot_dir(repo: RepoId, commit: str) -> Path:
    return repo_dir(repo) / "snapshots" / commit


def _atomic_write(target: Path, data: bytes) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=target.parent, prefix=".part-")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        os.replace(tmp, target)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _read_manifest(repo: RepoId, commit: str) -> dict[str, Any] | None:
    try:
        return json.loads((repo_dir(repo) / "manifests" / f"{commit}.json").read_text())
    except (OSError, ValueError):
        return None


def _cached_commit(repo: RepoId, revision: str) -> str | None:
    if revision == "main":
        try:
            commit = (repo_dir(repo) / "refs" / "main").read_text().strip()
        except OSError:
            return None
        return commit if COMMIT_RE.match(commit) else None
    folder = repo_dir(repo) / "manifests"
    matches = [p.stem for p in folder.glob(f"{revision}*.json")] if folder.is_dir() else []
    return matches[0] if len(matches) == 1 else None


def _validated(info: Any, repo: RepoId) -> dict[str, Any]:
    """Refuse a server answer that is malformed or names unsafe paths."""
    try:
        commit = info["commit"]
        if not COMMIT_RE.match(commit):
            raise ValueError(f"bad commit id {commit!r}")
        if info["repo_type"] not in ("generator", "network"):
            raise ValueError(f"bad repo type {info['repo_type']!r}")
        files = {}
        for f in info["files"]:
            path = check_path(f["path"])
            if not isinstance(f["sha256"], str) or len(f["sha256"]) != 64:
                raise ValueError(f"bad sha256 for {path}")
            files[path] = {"sha256": f["sha256"], "size": int(f["size"])}
    except (KeyError, TypeError, ValueError) as err:
        raise _http.HubError(f"the hub sent a malformed answer for {repo}: {err}") from None
    return {"commit": commit, "repo_type": info["repo_type"], "files": files}


def revision_info(repo: RepoId, revision: str) -> dict[str, Any]:
    """Resolve `revision` of `repo` to {commit, repo_type, files: {path: {sha256, size}}}."""
    if COMMIT_RE.match(revision):
        cached = _read_manifest(repo, revision)
        if cached is not None:
            return cached
    if _config.offline():
        return _from_cache(repo, revision, "CHEMART_HUB_OFFLINE is set")
    try:
        response = _http.request("GET", f"/api/repos/{repo.namespace}/{repo.name}/revision/{revision}")
    except _http.HubConnectionError as err:
        try:
            info = _from_cache(repo, revision, "the hub is unreachable")
        except _http.HubError:
            raise err from None
        warnings.warn(f"the Chemart Hub is unreachable; using the cached {repo}@{revision} "
                      f"({info['commit'][:12]})", stacklevel=4)
        return info
    info = _validated(response.json(), repo)
    _atomic_write(repo_dir(repo) / "manifests" / f"{info['commit']}.json", json.dumps(info).encode())
    if revision == "main":
        _atomic_write(repo_dir(repo) / "refs" / "main", info["commit"].encode())
    return info


def _from_cache(repo: RepoId, revision: str, why: str) -> dict[str, Any]:
    commit = _cached_commit(repo, revision)
    info = _read_manifest(repo, commit) if commit else None
    if info is None:
        raise _http.HubConnectionError(f"{why}, and {repo}@{revision} is not in the cache ({repo_dir(repo)})")
    return info


def fetch(repo: RepoId, info: dict[str, Any], paths: Iterable[str]) -> Path:
    """Make sure `paths` of this commit are on disk; return the snapshot folder."""
    commit = info["commit"]
    root = snapshot_dir(repo, commit)
    for path in paths:
        meta = info["files"].get(path)
        if meta is None:
            raise _http.RepoNotFoundError(f"{repo}@{commit[:12]} has no file {path!r}")
        target = root / check_path(path)
        if target.is_file():
            continue
        if _config.offline():
            raise _http.HubConnectionError(f"CHEMART_HUB_OFFLINE is set and {path} of {repo}@{commit[:12]} "
                                           "is not in the cache")
        response = _http.request("GET", f"/{repo.namespace}/{repo.name}/resolve/{commit}/{quote(path)}")
        digest = hashlib.sha256(response.body).hexdigest()
        if digest != meta["sha256"]:
            raise _http.HubError(f"{path} of {repo}@{commit[:12]} failed its integrity check "
                                 f"(sha256 {digest}, expected {meta['sha256']})")
        _atomic_write(target, response.body)
    return root


def trusted_before(repo: RepoId, commit: str) -> bool:
    """Has code from this commit already run on this machine? Records it if not."""
    record = repo_dir(repo) / "trusted.json"
    try:
        seen = set(json.loads(record.read_text()))
    except (OSError, ValueError):
        seen = set()
    if commit in seen:
        return True
    _atomic_write(record, json.dumps(sorted(seen | {commit})).encode())
    return False
