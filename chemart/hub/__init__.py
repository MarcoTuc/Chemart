"""The Chemart Hub client: share chemistries and networks, load other people's.

Hub ids are ``namespace/name``, optionally pinned with ``@revision``. They work
wherever a catalog id does:

    import chemart

    chemart.describe_chemistry("alice/my-chem")            # metadata only; runs nothing
    net = chemart.generate_network("alice/my-chem", seed=0,
                                   trust_remote_code=True, revision="3f2a9c1")
    net = chemart.load_network("bob/ecoli-core")            # a shared network: plain data
    net.push_to_hub("me/my-snapshot")                       # share a network
    chemart.hub.push_generator("./my-chem")                 # share a chemistry (code + entry)

A generator repo's code runs on your machine, so it needs
``trust_remote_code=True``; pin ``revision`` to the commit you read. Official
repos under ``chemart/`` are the built-in catalog and need no trust.

A hub is either a live server (``chemart-hub serve``) or a static site built
from a GitHub registry (``chemart-hub build-static``). Reading works the same
on both; on a static hub, sharing opens a pull request on the registry.

Configuration: CHEMART_HUB_URL, CHEMART_HOME, CHEMART_HUB_TOKEN,
CHEMART_HUB_OFFLINE (see `chemart.hub._config`).
"""

from __future__ import annotations

import getpass
import json
from pathlib import Path
from typing import Any

import yaml

from chemart.hub import _cache, _config, _format, _http, _pr, _remote_code, _static
from chemart.hub._http import (
    HubAuthError, HubConflictError, HubConnectionError, HubError,
    HubValidationError, RepoNotFoundError,
)
from chemart.hub._ids import RepoId, check_path, parse_repo_id
from chemart.hub._remote_code import RemoteChemistry
from chemart.network import Network

__all__ = [
    "HubAuthError", "HubConflictError", "HubConnectionError", "HubError",
    "HubValidationError", "RemoteChemistry", "RepoNotFoundError",
    "check", "create_repo", "delete_repo", "load_network", "login", "logout",
    "new", "push_generator", "push_network", "repo_info", "resolve_chemistry",
    "search", "snapshot_download", "upload_files", "whoami",
    "generator_files", "network_files",
]

#: Repos in this namespace may stand for built-in catalog entries.
OFFICIAL_NAMESPACE = "chemart"


# --------------------------------------------------------------------------
# Accounts
# --------------------------------------------------------------------------

def whoami(token: str | None = None) -> dict[str, Any]:
    """Who the hub thinks you are: {name, orgs, scope}.

    On a static hub your identity is your GitHub account (via the GitHub CLI)."""
    if _static.is_static():
        return {"name": _pr.github_login(), "orgs": [], "scope": "pull-request"}
    return _http.request("GET", "/api/whoami", token=token).json()


def login(token: str | None = None) -> str:
    """Save an API token for the configured hub (create one at /settings/tokens)."""
    url = _config.hub_url()
    if _static.is_static():
        raise HubAuthError(f"{url} is a static hub built from a GitHub registry: there is no "
                           "hub login. Sharing opens a pull request; sign in to GitHub with "
                           "`gh auth login`.")
    if token is None:
        token = getpass.getpass(f"Paste a token from {url}/settings/tokens: ").strip()
    if not token:
        raise HubAuthError("no token given")
    me = whoami(token)
    _config.save_token(token, url)
    return me["name"]


def logout() -> None:
    """Forget the saved token for the configured hub."""
    _config.save_token(None, _config.hub_url())


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------

def _repo(repo_id: str, revision: str | None = None) -> tuple[RepoId, str]:
    return parse_repo_id(repo_id, revision)


def repo_info(repo_id: str) -> dict[str, Any]:
    """The hub's metadata for a repo: type, head commit, title, facets, likes."""
    repo, _ = _repo(repo_id)
    if _static.is_static():
        return _static.repo_info(repo)
    return _http.request("GET", f"/api/repos/{repo.namespace}/{repo.name}").json()


def search(query: str | None = None, *, repo_type: str | None = None, family: str | None = None,
           kind: str | None = None, provides: list[str] | None = None, tag: str | None = None,
           author: str | None = None, sort: str = "trending", limit: int = 30) -> list[dict[str, Any]]:
    """Find repos. `provides` must all match, e.g. ["rate-constants", "energies"]."""
    from urllib.parse import urlencode

    if _static.is_static():
        return _static.search(query, repo_type=repo_type, family=family, kind=kind, provides=provides,
                              tag=tag, author=author, sort=sort, limit=limit)
    params = {"search": query, "repo_type": repo_type, "family": family, "kind": kind,
              "provides": ",".join(provides or []) or None, "tag": tag, "author": author,
              "sort": sort, "limit": limit}
    qs = urlencode({k: v for k, v in params.items() if v is not None})
    return _http.request("GET", f"/api/repos?{qs}").json()["repos"]


def snapshot_download(repo_id: str, revision: str | None = None) -> Path:
    """Download every file of a repo at a revision; return the local folder."""
    repo, rev = _repo(repo_id, revision)
    info = _cache.revision_info(repo, rev)
    return _cache.fetch(repo, info, info["files"])


def resolve_chemistry(chemistry_id: str, revision: str | None = None):
    """The catalog entry of a generator repo (what `generate_network` runs).

    Reads chemart.yaml only. For an official builtin-backed repo this is the
    built-in entry, so ``chemart/brusselator`` and ``brusselator`` agree.
    """
    repo, rev = _repo(chemistry_id, revision)
    info = _cache.revision_info(repo, rev)
    if info["repo_type"] != "generator":
        raise ValueError(f"{repo} is a network repo, not a chemistry: load it with "
                         f"chemart.load_network({str(repo)!r})")
    folder = _cache.fetch(repo, info, ["chemart.yaml"])
    hub, entry = _format.parse_chemart_yaml((folder / "chemart.yaml").read_bytes())
    if entry is None:
        raise HubError(f"{repo}@{info['commit'][:12]}: chemart.yaml has no chemistry entry")
    builtin = hub.get("builtin")
    if builtin and repo.namespace == OFFICIAL_NAMESPACE:
        from chemart import api

        try:
            return api._entry(builtin)
        except ValueError:
            raise HubError(f"{repo} is the built-in {builtin!r}, which this Chemart does not have; "
                           "upgrade chemart") from None
    return RemoteChemistry.from_entry(entry, repo=str(repo), commit=info["commit"],
                                      revision=rev, hub=hub, manifest=info)


def load_network(repo_id: str, revision: str | None = None) -> Network:
    """Load a shared reaction network (a network repo). It is plain data: no trust needed."""
    repo, rev = _repo(repo_id, revision)
    info = _cache.revision_info(repo, rev)
    if info["repo_type"] != "network":
        raise ValueError(f"{repo} is a chemistry (generator repo): run it with "
                         f"chemart.generate_network({str(repo)!r}, ...)")
    folder = _cache.fetch(repo, info, ["network.json"])
    return Network.from_dict(json.loads((folder / "network.json").read_text()))


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

def create_repo(repo_id: str, repo_type: str, *, exist_ok: bool = True) -> dict[str, Any]:
    repo, _ = _repo(repo_id)
    return _http.request("POST", "/api/repos", json_body={
        "name": repo.name, "namespace": repo.namespace, "repo_type": repo_type, "exist_ok": exist_ok,
    }).json()


def delete_repo(repo_id: str) -> None:
    """Delete a repo you own. The hub moves it to its archive: it is gone from
    the site at once, and a superadmin can restore it until it is purged."""
    repo, _ = _repo(repo_id)
    _http.request("DELETE", f"/api/repos/{repo.namespace}/{repo.name}")


def upload_files(repo_id: str, files: dict[str, bytes], *, repo_type: str,
                 message: str | None = None) -> dict[str, Any]:
    """Make `files` the new content of a repo (created if missing).

    The files are checked locally with the server's own rules first. Returns
    {commit, url, unchanged}; `unchanged` is True if nothing differed. On a
    static hub this opens a pull request on the registry instead and returns
    {pull_request, url, branch} (or, without the GitHub CLI, {folder, steps}).
    """
    import hashlib

    repo, _ = _repo(repo_id)
    try:
        _format.inspect(files, repo, repo_type, official_namespace=OFFICIAL_NAMESPACE)
    except _format.FormatError as err:
        raise HubValidationError(f"{repo}: these files do not form a valid {repo_type} repo",
                                 problems=err.problems) from None
    if _static.is_static():
        return _pr.submit(repo, files, message=message, meta=_static.meta())
    head = create_repo(str(repo), repo_type)["head"]
    operations = []
    for path, data in sorted(files.items()):
        sha = hashlib.sha256(data).hexdigest()
        _http.request("PUT", f"/api/repos/{repo.namespace}/{repo.name}/blobs/{sha}", data=data)
        operations.append({"op": "add", "path": path, "sha256": sha})
    return _http.request("POST", f"/api/repos/{repo.namespace}/{repo.name}/commit/main", json_body={
        "message": message or "Upload with chemart", "parent_commit": head,
        "operations": operations, "replace": True,
    }).json()


def _hub_yaml(block: dict[str, Any]) -> bytes:
    return yaml.safe_dump({"hub": block}, sort_keys=False, allow_unicode=True).encode()


def _network_readme(net: Network, repo: RepoId, title: str | None, description: str | None) -> str:
    """A starting card; the hub page shows the reactions and summary itself."""
    lines = [f"# {title or repo.name}", ""]
    if description:
        lines += [description, ""]
    if net.chemistry:
        how = f"Generated by `{net.chemistry}`"
        if net.seed is not None:
            how += f" with seed {net.seed}"
        lines += [how + (f" and parameters `{json.dumps(net.params)}`." if net.params else "."), ""]
    lines += [f"{len(net.species)} species, {len(net.reactions)} reactions ({net.status}).", "",
              "```python", "import chemart", "", f'net = chemart.load_network("{repo}")', "```", ""]
    return "\n".join(lines)


def network_files(net: Network, repo_id: str, *, readme: str | None = None, title: str | None = None,
                  description: str | None = None, license: str | None = None,
                  tags: list[str] | None = None) -> dict[str, bytes]:
    """The files of a network repo holding `net`."""
    repo, _ = _repo(repo_id)
    block: dict[str, Any] = {"repo_type": "network"}
    for key, value in (("title", title), ("description", description), ("license", license), ("tags", tags)):
        if value:
            block[key] = value
    return {
        "network.json": json.dumps(net.to_dict(), ensure_ascii=False).encode(),
        "chemart.yaml": _hub_yaml(block),
        "README.md": (readme or _network_readme(net, repo, title, description)).encode(),
    }


def push_network(net: Network, repo_id: str, *, message: str | None = None, readme: str | None = None,
                 title: str | None = None, description: str | None = None, license: str | None = None,
                 tags: list[str] | None = None) -> dict[str, Any]:
    """Share a reaction network as a network repo; also `Network.push_to_hub`."""
    files = network_files(net, repo_id, readme=readme, title=title, description=description,
                          license=license, tags=tags)
    return upload_files(repo_id, files, repo_type="network", message=message)


_SKIP_DIRS = {"__pycache__"}


def _collect(folder: Path) -> dict[str, bytes]:
    """The shareable files of a generator folder (one level of subfolders)."""
    files: dict[str, bytes] = {}
    for path in sorted(folder.rglob("*")):
        rel = path.relative_to(folder)
        if not path.is_file() or any(part.startswith(".") or part in _SKIP_DIRS for part in rel.parts):
            continue
        name = rel.as_posix()
        if name == "preview.json":
            continue                      # always regenerated
        try:
            check_path(name)
        except ValueError:
            continue                      # not a file type the hub stores
        files[name] = path.read_bytes()
    return files


def check(folder: str | Path) -> list[str]:
    """Run the hub's format rules and Chemart's contract on a generator folder.

    Returns the problems; an empty list means `push_generator` will be accepted.
    This imports and runs the folder's own code (it is yours).
    """
    from chemart import catalog, contract

    folder = Path(folder)
    files = _collect(folder)
    if "chemart.yaml" not in files:
        return [f"{folder} has no chemart.yaml (start one with `chemart new`)"]
    try:
        hub, entry = _format.parse_chemart_yaml(files["chemart.yaml"])
    except _format.FormatError as err:
        return err.problems
    if entry is None:
        return ["chemart.yaml has no entry under 'chemistries'"]
    problems = [f"chemart.yaml: {p}" for p in catalog.entry_problems(entry, hub=True)]
    if problems:
        return problems
    try:
        _remote_code.check_requirements(hub, str(folder))
        generate = _remote_code.local_generator(folder)
    except (ImportError, ValueError, SyntaxError) as err:
        return [f"cannot import generator.py: {type(err).__name__}: {err}"]
    return contract.problems(entry, generate)


def push_generator(folder: str | Path, repo_id: str | None = None, *, message: str | None = None,
                   check_contract: bool = True, seed: int = 0) -> dict[str, Any]:
    """Share the chemistry in `folder` (chemart.yaml + generator.py + helpers).

    Runs the contract checks on your machine, writes preview.json (the network
    at default parameters), and uploads. `repo_id` defaults to
    ``<you>/<entry id>``.
    """
    folder = Path(folder)
    if check_contract:
        problems = check(folder)
        if problems:
            raise HubValidationError(f"{folder} does not pass the checks", problems=problems)
    if repo_id is None:
        _, entry = _format.parse_chemart_yaml((folder / "chemart.yaml").read_bytes())
        repo_id = f"{_account_name()}/{entry.id}"
    files = generator_files(folder, repo_id, seed=seed)
    return upload_files(repo_id, files, repo_type="generator", message=message)


def generator_files(folder: str | Path, repo_id: str, *, seed: int = 0) -> dict[str, bytes]:
    """The files of a generator repo made from `folder`, with a fresh preview.json
    (the network at default parameters). Runs the folder's own code."""
    from chemart import api

    folder = Path(folder)
    files = _collect(folder)
    _, entry = _format.parse_chemart_yaml(files["chemart.yaml"])
    repo, _ = _repo(repo_id)
    generate = _remote_code.local_generator(folder)
    preview = api.run_generator(entry, generate, seed, {})
    preview.chemistry = str(repo)
    files["preview.json"] = json.dumps(preview.to_dict(), ensure_ascii=False).encode()
    return files


def _account_name() -> str:
    """Your name on the configured hub: its account, or your GitHub login for a static hub."""
    if _static.is_static():
        return _pr.github_login()
    return whoami()["name"]


def new(folder: str | Path, chemistry_id: str | None = None, name: str | None = None) -> list[Path]:
    """Write a working generator skeleton into `folder`."""
    from chemart.hub import _template

    folder = Path(folder)
    return _template.write(folder, chemistry_id or folder.resolve().name, name)
