"""Sharing on a static hub: a pull request on the registry it is built from.

    chemart push ./my-chem                      # opens a pull request
    net.push_to_hub("you/snapshot")             # the same, for a network

With the GitHub CLI (``gh``, signed in once with ``gh auth login``) the files
are committed on a new branch of your fork of the registry (or of the registry
itself, if you may push there) and a pull request is opened. Without it, the
files are written to a local folder laid out like the registry, and the steps
to open the pull request by hand are returned.

Once the pull request is merged, the registry's workflow rebuilds the site
and the repo is on the shelf.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from chemart.hub._http import HubAuthError, HubError
from chemart.hub._ids import RepoId, check_path

REPOS = "repos"


def _run(*args: str, cwd: str | Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    try:
        done = subprocess.run(list(args), cwd=cwd, capture_output=True, text=True)
    except FileNotFoundError:
        raise HubError(f"{args[0]} is not installed") from None
    if check and done.returncode != 0:
        raise HubError(f"{' '.join(args[:3])} failed: {(done.stderr or done.stdout).strip()}")
    return done


def gh_ready() -> bool:
    """Is the GitHub CLI installed and signed in?"""
    try:
        return subprocess.run(["gh", "auth", "status"], capture_output=True).returncode == 0
    except FileNotFoundError:
        return False


def github_login() -> str:
    """Your GitHub login, from the GitHub CLI."""
    if not gh_ready():
        raise HubAuthError("sharing on a static hub uses your GitHub account: install the GitHub CLI "
                           "(https://cli.github.com/) and run `gh auth login`, or give the repo id "
                           "explicitly as <your GitHub name>/<name>")
    return _run("gh", "api", "user", "--jq", ".login").stdout.strip()


def stage(repo: RepoId, files: dict[str, bytes], root: str | Path) -> Path:
    """Write `files` as repos/<ns>/<name>/ under `root`, replacing what was there."""
    folder = Path(root) / REPOS / repo.namespace / repo.name
    if folder.exists():
        shutil.rmtree(folder)
    for path, data in files.items():
        target = folder / check_path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    return folder


def _manual(repo: RepoId, files: dict[str, bytes], registry: str, branch: str) -> dict[str, Any]:
    root = Path.cwd() / "chemart-pull-request"
    folder = stage(repo, files, root)
    url = f"https://github.com/{registry}"
    steps = [
        f"Fork {url} on GitHub (the Fork button).",
        f"Copy {folder} into your fork as {REPOS}/{repo.namespace}/{repo.name}/ "
        f"on a new branch from {branch}, and commit it.",
        f"Open a pull request against {registry}:{branch}.",
        "Or install the GitHub CLI (https://cli.github.com/), run `gh auth login`, "
        "and push again: chemart then does all of this for you.",
    ]
    return {"folder": str(folder), "url": url, "steps": steps, "unchanged": False}


def submit(repo: RepoId, files: dict[str, bytes], *, message: str | None, meta: dict[str, Any]) -> dict[str, Any]:
    """Open a pull request that makes `files` the content of `repo` in the registry."""
    registry, branch = meta.get("registry"), meta.get("branch") or "main"
    if not registry or registry.count("/") != 1:
        raise HubError("this static hub does not say which GitHub repository it is built from "
                       "(hub.json has no 'registry')")
    if not gh_ready():
        return _manual(repo, files, registry, branch)

    login = github_login()
    owner, name = registry.split("/")
    can_push = _run("gh", "api", f"repos/{registry}", "--jq", ".permissions.push",
                    check=False).stdout.strip() == "true"
    target = registry
    if not can_push:
        _run("gh", "repo", "fork", registry, "--clone=false", "--remote=false", check=False)
        target = f"{login}/{name}"

    title = message or f"Share {repo}"
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    new_branch = f"chemart/{repo.namespace}-{repo.name}-{stamp}"
    helper = ["-c", "credential.helper=", "-c", "credential.helper=!gh auth git-credential"]
    with tempfile.TemporaryDirectory(prefix="chemart-pr-") as tmp:
        _run("git", *helper, "clone", "--quiet", "--depth", "1", "--branch", branch,
             f"https://github.com/{registry}.git", tmp)
        _run("git", "checkout", "--quiet", "-b", new_branch, cwd=tmp)
        stage(repo, files, tmp)
        _run("git", "add", "-A", f"{REPOS}/{repo.namespace}/{repo.name}", cwd=tmp)
        if not _run("git", "status", "--porcelain", cwd=tmp).stdout.strip():
            return {"unchanged": True, "url": f"https://github.com/{registry}/tree/{branch}/"
                                              f"{REPOS}/{repo.namespace}/{repo.name}"}
        _run("git", "-c", f"user.name={login}", "-c", f"user.email={login}@users.noreply.github.com",
             "commit", "--quiet", "-m", title, cwd=tmp)
        push_url = f"https://github.com/{target}.git"
        for attempt in range(6):                      # a fresh fork takes a moment to exist
            done = _run("git", *helper, "push", "--quiet", push_url, f"HEAD:refs/heads/{new_branch}",
                        cwd=tmp, check=False)
            if done.returncode == 0:
                break
            time.sleep(5)
        else:
            raise HubError(f"could not push to {target}: {done.stderr.strip()}")

    head = new_branch if target == registry else f"{login}:{new_branch}"
    body = (f"Shares `{repo}` on the Chemart Hub.\n\n"
            f"Files: {', '.join(f'`{p}`' for p in sorted(files))}\n\n"
            "Opened by `chemart`. The registry's checks validate the repo; a maintainer merges it.")
    pr = _run("gh", "pr", "create", "-R", registry, "--base", branch, "--head", head,
              "--title", title, "--body", body).stdout.strip().splitlines()[-1]
    return {"pull_request": pr, "url": pr, "branch": new_branch, "unchanged": False}


def describe(result: dict[str, Any]) -> str:
    """One line (or a few) for the CLI."""
    if result.get("pull_request"):
        return f"opened a pull request: {result['pull_request']}"
    if result.get("folder"):
        return "\n".join([f"wrote {result['folder']}; to share it:"] + [f"  {i}. {s}" for i, s in
                                                                         enumerate(result["steps"], 1)])
    if result.get("unchanged"):
        return f"no changes; the registry already has these files ({result['url']})"
    return json.dumps(result)
