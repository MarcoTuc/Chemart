"""The grammar of hub ids, revisions and file paths, shared by client and server.

    alice/my-chem              a repo; namespace is a user or an organisation
    alice/my-chem@3f2a9c1      the same, pinned to a commit (prefix of >= 7 hex digits)
    alice/my-chem@main         the moving head (the default)

Names follow the catalog's id style, so a generator repo's name is the id of
the chemistry it holds. Built-in catalog ids never contain "/" or "@".
"""

from __future__ import annotations

import re
from dataclasses import dataclass

NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
NAMESPACE_RE = re.compile(r"^[a-z0-9](?:[a-z0-9-]{0,37}[a-z0-9])?$")
COMMIT_RE = re.compile(r"^[0-9a-f]{64}$")
REVISION_RE = re.compile(r"^(?:main|[0-9a-f]{7,64})$")
MAX_NAME = 64

#: Names that are routes of the web site, so no user or org may take them.
RESERVED_NAMESPACES = frozenset({
    "api", "settings", "login", "logout", "signup", "new", "browse", "static",
    "docs", "admin", "search", "about", "help", "hub", "chemart",
})

#: One optional folder level; no dotfiles, no "..", portable characters only.
_SEGMENT = r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,99}"
PATH_RE = re.compile(rf"^(?:{_SEGMENT}/)?{_SEGMENT}$")
ALLOWED_SUFFIXES = (".py", ".yaml", ".yml", ".json", ".md", ".txt")
ALLOWED_NAMES = ("LICENSE",)


@dataclass(frozen=True)
class RepoId:
    namespace: str
    name: str

    def __str__(self) -> str:
        return f"{self.namespace}/{self.name}"

    @property
    def folder(self) -> str:
        """Cache folder name: '/' is not allowed in one path component."""
        return f"{self.namespace}--{self.name}"


def check_name(name: str) -> str:
    if not isinstance(name, str) or len(name) > MAX_NAME or not NAME_RE.match(name):
        raise ValueError(
            f"invalid repo name {name!r}: use lowercase letters, digits and single "
            f"hyphens (like a catalog id), at most {MAX_NAME} characters"
        )
    return name


def check_namespace(namespace: str) -> str:
    if not isinstance(namespace, str) or not NAMESPACE_RE.match(namespace):
        raise ValueError(
            f"invalid namespace {namespace!r}: use lowercase letters, digits and "
            "hyphens, at most 39 characters"
        )
    return namespace


def check_revision(revision: str) -> str:
    if not isinstance(revision, str) or not REVISION_RE.match(revision):
        raise ValueError(
            f"invalid revision {revision!r}: use 'main' or a commit id "
            "(at least its first 7 hex digits)"
        )
    return revision


def parse_repo_id(text: str, revision: str | None = None) -> tuple[RepoId, str]:
    """Split ``ns/name[@rev]`` into a RepoId and a revision (default 'main')."""
    if not isinstance(text, str):
        raise ValueError(f"a hub id must be a string, got {type(text).__name__}")
    body, sep, pinned = text.partition("@")
    if sep:
        if revision is not None and revision != pinned:
            raise ValueError(f"{text!r} pins revision {pinned!r} but revision={revision!r} was also given")
        revision = pinned
    namespace, slash, name = body.partition("/")
    if not slash or "/" in name:
        raise ValueError(f"invalid hub id {text!r}: expected 'namespace/name', e.g. 'alice/my-chem'")
    return RepoId(check_namespace(namespace), check_name(name)), check_revision(revision or "main")


def check_path(path: str) -> str:
    """A repo file path, or ValueError. Applied to uploads *and* to paths a
    server sends back, so a hostile server cannot write outside the cache."""
    if not isinstance(path, str) or not PATH_RE.match(path) or ".." in path:
        raise ValueError(
            f"invalid file path {path!r}: at most one folder level, no dotfiles, "
            "letters, digits, '_', '-' and '.' only"
        )
    base = path.rsplit("/", 1)[-1]
    if base not in ALLOWED_NAMES and not base.endswith(ALLOWED_SUFFIXES):
        raise ValueError(
            f"file type not allowed: {path!r} (allowed: {', '.join(ALLOWED_SUFFIXES)}, LICENSE)"
        )
    return path
