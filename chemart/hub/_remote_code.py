"""Chemistries whose code comes from the hub, and the gate in front of it.

A generator repo's code runs on *your* machine with *your* permissions, so
`generate_network` refuses to import it unless you pass
``trust_remote_code=True`` -- the same contract as `transformers`. Pin
``revision=`` to the commit you read, or the next push to 'main' runs too.

Loading: the snapshot folder becomes a synthetic package
``_chemart_hub_<commit>``, registered in sys.modules *before* its
``generator`` module is imported. That makes relative imports between the
repo's files work, and keeps modules that use @dataclass happy (dataclasses
look their module up in sys.modules while the class is being built).
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.machinery
import importlib.util
import sys
import warnings
from dataclasses import dataclass, field, fields
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

from chemart import catalog
from chemart.hub import _cache, _config, _format
from chemart.hub._ids import RepoId


@dataclass
class RemoteChemistry(catalog.Chemistry):
    """A catalog entry that lives in a hub repo at one commit."""

    repo: str = ""
    commit: str = ""
    revision: str = "main"
    hub: dict[str, Any] = field(default_factory=dict)
    manifest: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_entry(cls, entry: catalog.Chemistry, **extra: Any) -> RemoteChemistry:
        return cls(**{f.name: getattr(entry, f.name) for f in fields(catalog.Chemistry)}, **extra)

    @property
    def module(self) -> str:
        return f"{_package_name(self.commit)}.generator"

    @property
    def implemented(self) -> bool:
        return True

    @property
    def ref(self) -> str:
        """Provenance recorded in generated networks: enough to reproduce them."""
        return f"{self.repo}@{self.commit}"

    @property
    def repo_id(self) -> RepoId:
        namespace, name = self.repo.split("/")
        return RepoId(namespace, name)

    def hub_info(self) -> dict[str, Any]:
        return {
            "repo": self.repo,
            "commit": self.commit,
            "revision": self.revision,
            "url": f"{_config.hub_url()}/{self.repo}",
            "has_code": True,
            "license": self.hub.get("license"),
            "tags": self.hub.get("tags", []),
            "requires": self.hub.get("requires", []),
            "requires_chemart": self.hub.get("requires_chemart"),
            "files": sorted(self.manifest.get("files", {})),
        }

    def load_generator(self, trust_remote_code: bool = False) -> Callable:
        return load(self, trust_remote_code)


def _package_name(key: str) -> str:
    return f"_chemart_hub_{key[:24]}"


def import_folder(folder: Path, package: str, *, fresh: bool = False):
    """Import `folder/generator.py` as `<package>.generator` and return the module."""
    if fresh:
        for name in [m for m in sys.modules if m == package or m.startswith(package + ".")]:
            del sys.modules[name]
    if package not in sys.modules:
        spec = importlib.machinery.ModuleSpec(package, None, is_package=True)
        pkg = importlib.util.module_from_spec(spec)
        pkg.__path__ = [str(folder)]
        sys.modules[package] = pkg
    return importlib.import_module(f"{package}.generator")


def local_generator(folder: Path) -> Callable:
    """The `generate` of a chemistry being developed in `folder` (your own code)."""
    folder = Path(folder).resolve()
    package = "_chemart_local_" + hashlib.sha256(str(folder).encode()).hexdigest()[:16]
    module = import_folder(folder, package, fresh=True)
    if not callable(getattr(module, "generate", None)):
        raise ValueError(f"{folder / 'generator.py'} must define generate(p, rng)")
    return module.generate


def check_requirements(hub: dict[str, Any], what: str) -> None:
    spec = hub.get("requires_chemart")
    if spec:
        try:
            have = version("chemart")
        except PackageNotFoundError:
            have = None
        if have is not None and not _format.version_satisfies(have, spec):
            raise ImportError(f"{what} needs chemart {spec}, but chemart {have} is installed")
    missing = [m for m in hub.get("requires", []) if importlib.util.find_spec(m.split(".")[0]) is None]
    if missing:
        raise ImportError(
            f"{what} needs {', '.join(missing)}, which is not installed. Chemart never "
            f"installs anything for you: add it yourself (uv add {' '.join(missing)})."
        )


def load(c: RemoteChemistry, trust_remote_code: bool) -> Callable:
    url = f"{_config.hub_url()}/{c.repo}/blob/{c.commit}/generator.py"
    if not trust_remote_code:
        raise ValueError(
            f"{c.repo} runs its own Python code on your machine (generator.py at commit "
            f"{c.commit[:12]}). Read it first: {url}\n"
            f"If you trust it, pass trust_remote_code=True, and pin the version you read with "
            f"revision={c.commit[:12]!r} so later pushes cannot change what runs."
        )
    check_requirements(c.hub, c.repo)
    repo = c.repo_id
    info = c.manifest
    folder = _cache.fetch(repo, info, [p for p in info["files"] if p.endswith(".py")])
    # A commit id (or a unique prefix of one) pins the code; 'main' does not.
    first_run = not _cache.trusted_before(repo, c.commit)
    if first_run and c.revision == "main":
        warnings.warn(
            f"running code from {c.repo}@{c.revision}, which resolved to commit {c.commit[:12]}. "
            f"Make sure it holds nothing malicious ({url}); to keep running exactly this code, "
            f"pass revision={c.commit[:12]!r}.",
            stacklevel=5,  # load <- load_generator <- generator_for <- generate_network <- caller
        )
    module = import_folder(folder, _package_name(c.commit))
    generate = getattr(module, "generate", None)
    if not callable(generate):
        raise ImportError(f"{c.repo}@{c.commit[:12]}: generator.py defines no generate(p, rng)")
    return generate
