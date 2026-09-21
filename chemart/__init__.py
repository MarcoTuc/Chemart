"""Chemart: a mart of artificial chemistries.

    import chemart
    chemart.list_chemistries()
    chemart.describe_chemistry("matrix-chemistry")
    net = chemart.generate_network("matrix-chemistry", seed=0, N=4)

Chemistries shared on the Chemart Hub load the same way, by ``namespace/name``:

    net = chemart.generate_network("alice/my-chem", trust_remote_code=True, revision="3f2a9c1")
    net = chemart.load_network("alice/some-network")
    net.push_to_hub("alice/my-network")
"""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("chemart")
except PackageNotFoundError:  # running from a source tree that was never installed
    __version__ = "0.0.0+local"

_EXPORTS = {
    "call_tool": "chemart.api",
    "describe_chemistry": "chemart.api",
    "generate_network": "chemart.api",
    "list_chemistries": "chemart.api",
    "tool_definitions": "chemart.api",
    "load_network": "chemart.hub",
    "Network": "chemart.network",
    "Reaction": "chemart.network",
    "Species": "chemart.network",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    # Lazy, so `python -m chemart.catalog` does not import the catalog twice.
    if name in _EXPORTS:
        return getattr(import_module(_EXPORTS[name]), name)
    if name == "hub":
        return import_module("chemart.hub")
    raise AttributeError(f"module 'chemart' has no attribute {name!r}")
