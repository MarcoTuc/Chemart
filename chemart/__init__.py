"""Chemart: a mart of artificial chemistries.

    import chemart
    chemart.list_chemistries()
    chemart.describe_chemistry("matrix-chemistry")
    net = chemart.generate_network("matrix-chemistry", seed=0, N=4)
"""

from importlib import import_module

_EXPORTS = {
    "call_tool": "chemart.api",
    "describe_chemistry": "chemart.api",
    "generate_network": "chemart.api",
    "list_chemistries": "chemart.api",
    "tool_definitions": "chemart.api",
    "Network": "chemart.network",
    "Reaction": "chemart.network",
    "Species": "chemart.network",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    # Lazy, so `python -m chemart.catalog` does not import the catalog twice.
    if name in _EXPORTS:
        return getattr(import_module(_EXPORTS[name]), name)
    raise AttributeError(f"module 'chemart' has no attribute {name!r}")
