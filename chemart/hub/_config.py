"""Where the hub is, where the cache lives, and the limits both sides enforce.

    CHEMART_HUB_URL      the hub's base URL (default: the public hub on GitHub Pages)
    CHEMART_HOME         cache and credentials (default: ~/.cache/chemart)
    CHEMART_HUB_TOKEN    an API token; overrides the one saved by `chemart login`
    CHEMART_HUB_OFFLINE  1 = never touch the network, use the cache only
"""

from __future__ import annotations

import json
import os
from pathlib import Path

#: The public hub: a static site built from the MarcoTuc/chemart-hub registry.
#: Point CHEMART_HUB_URL at http://127.0.0.1:8000 for a hub you run with
#: `chemart-hub serve`.
DEFAULT_URL = "https://marcotuc.github.io/chemart-hub"

KB, MB = 1024, 1024 * 1024
#: Per-file caps by kind, and per-repo caps. The server enforces them; the
#: client checks them before uploading so errors come early.
MAX_TEXT_FILE = 1 * MB          # .py .yaml .yml .md .txt LICENSE
MAX_PREVIEW = 5 * MB            # preview.json of a generator repo
MAX_NETWORK = 50 * MB           # network.json of a network repo, and any other .json
MAX_FILES = 64
MAX_REPO = 100 * MB


def hub_url() -> str:
    return (os.environ.get("CHEMART_HUB_URL") or DEFAULT_URL).rstrip("/")


def home() -> Path:
    return Path(os.environ.get("CHEMART_HOME") or Path.home() / ".cache" / "chemart")


def offline() -> bool:
    return os.environ.get("CHEMART_HUB_OFFLINE", "").strip().lower() in {"1", "true", "yes", "on"}


def _tokens_file() -> Path:
    return home() / "tokens.json"


def _saved_tokens() -> dict[str, str]:
    try:
        data = json.loads(_tokens_file().read_text())
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def token(url: str | None = None) -> str | None:
    """The token for `url` (default: the configured hub), or None."""
    return os.environ.get("CHEMART_HUB_TOKEN") or _saved_tokens().get(url or hub_url())


def save_token(value: str | None, url: str | None = None) -> None:
    """Store (or with None, forget) the token for one hub; the file is 0600."""
    tokens = _saved_tokens()
    if value is None:
        tokens.pop(url or hub_url(), None)
    else:
        tokens[url or hub_url()] = value
    path = _tokens_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as f:
        json.dump(tokens, f, indent=2)
    os.replace(tmp, path)
