"""Server settings, from arguments or the environment.

    CHEMART_HUB_DATA         data folder: database + blob store (default ./hub-data)
    CHEMART_HUB_PUBLIC_URL   the URL people reach the hub at (for links and snippets)
    CHEMART_HUB_SECURE       1 = cookies only over https (set this behind TLS)
    CHEMART_HUB_SIGNUP       0 = no self-service sign-up
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _flag(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    return default if raw is None else raw.strip().lower() in {"1", "true", "yes", "on"}


@dataclass
class Settings:
    data_dir: Path
    public_url: str = "http://127.0.0.1:8000"
    #: The organisation whose repos may be backed by the built-in catalog.
    official_namespace: str = "chemart"
    secure_cookies: bool = False
    allow_signup: bool = True
    session_days: int = 30

    def __post_init__(self) -> None:
        self.data_dir = Path(self.data_dir)
        self.public_url = self.public_url.rstrip("/")

    @property
    def db_path(self) -> Path:
        return self.data_dir / "hub.sqlite3"

    @property
    def blobs_dir(self) -> Path:
        return self.data_dir / "blobs"

    @classmethod
    def from_env(cls, **overrides) -> Settings:
        values = {
            "data_dir": Path(os.environ.get("CHEMART_HUB_DATA", "hub-data")),
            "public_url": os.environ.get("CHEMART_HUB_PUBLIC_URL", "http://127.0.0.1:8000"),
            "secure_cookies": _flag("CHEMART_HUB_SECURE", False),
            "allow_signup": _flag("CHEMART_HUB_SIGNUP", True),
        }
        values.update({k: v for k, v in overrides.items() if v is not None})
        return cls(**values)
