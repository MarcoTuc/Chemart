"""HTTP to the hub over urllib (no third-party dependency), and its errors.

`transport` is the single I/O seam: tests replace it to talk to a fake hub.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable
from urllib.parse import urlsplit

from chemart.hub import _config


class HubError(RuntimeError):
    """Something went wrong talking to the Chemart Hub."""

    def __init__(self, message: str, status: int | None = None, problems: list[str] | None = None):
        self.status = status
        self.problems = list(problems or [])
        detail = "".join(f"\n  - {p}" for p in self.problems)
        super().__init__(message + detail)


class HubConnectionError(HubError):
    """The hub could not be reached."""


class HubAuthError(HubError):
    """Not logged in, bad token, or not allowed to do this."""


class RepoNotFoundError(HubError, LookupError):
    """No such repo, revision or file."""


class HubConflictError(HubError):
    """The repo moved on the server while you were pushing."""


class HubValidationError(HubError, ValueError):
    """The files do not form a valid repo; `problems` lists why."""


@dataclass
class Response:
    status: int
    body: bytes
    headers: dict[str, str] = field(default_factory=dict)

    def json(self) -> Any:
        return json.loads(self.body)


def _urllib(method: str, url: str, headers: dict[str, str], body: bytes | None, timeout: float) -> Response:
    req = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return Response(r.status, r.read(), {k.lower(): v for k, v in r.headers.items()})
    except urllib.error.HTTPError as err:
        return Response(err.code, err.read(), {k.lower(): v for k, v in (err.headers or {}).items()})


#: (method, url, headers, body, timeout) -> Response
transport: Callable[..., Response] = _urllib

_warned_plain_http: set[str] = set()


def _is_local(url: str) -> bool:
    host = urlsplit(url).hostname or ""
    return host in {"localhost", "127.0.0.1", "::1"} or host.endswith(".localhost")


def request(
    method: str,
    path: str,
    *,
    json_body: Any = None,
    data: bytes | None = None,
    auth: bool = True,
    token: str | None = None,
    timeout: float = 120,
) -> Response:
    """Send one request to the configured hub; raise a HubError subclass on failure."""
    base = _config.hub_url()
    url = base + path
    from chemart import __version__

    headers = {"User-Agent": f"chemart/{__version__}", "Accept": "application/json"}
    body = data
    if json_body is not None:
        body = json.dumps(json_body).encode()
        headers["Content-Type"] = "application/json"
    elif data is not None:
        headers["Content-Type"] = "application/octet-stream"
    secret = token or (_config.token(base) if auth else None)
    if secret:
        if urlsplit(base).scheme != "https" and not _is_local(base) and base not in _warned_plain_http:
            _warned_plain_http.add(base)
            warnings.warn(f"sending your Chemart Hub token over plain http to {base}", stacklevel=3)
        headers["Authorization"] = f"Bearer {secret}"
    try:
        response = transport(method, url, headers, body, timeout)
    except (urllib.error.URLError, OSError) as err:
        reason = getattr(err, "reason", err)
        raise HubConnectionError(
            f"cannot reach the Chemart Hub at {base} ({reason}). Is it running "
            "(`chemart-hub serve`)? Set CHEMART_HUB_URL to use another hub, or "
            "CHEMART_HUB_OFFLINE=1 to work from the cache."
        ) from None
    if response.status < 400:
        return response
    raise _error(response, method, path)


def _error(response: Response, method: str, path: str) -> HubError:
    try:
        payload = response.json()
        message, problems = payload.get("error") or "", payload.get("problems") or []
    except (ValueError, AttributeError):
        message, problems = response.body[:200].decode("utf-8", "replace"), []
    message = message or f"{method} {path} failed"
    status = response.status
    if status == 401:
        return HubAuthError(f"{message}. Run `chemart login` (or set CHEMART_HUB_TOKEN).", status)
    if status == 403:
        return HubAuthError(message, status)
    if status == 404:
        return RepoNotFoundError(message, status)
    if status == 409:
        return HubConflictError(message, status)
    if status == 422:
        return HubValidationError(message, status, problems)
    return HubError(f"{message} (HTTP {status})", status, problems)
