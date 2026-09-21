"""Turning repo content into safe HTML: markdown cards, code, previews.

README files are user content: markdown is rendered, then *everything* goes
through nh3's allow-list, so no script, handler attribute or javascript: URL
survives. Code is highlighted by Pygments and never interpreted.
"""

from __future__ import annotations

import json
import threading
from collections import OrderedDict
from datetime import UTC, datetime
from typing import Any

import markdown as _markdown
import nh3
from markupsafe import Markup
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import TextLexer, get_lexer_for_filename
from pygments.util import ClassNotFound

from chemart.api import params_schema
from chemart.hub import _format

_TAGS = {
    "a", "abbr", "b", "blockquote", "br", "code", "dd", "del", "details", "div", "dl", "dt",
    "em", "h1", "h2", "h3", "h4", "h5", "h6", "hr", "i", "img", "kbd", "li", "ol", "p", "pre",
    "s", "span", "strong", "sub", "summary", "sup", "table", "tbody", "td", "th", "thead", "tr", "ul",
}
_ATTRIBUTES = {
    "a": {"href", "title"},
    "img": {"src", "alt", "title", "width", "height"},
    "code": {"class"}, "span": {"class"}, "div": {"class"}, "pre": {"class"},
    "th": {"align"}, "td": {"align"}, "abbr": {"title"},
}
_FORMATTER = HtmlFormatter(cssclass="hl", wrapcode=True)


def sanitize(html: str) -> Markup:
    return Markup(nh3.clean(html, tags=_TAGS, attributes=_ATTRIBUTES,
                            url_schemes={"http", "https", "mailto"}, link_rel="nofollow ugc noopener"))


def markdown_html(text: str) -> Markup:
    html = _markdown.markdown(
        text,
        extensions=["fenced_code", "tables", "sane_lists", "codehilite"],
        extension_configs={"codehilite": {"css_class": "hl", "guess_lang": False}},
        output_format="html",
    )
    return sanitize(html)


def code_html(path: str, text: str) -> Markup:
    try:
        lexer = get_lexer_for_filename(path, stripnl=False)
    except ClassNotFound:
        lexer = TextLexer(stripnl=False)
    return Markup(highlight(text, lexer, _FORMATTER))


def pygments_css() -> str:
    light = HtmlFormatter(style="friendly").get_style_defs(".hl")
    dark = HtmlFormatter(style="monokai").get_style_defs(".hl")
    return (f"{light}\n.hl, .hl pre {{ background: transparent; }}\n"
            f"@media (prefers-color-scheme: dark) {{\n{dark}\n.hl, .hl pre {{ background: transparent; }}\n}}\n")


# --------------------------------------------------------------------------
# Repo views, cached per commit (commits never change)
# --------------------------------------------------------------------------

class ViewCache:
    def __init__(self, size: int = 256):
        self.size = size
        self._data: OrderedDict[str, dict[str, Any]] = OrderedDict()
        self._lock = threading.Lock()

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    def get(self, key: str, build) -> dict[str, Any]:
        with self._lock:
            if key in self._data:
                self._data.move_to_end(key)
                return self._data[key]
        value = build()
        with self._lock:
            self._data[key] = value
            while len(self._data) > self.size:
                self._data.popitem(last=False)
        return value


def _bounds(p) -> str:
    parts = []
    if p.choices:
        parts.append(" | ".join(map(str, p.choices)))
    if p.min is not None or p.max is not None:
        lo = "−∞" if p.min is None else f"{p.min:g}"
        hi = "∞" if p.max is None else f"{p.max:g}"
        parts.append(f"[{lo}, {hi}]")
    if p.range:
        parts.append(p.range)
    return "; ".join(parts)


def repo_view(card: _format.RepoCard, max_reactions: int = 40) -> dict[str, Any]:
    """Everything the repo page shows, derived from one snapshot's card."""
    entry = card.entry
    view: dict[str, Any] = {
        "readme_html": markdown_html(card.readme) if card.readme else None,
        "entry": entry,
        "imports": card.imports,
        "params": [],
        "preview": None,
        "hub": card.hub,
    }
    if entry is not None:
        schema = params_schema(entry)["properties"]
        view["params"] = [
            {"name": p.name, "type": p.type, "default": json.dumps(p.default, ensure_ascii=False),
             "bounds": _bounds(p), "role": p.role, "meaning": " ".join((p.meaning or "").split()),
             "schema": schema.get(p.name, {})}
            for p in entry.params
        ]
    net = card.network
    if net is not None:
        lines = net.to_text().splitlines()
        view["preview"] = {
            "summary": net.summary(),
            "species": len(net.species),
            "reactions_n": len(net.reactions),
            "reactions": lines[:max_reactions],
            "more": max(0, len(lines) - max_reactions),
            "provides": net.provides,
            "status": net.status,
            "chemistry": net.chemistry,
            "params": json.dumps(net.params, ensure_ascii=False) if net.params else "",
            "seed": net.seed,
        }
    return view


def ago(stamp: str | None) -> str:
    """'2024-05-01T10:00:00Z' -> '3 days ago'."""
    if not stamp:
        return ""
    try:
        then = datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=UTC)
    except ValueError:
        return stamp
    seconds = int((datetime.now(UTC) - then).total_seconds())
    for unit, size in (("year", 31536000), ("month", 2592000), ("day", 86400), ("hour", 3600), ("minute", 60)):
        if seconds >= size:
            n = seconds // size
            return f"{n} {unit}{'s' if n > 1 else ''} ago"
    return "just now"


def human_size(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB"):
        if size < 1024 or unit == "MB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{n} B"
