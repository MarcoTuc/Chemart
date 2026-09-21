"""The web site: server-rendered pages, no JavaScript build.

Pages read the same service layer as the API. Forms carry a CSRF token: the
session's own when logged in, otherwise a double-submit cookie. Fixed routes
(/browse, /login, ...) are registered before /{namespace}/{name}; those names
are reserved, so no user can shadow them.
"""

from __future__ import annotations

import hmac
import secrets
import sqlite3
import time
from collections import defaultdict, deque
from pathlib import Path
from typing import Any
from urllib.parse import urlencode

from fastapi import APIRouter, Depends, FastAPI, Form, Request
from fastapi.responses import HTMLResponse, PlainTextResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from chemart import catalog
from chemart.hub._ids import check_path

from chemart_hub import auth, curation, db, deps, render, service
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError, forbidden, not_found
from chemart_hub.storage import BlobStore

HERE = Path(__file__).resolve().parent.parent
templates = Jinja2Templates(directory=str(HERE / "templates"))
templates.env.filters["ago"] = render.ago
templates.env.filters["size"] = render.human_size
templates.env.globals["FAMILIES"] = sorted(catalog.FAMILIES)
templates.env.globals["PROVIDES"] = sorted(catalog.PROVIDES)
#: Static mode: the site is pre-rendered from a GitHub registry (chemart_hub.static),
#: so logins, likes and edit forms give way to "contribute by pull request".
templates.env.globals["STATIC"] = False
templates.env.globals["REGISTRY_URL"] = ""
templates.env.globals["REGISTRY_BRANCH"] = "main"

ANON_CSRF = "chemart_csrf"
PAGE = 24
router = APIRouter()

#: The storefront banner from the README, shown on the home page.
BANNER = """\
     °       H     H       O
    ∘  °      ╲   ╱       ╱ ╲
     °         C═C       H   H      ▗▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▖
    ° ∘       ╱   ╲                 ▐▓▓▓▓░░░░▓▓▓▓░░░░▓▓▓▓░░░░▓▓▓▓░░░░▓▓▌
    ┌─┐      H     H                ╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱
    │ │      ░░░░░░╗░░╗  ░░╗░░░░░░░╗███╗   ███╗ █████╗ ██████╗ ████████╗
    │ │     ░░╔════╝░░║  ░░║░░╔════╝████╗ ████║██╔══██╗██╔══██╗╚══██╔══╝
   ╱ ° ╲    ░░║     ░░░░░░░║░░░░░╗  ██╔████╔██║███████║██████╔╝   ██║
  ╱ ∘ ° ╲   ░░║     ░░╔══░░║░░╔══╝  ██║╚██╔╝██║██╔══██║██╔══██╗   ██║
 ╱≈≈≈≈≈≈≈╲  ╚░░░░░░╗░░║  ░░║░░░░░░░╗██║ ╚═╝ ██║██║  ██║██║  ██║   ██║
 ╰───────╯   ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝"""


# --------------------------------------------------------------------------
# Plumbing
# --------------------------------------------------------------------------

def _csrf(request: Request, p: auth.Principal | None) -> str:
    if p is not None and p.via == "session":
        return p.csrf
    return getattr(request.state, "anon_csrf", "")


def _check_csrf(request: Request, p: auth.Principal | None, sent: str) -> None:
    expected = p.csrf if (p is not None and p.via == "session") else request.cookies.get(ANON_CSRF)
    if not expected or not hmac.compare_digest(sent or "", expected):
        raise forbidden("the form expired or came from another site; reload the page and try again")


def _site(settings: Settings) -> dict[str, Any]:
    """Editable site texts; markdown ones also rendered (sanitised) as `<key>_html`."""
    c = db.connect(settings.db_path)
    try:
        texts = curation.site_texts(c)
    finally:
        c.close()
    for key, (_, _, kind, _) in curation.SITE_TEXT.items():
        if kind == "markdown":
            texts[f"{key}_html"] = render.markdown_html(texts[key]) if texts[key] else None
    return texts


def _page(request: Request, name: str, p: auth.Principal | None, status: int = 200, **ctx: Any) -> HTMLResponse:
    settings: Settings = request.app.state.settings
    ctx.update(me=p, csrf=_csrf(request, p), settings=settings, path=request.url.path,
               signup=settings.allow_signup, site=_site(settings))
    return templates.TemplateResponse(request, name, ctx, status_code=status)


def _safe_next(target: str | None) -> str:
    if target and target.startswith("/") and not target.startswith("//") and "\\" not in target:
        return target
    return "/"


def _login_redirect(request: Request) -> RedirectResponse:
    return RedirectResponse(f"/login?{urlencode({'next': request.url.path})}", status_code=303)


class _Throttle:
    """At most `limit` failures per key in `window` seconds (per process)."""

    def __init__(self, limit: int = 10, window: float = 600):
        self.limit, self.window = limit, window
        self.failures: dict[str, deque] = defaultdict(deque)

    def blocked(self, *keys: str) -> bool:
        now = time.monotonic()
        for key in keys:
            q = self.failures[key]
            while q and now - q[0] > self.window:
                q.popleft()
            if len(q) >= self.limit:
                return True
        return False

    def fail(self, *keys: str) -> None:
        for key in keys:
            self.failures[key].append(time.monotonic())


_throttle = _Throttle()


def _cards(c: sqlite3.Connection, rows, settings: Settings) -> list[dict[str, Any]]:
    return [service.repo_json(c, r, settings) for r in rows]


# --------------------------------------------------------------------------
# Home and browse
# --------------------------------------------------------------------------

@router.get("/", response_class=HTMLResponse)
def home(request: Request, c: sqlite3.Connection = Depends(deps.conn),
         p: auth.Principal | None = Depends(deps.principal),
         settings: Settings = Depends(deps.settings)) -> HTMLResponse:
    # Trending only means something once there are likes or downloads.
    trending = [r for r in service.search(c, sort="trending", limit=6)[0] if r["likes"] or r["downloads"]]
    featured, _ = service.search(c, featured=True, sort="updated", limit=12)
    recent, _ = service.search(c, sort="updated", limit=6)
    official, n_official = service.search(c, author=settings.official_namespace, sort="name", limit=12)
    counts = dict(c.execute(
        "SELECT repo_type, COUNT(*) FROM repos WHERE head IS NOT NULL AND archived_at IS NULL "
        "AND hidden = 0 GROUP BY repo_type").fetchall())
    return _page(request, "home.html", p, banner=BANNER, featured=_cards(c, featured, settings),
                 trending=_cards(c, trending, settings), recent=_cards(c, recent, settings),
                 official=_cards(c, official, settings), n_official=n_official,
                 n_generators=counts.get("generator", 0), n_networks=counts.get("network", 0))


@router.get("/browse", response_class=HTMLResponse)
def browse(request: Request, q: str = "", type: str = "", family: str = "", tag: str = "",
           fidelity: str = "", code: str = "", author: str = "", sort: str = "trending", page: int = 1,
           c: sqlite3.Connection = Depends(deps.conn),
           p: auth.Principal | None = Depends(deps.principal),
           settings: Settings = Depends(deps.settings)) -> HTMLResponse:
    provides = [v for v in request.query_params.getlist("provides") if v]
    page = max(1, page)
    rows, total = service.search(
        c, q=q or None, repo_type=type or None, family=family or None, tag=tag or None,
        fidelity=fidelity or None, author=author or None, provides=provides,
        has_code={"yes": True, "no": False}.get(code), sort=sort if sort in service.SORTS else "trending",
        limit=PAGE, offset=(page - 1) * PAGE, include_hidden=bool(p and p.is_admin),
    )
    filters = {"q": q, "type": type, "family": family, "tag": tag, "fidelity": fidelity,
               "code": code, "author": author, "sort": sort, "provides": provides}

    def link(**change: Any) -> str:
        merged = {**filters, **change}
        pairs = [(k, v) for k, v in merged.items() if k != "provides" and v not in ("", None, 1)]
        pairs += [("provides", v) for v in merged["provides"]]
        return "/browse" + (f"?{urlencode(pairs)}" if pairs else "")

    return _page(request, "browse.html", p, repos=_cards(c, rows, settings), total=total,
                 filters=filters, facets=service.facet_counts(c), page=page,
                 pages=max(1, (total + PAGE - 1) // PAGE), link=link, sorts=service.SORTS)


# --------------------------------------------------------------------------
# Accounts
# --------------------------------------------------------------------------

@router.get("/login", response_class=HTMLResponse)
def login_form(request: Request, next: str = "/",
               p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is not None and p.via == "session":
        return RedirectResponse(_safe_next(next), status_code=303)
    return _page(request, "login.html", p, next=_safe_next(next), error=None)


@router.post("/login", response_class=HTMLResponse)
def login(request: Request, username: str = Form(...), password: str = Form(...),
          csrf: str = Form(""), next: str = Form("/"),
          c: sqlite3.Connection = Depends(deps.conn),
          p: auth.Principal | None = Depends(deps.principal),
          settings: Settings = Depends(deps.settings)) -> Response:
    _check_csrf(request, None, csrf)
    name = username.strip().lower()
    ip = request.client.host if request.client else "?"
    if _throttle.blocked(f"ip:{ip}", f"user:{name}"):
        return _page(request, "login.html", p, 429, next=_safe_next(next),
                     error="Too many failed attempts. Wait a few minutes and try again.")
    user = service.account(c, name)
    if user is None or user["is_org"] or not auth.verify_password(password, user["password_hash"]):
        _throttle.fail(f"ip:{ip}", f"user:{name}")
        return _page(request, "login.html", p, 400, next=_safe_next(next),
                     error="Wrong username or password.", username=name)
    if user["disabled"]:
        return _page(request, "login.html", p, 403, next=_safe_next(next),
                     error="This account is suspended. Contact the hub's administrators.", username=name)
    response = RedirectResponse(_safe_next(next), status_code=303)
    _start_session(response, c, user["id"], settings)
    return response


def _start_session(response: Response, c: sqlite3.Connection, user_id: int, settings: Settings) -> None:
    secret = auth.create_session(c, user_id, settings.session_days)
    response.set_cookie(auth.SESSION_COOKIE, secret, max_age=settings.session_days * 86400,
                        httponly=True, samesite="lax", secure=settings.secure_cookies)


@router.get("/signup", response_class=HTMLResponse)
def signup_form(request: Request, p: auth.Principal | None = Depends(deps.principal),
                settings: Settings = Depends(deps.settings)) -> Response:
    if not settings.allow_signup:
        raise not_found("sign-up page")
    return _page(request, "signup.html", p, error=None)


@router.post("/signup", response_class=HTMLResponse)
def signup(request: Request, username: str = Form(...), password: str = Form(...),
           email: str = Form(""), csrf: str = Form(""),
           c: sqlite3.Connection = Depends(deps.conn),
           p: auth.Principal | None = Depends(deps.principal),
           settings: Settings = Depends(deps.settings)) -> Response:
    if not settings.allow_signup:
        raise not_found("sign-up page")
    _check_csrf(request, None, csrf)
    try:
        user = service.create_user(c, username, password, email.strip() or None)
    except ApiError as err:
        return _page(request, "signup.html", p, 400, error=err.message, username=username, email=email)
    response = RedirectResponse("/settings/tokens?welcome=1", status_code=303)
    _start_session(response, c, user["id"], settings)
    return response


@router.post("/logout")
def logout(request: Request, csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
           p: auth.Principal | None = Depends(deps.principal)) -> Response:
    _check_csrf(request, p, csrf)
    auth.end_session(c, request.cookies.get(auth.SESSION_COOKIE))
    response = RedirectResponse("/", status_code=303)
    response.delete_cookie(auth.SESSION_COOKIE)
    return response


@router.get("/settings/tokens", response_class=HTMLResponse)
def tokens_page(request: Request, welcome: int = 0, c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    return _page(request, "tokens.html", p, tokens=_tokens(c, p), new_token=None, welcome=bool(welcome))


def _tokens(c: sqlite3.Connection, p: auth.Principal) -> list[sqlite3.Row]:
    return c.execute("SELECT * FROM tokens WHERE user_id = ? ORDER BY revoked, created_at DESC", (p.id,)).fetchall()


@router.post("/settings/tokens", response_class=HTMLResponse)
def create_token(request: Request, name: str = Form("token"), scope: str = Form("write"),
                 csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                 p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    secret = auth.create_token(c, p.id, name.strip() or "token", scope)
    return _page(request, "tokens.html", p, tokens=_tokens(c, p), new_token=secret, welcome=False)


@router.post("/settings/tokens/{token_id}/revoke")
def revoke_token(request: Request, token_id: int, csrf: str = Form(""),
                 c: sqlite3.Connection = Depends(deps.conn),
                 p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    with c:
        c.execute("UPDATE tokens SET revoked = 1 WHERE id = ? AND user_id = ?", (token_id, p.id))
    return RedirectResponse("/settings/tokens", status_code=303)


@router.get("/new", response_class=HTMLResponse)
def new_page(request: Request, c: sqlite3.Connection = Depends(deps.conn),
             p: auth.Principal | None = Depends(deps.principal)) -> HTMLResponse:
    namespaces = [p.name, *auth.namespaces_of(c, p.id)] if p else []
    return _page(request, "new.html", p, namespaces=namespaces, error=None)


@router.post("/new", response_class=HTMLResponse)
def new_repo(request: Request, namespace: str = Form(...), name: str = Form(...),
             repo_type: str = Form(...), csrf: str = Form(""),
             c: sqlite3.Connection = Depends(deps.conn),
             p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    try:
        repo = service.create_repo(c, p, name.strip(), repo_type, namespace)
    except ApiError as err:
        namespaces = [p.name, *auth.namespaces_of(c, p.id)]
        return _page(request, "new.html", p, 400, namespaces=namespaces, error=err.message)
    return RedirectResponse(f"/{repo['namespace']}/{repo['name']}", status_code=303)


@router.get("/about", response_class=HTMLResponse)
def about(request: Request, p: auth.Principal | None = Depends(deps.principal)) -> HTMLResponse:
    return _page(request, "about.html", p)


@router.get("/settings/account", response_class=HTMLResponse)
def account_page(request: Request, p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    return _page(request, "account.html", p, error=None, done=request.query_params.get("done"))


@router.post("/settings/account/password", response_class=HTMLResponse)
def change_password(request: Request, current: str = Form(""), password: str = Form(...),
                    confirm: str = Form(...), csrf: str = Form(""),
                    c: sqlite3.Connection = Depends(deps.conn),
                    p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    try:
        if password != confirm:
            raise ApiError(422, "the two new passwords differ")
        if not auth.verify_password(current, p.user["password_hash"]):
            raise ApiError(403, "the current password is wrong")
        curation.set_password(c, p, p.name, password, current=current)
    except ApiError as err:
        return _page(request, "account.html", p, 400, error=err.message, done=None)
    return RedirectResponse("/settings/account?done=password", status_code=303)


def _profile_context(c: sqlite3.Connection, owner: sqlite3.Row) -> dict[str, Any]:
    members = c.execute(
        "SELECT a.name FROM org_members m JOIN accounts a ON a.id = m.user_id "
        "WHERE m.org_id = ? AND a.archived_at IS NULL ORDER BY a.name",
        (owner["id"],)).fetchall() if owner["is_org"] else []
    return {"owner": owner, "members": [m["name"] for m in members]}


@router.get("/settings/profile/{name}", response_class=HTMLResponse)
def profile_edit_page(request: Request, name: str, c: sqlite3.Connection = Depends(deps.conn),
                      p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    owner = service.account(c, name)
    if owner is None:
        raise not_found(f"account {name!r}")
    if not auth.can_write(c, p, owner):
        raise forbidden("you cannot edit this profile")
    return _page(request, "profile_edit.html", p, error=None, **_profile_context(c, owner))


@router.post("/settings/profile/{name}", response_class=HTMLResponse)
def profile_edit(request: Request, name: str, fullname: str = Form(""), bio: str = Form(""),
                 email: str = Form(""), csrf: str = Form(""),
                 c: sqlite3.Connection = Depends(deps.conn),
                 p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    curation.update_profile(c, p, name, fullname=fullname, bio=bio, email=email)
    return RedirectResponse(f"/{name.lower()}", status_code=303)


@router.post("/settings/profile/{name}/members", response_class=HTMLResponse)
def org_members(request: Request, name: str, user: str = Form(...), action: str = Form("add"),
                csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return _login_redirect(request)
    _check_csrf(request, p, csrf)
    curation.set_membership(c, p, name, user.strip(), action == "add")
    return RedirectResponse(f"/settings/profile/{name.lower()}", status_code=303)


# --------------------------------------------------------------------------
# Profiles and repos
# --------------------------------------------------------------------------

@router.get("/{namespace}", response_class=HTMLResponse)
def profile(request: Request, namespace: str, c: sqlite3.Connection = Depends(deps.conn),
            p: auth.Principal | None = Depends(deps.principal),
            settings: Settings = Depends(deps.settings)) -> HTMLResponse:
    owner = service.account(c, namespace)
    if owner is None:
        raise not_found(f"user or organisation {namespace!r}")
    mine = auth.can_write(c, p, owner)
    rows, total = service.search(c, author=owner["name"], sort="updated", limit=100,
                                 include_empty=mine, include_hidden=mine)
    liked = c.execute(
        "SELECT r.*, a.name AS namespace FROM likes l JOIN repos r ON r.id = l.repo_id "
        "JOIN accounts a ON a.id = r.owner_id WHERE l.user_id = ? AND r.head IS NOT NULL AND r.hidden = 0 "
        "AND r.archived_at IS NULL "
        "ORDER BY l.created_at DESC LIMIT 24", (owner["id"],)).fetchall()
    members = c.execute(
        "SELECT a.name FROM org_members m JOIN accounts a ON a.id = m.user_id "
        "WHERE m.org_id = ? AND a.archived_at IS NULL ORDER BY a.name",
        (owner["id"],)).fetchall() if owner["is_org"] else []
    return _page(request, "profile.html", p, owner=owner, repos=_cards(c, rows, settings), total=total,
                 liked=_cards(c, liked, settings), members=[m["name"] for m in members],
                 orgs=[] if owner["is_org"] else auth.namespaces_of(c, owner["id"]), can_edit=mine,
                 deleted=request.query_params.get("deleted") if mine else None,
                 bio_html=render.markdown_html(owner["bio"]) if owner["bio"] else None)


def _repo_context(request: Request, c: sqlite3.Connection, namespace: str, name: str,
                  p: auth.Principal | None) -> dict[str, Any]:
    settings: Settings = request.app.state.settings
    repo = service.get_repo(c, namespace, name)
    info = service.repo_json(c, repo, settings)
    owner = c.execute("SELECT * FROM accounts WHERE id = ?", (repo["owner_id"],)).fetchone()
    return {"repo": info, "row": repo, "liked": service.has_liked(c, p, repo),
            "can_write": auth.can_write(c, p, owner), "is_admin": bool(p and p.is_admin)}


@router.get("/{namespace}/{name}", response_class=HTMLResponse)
def repo_page(request: Request, namespace: str, name: str,
              c: sqlite3.Connection = Depends(deps.conn),
              p: auth.Principal | None = Depends(deps.principal),
              store: BlobStore = Depends(deps.store),
              settings: Settings = Depends(deps.settings)) -> HTMLResponse:
    ctx = _repo_context(request, c, namespace, name, p)
    repo = ctx["row"]
    view = None
    if repo["head"]:
        commit = service.resolve_revision(c, repo, repo["head"])
        cache: render.ViewCache = request.app.state.views
        view = cache.get(commit["id"], lambda: render.repo_view(service.card_for(c, store, repo, commit, settings)))
    return _page(request, "repo.html", p, tab="card", view=view, **ctx)


@router.get("/{namespace}/{name}/tree/{revision}", response_class=HTMLResponse)
def tree(request: Request, namespace: str, name: str, revision: str,
         c: sqlite3.Connection = Depends(deps.conn),
         p: auth.Principal | None = Depends(deps.principal)) -> HTMLResponse:
    ctx = _repo_context(request, c, namespace, name, p)
    commit = service.resolve_revision(c, ctx["row"], revision)
    files = sorted(service.manifest(commit).items())
    return _page(request, "tree.html", p, tab="files", commit=service.commit_json(c, commit),
                 files=files, revision=revision, **ctx)


@router.get("/{namespace}/{name}/blob/{revision}/{path:path}", response_class=HTMLResponse)
def blob(request: Request, namespace: str, name: str, revision: str, path: str,
         c: sqlite3.Connection = Depends(deps.conn),
         p: auth.Principal | None = Depends(deps.principal),
         store: BlobStore = Depends(deps.store)) -> HTMLResponse:
    ctx = _repo_context(request, c, namespace, name, p)
    commit = service.resolve_revision(c, ctx["row"], revision)
    try:
        check_path(path)
    except ValueError:
        raise not_found(f"file {path!r}") from None
    meta = service.manifest(commit).get(path)
    if meta is None:
        raise not_found(f"file {path!r}")
    data = store.get(meta["sha256"])
    too_big = len(data) > 512 * 1024
    text = None if too_big else data.decode("utf-8", "replace")
    rendered = None
    if text is not None:
        rendered = render.markdown_html(text) if path.endswith(".md") and request.query_params.get("plain") != "1" \
            else render.code_html(path, text)
    return _page(request, "blob.html", p, tab="files", commit=service.commit_json(c, commit),
                 file_path=path, meta=meta, rendered=rendered, too_big=too_big, revision=revision,
                 is_code=path.endswith(".py"), **ctx)


@router.get("/{namespace}/{name}/commits", response_class=HTMLResponse)
def commits(request: Request, namespace: str, name: str,
            c: sqlite3.Connection = Depends(deps.conn),
            p: auth.Principal | None = Depends(deps.principal)) -> HTMLResponse:
    ctx = _repo_context(request, c, namespace, name, p)
    return _page(request, "commits.html", p, tab="commits", history=service.history(c, ctx["row"]), **ctx)


@router.post("/{namespace}/{name}/like")
def like(request: Request, namespace: str, name: str, csrf: str = Form(""),
         c: sqlite3.Connection = Depends(deps.conn),
         p: auth.Principal | None = Depends(deps.principal)) -> Response:
    if p is None or p.via != "session":
        return RedirectResponse(f"/login?{urlencode({'next': f'/{namespace}/{name}'})}", status_code=303)
    _check_csrf(request, p, csrf)
    repo = service.get_repo(c, namespace, name)
    service.set_like(c, p, repo, not service.has_liked(c, p, repo))
    return RedirectResponse(f"/{repo['namespace']}/{repo['name']}", status_code=303)


# --------------------------------------------------------------------------
# Installation into the app
# --------------------------------------------------------------------------

def install(app: FastAPI) -> None:
    settings: Settings = app.state.settings
    app.state.views = render.ViewCache()
    css = render.pygments_css()

    @app.middleware("http")
    async def anonymous_csrf(request: Request, call_next):
        is_page = not request.url.path.startswith(("/api/", "/static/")) and "/resolve/" not in request.url.path
        token = request.cookies.get(ANON_CSRF) or secrets.token_urlsafe(24)
        request.state.anon_csrf = token
        response = await call_next(request)
        if is_page and ANON_CSRF not in request.cookies:
            response.set_cookie(ANON_CSRF, token, httponly=True, samesite="lax", secure=settings.secure_cookies)
        return response

    def render_error(request: Request, err: ApiError) -> Response:
        return templates.TemplateResponse(
            request, "error.html",
            {"status": err.status, "message": err.message, "problems": err.problems, "me": None,
             "csrf": "", "settings": settings, "path": request.url.path, "signup": settings.allow_signup,
             "site": _site(settings)},
            status_code=err.status,
        )

    app.state.render_error = render_error

    @app.get("/static/pygments.css", include_in_schema=False)
    def pygments() -> Response:
        return PlainTextResponse(css, media_type="text/css", headers={"Cache-Control": "public, max-age=3600"})

    app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")
    from chemart_hub.routes import admin  # /admin/... and /{ns}/{name}/edit, before the catch-alls

    app.include_router(admin.router)
    app.include_router(router)
