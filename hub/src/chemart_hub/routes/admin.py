"""Curation pages: the superadmin panel (/admin) and each repo's Edit tab.

The panel is for superadmins only. The Edit tab is open to anyone who may
push to the repo (its owner, members of its organisation, superadmins); the
curation controls on it (feature, hide, transfer) are superadmin-only. Every
form posts with a CSRF token, and every content change is a new commit.

Deleting never destroys anything: it archives (see `curation`). The archive
page restores items, or purges them after a typed confirmation.
"""

from __future__ import annotations

import sqlite3
from typing import Any
from urllib.parse import urlencode

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse, RedirectResponse, Response

from chemart.hub import _format

from chemart_hub import auth, curation, deps, service
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError, forbidden, invalid, not_found
from chemart_hub.routes.web import _cards, _check_csrf, _login_redirect, _page, _repo_context, _safe_next
from chemart_hub.storage import BlobStore

router = APIRouter()
PAGE = 50
EDITABLE_MAX = 512 * 1024


def _session(request: Request, p: auth.Principal | None) -> auth.Principal | None:
    return p if p is not None and p.via == "session" else None


def _admin_or_redirect(request: Request, p: auth.Principal | None) -> auth.Principal | Response:
    p = _session(request, p)
    if p is None:
        return _login_redirect(request)
    return curation.require_admin(p)


def _back(next_url: str, fallback: str) -> RedirectResponse:
    return RedirectResponse(_safe_next(next_url) if next_url else fallback, status_code=303)


# --------------------------------------------------------------------------
# The superadmin panel
# --------------------------------------------------------------------------

@router.get("/admin", response_class=HTMLResponse)
def overview(request: Request, c: sqlite3.Connection = Depends(deps.conn),
             p: auth.Principal | None = Depends(deps.principal),
             store: BlobStore = Depends(deps.store)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    count = lambda sql: c.execute(sql).fetchone()[0]  # noqa: E731
    live_acc, live_repo = "archived_at IS NULL", "archived_at IS NULL"
    stats = {
        "users": count(f"SELECT COUNT(*) FROM accounts WHERE is_org = 0 AND {live_acc}"),
        "organisations": count(f"SELECT COUNT(*) FROM accounts WHERE is_org = 1 AND {live_acc}"),
        "superadmins": count(f"SELECT COUNT(*) FROM accounts WHERE is_admin = 1 AND {live_acc}"),
        "suspended": count(f"SELECT COUNT(*) FROM accounts WHERE disabled = 1 AND {live_acc}"),
        "chemistries": count(f"SELECT COUNT(*) FROM repos WHERE repo_type = 'generator' AND head IS NOT NULL AND {live_repo}"),
        "networks": count(f"SELECT COUNT(*) FROM repos WHERE repo_type = 'network' AND head IS NOT NULL AND {live_repo}"),
        "empty repos": count(f"SELECT COUNT(*) FROM repos WHERE head IS NULL AND {live_repo}"),
        "featured": count(f"SELECT COUNT(*) FROM repos WHERE featured = 1 AND {live_repo}"),
        "hidden": count(f"SELECT COUNT(*) FROM repos WHERE hidden = 1 AND {live_repo}"),
        "in the archive": count("SELECT COUNT(*) FROM repos WHERE archived_at IS NOT NULL")
                          + count("SELECT COUNT(*) FROM accounts WHERE archived_at IS NOT NULL"),
        "commits": count("SELECT COUNT(*) FROM commits"),
        "downloads": count(f"SELECT COALESCE(SUM(downloads), 0) FROM repos WHERE {live_repo}"),
        "likes": count("SELECT COUNT(*) FROM likes"),
    }
    return _page(request, "admin_overview.html", p, tab="overview", stats=stats, log=curation.recent_log(c, 15))


@router.get("/admin/repos", response_class=HTMLResponse)
def admin_repos(request: Request, q: str = "", show: str = "all", page: int = 1,
                c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal),
                settings: Settings = Depends(deps.settings)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    where, args = ["r.archived_at IS NULL",
                   "(a.name || '/' || r.name || ' ' || COALESCE(r.title, '')) LIKE ?"], [f"%{q.strip()}%"]
    extra = {"hidden": "r.hidden = 1", "featured": "r.featured = 1", "empty": "r.head IS NULL",
             "code": "r.has_code = 1", "network": "r.repo_type = 'network'",
             "generator": "r.repo_type = 'generator'"}.get(show)
    if extra:
        where.append(extra)
    clause = " WHERE " + " AND ".join(where)
    total = c.execute(f"SELECT COUNT(*) FROM repos r JOIN accounts a ON a.id = r.owner_id{clause}", args).fetchone()[0]
    page = max(1, page)
    rows = c.execute(
        f"SELECT r.*, a.name AS namespace FROM repos r JOIN accounts a ON a.id = r.owner_id{clause} "
        f"ORDER BY r.updated_at DESC LIMIT ? OFFSET ?", [*args, PAGE, (page - 1) * PAGE]).fetchall()
    return _page(request, "admin_repos.html", p, tab="repos", repos=_cards(c, rows, settings), total=total,
                 q=q, show=show, page=page, pages=max(1, (total + PAGE - 1) // PAGE),
                 back=str(request.url.path) + ("?" + request.url.query if request.url.query else ""))


@router.get("/admin/users", response_class=HTMLResponse)
def admin_users(request: Request, q: str = "", c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    rows = c.execute(
        "SELECT a.*, (SELECT COUNT(*) FROM repos r WHERE r.owner_id = a.id AND r.archived_at IS NULL) AS n_repos, "
        "(SELECT COUNT(*) FROM tokens t WHERE t.user_id = a.id AND t.revoked = 0) AS n_tokens, "
        "(SELECT GROUP_CONCAT(u.name, ', ') FROM org_members m JOIN accounts u ON u.id = m.user_id "
        " WHERE m.org_id = a.id AND u.archived_at IS NULL) AS members "
        "FROM accounts a WHERE a.archived_at IS NULL AND "
        "(a.name LIKE ? OR COALESCE(a.email, '') LIKE ? OR COALESCE(a.fullname, '') LIKE ?) "
        "ORDER BY a.is_org, a.name LIMIT 500",
        (f"%{q.strip()}%",) * 3).fetchall()
    return _page(request, "admin_users.html", p, tab="users", accounts=rows, q=q,
                 back=str(request.url.path) + ("?" + request.url.query if request.url.query else ""))


@router.get("/admin/site", response_class=HTMLResponse)
def admin_site(request: Request, c: sqlite3.Connection = Depends(deps.conn),
               p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    stored = {r["key"]: r for r in c.execute("SELECT * FROM site_text")}
    fields = [{"key": k, "label": label, "default": default, "kind": kind, "help": help_,
               "value": curation.site_texts(c)[k], "edited": stored.get(k)}
              for k, (label, default, kind, help_) in curation.SITE_TEXT.items()]
    return _page(request, "admin_site.html", p, tab="site", fields=fields,
                 saved=request.query_params.get("saved"))


@router.post("/admin/site")
async def admin_site_save(request: Request, c: sqlite3.Connection = Depends(deps.conn),
                          p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    form = await request.form()
    _check_csrf(request, admin, str(form.get("csrf", "")))
    changed = 0
    for key in curation.SITE_TEXT:
        if key in form:
            changed += curation.set_site_text(c, admin, key, str(form[key]))
    return RedirectResponse(f"/admin/site?saved={changed}", status_code=303)


@router.get("/admin/log", response_class=HTMLResponse)
def admin_log(request: Request, c: sqlite3.Connection = Depends(deps.conn),
              p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    return _page(request, "admin_log.html", p, tab="log", log=curation.recent_log(c, 500))


# --- actions --------------------------------------------------------------

def _flag(value: str | None) -> bool | None:
    return None if value in (None, "") else value in ("1", "true", "on", "yes")


@router.post("/admin/repos/{namespace}/{name}/flags")
def repo_flags(request: Request, namespace: str, name: str, featured: str = Form(""), hidden: str = Form(""),
               csrf: str = Form(""), next: str = Form(""),
               c: sqlite3.Connection = Depends(deps.conn),
               p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    curation.set_flags(c, admin, service.get_repo(c, namespace, name),
                       featured=_flag(featured), hidden=_flag(hidden))
    return _back(next, f"/{namespace}/{name}/edit")


@router.post("/admin/repos/{namespace}/{name}/transfer")
def repo_transfer(request: Request, namespace: str, name: str, owner: str = Form(...),
                  csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                  p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    repo = curation.transfer(c, admin, service.get_repo(c, namespace, name), owner.strip())
    return RedirectResponse(f"/{repo['namespace']}/{repo['name']}/edit", status_code=303)


@router.post("/admin/users/{name}/superadmin")
def user_superadmin(request: Request, name: str, flag: str = Form(...), csrf: str = Form(""),
                    next: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                    p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    curation.set_admin(c, admin, name, bool(_flag(flag)))
    return _back(next, "/admin/users")


@router.post("/admin/users/{name}/suspend")
def user_suspend(request: Request, name: str, flag: str = Form(...), csrf: str = Form(""),
                 next: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                 p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    curation.set_disabled(c, admin, name, bool(_flag(flag)))
    return _back(next, "/admin/users")


@router.post("/admin/users/{name}/password")
def user_password(request: Request, name: str, password: str = Form(...), csrf: str = Form(""),
                  next: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                  p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    curation.set_password(c, admin, name, password)
    return _back(next, "/admin/users")


@router.post("/admin/users/{name}/delete")
def user_delete(request: Request, name: str, confirm: str = Form(""), csrf: str = Form(""),
                c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    if confirm.strip().lower() != name.lower():
        raise invalid(f"type the account name ({name}) to confirm the deletion")
    curation.delete_account(c, admin, name)
    return RedirectResponse("/admin/archive?" + urlencode({"archived": name.lower()}), status_code=303)


@router.post("/admin/orgs")
def org_create(request: Request, name: str = Form(...), owner: str = Form(...), csrf: str = Form(""),
               c: sqlite3.Connection = Depends(deps.conn),
               p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    org = curation.create_org(c, admin, name, owner.strip())
    return RedirectResponse(f"/settings/profile/{org['name']}", status_code=303)


# --------------------------------------------------------------------------
# A repo's Edit tab
# --------------------------------------------------------------------------

def _editor(request: Request, c: sqlite3.Connection, namespace: str, name: str,
            p: auth.Principal | None) -> tuple[dict[str, Any], auth.Principal] | Response:
    session = _session(request, p)
    if session is None:
        return _login_redirect(request)
    ctx = _repo_context(request, c, namespace, name, session)
    if not ctx["can_write"]:
        raise forbidden(f"you cannot edit {ctx['repo']['id']}")
    return ctx, session


def _hub_block(files: dict[str, bytes]) -> dict[str, Any]:
    if "chemart.yaml" not in files:
        return {}
    try:
        doc = _format.load_yaml(files["chemart.yaml"])
    except _format.FormatError:
        return {}
    hub = doc.get("hub") if isinstance(doc, dict) else None
    return hub if isinstance(hub, dict) else {}


def _edit_page(request: Request, c: sqlite3.Connection, store: BlobStore, ctx: dict[str, Any],
               p: auth.Principal, status: int = 200, **extra: Any) -> HTMLResponse:
    files = curation.head_files(c, store, ctx["row"])
    hub = _hub_block(files)
    meta = {"title": hub.get("title", ""), "description": hub.get("description", ""),
            "tags": ", ".join(hub.get("tags", []) if isinstance(hub.get("tags"), list) else []),
            "license": hub.get("license", "")}
    meta.update(extra.pop("meta", {}))
    listing = sorted((path, len(data)) for path, data in files.items())
    return _page(request, "repo_edit.html", p, status, tab="edit", files=listing, meta=meta,
                 error=extra.pop("error", None), problems=extra.pop("problems", []), **ctx, **extra)


@router.get("/{namespace}/{name}/edit", response_class=HTMLResponse)
def edit_page(request: Request, namespace: str, name: str, c: sqlite3.Connection = Depends(deps.conn),
              p: auth.Principal | None = Depends(deps.principal),
              store: BlobStore = Depends(deps.store)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    return _edit_page(request, c, store, ctx, session, saved=request.query_params.get("saved"))


@router.post("/{namespace}/{name}/edit/metadata", response_class=HTMLResponse)
def edit_metadata(request: Request, namespace: str, name: str, title: str = Form(""),
                  description: str = Form(""), tags: str = Form(""), license: str = Form(""),
                  csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                  p: auth.Principal | None = Depends(deps.principal),
                  store: BlobStore = Depends(deps.store),
                  settings: Settings = Depends(deps.settings)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    _check_csrf(request, session, csrf)
    if ctx["row"]["head"] is None:
        raise invalid("push something to this repo before editing its card")
    try:
        curation.update_metadata(c, store, settings, session, ctx["row"], title=title,
                                 description=description, tags=tags, license=license)
    except ApiError as err:
        return _edit_page(request, c, store, ctx, session, 400, error=err.message, problems=err.problems,
                          meta={"title": title, "description": description, "tags": tags, "license": license})
    return RedirectResponse(f"/{namespace}/{name}/edit?saved=card", status_code=303)


@router.get("/{namespace}/{name}/edit/file", response_class=HTMLResponse)
def file_editor(request: Request, namespace: str, name: str, path: str = "",
                c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal),
                store: BlobStore = Depends(deps.store)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    content = ""
    if path:
        files = curation.head_files(c, store, ctx["row"])
        if path not in files:
            raise not_found(f"file {path!r}")
        if len(files[path]) > EDITABLE_MAX:
            raise invalid(f"{path} is too large to edit in the browser; push it with chemart instead")
        content = files[path].decode("utf-8", "replace")
    return _page(request, "file_edit.html", session, tab="edit", path=path, original=path,
                 content=content, message="", error=None, problems=[], **ctx)


@router.post("/{namespace}/{name}/edit/file", response_class=HTMLResponse)
def file_save(request: Request, namespace: str, name: str, path: str = Form(...), original: str = Form(""),
              content: str = Form(""), message: str = Form(""), csrf: str = Form(""),
              c: sqlite3.Connection = Depends(deps.conn),
              p: auth.Principal | None = Depends(deps.principal),
              store: BlobStore = Depends(deps.store),
              settings: Settings = Depends(deps.settings)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    _check_csrf(request, session, csrf)
    path, original = path.strip(), original.strip()
    changes: dict[str, bytes | None] = {path: curation.text_from_form(content)}
    if original and original != path:
        changes[original] = None                      # a rename
    verb = "Create" if not original else ("Rename" if original != path else "Edit")
    try:
        result = curation.edit_files(c, store, settings, session, ctx["row"], changes,
                                     message.strip() or f"{verb} {path}")
    except ApiError as err:
        return _page(request, "file_edit.html", session, 400, tab="edit", path=path, original=original,
                     content=content, message=message, error=err.message, problems=err.problems, **ctx)
    if result.get("unchanged"):
        return RedirectResponse(f"/{namespace}/{name}/edit?saved=nothing", status_code=303)
    return RedirectResponse(f"/{namespace}/{name}/blob/main/{path}", status_code=303)


@router.post("/{namespace}/{name}/edit/file/delete", response_class=HTMLResponse)
def file_delete(request: Request, namespace: str, name: str, path: str = Form(...), csrf: str = Form(""),
                c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal),
                store: BlobStore = Depends(deps.store),
                settings: Settings = Depends(deps.settings)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    _check_csrf(request, session, csrf)
    try:
        curation.edit_files(c, store, settings, session, ctx["row"], {path: None}, f"Delete {path}")
    except ApiError as err:
        return _edit_page(request, c, store, ctx, session, 400, error=err.message, problems=err.problems)
    return RedirectResponse(f"/{namespace}/{name}/edit?saved=deleted", status_code=303)


@router.post("/{namespace}/{name}/edit/restore", response_class=HTMLResponse)
def restore(request: Request, namespace: str, name: str, commit: str = Form(...), csrf: str = Form(""),
            c: sqlite3.Connection = Depends(deps.conn),
            p: auth.Principal | None = Depends(deps.principal),
            store: BlobStore = Depends(deps.store),
            settings: Settings = Depends(deps.settings)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    _check_csrf(request, session, csrf)
    curation.restore(c, store, settings, session, ctx["row"], commit)
    return RedirectResponse(f"/{namespace}/{name}/commits", status_code=303)


@router.post("/{namespace}/{name}/edit/delete")
def delete_repo(request: Request, namespace: str, name: str, confirm: str = Form(""), csrf: str = Form(""),
                c: sqlite3.Connection = Depends(deps.conn),
                p: auth.Principal | None = Depends(deps.principal)) -> Response:
    found = _editor(request, c, namespace, name, p)
    if isinstance(found, Response):
        return found
    ctx, session = found
    _check_csrf(request, session, csrf)
    if confirm.strip() != ctx["repo"]["id"]:
        raise invalid(f"type {ctx['repo']['id']} to confirm the deletion")
    curation.delete_repo(c, session, ctx["row"])
    return RedirectResponse(f"/{ctx['repo']['namespace']}?" + urlencode({"deleted": ctx["repo"]["name"]}),
                            status_code=303)


# --------------------------------------------------------------------------
# The archive: deleted things wait here until a superadmin purges them
# --------------------------------------------------------------------------

@router.get("/admin/archive", response_class=HTMLResponse)
def archive_page(request: Request, c: sqlite3.Connection = Depends(deps.conn),
                 p: auth.Principal | None = Depends(deps.principal),
                 store: BlobStore = Depends(deps.store)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    q = request.query_params
    done = {k: q.get(k) for k in ("purged", "files", "freed", "restored", "skipped", "archived")
            if q.get(k) is not None}
    return _page(request, "admin_archive.html", p, tab="archive", archive=curation.archive(c, store), done=done)


@router.post("/admin/archive/repos/{repo_id}/restore")
def archive_restore_repo(request: Request, repo_id: int, csrf: str = Form(""),
                         c: sqlite3.Connection = Depends(deps.conn),
                         p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    repo = curation.restore_repo(c, admin, repo_id)
    return RedirectResponse("/admin/archive?" + urlencode({"restored": f"{repo['namespace']}/{repo['name']}"}),
                            status_code=303)


@router.post("/admin/archive/accounts/{account_id}/restore")
def archive_restore_account(request: Request, account_id: int, csrf: str = Form(""),
                            c: sqlite3.Connection = Depends(deps.conn),
                            p: auth.Principal | None = Depends(deps.principal)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    name = c.execute("SELECT name FROM accounts WHERE id = ?", (account_id,)).fetchone()
    skipped = curation.restore_account(c, admin, account_id)
    params = {"restored": name["name"] if name else f"#{account_id}"}
    if skipped:
        params["skipped"] = ", ".join(skipped)
    return RedirectResponse("/admin/archive?" + urlencode(params), status_code=303)


def _purge_target(c: sqlite3.Connection, store: BlobStore, kind: str, item_id: int | None) -> dict[str, Any]:
    """What a purge would remove, for the confirmation page."""
    archive = curation.archive(c, store)
    if kind == "all":
        return {"kind": "all", "id": "", "label": "everything in the archive",
                "repos": [r["id_text"] for r in archive["repos"]], "accounts": [a["name"] for a in archive["accounts"]],
                "freed": archive["freed_by_all"]}
    items = archive["repos"] if kind == "repo" else archive["accounts"] if kind == "account" else None
    if items is None:
        raise invalid("kind must be repo, account or all")
    item = next((i for i in items if i["id"] == item_id), None)
    if item is None:
        raise not_found(f"archived {kind} #{item_id}")
    if kind == "repo":
        return {"kind": "repo", "id": item_id, "label": item["id_text"], "repos": [item["id_text"]], "accounts": [],
                "freed": item["freed"], "commits": item["n_commits"]}
    repos = [f"{item['name']}/{r['name']}" for r in c.execute(
        "SELECT name FROM repos WHERE owner_id = ? ORDER BY name", (item_id,))]
    return {"kind": "account", "id": item_id, "label": item["name"], "repos": repos, "accounts": [item["name"]],
            "freed": item["freed"]}


@router.get("/admin/archive/purge", response_class=HTMLResponse)
def purge_confirm(request: Request, kind: str, id: int | None = None,
                  c: sqlite3.Connection = Depends(deps.conn),
                  p: auth.Principal | None = Depends(deps.principal),
                  store: BlobStore = Depends(deps.store)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    return _page(request, "admin_purge.html", p, tab="archive", target=_purge_target(c, store, kind, id))


@router.post("/admin/archive/purge")
def purge(request: Request, kind: str = Form(...), id: str = Form(""), confirm: str = Form(""),
          csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
          p: auth.Principal | None = Depends(deps.principal),
          store: BlobStore = Depends(deps.store)) -> Response:
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    item_id = int(id) if id.isdigit() else None
    target = _purge_target(c, store, kind, item_id)
    expected = "purge everything" if kind == "all" else target["label"]
    if confirm.strip() != expected:
        raise invalid(f"type “{expected}” to confirm the purge")
    if kind == "all":
        result = curation.purge_all(c, store, admin)
    elif kind == "repo":
        result = curation.purge_repo(c, store, admin, item_id)
    else:
        result = curation.purge_account(c, store, admin, item_id)
    request.app.state.views.clear()          # no rendered copy of a purged repo stays in memory
    return RedirectResponse("/admin/archive?" + urlencode(
        {"purged": target["label"], "files": result["files_removed"], "freed": result["bytes_freed"]}),
        status_code=303)


@router.post("/admin/archive/gc")
def storage_cleanup(request: Request, csrf: str = Form(""), c: sqlite3.Connection = Depends(deps.conn),
                    p: auth.Principal | None = Depends(deps.principal),
                    store: BlobStore = Depends(deps.store)) -> Response:
    """Remove orphaned files (from interrupted uploads); touches nothing in the archive."""
    admin = _admin_or_redirect(request, p)
    if isinstance(admin, Response):
        return admin
    _check_csrf(request, admin, csrf)
    result = curation.collect_garbage(c, store)
    curation.log(c, admin, "clean up storage", "blobs", f"{result['files_removed']} files, {result['bytes_freed']} bytes")
    return RedirectResponse("/admin/archive?" + urlencode(
        {"purged": "orphaned files", "files": result["files_removed"], "freed": result["bytes_freed"]}),
        status_code=303)
