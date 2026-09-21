"""The FastAPI application: `create_app(settings)`."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, Response

from chemart_hub import db
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError
from chemart_hub.routes import api, files
from chemart_hub.storage import BlobStore

#: Applied to every HTML page; raw files set their own stricter policy.
CSP = ("default-src 'self'; script-src 'self'; style-src 'self' https://fonts.googleapis.com; "
       "font-src 'self' https://fonts.gstatic.com; img-src 'self' https: data:; "
       "object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")


def _wants_json(request: Request) -> bool:
    path = request.url.path
    return path.startswith("/api/") or "/resolve/" in path or "text/html" not in request.headers.get("accept", "")


def create_app(settings: Settings | None = None) -> FastAPI:
    settings = settings or Settings.from_env()
    db.migrate(settings.db_path)
    app = FastAPI(title="Chemart Hub", docs_url="/api/docs", redoc_url=None, openapi_url="/api/openapi.json")
    app.state.settings = settings
    app.state.store = BlobStore(settings.blobs_dir)

    @app.exception_handler(ApiError)
    async def api_error(request: Request, err: ApiError) -> Response:
        if _wants_json(request) or not hasattr(app.state, "render_error"):
            return JSONResponse(err.body(), status_code=err.status)
        return app.state.render_error(request, err)

    @app.exception_handler(RequestValidationError)
    async def bad_request(request: Request, err: RequestValidationError) -> Response:
        problems = [f"{'.'.join(map(str, e['loc']))}: {e['msg']}" for e in err.errors()]
        return JSONResponse({"error": "malformed request", "problems": problems}, status_code=422)

    @app.middleware("http")
    async def security_headers(request: Request, call_next):
        response = await call_next(request)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "strict-origin-when-cross-origin")
        response.headers.setdefault("Content-Security-Policy", CSP)
        if settings.secure_cookies:
            response.headers.setdefault("Strict-Transport-Security", "max-age=31536000")
        return response

    app.include_router(api.router)
    app.include_router(files.router)
    from chemart_hub.routes import web

    web.install(app)
    return app
