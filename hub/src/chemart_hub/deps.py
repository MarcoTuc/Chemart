"""FastAPI dependencies shared by the API and web routes."""

from __future__ import annotations

import sqlite3
from typing import Iterator

from fastapi import Depends, Request

from chemart_hub import auth, db
from chemart_hub.config import Settings
from chemart_hub.storage import BlobStore


def settings(request: Request) -> Settings:
    return request.app.state.settings


def store(request: Request) -> BlobStore:
    return request.app.state.store


def conn(request: Request) -> Iterator[sqlite3.Connection]:
    c = db.connect(request.app.state.settings.db_path)
    try:
        yield c
    finally:
        c.close()


def principal(request: Request, c: sqlite3.Connection = Depends(conn)) -> auth.Principal | None:
    return auth.identify(c, request)


def writer(request: Request, p: auth.Principal | None = Depends(principal)) -> auth.Principal:
    return auth.require_writer(p, request)
