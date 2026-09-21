"""Content-addressed blob store and commit ids.

A blob lives at ``blobs/ab/cd/<sha256>`` and is never modified, so a file
shared by many commits (or many repos) is stored once. Blobs are only ever
removed by a purge (`curation.collect_garbage`), and only when no remaining
commit uses them. A commit id is the
sha256 of the canonical JSON of the commit record, so it pins the exact file
contents, the parent and the metadata.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import AsyncIterator

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class BlobTooLarge(ValueError):
    pass


class BlobMismatch(ValueError):
    pass


class BlobStore:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def path(self, sha: str) -> Path:
        if not SHA256_RE.match(sha):
            raise ValueError(f"not a sha256: {sha!r}")
        return self.root / sha[:2] / sha[2:4] / sha

    def has(self, sha: str) -> bool:
        return self.path(sha).is_file()

    def get(self, sha: str) -> bytes:
        return self.path(sha).read_bytes()

    def size(self, sha: str) -> int:
        return self.path(sha).stat().st_size

    def _finish(self, tmp: Path, sha: str, digest: str) -> None:
        if digest != sha:
            tmp.unlink(missing_ok=True)
            raise BlobMismatch(f"content hashes to {digest}, not {sha}")
        target = self.path(sha)
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(tmp, target)

    def touch(self, sha: str) -> None:
        """Mark a blob as just uploaded, so a concurrent purge leaves it alone."""
        os.utime(self.path(sha))

    def blobs(self):
        """Every stored blob: (sha, path)."""
        for path in self.root.glob("??/??/*"):
            if path.is_file() and SHA256_RE.match(path.name):
                yield path.name, path

    def partial_uploads(self):
        """Temporary files of uploads that never finished."""
        return [p for p in self.root.glob(".upload-*") if p.is_file()]

    def remove(self, sha: str) -> int:
        """Delete one blob from disk; return the bytes freed."""
        path = self.path(sha)
        try:
            size = path.stat().st_size
            path.unlink()
        except FileNotFoundError:
            return 0
        for parent in (path.parent, path.parent.parent):   # drop empty fan-out folders
            try:
                parent.rmdir()
            except OSError:
                break
        return size

    def put_bytes(self, data: bytes) -> str:
        sha = hashlib.sha256(data).hexdigest()
        if self.has(sha):
            self.touch(sha)
        else:
            fd, tmp = tempfile.mkstemp(dir=self.root, prefix=".upload-")
            with os.fdopen(fd, "wb") as f:
                f.write(data)
            self._finish(Path(tmp), sha, sha)
        return sha

    async def put_stream(self, sha: str, chunks: AsyncIterator[bytes], limit: int) -> int:
        """Store a streamed upload that must hash to `sha`; return its size."""
        self.path(sha)  # validates the name
        digest = hashlib.sha256()
        size = 0
        fd, tmp = tempfile.mkstemp(dir=self.root, prefix=".upload-")
        try:
            with os.fdopen(fd, "wb") as f:
                async for chunk in chunks:
                    size += len(chunk)
                    if size > limit:
                        raise BlobTooLarge(f"upload exceeds {limit} bytes")
                    digest.update(chunk)
                    f.write(chunk)
            self._finish(Path(tmp), sha, digest.hexdigest())
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
        return size


def commit_id(*, repo: str, parent: str | None, files: list[list], message: str,
              author: str, created_at: str) -> str:
    record = {
        "format": 1, "repo": repo, "parent": parent, "files": sorted(files),
        "message": message, "author": author, "created_at": created_at,
    }
    canonical = json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode()).hexdigest()
