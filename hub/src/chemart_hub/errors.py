"""The one error type routes raise; API paths answer it as JSON."""

from __future__ import annotations


class ApiError(Exception):
    def __init__(self, status: int, message: str, problems: list[str] | None = None):
        super().__init__(message)
        self.status = status
        self.message = message
        self.problems = problems or []

    def body(self) -> dict:
        out: dict = {"error": self.message}
        if self.problems:
            out["problems"] = self.problems
        return out


def not_found(what: str) -> ApiError:
    return ApiError(404, f"{what} not found")


def unauthorized(message: str = "log in or pass an API token (Authorization: Bearer chm_...)") -> ApiError:
    return ApiError(401, message)


def forbidden(message: str) -> ApiError:
    return ApiError(403, message)


def conflict(message: str) -> ApiError:
    return ApiError(409, message)


def invalid(message: str, problems: list[str] | None = None) -> ApiError:
    return ApiError(422, message, problems)
