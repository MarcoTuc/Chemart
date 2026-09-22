"""The simulation pit: pick a chemistry, run it, watch it, measure it.

    chemart-hub pit [--port 8765] [--no-browser]

A small local app, separate from the hub site: no database, no accounts. It
runs chemistries on this machine with native Python, so it binds to
127.0.0.1 only and is never part of the public static hub. Every request must
come to a localhost Host name (against DNS rebinding) and every POST must
carry the `X-Chemart-Pit` header, which a cross-site form or script cannot
send without a CORS preflight the pit never answers.

Runs stream as newline-delimited JSON (one message per line): the network,
then frames as they come, then `done`. The last few runs stay in memory and
can be downloaded as a whole trajectory.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterator

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

HERE = Path(__file__).resolve().parent
STATIC = HERE / "static"
PAGE = HERE / "templates" / "pit.html"
HEADER = "x-chemart-pit"
LOCAL_HOSTS = {"127.0.0.1", "localhost", "[::1]"}
CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; "
       "connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
#: Runs kept in memory for download.
KEEP_RUNS = 5
#: Species drawn per frame in a stream; the rest are summed as "(other)".
TOP = 12


def create_pit_app() -> FastAPI:
    app = FastAPI(title="Chemart pit", docs_url=None, redoc_url=None, openapi_url=None)
    runs: OrderedDict[str, dict] = OrderedDict()
    counter = itertools.count(1)

    @app.middleware("http")
    async def guard(request: Request, call_next):
        host = (request.headers.get("host") or "").rsplit(":", 1)[0]
        client = request.client.host if request.client else ""
        if host not in LOCAL_HOSTS or client not in ("127.0.0.1", "::1", "testclient"):
            return JSONResponse({"error": "the pit only serves this machine"}, status_code=403)
        if request.method == "POST" and request.headers.get(HEADER) != "1":
            return JSONResponse({"error": f"missing {HEADER} header"}, status_code=403)
        response = await call_next(request)
        response.headers.setdefault("Content-Security-Policy", CSP)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        return response

    app.mount("/static", StaticFiles(directory=STATIC), name="static")

    @app.get("/")
    def page():
        return FileResponse(PAGE, media_type="text/html")

    @app.get("/api/chemistries")
    def chemistries():
        return _chemistries()

    @app.get("/api/chemistry/{cid}")
    def chemistry(cid: str):
        import chemart
        from chemart import measures

        try:
            info = chemart.describe_chemistry(cid)
        except ValueError as err:
            raise HTTPException(404, str(err)) from None
        info["measures"] = measures.describe()
        return info

    @app.post("/api/run")
    async def run(request: Request):
        body = await request.json()
        run_id = str(next(counter))

        def stream() -> Iterator[bytes]:
            for message in _run(body, run_id, runs):
                yield _line(message)

        return StreamingResponse(stream(), media_type="application/x-ndjson")

    @app.get("/api/runs/{run_id}.json")
    def download(run_id: str):
        if run_id not in runs:
            raise HTTPException(404, "that run is no longer kept")
        return JSONResponse(runs[run_id], headers={
            "Content-Disposition": f'attachment; filename="chemart-run-{run_id}.json"'})

    @app.post("/api/measure")
    async def measure(request: Request):
        body = await request.json()
        try:
            net = _network(body)
            return _measures(net, body.get("names"), body.get("cost", "cheap"))
        except Exception as err:  # the page shows the message
            return JSONResponse({"error": _message(err)}, status_code=400)

    @app.post("/api/sweep")
    async def sweep(request: Request):
        body = await request.json()

        def stream() -> Iterator[bytes]:
            for message in _sweep(body):
                yield _line(message)

        return StreamingResponse(stream(), media_type="application/x-ndjson")

    return app


# --------------------------------------------------------------------------
# What the routes do
# --------------------------------------------------------------------------

_CHEMISTRIES: list[dict] | None = None


def _chemistries() -> list[dict]:
    """Every catalogued chemistry with its type and faces (computed once)."""
    global _CHEMISTRIES
    if _CHEMISTRIES is None:
        import chemart
        from chemart import api

        out = []
        for row in chemart.list_chemistries():
            c = api._entry(row["id"])
            try:
                faces = api.faces(c)
            except Exception:  # a module that fails to import is shown without faces
                faces = []
            out.append({"id": c.id, "name": c.name, "type": c.type, "faces": faces,
                        "clock": c.clock, "family": c.family, "provides": c.provides})
        _CHEMISTRIES = out
    return _CHEMISTRIES


def _network(body: dict):
    import chemart

    return chemart.generate_network(body["chemistry"], body.get("seed"), **(body.get("params") or {}))


def _measures(net, names=None, cost="cheap") -> dict:
    from chemart import measures

    values = measures.measure(net, names or None, cost=cost)
    reasons = {k: v for k, v in measures.applicable(net, names or None, cost=cost).items() if v}
    return {"values": values, "reasons": reasons}


def _summary(net) -> dict:
    return {"chemistry": net.chemistry, "species": len(net.species), "reactions": len(net.reactions),
            "status": net.status, "provides": net.provides, "text": net.summary()}


def _top(state: dict, keep: list[str] | None = None, k: int = TOP) -> dict:
    """The k most abundant species (or the `keep` list) and the rest as (other)."""
    names = keep if keep is not None else sorted(state, key=state.get, reverse=True)[:k]
    out = {s: state.get(s, 0.0) for s in names}
    rest = sum(v for s, v in state.items() if s not in out)
    if rest:
        out["(other)"] = rest
    return out


def _run(body: dict, run_id: str, runs: OrderedDict) -> Iterator[dict]:
    from chemart import api, measures, simulate
    from chemart.trajectory import Frame, Trajectory

    method = body.get("method", "ode")
    tracked = body.get("measures") or ["richness", "shannon", "dominance", "population"]
    try:
        if method in ("ode", "ssa"):
            net = _network(body)
            yield {"type": "network", "summary": _summary(net), "measures": _measures(net)}
            options = dict(rates=_spec(body.get("rates")), x0=_spec(body.get("x0")),
                           points=int(body.get("points", 200)), seed=body.get("seed"),
                           fill_only=bool(body.get("fill_only", False)))
            if method == "ode":
                traj = simulate.ode(net, float(body.get("t_end", 40)), **options)
            else:
                traj = simulate.ssa(net, float(body.get("t_end", 40)), volume=float(body.get("volume", 1.0)),
                                    **options)
            peak: dict[str, float] = {}
            for f in traj.frames:
                for s, v in f.state.items():
                    peak[s] = max(peak.get(s, 0.0), abs(v))
            keep = sorted(peak, key=peak.get, reverse=True)[: int(body.get("top", TOP))]
            for row, f in zip(measures.track(traj.frames, [m for m in tracked if _is_state(m)]), traj.frames):
                yield {"type": "frame", "t": f.t, "state": _top(f.state, keep),
                       "measures": {k: v for k, v in row.items() if k != "t"}}
        elif method == "evolve":
            c = api._entry(body["chemistry"])
            run = api.run_evolver(c, api.evolver_for(c), body.get("seed"), body.get("params") or {},
                                  int(body.get("every", 1)))
            frames: list[Frame] = []
            window = body.get("window", 1)
            chosen = [m for m in tracked if m in measures.REGISTRY]
            yield {"type": "start", "clock": c.clock}
            while True:
                try:
                    frame = next(run)
                except StopIteration as stop:
                    net = stop.value
                    break
                frames.append(frame)
                row = next(measures.track([frame], [m for m in chosen if _is_state(m)]))
                if window and any(not _is_state(m) for m in chosen):
                    recent = frames[-int(window):]
                    fired = [e for f in recent for e in f.fired]
                    if fired:
                        sub = measures.fired_network(fired)
                        row.update(measures.measure(sub, [m for m in chosen if not _is_state(m)
                                                          and measures.REGISTRY[m].input == "network"]))
                yield {"type": "frame", "t": frame.t, "state": _top(frame.state, k=int(body.get("top", TOP))),
                       "observables": frame.observables,
                       "measures": {k: v for k, v in row.items() if k != "t"}}
            traj = Trajectory(network=net, frames=frames, method="evolve", clock=c.clock or "steps",
                              settings={"seed": body.get("seed"), "params": net.params})
            yield {"type": "network", "summary": _summary(net), "measures": _measures(net)}
        else:
            raise ValueError(f"method must be ode, ssa or evolve, got {method!r}")
    except Exception as err:
        yield {"type": "error", "message": _message(err)}
        return
    runs[run_id] = traj.to_dict()
    while len(runs) > KEEP_RUNS:
        runs.popitem(last=False)
    yield {"type": "done", "run": run_id, "frames": len(traj.frames),
           "trajectory_measures": measures.measure(traj, [m for m in measures.REGISTRY
                                                          if measures.REGISTRY[m].input == "trajectory"])}


def _sweep(body: dict) -> Iterator[dict]:
    """One row per value of the swept argument and seed, as they are measured."""
    import chemart
    from chemart import measures

    try:
        name, values = body["param"], body["values"]
        seeds = body.get("seeds", [0, 1, 2])
        fixed = body.get("params") or {}
        fixed.pop(name, None)
        for value in values:
            for seed in seeds:
                net = chemart.generate_network(body["chemistry"], seed, **fixed, **{name: value})
                row = {name: value, "seed": seed}
                for k, v in measures.measure(net, body.get("names") or None).items():
                    if isinstance(v, dict):
                        row.update({f"{k}.{kk}": vv for kk, vv in v.items()})
                    else:
                        row[k] = v
                yield {"type": "row", "row": row}
    except Exception as err:
        yield {"type": "error", "message": _message(err)}
        return
    yield {"type": "done"}


def _is_state(name: str) -> bool:
    from chemart.measures import REGISTRY

    return name in REGISTRY and REGISTRY[name].input == "state"


def _spec(value: Any) -> Any:
    """Rate and state specs arrive as JSON; file contents are sent as tables."""
    return None if value in (None, "", {}) else value


def _message(err: Exception) -> str:
    return f"{type(err).__name__}: {err}"


def _line(message: dict) -> bytes:
    """One NDJSON line; non-finite floats become null so the browser can parse it."""
    return (json.dumps(_clean(message), ensure_ascii=False) + "\n").encode()


def _clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_clean(v) for v in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()                       # numpy scalars
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)
