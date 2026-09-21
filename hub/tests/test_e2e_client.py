"""End to end: the real `chemart` client over real HTTP to a real server.

A uvicorn server runs in a thread on a free port, with a temporary data
folder; CHEMART_HUB_URL and CHEMART_HOME point the client at it.
"""

from __future__ import annotations

import json
import socket
import textwrap
import threading
import time
from pathlib import Path

import pytest
import uvicorn

import chemart
from chemart import api, hub
from chemart.cli import main as cli_main
from chemart.hub import _format
from chemart_hub import auth, db, service
from chemart_hub.app import create_app
from chemart_hub.config import Settings

from conftest import TINY_GENERATOR, TINY_YAML, network_files, tiny_files


class LiveHub:
    def __init__(self, url: str, settings: Settings, monkeypatch):
        self.url = url
        self.settings = settings
        self.monkeypatch = monkeypatch
        self.tokens: dict[str, str] = {}

    def user(self, name: str, *, login: bool = True, scope: str = "write") -> str:
        c = db.connect(self.settings.db_path)
        try:
            row = service.account(c, name) or service.create_user(c, name, "correct horse battery")
            self.tokens[name] = auth.create_token(c, row["id"], "test", scope)
        finally:
            c.close()
        if login:
            self.as_user(name)
        return self.tokens[name]

    def as_user(self, name: str | None) -> None:
        if name is None:
            self.monkeypatch.delenv("CHEMART_HUB_TOKEN", raising=False)
        else:
            self.monkeypatch.setenv("CHEMART_HUB_TOKEN", self.tokens[name])


@pytest.fixture
def live(tmp_path: Path, monkeypatch):
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    url = f"http://127.0.0.1:{port}"
    settings = Settings(tmp_path / "data", public_url=url)
    server = uvicorn.Server(uvicorn.Config(create_app(settings), log_level="warning", lifespan="off"))
    thread = threading.Thread(target=server.run, kwargs={"sockets": [sock]}, daemon=True)
    thread.start()
    deadline = time.time() + 10
    while not server.started:
        assert time.time() < deadline, "server did not start"
        time.sleep(0.02)
    monkeypatch.setenv("CHEMART_HUB_URL", url)
    monkeypatch.setenv("CHEMART_HOME", str(tmp_path / "home"))
    monkeypatch.delenv("CHEMART_HUB_TOKEN", raising=False)
    monkeypatch.delenv("CHEMART_HUB_OFFLINE", raising=False)
    live = LiveHub(url, settings, monkeypatch)
    live.server = server
    yield live
    server.should_exit = True
    thread.join(timeout=10)
    sock.close()


def write_tiny(folder: Path, generator: str = TINY_GENERATOR, yaml_text: str = TINY_YAML) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "chemart.yaml").write_text(yaml_text)
    (folder / "generator.py").write_text(generator)
    (folder / "README.md").write_text("# Tiny chem\n")
    return folder


# --------------------------------------------------------------------------

def test_network_round_trip(live):
    live.user("bob")
    net = chemart.generate_network("brusselator", seed=3)
    result = net.push_to_hub("bob/bruss", title="Brusselator", tags=["oscillator"])
    assert result["url"] == f"{live.url}/bob/bruss" and not result["unchanged"]
    live.as_user(None)                                       # reading needs no login
    assert chemart.load_network("bob/bruss") == net
    assert net.push_to_hub  # still a method
    live.as_user("bob")
    assert net.push_to_hub("bob/bruss", title="Brusselator", tags=["oscillator"])["unchanged"] is True


def test_generator_gate_and_provenance(live, tmp_path):
    live.user("alice")
    result = hub.push_generator(write_tiny(tmp_path / "tiny-chem"))
    commit = result["commit"]
    assert result["url"] == f"{live.url}/alice/tiny-chem"

    info = chemart.describe_chemistry("alice/tiny-chem")     # metadata only: no trust needed
    assert info["name"] == "Tiny autocatalytic chemistry"
    assert info["hub"]["commit"] == commit and info["hub"]["has_code"] is True
    assert info["params"]["properties"]["k"]["default"] == 1.0

    with pytest.raises(ValueError, match="trust_remote_code=True") as err:
        chemart.generate_network("alice/tiny-chem", seed=0)
    assert commit[:12] in str(err.value)

    with pytest.warns(UserWarning, match="pass revision="):
        net = chemart.generate_network("alice/tiny-chem", seed=0, trust_remote_code=True, k=2.0)
    assert net.chemistry == f"alice/tiny-chem@{commit}"
    assert net.params == {"k": 2.0, "d": 0.1}
    assert net.reactions[0].rate == {"law": "mass-action", "k": 2.0}

    again = chemart.generate_network(net.chemistry, seed=net.seed, trust_remote_code=True, **net.params)
    assert again == net                                      # provenance reproduces the network

    _, entry = _format.parse_chemart_yaml(TINY_YAML.encode())
    from conftest import tiny_generate
    local = api.run_generator(entry, tiny_generate, 0, {"k": 2.0})
    assert local.to_dict() | {"chemistry": ""} == net.to_dict() | {"chemistry": ""}


def test_revisions_pin_code(live, tmp_path):
    live.user("alice")
    folder = write_tiny(tmp_path / "tiny-chem")
    first = hub.push_generator(folder)["commit"]
    write_tiny(folder, TINY_GENERATOR.replace('("X -> ", p.d)', '("X -> ", 2 * p.d)'))
    second = hub.push_generator(folder)["commit"]
    assert first != second

    old = chemart.generate_network(f"alice/tiny-chem@{first[:7]}", trust_remote_code=True)
    new = chemart.generate_network("alice/tiny-chem", revision=second, trust_remote_code=True)
    head = chemart.generate_network("alice/tiny-chem", trust_remote_code=True)
    assert old.reactions[1].rate["k"] == 0.1 and new.reactions[1].rate["k"] == 0.2
    assert head == new
    assert old.chemistry == f"alice/tiny-chem@{first}"


def test_offline_and_unreachable_hub_use_the_cache(live, monkeypatch):
    live.user("bob")
    net = chemart.generate_network("brusselator", seed=1)
    net.push_to_hub("bob/bruss")
    assert chemart.load_network("bob/bruss") == net          # now cached

    monkeypatch.setenv("CHEMART_HUB_OFFLINE", "1")
    assert chemart.load_network("bob/bruss") == net
    with pytest.raises(hub.HubConnectionError, match="not in the cache"):
        chemart.load_network("bob/never-fetched")
    monkeypatch.delenv("CHEMART_HUB_OFFLINE")

    live.server.should_exit = True
    time.sleep(0.5)
    with pytest.warns(UserWarning, match="unreachable"):
        assert chemart.load_network("bob/bruss") == net
    with pytest.raises(hub.HubConnectionError, match="cannot reach"):
        hub.repo_info("bob/bruss")


def test_uploaded_code_never_runs_until_trusted(live, tmp_path):
    """The canary: pushing and describing never import the code."""
    live.user("alice")
    sentinel = tmp_path / "canary-ran"
    evil = TINY_GENERATOR + textwrap.dedent(f"""
        import pathlib
        pathlib.Path({str(sentinel)!r}).write_text("ran")
        """)
    files = tiny_files(**{"generator.py": evil.encode()})
    hub.upload_files("alice/tiny-chem", files, repo_type="generator")
    chemart.describe_chemistry("alice/tiny-chem")
    hub.snapshot_download("alice/tiny-chem")
    with pytest.raises(ValueError):
        chemart.generate_network("alice/tiny-chem")
    assert not sentinel.exists(), "uploaded code ran without trust_remote_code"

    chemart.generate_network("alice/tiny-chem", trust_remote_code=True)
    assert sentinel.exists()                                 # the gate was the only thing in the way


def test_helpers_dataclasses_and_requirements(live, tmp_path):
    live.user("alice")
    folder = write_tiny(tmp_path / "tiny-chem", textwrap.dedent('''\
        from dataclasses import dataclass

        from chemart.helpers.explicit import network
        from .rules import REACTIONS


        @dataclass
        class Rate:
            k: float


        def generate(p, rng):
            return network([(r, Rate(p.k if i == 0 else p.d).k) for i, r in enumerate(REACTIONS)],
                           species=["F", "X"])
        '''))
    (folder / "rules.py").write_text('REACTIONS = ["F + X -> 2 X", "X -> "]\n')
    (folder / "notes.bin").write_bytes(b"\x00")            # not a hub file type: skipped
    hub.push_generator(folder)
    net = chemart.generate_network("alice/tiny-chem", trust_remote_code=True)
    assert len(net.reactions) == 2

    write_tiny(folder, yaml_text=TINY_YAML.replace("  license: MIT", "  license: MIT\n  requires: [no_such_module_xyz]"))
    with pytest.raises(hub.HubValidationError, match="no_such_module_xyz"):
        hub.push_generator(folder)                           # the local check imports nothing missing


def test_type_confusion_is_explained(live, tmp_path):
    live.user("alice")
    chemart.generate_network("brusselator").push_to_hub("alice/snap")
    hub.push_generator(write_tiny(tmp_path / "tiny-chem"))
    with pytest.raises(ValueError, match="load_network"):
        chemart.generate_network("alice/snap")
    with pytest.raises(ValueError, match="generate_network"):
        chemart.load_network("alice/tiny-chem")


def test_permissions_and_conflicts(live, tmp_path):
    live.user("alice")
    live.user("mallory")
    live.as_user("alice")
    hub.push_generator(write_tiny(tmp_path / "tiny-chem"))
    live.as_user("mallory")
    with pytest.raises(hub.HubAuthError):
        hub.push_generator(tmp_path / "tiny-chem", "alice/tiny-chem")
    live.as_user(None)
    with pytest.raises(hub.HubAuthError, match="chemart login"):
        chemart.generate_network("brusselator").push_to_hub("alice/snap")


def test_invalid_push_is_refused_before_upload(live, tmp_path):
    live.user("alice")
    folder = write_tiny(tmp_path / "tiny-chem", yaml_text=TINY_YAML.replace("fidelity: original", "fidelity: vibes"))
    with pytest.raises(hub.HubValidationError, match="fidelity"):
        hub.push_generator(folder)
    with pytest.raises(hub.RepoNotFoundError):
        hub.repo_info("alice/tiny-chem")                     # nothing was created


def test_cli_flow(live, tmp_path, capsys):
    token = live.user("carol", login=False)
    assert cli_main(["login", "--token", token]) == 0
    assert "as carol" in capsys.readouterr().out
    assert cli_main(["whoami"]) == 0 and "carol" in capsys.readouterr().out

    folder = tmp_path / "my-chem"
    assert cli_main(["new", str(folder)]) == 0
    assert cli_main(["check", str(folder)]) == 0, capsys.readouterr().out
    assert "ready to push" in capsys.readouterr().out
    assert cli_main(["push", str(folder)]) == 0
    assert f"{live.url}/carol/my-chem" in capsys.readouterr().out

    assert cli_main(["generate", "carol/my-chem"]) == 2
    assert "trust_remote_code" in capsys.readouterr().err
    assert cli_main(["generate", "carol/my-chem", "--trust-remote-code", "--format", "text"]) == 0
    assert "F + X -> 2 X" in capsys.readouterr().out

    net_file = tmp_path / "net.json"
    net_file.write_text(json.dumps(chemart.generate_network("brusselator").to_dict()))
    assert cli_main(["push-network", str(net_file), "carol/bruss", "--title", "Bruss", "--tag", "oscillator"]) == 0
    capsys.readouterr()
    assert cli_main(["search", "--tag", "oscillator"]) == 0
    assert "carol/bruss" in capsys.readouterr().out
    assert cli_main(["download", "carol/bruss"]) == 0
    assert Path(capsys.readouterr().out.strip(), "network.json").is_file()

    assert cli_main(["logout"]) == 0
    capsys.readouterr()
    assert cli_main(["whoami"]) == 3


def test_official_builtin_backed_repos(live, capsys):
    from chemart_hub.seed import seed

    live.user("curator")
    c = db.connect(live.settings.db_path)
    service.create_org(c, "chemart", "curator", allow_reserved=True)
    c.close()
    assert seed(only=["brusselator", "matrix-chemistry"]) == 0
    assert "2 of 2 official repos updated" in capsys.readouterr().out
    assert seed(only=["brusselator"]) == 0
    assert "unchanged" in capsys.readouterr().out

    live.as_user(None)
    official = chemart.generate_network("chemart/brusselator", seed=4)      # no trust needed
    assert official == chemart.generate_network("brusselator", seed=4)
    assert chemart.describe_chemistry("chemart/matrix-chemistry")["id"] == "matrix-chemistry"
    info = hub.repo_info("chemart/brusselator")
    assert info["builtin"] == "brusselator" and info["has_code"] is False
    assert "banzhaf-yamamoto" in info["tags"]

    # only the official org may point at built-ins
    live.user("mallory")
    from chemart_hub.seed import official_files
    with pytest.raises(hub.HubValidationError, match="reserved for the official"):
        hub.upload_files("mallory/brusselator", official_files("brusselator"), repo_type="generator")
