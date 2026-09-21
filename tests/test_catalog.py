"""The catalog is data we rely on, so keep it honest."""
from chemart.catalog import PROVIDES, load, validate


def test_catalog_is_valid():
    assert validate(load()) == []


def test_ids_are_unique_and_slugs():
    ids = [c.id for c in load()]
    assert len(ids) == len(set(ids))
    assert all(i == i.lower() and " " not in i for i in ids)


def test_every_generator_declares_capabilities():
    for c in load():
        if c.implemented or (c.kind == "generator" and c.generator_ready in {"yes", "partial"}):
            assert set(c.provides) <= PROVIDES
            assert c.provides, c.id


def test_every_entry_has_documented_params_or_says_why():
    """A parameterless chemistry is suspicious; only analyses may have none."""
    for c in load():
        if c.kind == "generator" and (c.implemented or c.generator_ready == "yes"):
            assert c.params, f"{c.id} is ready or implemented but has no parameters"


def test_archive_groups_are_known():
    from dataclasses import replace

    entries = load()
    first = entries[0]
    assert any("unknown archive group" in p for p in validate([replace(first, archived="attic")] + entries[1:]))
    assert validate([replace(first, archived="artificial-life")] + entries[1:]) == []


def test_archived_entries_leave_the_listings_but_still_run():
    import chemart
    from chemart.catalog import active

    entries = load()
    archived = {c.id for c in entries if c.archived}
    assert archived and len(active(entries)) == len(entries) - len(archived)
    assert not archived & {row["id"] for row in chemart.list_chemistries()}
    full = chemart.list_chemistries(include_archived=True)
    assert {row["id"] for row in full if row["archived"]} == archived
    tools = {t["name"]: t for t in chemart.tool_definitions()}
    assert not archived & set(tools["generate_network"]["input_schema"]["properties"]["chemistry"]["enum"])
    assert chemart.describe_chemistry("tierra")["archived"] == "artificial-life"
    assert chemart.generate_network("lotka-volterra", seed=0).reactions   # archived, still runnable


def test_catalog_index_is_current():
    from pathlib import Path

    from chemart.catalog import ROOT, render_index

    generated = render_index(load())
    on_disk = Path(ROOT / "docs" / "CATALOG.md").read_text()
    assert generated == on_disk, "run `python -m chemart.catalog index`"


# --- single entries (the hub validates uploads with the same functions) ----
import copy

import pytest
import yaml

from chemart.catalog import CATALOG_DIR, entry_problems, parse_entry


def _raw(cid="brusselator"):
    return copy.deepcopy(yaml.safe_load((CATALOG_DIR / f"{cid}.yaml").read_text())["chemistries"][0])


def test_parse_entry_round_trips_a_catalog_file():
    c = parse_entry(_raw(), "brusselator.yaml")
    assert c.id == "brusselator" and c.params[0].name == "a"
    assert entry_problems(c) == []


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda r: r.update(colour="teal"), "unknown field.*colour"),
        (lambda r: r.pop("family"), "missing required field.*family"),
        (lambda r: r["params"][0].update(unit="mM"), r"params\[0\].*unknown field.*unit"),
        (lambda r: r["params"][0].pop("role"), r"params\[0\].*missing.*role"),
        (lambda r: r.update(params="N=4"), "params must be a list"),
    ],
)
def test_parse_entry_errors_are_readable(mutate, message):
    raw = _raw()
    mutate(raw)
    with pytest.raises(ValueError, match=message):
        parse_entry(raw)


def test_hub_rules_differ_from_the_builtin_catalog():
    raw = _raw()
    raw.pop("book")
    raw["fidelity"] = "original"
    c = parse_entry(raw)
    assert "missing book section" in entry_problems(c)
    assert any("fidelity" in p for p in entry_problems(c))
    assert entry_problems(c, hub=True) == []


def test_reserved_parameter_names():
    raw = _raw()
    raw["params"][0]["name"] = "seed"
    problems = entry_problems(parse_entry(raw), hub=True)
    assert any("'seed' is reserved" in p for p in problems)


def test_load_refuses_an_empty_catalog(tmp_path):
    with pytest.raises(RuntimeError, match="no catalog entries"):
        load(tmp_path)
