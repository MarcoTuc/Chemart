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


def test_every_implemented_entry_says_what_it_is():
    from chemart.catalog import TYPES

    missing = [c.id for c in load() if c.implemented and c.type not in TYPES]
    assert not missing, f"set type: given | generator | gas for {missing}"


def test_the_type_table_matches_the_catalog():
    """catalog/TYPES.md gives the reasoning; it must agree with the YAML."""
    import re

    from chemart.catalog import ROOT

    rows = {}
    for line in (ROOT / "catalog" / "TYPES.md").read_text().splitlines():
        m = re.match(r"\| `([a-z0-9-]+)`[^|]* \| (given|generator|gas)\b", line)
        if m:
            rows[m.group(1)] = m.group(2)
    assert rows == {c.id: c.type for c in load()}


def test_describe_reports_the_type():
    import chemart

    assert chemart.describe_chemistry("brusselator")["type"] == "given"
    assert chemart.describe_chemistry("kauffman-autocatalytic-sets")["type"] == "generator"
    assert chemart.describe_chemistry("alchemy")["type"] == "gas"
    assert {row["type"] for row in chemart.list_chemistries()} == {"given", "generator", "gas"}


def test_a_given_network_needs_no_closure():
    """A given network is written down whole; only a computed one can be cut short."""
    import chemart

    for c in load():
        if c.type == "given" and c.implemented:
            assert chemart.generate_network(c.id, seed=1).status != "truncated", c.id


def test_face_and_every_are_checked():
    from dataclasses import replace

    from chemart.catalog import entry_problems

    c = next(c for c in load() if c.id == "brusselator")
    p = c.params[0]
    assert any("unknown face" in m for m in entry_problems(replace(c, params=[replace(p, face="closure")])))
    assert any("reserved" in m for m in entry_problems(replace(c, params=[replace(p, name="every")])))
    assert any("unknown type" in m for m in entry_problems(replace(c, type="given-ish")))


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
