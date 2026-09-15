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


def test_catalog_index_is_current():
    from pathlib import Path

    from chemart.catalog import ROOT, render_index

    generated = render_index(load())
    on_disk = Path(ROOT / "docs" / "CATALOG.md").read_text()
    assert generated == on_disk, "run `python -m chemart.catalog index`"
