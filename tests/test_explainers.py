"""The hand-written prose of the catalog pages: catalog/explainers/<id>.md.

Every explainer must belong to a catalog entry and have the sections the page
generator expects, so a malformed one fails here rather than in the docs build.
"""

import importlib.util
from pathlib import Path

import pytest

from chemart.catalog import load

REPO = Path(__file__).resolve().parent.parent
EXPLAINERS = REPO / "catalog" / "explainers"

_spec = importlib.util.spec_from_file_location("gen_catalog_pages", REPO / "tools" / "gen_catalog_pages.py")
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)

IDS = {c.id for c in load()}
FILES = sorted(p for p in EXPLAINERS.glob("*.md") if p.name != "README.md")


def test_every_explainer_belongs_to_an_entry():
    assert {p.stem for p in FILES} <= IDS


def test_every_catalog_chemistry_has_an_explainer():
    """Archived entries may go without one; the chemistry catalog may not."""
    from chemart.catalog import active

    missing = sorted({c.id for c in active(load())} - {p.stem for p in FILES})
    assert not missing, f"write catalog/explainers/<id>.md for: {', '.join(missing)}"


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.stem)
def test_explainer_has_the_page_sections(path):
    sections = gen.load_explainer(path.stem)
    required = [s for s, needed in gen.SECTIONS.items() if needed]
    assert all(sections[s] for s in required)


GOOD = "## Introduction\nA.\n## How it works\nB.\n## Using it\nC.\n## Results\nD.\n"


@pytest.mark.parametrize("text, message", [
    ("# Title\n" + GOOD, "no '# ' title"),
    ("Preamble.\n" + GOOD, "text before the first"),
    (GOOD + "## Notes\nE.\n", "unknown section"),
    (GOOD + "## Results\nE.\n", "appears twice"),
    (GOOD.replace("## Results\nD.\n", ""), "missing section"),
    ("## How it works\nB.\n## Introduction\nA.\n## Using it\nC.\n## Results\nD.\n", "order"),
])
def test_malformed_explainers_are_refused(tmp_path, text, message):
    (tmp_path / "x.md").write_text(text)
    with pytest.raises(ValueError, match=message):
        gen.load_explainer("x", tmp_path)


def test_headings_inside_code_are_not_sections(tmp_path):
    text = GOOD.replace("C.\n", "```python\n## not a section\n# nor a title\n```\n")
    (tmp_path / "x.md").write_text(text)
    assert "## not a section" in gen.load_explainer("x", tmp_path)["Using it"]
