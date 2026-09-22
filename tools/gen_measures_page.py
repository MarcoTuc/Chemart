"""Regenerate the measure tables of docs/guide/measures.md from the registry.

    uv run python tools/gen_measures_page.py           # rewrite the tables
    uv run python tools/gen_measures_page.py --check   # exit 1 if they are stale

Each section of the page has a block between `<!-- measures X -->` and
`<!-- /measures -->`; the block is rebuilt from `chemart.measures.REGISTRY`,
so the page cannot drift from the code.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from chemart.measures import REGISTRY  # noqa: E402

PAGE = ROOT / "docs" / "guide" / "measures.md"

def table(section: str) -> str:
    rows = ["| measure | meaning | needs | input | cost |", "|---|---|---|---|---|"]
    for m in REGISTRY.values():
        if m.section == section:
            limit = f", up to {m.limit} nodes" if m.limit else ""
            rows.append(f"| `{m.name}` | {m.doc} | {' '.join(sorted(m.needs))} | {m.input} | {m.cost}{limit} |")
    return "\n".join(rows)


BLOCK = re.compile(r"<!-- measures ([A-J]) -->\n.*?<!-- /measures -->", re.S)


def render(text: str) -> str:
    return BLOCK.sub(lambda m: f"<!-- measures {m.group(1)} -->\n{table(m.group(1))}\n<!-- /measures -->", text)


def main(argv: list[str]) -> int:
    text = PAGE.read_text()
    new = render(text)
    if "--check" in argv:
        if new != text:
            print("docs/guide/measures.md is stale: run tools/gen_measures_page.py")
            return 1
        return 0
    PAGE.write_text(new)
    print(f"wrote {PAGE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
