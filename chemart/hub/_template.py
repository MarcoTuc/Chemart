"""The skeleton `chemart new` writes: a working generator repo to edit."""

from __future__ import annotations

from pathlib import Path

CHEMART_YAML = """\
# The chemistry's catalog entry: the same schema as Chemart's built-in
# catalog (catalog/SCHEMA.md), plus a `hub:` block. Edit freely, then
#   chemart check .      runs the contract checks locally
#   chemart push .       publishes it
hub:
  repo_type: generator
  license: MIT
  tags: [autocatalysis]
  # requires: [networkx]          # extra modules your code imports (checked, never installed)
chemistries:
- id: {id}
  name: {name}
  intuition: >
    Say in plain words what the molecules stand for, what a reaction does to
    them, and why this chemistry is interesting. This is the first thing
    people read on the hub.
  origin: "you, {year}"
  family: core             # core, rewriting, automata, bio-inspired, origin-of-life, network, spatial, ...
  kind: generator
  constructive: false      # true if new species appear while it runs
  S: {{repr: "F (food), X (replicator)", definition: explicit}}
  R: {{definition: explicit, arity: [1, 2], scheme: "F + X -> 2X ; X -> 0"}}
  A: {{reactor: [ode, ssa], dilution: "none"}}
  params:
    - {{name: k, type: float, default: 1.0, min: 0, role: kinetic, meaning: "replication rate constant"}}
    - {{name: d, type: float, default: 0.1, min: 0, role: kinetic, meaning: "decay rate constant"}}
  provides: [topology, stoichiometry, catalysts, rate-constants, initial-state]
  fidelity: original       # or book / book+decisions / reconstructed, for published chemistries
"""

GENERATOR_PY = '''\
"""{name}.

generate(p, rng) receives the validated parameters as attributes of `p` and a
numpy Generator `rng` (draw all randomness from it), and returns a Network.
Helpers from Chemart (chemart.expand, chemart.soup, chemart.helpers) are
available; other files of this repo import relatively: `from .rules import x`.
"""

from chemart.helpers.explicit import network


def generate(p, rng):
    return network(
        [
            ("F + X -> 2 X", p.k),
            ("X -> ", p.d),
        ],
        species=["F", "X"],
        initial_state={{"F": 1.0, "X": 0.01}},
    )
'''

README_MD = """\
# {name}

What this chemistry is, where it comes from, and what to look for when you
run it.

```python
import chemart

net = chemart.generate_network("<you>/{id}", seed=0, trust_remote_code=True)
print(net.to_text())
```
"""


def write(folder: Path, chemistry_id: str, name: str | None = None) -> list[Path]:
    """Create the skeleton in `folder`; refuse to overwrite existing files."""
    from datetime import date

    from chemart.hub._ids import check_name

    check_name(chemistry_id)
    name = name or chemistry_id.replace("-", " ").capitalize()
    values = {"id": chemistry_id, "name": name, "year": date.today().year}
    files = {
        "chemart.yaml": CHEMART_YAML.format(**values),
        "generator.py": GENERATOR_PY.format(**values),
        "README.md": README_MD.format(**values),
    }
    folder = Path(folder)
    clashes = [p for p in files if (folder / p).exists()]
    if clashes:
        raise FileExistsError(f"{folder} already has {', '.join(clashes)}")
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for path, text in files.items():
        (folder / path).write_text(text)
        written.append(folder / path)
    return written
