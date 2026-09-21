# Examples

| notebook | what it shows |
|---|---|
| [`simulating-chemistries.ipynb`](simulating-chemistries.ipynb) | Simulating a **given** network (the Brusselator: its record turned into rate equations and integrated, below and above the oscillation threshold) and a **generator** (the prime-number chemistry: running the chemistry *is* the simulation, and each run generates its own network) |

Open them with the optional `notebooks` dependency group:

```bash
uv sync --group notebooks
uv run --group notebooks jupyter lab examples/
```

The notebooks are written and executed by `build_notebooks.py`, so they ship
with real outputs; rerun it after changing one:

```bash
uv run --group notebooks python examples/build_notebooks.py
```
