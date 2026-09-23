# Examples

| notebook | what it shows |
|---|---|
| [`simulating-chemistries.ipynb`](simulating-chemistries.ipynb) | The three types in practice: a **given** network simulated with `chemart.simulate` (the Brusselator, by rate equations below and above its oscillation threshold, and by stochastic paths at two volumes); two **gases** evolved with `chemart.evolve` (the prime-number soup turning prime, and AlChemy with measures followed in evolutionary time by `measures.over`); and networks of all three types compared with `chemart.measure` |

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
