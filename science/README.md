# science/

This is research, not platform. The folder holds experiments on the neutral networks of molecules, meaning how much a one-symbol or one-atom change alters what a molecule does. It covers two artificial chemistries (AlChemy, Stringmol) and two real reaction networks (MCM, Rhea). The motivation and the literature behind it are in `../research/measures-literature.md`. The plan is in [`PLAN.md`](PLAN.md).

## How it relates to chemart

- **It uses chemart but never changes it.** Code here imports `chemart`, for example Stringmol's `Machine.react` and toychem's RDKit helpers, but no file under `chemart/` is modified for this work. If an experiment needs something chemart lacks, it is built here.
- **External engines are vendored here at a pinned commit.** Each one has a `VENDOR.md` recording its source, commit and licence. The AlChemy engine is the Mathis group's Rust reimplementation, `AgentElement/functional-supercollider`. Both it and the Stringmol C++ fallback are GPL-3.0, so code derived from them is GPL-3.0 too.

## Layout

This is the planned layout. Folders are created when their phase starts.

```
science/
  README.md  PLAN.md  PREREGISTRATION.md (per phase)
  common/      shared Python: phenotypes, locality ratio, statistics, IO
  alchemy/     engine/ (vendored Rust, pinned), mut/ (our crate), *.py
  stringmol/   *.py (chemart oracle); C++ fallback if needed
  realchem/    mcm/, rhea/, mcgillen/ loaders and analyses
  notes/       results notes, one per phase
  results/     small summary tables (tracked)
  data/        downloaded datasets (git-ignored)
  runs/        raw run outputs (git-ignored)
```

## Conventions

- **Python.** Use uv and the repository's `./.venv`. Run everything with `uv run` and never with pip or a system Python. Dependencies needed only here go in a `science` dependency group in the root `pyproject.toml`, added in Phase 1.
- **Rust.** Build with cargo. Build outputs (`target/`) are git-ignored.
- **Every run writes a manifest** recording the engine commit, config, seed and command, so that any result can be regenerated.
- **Predictions are registered before confirmatory runs.** They go in `PREREGISTRATION.md`, and anything decided afterwards is labelled exploratory.
- **Raw data and large outputs stay out of git.** Downloaded papers live in `../research/library/`, which is covered by the repository's `*.pdf` rule.
