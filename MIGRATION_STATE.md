# EvolvePro CLI migration state

This file records the current migration state so work can resume after context compaction.

## Branch and baseline

- Active branch: `cli-migration`.
- Frozen baseline tag: `baseline-v0` -> `e54fd1e` (`init flake`).
- Do not move or delete `baseline-v0`.
- The original BenM/TphR migration plan is obsolete for this repository state. This repo is currently the upstream-style fork with Colab/DMS examples.

## Ground truth for regression tests

Regression tests are based on checked-in upstream examples:

- `README.md` workflow description.
- `colab/EvolvePro.ipynb` Kelsic tutorial.
- `colab/rounds_data/kelsic_Round1.xlsx` through `kelsic_Round5.xlsx`.
- `data/dms/activity/*` and `data/dms/wt_fasta/*`.
- `scripts/{process,dms,exp,plot}/README.md` and script entrypoints.

See:

- `tests/UPSTREAM_EXAMPLES.md` for upstream workflow audit and deletion proposal.
- `tests/REACHABLE.md` for reachable-code audit and CLI replacement gates.
- `docs/WORKSPACE_CONTRACT.md` for the planned file-based workspace, round, lock, and manifest contract.

## Completed commits on `cli-migration`

```text
79ab014 docs: audit upstream examples for regression baseline
03b74c7 build: declare upstream runtime and test dependencies
56342ee test: lock upstream example behavior
9350fef chore: apply formatter and ruff cleanup
b1a141c test: capture baseline file outputs
a40770c docs: record baseline coverage and deletion scope
15af389 docs: audit reachable code before migration
14fb720 refactor: lazy import torch for pt embeddings
214160c feat: write reproducible round dump artifacts
```

## Current validation command

Run after every migration commit:

```bash
nix develop -c uvx ruff check .
nix develop -c uv run pytest tests/ -v
```

Current expected result:

```text
29 passed, 12 warnings
```

The warnings are known and not part of the behavioral contract:

- `pandas` `FutureWarning` from `evolvepro/src/model.py` `pd.concat([df_train, df_test])`.
- `sklearn_extra` `distutils` `DeprecationWarning`.

Do not fix these inside unrelated migration commits.

## Stage 0 baseline coverage

The Stage 0 tests lock both function behavior and file outputs.

- Process workflow:
  - `generate_wt`.
  - `generate_single_aa_mutants`.
  - `suggest_initial_mutants`.
  - Byte-for-byte Kelsic FASTA outputs.
- Colab Kelsic experimental workflow:
  - `load_experimental_data`.
  - `create_iteration_dataframes`.
  - `read_exp_data`.
  - `evolve_experimental` Round1 with deterministic synthetic embeddings.
  - Byte-for-byte existing CSV outputs:
    - `iteration.csv`
    - `this_round_variants.csv`
    - `df_test.csv`
    - `df_sorted_all.csv`
- Model workflow:
  - `first_round(random)`.
  - `first_round(diverse_medoids)`.
  - `top_layer(ridge)`.
  - `top_layer(randomforest)`.
  - `top_layer(knn)`.
- DMS workflow:
  - Small `directed_evolution_simulation` golden.
  - `process_dataset` on `sarscov2_preprocessed.csv` with `cov2_S_WT.fasta`.
- Script contracts:
  - `scripts/dms/dms_main.py` parser.
  - static `scripts/exp/*.py` entrypoint calls.

## Dependency and Nix notes

- Runtime dependencies are declared in `pyproject.toml` and locked in `uv.lock`.
- `joblib` and `pyarrow` are runtime dependencies because dump artifacts write `model.pkl` and `predictions.parquet`.
- `torch` is intentionally not a required runtime dependency for core CSV/Excel workflows.
- `evolvepro/src/data.py` now lazy-imports `torch` only for `.pt` embedding files.
- Keep CLI/subprocess implementations from overwriting `LD_LIBRARY_PATH`; Nix devShell sets it for native libraries and host CUDA driver paths.
- PLM/GPU/model-download paths should remain optional and smoke-tested only in normal CI.

## Dump artifact implementation

Implemented in `214160c`.

New files:

- `evolvepro/__init__.py`
- `evolvepro/io/__init__.py`
- `evolvepro/io/dump.py`
- `evolvepro/io/hashing.py`
- `tests/test_dump_artifacts.py`

`top_layer(..., experimental=True)` now returns:

```python
this_round_variants, df_test, df_sorted_all, metrics, model
```

Updated callers:

- `evolve_experimental`
- `evolve_experimental_multi`

Each experimental round output directory now contains legacy CSVs plus:

- `model.pkl`
- `predictions.parquet`
- `metrics.json`
- `inputs.json`
- `version.txt`

Existing CSV outputs remain byte-identical to the Stage 0 baseline tests.

## Deletion policy

No code deletion has been performed.

Do not delete legacy scripts or example data until equivalent CLI commands exist and tests compare CLI outputs against the locked baseline.

Defer deletion of:

- `scripts/process/*`
- `scripts/dms/*`
- `scripts/exp/*`
- `scripts/plot/*`
- `scripts/plm/*`
- `colab/*`
- `data/*`
- `images/*` while README references them

## Next stage: workspace/config/manifest abstraction

Start here after context compaction. Read `docs/WORKSPACE_CONTRACT.md` first.

Goal: introduce file-based TOML config, round lock files, and manifests for the current upstream-style experimental workflow without removing existing function/script behavior.

Key decision: do not use SQLite as source of truth. Round state is authoritative in `round.toml`, `round.lock.json`, and `manifest.json`. A future SQLite index may only be a rebuildable `.evolvepro/` cache.

Suggested files:

- `evolvepro/io/schemas.py`
- `evolvepro/io/config.py`
- `evolvepro/io/manifest.py`
- `evolvepro/io/canonical.py`
- `tests/fixtures/kelsic_round1.toml`
- `tests/test_config.py`
- `tests/test_manifest.py`

Suggested config model fields for the current Kelsic experimental path:

```python
@dataclass
class ExperimentalRoundConfig:
    protein_name: str
    round_name: str
    embeddings_base_path: Path
    embeddings_file_name: str
    round_base_path: Path
    round_file_names: list[str]
    wt_fasta_path: Path
    rename_WT: bool = False
    number_of_variants: int = 12
    output_dir: Path | None = None
```

Also plan for multi-round config later:

```python
@dataclass
class ExperimentalMultiRoundConfig:
    protein_name: str
    round_name: str
    embeddings_base_path: Path
    embeddings_file_names: list[str]
    round_base_path: Path
    round_file_names_single: list[str]
    round_file_names_multi: list[str]
    wt_fasta_path: Path
    rename_WT: bool = False
    number_of_variants: int = 12
    output_dir: Path | None = None
```

Tests to add first:

1. TOML loader validates required fields and path/string/list types.
2. Invalid TOML raises a clear exception.
3. Canonical JSON/TOML writes are deterministic and stable.
4. Raw file SHA256 and table semantic hashes are reproducible.
5. Round lock generation captures resolved config, input hashes, previous manifest hash, and model parameters.
6. Manifest validation stops on input or previous-manifest hash mismatch.
7. TOML-driven `evolve_experimental` produces the same legacy CSV outputs and dump artifacts as direct function invocation using the existing Kelsic synthetic embedding fixture.

Important: there is no `data/<enzyme>/config.py` legacy adapter in this repo. The adapter should target current function arguments and later the Typer CLI, not the obsolete BenM/TphR layout.
