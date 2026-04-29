# Upstream examples audit

This repository currently matches the upstream-style EvolvePro layout plus the local Nix/uv shell. The CLI migration regression baseline must therefore be derived from the checked-in upstream examples, not from downstream BenM/TphR experiment data.

## Source baseline

- Frozen reference tag: `baseline-v0` at `e54fd1e`.
- Public workflow description: `README.md`.
- Tutorial workflow: `colab/EvolvePro.ipynb`.
- Example data shipped in the repo:
  - `colab/rounds_data/kelsic_Round1.xlsx` through `kelsic_Round5.xlsx`.
  - `data/dms/activity/*`.
  - `data/dms/wt_fasta/*`.
  - `data/exp/Evolution_Summary_Data.xlsx`.
- Existing script entrypoints:
  - `scripts/process/{exp_process.py,dms_process.py}`.
  - `scripts/dms/dms_main.py`.
  - `scripts/exp/{t7_pol.py,mlv.py,bxb1.py}`.
  - `scripts/plot/{dms.py,exp.py}`.

## README workflow map

| README step | Existing entrypoints | Checked-in data | Regression plan |
| --- | --- | --- | --- |
| Process | `evolvepro.src.process.generate_wt`, `generate_single_aa_mutants`, `generate_n_mutant_combinations`, `suggest_initial_mutants` | Kelsic WT sequence from `colab/EvolvePro.ipynb`; DMS WT FASTA files | Golden tests for FASTA generation and seeded initial mutant suggestions. |
| PLM | `evolvepro/plm/*/extract.py`; `scripts/plm/*.sh` | No small checked-in embeddings; extractors require model downloads/GPU | Static/import smoke only. Do not run model downloads in regression tests. |
| Run EVOLVEpro: DMS | `evolvepro.src.evolve.directed_evolution_simulation`, `grid_search`; `scripts/dms/dms_main.py` | `data/dms/activity/*`, `data/dms/wt_fasta/*` | Golden tests for a small deterministic simulation and smoke/parser tests for `dms_main.py`; lightweight DMS data processing smoke where runtime is acceptable. |
| Run EVOLVEpro: Experimental | `evolvepro.src.evolve.evolve_experimental`, `evolve_experimental_multi`; `evolvepro.src.data.load_experimental_data`, `create_iteration_dataframes` | `colab/rounds_data/kelsic_Round*.xlsx`; Kelsic WT sequence from notebook | Golden tests for Kelsic round loading, iteration/label dataframe construction, and a small experimental round prediction path using synthetic embeddings aligned to Kelsic variants. |
| Plot | `evolvepro.src.plot.read_exp_data`, `process_dataframe`, `filter_dataframe`, plotting helpers; `scripts/plot/*.py` | Colab Kelsic rounds; generated DMS/experiment summaries | Golden tests for non-plotting dataframe transforms; plot rendering is smoke-only or deferred. |

## Canonical tutorial: Colab Kelsic workflow

`colab/EvolvePro.ipynb` documents the tutorial sequence:

1. Generate a WT FASTA for Kelsic using the sequence:
   `MAKEDNIEMQGTVLETLPNTMFRVELENGHVVTAHISGKMRKNYIRILTGDKVTVELTPYDLSKGRIVFRSR`.
2. Generate all single amino-acid mutants with `generate_single_aa_mutants`.
3. Suggest 12 initial mutants with `suggest_initial_mutants(..., random_seed=42)`.
4. Extract ESM embeddings with `evolvepro/plm/esm/extract.py`.
5. Run `evolve_experimental` for `Round1` through `Round5` using `colab/rounds_data/kelsic_Round*.xlsx`.
6. Read and plot experimental data with `read_exp_data` and `plot_variants_by_iteration`.

Regression tests should execute steps 1, 2, 3, 5 data-loading pieces, and 6 dataframe-loading pieces. Step 4 is excluded from normal tests because it downloads/runs a PLM.

## Reachable code for migration

Do not remove these paths before regression coverage is in place:

- `evolvepro/src/model.py`
  - `first_round` strategies: `random`, `diverse_medoids`, `explicit_variants`.
  - `top_layer` regression types used by examples or regression tests: `ridge`, `randomforest`, `knn`.
- `evolvepro/src/evolve.py`
  - `directed_evolution_simulation`.
  - `grid_search`.
  - `evolve_experimental`.
  - `evolve_experimental_multi`.
- `evolvepro/src/data.py`
  - `load_dms_data` and `process_pt_embeddings`.
  - `load_experimental_embeddings`.
  - `load_experimental_data`.
  - `create_iteration_dataframes`.
- `evolvepro/src/process.py`
  - DMS processing helpers used by `scripts/process/dms_process.py`.
  - experimental FASTA helpers used by `scripts/process/exp_process.py` and the Colab notebook.
- `evolvepro/src/plot.py`
  - `read_exp_data` and dataframe aggregation/filter helpers.
- `scripts/dms/dms_main.py` parser contract.
- `scripts/plm/*` and `evolvepro/plm/*` as static entrypoints until the CLI has equivalent wrappers.

## Stage 0 regression matrix

| Test module | Coverage | Expected data policy |
| --- | --- | --- |
| `tests/test_process_examples.py` | Kelsic WT FASTA, Kelsic single-AA mutant subset, seeded suggestions | `pytest-regressions` golden generated from `baseline-v0`. |
| `tests/test_colab_kelsic_data.py` | Kelsic Excel loading, iteration/label construction, plot reader | Golden summaries/dataframes generated from checked-in Colab data. |
| `tests/test_model_regression.py` | `first_round` and `top_layer` behavior on deterministic synthetic data | Golden selected variants, metrics, and predictions. Do not recompute expected values with duplicate sklearn logic. |
| `tests/test_dms_simulation.py` | Small deterministic `directed_evolution_simulation` run | Golden metrics dataframe. |
| `tests/test_dms_data_examples.py` | Lightweight DMS checked-in data processing/load path | Golden schema/summary. Skip expensive full benchmark runs. |
| `tests/test_scripts_smoke.py` | `scripts/dms/dms_main.py` parser and static script contracts | No script execution with absolute lab paths. |
| `tests/test_baseline_file_outputs.py` | Existing file outputs from `evolve_experimental`, `generate_wt`, and `generate_single_aa_mutants` | Byte-for-byte file regression for CSV/FASTA artifacts used as CLI migration targets. |

## Locked Stage 0 baseline

The Stage 0 baseline currently consists of 22 tests and checked-in `pytest-regressions` golden files. It intentionally covers both function-level behavior and file-level outputs before any CLI refactor.

- Process workflow:
  - `tests/test_process_examples.py` locks FASTA metadata, mutant IDs, and seeded suggestion stdout.
  - `tests/test_baseline_file_outputs.py` locks generated Kelsic WT and mutant FASTA files byte-for-byte.
- Colab Kelsic experimental workflow:
  - `tests/test_colab_kelsic_data.py` locks Kelsic Excel loading, iteration/label construction, `read_exp_data`, and an `evolve_experimental` Round1 run with deterministic synthetic embeddings.
  - `tests/test_baseline_file_outputs.py` locks the existing `evolve_experimental` output files: `iteration.csv`, `this_round_variants.csv`, `df_test.csv`, and `df_sorted_all.csv`.
- Model workflow:
  - `tests/test_model_regression.py` locks `first_round` (`random`, `diverse_medoids`) and `top_layer` (`ridge`, `randomforest`, `knn`) on deterministic synthetic data.
- DMS workflow:
  - `tests/test_dms_simulation.py` locks a small `directed_evolution_simulation` metrics dataframe.
  - `tests/test_dms_data_examples.py` locks processing of `data/dms/activity/sarscov2_preprocessed.csv` with `data/dms/wt_fasta/cov2_S_WT.fasta`.
- Script contracts:
  - `tests/test_scripts_smoke.py` locks the `scripts/dms/dms_main.py` parser and statically records the experimental script entrypoint calls.

Run before and after every migration commit:

```bash
nix develop -c uvx ruff check .
nix develop -c uv run pytest tests/ -v
```

Current expected result: `22 passed`. The warning count is not part of the behavioral contract; warnings are handled in separate cleanup commits.

## Stage 0b deletion proposal

Deletion must be conservative: remove only paths that are neither documented by the upstream examples nor covered by the locked Stage 0 baseline. Anything used by README, the Colab tutorial, or scripts remains until the CLI has equivalent coverage.

### Safe to remove now

No Python package code is safe to remove immediately. The baseline intentionally keeps the upstream public surface broad until the CLI wrappers exist.

The only candidates for immediate deletion are non-runtime presentation artifacts if repository size cleanup is desired:

- `images/*.gif`, `images/*.png` after README is rewritten and no longer references them.

Do not delete these in Stage 0b if README still references them.

### Safe to remove only after CLI replacement exists

These script entrypoints are legacy wrappers, but they document the current public workflow. Delete them only after Typer commands provide equivalent behavior and tests compare against the locked baseline outputs.

- `scripts/process/exp_process.py`
- `scripts/process/dms_process.py`
- `scripts/dms/dms_main.py`
- `scripts/exp/{t7_pol.py,mlv.py,bxb1.py}`
- `scripts/plot/{dms.py,exp.py}`
- `scripts/plm/*.sh`

Suggested deletion stage: after CLI commands for `process`, `fit-round`, `dms-grid`, `plot`, and `embed` exist and `tests/test_scripts_smoke.py` has been replaced by CLI tests.

### Defer deletion

These paths are upstream example data or tutorial references and should remain until the final library-only cleanup stage.

- `colab/EvolvePro.ipynb`
- `colab/rounds_data/kelsic_Round*.xlsx`
- `data/dms/activity/*`
- `data/dms/wt_fasta/*`
- `data/exp/Evolution_Summary_Data.xlsx`

Reason: Stage 0 uses these as canonical public fixtures. Move/delete only after the CLI package has fixture copies or a documented external experiment-data strategy.

### Code paths to keep through Stage 3

Keep these through dump/config/CLI wrapping because they are either covered by tests or are needed by upstream workflows:

- `evolvepro/src/model.py`
  - `first_round`: `random`, `diverse_medoids`, `explicit_variants`.
  - `top_layer`: at least `ridge`, `randomforest`, `knn`.
  - Other regression branches (`lasso`, `elasticnet`, `linear`, `neuralnet`, `gradientboosting`, `gp`) are not Stage-0-covered; remove only after confirming no documented script uses them, or add tests first.
- `evolvepro/src/evolve.py`
  - `directed_evolution_simulation`, `grid_search`, `evolve_experimental`, `evolve_experimental_multi`.
- `evolvepro/src/data.py`
  - CSV/Excel loaders, DMS loader, `.pt` loader support. Move `torch` to lazy import before deciding whether `.pt` support is optional.
- `evolvepro/src/process.py`
  - FASTA generation, mutant suggestion, DMS processing helpers.
- `evolvepro/src/plot.py`
  - `read_exp_data` and dataframe helpers.
- `evolvepro/plm/*`
  - Keep as optional PLM/embed backend until `evolvepro embed` exists.

### Recommended Stage 0b action

Do not delete code yet. Instead, make Stage 0b a reachable-code audit commit that records:

1. which CLI commands will replace each legacy script,
2. which uncovered regression branches are candidates for later removal,
3. which data/notebook assets are test fixtures vs. documentation-only assets,
4. where `torch` must become a lazy optional dependency.

After that audit, start Stage 1 dump infrastructure on top of the now-locked 22-test baseline.

## Exclusions from regression tests

- Network access.
- GPU execution.
- PLM model downloads.
- Full notebook execution.
- Full DMS benchmark/grid search sweeps.
- Experimental scripts with absolute `/orcd/...` paths; test them statically unless they are refactored to safe functions.

## Commit boundaries

Stage 0 test commits may change only:

- `tests/`.
- `pyproject.toml`.
- `uv.lock`.
- `.gitignore` if needed for generated regression artifacts.

They must not change library code, scripts, checked-in data, notebook files, or Nix formatting.
