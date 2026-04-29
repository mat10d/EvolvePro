# Reachable code audit for CLI migration

This audit is the Stage 0b replacement for the original deletion-first plan. The current upstream fork exposes a broad example surface through README, Colab, checked-in data, and scripts. Deletion is therefore deferred until equivalent CLI commands exist and are tested against the Stage 0 golden baseline.

## Baseline command

Run before and after every migration commit:

```bash
nix develop -c uvx ruff check .
nix develop -c uv run pytest tests/ -v
```

Current baseline: `22 passed`.

## Legacy entrypoint to CLI mapping

| Legacy entrypoint | Current purpose | Future CLI command | Deletion gate |
| --- | --- | --- | --- |
| `scripts/process/exp_process.py` | Generate WT FASTA, single-AA mutants, and n-mutant combinations for experimental workflows. | `evolvepro init` / `evolvepro mutants` | Delete after CLI generates byte-identical Kelsic FASTA outputs covered by `tests/test_baseline_file_outputs.py`. |
| `scripts/process/dms_process.py` | Process checked-in DMS source files and generate labels/FASTA. | `evolvepro process-dms` | Delete after CLI reproduces a checked-in DMS processing golden, at least the SARS-CoV-2 smoke path. |
| `scripts/dms/dms_main.py` | Argparse wrapper over `evolvepro.src.evolve.grid_search`. | `evolvepro dms-grid` | Delete after Typer parser tests replace `tests/test_scripts_smoke.py::test_dms_main_parser_contract`. |
| `scripts/exp/{t7_pol.py,mlv.py,bxb1.py}` | Absolute-path examples for experimental evolution campaigns. | `evolvepro fit-round --config ...` and `evolvepro fit-multi-round --config ...` | Delete after CLI config examples cover single- and multi-mutant experimental paths. |
| `scripts/plot/{dms.py,exp.py}` | Absolute-path plot examples. | `evolvepro plot dms` / `evolvepro plot exp` | Delete after dataframe transform tests and CLI smoke tests cover plot inputs without requiring display/GPU. |
| `scripts/plm/*.sh` | Shell wrappers around PLM embedding extractors. | `evolvepro embed` | Delete after PLM/embed wrappers exist; normal CI should remain CPU/no-download smoke only. |

## Library API reachability

### `evolvepro/src/model.py`

Covered by Stage 0:

- `first_round(first_round_strategy="random")`.
- `first_round(first_round_strategy="diverse_medoids")`.
- `top_layer(regression_type="ridge")`.
- `top_layer(regression_type="randomforest")`.
- `top_layer(regression_type="knn")`.
- `top_layer(..., experimental=True)` indirectly through `evolve_experimental`.

Reachable but not directly covered:

- `first_round(first_round_strategy="explicit_variants")`.
  - Keep until CLI/config has explicit initial-candidate support or add a direct regression test first.

Uncovered regression branches:

- `lasso`.
- `elasticnet`.
- `linear`.
- `neuralnet`.
- `gradientboosting`.
- `gp`.

Deletion decision:

- Do not remove these branches before Stage 3 because `scripts/dms/esm2_15B_full.sh` documents a full DMS sweep using all regression types.
- Later options:
  1. keep them and add focused tests,
  2. move them behind an unsupported/legacy flag,
  3. remove them only after CLI scope explicitly drops the full upstream DMS sweep.

### `evolvepro/src/evolve.py`

Covered by Stage 0:

- `directed_evolution_simulation` with `ridge`, `topn`, `random` first round.
- `evolve_experimental` with deterministic synthetic embeddings and Colab Kelsic Round1 data.

Reachable but not fully covered:

- `grid_search`.
- `evolve_experimental_multi`.
- `directed_evolution_simulation` learning strategies: `dist`, `random`, `topn2bottomn2`.

Deletion decision:

- Keep all through CLI wrapping. `grid_search` is the documented DMS entrypoint and `evolve_experimental_multi` is documented by `scripts/exp/README.md` and `scripts/exp/t7_pol.py`.
- Add CLI-level tests before deleting legacy scripts.

### `evolvepro/src/data.py`

Covered by Stage 0:

- `load_experimental_data`.
- `create_iteration_dataframes`.
- `add_missing_variants` through `create_iteration_dataframes`.
- `process_variant` through Kelsic Excel loading.

Reachable but not fully covered:

- `load_experimental_embeddings` through `evolve_experimental`.
- `load_dms_data` CSV path through `grid_search`, not directly tested yet.
- `.pt` embedding path via `torch.load` and `process_pt_embeddings`.

Required migration change:

- Move `import torch` from module import time into the `.pt` branch of `load_dms_data`.
- Keep `.pt` support optional so `evolvepro --help`, CSV/Excel processing, and DMS CSV workflows work without installing torch.

Deletion decision:

- Do not remove `.pt` support until the PLM/embed packaging strategy is settled.

### `evolvepro/src/process.py`

Covered by Stage 0:

- `generate_wt`.
- `generate_single_aa_mutants`.
- `suggest_initial_mutants`.
- `process_dataset` on a checked-in SARS-CoV-2 DMS file.
- `generate_mutation_fasta` through `process_dataset`.

Reachable but not fully covered:

- `generate_n_mutant_combinations`.
- DMS preprocessors: `preprocess_cas12f`, `preprocess_zikv_E`, `preprocess_cov2_S`.
- Plot helpers in this module: `plot_mutations_per_position`, `plot_histogram_of_readout`.
- `markin_custom_cutoff`.

Deletion decision:

- Keep through CLI wrapping. These are documented by `scripts/process/README.md` and `scripts/process/dms_process.py`.
- Prefer adding small tests before refactoring/deleting individual helpers.

### `evolvepro/src/plot.py`

Covered by Stage 0:

- `read_exp_data` through Colab Kelsic rounds.

Reachable but not fully covered:

- `read_dms_data`.
- `process_dataframe`.
- `filter_dataframe`.
- `save_dataframe`.
- `apply_labels`.
- `load_external_data`.
- `concatenate_dataframes`.
- plotting functions.

Deletion decision:

- Keep through CLI wrapping because `scripts/plot/README.md` documents these functions.
- Plot rendering can stay smoke-only; dataframe transform helpers should get unit tests when `evolvepro plot` is implemented.

### `evolvepro/plm/*` and `scripts/plm/*`

Covered by Stage 0:

- Static reachability only; no model download/GPU execution.

Deletion decision:

- Keep until `evolvepro embed` exists and supports the documented PLM backends, or until project scope explicitly reduces supported backends.
- PLM dependencies, including torch, should be optional and lazy-loaded.

## Data and notebook assets

### Test fixtures, keep until final cleanup

- `colab/rounds_data/kelsic_Round*.xlsx`.
- `data/dms/activity/sarscov2_preprocessed.csv`.
- `data/dms/wt_fasta/cov2_S_WT.fasta`.
- `data/dms/activity/*` and `data/dms/wt_fasta/*` more broadly until DMS examples are externalized.

### Documentation/tutorial assets, defer deletion

- `colab/EvolvePro.ipynb`.
- `images/*` while README still references them.
- `data/exp/Evolution_Summary_Data.xlsx` until README/docs are rewritten.

## Stage 0b decision

No code deletion in Stage 0b.

Rationale:

1. The upstream examples intentionally expose DMS, experimental, plotting, and PLM workflows.
2. The new 22-test baseline locks enough behavior to start refactoring, but not enough to justify deleting whole workflow branches.
3. Removing code now would silently narrow upstream compatibility before the CLI replacement exists.

Next implementation stages:

1. Make `torch` lazy/optional so basic imports and CLI help work without the PLM stack.
2. Add standardized round dump artifacts to `evolve_experimental`/`top_layer` paths without changing existing CSV outputs.
3. Introduce config schemas and Typer wrappers.
4. Delete legacy scripts only after CLI tests replace their contracts.
