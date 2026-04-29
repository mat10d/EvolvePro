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
