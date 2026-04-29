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
- `docs/CLI_INTERFACE.md` for the target user-facing CLI (`init`, `evolve`, `validate`, `plot`, optional `dms`).
- `docs/WORKSPACE_CONTRACT.md` for the planned file-based workspace, round, lock, and manifest contract.
- `docs/INTERACTIVE_CLI.md` for the deferred rich-click interactive CLI plan.
- `docs/SCRIPT_MIGRATION.md` for the script-to-CLI replacement map and compatibility command table.

## Completed commits on `main`

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
3b7f3a5 docs: record migration state after dump artifacts
1dde955 docs: define file-based workspace contract
035c4ab feat: add file-based config and lock primitives
79559ef feat: write round lock and manifest for config runs
d514b91 feat: plumb model seed and parameters
5ba7254 feat: materialize planned round layout
f95cbce feat: add workspace initialization primitives
c835ab9 feat: validate completed round manifests
7658092 feat: add minimal workspace CLI
872b3fc fix: install evolvepro console script in uv env
3d597bb test: lock experimental CLI workflow
80bcb5a feat: add process FASTA CLI commands
0bc484f feat: add DMS grid search CLI
221cc5c feat: add process dataset CLI
d4a5f7c docs: map scripts to CLI replacements
d1ef6c1 docs: define target CLI interface
4b86990 docs: limit xlsx to compatibility paths
4aaccfb feat: start workspace evolve state machine
ed8c0c8 feat: advance workspace evolve rounds
1428ef9 feat: train workspace evolve from embeddings
5beacd9 docs: make embedding workflow gpu first
79ecfde feat: hash workspace embedding inputs
2df0da7 feat: add gpu-first embedding config
a6570b9 feat: scaffold embedding backend api
```

## Current validation command

Run after every migration commit:

```bash
nix develop -c uvx ruff check .
nix develop -c uv run pytest tests/ -v
```

Current expected result:

```text
86 passed, 26 warnings
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
- `torch` is intentionally not a required runtime dependency for core CSV workflows.
- `evolvepro/src/data.py` now lazy-imports `torch` only for `.pt` embedding files.
- Keep CLI/subprocess implementations from overwriting `LD_LIBRARY_PATH`; Nix devShell sets it for native libraries and host CUDA driver paths.
- Target workflow is GPU-first for PLM embedding extraction, but GPU/model-download paths remain optional or smoke-tested in normal CI.

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

## Current workspace/config/manifest state

Read `docs/WORKSPACE_CONTRACT.md` first after context compaction.

Implemented in `035c4ab` and `79559ef`:

- deterministic canonical JSON bytes and JSON SHA256
- table semantic hashes with explicit sort columns
- TOML loading into typed config dataclasses
- deterministic TOML writing
- config-driven `evolve_experimental` adapter
- round lock construction with raw input hashes, model config, config hash, and previous manifest hash
- round lock validation that stops on input or previous-manifest hash mismatch
- config-driven runs write `round.lock.json` and `manifest.json` into the legacy round output directory
- manifest captures lock hash, raw input hashes, model checkpoint hash, prediction file hash, prediction semantic hash, dump file hashes, and legacy CSV hashes
- config-driven runs pass model seed and model parameters into training; legacy calls without these arguments keep upstream defaults
- config-driven runs materialize planned layout copies under `legacy/`, `model/`, and `predictions/` while keeping root-level legacy/dump outputs intact
- `model/params.json` persists resolved model parameters from the trained estimator
- workspace primitives initialize `evolvepro.toml`, `evolvepro.lock`, `.evolvepro/workspace.toml`, `inputs/`, `outputs/`, `rounds/`, `.evolvepro/cache/`, `.evolvepro/locks/`, and `.evolvepro/cas/sha256/`
- workspace discovery walks parent directories from an arbitrary path
- completed round validation checks `round.lock.json`, `manifest.json`, input hashes, output hashes, prediction semantic hash, and previous manifest hash chain
- minimal non-interactive CLI exists with `evolvepro init`, `evolvepro evolve`, `evolvepro run --config`, and `evolvepro validate`
- `uv sync --frozen` installs the project editable and exposes `.venv/bin/evolvepro`; `setup.py` is now a shim over `pyproject.toml`
- README documents the experimental CLI workspace workflow
- subprocess-level CLI workflow test verifies `evolvepro run --config` and `evolvepro validate` against config API output
- process FASTA generation has CLI coverage via `evolvepro process generate-wt` and `evolvepro process generate-single-aa-mutants`
- DMS grid search has CLI coverage via `evolvepro dms grid-search` and README documentation
- DMS/process dataset preparation has CLI coverage via `evolvepro process dataset`
- `evolvepro init` accepts target setup flags `--set-default`, `--reconfigure`, and `--yes`; `--force` remains as compatibility overwrite
- `evolvepro evolve` starts the target workspace state machine by generating `round_000/candidates.csv` from the configured WT FASTA, snapshotting it into CAS, writing read-only `round.lock` and `manifest.json`, and updating `evolvepro.lock`
- `evolvepro evolve` now ingests partially filled candidate CSV files, records measured/unmeasured counts, snapshots the edited handoff into CAS, and writes the next round candidate handoff
- CSV candidate handoff helpers enforce exact `variant,activity` columns, numeric filled activity values, duplicate rejection, and blank activity as unmeasured
- CAS helpers write content-addressed read-only snapshots under `.evolvepro/cas/sha256/`
- workspace validation now checks `evolvepro.lock`, `round.lock`, manifest hashes, previous-manifest links, and CAS snapshots while allowing mutable handoff files to change
- when `inputs/embeddings.csv` exists, `evolvepro evolve` now trains the configured top-layer model on accumulated measured rows, ranks unmeasured embedded variants, and writes model/training/prediction artifacts
- core top-layer training is CPU-bound scikit-learn/XGBoost style regression; GPU is not used in the train step
- target workflow is GPU-first for PLM embedding extraction, but embedding extraction is still not part of `evolve`; current implementation consumes precomputed embedding CSV files
- embedding CSV inputs used by `evolvepro evolve` are recorded in `round.lock`/`manifest.json` with raw file hash and semantic table hash; workspace validation rejects embedding file mutation
- `EmbeddingConfig` is part of typed config loading/writing and default workspaces now include GPU-first embedding settings without importing torch
- embedding backend API scaffold exists with `EmbeddingRequest`, `EmbeddingResult`, and `EmbeddingBackend` protocol; device, FASTA, normalization, fake backend, and ESM backend can now branch independently

Key decisions:

- Do not use SQLite as source of truth. A future SQLite index may only be a rebuildable `.evolvepro/` cache.
- Target user-facing commands are `evolvepro init`, `evolvepro evolve`, `evolvepro validate`, `evolvepro plot`, plus optional `evolvepro dms ...`.
- A separate public `configure` command is not planned; `init` owns setup and `init --reconfigure` owns overrides.
- One workspace manages one protein.
- User handoff files are CSV-only with `variant,activity` columns; XLSX support is limited to legacy regression/compatibility paths and must normalize to CSV before target provenance.
- `evolve` owns initial candidate generation, measured activity ingest, model training, prediction, checkpointing, implicit recursive validation, and state updates.
- Mutable candidate CSV files are handoff files; immutable provenance comes from CAS snapshots, `round.lock`, `manifest.json`, and `evolvepro.lock`.
- Use `.lock` filenames (`evolvepro.lock`, `round.lock`) for target lock files, even if internal encoding remains deterministic JSON.

Implemented/planning files:

- `docs/GPU_FIRST_EMBEDDING_PLAN.md`
- `evolvepro/embeddings/__init__.py`
- `evolvepro/embeddings/schemas.py`
- `evolvepro/embeddings/backends/__init__.py`
- `evolvepro/embeddings/backends/base.py`
- `evolvepro/io/schemas.py`
- `evolvepro/io/config.py`
- `evolvepro/io/manifest.py`
- `evolvepro/io/canonical.py`
- `tests/fixtures/kelsic_round1.toml`
- `tests/test_canonical.py`
- `tests/test_config.py`
- `tests/test_manifest.py`
- `evolvepro/io/cas.py`
- `evolvepro/io/handoff.py`
- `evolvepro/workflows.py`
- `evolvepro/validation.py`
- `tests/test_handoff_and_cas.py`
- `tests/test_evolve_workspace.py`
- `tests/test_workspace_validate.py`

Compatibility output policy is documented in `docs/WORKSPACE_CONTRACT.md`:

- root-level outputs are migration compatibility artifacts and must remain until all legacy scripts/examples have CLI replacements with regression coverage
- planned layout is the forward-looking CLI contract
- `manifest.json` records both compatibility and planned-layout outputs

Interactive CLI plan is documented in `docs/INTERACTIVE_CLI.md`. Keep `argparse` for the current compatibility contract; defer `rich-click` until the target workflow state machine is stable.

Next implementation stage:

1. Implement M2 parallel branches from `docs/GPU_FIRST_EMBEDDING_PLAN.md`:
   - device resolver with lazy torch import
   - canonical embedding FASTA generation
   - embedding CSV normalization and CAS helpers
   - deterministic fake backend without torch/GPU/network in CI
   - ESM backend wrapper smoke-tested with monkeypatching only
2. Implement M3: integrate embedding preparation into `evolvepro evolve`, reuse valid embeddings, and record embedding artifacts in `round.lock`/`manifest.json`.
3. Implement M4: finalize ESM extraction with GPU-first device resolution and explicit CPU fallback only.
4. Implement M5: CLI polish and compatibility bridge for legacy `scripts/plm/*`.
5. Persist normalized observation/training tables as artifacts with semantic hashes in both lock and manifest payloads.
6. Add stronger workspace validation coverage for manifest mutation, lock mutation, and previous-manifest chain breakage.

Important: there is no `data/<enzyme>/config.py` legacy adapter in this repo. The adapter targets current function arguments and the evolving argparse CLI, not the obsolete BenM/TphR layout.
