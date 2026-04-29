# Script migration map

This document maps the checked-in `scripts/` surface to the new CLI. It is a migration aid, not a deletion plan. Do not delete scripts until their replacement command has regression coverage and public docs point to the CLI.

## Current CLI interface

The current CLI is intentionally non-interactive and provisional. It is implemented with `argparse`; `rich-click` interactive UX is deferred. The target user-facing interface is documented in `docs/CLI_INTERFACE.md` and is centered on `init`, `evolve`, `validate`, and `plot`.

| Command | Purpose | Status | Tests |
|---|---|---|---|
| `evolvepro init [DIR] [--force]` | Initialize a file-based workspace. | Implemented | `tests/test_cli.py`, `tests/test_workspace.py` |
| `evolvepro run [--config evolvepro.toml]` | Run an experimental round from TOML. | Implemented | `tests/test_cli.py`, `tests/test_cli_workflows.py` |
| `evolvepro validate [ROUND_DIR]` | Validate a completed round's lock/manifest/hash chain. | Implemented | `tests/test_cli.py`, `tests/test_config.py` |
| `evolvepro process generate-wt --sequence SEQ --output FILE` | Write WT FASTA. | Implemented | `tests/test_cli_process.py` |
| `evolvepro process generate-single-aa-mutants --wt-fasta FILE --output FILE [--positions 1,2,3]` | Write single-AA mutant FASTA. | Implemented | `tests/test_cli_process.py` |
| `evolvepro process dataset --file-path FILE --wt-fasta FILE --dataset-name NAME --activity-column COL --cutoff-value FLOAT --output-dir DIR [--sheet-name NAME] [--cutoff-rule RULE] [--cutoff-percentiles ...] [--aa-shift INT] [--drop-columns]` | Process a DMS/CSV/Excel activity dataset into labels CSV and FASTA. | Implemented | `tests/test_cli_process_dataset.py` |
| `evolvepro dms grid-search ...` | Run DMS simulation grid search. | Implemented | `tests/test_cli_dms.py` |

## Interface decision

The current nested command shape is too low-level for primary UX:

```bash
evolvepro process generate-single-aa-mutants ...
evolvepro process dataset ...
evolvepro dms grid-search ...
```

These commands mirror upstream script directories. They remain compatibility scaffolding, not the target interface.

Target primary UX:

```bash
evolvepro init [DIR]
evolvepro evolve
evolvepro validate
evolvepro plot
```

`evolve` absorbs initial candidate generation, measured activity ingest, model training, prediction, checkpointing, artifact writing, lock updates, and implicit validation. CSV handoff files are mutable user inputs; immutable provenance is maintained through CAS snapshots, `round.lock`, `manifest.json`, and `evolvepro.lock`.

Recommendation for now:

- Keep existing commands stable while target commands are built.
- Do not add more primary `process ...` commands.
- Treat `process ...` as low-level compatibility utilities.
- Migrate README and scripts toward `init`/`evolve`/`validate`/`plot` once implemented.

## Migration/deletion gate

A script is not deletable until all are true:

1. CLI replacement exists.
2. CLI replacement has subprocess-level regression coverage.
3. README or script docs point users to the CLI replacement.
4. Existing legacy outputs remain covered by tests or compatibility policy.
5. `tests/REACHABLE.md` is updated with the deletion decision.

Until then, scripts are compatibility examples.

## `scripts/process` map

| Script | Role | Current CLI replacement | Coverage | Delete now? |
|---|---|---|---|---|
| `scripts/process/exp_process.py` | Example calls for WT FASTA, single-AA mutants, initial mutant suggestions, n-mutant combinations. | Partial: `evolvepro process generate-wt`, `evolvepro process generate-single-aa-mutants`. | CLI tests for FASTA commands. Existing function tests for `suggest_initial_mutants`. | No |
| `scripts/process/dms_process.py` | Hard-coded DMS dataset preprocessing examples. | Partial: `evolvepro process dataset`. | CLI test for SARS-CoV-2 dataset. Existing DMS data test. | No |

Missing process replacements:

- `suggest_initial_mutants`
- `generate_n_mutant_combinations`
- dataset-specific helpers such as `preprocess_cas12f`
- plotting helpers in `evolvepro/src/process.py`

## `scripts/dms` map

| Script | Role | Current CLI replacement | Coverage | Delete now? |
|---|---|---|---|---|
| `scripts/dms/dms_main.py` | argparse wrapper around `grid_search`. | `evolvepro dms grid-search`. | CLI test compares output to API. Existing parser smoke test. | No |
| `scripts/dms/*.sh` | SLURM/grid examples for specific embedding/model combinations. | Partial: command body maps to `evolvepro dms grid-search`; scheduling remains external. | No direct CLI coverage for each shell script. | No |

SLURM shell scripts should likely remain as examples until a separate batch/run-template story exists.

## `scripts/exp` map

| Script | Role | Current CLI replacement | Coverage | Delete now? |
|---|---|---|---|---|
| `scripts/exp/t7_pol.py` | Hard-coded experimental round example. | Conceptual replacement: `evolvepro run --config`. | Static script contract; Kelsic config-driven CLI workflow test. | No |
| `scripts/exp/bxb1.py` | Hard-coded experimental round example. | Conceptual replacement: `evolvepro run --config`. | Static script contract only. | No |
| `scripts/exp/mlv.py` | Hard-coded experimental round example. | Conceptual replacement: `evolvepro run --config`. | Static script contract only. | No |

The better replacement for `scripts/exp` is probably config templates, not one CLI command per script.

Potential future command:

```bash
evolvepro config experimental-round --template kelsic --output evolvepro.toml
```

## `scripts/plot` map

| Script | Role | Current CLI replacement | Coverage | Delete now? |
|---|---|---|---|---|
| `scripts/plot/dms.py` | DMS plotting example. | None. | Static/docs only. | No |
| `scripts/plot/exp.py` | Experimental plotting example. | None. | Static/docs only. | No |

Plot CLI should wait until output layout stabilizes.

Potential future commands:

```bash
evolvepro plot dms --input ... --output ...
evolvepro plot experimental --round-dir ... --output ...
```

## `scripts/plm` map

| Script | Role | Current CLI replacement | Coverage | Delete now? |
|---|---|---|---|---|
| `scripts/plm/*.sh` | Embedding extraction examples for PLM backends. | None. | Smoke/static only through documentation. | No |

PLM commands are intentionally deferred because they involve optional torch/GPU/model-download dependencies. They should remain outside the core CSV handoff workflow contract until a separate embedding manifest is designed.

Potential future command group:

```bash
evolvepro embeddings extract ...
evolvepro embeddings validate ...
```

## Recommended next work

1. Implement the target `evolvepro evolve` state machine before adding more low-level `process ...` commands.
2. Keep existing low-level commands as compatibility coverage until their behavior is absorbed by `evolve`.
3. For `scripts/exp`, migrate to workspace-driven `evolve` rather than script-specific commands.
4. Keep PLM and plot workflows deferred until the core workspace state machine is stable.
