# Workspace, round, and manifest contract

This contract defines the target CLI/workspace shape for the migration. It is intentionally file-based. SQLite or another database must not be the source of truth.

## Decision

Use per-workspace and per-round TOML/JSON files as the authoritative state:

- `evolvepro.toml` records workspace-level defaults and interactive choices.
- `rounds/<round-id>/round.toml` records user-editable round intent.
- `rounds/<round-id>/round.lock.json` records resolved immutable inputs for one run.
- `rounds/<round-id>/manifest.json` records outputs, hashes, and round lineage.

Do not introduce SQLite as primary state. If a future interactive UI needs faster listing/searching, an index database may be added only as a rebuildable cache under `.evolvepro/`.

## Rationale

File state is preferred because it provides:

- Human-readable review via Git diffs.
- Easy archival with experiment artifacts.
- Explicit reproducibility through config, lock files, manifests, and hashes.
- Stable idempotent writes through canonical TOML/JSON serialization.
- Round-level corruption isolation.
- Simple CI and golden-file regression testing.

SQLite primary state is rejected because it makes review, archival, idempotency, and schema migration harder for this experiment workflow.

## Target workspace layout

```text
workspace/
  evolvepro.toml

  inputs/
    wt.fasta
    embeddings.csv

  rounds/
    round-001/
      round.toml
      round.lock.json
      manifest.json

      inputs/
        activity.raw.xlsx
        activity.normalized.parquet

      train/
        training_set.parquet
        training_manifest.json

      model/
        model.pkl
        params.json

      predictions/
        predictions.parquet
        candidates.csv

      legacy/
        iteration.csv
        this_round_variants.csv
        df_test.csv
        df_sorted_all.csv

    round-002/
      round.toml
      round.lock.json
      manifest.json
      ...

  .evolvepro/
    cache/
    locks/
```

## Command model

Planned CLI flow:

```bash
evolvepro init [DIR]
evolvepro
```

`evolvepro` without subcommands should default to interactive mode:

1. Detect or create a workspace.
2. Let the user choose the command/workflow.
3. Let the user choose input files and round settings.
4. Write deterministic TOML.
5. Run from TOML after confirmation.

Non-interactive mode must be supported:

```bash
evolvepro run --config evolvepro.toml
```

Interactive mode is only a config generator/editor. Execution must always resolve through TOML and lock/manifest files.

## File roles

### `evolvepro.toml`

Workspace-level user-facing config. It is editable and suitable for version control.

It should contain defaults such as:

- schema version
- workflow
- output locations
- model defaults
- path defaults

### `round.toml`

Round-level user intent. It is editable and suitable for version control.

Example shape:

```toml
schema_version = 1
round = "round-001"
workflow = "experimental_round"

[inputs]
activity = "inputs/activity.raw.xlsx"

[model]
regression_type = "randomforest"
seed = 0
```

Paths are interpreted relative to the TOML file location unless documented otherwise.

### `round.lock.json`

Resolved immutable run inputs. It is generated immediately before a run and should not be edited by users.

It must include:

- schema version
- resolved config hash
- resolved relative paths
- raw input byte hashes
- previous round manifest hash
- resolved model parameters with defaults
- package/environment metadata required for reproducibility

The lock file answers: "What exact inputs and settings were used for this run?"

### `manifest.json`

Run result ledger. It is generated after a successful run and should not be edited by users.

It must include:

- lock file hash
- raw input hashes
- normalized input semantic hashes
- training table semantic hash
- model parameter semantic hash
- model checkpoint byte hash
- prediction semantic hash
- legacy output file hashes
- previous round manifest hash

The manifest answers: "What artifacts were produced, and how do they link to the previous round?"

## Round progression

EvolvePro is a round-based workflow:

1. A human measures activity for the current candidates.
2. The activity file is added to the next round directory.
3. EvolvePro canonicalizes/normalizes the activity input.
4. The model trains on accumulated measured rounds.
5. The model predicts candidates for the next round.
6. A manifest records the output and links to the previous round manifest.

The next round must verify the previous round manifest hash before running.

## Hash policy

Hash validation is required before advancing rounds.

### Raw inputs

Use byte SHA256 for raw user inputs:

- FASTA
- CSV
- XLSX
- embedding files
- explicit config files

A hash mismatch means the input changed. The run must stop unless the user explicitly accepts the change through a future command such as `--accept-input-changes`.

### Normalized tables and predictions

Use semantic hashes rather than Parquet/CSV byte hashes for tables.

Canonicalization must define:

- column order
- row sort order
- dtype handling
- float formatting/precision
- null representation
- UTF-8 encoding

This avoids false differences from Parquet metadata, writer versions, row group layout, or CSV formatting.

### Model checkpoint

Record both:

- `model.pkl` byte SHA256
- `params.json` semantic SHA256

The checkpoint byte hash is useful, but long-term reproducibility should depend primarily on:

1. resolved model parameters
2. training table semantic hash
3. prediction semantic hash
4. environment/package metadata

## Determinism policy

A run is deterministic if the following are unchanged:

- TOML config semantic hash
- raw input byte hashes
- previous round manifest hash
- resolved model parameters
- package/runtime environment within supported bounds

Model training must expose and persist all stochastic parameters. Random forest must use explicit `random_state` and deterministic defaults for config-driven runs. Prefer `n_jobs = 1` for strict determinism.

GPU embedding generation is outside the first deterministic run contract. The first CLI run contract consumes precomputed embedding files and verifies their hashes. A later embedding-preparation workflow must write its own manifest with model name, model revision, dtype, device, torch/CUDA metadata, and seed.

## SQLite policy

Do not implement SQLite in the next stages.

A future SQLite index is allowed only if all of the following are true:

- It lives under `.evolvepro/`.
- It is explicitly documented as cache/index state.
- It can be deleted and fully rebuilt from TOML/JSON/manifest files.
- It is never required to reproduce a run.
