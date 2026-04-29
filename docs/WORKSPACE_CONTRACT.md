# Workspace, round, and manifest contract

This contract defines the target CLI/workspace shape for the migration. It is intentionally file-based. SQLite or another database must not be the source of truth.

## Decision

Use per-workspace config, lock/index files, per-round manifests, and CAS snapshots as the authoritative state:

- `evolvepro.toml` records stable workspace policy for one protein.
- `evolvepro.lock` records mutable workspace index state and the latest round pointer.
- `rounds/<round-id>/round.lock` records resolved immutable inputs for one round transition.
- `rounds/<round-id>/manifest.json` records outputs, hashes, and round lineage.
- `.evolvepro/cas/sha256/...` stores immutable snapshots of every input/output artifact used for validation.

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
  evolvepro.lock

  inputs/
    embeddings.csv

  rounds/
    round_000/
      candidates.csv          # mutable user handoff
      round.lock              # immutable resolved snapshot
      manifest.json           # immutable artifact ledger

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

    round_001/
      candidates.csv
      round.lock
      manifest.json
      ...

  .evolvepro/
    cas/
      sha256/
    cache/
```

## Command model

Target CLI flow:

```bash
evolvepro init [DIR]
evolvepro evolve
evolvepro validate
evolvepro plot
```

`evolvepro init` owns workspace setup and interactive configuration. A separate public `configure` command is not part of the target UX.

`evolvepro evolve` advances the workspace state machine:

1. Validate previous rounds recursively.
2. If no round exists, generate the initial candidate CSV.
3. If a candidate CSV has measured activity values, ingest measured rows and record blank rows as unmeasured.
4. Train/update the model when enough measured data exists.
5. Predict the next candidates.
6. Write CAS snapshots, `round.lock`, `manifest.json`, and update `evolvepro.lock`.

Current `evolvepro run --config` remains a compatibility scaffold until this target flow is implemented.

## File roles

### `evolvepro.toml`

Workspace-level user-facing config. It is editable and suitable for version control.

It should contain defaults such as:

- schema version
- workflow
- output locations
- model defaults
- path defaults

### `evolvepro.lock`

Workspace-level mutable index. It is generated and updated by EvolvePro.

It should contain:

- schema version
- workspace id
- protein id/name
- latest round pointer
- list of round manifest paths and hashes
- global CAS/index metadata required to validate the workspace

It should not contain per-round observations or experiment outcomes. Those belong in round artifacts and manifests.

### `round.lock`

Resolved immutable round-transition inputs. It is generated immediately before a state transition and should not be edited by users.

It must include:

- schema version
- resolved config hash
- resolved relative paths
- raw input byte hashes
- previous round manifest hash
- resolved model parameters with defaults
- package/environment metadata required for reproducibility

The lock file answers: "What exact inputs and settings were used for this round transition?"

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

1. `evolvepro evolve` writes `rounds/round_N/candidates.csv` with `variant,activity` columns.
2. A human fills measured activity values directly into that CSV.
3. Blank activity cells are treated as unmeasured candidates.
4. The next `evolvepro evolve` snapshots the edited CSV into CAS.
5. EvolvePro canonicalizes/normalizes measured activity rows.
6. The model trains on accumulated measured rounds.
7. The model predicts candidates for the next round.
8. A manifest records every input/output artifact hash and links to the previous round manifest.

Every `evolve` invocation must verify previous round manifest hashes before advancing.

## CSV handoff policy

Primary user handoff files are CSV only. The supported schema is exactly:

```csv
variant,activity
A12V,
G18D,0.42
```

Rules:

- exactly two columns: `variant` and `activity`
- UTF-8 comma-delimited CSV
- duplicate variants are invalid
- blank activity means unmeasured
- non-blank activity must parse as numeric
- partial measurement is normal state

XLSX is not part of the primary handoff contract. XLSX remains allowed only for legacy regression tests and upstream compatibility fixtures during migration. If legacy XLSX import remains available, it must immediately normalize to canonical CSV and validation must proceed from the CSV artifact.

## Mutable handoff and immutable artifacts

`rounds/<round-id>/candidates.csv` is writable by the user. It is a handoff file, not the immutable provenance record.

Before a file participates in validation, EvolvePro snapshots it into `.evolvepro/cas/sha256/...` and records the CAS hash in `round.lock` and `manifest.json`.

Generated immutable files should be written atomically and then made read-only where the filesystem supports it:

- CAS snapshots
- `round.lock`
- `manifest.json`

Read-only permissions are guardrails against accidental edits. Hash validation is the integrity mechanism.

## Hash policy

Hash validation is required before advancing rounds.

### Raw inputs

Use byte SHA256 for raw user inputs:

- CSV handoff files
- embedding files
- explicit config files

A hash mismatch means an immutable snapshot or previously recorded input changed. The run must stop. Mutable handoff CSV files are snapshotted as new input artifacts during the next `evolve` transition.

### Normalized tables and predictions

Use semantic hashes in addition to byte hashes for CSV handoff files and normalized tables.

Canonicalization must define:

- column order
- row sort order
- dtype handling
- float formatting/precision
- null representation
- UTF-8 encoding

This avoids false differences from Parquet metadata, writer versions, row group layout, or allowed CSV formatting differences.

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

GPU embedding generation is target first-class workflow state. The target `evolve` state machine reuses precomputed embedding files when present and can prepare missing embeddings through configured backends. Generated embedding provenance records backend, model name, model revision, dtype, device metadata, seed, input FASTA hash, output embedding hash, metadata hash, CAS paths, and embedding config semantic hash. CPU/no-PLM paths remain compatibility and test paths, not the primary UX.

## Compatibility output policy

During migration, config-driven runs write both the historical root-level outputs and the planned round layout.

Historical root-level outputs are compatibility artifacts:

```text
rounds/<round-id>/iteration.csv
rounds/<round-id>/this_round_variants.csv
rounds/<round-id>/df_test.csv
rounds/<round-id>/df_sorted_all.csv
rounds/<round-id>/model.pkl
rounds/<round-id>/predictions.parquet
rounds/<round-id>/metrics.json
rounds/<round-id>/inputs.json
rounds/<round-id>/version.txt
```

Planned layout outputs are the forward-looking CLI contract:

```text
rounds/<round-id>/legacy/iteration.csv
rounds/<round-id>/legacy/this_round_variants.csv
rounds/<round-id>/legacy/df_test.csv
rounds/<round-id>/legacy/df_sorted_all.csv
rounds/<round-id>/model/model.pkl
rounds/<round-id>/model/params.json
rounds/<round-id>/predictions/predictions.parquet
rounds/<round-id>/predictions/candidates.csv
```

Compatibility rules:

- Root-level CSV outputs remain byte-identical to the Stage 0 golden baseline.
- Planned-layout files are copies or deterministic derivatives of root-level compatibility files.
- `manifest.json` records hashes for both compatibility outputs and planned-layout outputs.
- Root-level outputs must remain until all legacy scripts and public examples have CLI replacements with regression coverage.
- After CLI migration is complete, root-level outputs may be deprecated only through a documented compatibility window.
- The planned layout is the preferred target for new CLI and interactive UI code.

## SQLite policy

Do not implement SQLite in the next stages.

A future SQLite index is allowed only if all of the following are true:

- It lives under `.evolvepro/`.
- It is explicitly documented as cache/index state.
- It can be deleted and fully rebuilt from TOML/JSON/manifest files.
- It is never required to reproduce a run.
