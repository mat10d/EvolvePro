# GPU-first embedding implementation plan

This is the next migration phase after workspace model-backed `evolve` and embedding input hashes.

## Goal

Make PLM embedding extraction a first-class, GPU-first workflow phase owned by `evolvepro evolve`.

Current state:

- `evolvepro evolve` trains from precomputed `inputs/embeddings.csv` when present.
- The top-layer regression model runs on CPU.
- Legacy PLM scripts can extract embeddings, but they are outside the workspace state machine.

Target state:

- `evolvepro evolve` ensures embeddings exist before model-backed candidate ranking.
- GPU is preferred and explicit for PLM extraction.
- CPU/no-PLM remains a compatibility/test path, not the primary UX.
- Every FASTA, embedding file, parameter set, device decision, and environment detail is captured in lock/manifest/CAS.

## Non-goals for this phase

- Do not make torch a mandatory dependency for core CSV/help/validation workflows.
- Do not run model downloads, GPU jobs, or full PLM extraction in normal CI.
- Do not replace the CPU top-layer regression path unless a future backend explicitly supports GPU training.
- Do not delete legacy `scripts/plm/*` before CLI parity and regression coverage exist.

## User-facing behavior

Primary workflow stays high-level:

```bash
evolvepro init my-protein
evolvepro evolve
```

When embeddings are missing or stale, `evolve` should prepare them before training. In GPU-first mode:

1. Resolve configured PLM backend/model/device.
2. Generate a canonical FASTA for the candidate universe or current embedding target set.
3. Snapshot FASTA into CAS.
4. Run embedding extraction on GPU.
5. Normalize embedding output to canonical CSV.
6. Snapshot embedding CSV into CAS.
7. Record metadata in `round.lock` and `manifest.json`.
8. Continue CPU top-layer training/ranking.

Failure policy:

- `device = "auto"` prefers CUDA.
- If CUDA is unavailable, fail with an actionable message unless explicit CPU compatibility mode is enabled.
- If torch/backend dependencies are missing, fail only when embedding extraction is required.
- Existing valid embeddings may be reused without importing torch.

## Config additions

Add an embedding section to `evolvepro.toml`:

```toml
[embedding]
backend = "esm"
model_name = "esm2_t33_650M_UR50D"
# model_revision is omitted unless a backend-specific revision is pinned.
representation = "mean"
repr_layers = [-1]
toks_per_batch = 512
device = "auto"         # auto | cuda | cuda:N | cpu
dtype = "float32"
allow_download = false
allow_cpu_fallback = false
cache_dir = ".evolvepro/cache/plm"
```

Compatibility/test config may use:

```toml
[embedding]
backend = "none"        # skip extraction; deterministic fallback remains allowed
```

or:

```toml
[embedding]
backend = "fake"        # deterministic test backend, never public default
```

## Internal API

Introduce an embedding package instead of wiring legacy scripts directly into the state machine. The package scaffold exists; device, FASTA, normalization, and backend implementations are filled in by later milestones.

```text
evolvepro/embeddings/
  __init__.py
  schemas.py
  device.py
  fasta.py
  normalize.py
  backends/
    base.py
    esm.py
    fake.py
```

Core dataclasses:

```python
@dataclass(frozen=True)
class EmbeddingRequest:
    fasta_path: Path
    output_dir: Path
    backend: str
    model_name: str
    model_revision: str | None
    repr_layers: list[int]
    representation: str
    toks_per_batch: int
    device: str
    dtype: str
    seed: int

@dataclass(frozen=True)
class EmbeddingResult:
    embeddings_csv: Path
    backend_metadata: dict[str, object]
    device_metadata: dict[str, object]
```

Backend contract:

```python
class EmbeddingBackend(Protocol):
    def extract(self, request: EmbeddingRequest) -> EmbeddingResult: ...
```

The ESM backend should wrap/import `evolvepro/plm/esm/extract.py` logic initially. Later backends can be added without changing the workspace state machine.

## Device resolution

Add `evolvepro/embeddings/device.py`:

- Lazy import torch.
- `resolve_device("auto", allow_cpu_fallback=False)`:
  - return `cuda:0` if available
  - otherwise raise `EmbeddingDeviceError`
- `resolve_device("cuda")`:
  - return `cuda:0` if available
  - otherwise raise
- `resolve_device("cuda:N")`:
  - validate device index
- `resolve_device("cpu")`:
  - allowed only when `allow_cpu_fallback = true` or backend is test/compatibility

Record:

```json
{
  "requested_device": "auto",
  "resolved_device": "cuda:0",
  "cuda_available": true,
  "cuda_device_count": 1,
  "cuda_device_name": "...",
  "torch_version": "...",
  "cuda_version": "..."
}
```

## Artifact contract

Embedding preparation writes under the round that first needs the embeddings:

```text
rounds/round_001/
  embeddings/
    candidates.fasta
    embeddings.csv
    metadata.json
```

All three files are recorded in CAS and manifest outputs. If an existing embedding CSV is reused from `inputs/embeddings.csv`, it remains a manifest input with raw and semantic hashes.

Manifest example:

```json
{
  "inputs": {
    "embedding_fasta": {
      "path": "rounds/round_001/embeddings/candidates.fasta",
      "sha256": "...",
      "cas_path": ".evolvepro/cas/sha256/..."
    },
    "embedding_config": {
      "semantic_sha256": "..."
    }
  },
  "outputs": {
    "embeddings": {
      "path": "rounds/round_001/embeddings/embeddings.csv",
      "sha256": "...",
      "semantic_sha256": "...",
      "cas_path": ".evolvepro/cas/sha256/..."
    },
    "embedding_metadata": {
      "path": "rounds/round_001/embeddings/metadata.json",
      "sha256": "..."
    }
  }
}
```

## Milestones

### M1: Schema and docs

- Add `EmbeddingConfig` to config schemas and TOML loader/writer.
- Default workspace config uses GPU-first ESM settings but does not import torch during `init`.
- Document device/fallback behavior.
- Tests: config round-trip and lazy optional dependency checks.

### M2: FASTA and deterministic fake backend

- Generate canonical candidate FASTA from configured WT sequence and mutation policy.
- Add deterministic `fake` backend for CI.
- Normalize backend output to canonical embedding CSV.
- Snapshot FASTA/CSV/metadata into CAS.
- Tests: no torch/GPU/network required.

### M3: Workspace state integration

Status: implemented for the current workspace state machine. `evolvepro evolve` prepares missing embeddings through configured backends and has CI coverage for the deterministic `fake` backend. Existing valid `inputs/embeddings.csv` files are reused without extraction. Generated FASTA, embedding CSV, metadata, backend/device metadata, and embedding config hashes are recorded in `round.lock` and `manifest.json`.

Validation rejects generated FASTA/CSV/metadata mutation, generated embedding CAS mutation, and embedding config payload/hash mismatch.

### M4: ESM GPU backend

- Wrap current ESM extraction logic behind the backend interface.
- Implement GPU-first device resolution.
- Support explicit CPU fallback only when configured.
- Add smoke tests that mock/monkeypatch torch availability and backend calls; no real model download.

### M5: CLI polish and compatibility bridge

- Add targeted flags only if needed:
  - `evolvepro evolve --device cuda:0`
  - `evolvepro evolve --allow-cpu-embedding`
- Keep config as source of truth; flags are run-local overrides recorded in the lock.
- Update `scripts/plm/README.md` to point to target CLI workflow once parity exists.

## Acceptance criteria

- `evolvepro init --set-default --yes` writes an embedding config with GPU-first defaults.
- `evolvepro evolve` can prepare embeddings via fake backend in CI and continue training.
- A missing CUDA device fails clearly in GPU-first ESM mode.
- Existing precomputed embeddings are reused and validated by hash.
- All generated embedding artifacts are in CAS and manifest/lock payloads.
- Normal tests pass without torch, GPU, network, or model downloads.
