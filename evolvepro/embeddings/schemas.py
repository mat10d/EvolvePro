from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


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
    allow_download: bool = False
    allow_cpu_fallback: bool = False
    cache_dir: Path | None = None


@dataclass(frozen=True)
class EmbeddingResult:
    embeddings_csv: Path
    metadata_json: Path
    backend_metadata: dict[str, Any] = field(default_factory=dict)
    device_metadata: dict[str, Any] = field(default_factory=dict)
    artifacts: dict[str, Path] = field(default_factory=dict)
