from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from evolvepro.embeddings.backends.registry import get_embedding_backend
from evolvepro.embeddings.fasta import read_single_fasta_sequence, write_embedding_fasta
from evolvepro.embeddings.normalize import snapshot_embedding_csv, write_embedding_csv
from evolvepro.embeddings.schemas import EmbeddingRequest
from evolvepro.io.cas import snapshot_file
from evolvepro.io.canonical import canonical_json_sha256
from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import EvolveProConfig
from evolvepro.workspace import WORKSPACE_DIR


@dataclass(frozen=True)
class PreparedEmbeddings:
    embeddings_csv: Path
    inputs: dict[str, object]
    outputs: dict[str, object]


def embedding_config_payload(config: EvolveProConfig) -> dict[str, Any]:
    return asdict(config.embedding)


def embedding_config_sha256(config: EvolveProConfig) -> str:
    return canonical_json_sha256(embedding_config_payload(config))


def _relative_path(path: Path, workspace: Path) -> str:
    try:
        return path.resolve().relative_to(workspace.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _snapshot_artifact(path: Path, workspace: Path) -> dict[str, object]:
    snapshot = snapshot_file(path, workspace / WORKSPACE_DIR / "cas" / "sha256")
    return {
        "path": _relative_path(path, workspace),
        "sha256": sha256_file(path),
        "size": path.stat().st_size,
        "cas_path": _relative_path(snapshot.path, workspace),
    }


def prepare_round_embeddings(
    workspace: Path,
    round_dir: Path,
    config: EvolveProConfig,
    variants: list[str],
) -> PreparedEmbeddings:
    if config.embedding.backend == "none":
        raise ValueError("embedding backend 'none' cannot prepare embeddings")

    output_dir = round_dir / "embeddings"
    fasta_path = output_dir / "candidates.fasta"
    _, wt_sequence = read_single_fasta_sequence(config.experiment.wt_fasta_path)
    write_embedding_fasta(fasta_path, wt_sequence, variants)

    request = EmbeddingRequest(
        fasta_path=fasta_path,
        output_dir=output_dir,
        backend=config.embedding.backend,
        model_name=config.embedding.model_name,
        model_revision=config.embedding.model_revision,
        repr_layers=config.embedding.repr_layers,
        representation=config.embedding.representation,
        toks_per_batch=config.embedding.toks_per_batch,
        device=config.embedding.device,
        dtype=config.embedding.dtype,
        seed=config.run.seed,
        allow_download=config.embedding.allow_download,
        allow_cpu_fallback=config.embedding.allow_cpu_fallback,
        cache_dir=config.embedding.cache_dir,
    )
    result = get_embedding_backend(config.embedding.backend).extract(request)
    write_embedding_csv(result.embeddings_csv, pd.read_csv(result.embeddings_csv))

    inputs: dict[str, object] = {
        "embedding_config": {
            "semantic_sha256": embedding_config_sha256(config),
            "payload": embedding_config_payload(config),
        }
    }
    outputs: dict[str, object] = {
        "embedding_preparation": {
            "fasta": _snapshot_artifact(fasta_path, workspace),
            "embeddings": snapshot_embedding_csv(result.embeddings_csv, workspace),
            "metadata": _snapshot_artifact(result.metadata_json, workspace),
            "backend_metadata": result.backend_metadata,
            "device_metadata": result.device_metadata,
        }
    }
    return PreparedEmbeddings(
        embeddings_csv=result.embeddings_csv,
        inputs=inputs,
        outputs=outputs,
    )
