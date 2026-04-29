from __future__ import annotations

from pathlib import Path

from evolvepro.embeddings import EmbeddingRequest, EmbeddingResult
from evolvepro.embeddings.backends.base import EmbeddingBackend


class RecordingBackend:
    def extract(self, request: EmbeddingRequest) -> EmbeddingResult:
        return EmbeddingResult(
            embeddings_csv=request.output_dir / "embeddings.csv",
            metadata_json=request.output_dir / "metadata.json",
            backend_metadata={"backend": request.backend},
            device_metadata={"requested_device": request.device},
            artifacts={"fasta": request.fasta_path},
        )


def test_embedding_backend_protocol_scaffold():
    backend: EmbeddingBackend = RecordingBackend()
    request = EmbeddingRequest(
        fasta_path=Path("candidates.fasta"),
        output_dir=Path("embeddings"),
        backend="fake",
        model_name="fake-model",
        model_revision=None,
        repr_layers=[-1],
        representation="mean",
        toks_per_batch=512,
        device="auto",
        dtype="float32",
        seed=0,
        allow_download=False,
        allow_cpu_fallback=False,
        cache_dir=Path(".evolvepro/cache/plm"),
    )

    result = backend.extract(request)

    assert result.embeddings_csv == Path("embeddings/embeddings.csv")
    assert result.metadata_json == Path("embeddings/metadata.json")
    assert result.backend_metadata == {"backend": "fake"}
    assert result.device_metadata == {"requested_device": "auto"}
    assert result.artifacts == {"fasta": Path("candidates.fasta")}
