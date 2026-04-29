from __future__ import annotations

from evolvepro.embeddings.backends.base import EmbeddingBackend
from evolvepro.embeddings.backends.esm import EsmEmbeddingBackend
from evolvepro.embeddings.backends.fake import FakeEmbeddingBackend


def get_embedding_backend(name: str) -> EmbeddingBackend:
    normalized = name.strip().lower()
    if normalized == "fake":
        return FakeEmbeddingBackend()
    if normalized == "esm":
        return EsmEmbeddingBackend()
    raise ValueError(f"unsupported embedding backend: {name}")
