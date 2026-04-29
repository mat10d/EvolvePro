from __future__ import annotations

from typing import Protocol

from evolvepro.embeddings.schemas import EmbeddingRequest, EmbeddingResult


class EmbeddingBackend(Protocol):
    def extract(self, request: EmbeddingRequest) -> EmbeddingResult: ...
