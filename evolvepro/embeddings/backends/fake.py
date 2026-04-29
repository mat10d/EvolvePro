from __future__ import annotations

import csv
import hashlib
import json
import struct
from pathlib import Path
from typing import TextIO

from evolvepro.embeddings.schemas import EmbeddingRequest, EmbeddingResult

_EMBEDDING_DIMENSIONS = 4
_HASH_INTEGER_RANGE = 1 << 32
_FLOAT_SCALE = 10_000_000_000


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seen_labels: set[str] = set()
    current_label: str | None = None
    current_sequence: list[str] = []

    def finish_record() -> None:
        nonlocal current_label, current_sequence
        if current_label is None:
            return
        sequence = "".join(current_sequence)
        if not sequence:
            raise ValueError(f"FASTA record '{current_label}' has an empty sequence")
        if current_label in seen_labels:
            raise ValueError(f"Duplicate FASTA label: {current_label}")
        seen_labels.add(current_label)
        records.append((current_label, sequence))
        current_label = None
        current_sequence = []

    with path.open("r", encoding="utf-8") as fasta_file:
        for line_number, raw_line in enumerate(fasta_file, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                label = line[1:].strip()
                if not label:
                    raise ValueError(f"FASTA label is empty at line {line_number}")
                current_label = label
                current_sequence = []
                continue
            if current_label is None:
                raise ValueError(
                    f"FASTA sequence data appears before the first label at line {line_number}"
                )
            current_sequence.append(line)

    finish_record()
    if not records:
        raise ValueError("FASTA file contains no records")
    return records


def _embedding_values(seed: int, label: str, sequence: str) -> list[str]:
    payload = (
        f"fake-embedding-v1\nseed={seed}\nlabel={label}\nsequence={sequence}\n"
    ).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    values = struct.unpack(">4I", digest[: _EMBEDDING_DIMENSIONS * 4])
    return [f"0.{value * _FLOAT_SCALE // _HASH_INTEGER_RANGE:010d}" for value in values]


def _write_embedding_csv(path: Path, records: list[tuple[str, str]], seed: int) -> None:
    with path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.writer(csv_file, lineterminator="\n")
        writer.writerow(
            ["variant", *[f"feature_{index}" for index in range(_EMBEDDING_DIMENSIONS)]]
        )
        for label, sequence in sorted(records, key=lambda record: record[0]):
            writer.writerow([label, *_embedding_values(seed, label, sequence)])


def _write_metadata_json(path: Path, metadata: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as metadata_file:
        _dump_stable_json(metadata, metadata_file)


def _dump_stable_json(data: dict[str, object], file: TextIO) -> None:
    json.dump(data, file, indent=2, sort_keys=True)
    file.write("\n")


class FakeEmbeddingBackend:
    def extract(self, request: EmbeddingRequest) -> EmbeddingResult:
        records = _parse_fasta(request.fasta_path)
        request.output_dir.mkdir(parents=True, exist_ok=True)

        embeddings_csv = request.output_dir / "embeddings.csv"
        metadata_json = request.output_dir / "metadata.json"
        metadata: dict[str, object] = {
            "backend": "fake",
            "device": "fake",
            "dtype": "float32",
            "model_name": request.model_name,
            "model_revision": request.model_revision,
            "repr_layers": request.repr_layers,
            "representation": request.representation,
            "seed": request.seed,
        }

        _write_embedding_csv(embeddings_csv, records, request.seed)
        _write_metadata_json(metadata_json, metadata)

        return EmbeddingResult(
            embeddings_csv=embeddings_csv,
            metadata_json=metadata_json,
            backend_metadata=metadata,
            device_metadata={"device": "fake", "dtype": "float32"},
        )
