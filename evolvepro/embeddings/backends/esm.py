from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
from typing import Any

from evolvepro.embeddings.schemas import EmbeddingRequest, EmbeddingResult


_CANONICAL_EMBEDDINGS_CSV = "embeddings.csv"
_METADATA_JSON = "metadata.json"
_LEGACY_TRUNCATION_SEQ_LENGTH = 1022


class EsmEmbeddingBackend:
    def extract(self, request: EmbeddingRequest) -> EmbeddingResult:
        output_dir = request.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        pt_output_dir = output_dir / "pt"
        embeddings_csv = output_dir / _CANONICAL_EMBEDDINGS_CSV
        metadata_json = output_dir / _METADATA_JSON

        args = _build_legacy_args(request, pt_output_dir)
        _run_legacy_esm(args)
        if not embeddings_csv.exists():
            _concatenate_legacy_outputs(pt_output_dir, embeddings_csv)
        if not embeddings_csv.exists():
            raise FileNotFoundError(
                f"ESM extraction did not produce canonical CSV: {embeddings_csv}"
            )

        backend_metadata = _backend_metadata(request)
        device_metadata = _device_metadata(request, args.nogpu)
        _write_metadata(
            metadata_json=metadata_json,
            request=request,
            backend_metadata=backend_metadata,
            device_metadata=device_metadata,
            artifacts={
                "pt_dir": pt_output_dir,
                "embeddings_csv": embeddings_csv,
            },
        )

        return EmbeddingResult(
            embeddings_csv=embeddings_csv,
            metadata_json=metadata_json,
            backend_metadata=backend_metadata,
            device_metadata=device_metadata,
            artifacts={"pt_dir": pt_output_dir},
        )


def _build_legacy_args(request: EmbeddingRequest, pt_output_dir: Path) -> Namespace:
    return Namespace(
        model_location=request.model_name,
        fasta_file=request.fasta_path,
        output_dir=pt_output_dir,
        toks_per_batch=request.toks_per_batch,
        repr_layers=request.repr_layers,
        include=[request.representation],
        truncation_seq_length=_LEGACY_TRUNCATION_SEQ_LENGTH,
        nogpu=_nogpu_for_device(request.device),
        concatenate_dir=request.output_dir,
    )


def _nogpu_for_device(device: str) -> bool:
    normalized = device.lower()
    if normalized == "cpu":
        return True
    if normalized == "auto" or normalized.startswith("cuda"):
        return False
    return False


def _run_legacy_esm(args: object) -> None:
    from evolvepro.plm.esm import extract as legacy_extract

    legacy_extract.run(args)


def _concatenate_legacy_outputs(output_dir: Path, output_csv: Path) -> None:
    from evolvepro.plm.esm import extract as legacy_extract

    legacy_extract.concatenate_files(output_dir, output_csv)


def _backend_metadata(request: EmbeddingRequest) -> dict[str, Any]:
    return {
        "backend": request.backend,
        "legacy_wrapper": "evolvepro.plm.esm.extract",
        "model_name": request.model_name,
        "model_revision": request.model_revision,
        "representation": request.representation,
        "repr_layers": request.repr_layers,
    }


def _device_metadata(request: EmbeddingRequest, nogpu: bool) -> dict[str, Any]:
    return {
        "requested_device": request.device,
        "nogpu": nogpu,
        "dtype": request.dtype,
        "seed": request.seed,
    }


def _write_metadata(
    *,
    metadata_json: Path,
    request: EmbeddingRequest,
    backend_metadata: dict[str, Any],
    device_metadata: dict[str, Any],
    artifacts: dict[str, Path],
) -> None:
    payload = {
        "backend": backend_metadata,
        "model": {
            "name": request.model_name,
            "revision": request.model_revision,
        },
        "representation": {
            "name": request.representation,
            "repr_layers": request.repr_layers,
        },
        "device": device_metadata,
        "dtype": request.dtype,
        "seed": request.seed,
        "artifacts": {key: str(value) for key, value in sorted(artifacts.items())},
    }
    metadata_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
