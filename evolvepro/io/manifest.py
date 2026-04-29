from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

from evolvepro.io.canonical import canonical_json_sha256, write_json_atomic
from evolvepro.io.config import config_sha256
from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import EvolveProConfig, FileHash, RoundLock


def manifest_sha256(path: Path) -> str:
    return sha256_file(path)


def _input_paths(config: EvolveProConfig) -> list[Path]:
    experiment = config.experiment
    return [
        experiment.embeddings_base_path / experiment.embeddings_file_name,
        *[experiment.round_base_path / name for name in experiment.round_file_names],
        experiment.wt_fasta_path,
    ]


def _file_hash(path: Path) -> FileHash:
    return FileHash(path=str(path), sha256=sha256_file(path), size=path.stat().st_size)


def build_round_lock(config: EvolveProConfig) -> RoundLock:
    previous_path = config.run.previous_manifest
    previous_sha256 = manifest_sha256(previous_path) if previous_path else None
    return RoundLock(
        schema_version=1,
        workflow=config.workflow,
        config_sha256=config_sha256(config),
        inputs=[_file_hash(path) for path in _input_paths(config)],
        model=config.model,
        previous_manifest_path=str(previous_path) if previous_path else None,
        previous_manifest_sha256=previous_sha256,
    )


def round_lock_payload(lock: RoundLock) -> dict[str, Any]:
    return asdict(lock)


def round_lock_sha256(lock: RoundLock) -> str:
    return canonical_json_sha256(round_lock_payload(lock))


def write_round_lock(path: Path, lock: RoundLock) -> None:
    write_json_atomic(path, round_lock_payload(lock))


def validate_previous_manifest(expected_sha256: str, path: Path) -> None:
    actual = manifest_sha256(path)
    if actual != expected_sha256:
        raise ValueError(
            f"previous manifest hash mismatch: expected {expected_sha256}, got {actual}"
        )


def validate_round_lock(lock: RoundLock) -> None:
    for entry in lock.inputs:
        path = Path(entry.path)
        actual = sha256_file(path)
        if actual != entry.sha256:
            raise ValueError(
                f"input hash mismatch for {entry.path}: expected {entry.sha256}, got {actual}"
            )
        actual_size = path.stat().st_size
        if actual_size != entry.size:
            raise ValueError(
                f"input size mismatch for {entry.path}: expected {entry.size}, got {actual_size}"
            )
    if lock.previous_manifest_path is not None:
        if lock.previous_manifest_sha256 is None:
            raise ValueError("previous manifest path set without hash")
        validate_previous_manifest(
            lock.previous_manifest_sha256, Path(lock.previous_manifest_path)
        )
