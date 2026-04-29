from __future__ import annotations

import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import joblib
import pandas as pd

from evolvepro.io.canonical import (
    canonical_json_sha256,
    semantic_hash_table,
    write_json_atomic,
)
from evolvepro.io.config import config_sha256
from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import EvolveProConfig, FileHash, RoundLock, RoundManifest


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


LEGACY_NAMES = [
    "iteration.csv",
    "this_round_variants.csv",
    "df_test.csv",
    "df_sorted_all.csv",
]
DUMP_NAMES = ["metrics.json", "inputs.json", "version.txt"]


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, str | int | float | bool):
        return value
    if hasattr(value, "item"):
        return _jsonable(value.item())
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_jsonable(item) for item in value]
    return repr(value)


def _copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def materialize_round_layout(round_dir: Path) -> None:
    for name in LEGACY_NAMES:
        source = round_dir / name
        if source.exists():
            _copy_file(source, round_dir / "legacy" / name)

    _copy_file(round_dir / "model.pkl", round_dir / "model" / "model.pkl")
    model = joblib.load(round_dir / "model.pkl")
    write_json_atomic(
        round_dir / "model" / "params.json",
        _jsonable(model.get_params(deep=False)),
    )

    _copy_file(
        round_dir / "predictions.parquet",
        round_dir / "predictions" / "predictions.parquet",
    )
    _copy_file(round_dir / "df_test.csv", round_dir / "predictions" / "candidates.csv")


def build_round_manifest(
    lock: RoundLock, *, lock_path: Path, round_dir: Path
) -> RoundManifest:
    predictions = pd.read_parquet(round_dir / "predictions.parquet")
    outputs: dict[str, object] = {
        "legacy_files": [
            asdict(_file_hash(round_dir / name))
            for name in LEGACY_NAMES
            if (round_dir / name).exists()
        ],
        "dump_files": [
            asdict(_file_hash(round_dir / name))
            for name in DUMP_NAMES
            if (round_dir / name).exists()
        ],
        "model_checkpoint": asdict(_file_hash(round_dir / "model.pkl")),
        "predictions_file": asdict(_file_hash(round_dir / "predictions.parquet")),
        "predictions_semantic_sha256": semantic_hash_table(
            predictions, sort_by=["variant"]
        ),
    }
    if (round_dir / "model" / "params.json").exists():
        outputs["planned_layout"] = {
            "legacy_files": [
                asdict(_file_hash(round_dir / "legacy" / name))
                for name in LEGACY_NAMES
                if (round_dir / "legacy" / name).exists()
            ],
            "model_checkpoint": asdict(_file_hash(round_dir / "model" / "model.pkl")),
            "model_params": asdict(_file_hash(round_dir / "model" / "params.json")),
            "predictions_file": asdict(
                _file_hash(round_dir / "predictions" / "predictions.parquet")
            ),
            "candidates_file": asdict(
                _file_hash(round_dir / "predictions" / "candidates.csv")
            ),
        }
    return RoundManifest(
        schema_version=1,
        workflow=lock.workflow,
        lock_sha256=sha256_file(lock_path),
        inputs=lock.inputs,
        outputs=outputs,
        previous_manifest_path=lock.previous_manifest_path,
        previous_manifest_sha256=lock.previous_manifest_sha256,
    )


def round_manifest_payload(manifest: RoundManifest) -> dict[str, Any]:
    return asdict(manifest)


def write_round_manifest(path: Path, manifest: RoundManifest) -> None:
    write_json_atomic(path, round_manifest_payload(manifest))


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
