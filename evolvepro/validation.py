from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from evolvepro.io.canonical import canonical_json_sha256
from evolvepro.io.hashing import sha256_file
from evolvepro.workspace import LOCK_FILE, find_workspace


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _validate_cas_entry(workspace: Path, entry: dict[str, Any]) -> None:
    cas_path = workspace / entry["cas_path"]
    actual = sha256_file(cas_path)
    if actual != entry["sha256"]:
        raise ValueError(
            f"CAS hash mismatch for {cas_path}: expected {entry['sha256']}, got {actual}"
        )
    actual_size = cas_path.stat().st_size
    if actual_size != entry["size"]:
        raise ValueError(
            f"CAS size mismatch for {cas_path}: expected {entry['size']}, got {actual_size}"
        )


def _validate_file_entry(entry: dict[str, Any], workspace: Path) -> None:
    path = Path(entry["path"])
    if not path.is_absolute():
        path = workspace / path
    actual = sha256_file(path)
    if actual != entry["sha256"]:
        raise ValueError(
            f"file hash mismatch for {path}: expected {entry['sha256']}, got {actual}"
        )
    actual_size = path.stat().st_size
    if actual_size != entry["size"]:
        raise ValueError(
            f"file size mismatch for {path}: expected {entry['size']}, got {actual_size}"
        )


def _walk_hash_entries(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, dict):
        entries: list[dict[str, Any]] = []
        if "cas_path" in value and {"sha256", "size"} <= value.keys():
            entries.append(value)
        elif {"path", "sha256", "size"} <= value.keys():
            entries.append(value)
        for item in value.values():
            entries.extend(_walk_hash_entries(item))
        return entries
    if isinstance(value, list):
        entries = []
        for item in value:
            entries.extend(_walk_hash_entries(item))
        return entries
    return []


def _walk_semantic_payload_entries(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, dict):
        entries: list[dict[str, Any]] = []
        if {"payload", "semantic_sha256"} <= value.keys():
            entries.append(value)
        for item in value.values():
            entries.extend(_walk_semantic_payload_entries(item))
        return entries
    if isinstance(value, list):
        entries = []
        for item in value:
            entries.extend(_walk_semantic_payload_entries(item))
        return entries
    return []


def _validate_semantic_payload_entry(entry: dict[str, Any]) -> None:
    actual = canonical_json_sha256(entry["payload"])
    if actual != entry["semantic_sha256"]:
        raise ValueError(
            f"semantic hash mismatch: expected {entry['semantic_sha256']}, got {actual}"
        )


def _validate_manifest(workspace: Path, manifest_path: Path) -> None:
    manifest = _load_json(manifest_path)
    lock_path = manifest_path.parent / "round.lock"
    actual_lock_hash = sha256_file(lock_path)
    if manifest["lock_sha256"] != actual_lock_hash:
        raise ValueError(
            "round lock hash mismatch: "
            f"expected {manifest['lock_sha256']}, got {actual_lock_hash}"
        )

    for entry in _walk_hash_entries(manifest):
        if "cas_path" in entry:
            _validate_cas_entry(workspace, entry)
            if not entry.get("mutable", False):
                _validate_file_entry(entry, workspace)
        else:
            _validate_file_entry(entry, workspace)

    for entry in _walk_semantic_payload_entries(manifest):
        _validate_semantic_payload_entry(entry)

    previous_path = manifest.get("previous_manifest_path")
    previous_sha256 = manifest.get("previous_manifest_sha256")
    if previous_path is not None:
        actual = sha256_file(workspace / previous_path)
        if actual != previous_sha256:
            raise ValueError(
                "previous manifest hash mismatch: "
                f"expected {previous_sha256}, got {actual}"
            )


def validate_workspace(start: Path) -> None:
    workspace = find_workspace(start)
    if workspace is None:
        raise ValueError("not inside an EvolvePro workspace")

    lock_path = workspace / LOCK_FILE
    workspace_lock = _load_json(lock_path)
    for entry in workspace_lock["rounds"]:
        manifest_path = workspace / entry["manifest"]
        actual_manifest_hash = sha256_file(manifest_path)
        if actual_manifest_hash != entry["manifest_sha256"]:
            raise ValueError(
                "workspace manifest hash mismatch: "
                f"expected {entry['manifest_sha256']}, got {actual_manifest_hash}"
            )
        _validate_manifest(workspace, manifest_path)
