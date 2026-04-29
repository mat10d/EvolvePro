from __future__ import annotations

import json
import stat
from dataclasses import asdict
from pathlib import Path

from Bio import SeqIO

from evolvepro.io.cas import snapshot_file
from evolvepro.io.canonical import semantic_hash_table, write_json_atomic
from evolvepro.io.config import config_sha256, load_config
from evolvepro.io.handoff import read_candidate_handoff, write_candidate_handoff
from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import FileHash
from evolvepro.workspace import CONFIG_FILE, LOCK_FILE, WORKSPACE_DIR, find_workspace

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


class UserActionRequired(RuntimeError):
    pass


def _relative_to_workspace(path: Path, workspace: Path) -> str:
    return path.relative_to(workspace).as_posix()


def _make_read_only(path: Path) -> None:
    path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


def _write_read_only_json(path: Path, payload: dict[str, object]) -> None:
    write_json_atomic(path, payload)
    _make_read_only(path)


def _file_hash(path: Path) -> FileHash:
    return FileHash(path=str(path), sha256=sha256_file(path), size=path.stat().st_size)


def _artifact_entry(
    snapshot_path: Path, file_hash: FileHash, workspace: Path
) -> dict[str, object]:
    return {
        **asdict(file_hash),
        "cas_path": _relative_to_workspace(snapshot_path, workspace),
    }


def _handoff_artifact_entry(path: Path, workspace: Path) -> dict[str, object]:
    snapshot = snapshot_file(path, workspace / WORKSPACE_DIR / "cas" / "sha256")
    handoff = read_candidate_handoff(path)
    return {
        **_artifact_entry(snapshot.path, snapshot.file_hash, workspace),
        "semantic_sha256": semantic_hash_table(handoff.table, sort_by=["variant"]),
    }


def _single_mutant_variants(
    sequence: str, limit: int, *, exclude: set[str] | None = None
) -> list[str]:
    excluded = exclude or set()
    variants: list[str] = []
    for index, wt_aa in enumerate(sequence, start=1):
        for mutant_aa in AA_ALPHABET:
            if mutant_aa == wt_aa:
                continue
            variant = f"{wt_aa}{index}{mutant_aa}"
            if variant in excluded:
                continue
            variants.append(variant)
            if len(variants) == limit:
                return variants
    return variants


def _load_workspace_lock(workspace: Path) -> dict[str, object]:
    return json.loads((workspace / LOCK_FILE).read_text())


def _write_workspace_lock(workspace: Path, payload: dict[str, object]) -> None:
    write_json_atomic(workspace / LOCK_FILE, payload)


def _next_round_name(round_name: str) -> str:
    prefix, _, suffix = round_name.rpartition("_")
    if prefix != "round" or not suffix.isdigit():
        raise ValueError(f"unsupported round name in workspace lock: {round_name}")
    return f"round_{int(suffix) + 1:03d}"


def _workspace_round_entry(
    workspace: Path, workspace_lock: dict[str, object], round_name: str
) -> dict[str, object]:
    for entry in workspace_lock["rounds"]:
        if entry["round"] == round_name:
            return entry
    raise ValueError(f"round missing from workspace lock: {round_name}")


def _validate_manifest_entry(workspace: Path, entry: dict[str, object]) -> None:
    manifest_path = workspace / str(entry["manifest"])
    actual = sha256_file(manifest_path)
    if actual != entry["manifest_sha256"]:
        raise ValueError(
            "workspace manifest hash mismatch: "
            f"expected {entry['manifest_sha256']}, got {actual}"
        )


def _used_variants(workspace: Path, workspace_lock: dict[str, object]) -> set[str]:
    variants: set[str] = set()
    for entry in workspace_lock["rounds"]:
        candidates_path = workspace / "rounds" / str(entry["round"]) / "candidates.csv"
        if candidates_path.exists():
            variants.update(
                read_candidate_handoff(candidates_path).table["variant"].tolist()
            )
    return variants


def _initial_candidate_round(workspace: Path) -> Path:
    config = load_config(workspace / CONFIG_FILE)
    round_name = "round_000"
    round_dir = workspace / "rounds" / round_name
    candidates_path = round_dir / "candidates.csv"

    wt_record = SeqIO.read(config.experiment.wt_fasta_path, "fasta")
    variants = _single_mutant_variants(
        str(wt_record.seq), config.experiment.number_of_variants
    )
    write_candidate_handoff(candidates_path, variants)

    candidate_artifact = _handoff_artifact_entry(candidates_path, workspace)
    wt_hash = _file_hash(config.experiment.wt_fasta_path)

    lock_path = round_dir / "round.lock"
    lock_payload: dict[str, object] = {
        "schema_version": 1,
        "workflow": config.workflow,
        "step": "initial_candidate_generation",
        "round": round_name,
        "config_sha256": config_sha256(config),
        "inputs": {"wt_fasta": asdict(wt_hash)},
        "outputs": {"candidate_handoff": candidate_artifact},
    }
    _write_read_only_json(lock_path, lock_payload)

    manifest_path = round_dir / "manifest.json"
    manifest_payload: dict[str, object] = {
        "schema_version": 1,
        "workflow": config.workflow,
        "step": "initial_candidate_generation",
        "round": round_name,
        "lock_sha256": sha256_file(lock_path),
        "inputs": lock_payload["inputs"],
        "outputs": lock_payload["outputs"],
        "previous_manifest_path": None,
        "previous_manifest_sha256": None,
    }
    _write_read_only_json(manifest_path, manifest_payload)

    workspace_lock = _load_workspace_lock(workspace)
    workspace_lock["current_round"] = round_name
    workspace_lock["rounds"] = [
        {
            "round": round_name,
            "manifest": _relative_to_workspace(manifest_path, workspace),
            "manifest_sha256": sha256_file(manifest_path),
        }
    ]
    _write_workspace_lock(workspace, workspace_lock)
    return candidates_path


def _advance_from_activity(workspace: Path, current_round: str) -> Path:
    config = load_config(workspace / CONFIG_FILE)
    workspace_lock = _load_workspace_lock(workspace)
    previous_entry = _workspace_round_entry(workspace, workspace_lock, current_round)
    _validate_manifest_entry(workspace, previous_entry)

    activity_path = workspace / "rounds" / current_round / "candidates.csv"
    handoff = read_candidate_handoff(activity_path)
    if handoff.measured.empty:
        raise UserActionRequired(
            f"fill the activity column in {activity_path} before running evolve again"
        )

    round_name = _next_round_name(current_round)
    round_dir = workspace / "rounds" / round_name
    candidates_path = round_dir / "candidates.csv"

    wt_record = SeqIO.read(config.experiment.wt_fasta_path, "fasta")
    variants = _single_mutant_variants(
        str(wt_record.seq),
        config.experiment.number_of_variants,
        exclude=_used_variants(workspace, workspace_lock),
    )
    write_candidate_handoff(candidates_path, variants)

    activity_artifact = _handoff_artifact_entry(activity_path, workspace)
    candidate_artifact = _handoff_artifact_entry(candidates_path, workspace)
    previous_manifest_path = workspace / str(previous_entry["manifest"])

    lock_path = round_dir / "round.lock"
    lock_payload: dict[str, object] = {
        "schema_version": 1,
        "workflow": config.workflow,
        "step": "activity_ingest_and_candidate_generation",
        "round": round_name,
        "config_sha256": config_sha256(config),
        "inputs": {
            "activity_handoff": activity_artifact,
            "previous_manifest": {
                "path": str(previous_entry["manifest"]),
                "sha256": previous_entry["manifest_sha256"],
            },
        },
        "outputs": {"candidate_handoff": candidate_artifact},
        "observations": {
            "measured_count": int(len(handoff.measured)),
            "unmeasured_count": int(len(handoff.unmeasured)),
        },
    }
    _write_read_only_json(lock_path, lock_payload)

    manifest_path = round_dir / "manifest.json"
    manifest_payload: dict[str, object] = {
        "schema_version": 1,
        "workflow": config.workflow,
        "step": "activity_ingest_and_candidate_generation",
        "round": round_name,
        "lock_sha256": sha256_file(lock_path),
        "inputs": lock_payload["inputs"],
        "outputs": lock_payload["outputs"],
        "observations": lock_payload["observations"],
        "previous_manifest_path": str(previous_entry["manifest"]),
        "previous_manifest_sha256": sha256_file(previous_manifest_path),
    }
    _write_read_only_json(manifest_path, manifest_payload)

    workspace_lock["current_round"] = round_name
    workspace_lock["rounds"].append(
        {
            "round": round_name,
            "manifest": _relative_to_workspace(manifest_path, workspace),
            "manifest_sha256": sha256_file(manifest_path),
        }
    )
    _write_workspace_lock(workspace, workspace_lock)
    return candidates_path


def evolve_workspace(start: Path) -> Path:
    workspace = find_workspace(start)
    if workspace is None:
        raise ValueError(
            "not inside an EvolvePro workspace; run `evolvepro init` first"
        )

    workspace_lock = _load_workspace_lock(workspace)
    current_round = workspace_lock.get("current_round")
    if current_round is None:
        return _initial_candidate_round(workspace)

    return _advance_from_activity(workspace, str(current_round))
