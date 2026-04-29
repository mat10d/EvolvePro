from __future__ import annotations

import json
import stat
from dataclasses import asdict
from pathlib import Path

from Bio import SeqIO

from evolvepro.io.cas import snapshot_file
from evolvepro.io.canonical import write_json_atomic
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


def _single_mutant_variants(sequence: str, limit: int) -> list[str]:
    variants: list[str] = []
    for index, wt_aa in enumerate(sequence, start=1):
        for mutant_aa in AA_ALPHABET:
            if mutant_aa == wt_aa:
                continue
            variants.append(f"{wt_aa}{index}{mutant_aa}")
            if len(variants) == limit:
                return variants
    return variants


def _load_workspace_lock(workspace: Path) -> dict[str, object]:
    return json.loads((workspace / LOCK_FILE).read_text())


def _write_workspace_lock(workspace: Path, payload: dict[str, object]) -> None:
    write_json_atomic(workspace / LOCK_FILE, payload)


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

    candidate_snapshot = snapshot_file(
        candidates_path, workspace / WORKSPACE_DIR / "cas" / "sha256"
    )
    wt_hash = _file_hash(config.experiment.wt_fasta_path)

    lock_path = round_dir / "round.lock"
    lock_payload: dict[str, object] = {
        "schema_version": 1,
        "workflow": config.workflow,
        "step": "initial_candidate_generation",
        "round": round_name,
        "config_sha256": config_sha256(config),
        "inputs": {"wt_fasta": asdict(wt_hash)},
        "outputs": {
            "candidate_handoff": {
                **asdict(candidate_snapshot.file_hash),
                "cas_path": _relative_to_workspace(candidate_snapshot.path, workspace),
            }
        },
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

    candidates_path = workspace / "rounds" / str(current_round) / "candidates.csv"
    handoff = read_candidate_handoff(candidates_path)
    if handoff.measured.empty:
        raise UserActionRequired(
            f"fill the activity column in {candidates_path} before running evolve again"
        )
    raise UserActionRequired("measured activity ingest is not implemented yet")
