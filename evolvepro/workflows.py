from __future__ import annotations

import json
import stat
from dataclasses import asdict
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO

from evolvepro.io.cas import snapshot_file
from evolvepro.io.canonical import semantic_hash_table, write_json_atomic
from evolvepro.io.config import config_sha256, load_config
from evolvepro.io.handoff import read_candidate_handoff, write_candidate_handoff
from evolvepro.io.hashing import sha256_file
from evolvepro.io.schemas import FileHash
from evolvepro.src.model import top_layer
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


def _round_index(round_name: str) -> int:
    return int(round_name.rsplit("_", maxsplit=1)[1])


def _measured_observations(
    workspace: Path, workspace_lock: dict[str, object]
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for entry in workspace_lock["rounds"]:
        round_name = str(entry["round"])
        candidates_path = workspace / "rounds" / round_name / "candidates.csv"
        if not candidates_path.exists():
            continue
        measured = read_candidate_handoff(candidates_path).measured.copy()
        if measured.empty:
            continue
        measured["iteration"] = _round_index(round_name)
        frames.append(measured)
    if not frames:
        return pd.DataFrame(columns=["variant", "activity", "iteration"])
    observations = pd.concat(frames, ignore_index=True)
    duplicated = observations["variant"].duplicated(keep=False)
    if duplicated.any():
        duplicates = sorted(set(observations.loc[duplicated, "variant"].tolist()))
        raise ValueError(f"duplicate measured variants across rounds: {duplicates}")
    return observations


def _labels_for_embeddings(
    embeddings: pd.DataFrame, observations: pd.DataFrame
) -> pd.DataFrame:
    labels = pd.DataFrame({"variant": embeddings.index.astype(str)})
    labels = labels.merge(observations, on="variant", how="left")
    measured = labels["activity"].notna()
    labels["activity_binary"] = np.nan
    labels.loc[measured, "activity_binary"] = (
        labels.loc[measured, "activity"].astype(float) >= 1
    ).astype(int)
    labels["activity_scaled"] = np.nan
    if measured.any():
        activity = labels.loc[measured, "activity"].astype(float)
        min_activity = activity.min()
        max_activity = activity.max()
        if max_activity == min_activity:
            labels.loc[measured, "activity_scaled"] = 1.0
        else:
            labels.loc[measured, "activity_scaled"] = (activity - min_activity) / (
                max_activity - min_activity
            )
    return labels


def _write_table_artifact(
    path: Path, table: pd.DataFrame, workspace: Path, *, sort_by: list[str]
) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    snapshot = snapshot_file(path, workspace / WORKSPACE_DIR / "cas" / "sha256")
    return {
        **_artifact_entry(snapshot.path, snapshot.file_hash, workspace),
        "semantic_sha256": semantic_hash_table(table, sort_by=sort_by),
    }


def _rank_with_model(
    workspace: Path,
    round_dir: Path,
    workspace_lock: dict[str, object],
    exclude: set[str],
) -> tuple[list[str], dict[str, object]] | None:
    config = load_config(workspace / CONFIG_FILE)
    embeddings_path = (
        config.experiment.embeddings_base_path / config.experiment.embeddings_file_name
    )
    if not embeddings_path.is_file():
        return None

    embeddings = pd.read_csv(embeddings_path, index_col=0)
    observations = _measured_observations(workspace, workspace_lock)
    if observations.empty:
        return None

    labels = _labels_for_embeddings(embeddings, observations)
    train_iterations = sorted(
        labels.loc[labels["iteration"].notna(), "iteration"].unique().tolist()
    )
    _, df_test, df_sorted_all, metrics, model = top_layer(
        iter_train=train_iterations,
        iter_test=None,
        embeddings_pd=embeddings,
        labels_pd=labels,
        measured_var="activity",
        regression_type=config.model.regression_type,
        experimental=True,
        model_seed=config.model.seed,
        model_params=config.model.params,
    )
    ranked = df_test.sort_values(by=["y_pred"], ascending=False)
    variants = [
        variant
        for variant in ranked["variant"].astype(str).tolist()
        if variant not in exclude
    ][: config.experiment.number_of_variants]

    model_dir = round_dir / "model"
    model_dir.mkdir(parents=True, exist_ok=True)
    model_path = model_dir / "model.pkl"
    joblib.dump(model, model_path)
    params_path = model_dir / "params.json"
    write_json_atomic(params_path, model.get_params(deep=False))

    predictions_dir = round_dir / "predictions"
    predictions_dir.mkdir(parents=True, exist_ok=True)
    predictions_path = predictions_dir / "predictions.parquet"
    df_sorted_all.to_parquet(predictions_path, index=False)

    training_artifact = _write_table_artifact(
        round_dir / "training" / "training_set.csv",
        observations,
        workspace,
        sort_by=["variant"],
    )
    outputs = {
        "training_table": training_artifact,
        "model_checkpoint": asdict(_file_hash(model_path)),
        "model_params": asdict(_file_hash(params_path)),
        "predictions": {
            **asdict(_file_hash(predictions_path)),
            "semantic_sha256": semantic_hash_table(df_sorted_all, sort_by=["variant"]),
        },
        "metrics": metrics,
    }
    return variants, outputs


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

    used_variants = _used_variants(workspace, workspace_lock)
    model_result = _rank_with_model(workspace, round_dir, workspace_lock, used_variants)
    model_outputs: dict[str, object] = {}
    if model_result is None:
        wt_record = SeqIO.read(config.experiment.wt_fasta_path, "fasta")
        variants = _single_mutant_variants(
            str(wt_record.seq),
            config.experiment.number_of_variants,
            exclude=used_variants,
        )
    else:
        variants, model_outputs = model_result
    write_candidate_handoff(candidates_path, variants)

    activity_artifact = _handoff_artifact_entry(activity_path, workspace)
    candidate_artifact = _handoff_artifact_entry(candidates_path, workspace)
    outputs: dict[str, object] = {"candidate_handoff": candidate_artifact}
    outputs.update(model_outputs)
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
        "outputs": outputs,
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
