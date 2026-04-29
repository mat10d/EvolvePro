from __future__ import annotations

import json
import os
import stat

import pandas as pd
import pytest

from evolvepro.cli import main
from evolvepro.io.hashing import sha256_file
from evolvepro.validation import validate_workspace
from evolvepro.workspace import init_workspace


def _init_small_workspace(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text()
        .replace("number_of_variants = 12", "number_of_variants = 2")
        .replace('backend = "esm"', 'backend = "none"')
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")


def test_validate_workspace_accepts_mutated_handoff_when_cas_is_intact(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    round0_candidates = tmp_path / "rounds" / "round_000" / "candidates.csv"
    round0_candidates.write_text("variant,activity\nA1C,0.5\nA1D,\n")
    assert main(["evolve"], cwd=tmp_path) == 0

    validate_workspace(tmp_path)
    assert main(["validate"], cwd=tmp_path) == 0


def test_validate_workspace_rejects_cas_mutation(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    manifest = json.loads(
        (tmp_path / "rounds" / "round_000" / "manifest.json").read_text()
    )
    cas_path = tmp_path / manifest["outputs"]["candidate_handoff"]["cas_path"]
    os.chmod(cas_path, stat.S_IRUSR | stat.S_IWUSR)
    cas_path.write_text("variant,activity\nBROKEN,\n")

    with pytest.raises(ValueError, match="CAS hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_manifest_mutation(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    manifest_path = tmp_path / "rounds" / "round_000" / "manifest.json"
    os.chmod(manifest_path, stat.S_IRUSR | stat.S_IWUSR)
    manifest = json.loads(manifest_path.read_text())
    manifest["step"] = "broken"
    manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n")

    with pytest.raises(ValueError, match="workspace manifest hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_round_lock_mutation(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    lock_path = tmp_path / "rounds" / "round_000" / "round.lock"
    os.chmod(lock_path, stat.S_IRUSR | stat.S_IWUSR)
    lock = json.loads(lock_path.read_text())
    lock["step"] = "broken"
    lock_path.write_text(json.dumps(lock, sort_keys=True) + "\n")

    with pytest.raises(ValueError, match="round lock hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_previous_manifest_chain_break(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    round0_manifest_path = tmp_path / "rounds" / "round_000" / "manifest.json"
    (tmp_path / "rounds" / "round_000" / "candidates.csv").write_text(
        "variant,activity\nA1C,0.0\nA1D,1.0\n"
    )
    assert main(["evolve"], cwd=tmp_path) == 0

    os.chmod(round0_manifest_path, stat.S_IRUSR | stat.S_IWUSR)
    round0_manifest = json.loads(round0_manifest_path.read_text())
    round0_manifest["step"] = "rewritten"
    round0_manifest_path.write_text(json.dumps(round0_manifest, sort_keys=True) + "\n")

    workspace_lock_path = tmp_path / "evolvepro.lock"
    workspace_lock = json.loads(workspace_lock_path.read_text())
    workspace_lock["rounds"][0]["manifest_sha256"] = sha256_file(round0_manifest_path)
    workspace_lock_path.write_text(json.dumps(workspace_lock, sort_keys=True) + "\n")

    with pytest.raises(ValueError, match="previous manifest hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_observation_artifact_mutation(tmp_path):
    _init_small_workspace(tmp_path)
    assert main(["evolve"], cwd=tmp_path) == 0
    (tmp_path / "rounds" / "round_000" / "candidates.csv").write_text(
        "variant,activity\nA1C,0.0\nA1D,1.0\n"
    )
    assert main(["evolve"], cwd=tmp_path) == 0

    (
        tmp_path / "rounds" / "round_001" / "observations" / "activity_observations.csv"
    ).write_text("variant,activity\nA1C,999\n")

    with pytest.raises(ValueError, match="file hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_embedding_input_mutation(tmp_path):
    _init_small_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text()
        .replace('regression_type = "randomforest"', 'regression_type = "linear"')
        .replace("n_jobs = 1\n", "")
    )
    embeddings_path = tmp_path / "inputs" / "embeddings.csv"
    pd.DataFrame({"feature": [0.0, 1.0, 0.9]}, index=["A1C", "A1D", "A1E"]).to_csv(
        embeddings_path, index_label="variant"
    )
    assert main(["evolve"], cwd=tmp_path) == 0
    (tmp_path / "rounds" / "round_000" / "candidates.csv").write_text(
        "variant,activity\nA1C,0.0\nA1D,1.0\n"
    )
    assert main(["evolve"], cwd=tmp_path) == 0

    embeddings_path.write_text("variant,feature\nA1C,999\nA1D,1\nA1E,0.9\n")

    with pytest.raises(ValueError, match="file hash mismatch"):
        validate_workspace(tmp_path)


def _workspace_with_generated_embeddings(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text()
        .replace("number_of_variants = 12", "number_of_variants = 2")
        .replace('regression_type = "randomforest"', 'regression_type = "linear"')
        .replace("n_jobs = 1\n", "")
        .replace('backend = "esm"', 'backend = "fake"')
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")
    assert main(["evolve"], cwd=tmp_path) == 0
    (tmp_path / "rounds" / "round_000" / "candidates.csv").write_text(
        "variant,activity\nA1C,0.0\nA1D,1.0\n"
    )
    assert main(["evolve"], cwd=tmp_path) == 0
    return tmp_path / "rounds" / "round_001"


@pytest.mark.parametrize(
    ("relative_path", "replacement"),
    [
        ("embeddings/candidates.fasta", ">BROKEN\nAC\n"),
        ("embeddings/embeddings.csv", "variant,feature_0\nA1C,999\n"),
        ("embeddings/metadata.json", '{"backend":"broken"}\n'),
    ],
)
def test_validate_workspace_rejects_generated_embedding_mutation(
    tmp_path, relative_path, replacement
):
    round1 = _workspace_with_generated_embeddings(tmp_path)

    round1.joinpath(relative_path).write_text(replacement)

    with pytest.raises(ValueError, match="file hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_training_artifact_mutation(tmp_path):
    round1 = _workspace_with_generated_embeddings(tmp_path)
    round1.joinpath("training", "training_set.csv").write_text(
        "variant,activity,iteration\nA1C,999,0\n"
    )

    with pytest.raises(ValueError, match="file hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_generated_embedding_cas_mutation(tmp_path):
    round1 = _workspace_with_generated_embeddings(tmp_path)
    manifest = json.loads(round1.joinpath("manifest.json").read_text())
    cas_path = (
        tmp_path
        / manifest["outputs"]["embedding_preparation"]["embeddings"]["cas_path"]
    )
    os.chmod(cas_path, stat.S_IRUSR | stat.S_IWUSR)
    cas_path.write_text("variant,feature_0\nA1C,999\n")

    with pytest.raises(ValueError, match="CAS hash mismatch"):
        validate_workspace(tmp_path)


def test_validate_workspace_rejects_embedding_config_payload_mismatch(tmp_path):
    round1 = _workspace_with_generated_embeddings(tmp_path)
    manifest_path = round1 / "manifest.json"
    os.chmod(manifest_path, stat.S_IRUSR | stat.S_IWUSR)
    manifest = json.loads(manifest_path.read_text())
    manifest["inputs"]["embedding_config"]["payload"]["backend"] = "esm"
    manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n")

    lock_path = tmp_path / "evolvepro.lock"
    workspace_lock = json.loads(lock_path.read_text())
    workspace_lock["rounds"][-1]["manifest_sha256"] = sha256_file(manifest_path)
    lock_path.write_text(json.dumps(workspace_lock, sort_keys=True) + "\n")

    with pytest.raises(ValueError, match="semantic hash mismatch"):
        validate_workspace(tmp_path)
