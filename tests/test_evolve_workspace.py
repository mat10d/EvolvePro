from __future__ import annotations

import json
import stat

import pandas as pd

from evolvepro.cli import main
from evolvepro.io.hashing import sha256_file
from evolvepro.workspace import init_workspace


def test_evolve_initial_round_writes_candidate_handoff_and_manifest(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text().replace(
            "number_of_variants = 12", "number_of_variants = 3"
        )
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")

    assert main(["evolve"], cwd=tmp_path) == 0

    round_dir = tmp_path / "rounds" / "round_000"
    candidates = round_dir / "candidates.csv"
    assert candidates.read_text() == "variant,activity\nA1C,\nA1D,\nA1E,\n"
    assert (round_dir / "round.lock").exists()
    assert (round_dir / "manifest.json").exists()

    workspace_lock = json.loads((tmp_path / "evolvepro.lock").read_text())
    assert workspace_lock["current_round"] == "round_000"
    assert workspace_lock["rounds"] == [
        {
            "round": "round_000",
            "manifest": "rounds/round_000/manifest.json",
            "manifest_sha256": sha256_file(round_dir / "manifest.json"),
        }
    ]

    manifest = json.loads((round_dir / "manifest.json").read_text())
    snapshot_path = tmp_path / manifest["outputs"]["candidate_handoff"]["cas_path"]
    assert snapshot_path.read_bytes() == candidates.read_bytes()
    assert snapshot_path.stat().st_mode & stat.S_IWUSR == 0
    assert round_dir.joinpath("round.lock").stat().st_mode & stat.S_IWUSR == 0
    assert round_dir.joinpath("manifest.json").stat().st_mode & stat.S_IWUSR == 0


def test_evolve_stops_when_latest_candidates_have_no_activity(tmp_path, capsys):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text().replace(
            "number_of_variants = 12", "number_of_variants = 1"
        )
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")
    assert main(["evolve"], cwd=tmp_path) == 0

    assert main(["evolve"], cwd=tmp_path) == 1

    assert "fill the activity column" in capsys.readouterr().out


def test_evolve_ingests_partial_activity_and_writes_next_round(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text()
        .replace("number_of_variants = 12", "number_of_variants = 2")
        .replace('backend = "esm"', 'backend = "none"')
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")
    assert main(["evolve"], cwd=tmp_path) == 0
    round0 = tmp_path / "rounds" / "round_000"
    round0.joinpath("candidates.csv").write_text("variant,activity\nA1C,0.5\nA1D,\n")

    assert main(["evolve"], cwd=tmp_path) == 0

    round1 = tmp_path / "rounds" / "round_001"
    assert round1.joinpath("candidates.csv").read_text() == (
        "variant,activity\nA1E,\nA1F,\n"
    )
    manifest = json.loads(round1.joinpath("manifest.json").read_text())
    assert manifest["step"] == "activity_ingest_and_candidate_generation"
    assert manifest["observations"] == {
        "measured_count": 1,
        "unmeasured_count": 1,
    }
    assert manifest["previous_manifest_path"] == "rounds/round_000/manifest.json"
    assert manifest["previous_manifest_sha256"] == sha256_file(round0 / "manifest.json")
    assert "semantic_sha256" in manifest["inputs"]["activity_handoff"]
    activity_snapshot = tmp_path / manifest["inputs"]["activity_handoff"]["cas_path"]
    assert activity_snapshot.read_text() == "variant,activity\nA1C,0.5\nA1D,\n"

    workspace_lock = json.loads((tmp_path / "evolvepro.lock").read_text())
    assert workspace_lock["current_round"] == "round_001"
    assert [entry["round"] for entry in workspace_lock["rounds"]] == [
        "round_000",
        "round_001",
    ]


def test_evolve_prepares_fake_embeddings_when_missing(tmp_path):
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
    round0 = tmp_path / "rounds" / "round_000"
    round0.joinpath("candidates.csv").write_text("variant,activity\nA1C,0.0\nA1D,1.0\n")

    assert main(["evolve"], cwd=tmp_path) == 0

    round1 = tmp_path / "rounds" / "round_001"
    assert round1.joinpath("embeddings", "candidates.fasta").exists()
    assert round1.joinpath("embeddings", "embeddings.csv").exists()
    assert round1.joinpath("embeddings", "metadata.json").exists()
    assert round1.joinpath("model", "model.pkl").exists()
    assert round1.joinpath("training", "training_set.csv").exists()
    assert round1.joinpath("predictions", "predictions.parquet").exists()
    manifest = json.loads(round1.joinpath("manifest.json").read_text())
    assert "embedding_config" in manifest["inputs"]
    assert "semantic_sha256" in manifest["inputs"]["embedding_config"]
    embedding_outputs = manifest["outputs"]["embedding_preparation"]
    assert {"fasta", "embeddings", "metadata"} <= embedding_outputs.keys()
    assert (tmp_path / embedding_outputs["fasta"]["cas_path"]).exists()
    assert (tmp_path / embedding_outputs["embeddings"]["cas_path"]).exists()
    assert (tmp_path / embedding_outputs["metadata"]["cas_path"]).exists()


def test_evolve_uses_embeddings_and_model_when_available(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        config_path.read_text()
        .replace("number_of_variants = 12", "number_of_variants = 2")
        .replace('regression_type = "randomforest"', 'regression_type = "linear"')
        .replace("n_jobs = 1\n", "")
    )
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nAC\n")
    pd.DataFrame(
        {"feature": [0.0, 1.0, 0.9, 0.2, 0.8]},
        index=["A1C", "A1D", "A1E", "A1F", "A1G"],
    ).to_csv(tmp_path / "inputs" / "embeddings.csv", index_label="variant")

    assert main(["evolve"], cwd=tmp_path) == 0
    round0 = tmp_path / "rounds" / "round_000"
    round0.joinpath("candidates.csv").write_text("variant,activity\nA1C,0.0\nA1D,1.0\n")

    assert main(["evolve"], cwd=tmp_path) == 0

    round1 = tmp_path / "rounds" / "round_001"
    assert round1.joinpath("candidates.csv").read_text() == (
        "variant,activity\nA1E,\nA1G,\n"
    )
    assert round1.joinpath("model", "model.pkl").exists()
    assert round1.joinpath("training", "training_set.csv").exists()
    assert round1.joinpath("predictions", "predictions.parquet").exists()
    manifest = json.loads(round1.joinpath("manifest.json").read_text())
    assert "embeddings" in manifest["inputs"]
    assert "semantic_sha256" in manifest["inputs"]["embeddings"]
    assert "training_table" in manifest["outputs"]
    assert "predictions" in manifest["outputs"]
    assert "model_checkpoint" in manifest["outputs"]
