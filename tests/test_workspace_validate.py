from __future__ import annotations

import json
import os
import stat

import pandas as pd
import pytest

from evolvepro.cli import main
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


def test_validate_workspace_rejects_generated_embedding_mutation(tmp_path):
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

    (tmp_path / "rounds" / "round_001" / "embeddings" / "embeddings.csv").write_text(
        "variant,feature_0\nA1C,999\n"
    )

    with pytest.raises(ValueError, match="file hash mismatch"):
        validate_workspace(tmp_path)
