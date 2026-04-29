from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest

from evolvepro.io.config import (
    config_sha256,
    evolve_experimental_from_config,
    load_config,
    write_config,
)
from evolvepro.io.hashing import sha256_file
from evolvepro.src.data import load_experimental_data
from evolvepro.src.evolve import evolve_experimental

ROUND_BASE = Path("colab") / "rounds_data"
LEGACY_OUTPUTS = [
    "iteration.csv",
    "this_round_variants.csv",
    "df_test.csv",
    "df_sorted_all.csv",
]
DUMP_OUTPUTS = [
    "model.pkl",
    "predictions.parquet",
    "metrics.json",
    "inputs.json",
    "version.txt",
]


def _make_embeddings(tmp_path: Path, repo_root: Path, kelsic_wt_fasta: Path) -> Path:
    round1 = load_experimental_data(
        str(repo_root / ROUND_BASE), "kelsic_Round1.xlsx", str(kelsic_wt_fasta)
    )
    variants = sorted(
        round1["updated_variant"].tolist() + [f"UNMEASURED_{i}" for i in range(24)]
    )
    rng = np.random.default_rng(123)
    embeddings = pd.DataFrame(
        rng.normal(size=(len(variants), 10)),
        index=variants,
        columns=[f"embedding_{i}" for i in range(10)],
    )
    embeddings_file = tmp_path / "kelsic_synthetic.csv"
    embeddings.to_csv(embeddings_file)
    return embeddings_file


def _write_round_config(
    tmp_path: Path, repo_root: Path, kelsic_wt_fasta: Path, output_dir: Path
) -> Path:
    embeddings_file = _make_embeddings(tmp_path, repo_root, kelsic_wt_fasta)
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        f"""
schema_version = 1
workflow = "experimental_round"

[run]
seed = 0
output_dir = "{output_dir.as_posix()}"
overwrite = false

[experiment]
protein_name = "kelsic"
round_name = "Round1"
rename_WT = false
number_of_variants = 12

[inputs]
embeddings_base_path = "{tmp_path.as_posix()}"
embeddings_file_name = "{embeddings_file.name}"
round_base_path = "{(repo_root / ROUND_BASE).as_posix()}"
round_file_names = ["kelsic_Round1.xlsx"]
wt_fasta_path = "{kelsic_wt_fasta.as_posix()}"

[model]
regression_type = "randomforest"
seed = 0
""".lstrip()
    )
    return config_path


def test_load_config_resolves_paths_relative_to_toml(tmp_path):
    (tmp_path / "inputs").mkdir()
    (tmp_path / "rounds").mkdir()
    (tmp_path / "inputs" / "embeddings.csv").write_text("variant,e0\nA,1.0\n")
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nACD\n")
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        """
schema_version = 1
workflow = "experimental_round"

[run]
seed = 7
output_dir = "outputs"
overwrite = false

[experiment]
protein_name = "protein"
round_name = "round-001"
number_of_variants = 3

[inputs]
embeddings_base_path = "inputs"
embeddings_file_name = "embeddings.csv"
round_base_path = "rounds"
round_file_names = ["activity.xlsx"]
wt_fasta_path = "inputs/wt.fasta"

[model]
regression_type = "randomforest"
seed = 7
""".lstrip()
    )

    config = load_config(config_path)

    assert config.schema_version == 1
    assert config.workflow == "experimental_round"
    assert config.run.output_dir == tmp_path / "outputs"
    assert config.run.seed == 7
    assert config.experiment.embeddings_base_path == tmp_path / "inputs"
    assert config.experiment.round_file_names == ["activity.xlsx"]
    assert config.model.regression_type == "randomforest"
    assert config.model.seed == 7


def test_load_config_missing_required_field_raises_clear_error(tmp_path):
    config_path = tmp_path / "bad.toml"
    config_path.write_text('schema_version = 1\nworkflow = "experimental_round"\n')

    with pytest.raises(ValueError, match="missing required section: run"):
        load_config(config_path)


def test_load_config_wrong_type_raises_clear_error(tmp_path):
    config_path = tmp_path / "bad.toml"
    config_path.write_text(
        """
schema_version = 1
workflow = "experimental_round"

[run]
seed = "zero"
output_dir = "outputs"

[experiment]
protein_name = "protein"
round_name = "round-001"

[inputs]
embeddings_base_path = "inputs"
embeddings_file_name = "embeddings.csv"
round_base_path = "rounds"
round_file_names = "activity.xlsx"
wt_fasta_path = "inputs/wt.fasta"

[model]
regression_type = "randomforest"
""".lstrip()
    )

    with pytest.raises(ValueError, match="run.seed must be int"):
        load_config(config_path)


def test_write_config_is_deterministic(tmp_path):
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        """
schema_version = 1
workflow = "experimental_round"

[run]
output_dir = "outputs"
seed = 0

[experiment]
round_name = "round-001"
protein_name = "protein"

[inputs]
round_file_names = ["activity.xlsx"]
round_base_path = "rounds"
embeddings_file_name = "embeddings.csv"
embeddings_base_path = "inputs"
wt_fasta_path = "inputs/wt.fasta"

[model]
seed = 0
regression_type = "randomforest"
""".lstrip()
    )
    config = load_config(config_path)
    first = tmp_path / "first.toml"
    second = tmp_path / "second.toml"

    write_config(config, first)
    write_config(config, second)

    assert first.read_bytes() == second.read_bytes()


def test_config_driven_experimental_round_matches_direct_call(
    tmp_path, repo_root, kelsic_wt_fasta
):
    config_output = tmp_path / "config-output"
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, config_output
    )
    config = load_config(config_path)

    direct_output = tmp_path / "direct-output"
    direct_result = evolve_experimental(
        protein_name=config.experiment.protein_name,
        round_name=config.experiment.round_name,
        embeddings_base_path=str(config.experiment.embeddings_base_path),
        embeddings_file_name=config.experiment.embeddings_file_name,
        round_base_path=str(config.experiment.round_base_path),
        round_file_names=config.experiment.round_file_names,
        wt_fasta_path=str(config.experiment.wt_fasta_path),
        rename_WT=config.experiment.rename_WT,
        number_of_variants=config.experiment.number_of_variants,
        output_dir=str(direct_output),
        regression_type=config.model.regression_type,
        model_seed=config.model.seed,
        model_params=config.model.params,
    )
    config_result = evolve_experimental_from_config(config)

    pd.testing.assert_series_equal(
        config_result[0].reset_index(drop=True), direct_result[0].reset_index(drop=True)
    )
    for actual, expected in zip(config_result[1:], direct_result[1:], strict=True):
        pd.testing.assert_frame_equal(
            actual.reset_index(drop=True), expected.reset_index(drop=True)
        )

    direct_round = direct_output / "kelsic" / "Round1"
    config_round = config_output / "kelsic" / "Round1"
    for name in LEGACY_OUTPUTS:
        assert (config_round / name).read_bytes() == (direct_round / name).read_bytes()
    for name in DUMP_OUTPUTS:
        assert (config_round / name).exists()


def test_config_driven_run_writes_round_lock_and_manifest(
    tmp_path, repo_root, kelsic_wt_fasta
):
    config_output = tmp_path / "config-output"
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, config_output
    )
    config = load_config(config_path)

    evolve_experimental_from_config(config)

    round_dir = config_output / "kelsic" / "Round1"
    lock_path = round_dir / "round.lock.json"
    manifest_path = round_dir / "manifest.json"
    assert lock_path.exists()
    assert manifest_path.exists()

    lock = json.loads(lock_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    assert lock["config_sha256"] == config_sha256(config)
    assert manifest["lock_sha256"] == sha256_file(lock_path)
    assert manifest["previous_manifest_sha256"] is None
    assert manifest["outputs"]["model_checkpoint"]["path"].endswith("model.pkl")
    assert len(manifest["outputs"]["predictions_semantic_sha256"]) == 64
    assert "df_sorted_all.csv" in {
        Path(entry["path"]).name for entry in manifest["outputs"]["legacy_files"]
    }


def test_config_driven_run_writes_planned_round_layout(
    tmp_path, repo_root, kelsic_wt_fasta
):
    config_output = tmp_path / "config-output"
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, config_output
    )
    config = load_config(config_path)

    evolve_experimental_from_config(config)

    round_dir = config_output / "kelsic" / "Round1"
    assert (round_dir / "legacy" / "df_sorted_all.csv").read_bytes() == (
        round_dir / "df_sorted_all.csv"
    ).read_bytes()
    assert (round_dir / "legacy" / "iteration.csv").read_bytes() == (
        round_dir / "iteration.csv"
    ).read_bytes()
    assert sha256_file(round_dir / "model" / "model.pkl") == sha256_file(
        round_dir / "model.pkl"
    )
    assert sha256_file(
        round_dir / "predictions" / "predictions.parquet"
    ) == sha256_file(round_dir / "predictions.parquet")
    assert (round_dir / "predictions" / "candidates.csv").read_bytes() == (
        round_dir / "df_test.csv"
    ).read_bytes()

    params = json.loads((round_dir / "model" / "params.json").read_text())
    manifest = json.loads((round_dir / "manifest.json").read_text())
    assert params["random_state"] == 0
    assert params["n_estimators"] == 100
    assert manifest["outputs"]["planned_layout"]["model_params"]["path"].endswith(
        "model/params.json"
    )


def test_config_model_seed_and_params_control_training(
    tmp_path, repo_root, kelsic_wt_fasta
):
    config_output = tmp_path / "config-output"
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, config_output
    )
    config_path.write_text(
        config_path.read_text().replace(
            '[model]\nregression_type = "randomforest"\nseed = 0\n',
            '[model]\nregression_type = "randomforest"\nseed = 7\nn_estimators = 11\nn_jobs = 1\n',
        )
    )
    config = load_config(config_path)

    evolve_experimental_from_config(config)

    round_dir = config_output / "kelsic" / "Round1"
    model = joblib.load(round_dir / "model.pkl")
    lock = json.loads((round_dir / "round.lock.json").read_text())
    assert model.random_state == 7
    assert model.n_estimators == 11
    assert model.n_jobs == 1
    assert lock["model"]["seed"] == 7
    assert lock["model"]["params"] == {"n_estimators": 11, "n_jobs": 1}
