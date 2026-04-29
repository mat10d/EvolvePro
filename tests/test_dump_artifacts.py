from __future__ import annotations

import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from evolvepro.src.data import load_experimental_data
from evolvepro.src.evolve import evolve_experimental

ROUND_BASE = Path("colab") / "rounds_data"
DUMP_FILES = {
    "model.pkl",
    "predictions.parquet",
    "metrics.json",
    "inputs.json",
    "version.txt",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _run_round(
    tmp_path: Path, repo_root: Path, kelsic_wt_fasta: Path
) -> tuple[Path, Path]:
    embeddings_file = _make_embeddings(tmp_path, repo_root, kelsic_wt_fasta)
    output_root = tmp_path / "output"
    evolve_experimental(
        protein_name="kelsic",
        round_name="Round1",
        embeddings_base_path=str(tmp_path),
        embeddings_file_name=embeddings_file.name,
        round_base_path=str(repo_root / ROUND_BASE),
        round_file_names=["kelsic_Round1.xlsx"],
        wt_fasta_path=str(kelsic_wt_fasta),
        rename_WT=False,
        number_of_variants=12,
        output_dir=str(output_root),
    )
    return output_root / "kelsic" / "Round1", embeddings_file


def test_round_dump_files_exist(tmp_path, repo_root, kelsic_wt_fasta):
    output_dir, _ = _run_round(tmp_path, repo_root, kelsic_wt_fasta)

    assert DUMP_FILES <= {path.name for path in output_dir.iterdir()}


def test_round_dump_metrics_schema(
    tmp_path, repo_root, kelsic_wt_fasta, data_regression
):
    output_dir, _ = _run_round(tmp_path, repo_root, kelsic_wt_fasta)
    metrics = json.loads((output_dir / "metrics.json").read_text())

    data_regression.check(metrics)


def test_round_dump_predictions_parquet(tmp_path, repo_root, kelsic_wt_fasta):
    output_dir, _ = _run_round(tmp_path, repo_root, kelsic_wt_fasta)
    predictions = pd.read_parquet(output_dir / "predictions.parquet")

    assert predictions.columns.tolist() == [
        "variant",
        "y_pred",
        "y_actual",
        "y_actual_scaled",
        "y_actual_binary",
        "dist_metric",
        "std_predictions",
    ]
    assert len(predictions) == 36
    assert predictions["y_pred"].is_monotonic_decreasing


def test_round_dump_model_loadable(tmp_path, repo_root, kelsic_wt_fasta):
    output_dir, embeddings_file = _run_round(tmp_path, repo_root, kelsic_wt_fasta)
    model = joblib.load(output_dir / "model.pkl")
    predictions = pd.read_parquet(output_dir / "predictions.parquet")
    embeddings = pd.read_csv(embeddings_file, index_col=0).sort_index()

    actual = pd.Series(model.predict(embeddings), index=embeddings.index, name="actual")
    expected = predictions.set_index("variant")["y_pred"]

    pd.testing.assert_series_equal(
        actual.loc[expected.index].round(10),
        expected.round(10).rename("actual"),
        check_names=True,
    )


def test_round_dump_inputs_hashes(tmp_path, repo_root, kelsic_wt_fasta):
    output_dir, embeddings_file = _run_round(tmp_path, repo_root, kelsic_wt_fasta)
    inputs = json.loads((output_dir / "inputs.json").read_text())
    by_name = {Path(entry["path"]).name: entry for entry in inputs["files"]}

    for path in [
        embeddings_file,
        kelsic_wt_fasta,
        repo_root / ROUND_BASE / "kelsic_Round1.xlsx",
    ]:
        entry = by_name[path.name]
        assert entry["sha256"] == _sha256(path)
        assert entry["size"] == path.stat().st_size


def test_round_dump_version_txt(tmp_path, repo_root, kelsic_wt_fasta):
    output_dir, _ = _run_round(tmp_path, repo_root, kelsic_wt_fasta)
    version = (output_dir / "version.txt").read_text().splitlines()

    assert version[0].startswith("evolvepro==")
    assert version[1].startswith("git_commit=")
