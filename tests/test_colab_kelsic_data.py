from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from evolvepro.src.data import create_iteration_dataframes, load_experimental_data
from evolvepro.src.evolve import evolve_experimental
from evolvepro.src.plot import read_exp_data

ROUND_FILES = [f"kelsic_Round{i}.xlsx" for i in range(1, 6)]


def _round_base(repo_root: Path) -> Path:
    return repo_root / "colab" / "rounds_data"


def _round_summary(df: pd.DataFrame) -> dict:
    return {
        "shape": list(df.shape),
        "columns": df.columns.tolist(),
        "head": df.head(5).to_dict(orient="records"),
        "tail": df.tail(5).to_dict(orient="records"),
    }


def _load_rounds(repo_root: Path, kelsic_wt_fasta: Path) -> list[pd.DataFrame]:
    base = _round_base(repo_root)
    return [
        load_experimental_data(str(base), round_file, str(kelsic_wt_fasta))
        for round_file in ROUND_FILES
    ]


def test_load_kelsic_round1(repo_root, kelsic_wt_fasta, data_regression):
    df = load_experimental_data(
        str(_round_base(repo_root)), "kelsic_Round1.xlsx", str(kelsic_wt_fasta)
    )

    data_regression.check(_round_summary(df))


def test_create_iteration_dataframes_kelsic_rounds(
    repo_root, kelsic_wt_fasta, data_regression, dataframe_regression
):
    rounds = _load_rounds(repo_root, kelsic_wt_fasta)
    measured_variants = sorted(
        {variant for df in rounds for variant in df["updated_variant"]}
    )
    expected_variants = measured_variants + ["UNMEASURED_A1C", "UNMEASURED_A1D"]

    iteration, labels = create_iteration_dataframes(rounds, expected_variants)
    assert iteration is not None
    assert labels is not None

    data_regression.check(
        {
            "iteration_shape": list(iteration.shape),
            "labels_shape": list(labels.shape),
            "iteration_counts": {
                str(key): int(value)
                for key, value in iteration["iteration"]
                .value_counts()
                .sort_index()
                .items()
            },
            "labels_iteration_counts": {
                str(key): int(value)
                for key, value in labels["iteration"]
                .value_counts(dropna=False)
                .sort_index()
                .items()
            },
            "activity_min": round(float(labels["activity"].min()), 6),
            "activity_max": round(float(labels["activity"].max()), 6),
        }
    )
    dataframe_regression.check(
        labels.head(12).round(6).reset_index(drop=True), basename="kelsic_labels_head"
    )
    dataframe_regression.check(
        labels.tail(8).round(6).reset_index(drop=True), basename="kelsic_labels_tail"
    )


def test_read_exp_data_kelsic_rounds(repo_root, kelsic_wt_fasta, data_regression):
    df = read_exp_data(str(_round_base(repo_root)), ROUND_FILES, str(kelsic_wt_fasta))

    data_regression.check(
        {
            "shape": list(df.shape),
            "columns": df.columns.tolist(),
            "iteration_counts": {
                str(key): int(value)
                for key, value in df["iteration"].value_counts().sort_index().items()
            },
            "activity_by_iteration": {
                str(key): round(float(value), 6)
                for key, value in df.groupby("iteration")["activity"].mean().items()
            },
            "head": df.head(5).round(6).to_dict(orient="records"),
            "tail": df.tail(5).round(6).to_dict(orient="records"),
        }
    )


def test_evolve_experimental_kelsic_round1_synthetic_embeddings(
    tmp_path, repo_root, kelsic_wt_fasta, data_regression
):
    round1 = load_experimental_data(
        str(_round_base(repo_root)), "kelsic_Round1.xlsx", str(kelsic_wt_fasta)
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
    embeddings_dir = tmp_path / "embeddings"
    embeddings_dir.mkdir()
    embeddings.to_csv(embeddings_dir / "kelsic_synthetic.csv")

    this_round_variants, df_test, df_sorted_all = evolve_experimental(
        protein_name="kelsic",
        round_name="Round1",
        embeddings_base_path=str(embeddings_dir),
        embeddings_file_name="kelsic_synthetic.csv",
        round_base_path=str(_round_base(repo_root)),
        round_file_names=["kelsic_Round1.xlsx"],
        wt_fasta_path=str(kelsic_wt_fasta),
        rename_WT=False,
        number_of_variants=12,
        output_dir=str(tmp_path / "output"),
    )

    data_regression.check(
        {
            "this_round_variants": this_round_variants.tolist(),
            "df_test_top12": df_test.sort_values("y_pred", ascending=False)
            .head(12)
            .round(6)
            .to_dict(orient="records"),
            "df_sorted_all_top12": df_sorted_all.head(12)
            .round(6)
            .to_dict(orient="records"),
        }
    )
