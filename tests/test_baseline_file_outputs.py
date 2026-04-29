from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from evolvepro.src.data import load_experimental_data
from evolvepro.src.evolve import evolve_experimental
from evolvepro.src.process import generate_single_aa_mutants, generate_wt

ROUND_BASE = Path("colab") / "rounds_data"
EVOLVE_EXPERIMENTAL_ARTIFACTS = [
    "iteration.csv",
    "this_round_variants.csv",
    "df_test.csv",
    "df_sorted_all.csv",
]


def _make_kelsic_synthetic_embeddings(
    tmp_path: Path, repo_root: Path, kelsic_wt_fasta: Path
) -> Path:
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
    embeddings_dir = tmp_path / "embeddings"
    embeddings_dir.mkdir()
    embeddings.to_csv(embeddings_dir / "kelsic_synthetic.csv")
    return embeddings_dir


def _run_evolve_experimental_round1(
    tmp_path: Path, repo_root: Path, kelsic_wt_fasta: Path
) -> Path:
    embeddings_dir = _make_kelsic_synthetic_embeddings(
        tmp_path, repo_root, kelsic_wt_fasta
    )
    output_root = tmp_path / "output"
    evolve_experimental(
        protein_name="kelsic",
        round_name="Round1",
        embeddings_base_path=str(embeddings_dir),
        embeddings_file_name="kelsic_synthetic.csv",
        round_base_path=str(repo_root / ROUND_BASE),
        round_file_names=["kelsic_Round1.xlsx"],
        wt_fasta_path=str(kelsic_wt_fasta),
        rename_WT=False,
        number_of_variants=12,
        output_dir=str(output_root),
    )
    return output_root / "kelsic" / "Round1"


@pytest.mark.parametrize("artifact", EVOLVE_EXPERIMENTAL_ARTIFACTS)
def test_evolve_experimental_round1_file_outputs(
    artifact, tmp_path, repo_root, kelsic_wt_fasta, file_regression
):
    output_dir = _run_evolve_experimental_round1(tmp_path, repo_root, kelsic_wt_fasta)
    artifact_path = output_dir / artifact

    assert artifact_path.exists()
    file_regression.check(
        artifact_path.read_text(),
        extension=".csv",
        basename=f"evolve_experimental_round1_{artifact_path.stem}",
    )


def test_generate_wt_kelsic_fasta_file(tmp_path, kelsic_wt_sequence, file_regression):
    output = tmp_path / "kelsic_WT.fasta"
    generate_wt(kelsic_wt_sequence, str(output))

    file_regression.check(
        output.read_text(), extension=".fasta", basename="kelsic_wt_fasta"
    )


def test_generate_single_aa_mutants_kelsic_subset_fasta_file(
    tmp_path, kelsic_wt_fasta, file_regression
):
    output = tmp_path / "kelsic_subset.fasta"
    generate_single_aa_mutants(str(kelsic_wt_fasta), str(output), positions=[1, 2, 3])

    file_regression.check(
        output.read_text(), extension=".fasta", basename="kelsic_subset_mutants_fasta"
    )
