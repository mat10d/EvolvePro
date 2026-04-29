from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

from evolvepro.src.evolve import grid_search
from tests.fixtures.synthetic import add_activity_metadata, make_embeddings


def _write_dms_inputs(base: Path) -> tuple[Path, Path]:
    dataset_name = "synthetic_dms"
    model_name = "synthetic_model"
    embeddings_dir = base / "embeddings"
    labels_dir = base / "labels"
    embeddings_dir.mkdir()
    labels_dir.mkdir()

    embeddings = make_embeddings(n_variants=10, n_dim=4, seed=17, include_wt=True)
    variants = sorted(embeddings.index.tolist())
    labels = pd.DataFrame(
        {
            "variant": variants,
            "activity": [1.0 + 0.1 * index for index in range(len(variants))],
        }
    )
    labels = add_activity_metadata(labels)

    embeddings.to_csv(embeddings_dir / f"{dataset_name}_{model_name}.csv")
    labels.to_csv(labels_dir / f"{dataset_name}_labels.csv", index=False)
    return embeddings_dir, labels_dir


def _run_evolvepro(args: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["evolvepro", *args],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def test_cli_dms_grid_search_matches_api(tmp_path):
    dataset_name = "synthetic_dms"
    model_name = "synthetic_model"
    experiment_name = "smoke"
    embeddings_dir, labels_dir = _write_dms_inputs(tmp_path)
    api_output = tmp_path / "api-output"
    cli_output = tmp_path / "cli-output"

    grid_search(
        dataset_name=dataset_name,
        experiment_name=experiment_name,
        model_name=model_name,
        embeddings_path=str(embeddings_dir),
        labels_path=str(labels_dir),
        num_simulations=1,
        num_iterations=[1],
        measured_var=["activity"],
        learning_strategies=["topn"],
        num_mutants_per_round=[2],
        num_final_round_mutants=2,
        first_round_strategies=["random"],
        embedding_types=["embeddings"],
        pca_components=None,
        regression_types=["ridge"],
        embeddings_file_type="csv",
        output_dir=str(api_output),
        embeddings_type_pt=None,
    )

    run = _run_evolvepro(
        [
            "dms",
            "grid-search",
            "--dataset-name",
            dataset_name,
            "--experiment-name",
            experiment_name,
            "--model-name",
            model_name,
            "--embeddings-path",
            str(embeddings_dir),
            "--labels-path",
            str(labels_dir),
            "--num-simulations",
            "1",
            "--num-iterations",
            "1",
            "--measured-var",
            "activity",
            "--learning-strategies",
            "topn",
            "--num-mutants-per-round",
            "2",
            "--num-final-round-mutants",
            "2",
            "--first-round-strategies",
            "random",
            "--embedding-types",
            "embeddings",
            "--regression-types",
            "ridge",
            "--embeddings-file-type",
            "csv",
            "--output-dir",
            str(cli_output),
        ],
        cwd=tmp_path,
    )

    output_name = f"{dataset_name}_{model_name}_{experiment_name}.csv"
    assert "completed DMS grid search" in run.stdout
    assert (cli_output / output_name).read_bytes() == (
        api_output / output_name
    ).read_bytes()
