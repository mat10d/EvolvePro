from __future__ import annotations

import subprocess
from pathlib import Path

from evolvepro.src.process import process_dataset


def _run_evolvepro(args: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["evolvepro", *args],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def test_cli_process_dataset_matches_api(tmp_path, repo_root):
    source = repo_root / "data" / "dms" / "activity" / "sarscov2_preprocessed.csv"
    wt_fasta = repo_root / "data" / "dms" / "wt_fasta" / "cov2_S_WT.fasta"
    api_output = tmp_path / "api-output"
    cli_output = tmp_path / "cli-output"

    process_dataset(
        file_path=str(source),
        wt_fasta_path=str(wt_fasta),
        dataset_name="sarscov2",
        activity_column="mut_escape",
        cutoff_value=0.1,
        output_dir=str(api_output),
        cutoff_rule="greater_than",
        cutoff_percentiles=[90, 95],
        drop_columns=True,
    )

    run = _run_evolvepro(
        [
            "process",
            "dataset",
            "--file-path",
            str(source),
            "--wt-fasta",
            str(wt_fasta),
            "--dataset-name",
            "sarscov2",
            "--activity-column",
            "mut_escape",
            "--cutoff-value",
            "0.1",
            "--output-dir",
            str(cli_output),
            "--cutoff-rule",
            "greater_than",
            "--cutoff-percentiles",
            "90",
            "95",
            "--drop-columns",
        ],
        cwd=tmp_path,
    )

    assert "processed dataset" in run.stdout
    for name in ["sarscov2.fasta", "sarscov2_labels.csv"]:
        assert (cli_output / name).read_bytes() == (api_output / name).read_bytes()
