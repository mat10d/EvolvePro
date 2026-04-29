from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

from evolvepro.io.config import evolve_experimental_from_config, load_config
from tests.test_config import _write_round_config

LEGACY_OUTPUTS = [
    "iteration.csv",
    "this_round_variants.csv",
    "df_test.csv",
    "df_sorted_all.csv",
]


def _run_evolvepro(args: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["evolvepro", *args],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def test_experimental_round_cli_workflow_matches_config_api(
    tmp_path, repo_root, kelsic_wt_fasta
):
    cli_output = tmp_path / "cli-output"
    cli_config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, cli_output
    )

    direct_output = tmp_path / "direct-output"
    direct_config_path = tmp_path / "direct.toml"
    direct_config_path.write_text(
        cli_config_path.read_text().replace(
            cli_output.as_posix(), direct_output.as_posix()
        )
    )
    evolve_experimental_from_config(load_config(direct_config_path))

    run = _run_evolvepro(["run", "--config", cli_config_path.name], cwd=tmp_path)
    cli_round = cli_output / "kelsic" / "Round1"
    direct_round = direct_output / "kelsic" / "Round1"
    validate = _run_evolvepro(["validate", str(cli_round)], cwd=tmp_path)

    assert "completed run from config" in run.stdout
    assert "validated round" in validate.stdout
    for name in LEGACY_OUTPUTS:
        assert (cli_round / name).read_bytes() == (direct_round / name).read_bytes()
        assert (cli_round / "legacy" / name).read_bytes() == (
            cli_round / name
        ).read_bytes()

    pd.testing.assert_frame_equal(
        pd.read_parquet(cli_round / "predictions" / "predictions.parquet"),
        pd.read_parquet(direct_round / "predictions" / "predictions.parquet"),
    )
    assert (cli_round / "round.lock.json").exists()
    assert (cli_round / "manifest.json").exists()
    assert (cli_round / "model" / "params.json").exists()
