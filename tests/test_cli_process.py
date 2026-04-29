from __future__ import annotations

import subprocess
from pathlib import Path

from evolvepro.src.process import generate_single_aa_mutants, generate_wt


def _run_evolvepro(args: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["evolvepro", *args],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def test_cli_process_generate_wt_matches_api(tmp_path, kelsic_wt_sequence):
    cli_output = tmp_path / "cli_WT.fasta"
    api_output = tmp_path / "api_WT.fasta"

    run = _run_evolvepro(
        [
            "process",
            "generate-wt",
            "--sequence",
            kelsic_wt_sequence,
            "--output",
            str(cli_output),
        ],
        cwd=tmp_path,
    )
    generate_wt(kelsic_wt_sequence, str(api_output))

    assert "wrote WT FASTA" in run.stdout
    assert cli_output.read_bytes() == api_output.read_bytes()


def test_cli_process_generate_single_aa_mutants_matches_api(tmp_path, kelsic_wt_fasta):
    cli_output = tmp_path / "cli_mutants.fasta"
    api_output = tmp_path / "api_mutants.fasta"

    run = _run_evolvepro(
        [
            "process",
            "generate-single-aa-mutants",
            "--wt-fasta",
            str(kelsic_wt_fasta),
            "--output",
            str(cli_output),
            "--positions",
            "1,2,3",
        ],
        cwd=tmp_path,
    )
    generate_single_aa_mutants(
        str(kelsic_wt_fasta), str(api_output), positions=[1, 2, 3]
    )

    assert "wrote single amino acid mutant FASTA" in run.stdout
    assert cli_output.read_bytes() == api_output.read_bytes()
