from __future__ import annotations

import json
import stat

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
