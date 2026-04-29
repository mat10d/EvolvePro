from __future__ import annotations

import builtins
from importlib.metadata import entry_points

from evolvepro.cli import main
from tests.test_config import _write_round_config


def test_console_script_entrypoint_registered():
    scripts = entry_points(group="console_scripts")
    assert any(script.name == "evolvepro" for script in scripts)


def test_cli_init_creates_workspace(tmp_path, capsys):
    assert main(["init", str(tmp_path)]) == 0

    assert (tmp_path / "evolvepro.toml").exists()
    assert (tmp_path / "evolvepro.lock").exists()
    assert (tmp_path / ".evolvepro" / "workspace.toml").exists()
    assert (tmp_path / ".evolvepro" / "cas" / "sha256").is_dir()
    assert "initialized workspace" in capsys.readouterr().out


def test_cli_init_reconfigure_requires_confirmation(tmp_path, monkeypatch):
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text("user config\n")
    monkeypatch.setattr(builtins, "input", lambda _prompt: "n")

    assert main(["init", str(tmp_path), "--reconfigure"]) == 1

    assert config_path.read_text() == "user config\n"


def test_cli_init_reconfigure_yes_overwrites_config(tmp_path):
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text("user config\n")

    assert main(["init", str(tmp_path), "--reconfigure", "--yes"]) == 0

    assert config_path.read_text().startswith("schema_version = 1\n")


def test_cli_init_set_default_is_non_interactive(tmp_path):
    assert main(["init", str(tmp_path), "--set-default", "--yes"]) == 0

    assert (tmp_path / "evolvepro.toml").exists()


def test_cli_run_config_writes_round_manifest(tmp_path, repo_root, kelsic_wt_fasta):
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, tmp_path / "output"
    )

    assert main(["run", "--config", str(config_path)]) == 0

    assert (tmp_path / "output" / "kelsic" / "Round1" / "manifest.json").exists()
    assert (tmp_path / "output" / "kelsic" / "Round1" / "round.lock.json").exists()


def test_cli_validate_completed_round(tmp_path, repo_root, kelsic_wt_fasta, capsys):
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, tmp_path / "output"
    )
    assert main(["run", "--config", str(config_path)]) == 0
    round_dir = tmp_path / "output" / "kelsic" / "Round1"

    assert main(["validate", str(round_dir)]) == 0

    assert "validated round" in capsys.readouterr().out


def test_cli_no_args_prints_help(capsys):
    assert main([]) == 0

    output = capsys.readouterr().out
    assert "usage:" in output
    assert "init" in output
    assert "run" in output
    assert "validate" in output


def test_cli_run_uses_workspace_config_by_default(tmp_path, repo_root, kelsic_wt_fasta):
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, tmp_path / "output"
    )
    workspace_config = tmp_path / "evolvepro.toml"
    workspace_config.write_text(config_path.read_text())

    assert main(["run"], cwd=tmp_path) == 0

    assert (tmp_path / "output" / "kelsic" / "Round1" / "manifest.json").exists()


def test_cli_validate_uses_cwd_when_round_dir_omitted(
    tmp_path, repo_root, kelsic_wt_fasta
):
    config_path = _write_round_config(
        tmp_path, repo_root, kelsic_wt_fasta, tmp_path / "output"
    )
    assert main(["run", "--config", str(config_path)]) == 0
    round_dir = tmp_path / "output" / "kelsic" / "Round1"

    assert main(["validate"], cwd=round_dir) == 0
