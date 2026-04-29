from __future__ import annotations

import json

import pytest

from evolvepro.workspace import find_workspace, init_workspace


def test_init_workspace_creates_file_scaffold(tmp_path):
    init_workspace(tmp_path)

    assert (tmp_path / "evolvepro.toml").exists()
    assert (tmp_path / ".evolvepro" / "workspace.toml").exists()
    assert (tmp_path / "inputs").is_dir()
    assert (tmp_path / "outputs").is_dir()
    assert (tmp_path / "rounds").is_dir()
    assert (tmp_path / ".evolvepro" / "cache").is_dir()
    assert (tmp_path / ".evolvepro" / "locks").is_dir()
    assert (tmp_path / ".evolvepro" / "cas" / "sha256").is_dir()
    assert (tmp_path / "evolvepro.lock").exists()

    workspace_lock = json.loads((tmp_path / "evolvepro.lock").read_text())
    assert workspace_lock == {
        "current_round": None,
        "lock_version": 1,
        "rounds": [],
        "schema_version": 1,
        "workspace": {"config": "evolvepro.toml"},
    }

    workspace = (tmp_path / ".evolvepro" / "workspace.toml").read_text()
    assert workspace == (
        'schema_version = 1\nworkspace_version = 1\nconfig = "evolvepro.toml"\n'
    )


def test_init_workspace_is_idempotent_without_force(tmp_path):
    init_workspace(tmp_path)
    config_path = tmp_path / "evolvepro.toml"
    original = config_path.read_bytes()

    init_workspace(tmp_path)

    assert config_path.read_bytes() == original


def test_init_workspace_does_not_overwrite_workspace_lock_without_force(tmp_path):
    init_workspace(tmp_path)
    lock_path = tmp_path / "evolvepro.lock"
    lock_path.write_text('{"user":"state"}\n')

    init_workspace(tmp_path)

    assert lock_path.read_text() == '{"user":"state"}\n'


def test_init_workspace_refuses_to_overwrite_user_config_without_force(tmp_path):
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text("user config\n")

    init_workspace(tmp_path)

    assert config_path.read_text() == "user config\n"


def test_init_workspace_force_overwrites_generated_files(tmp_path):
    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text("user config\n")

    init_workspace(tmp_path, force=True)

    assert config_path.read_text().startswith("schema_version = 1\n")


def test_find_workspace_walks_up_to_initialized_directory(tmp_path):
    init_workspace(tmp_path)
    nested = tmp_path / "rounds" / "round-001" / "inputs"
    nested.mkdir(parents=True)

    assert find_workspace(nested) == tmp_path


def test_find_workspace_returns_none_when_not_initialized(tmp_path):
    assert find_workspace(tmp_path) is None


def test_init_workspace_rejects_file_path(tmp_path):
    path = tmp_path / "not-a-directory"
    path.write_text("content\n")

    with pytest.raises(ValueError, match="workspace path is not a directory"):
        init_workspace(path)
