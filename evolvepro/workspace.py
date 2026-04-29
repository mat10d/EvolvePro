from __future__ import annotations

from pathlib import Path

from evolvepro.io.canonical import write_json_atomic

WORKSPACE_DIR = ".evolvepro"
CONFIG_FILE = "evolvepro.toml"
LOCK_FILE = "evolvepro.lock"
WORKSPACE_FILE = "workspace.toml"


def _default_config() -> str:
    return (
        "schema_version = 1\n"
        'workflow = "experimental_round"\n'
        "\n"
        "[run]\n"
        "seed = 0\n"
        'output_dir = "outputs"\n'
        "overwrite = false\n"
        "\n"
        "[experiment]\n"
        'protein_name = "example"\n'
        'round_name = "round-001"\n'
        "rename_WT = false\n"
        "number_of_variants = 12\n"
        "\n"
        "[inputs]\n"
        'embeddings_base_path = "inputs"\n'
        'embeddings_file_name = "embeddings.csv"\n'
        'round_base_path = "rounds/round-001/inputs"\n'
        'round_file_names = ["activity.xlsx"]\n'
        'wt_fasta_path = "inputs/wt.fasta"\n'
        "\n"
        "[model]\n"
        'regression_type = "randomforest"\n'
        "seed = 0\n"
        "n_jobs = 1\n"
    )


def _workspace_metadata() -> str:
    return 'schema_version = 1\nworkspace_version = 1\nconfig = "evolvepro.toml"\n'


def _workspace_lock() -> dict[str, object]:
    return {
        "schema_version": 1,
        "lock_version": 1,
        "workspace": {"config": CONFIG_FILE},
        "current_round": None,
        "rounds": [],
    }


def _write_if_needed(path: Path, content: str, *, force: bool) -> None:
    if path.exists() and not force:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _write_json_if_needed(
    path: Path, content: dict[str, object], *, force: bool
) -> None:
    if path.exists() and not force:
        return
    write_json_atomic(path, content)


def init_workspace(path: Path, *, force: bool = False) -> Path:
    workspace = path.expanduser()
    if workspace.exists() and not workspace.is_dir():
        raise ValueError(f"workspace path is not a directory: {workspace}")
    workspace.mkdir(parents=True, exist_ok=True)

    for directory in [
        workspace / "inputs",
        workspace / "outputs",
        workspace / "rounds",
        workspace / WORKSPACE_DIR / "cache",
        workspace / WORKSPACE_DIR / "locks",
        workspace / WORKSPACE_DIR / "cas" / "sha256",
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    _write_if_needed(workspace / CONFIG_FILE, _default_config(), force=force)
    _write_if_needed(
        workspace / WORKSPACE_DIR / WORKSPACE_FILE,
        _workspace_metadata(),
        force=force,
    )
    _write_json_if_needed(workspace / LOCK_FILE, _workspace_lock(), force=force)
    return workspace


def find_workspace(start: Path) -> Path | None:
    current = start.expanduser().resolve()
    if current.is_file():
        current = current.parent

    for candidate in [current, *current.parents]:
        if (candidate / WORKSPACE_DIR / WORKSPACE_FILE).is_file() and (
            candidate / CONFIG_FILE
        ).is_file():
            return candidate
    return None
