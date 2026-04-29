from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from evolvepro.io.config import load_config
from evolvepro.io.hashing import sha256_file
from evolvepro.io.manifest import (
    build_round_lock,
    manifest_sha256,
    validate_previous_manifest,
    validate_round_lock,
    write_round_lock,
)


def _write_config_tree(tmp_path: Path, previous_manifest: Path | None = None) -> Path:
    (tmp_path / "inputs").mkdir()
    (tmp_path / "rounds").mkdir()
    (tmp_path / "outputs").mkdir()
    (tmp_path / "inputs" / "embeddings.csv").write_text("variant,e0\nA,1.0\n")
    (tmp_path / "inputs" / "wt.fasta").write_text(">WT\nACD\n")
    (tmp_path / "rounds" / "activity.csv").write_text("variant,activity\nA,1.0\n")

    previous_line = ""
    if previous_manifest is not None:
        previous_line = f'previous_manifest = "{previous_manifest.as_posix()}"\n'

    config_path = tmp_path / "evolvepro.toml"
    config_path.write_text(
        f"""
schema_version = 1
workflow = "experimental_round"

[run]
seed = 0
output_dir = "outputs"
overwrite = false
{previous_line}
[experiment]
protein_name = "protein"
round_name = "round-001"
number_of_variants = 3

[inputs]
embeddings_base_path = "inputs"
embeddings_file_name = "embeddings.csv"
round_base_path = "rounds"
round_file_names = ["activity.csv"]
wt_fasta_path = "inputs/wt.fasta"

[model]
regression_type = "randomforest"
seed = 0
""".lstrip()
    )
    return config_path


def test_round_lock_captures_input_hashes_and_model_params(tmp_path):
    config = load_config(_write_config_tree(tmp_path))

    lock = build_round_lock(config)

    by_name = {Path(entry.path).name: entry for entry in lock.inputs}
    for name in ["embeddings.csv", "activity.csv", "wt.fasta"]:
        path = tmp_path / "inputs" / name
        if name == "activity.csv":
            path = tmp_path / "rounds" / name
        assert by_name[name].sha256 == sha256_file(path)
        assert by_name[name].size == path.stat().st_size
    assert lock.model.regression_type == "randomforest"
    assert lock.model.seed == 0
    assert lock.previous_manifest_sha256 is None


def test_write_round_lock_is_deterministic(tmp_path):
    config = load_config(_write_config_tree(tmp_path))
    lock = build_round_lock(config)
    first = tmp_path / "first.lock.json"
    second = tmp_path / "second.lock.json"

    write_round_lock(first, lock)
    write_round_lock(second, lock)

    assert first.read_bytes() == second.read_bytes()
    assert json.loads(first.read_text())["schema_version"] == 1


def test_validate_round_lock_stops_on_input_hash_mismatch(tmp_path):
    config = load_config(_write_config_tree(tmp_path))
    lock = build_round_lock(config)
    (tmp_path / "inputs" / "embeddings.csv").write_text("variant,e0\nA,2.0\n")

    with pytest.raises(ValueError, match="input hash mismatch"):
        validate_round_lock(lock)


def test_validate_previous_manifest_stops_on_hash_mismatch(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text('{"schema_version":1}\n')

    with pytest.raises(ValueError, match="previous manifest hash mismatch"):
        validate_previous_manifest("0" * 64, manifest)


def test_round_lock_records_previous_manifest_hash(tmp_path):
    previous_manifest = tmp_path / "previous-manifest.json"
    previous_manifest.write_text('{"schema_version":1,"round":"round-000"}\n')
    config = load_config(_write_config_tree(tmp_path, previous_manifest))

    lock = build_round_lock(config)

    assert lock.previous_manifest_path == str(previous_manifest)
    assert lock.previous_manifest_sha256 == manifest_sha256(previous_manifest)


def test_semantic_table_hash_can_be_used_for_manifest_payload():
    from evolvepro.io.canonical import semantic_hash_table

    table = pd.DataFrame({"variant": ["B", "A"], "score": [0.2, 0.1]})

    assert len(semantic_hash_table(table, sort_by=["variant"])) == 64
