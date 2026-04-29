from __future__ import annotations

import stat

import pytest

from evolvepro.io.cas import snapshot_file
from evolvepro.io.handoff import read_candidate_handoff, write_candidate_handoff
from evolvepro.io.hashing import sha256_file


def test_candidate_handoff_round_trips_partial_measurements(tmp_path):
    path = tmp_path / "candidates.csv"
    write_candidate_handoff(path, ["A12V", "G18D", "L22F"])
    path.write_text("variant,activity\nA12V,0.5\nG18D,\nL22F,1.25\n")

    handoff = read_candidate_handoff(path)

    assert handoff.table.to_dict(orient="records") == [
        {"variant": "A12V", "activity": 0.5},
        {"variant": "G18D", "activity": None},
        {"variant": "L22F", "activity": 1.25},
    ]
    assert handoff.measured["variant"].tolist() == ["A12V", "L22F"]
    assert handoff.unmeasured["variant"].tolist() == ["G18D"]


def test_candidate_handoff_requires_exact_columns(tmp_path):
    path = tmp_path / "candidates.csv"
    path.write_text("variant,activity,notes\nA12V,0.5,ok\n")

    with pytest.raises(ValueError, match="exactly two columns"):
        read_candidate_handoff(path)


def test_candidate_handoff_rejects_duplicate_variants(tmp_path):
    path = tmp_path / "candidates.csv"
    path.write_text("variant,activity\nA12V,0.5\nA12V,0.6\n")

    with pytest.raises(ValueError, match="duplicate variants"):
        read_candidate_handoff(path)


def test_candidate_handoff_rejects_non_numeric_activity(tmp_path):
    path = tmp_path / "candidates.csv"
    path.write_text("variant,activity\nA12V,strong\n")

    with pytest.raises(ValueError, match="activity must be numeric"):
        read_candidate_handoff(path)


def test_snapshot_file_writes_content_addressed_read_only_copy(tmp_path):
    source = tmp_path / "candidates.csv"
    source.write_text("variant,activity\nA12V,0.5\n")
    cas_root = tmp_path / ".evolvepro" / "cas" / "sha256"

    snapshot = snapshot_file(source, cas_root)

    assert snapshot.file_hash.sha256 == sha256_file(source)
    assert snapshot.path.read_bytes() == source.read_bytes()
    assert (
        snapshot.path
        == cas_root
        / snapshot.file_hash.sha256[:2]
        / snapshot.file_hash.sha256[2:4]
        / snapshot.file_hash.sha256
    )
    assert snapshot.path.stat().st_mode & stat.S_IWUSR == 0

    source.write_text("variant,activity\nA12V,0.7\n")
    assert snapshot.path.read_text() == "variant,activity\nA12V,0.5\n"
