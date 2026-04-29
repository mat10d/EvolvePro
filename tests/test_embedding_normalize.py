from __future__ import annotations

import stat

import pandas as pd
import pytest

from evolvepro.embeddings.normalize import (
    embedding_table_semantic_sha256,
    normalize_embedding_table,
    snapshot_embedding_csv,
    write_embedding_csv,
)
from evolvepro.io.hashing import sha256_file


def test_index_based_dataframe_normalizes_to_variant_column():
    table = pd.DataFrame({"feat_b": [0.2, 0.1], "feat_a": [2, 1]}, index=["B2A", "A1C"])

    normalized = normalize_embedding_table(table)

    assert normalized.columns.tolist() == ["variant", "feat_b", "feat_a"]
    assert normalized.to_dict(orient="records") == [
        {"variant": "A1C", "feat_b": 0.1, "feat_a": 1},
        {"variant": "B2A", "feat_b": 0.2, "feat_a": 2},
    ]
    assert normalized.index.tolist() == [0, 1]


def test_variant_column_dataframe_sorts_deterministically_preserving_feature_order():
    table = pd.DataFrame(
        {
            "feat_z": [2, 1],
            "variant": ["B2A", "A1C"],
            "feat_a": [0.2, 0.1],
        }
    )

    normalized = normalize_embedding_table(table)

    assert normalized.columns.tolist() == ["variant", "feat_z", "feat_a"]
    assert normalized["variant"].tolist() == ["A1C", "B2A"]
    assert normalized["feat_z"].tolist() == [1, 2]


def test_semantic_hash_is_row_order_independent():
    first = pd.DataFrame(
        {"variant": ["B2A", "A1C"], "feat_0": [0.2, 0.1], "feat_1": [2, 1]}
    )
    second = pd.DataFrame(
        {"feat_0": [0.1, 0.2], "feat_1": [1, 2]}, index=["A1C", "B2A"]
    )

    assert embedding_table_semantic_sha256(first) == embedding_table_semantic_sha256(
        second
    )


def test_duplicate_variants_raise():
    table = pd.DataFrame({"feat_0": [0.1, 0.2]}, index=["A1C", "A1C"])

    with pytest.raises(ValueError, match="duplicate variants"):
        normalize_embedding_table(table)


def test_blank_variants_raise():
    table = pd.DataFrame({"variant": ["A1C", " "], "feat_0": [0.1, 0.2]})

    with pytest.raises(ValueError, match="blank variants"):
        normalize_embedding_table(table)


def test_non_numeric_feature_values_raise():
    table = pd.DataFrame({"variant": ["A1C"], "feat_0": ["not-a-number"]})

    with pytest.raises(ValueError, match="feature column 'feat_0' must be numeric"):
        normalize_embedding_table(table)


def test_no_feature_columns_raise():
    table = pd.DataFrame({"variant": ["A1C"]})

    with pytest.raises(ValueError, match="at least one feature column"):
        normalize_embedding_table(table)


def test_write_embedding_csv_produces_byte_stable_output(tmp_path):
    path = tmp_path / "nested" / "embeddings.csv"
    table = pd.DataFrame(
        {"variant": ["B2A", "A1C"], "feat_0": [2.0, 1.25], "feat_1": [4, 3]}
    )

    returned = write_embedding_csv(path, table)

    assert returned == path
    assert path.read_bytes() == (b"variant,feat_0,feat_1\nA1C,1.25,3\nB2A,2.0,4\n")


def test_snapshot_embedding_csv_writes_cas_entry_and_manifest_payload(tmp_path):
    source = tmp_path / "embeddings.csv"
    write_embedding_csv(
        source,
        pd.DataFrame({"variant": ["B2A", "A1C"], "feat_0": [0.2, 0.1]}),
    )
    expected_sha256 = sha256_file(source)
    expected_semantic_sha256 = embedding_table_semantic_sha256(pd.read_csv(source))

    artifact = snapshot_embedding_csv(source, tmp_path)

    assert artifact["path"] == "embeddings.csv"
    assert artifact["sha256"] == expected_sha256
    assert artifact["size"] == source.stat().st_size
    assert artifact["semantic_sha256"] == expected_semantic_sha256
    assert artifact["cas_path"] == (
        f".evolvepro/cas/sha256/{expected_sha256[:2]}/"
        f"{expected_sha256[2:4]}/{expected_sha256}"
    )
    cas_path = tmp_path / str(artifact["cas_path"])
    assert cas_path.read_bytes() == source.read_bytes()
    assert cas_path.stat().st_mode & stat.S_IWUSR == 0
