from __future__ import annotations

import pandas as pd

from evolvepro.io.canonical import canonical_json_bytes, semantic_hash_table


def test_canonical_json_bytes_are_key_order_independent():
    left = {"b": [2, 1], "a": {"d": False, "c": "value"}}
    right = {"a": {"c": "value", "d": False}, "b": [2, 1]}

    assert canonical_json_bytes(left) == canonical_json_bytes(right)
    assert canonical_json_bytes(left) == b'{"a":{"c":"value","d":false},"b":[2,1]}\n'


def test_semantic_hash_table_is_row_order_independent_when_sorted():
    left = pd.DataFrame(
        {
            "variant": ["B", "A"],
            "score": [1.25, 0.5],
            "called": [True, False],
        }
    )
    right = left.iloc[[1, 0]].reset_index(drop=True)

    assert semantic_hash_table(left, sort_by=["variant"]) == semantic_hash_table(
        right, sort_by=["variant"]
    )


def test_semantic_hash_table_changes_when_values_change():
    original = pd.DataFrame({"variant": ["A"], "score": [0.5]})
    changed = pd.DataFrame({"variant": ["A"], "score": [0.5001]})

    assert semantic_hash_table(original, sort_by=["variant"]) != semantic_hash_table(
        changed, sort_by=["variant"]
    )
