from __future__ import annotations

import hashlib
import json
import tempfile
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


def _canonical_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _canonical_value(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _canonical_value(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_canonical_value(item) for item in value]
    if isinstance(value, Path):
        return value.as_posix()
    if hasattr(value, "item"):
        return _canonical_value(value.item())
    return value


def canonical_json_bytes(data: Mapping[str, Any]) -> bytes:
    payload = _canonical_value(data)
    return (
        json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def canonical_json_sha256(data: Mapping[str, Any]) -> str:
    return hashlib.sha256(canonical_json_bytes(data)).hexdigest()


def write_json_atomic(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=path.parent, delete=False) as handle:
        handle.write(canonical_json_bytes(data))
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def _table_value(value: Any) -> Any:
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float):
        return format(value, ".17g")
    return value


def semantic_hash_table(df: pd.DataFrame, *, sort_by: list[str]) -> str:
    missing = [column for column in sort_by if column not in df.columns]
    if missing:
        raise ValueError(f"sort columns missing from table: {missing}")

    canonical = df.sort_values(sort_by, kind="mergesort").reset_index(drop=True)
    payload = {
        "columns": [str(column) for column in canonical.columns],
        "rows": [
            [_table_value(value) for value in row]
            for row in canonical.itertuples(index=False, name=None)
        ],
    }
    return canonical_json_sha256(payload)
