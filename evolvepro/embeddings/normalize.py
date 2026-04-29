from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from evolvepro.io.cas import snapshot_file
from evolvepro.io.canonical import semantic_hash_table
from evolvepro.io.hashing import sha256_file
from evolvepro.workspace import WORKSPACE_DIR


def _relative_path(path: Path, workspace: Path) -> str:
    try:
        return path.resolve().relative_to(workspace.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _variant_strings(variants: pd.Series) -> pd.Series:
    if variants.isna().any():
        raise ValueError("embedding table contains missing variants")

    normalized = variants.astype(str)
    if normalized.str.strip().eq("").any():
        raise ValueError("embedding table contains blank variants")
    return normalized


def _duplicate_variants(variants: pd.Series) -> list[str]:
    duplicated = variants.duplicated(keep=False)
    if not duplicated.any():
        return []
    return sorted(set(variants.loc[duplicated].tolist()))


def _numeric_feature_column(column: pd.Series) -> pd.Series:
    try:
        numeric = pd.to_numeric(column, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"feature column '{column.name}' must be numeric") from exc

    if pd.api.types.is_complex_dtype(numeric.dtype):
        raise ValueError(f"feature column '{column.name}' must be real-valued numeric")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(
            f"feature column '{column.name}' must contain finite numeric values"
        )
    return numeric


def normalize_embedding_table(table: pd.DataFrame) -> pd.DataFrame:
    variant_column_count = list(table.columns).count("variant")
    if variant_column_count > 1:
        raise ValueError("embedding table must contain exactly one variant column")

    if variant_column_count == 1:
        variants = table["variant"]
        feature_columns = [column for column in table.columns if column != "variant"]
    else:
        variants = pd.Series(table.index, index=table.index, name="variant")
        feature_columns = list(table.columns)

    if not feature_columns:
        raise ValueError("embedding table must contain at least one feature column")

    normalized_variants = _variant_strings(variants)
    duplicates = _duplicate_variants(normalized_variants)
    if duplicates:
        raise ValueError(f"duplicate variants in embedding table: {duplicates}")

    normalized = pd.DataFrame({"variant": normalized_variants})
    for column in feature_columns:
        normalized[column] = _numeric_feature_column(table[column])

    return normalized.sort_values("variant", kind="mergesort").reset_index(drop=True)


def write_embedding_csv(path: Path, table: pd.DataFrame) -> Path:
    normalized = normalize_embedding_table(table)
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
    return path


def embedding_table_semantic_sha256(table: pd.DataFrame) -> str:
    normalized = normalize_embedding_table(table)
    return semantic_hash_table(normalized, sort_by=["variant"])


def embedding_file_artifact(
    path: Path, *, workspace: Path | None = None
) -> dict[str, object]:
    table = pd.read_csv(path)
    artifact: dict[str, object] = {
        "path": _relative_path(path, workspace)
        if workspace is not None
        else path.as_posix(),
        "sha256": sha256_file(path),
        "size": path.stat().st_size,
        "semantic_sha256": embedding_table_semantic_sha256(table),
    }
    return artifact


def snapshot_embedding_csv(path: Path, workspace: Path) -> dict[str, object]:
    snapshot = snapshot_file(path, workspace / WORKSPACE_DIR / "cas" / "sha256")
    artifact = embedding_file_artifact(path, workspace=workspace)
    artifact["cas_path"] = _relative_path(snapshot.path, workspace)
    return artifact
