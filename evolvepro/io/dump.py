from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import joblib
import pandas as pd

from evolvepro.io.hashing import sha256_file

PREDICTION_COLUMNS = [
    "variant",
    "y_pred",
    "y_actual",
    "y_actual_scaled",
    "y_actual_binary",
    "dist_metric",
    "std_predictions",
]


def _json_default(value: Any) -> Any:
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _write_text_atomic(path: Path, content: str) -> None:
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(content)
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def _write_bytes_atomic(path: Path, writer) -> None:
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        tmp_path = Path(handle.name)
    try:
        writer(tmp_path)
        tmp_path.replace(path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def dump_round(
    output_dir: Path,
    *,
    model: Any,
    predictions: pd.DataFrame,
    metrics: dict,
    inputs: list[Path],
    config: dict,
    evolvepro_version: str,
) -> None:
    """Write reproducibility artifacts for one round without mutating legacy CSVs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    prediction_frame = predictions.loc[:, PREDICTION_COLUMNS].reset_index(drop=True)
    _write_bytes_atomic(output_dir / "model.pkl", lambda path: joblib.dump(model, path))
    _write_bytes_atomic(
        output_dir / "predictions.parquet",
        lambda path: prediction_frame.to_parquet(path, index=False),
    )

    metrics_payload = {**metrics, "config": config}
    _write_text_atomic(
        output_dir / "metrics.json",
        json.dumps(metrics_payload, indent=2, sort_keys=True, default=_json_default)
        + "\n",
    )

    input_payload = {
        "files": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "size": path.stat().st_size,
            }
            for path in inputs
        ]
    }
    _write_text_atomic(
        output_dir / "inputs.json",
        json.dumps(input_payload, indent=2, sort_keys=True) + "\n",
    )

    _write_text_atomic(
        output_dir / "version.txt",
        f"evolvepro=={evolvepro_version}\ngit_commit={_git_commit()}\n",
    )
