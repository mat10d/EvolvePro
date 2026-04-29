"""Deterministic fixture builders for Stage-0 regression tests.

RNG contract:
- embeddings use ``numpy.random.default_rng(seed).normal`` and are float64.
- labels use the same RNG API for small additive noise; the activity signal is a
  deterministic function of the variant index.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def make_embeddings(
    n_variants: int, n_dim: int, seed: int, include_wt: bool = True
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    values = rng.normal(size=(n_variants, n_dim)).astype("float64")
    variants = [f"variant_{i}" for i in range(n_variants)]
    columns = [f"embedding_{i}" for i in range(n_dim)]
    df = pd.DataFrame(values, index=variants, columns=columns)
    if include_wt:
        df.loc["WT"] = np.zeros(n_dim, dtype="float64")
    return df


def make_labels(
    variants: list[str], seed: int, include_metadata: bool = True
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x = np.arange(len(variants), dtype="float64")
    noise = rng.normal(scale=0.03, size=len(variants)).astype("float64")
    activity = (1.0 + 0.015 * x + 0.2 * np.sin(x / 4.0) + noise).astype("float64")
    labels = pd.DataFrame({"variant": variants, "activity": activity})
    if include_metadata:
        labels = add_activity_metadata(labels)
    return labels


def add_activity_metadata(labels: pd.DataFrame) -> pd.DataFrame:
    labels = labels.copy()
    activity = labels["activity"].astype("float64")
    labels["activity_binary"] = (activity >= 1.0).astype("int64")
    labels["activity_scaled"] = (
        (activity - float(activity.min()))
        / (float(activity.max()) - float(activity.min()))
    ).astype("float64")
    return labels
