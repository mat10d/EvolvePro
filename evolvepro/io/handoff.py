from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class CandidateHandoff:
    table: pd.DataFrame
    measured: pd.DataFrame
    unmeasured: pd.DataFrame


def write_candidate_handoff(path: Path, variants: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant", "activity"])
        for variant in variants:
            writer.writerow([variant, ""])


def _parse_activity(value: str, *, variant: str) -> float | None:
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(
            f"activity must be numeric for variant {variant}: {value}"
        ) from error


def read_candidate_handoff(path: Path) -> CandidateHandoff:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["variant", "activity"]:
            raise ValueError(
                "candidate handoff CSV must have exactly two columns: variant,activity"
            )

        variants: list[str] = []
        activities: list[float | None] = []
        for row in reader:
            variant = row["variant"].strip()
            if variant == "":
                raise ValueError("variant must not be blank")
            variants.append(variant)
            activities.append(_parse_activity(row["activity"], variant=variant))

    duplicated = pd.Series(variants).duplicated(keep=False)
    if duplicated.any():
        duplicates = sorted(set(pd.Series(variants)[duplicated].tolist()))
        raise ValueError(f"candidate handoff contains duplicate variants: {duplicates}")

    table = pd.DataFrame(
        {
            "variant": pd.Series(variants, dtype=object),
            "activity": pd.Series(activities, dtype=object),
        }
    )
    measured_mask = table["activity"].notna()
    return CandidateHandoff(
        table=table,
        measured=table[measured_mask].reset_index(drop=True),
        unmeasured=table[~measured_mask].reset_index(drop=True),
    )
