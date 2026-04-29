from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class RunConfig:
    seed: int
    output_dir: Path
    overwrite: bool = False
    previous_manifest: Path | None = None


@dataclass(frozen=True)
class ExperimentalRoundConfig:
    protein_name: str
    round_name: str
    embeddings_base_path: Path
    embeddings_file_name: str
    round_base_path: Path
    round_file_names: list[str]
    wt_fasta_path: Path
    rename_WT: bool = False
    number_of_variants: int = 12


@dataclass(frozen=True)
class ModelConfig:
    regression_type: str = "randomforest"
    seed: int = 0
    params: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class EvolveProConfig:
    schema_version: int
    workflow: Literal["experimental_round"]
    run: RunConfig
    experiment: ExperimentalRoundConfig
    model: ModelConfig
    config_path: Path | None = None


@dataclass(frozen=True)
class FileHash:
    path: str
    sha256: str
    size: int


@dataclass(frozen=True)
class RoundLock:
    schema_version: int
    workflow: Literal["experimental_round"]
    config_sha256: str
    inputs: list[FileHash]
    model: ModelConfig
    previous_manifest_path: str | None = None
    previous_manifest_sha256: str | None = None
