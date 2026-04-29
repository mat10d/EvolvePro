from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

import pandas as pd
import tomli

from evolvepro.io.canonical import canonical_json_sha256
from evolvepro.io.schemas import (
    EvolveProConfig,
    ExperimentalRoundConfig,
    ModelConfig,
    RunConfig,
)
from evolvepro.src.evolve import evolve_experimental


def _require_section(data: dict[str, Any], name: str) -> dict[str, Any]:
    value = data.get(name)
    if value is None:
        raise ValueError(f"missing required section: {name}")
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a table")
    return value


def _require(data: dict[str, Any], key: str, expected_type: type, prefix: str) -> Any:
    if key not in data:
        raise ValueError(f"missing required field: {prefix}.{key}")
    value = data[key]
    if not isinstance(value, expected_type):
        raise ValueError(f"{prefix}.{key} must be {expected_type.__name__}")
    return value


def _optional_bool(data: dict[str, Any], key: str, default: bool, prefix: str) -> bool:
    value = data.get(key, default)
    if not isinstance(value, bool):
        raise ValueError(f"{prefix}.{key} must be bool")
    return value


def _optional_int(data: dict[str, Any], key: str, default: int, prefix: str) -> int:
    value = data.get(key, default)
    if not isinstance(value, int):
        raise ValueError(f"{prefix}.{key} must be int")
    return value


def _path(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return base / path


def _require_string_list(data: dict[str, Any], key: str, prefix: str) -> list[str]:
    if key not in data:
        raise ValueError(f"missing required field: {prefix}.{key}")
    value = data[key]
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{prefix}.{key} must be list[str]")
    return value


def load_config(path: Path) -> EvolveProConfig:
    config_path = path.resolve()
    with config_path.open("rb") as handle:
        data = tomli.load(handle)

    schema_version = _require(data, "schema_version", int, "root")
    workflow = _require(data, "workflow", str, "root")
    if workflow != "experimental_round":
        raise ValueError(f"unsupported workflow: {workflow}")

    run_data = _require_section(data, "run")
    experiment_data = _require_section(data, "experiment")
    inputs_data = _require_section(data, "inputs")
    model_data = _require_section(data, "model")
    base = config_path.parent

    run = RunConfig(
        seed=_require(run_data, "seed", int, "run"),
        output_dir=_path(base, _require(run_data, "output_dir", str, "run")),
        overwrite=_optional_bool(run_data, "overwrite", False, "run"),
        previous_manifest=(
            _path(base, run_data["previous_manifest"])
            if "previous_manifest" in run_data
            else None
        ),
    )

    experiment = ExperimentalRoundConfig(
        protein_name=_require(experiment_data, "protein_name", str, "experiment"),
        round_name=_require(experiment_data, "round_name", str, "experiment"),
        embeddings_base_path=_path(
            base, _require(inputs_data, "embeddings_base_path", str, "inputs")
        ),
        embeddings_file_name=_require(
            inputs_data, "embeddings_file_name", str, "inputs"
        ),
        round_base_path=_path(
            base, _require(inputs_data, "round_base_path", str, "inputs")
        ),
        round_file_names=_require_string_list(
            inputs_data, "round_file_names", "inputs"
        ),
        wt_fasta_path=_path(
            base, _require(inputs_data, "wt_fasta_path", str, "inputs")
        ),
        rename_WT=_optional_bool(experiment_data, "rename_WT", False, "experiment"),
        number_of_variants=_optional_int(
            experiment_data, "number_of_variants", 12, "experiment"
        ),
    )

    model = ModelConfig(
        regression_type=model_data.get("regression_type", "randomforest"),
        seed=_optional_int(model_data, "seed", run.seed, "model"),
        params={
            key: value
            for key, value in model_data.items()
            if key not in {"regression_type", "seed"}
        },
    )
    if not isinstance(model.regression_type, str):
        raise ValueError("model.regression_type must be str")

    return EvolveProConfig(
        schema_version=schema_version,
        workflow="experimental_round",
        run=run,
        experiment=experiment,
        model=model,
        config_path=config_path,
    )


def config_payload(config: EvolveProConfig) -> dict[str, Any]:
    payload = asdict(replace(config, config_path=None))
    payload.pop("config_path", None)
    return payload


def config_sha256(config: EvolveProConfig) -> str:
    return canonical_json_sha256(config_payload(config))


def _toml_string(value: str) -> str:
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _toml_path(value: Path) -> str:
    return _toml_string(value.as_posix())


def _toml_list(values: list[str]) -> str:
    return "[" + ", ".join(_toml_string(value) for value in values) + "]"


def _toml_bool(value: bool) -> str:
    return "true" if value else "false"


def write_config(config: EvolveProConfig, path: Path) -> None:
    lines = [
        f"schema_version = {config.schema_version}",
        f"workflow = {_toml_string(config.workflow)}",
        "",
        "[run]",
        f"seed = {config.run.seed}",
        f"output_dir = {_toml_path(config.run.output_dir)}",
        f"overwrite = {_toml_bool(config.run.overwrite)}",
    ]
    if config.run.previous_manifest is not None:
        lines.append(f"previous_manifest = {_toml_path(config.run.previous_manifest)}")
    lines.extend(
        [
            "",
            "[experiment]",
            f"protein_name = {_toml_string(config.experiment.protein_name)}",
            f"round_name = {_toml_string(config.experiment.round_name)}",
            f"rename_WT = {_toml_bool(config.experiment.rename_WT)}",
            f"number_of_variants = {config.experiment.number_of_variants}",
            "",
            "[inputs]",
            f"embeddings_base_path = {_toml_path(config.experiment.embeddings_base_path)}",
            f"embeddings_file_name = {_toml_string(config.experiment.embeddings_file_name)}",
            f"round_base_path = {_toml_path(config.experiment.round_base_path)}",
            f"round_file_names = {_toml_list(config.experiment.round_file_names)}",
            f"wt_fasta_path = {_toml_path(config.experiment.wt_fasta_path)}",
            "",
            "[model]",
            f"regression_type = {_toml_string(config.model.regression_type)}",
            f"seed = {config.model.seed}",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def evolve_experimental_from_config(
    config: EvolveProConfig,
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    experiment = config.experiment
    return evolve_experimental(
        protein_name=experiment.protein_name,
        round_name=experiment.round_name,
        embeddings_base_path=str(experiment.embeddings_base_path),
        embeddings_file_name=experiment.embeddings_file_name,
        round_base_path=str(experiment.round_base_path),
        round_file_names=experiment.round_file_names,
        wt_fasta_path=str(experiment.wt_fasta_path),
        rename_WT=experiment.rename_WT,
        number_of_variants=experiment.number_of_variants,
        output_dir=str(config.run.output_dir),
    )
