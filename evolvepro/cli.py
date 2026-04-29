from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from evolvepro.io.config import evolve_experimental_from_config, load_config
from evolvepro.io.manifest import validate_round_directory
from evolvepro.src.evolve import grid_search
from evolvepro.src.process import (
    generate_single_aa_mutants,
    generate_wt,
    process_dataset,
)
from evolvepro.workspace import CONFIG_FILE, init_workspace


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="evolvepro",
        description="File-based EvolvePro workspace CLI.",
    )
    subparsers = parser.add_subparsers(dest="command")

    init_parser = subparsers.add_parser("init", help="initialize a workspace")
    init_parser.add_argument("directory", nargs="?", default=".")
    init_parser.add_argument("--force", action="store_true")

    run_parser = subparsers.add_parser("run", help="run from a TOML config")
    run_parser.add_argument("--config", default=CONFIG_FILE)

    validate_parser = subparsers.add_parser(
        "validate", help="validate a completed round directory"
    )
    validate_parser.add_argument("round_dir", nargs="?", default=".")

    process_parser = subparsers.add_parser("process", help="process sequence inputs")
    process_subparsers = process_parser.add_subparsers(dest="process_command")

    wt_parser = process_subparsers.add_parser("generate-wt", help="write WT FASTA")
    wt_parser.add_argument("--sequence", required=True)
    wt_parser.add_argument("--output", required=True)

    mutants_parser = process_subparsers.add_parser(
        "generate-single-aa-mutants", help="write single amino acid mutant FASTA"
    )
    mutants_parser.add_argument("--wt-fasta", required=True)
    mutants_parser.add_argument("--output", required=True)
    mutants_parser.add_argument("--positions")

    dataset_parser = process_subparsers.add_parser(
        "dataset", help="process a CSV or Excel activity dataset"
    )
    dataset_parser.add_argument("--file-path", required=True)
    dataset_parser.add_argument("--wt-fasta", required=True)
    dataset_parser.add_argument("--dataset-name", required=True)
    dataset_parser.add_argument("--activity-column", required=True)
    dataset_parser.add_argument("--cutoff-value", required=True, type=float)
    dataset_parser.add_argument("--output-dir", required=True)
    dataset_parser.add_argument("--sheet-name")
    dataset_parser.add_argument("--cutoff-rule", default="greater_than")
    dataset_parser.add_argument("--cutoff-percentiles", nargs="*")
    dataset_parser.add_argument("--aa-shift", type=int)
    dataset_parser.add_argument("--drop-columns", action="store_true")

    dms_parser = subparsers.add_parser("dms", help="run DMS workflows")
    dms_subparsers = dms_parser.add_subparsers(dest="dms_command")
    grid_parser = dms_subparsers.add_parser("grid-search", help="run a DMS grid search")
    grid_parser.add_argument("--dataset-name", required=True)
    grid_parser.add_argument("--experiment-name", required=True)
    grid_parser.add_argument("--model-name", required=True)
    grid_parser.add_argument("--embeddings-path", required=True)
    grid_parser.add_argument("--labels-path", required=True)
    grid_parser.add_argument("--num-simulations", required=True, type=int)
    grid_parser.add_argument("--num-iterations", required=True, type=int, nargs="+")
    grid_parser.add_argument("--measured-var", required=True, nargs="+")
    grid_parser.add_argument("--learning-strategies", required=True, nargs="+")
    grid_parser.add_argument(
        "--num-mutants-per-round", required=True, type=int, nargs="+"
    )
    grid_parser.add_argument("--num-final-round-mutants", required=True, type=int)
    grid_parser.add_argument("--first-round-strategies", required=True, nargs="+")
    grid_parser.add_argument("--embedding-types", required=True, nargs="+")
    grid_parser.add_argument("--pca-components", type=int, nargs="*")
    grid_parser.add_argument("--regression-types", required=True, nargs="+")
    grid_parser.add_argument("--embeddings-file-type", required=True)
    grid_parser.add_argument("--output-dir", required=True)
    grid_parser.add_argument("--embeddings-type-pt")

    return parser


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return base / path


def _parse_positions(value: str | None) -> list[int] | None:
    if value is None or value == "":
        return None
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def _parse_number(value: str) -> int | float:
    number = float(value)
    if number.is_integer():
        return int(number)
    return number


def _parse_number_list(values: list[str] | None) -> list[int | float] | None:
    if values is None:
        return None
    return [_parse_number(value) for value in values]


def main(argv: Sequence[str] | None = None, *, cwd: Path | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    base = (cwd or Path.cwd()).expanduser()

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "init":
        workspace = _resolve_path(base, args.directory)
        init_workspace(workspace, force=args.force)
        print(f"initialized workspace: {workspace}")
        return 0

    if args.command == "run":
        config_path = _resolve_path(base, args.config)
        config = load_config(config_path)
        evolve_experimental_from_config(config)
        print(f"completed run from config: {config_path}")
        return 0

    if args.command == "validate":
        round_dir = _resolve_path(base, args.round_dir)
        validate_round_directory(round_dir)
        print(f"validated round: {round_dir}")
        return 0

    if args.command == "process":
        if args.process_command == "generate-wt":
            output = _resolve_path(base, args.output)
            generate_wt(args.sequence, str(output))
            print(f"wrote WT FASTA: {output}")
            return 0
        if args.process_command == "generate-single-aa-mutants":
            wt_fasta = _resolve_path(base, args.wt_fasta)
            output = _resolve_path(base, args.output)
            generate_single_aa_mutants(
                str(wt_fasta), str(output), positions=_parse_positions(args.positions)
            )
            print(f"wrote single amino acid mutant FASTA: {output}")
            return 0
        if args.process_command == "dataset":
            file_path = _resolve_path(base, args.file_path)
            wt_fasta = _resolve_path(base, args.wt_fasta)
            output_dir = _resolve_path(base, args.output_dir)
            process_dataset(
                file_path=str(file_path),
                wt_fasta_path=str(wt_fasta),
                dataset_name=args.dataset_name,
                activity_column=args.activity_column,
                cutoff_value=args.cutoff_value,
                output_dir=str(output_dir),
                sheet_name=args.sheet_name,
                cutoff_rule=args.cutoff_rule,
                cutoff_percentiles=_parse_number_list(args.cutoff_percentiles),
                AA_shift=args.aa_shift,
                drop_columns=args.drop_columns,
            )
            print(f"processed dataset: {output_dir}")
            return 0
        parser.error("process requires a subcommand")

    if args.command == "dms":
        if args.dms_command == "grid-search":
            embeddings_path = _resolve_path(base, args.embeddings_path)
            labels_path = _resolve_path(base, args.labels_path)
            output_dir = _resolve_path(base, args.output_dir)
            grid_search(
                dataset_name=args.dataset_name,
                experiment_name=args.experiment_name,
                model_name=args.model_name,
                embeddings_path=str(embeddings_path),
                labels_path=str(labels_path),
                num_simulations=args.num_simulations,
                num_iterations=args.num_iterations,
                measured_var=args.measured_var,
                learning_strategies=args.learning_strategies,
                num_mutants_per_round=args.num_mutants_per_round,
                num_final_round_mutants=args.num_final_round_mutants,
                first_round_strategies=args.first_round_strategies,
                embedding_types=args.embedding_types,
                pca_components=args.pca_components,
                regression_types=args.regression_types,
                embeddings_file_type=args.embeddings_file_type,
                output_dir=str(output_dir),
                embeddings_type_pt=args.embeddings_type_pt,
            )
            print(f"completed DMS grid search: {output_dir}")
            return 0
        parser.error("dms requires a subcommand")

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
