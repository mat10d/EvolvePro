from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from evolvepro.io.config import evolve_experimental_from_config, load_config
from evolvepro.io.manifest import validate_round_directory
from evolvepro.src.process import generate_single_aa_mutants, generate_wt
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
        parser.error("process requires a subcommand")

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
