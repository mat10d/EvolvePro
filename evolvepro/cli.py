from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from evolvepro.io.config import evolve_experimental_from_config, load_config
from evolvepro.io.manifest import validate_round_directory
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

    return parser


def main(argv: Sequence[str] | None = None, *, cwd: Path | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    base = (cwd or Path.cwd()).expanduser()

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "init":
        workspace = Path(args.directory).expanduser()
        if not workspace.is_absolute():
            workspace = base / workspace
        init_workspace(workspace, force=args.force)
        print(f"initialized workspace: {workspace}")
        return 0

    if args.command == "run":
        config_path = Path(args.config).expanduser()
        if not config_path.is_absolute():
            config_path = base / config_path
        config = load_config(config_path)
        evolve_experimental_from_config(config)
        print(f"completed run from config: {config_path}")
        return 0

    if args.command == "validate":
        round_dir = Path(args.round_dir).expanduser()
        if not round_dir.is_absolute():
            round_dir = base / round_dir
        validate_round_directory(round_dir)
        print(f"validated round: {round_dir}")
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
