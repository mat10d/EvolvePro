from __future__ import annotations

import ast
import importlib.util
from pathlib import Path


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dms_main_parser_contract(repo_root, data_regression):
    module = _load_module(repo_root / "scripts" / "dms" / "dms_main.py", "dms_main")
    parser = module.create_parser()
    option_dests = sorted(
        action.dest
        for action in parser._actions
        if action.dest != "help"  # noqa: SLF001
    )

    args = parser.parse_args(
        [
            "--dataset_name",
            "kelsic",
            "--experiment_name",
            "smoke",
            "--model_name",
            "esm_smoke",
            "--embeddings_path",
            "embeddings",
            "--labels_path",
            "labels",
            "--num_simulations",
            "1",
            "--num_iterations",
            "2",
            "--measured_var",
            "activity",
            "--learning_strategies",
            "topn",
            "--num_mutants_per_round",
            "8",
            "--num_final_round_mutants",
            "8",
            "--first_round_strategies",
            "random",
            "--embedding_types",
            "embeddings",
            "--regression_types",
            "ridge",
            "--embeddings_file_type",
            "csv",
            "--output_dir",
            "out",
        ]
    )

    data_regression.check(
        {
            "option_dests": option_dests,
            "parsed": vars(args),
        }
    )


def test_experimental_scripts_static_contract(repo_root, data_regression):
    scripts = sorted((repo_root / "scripts" / "exp").glob("*.py"))
    payload = {}
    for script in scripts:
        tree = ast.parse(script.read_text())
        calls = sorted(
            {
                node.func.id
                for node in ast.walk(tree)
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            }
        )
        imports = [
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and node.module == "evolvepro.src.evolve"
            for alias in node.names
        ]
        payload[script.name] = {"imports": imports, "top_level_calls": calls}

    data_regression.check(payload)
