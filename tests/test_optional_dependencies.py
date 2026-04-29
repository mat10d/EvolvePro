from __future__ import annotations

import importlib
import sys

from evolvepro.workspace import init_workspace


def test_data_module_import_does_not_require_torch(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    sys.modules.pop("evolvepro.src.data", None)

    module = importlib.import_module("evolvepro.src.data")

    assert module.load_dms_data is not None
    assert module.load_experimental_data is not None


def test_workspace_config_load_does_not_require_torch(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    sys.modules.pop("evolvepro.io.config", None)
    init_workspace(tmp_path)

    module = importlib.import_module("evolvepro.io.config")
    config = module.load_config(tmp_path / "evolvepro.toml")

    assert config.embedding.backend == "esm"
    assert config.embedding.device == "auto"
