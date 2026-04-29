from __future__ import annotations

import importlib
import sys


def test_data_module_import_does_not_require_torch(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    sys.modules.pop("evolvepro.src.data", None)

    module = importlib.import_module("evolvepro.src.data")

    assert module.load_dms_data is not None
    assert module.load_experimental_data is not None
