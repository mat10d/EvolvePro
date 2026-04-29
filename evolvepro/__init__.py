from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("evolvepro")
except PackageNotFoundError:  # pragma: no cover - editable metadata is present in tests
    __version__ = "0.0.0"
