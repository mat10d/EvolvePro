from __future__ import annotations

import importlib
import sys
from types import SimpleNamespace

import pytest


def _fake_torch(
    *,
    cuda_available: bool,
    device_names: list[str] | None = None,
    torch_version: str = "2.5.0",
    cuda_version: str | None = "12.4",
):
    names = device_names or []

    class FakeCuda:
        @staticmethod
        def is_available() -> bool:
            return cuda_available

        @staticmethod
        def device_count() -> int:
            return len(names)

        @staticmethod
        def get_device_name(index: int) -> str:
            return names[index]

    return SimpleNamespace(
        __version__=torch_version,
        version=SimpleNamespace(cuda=cuda_version),
        cuda=FakeCuda,
    )


def test_importing_embedding_device_does_not_require_torch(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    sys.modules.pop("evolvepro.embeddings.device", None)

    module = importlib.import_module("evolvepro.embeddings.device")

    assert module.EmbeddingDeviceError is not None
    assert module.resolve_device is not None


def test_auto_with_cuda_available_resolves_cuda_zero(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=True, device_names=["GPU 0", "GPU 1"])
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    metadata = device.resolve_device("auto")

    assert metadata == {
        "requested_device": "auto",
        "resolved_device": "cuda:0",
        "cuda_available": True,
        "cuda_device_count": 2,
        "cuda_device_name": "GPU 0",
        "torch_version": "2.5.0",
        "cuda_version": "12.4",
    }


def test_auto_with_no_cuda_and_no_fallback_raises(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=False)
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    with pytest.raises(device.EmbeddingDeviceError, match="CUDA is not available"):
        device.resolve_device("auto")


def test_auto_with_no_cuda_and_fallback_resolves_cpu(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=False, cuda_version=None)
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    metadata = device.resolve_device("auto", allow_cpu_fallback=True)

    assert metadata == {
        "requested_device": "auto",
        "resolved_device": "cpu",
        "cuda_available": False,
        "cuda_device_count": 0,
        "cuda_device_name": None,
        "torch_version": "2.5.0",
        "cuda_version": None,
    }


def test_cuda_with_no_cuda_raises(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=False)
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    with pytest.raises(device.EmbeddingDeviceError, match="CUDA is not available"):
        device.resolve_device("cuda")


def test_cuda_index_validates_device_count_and_returns_name(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=True, device_names=["GPU 0", "GPU 1"])
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    metadata = device.resolve_device("cuda:1")

    assert metadata["resolved_device"] == "cuda:1"
    assert metadata["cuda_device_count"] == 2
    assert metadata["cuda_device_name"] == "GPU 1"


def test_out_of_range_cuda_index_raises(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=True, device_names=["GPU 0"])
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    with pytest.raises(device.EmbeddingDeviceError, match="out of range"):
        device.resolve_device("cuda:1")


def test_explicit_cpu_without_fallback_raises(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=True, device_names=["GPU 0"])
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    with pytest.raises(device.EmbeddingDeviceError, match="CPU embedding extraction"):
        device.resolve_device("cpu")


def test_explicit_cpu_with_fallback_resolves_cpu(monkeypatch):
    from evolvepro.embeddings import device

    fake_torch = _fake_torch(cuda_available=True, device_names=["GPU 0"])
    monkeypatch.setattr(device, "_load_torch", lambda: fake_torch)

    metadata = device.resolve_device("cpu", allow_cpu_fallback=True)

    assert metadata["requested_device"] == "cpu"
    assert metadata["resolved_device"] == "cpu"
    assert metadata["cuda_available"] is True
    assert metadata["cuda_device_count"] == 1
    assert metadata["cuda_device_name"] is None


def test_unknown_device_raises(monkeypatch):
    from evolvepro.embeddings import device

    monkeypatch.setattr(
        device,
        "_load_torch",
        lambda: pytest.fail("unknown devices should fail before loading torch"),
    )

    with pytest.raises(
        device.EmbeddingDeviceError, match="Unsupported embedding device"
    ):
        device.resolve_device("mps")
