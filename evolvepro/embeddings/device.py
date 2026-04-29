from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any


class EmbeddingDeviceError(RuntimeError):
    """Raised when an embedding device cannot be resolved."""


def _load_torch() -> ModuleType:
    try:
        return importlib.import_module("torch")
    except ImportError as exc:
        raise EmbeddingDeviceError(
            "Torch is required to resolve embedding devices. Install torch with CUDA "
            "support or reuse existing embeddings."
        ) from exc


def _cuda_device_name(torch: Any, index: int) -> str | None:
    name = torch.cuda.get_device_name(index)
    return str(name) if name is not None else None


def _optional_string(value: object) -> str | None:
    return None if value is None else str(value)


def _device_metadata(
    torch: Any,
    *,
    requested_device: str,
    resolved_device: str,
    cuda_available: bool,
    cuda_device_count: int,
    cuda_device_name: str | None,
) -> dict[str, object]:
    return {
        "requested_device": requested_device,
        "resolved_device": resolved_device,
        "cuda_available": bool(cuda_available),
        "cuda_device_count": int(cuda_device_count),
        "cuda_device_name": cuda_device_name,
        "torch_version": _optional_string(getattr(torch, "__version__", None)),
        "cuda_version": _optional_string(
            getattr(getattr(torch, "version", None), "cuda", None)
        ),
    }


def _parse_requested_device(requested_device: str) -> int | None:
    if requested_device == "cuda":
        return 0
    if requested_device.startswith("cuda:"):
        index_text = requested_device.removeprefix("cuda:")
        if not index_text.isdecimal():
            raise EmbeddingDeviceError(
                f"Unsupported embedding device '{requested_device}'. Expected one of: "
                "auto, cuda, cuda:N, cpu."
            )
        return int(index_text)
    return None


def resolve_device(
    requested_device: str, *, allow_cpu_fallback: bool = False
) -> dict[str, object]:
    if requested_device not in {
        "auto",
        "cuda",
        "cpu",
    } and not requested_device.startswith("cuda:"):
        raise EmbeddingDeviceError(
            f"Unsupported embedding device '{requested_device}'. Expected one of: "
            "auto, cuda, cuda:N, cpu."
        )

    torch = _load_torch()
    cuda_available = bool(torch.cuda.is_available())
    cuda_device_count = int(torch.cuda.device_count())

    if requested_device == "auto":
        if cuda_available and cuda_device_count > 0:
            return _device_metadata(
                torch,
                requested_device=requested_device,
                resolved_device="cuda:0",
                cuda_available=cuda_available,
                cuda_device_count=cuda_device_count,
                cuda_device_name=_cuda_device_name(torch, 0),
            )
        if allow_cpu_fallback:
            return _device_metadata(
                torch,
                requested_device=requested_device,
                resolved_device="cpu",
                cuda_available=cuda_available,
                cuda_device_count=cuda_device_count,
                cuda_device_name=None,
            )
        raise EmbeddingDeviceError(
            "CUDA is not available for GPU-first embedding extraction. Set "
            "allow_cpu_fallback=True to use CPU compatibility mode."
        )

    if requested_device == "cpu":
        if not allow_cpu_fallback:
            raise EmbeddingDeviceError(
                "CPU embedding extraction is compatibility mode. Set "
                "allow_cpu_fallback=True to request cpu explicitly."
            )
        return _device_metadata(
            torch,
            requested_device=requested_device,
            resolved_device="cpu",
            cuda_available=cuda_available,
            cuda_device_count=cuda_device_count,
            cuda_device_name=None,
        )

    cuda_index = _parse_requested_device(requested_device)
    if cuda_index is None:
        raise EmbeddingDeviceError(
            f"Unsupported embedding device '{requested_device}'. Expected one of: "
            "auto, cuda, cuda:N, cpu."
        )
    if not cuda_available:
        raise EmbeddingDeviceError(
            f"CUDA is not available for requested embedding device '{requested_device}'."
        )
    if cuda_index >= cuda_device_count:
        raise EmbeddingDeviceError(
            f"Requested embedding device '{requested_device}' is out of range; "
            f"CUDA device count is {cuda_device_count}."
        )

    return _device_metadata(
        torch,
        requested_device=requested_device,
        resolved_device=f"cuda:{cuda_index}",
        cuda_available=cuda_available,
        cuda_device_count=cuda_device_count,
        cuda_device_name=_cuda_device_name(torch, cuda_index),
    )
