from __future__ import annotations

import builtins
import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest

from evolvepro.embeddings.schemas import EmbeddingRequest


MODULE_NAME = "evolvepro.embeddings.backends.esm"


def test_importing_esm_backend_does_not_require_torch_or_esm(
    monkeypatch: pytest.MonkeyPatch,
):
    sys.modules.pop(MODULE_NAME, None)
    real_import = builtins.__import__

    def guarded_import(
        name: str,
        globals_: dict[str, Any] | None = None,
        locals_: dict[str, Any] | None = None,
        fromlist: tuple[str, ...] = (),
        level: int = 0,
    ) -> Any:
        if (
            name == "torch"
            or name.startswith("torch.")
            or name == "esm"
            or name.startswith("esm.")
        ):
            raise AssertionError(f"unexpected eager import: {name}")
        return real_import(name, globals_, locals_, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    module = importlib.import_module(MODULE_NAME)

    assert hasattr(module, "EsmEmbeddingBackend")


def test_extract_passes_expected_args_to_legacy_runner(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    esm_backend = importlib.import_module(MODULE_NAME)
    calls = _install_fake_legacy(monkeypatch, esm_backend)
    request = _request(tmp_path, device="cuda:0")

    result = esm_backend.EsmEmbeddingBackend().extract(request)

    assert len(calls) == 1
    args = calls[0]
    assert args.model_location == "esm2_t6_8M_UR50D"
    assert args.fasta_file == tmp_path / "input.fasta"
    assert args.output_dir == tmp_path / "embeddings" / "pt"
    assert args.toks_per_batch == 128
    assert args.repr_layers == [-1]
    assert args.include == ["mean"]
    assert args.truncation_seq_length == 1022
    assert args.nogpu is False
    assert args.concatenate_dir == tmp_path / "embeddings"
    assert result.embeddings_csv == tmp_path / "embeddings" / "embeddings.csv"
    assert result.metadata_json == tmp_path / "embeddings" / "metadata.json"
    assert result.backend_metadata["backend"] == "esm"
    assert result.device_metadata["requested_device"] == "cuda:0"
    assert result.device_metadata["nogpu"] is False
    assert result.artifacts == {"pt_dir": tmp_path / "embeddings" / "pt"}


@pytest.mark.parametrize(
    ("device", "expected_nogpu"),
    [
        ("cuda", False),
        ("cuda:0", False),
        ("auto", False),
        ("cpu", True),
    ],
)
def test_extract_maps_request_device_to_legacy_nogpu(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    device: str,
    expected_nogpu: bool,
):
    esm_backend = importlib.import_module(MODULE_NAME)
    calls = _install_fake_legacy(monkeypatch, esm_backend)

    esm_backend.EsmEmbeddingBackend().extract(_request(tmp_path, device=device))

    assert calls[0].nogpu is expected_nogpu


def test_metadata_json_is_written_and_deterministic(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    esm_backend = importlib.import_module(MODULE_NAME)
    _install_fake_legacy(monkeypatch, esm_backend)
    request = _request(tmp_path, device="cpu")
    backend = esm_backend.EsmEmbeddingBackend()

    result = backend.extract(request)
    first_metadata = result.metadata_json.read_text(encoding="utf-8")
    backend.extract(request)
    second_metadata = result.metadata_json.read_text(encoding="utf-8")

    assert first_metadata == second_metadata
    metadata = json.loads(first_metadata)
    assert metadata == {
        "artifacts": {
            "embeddings_csv": str(tmp_path / "embeddings" / "embeddings.csv"),
            "pt_dir": str(tmp_path / "embeddings" / "pt"),
        },
        "backend": {
            "backend": "esm",
            "legacy_wrapper": "evolvepro.plm.esm.extract",
            "model_name": "esm2_t6_8M_UR50D",
            "model_revision": "main",
            "repr_layers": [-1],
            "representation": "mean",
        },
        "device": {
            "dtype": "float32",
            "nogpu": True,
            "requested_device": "cpu",
            "seed": 7,
        },
        "dtype": "float32",
        "model": {
            "name": "esm2_t6_8M_UR50D",
            "revision": "main",
        },
        "representation": {
            "name": "mean",
            "repr_layers": [-1],
        },
        "seed": 7,
    }


def _install_fake_legacy(
    monkeypatch: pytest.MonkeyPatch, esm_backend: Any
) -> list[Any]:
    calls: list[Any] = []

    def fake_run(args: Any) -> None:
        calls.append(args)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "seq1.pt").write_text("fake legacy tensor", encoding="utf-8")

    def fake_concatenate(output_dir: Path, output_csv: Path) -> None:
        assert output_dir.exists()
        output_csv.write_text(",0\nseq1,1.0\n", encoding="utf-8")

    monkeypatch.setattr(esm_backend, "_run_legacy_esm", fake_run)
    monkeypatch.setattr(esm_backend, "_concatenate_legacy_outputs", fake_concatenate)
    return calls


def _request(tmp_path: Path, *, device: str) -> EmbeddingRequest:
    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text(">seq1\nM\n", encoding="utf-8")
    return EmbeddingRequest(
        fasta_path=fasta_path,
        output_dir=tmp_path / "embeddings",
        backend="esm",
        model_name="esm2_t6_8M_UR50D",
        model_revision="main",
        repr_layers=[-1],
        representation="mean",
        toks_per_batch=128,
        device=device,
        dtype="float32",
        seed=7,
        allow_download=False,
        allow_cpu_fallback=False,
        cache_dir=tmp_path / "cache",
    )
