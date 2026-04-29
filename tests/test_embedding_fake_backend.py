from __future__ import annotations

import builtins
import csv
import importlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from evolvepro.embeddings import EmbeddingRequest


def _write_fasta(path: Path, content: str = ">var_b\nMTEYK\n>var_a\nGAGGVG\n") -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def _request(
    tmp_path: Path, *, seed: int = 0, output_name: str = "out"
) -> EmbeddingRequest:
    return EmbeddingRequest(
        fasta_path=_write_fasta(tmp_path / "input.fasta"),
        output_dir=tmp_path / output_name,
        backend="fake",
        model_name="fake-model",
        model_revision=None,
        repr_layers=[-1],
        representation="mean",
        toks_per_batch=512,
        device="auto",
        dtype="float32",
        seed=seed,
        allow_download=False,
        allow_cpu_fallback=False,
        cache_dir=None,
    )


def _extract(request: EmbeddingRequest):
    from evolvepro.embeddings.backends.fake import FakeEmbeddingBackend

    return FakeEmbeddingBackend().extract(request)


def test_extract_writes_embeddings_and_metadata(tmp_path: Path):
    request = _request(tmp_path)

    result = _extract(request)

    assert result.embeddings_csv == request.output_dir / "embeddings.csv"
    assert result.metadata_json == request.output_dir / "metadata.json"
    assert result.embeddings_csv.is_file()
    assert result.metadata_json.is_file()
    metadata = json.loads(result.metadata_json.read_text(encoding="utf-8"))
    assert metadata == {
        "backend": "fake",
        "device": "fake",
        "dtype": "float32",
        "model_name": "fake-model",
        "model_revision": None,
        "repr_layers": [-1],
        "representation": "mean",
        "seed": 0,
    }
    assert result.backend_metadata == metadata
    assert result.device_metadata == {"device": "fake", "dtype": "float32"}


def test_output_is_byte_stable_across_repeated_runs(tmp_path: Path):
    request = _request(tmp_path)

    first = _extract(request)
    first_csv = first.embeddings_csv.read_bytes()
    first_metadata = first.metadata_json.read_bytes()
    second = _extract(request)

    assert second.embeddings_csv.read_bytes() == first_csv
    assert second.metadata_json.read_bytes() == first_metadata


def test_output_changes_when_seed_changes(tmp_path: Path):
    seed_zero = _request(tmp_path, seed=0, output_name="seed-zero")
    seed_one = replace(seed_zero, seed=1, output_dir=tmp_path / "seed-one")

    seed_zero_result = _extract(seed_zero)
    seed_one_result = _extract(seed_one)

    assert (
        seed_zero_result.embeddings_csv.read_bytes()
        != seed_one_result.embeddings_csv.read_bytes()
    )
    assert (
        seed_zero_result.metadata_json.read_bytes()
        != seed_one_result.metadata_json.read_bytes()
    )


def test_records_requested_dtype_in_metadata(tmp_path: Path):
    request = replace(_request(tmp_path), dtype="float16")

    result = _extract(request)

    metadata = json.loads(result.metadata_json.read_text(encoding="utf-8"))
    assert metadata["dtype"] == "float16"
    assert result.device_metadata == {"device": "fake", "dtype": "float16"}


def test_csv_has_variant_and_numeric_feature_columns(tmp_path: Path):
    result = _extract(_request(tmp_path))

    with result.embeddings_csv.open("r", encoding="utf-8", newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        rows = list(reader)

    assert reader.fieldnames == [
        "variant",
        "feature_0",
        "feature_1",
        "feature_2",
        "feature_3",
    ]
    assert [row["variant"] for row in rows] == ["var_a", "var_b"]
    for row in rows:
        for feature_name in reader.fieldnames[1:]:
            assert isinstance(float(row[feature_name]), float)


def test_duplicate_fasta_labels_raise(tmp_path: Path):
    request = _request(tmp_path)
    request.fasta_path.write_text(">dup\nAAAA\n>dup\nCCCC\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Duplicate FASTA label: dup"):
        _extract(request)


def test_importing_fake_backend_does_not_require_torch(monkeypatch: pytest.MonkeyPatch):
    sys.modules.pop("evolvepro.embeddings.backends.fake", None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "torch" or name.startswith("torch."):
            raise AssertionError("fake backend must not import torch")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    module = importlib.import_module("evolvepro.embeddings.backends.fake")

    assert module.FakeEmbeddingBackend is not None
