from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tests.fixtures.synthetic import make_embeddings, make_labels


KELSIC_WT_SEQUENCE = (
    "MAKEDNIEMQGTVLETLPNTMFRVELENGHVVTAHISGKMRKNYIRILTGDKVTVELTPYDLSKGRIVFRSR"
)


@pytest.fixture
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture
def golden_dir(repo_root: Path) -> Path:
    return repo_root / "tests" / "fixtures" / "golden"


@pytest.fixture
def kelsic_wt_sequence() -> str:
    return KELSIC_WT_SEQUENCE


@pytest.fixture
def kelsic_wt_fasta(tmp_path: Path, kelsic_wt_sequence: str) -> Path:
    from evolvepro.src.process import generate_wt

    output = tmp_path / "kelsic_WT.fasta"
    generate_wt(kelsic_wt_sequence, str(output))
    return output


@pytest.fixture
def synthetic_embeddings():
    def _make(
        n_variants: int, n_dim: int, seed: int, include_wt: bool = True
    ) -> pd.DataFrame:
        return make_embeddings(
            n_variants=n_variants, n_dim=n_dim, seed=seed, include_wt=include_wt
        )

    return _make


@pytest.fixture
def synthetic_labels():
    def _make(
        variants: list[str], seed: int, include_metadata: bool = True
    ) -> pd.DataFrame:
        return make_labels(
            variants=variants, seed=seed, include_metadata=include_metadata
        )

    return _make
