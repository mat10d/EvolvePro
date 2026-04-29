from __future__ import annotations

import pytest

from evolvepro.embeddings.fasta import (
    apply_single_substitution,
    read_single_fasta_sequence,
    write_embedding_fasta,
)


def test_read_single_fasta_sequence_joins_wrapped_lines(tmp_path):
    fasta_path = tmp_path / "wt.fasta"
    fasta_path.write_text("  \n>WT\n AC\nDE \n\nFG\n", encoding="utf-8")

    record_id, sequence = read_single_fasta_sequence(fasta_path)

    assert record_id == "WT"
    assert sequence == "ACDEFG"


def test_read_single_fasta_sequence_rejects_multiple_records(tmp_path):
    fasta_path = tmp_path / "wt.fasta"
    fasta_path.write_text(">one\nAC\n>two\nDE\n", encoding="utf-8")

    with pytest.raises(ValueError, match="multiple FASTA records"):
        read_single_fasta_sequence(fasta_path)


def test_write_embedding_fasta_is_byte_stable_and_sorted(tmp_path):
    fasta_path = tmp_path / "nested" / "embedding.fasta"

    result = write_embedding_fasta(fasta_path, "AC", ["A1D", "A1C"])

    assert result == fasta_path
    assert fasta_path.read_bytes() == b">A1C\nCC\n>A1D\nDC\n"

    write_embedding_fasta(fasta_path, "AC", ["A1C", "A1D"])
    assert fasta_path.read_bytes() == b">A1C\nCC\n>A1D\nDC\n"


def test_apply_single_substitution_uses_one_based_positions():
    assert apply_single_substitution("MAG", "A2D") == "MDG"


def test_apply_single_substitution_rejects_mismatched_wt_residue():
    with pytest.raises(ValueError, match="expected WT residue A at position 2"):
        apply_single_substitution("MG", "A2D")


def test_write_embedding_fasta_rejects_duplicate_variants(tmp_path):
    with pytest.raises(ValueError, match="duplicate variants"):
        write_embedding_fasta(tmp_path / "embedding.fasta", "AC", ["A1C", "A1C"])


def test_apply_single_substitution_rejects_malformed_variant():
    with pytest.raises(ValueError, match="malformed variant"):
        apply_single_substitution("AC", "A0C")
