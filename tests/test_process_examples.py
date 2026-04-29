from __future__ import annotations

from Bio import SeqIO

from evolvepro.src.process import (
    generate_single_aa_mutants,
    generate_wt,
    suggest_initial_mutants,
)


def test_generate_wt_kelsic(tmp_path, kelsic_wt_sequence, data_regression):
    output = tmp_path / "kelsic_WT.fasta"
    generate_wt(kelsic_wt_sequence, str(output))
    record = SeqIO.read(output, "fasta")

    data_regression.check(
        {
            "id": record.id,
            "description": record.description,
            "sequence_length": len(record.seq),
            "sequence": str(record.seq),
        }
    )


def test_generate_single_aa_mutants_kelsic_subset(
    tmp_path, kelsic_wt_fasta, data_regression
):
    output = tmp_path / "kelsic_subset.fasta"
    generate_single_aa_mutants(str(kelsic_wt_fasta), str(output), positions=[1, 2, 3])
    records = list(SeqIO.parse(output, "fasta"))

    data_regression.check(
        {
            "count": len(records),
            "ids": [record.id for record in records],
            "first_sequences": {record.id: str(record.seq) for record in records[:5]},
            "last_sequences": {record.id: str(record.seq) for record in records[-5:]},
        }
    )


def test_suggest_initial_mutants_seeded(
    tmp_path, kelsic_wt_fasta, capsys, data_regression
):
    mutants = tmp_path / "kelsic_subset.fasta"
    generate_single_aa_mutants(str(kelsic_wt_fasta), str(mutants), positions=[1, 2, 3])
    capsys.readouterr()

    suggest_initial_mutants(str(mutants), num_mutants=12, random_seed=42)
    output = capsys.readouterr().out.strip().splitlines()

    data_regression.check({"stdout": output})
