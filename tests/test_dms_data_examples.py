from __future__ import annotations

from Bio import SeqIO

from evolvepro.src.process import process_dataset


def test_process_sarscov2_preprocessed_dms(tmp_path, repo_root, data_regression):
    output_dir = tmp_path / "dms"
    df, fractions = process_dataset(
        file_path=str(
            repo_root / "data" / "dms" / "activity" / "sarscov2_preprocessed.csv"
        ),
        wt_fasta_path=str(repo_root / "data" / "dms" / "wt_fasta" / "cov2_S_WT.fasta"),
        dataset_name="sarscov2_test",
        activity_column="mut_escape",
        cutoff_value=0.1,
        output_dir=str(output_dir),
        cutoff_rule="greater_than",
        drop_columns=True,
    )
    fasta_records = list(SeqIO.parse(output_dir / "sarscov2_test.fasta", "fasta"))

    data_regression.check(
        {
            "shape": list(df.shape),
            "columns": df.columns.tolist(),
            "fractions": [round(float(value), 8) for value in fractions],
            "activity_min": round(float(df["activity"].min()), 8),
            "activity_max": round(float(df["activity"].max()), 8),
            "activity_binary_counts": {
                str(key): int(value)
                for key, value in df["activity_binary"]
                .value_counts()
                .sort_index()
                .items()
            },
            "head": df.head(5).round(8).to_dict(orient="records"),
            "tail": df.tail(5).round(8).to_dict(orient="records"),
            "fasta_count": len(fasta_records),
            "fasta_head_ids": [record.id for record in fasta_records[:5]],
        }
    )
