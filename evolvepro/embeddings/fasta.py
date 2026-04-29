from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

_SINGLE_SUBSTITUTION_RE = re.compile(r"^([A-Z])([1-9][0-9]*)([A-Z])$")


def read_single_fasta_sequence(path: Path) -> tuple[str, str]:
    record_id: str | None = None
    sequence_parts: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if line == "":
                continue
            if line.startswith(">"):
                if record_id is not None:
                    raise ValueError(f"multiple FASTA records in {path}")
                header = line[1:].strip()
                if header == "":
                    raise ValueError(f"FASTA record id must not be empty in {path}")
                record_id = header.split()[0]
                continue
            if record_id is None:
                raise ValueError(
                    f"sequence data before FASTA header in {path} at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))

    if record_id is None:
        raise ValueError(f"empty FASTA: {path}")

    sequence = "".join(sequence_parts)
    if sequence == "":
        raise ValueError(f"empty sequence for FASTA record {record_id!r} in {path}")

    return record_id, sequence


def apply_single_substitution(wt_sequence: str, variant: str) -> str:
    match = _SINGLE_SUBSTITUTION_RE.fullmatch(variant)
    if match is None:
        raise ValueError(
            f"malformed variant {variant!r}; expected single substitution like A1C"
        )

    wt_aa, position_text, mutant_aa = match.groups()
    position = int(position_text)
    index = position - 1

    if index >= len(wt_sequence):
        raise ValueError(
            f"variant {variant!r} position {position} is out of range for WT "
            f"sequence length {len(wt_sequence)}"
        )

    if wt_aa == mutant_aa:
        raise ValueError(
            f"same-aa substitution is not a mutation candidate: {variant!r}"
        )

    actual_wt_aa = wt_sequence[index]
    if actual_wt_aa != wt_aa:
        raise ValueError(
            f"variant {variant!r} expected WT residue {wt_aa} at position "
            f"{position}, found {actual_wt_aa}"
        )

    return wt_sequence[:index] + mutant_aa + wt_sequence[index + 1 :]


def embedding_fasta_records(
    wt_sequence: str, variants: list[str]
) -> list[tuple[str, str]]:
    counts = Counter(variants)
    duplicates = sorted(variant for variant, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(f"duplicate variants: {duplicates}")

    return [
        (variant, apply_single_substitution(wt_sequence, variant))
        for variant in sorted(variants)
    ]


def write_embedding_fasta(path: Path, wt_sequence: str, variants: list[str]) -> Path:
    records = embedding_fasta_records(wt_sequence, variants)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n{sequence}\n")

    return path
