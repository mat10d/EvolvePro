import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta(file_path):
    """Read a FASTA file and return a dictionary of sequences."""
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}

def compare_fasta_files(dir1, dir2):
    """Compare FASTA files in two directories."""
    fasta_files = [f for f in os.listdir(dir1) if f.endswith('.fasta') or f.endswith('.fa')]

    for fasta_file in fasta_files:
        file1_path = os.path.join(dir1, fasta_file)
        file2_path = os.path.join(dir2, fasta_file)

        if not os.path.exists(file2_path):
            print(f"File {fasta_file} not found in {dir2}")
            continue

        sequences1 = read_fasta(file1_path)
        sequences2 = read_fasta(file2_path)

        compare_sequences(fasta_file, sequences1, sequences2)

def compare_sequences(file_name, sequences1, sequences2):
    """Compare two sets of sequences."""
    if sequences1 == sequences2:
        print(f"{file_name}: Identical")
        return

    print(f"{file_name}: Different")
    
    # Compare sequence IDs
    ids1 = set(sequences1.keys())
    ids2 = set(sequences2.keys())
    
    if ids1 != ids2:
        print("  Sequence IDs are different")
        print(f"  Sequences only in file1: {ids1 - ids2}")
        print(f"  Sequences only in file2: {ids2 - ids1}")
    
    # Compare sequences for common IDs
    common_ids = ids1.intersection(ids2)
    for seq_id in common_ids:
        if sequences1[seq_id] != sequences2[seq_id]:
            print(f"  Sequence '{seq_id}' differs:")
            print(f"    File1: {sequences1[seq_id][:50]}...")  # Show first 50 characters
            print(f"    File2: {sequences2[seq_id][:50]}...")

# Usage
dir1 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp/'
dir2 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/data_processing/output/'

compare_fasta_files(dir1, dir2)
