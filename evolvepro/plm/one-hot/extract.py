#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import os
from typing import List, Union

def create_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser."""
    parser = argparse.ArgumentParser(description="Encode protein sequences from FASTA files.")
    parser.add_argument("fasta_file", help="Path to a FASTA file containing protein sequences")
    parser.add_argument("--method", choices=['one_hot', 'integer'], required=True, help="Encoding method")
    parser.add_argument("--results_path", default="results", help="Path to save results")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    return parser

def integer_encode(sequence: str) -> List[int]:
    """Encode a protein sequence as integers."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_to_int = {aa: i for i, aa in enumerate(amino_acids)}
    return [aa_to_int[aa] for aa in sequence if aa in amino_acids]

def one_hot_encode(sequence: str) -> np.ndarray:
    """Encode a protein sequence as a one-hot array."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    one_hot = np.zeros((len(sequence), len(amino_acids)))
    for i, aa in enumerate(sequence):
        if aa in amino_acids:
            one_hot[i, aa_to_index[aa]] = 1
    return one_hot

def encode_sequences_from_fasta(fasta_file: str, method: str, results_path: str, verbose: bool = False) -> pd.DataFrame:
    """
    Encode sequences from a FASTA file using the specified method.
    
    Args:
        fasta_file (str): Path to the input FASTA file.
        method (str): Encoding method ('one_hot' or 'integer').
        results_path (str): Path to save the results.
        verbose (bool): Whether to print verbose output.
    
    Returns:
        pd.DataFrame: Encoded sequences as a DataFrame.
    """
    dataset_name = os.path.splitext(os.path.basename(fasta_file))[0]
    print(f"Processing dataset: {dataset_name}")
    encoded_sequences = []
    row_names = []
    
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            row_names.append(record.id)
            if method == 'one_hot':
                encoded_sequence = one_hot_encode(sequence).flatten()
            elif method == 'integer':
                encoded_sequence = integer_encode(sequence)
            else:
                raise ValueError(f"Unknown encoding method: {method}")
            encoded_sequences.append(encoded_sequence)
            if verbose:
                print(f"Encoded sequence: {record.id}")
    
    encoded_sequences = pd.DataFrame(np.array(encoded_sequences), index=row_names)
    print(f"Encoded shape: {encoded_sequences.shape}")
    
    output_file = os.path.join(results_path, f"{dataset_name}_{method}_encoded.csv")
    encoded_sequences.to_csv(output_file)
    print(f"Saved encoded sequences to: {output_file}")
    return encoded_sequences

def main():
    parser = create_parser()
    args = parser.parse_args()
    os.makedirs(args.results_path, exist_ok=True)
    encode_sequences_from_fasta(args.fasta_file, args.method, args.results_path, args.verbose)

if __name__ == "__main__":
    main()