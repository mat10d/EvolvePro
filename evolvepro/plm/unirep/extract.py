import argparse
import time
from pathlib import Path
import pandas as pd
from Bio import SeqIO

from jax_unirep import get_reps

def create_arg_parser():
    """Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(description="Arguments to extract UniRep embeddings from fasta file and save as CSV.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to input FASTA file.")
    parser.add_argument("-o", "--output", required=True, type=str, help="Path to output CSV file.")
    return parser

def read_fasta(fasta_path):
    sequences = {}
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_id = record.id
            sequence = str(record.seq).upper()
            sequences[sequence_id] = sequence
    return sequences

def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    
    seq_path = Path(args.input)
    csv_path = Path(args.output)
    
    # Read FASTA file
    seqs = read_fasta(seq_path)
    sequences = list(seqs.values())
    sequence_ids = list(seqs.keys())
    
    # Get UniRep embeddings
    h_avg, h_final, c_final = get_reps(sequences)

    # Print information about the last embedded protein
    last_protein_id = sequence_ids[-1]
    last_protein_length = len(seqs[last_protein_id])
    last_global_representation_shape = h_avg[-1].shape

    print(f"Embedded protein {last_protein_id} with length {last_protein_length} to emb. of global shape: {last_global_representation_shape}")
    
    # Create a dictionary to hold the embeddings
    embeddings_dict = {seq_id: emb for seq_id, emb in zip(sequence_ids, h_avg)}
    
    # Convert embeddings to DataFrame
    df = pd.DataFrame.from_dict(embeddings_dict, orient='index')
    
    # Save embeddings to CSV
    df.to_csv(csv_path, index_label='sequence_id')
    print(f"Embeddings saved to {csv_path}")

if __name__ == '__main__':
    main()