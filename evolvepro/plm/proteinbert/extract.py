import argparse
import time
from pathlib import Path
import torch
import pandas as pd
import os
import sys
from Bio import SeqIO

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../../.."))
sys.path.append(project_root)

# Create a models directory next to extract.py
model_path = os.path.join(current_dir, "models")
os.makedirs(model_path, exist_ok=True)

from external.proteinbert.proteinbert.existing_model_loading import load_pretrained_model
from external.proteinbert.proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

def create_arg_parser():
    """Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(description="Arguments to extract ProteinBERT global embeddings from fasta file and save as CSV.")
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
    batch_size = 50
    
    # Read FASTA file
    seqs = read_fasta(seq_path)
    sequences = list(seqs.values())
    sequence_ids = list(seqs.keys())
    max_seq_len = max(len(seq) for seq in sequences)
    
    # Calculate seq_len based on the condition
    seq_len = 512 if max_seq_len <= 512 else max_seq_len + 2
    
    pretrained_model_generator, input_encoder = load_pretrained_model(local_model_dump_dir=model_path, validate_downloading=False)
    model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))
    encoded_x = input_encoder.encode_X(sequences, seq_len)
    local_representations, global_representations = model.predict(encoded_x, batch_size=batch_size)

    # Print information about the last embedded protein
    last_protein_id = sequence_ids[-1]
    last_protein_length = len(seqs[last_protein_id])
    last_global_representation_shape = global_representations[-1].shape
    print(f"Embedded protein {last_protein_id} with length {last_protein_length} to emb. of global shape: {last_global_representation_shape}")
    
    # Create a dictionary to hold the embeddings
    embeddings_dict = {seq_id: emb.flatten() for seq_id, emb in zip(sequence_ids, global_representations)}
    
    # Convert embeddings to DataFrame
    df = pd.DataFrame.from_dict(embeddings_dict, orient='index')
    
    # Save embeddings to CSV
    df.to_csv(csv_path, index_label='sequence_id')
    print(f"Embeddings saved to {csv_path}")

if __name__ == '__main__':
    main()