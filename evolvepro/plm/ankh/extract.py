import ankh
import torch
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO

# adapted from https://huggingface.co/ElnaggarLab/ankh2-large

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

def create_arg_parser():
    """Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(description="Arguments to extract protein_bert global embeddings from fasta file.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to input FASTA file.")
    parser.add_argument("-o", "--output", required=True, type=str, help="Path to output CSV file.")
    parser.add_argument("--model", choices=['large', 'base'], default='large', help="Choose the model type (large/base)")
    parser.add_argument("--batch_size", type=int, default=50, help="Batch size for processing sequences.")
    return parser

def read_fasta(fasta_path):
    sequences = {}
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_id = record.id
            sequence = str(record.seq).upper()
            sequences[sequence_id] = sequence
    return sequences

def batch_iterable(iterable, batch_size):
    """Yields batches of a given size from an iterable."""
    for i in range(0, len(iterable), batch_size):
        yield iterable[i:i + batch_size]

def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    
    seq_path = Path(args.input)
    csv_path = Path(args.output)
    
    # Read FASTA file
    seqs = read_fasta(seq_path)

    # Prepare input for the model
    sequences = list(seqs.values())
    sequence_ids = list(seqs.keys())

    # Load the selected model
    if args.model == 'large':
        model, tokenizer = ankh.load_large_model()
    elif args.model == 'base':
        model, tokenizer = ankh.load_base_model()
    else:
        raise ValueError("Invalid model type. Choose 'large' or 'base'.")

    model.to(device)
    model.eval()

    if device.type == 'cuda':
        print(f"Using GPU: {torch.cuda.get_device_name(0)}")

    embeddings_dict = {}

    for batch_seqs, batch_ids in zip(batch_iterable(sequences, args.batch_size), batch_iterable(sequence_ids, args.batch_size)):
        ids = tokenizer.batch_encode_plus(batch_seqs, add_special_tokens=True, padding="longest")
        input_ids = torch.tensor(ids['input_ids']).to(device)
        attention_mask = torch.tensor(ids['attention_mask']).to(device)
        
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask)

        last_hidden_states = embedding.last_hidden_state
        batch_embeddings = torch.mean(last_hidden_states, dim=1).cpu().numpy()

        for sequence_id, emb in zip(batch_ids, batch_embeddings):
            embeddings_dict[sequence_id] = emb

    df = pd.DataFrame.from_dict(embeddings_dict, orient='index')
    df.to_csv(csv_path)
    print(f"Embeddings saved to {csv_path}")

if __name__ == '__main__':
    main()