#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:33:22 2020
Modified to output CSV instead of HDF5

@author: mheinzinger
"""

import argparse
import time
from pathlib import Path
import pandas as pd

import torch
from transformers import T5EncoderModel, T5Tokenizer

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

def get_T5_model(model_dir, transformer_link = "Rostlab/prot_t5_xl_half_uniref50-enc"):
    print(f"Loading: {transformer_link}")
    if model_dir is not None:
        print("##########################")
        print(f"Loading cached model from: {model_dir}")
        print("##########################")
    model = T5EncoderModel.from_pretrained(transformer_link, cache_dir=model_dir)
    if device == torch.device("cpu"):
        print("Casting model to full precision for running on CPU ...")
        model.to(torch.float32)

    model = model.to(device)
    model = model.eval()
    vocab = T5Tokenizer.from_pretrained(transformer_link, do_lower_case=False)
    return model, vocab

def read_fasta(fasta_path):
    sequences = {}
    with open(fasta_path, 'r') as fasta_f:
        for line in fasta_f:
            if line.startswith('>'):
                uniprot_id = line.replace('>', '').strip()
                uniprot_id = uniprot_id.replace("/","_").replace(".","_")
                sequences[uniprot_id] = ''
            else:
                sequences[uniprot_id] += ''.join(line.split()).upper().replace("-","")
    return sequences

def get_embeddings(seq_path, emb_path, model_dir, per_protein, max_residues=4000, max_seq_len=1000, max_batch=50):
    seq_dict = read_fasta(seq_path)
    model, vocab = get_T5_model(model_dir)

    print('########################################')
    print(f'Example sequence: {next(iter(seq_dict.keys()))}\n{next(iter(seq_dict.values()))}')
    print('########################################')
    print(f'Total number of sequences: {len(seq_dict)}')

    avg_length = sum(len(seq) for seq in seq_dict.values()) / len(seq_dict)
    n_long = sum(1 for seq in seq_dict.values() if len(seq) > max_seq_len)
    seq_dict = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    
    print(f"Average sequence length: {avg_length}")
    print(f"Number of sequences >{max_seq_len}: {n_long}")
    
    start = time.time()
    batch = []
    emb_dict = {}

    for seq_idx, (pdb_id, seq) in enumerate(seq_dict, 1):
        seq = seq.replace('U','X').replace('Z','X').replace('O','X')
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id, seq, seq_len))

        n_res_batch = sum(s_len for _, _, s_len in batch) + seq_len 
        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(seq_dict) or seq_len > max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = []

            token_encoding = vocab.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)
            
            try:
                with torch.no_grad():
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                print(f"RuntimeError during embedding for {pdb_id} (L={seq_len}). Try lowering batch size. " +
                      "If single sequence processing does not work, you need more vRAM to process your protein.")
                continue
            
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                emb = embedding_repr.last_hidden_state[batch_idx, :s_len]
                
                if per_protein:
                    emb = emb.mean(dim=0)
            
                if len(emb_dict) == 0:
                    print(f"Embedded protein {identifier} with length {s_len} to emb. of shape: {emb.shape}")

                emb_dict[identifier] = emb.detach().cpu().numpy().squeeze()

    end = time.time()
    
    # Convert emb_dict to DataFrame and save as CSV
    df = pd.DataFrame.from_dict(emb_dict, orient='index')
    df.to_csv(emb_path, index_label='sequence_id')

    print('\n############# STATS #############')
    print(f'Total number of embeddings: {len(emb_dict)}')
    print(f'Total time: {end-start:.2f}[s]; time/prot: {(end-start)/len(emb_dict):.4f}[s]; avg. len= {avg_length:.2f}')
    return True

def create_arg_parser():
    parser = argparse.ArgumentParser(description='t5_embedder.py creates T5 embeddings for a given text file containing sequence(s) in FASTA-format.')
    parser.add_argument('-i', '--input', required=True, type=str, help='A path to a fasta-formatted text file containing protein sequence(s).')
    parser.add_argument('-o', '--output', required=True, type=str, help='A path for saving the created embeddings as CSV file.')
    parser.add_argument('--model', required=False, type=str, default=None, help='A path to a directory holding the checkpoint for a pre-trained model')
    parser.add_argument('--per_protein', type=int, default=0, help="Whether to return per-residue embeddings (0: default) or the mean-pooled per-protein representation (1).")
    return parser

def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    
    seq_path = Path(args.input)
    emb_path = Path(args.output)
    model_dir = Path(args.model) if args.model is not None else None

    per_protein = bool(int(args.per_protein))
    
    get_embeddings(seq_path, emb_path, model_dir, per_protein=per_protein)

if __name__ == '__main__':
    main()