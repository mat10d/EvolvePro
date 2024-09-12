import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
import pandas as pd

def generate_single_aa_mutants_fasta(wt_sequence, output_file):
    """
    Generate a FASTA file containing all possible single amino acid mutations.

    Args:
        wt_sequence (str): Wild-type protein sequence.
        output_file (str): Path to the output FASTA file.
    """
    aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    records = [SeqRecord(Seq(wt_sequence), id="WT", description="Wild-type sequence")]
    
    for i, wt_aa in enumerate(wt_sequence):
        for mutant_aa in aa_alphabet:
            if mutant_aa != wt_aa:
                mutant_sequence = wt_sequence[:i] + mutant_aa + wt_sequence[i+1:]
                variant = f'{wt_aa}{i+1}{mutant_aa}'
                record = SeqRecord(Seq(mutant_sequence), id=variant, description="")
                records.append(record)

    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")
    
    print(f"Number of records: {len(records)}")

def generate_n_mutant_combinations_fasta(wt_sequence, output_file, mutants, n):
    """
    Generate a FASTA file containing combinations of n mutations.

    Args:
        wt_sequence (str): Wild-type protein sequence.
        output_file (str): Path to the output FASTA file.
        mutants (pd.DataFrame): DataFrame containing variant information.
        n (int): Number of mutations to combine.
    """
    records = [SeqRecord(Seq(wt_sequence), id="WT", description="Wild-type sequence")]
    mutant_combinations = list(combinations(mutants['variant'], n))

    for combination in mutant_combinations:
        positions = set()
        valid_combination = True
        mutant_sequence = wt_sequence
        variant = ""

        for mutant in combination:
            wt_aa, position, mutant_aa = mutant[0], mutant[1:-1], mutant[-1]
            i = int(position) - 1

            if i in positions:
                print(f"Invalid combination: {combination}")
                valid_combination = False
                break

            positions.add(i)
            mutant_sequence = mutant_sequence[:i] + mutant_aa + mutant_sequence[i + 1:]
            variant += f'{wt_aa}{position}{mutant_aa}_'

        if valid_combination:
            print(f"Combination: {combination}")
            record = SeqRecord(Seq(mutant_sequence), id=variant.rstrip('_'), description="")
            records.append(record)

    print(f"Number of combinations: {len(mutant_combinations)}")
    print(f"Number of records: {len(records)}")
    
    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")

def process_mutations(wt_sequence, output_dir, mutants=None, n=None):
    """
    Process mutations and generate FASTA files.

    Args:
        wt_sequence (str): Wild-type protein sequence.
        output_dir (str): Directory to save output files.
        mutants (pd.DataFrame, optional): DataFrame containing variant information for n-mutant combinations.
        n (int, optional): Number of mutations to combine if using n-mutant combinations.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate single amino acid mutants
    single_mutant_file = os.path.join(output_dir, "single_mutants.fasta")
    generate_single_aa_mutants_fasta(wt_sequence, single_mutant_file)
    
    # Generate n-mutant combinations if mutants and n are provided
    if mutants is not None and n is not None:
        n_mutant_file = os.path.join(output_dir, f"{n}_mutant_combinations.fasta")
        generate_n_mutant_combinations_fasta(wt_sequence, n_mutant_file, mutants, n)

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from io import StringIO
# import os
# import xlrd
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from itertools import combinations
# import warnings

# def generate_single_aa_mutants_fasta(wt_sequence, output_file):
#     aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"

#     records = []
    
#     # Add the wild-type sequence as the first record
#     wt_record = SeqRecord(Seq(wt_sequence), id="WT", description="Wild-type sequence")
#     records.append(wt_record)
    
#     for i, wt_aa in enumerate(wt_sequence):
#         for mutant_aa in aa_alphabet:
#             if mutant_aa != wt_aa:
#                 mutant_sequence = wt_sequence[:i] + mutant_aa + wt_sequence[i+1:]
#                 variant = f'{wt_aa}{i+1}{mutant_aa}'
#                 record = SeqRecord(Seq(mutant_sequence), id=variant, description="")
#                 records.append(record)

#     with open(output_file, "w") as handle:
#         SeqIO.write(records, handle, "fasta")
    
#     # Print the number of records
#     num_records = len(records)
#     print(f"Number of records: {num_records}")

# def generate_n_mutant_combinations_fasta(wt_sequence, output_file, mutants, n):
#     records = []

#     # Add the wild-type sequence as the first record
#     wt_record = SeqRecord(Seq(wt_sequence), id="WT", description="Wild-type sequence")
#     records.append(wt_record)

#     mutant_combinations = list(combinations(mutants['variant'], n))

#     for combination in mutant_combinations:
#         # rest of the code remains the same
#         positions = set()
#         valid_combination = True
#         mutant_sequence = wt_sequence
#         variant = ""

#         for mutant in combination:
#             wt_aa = mutant[0]
#             position = mutant[1:-1]  # Extract position from the middle of the string
#             mutant_aa = mutant[-1]            
            
#             i = int(position) - 1  # Convert position to 0-based index
#             if i in positions:
#                 # Position is already used in this combination
#                 print(f"Invalid combination: {combination}")
#                 valid_combination = False
#                 break

#             positions.add(i)
#             mutant_sequence = mutant_sequence[:i] + mutant_aa + mutant_sequence[i + 1:]
#             variant += f'{wt_aa}{position}{mutant_aa}_'

#         if valid_combination:
#             # print the combination
#             print(f"Combination: {combination}")
#             record = SeqRecord(Seq(mutant_sequence), id=variant.rstrip('_'), description="")
#             records.append(record)

#     # print number of records
#     print(f"Number of combinations: {len(mutant_combinations)}")
#     print(f"Number of records: {len(records)}")
#     with open(output_file, "w") as handle:
#         SeqIO.write(records, handle, "fasta")

