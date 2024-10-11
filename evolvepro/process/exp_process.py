import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
import pandas as pd

def generate_wt(wt_sequence, output_file):
    """
    Generate a FASTA file containing the wild-type protein sequence.

    Args:
        wt_sequence (str): Wild-type protein sequence.
        output_file (str): Path to the output FASTA file.
    """
    record = SeqRecord(Seq(wt_sequence), id="WT", description="Wild-type sequence")
    with open(output_file, "w") as handle:
        SeqIO.write(record, handle, "fasta")

def generate_single_aa_mutants(wt_fasta, output_file):
    """
    Generate a FASTA file containing all possible single amino acid mutations.

    Args:
        wt_fasta (str): Path to the FASTA file containing the wild-type protein sequence.
        output_file (str): Path to the output FASTA file.
    """
    aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    wt_sequence = SeqIO.read(wt_fasta, "fasta").seq
    records = [SeqRecord(wt_sequence, id="WT", description="Wild-type sequence")]

    for i, wt_aa in enumerate(wt_sequence):
        for mutant_aa in aa_alphabet:
            if mutant_aa != wt_aa:
                mutant_sequence = wt_sequence[:i] + mutant_aa + wt_sequence[i+1:]
                variant = f'{wt_aa}{i+1}{mutant_aa}'
                record = SeqRecord(Seq(mutant_sequence), id=variant, description="")
                records.append(record)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")
    
    print(f"Number of mutants: {len(records)}")

def generate_n_mutant_combinations(wt_fasta, mutant_file, n, output_file, threshold=1):
    """
    Generate a FASTA file containing combinations of n mutations, filtered by a threshold.

    Args:
        wt_fasta (str): Path to the FASTA file containing the wild-type protein sequence.
        mutant_file (str): Path to the Excel file containing mutant information.
        n (int): Number of mutations to combine.
        output_file (str): Path to the output FASTA file.
        threshold (float): Minimum value for including a mutant (default: 1).
    """
    # Read wild-type sequence
    wt_sequence = str(SeqIO.read(wt_fasta, "fasta").seq)
    
    # Read and process mutant data
    mutants = pd.read_excel(mutant_file, header=None)
    mutants = mutants[mutants[1] > threshold]
    mutants[['position', 'mutant_aa']] = mutants[0].str.extract('(\d+)([A-Z]+)', expand=True)
    mutants['wt_aa'] = mutants.apply(lambda row: wt_sequence[int(row['position'])-1], axis=1)
    mutants['variant'] = mutants['wt_aa'] + mutants['position'] + mutants['mutant_aa']
    
    # Generate combinations
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
                valid_combination = False
                break

            positions.add(i)
            mutant_sequence = mutant_sequence[:i] + mutant_aa + mutant_sequence[i + 1:]
            variant += f'{wt_aa}{position}{mutant_aa}_'

        if valid_combination:
            record = SeqRecord(Seq(mutant_sequence), id=variant.rstrip('_'), description="")
            records.append(record)

    print(f"Number of mutant combinations: {len(mutant_combinations)}")
    print(f"Number of valid mutant combinations: {len(records)}")
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as handle:
        SeqIO.write(records, handle, "fasta")