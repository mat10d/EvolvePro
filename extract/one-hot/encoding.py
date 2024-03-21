from Bio import SeqIO
import numpy as np
import pandas as pd

encoding_studies = ["brenan", "jones", "stiffler", "haddox", "doud", "giacomelli", "kelsic", "lee", "markin", "cas12f", "cov2_S", "zikv_E"]

fasta_path="directed_evolution/data_processing/output/"
results_path="directed_evolution/extract/one-hot/results/"

def integer_encode(sequence):
    # Define mapping of amino acids to integers
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_to_int = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Integer encode the sequence
    integer_encoded_sequence = [aa_to_int[aa] for aa in sequence if aa in amino_acids]
    
    return integer_encoded_sequence

def one_hot_encode(sequence):
    # Define mapping of amino acids to indices
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Initialize one-hot encoded array
    one_hot = np.zeros((len(sequence), len(amino_acids)))
    
    # Fill in the one-hot encoded array
    for i, aa in enumerate(sequence):
        if aa in amino_acids:
            one_hot[i, aa_to_index[aa]] = 1
    
    return one_hot

def encode_sequences_from_fasta(fasta_file, method):
    dataset_name = fasta_file.split('/')[-1].split('.')[0]
    print(dataset_name)
    encoded_sequences = []
    row_names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        row_names.append(record.id)  # Store the title of each sequence
        if method == 'one_hot':
            encoded_sequence = one_hot_encode(sequence)
            # convert this from 2D array to 1D vector
            encoded_sequence = encoded_sequence.flatten()
        elif method == 'integer':
            encoded_sequence = integer_encode(sequence)
        encoded_sequences.append(encoded_sequence)
    # convert the list of arrays to a 2D array
    encoded_sequences = np.array(encoded_sequences)
    # set the row names of encoded_sequences as row_names
    encoded_sequences = pd.DataFrame(encoded_sequences, index=row_names)
    # print shape of encoded_sequences
    print(encoded_sequences.shape)
    # save the encoded_sequences to a csv file
    encoded_sequences.to_csv(f'{results_path}/csvs/{dataset_name}_{method}_encoded.csv')
    return encoded_sequences

# Loop over encoding methods and studies
for method in ['one_hot', 'integer']:
    for study in encoding_studies:
        fasta_file = f"{fasta_path}{study}.fasta"
        encode_sequences_from_fasta(fasta_file, method)

