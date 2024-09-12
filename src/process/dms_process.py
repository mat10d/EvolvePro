import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

def process_dataset(file_path, wt_fasta_path, dataset_name, fitness_column, cutoff_value, output_dir, sheet_name=None, cutoff_rule='greater_than', cutoff_percentiles=None, cutoff_function = None, AA_shift=None, drop_columns=False):
    """
    Process a dataset from an Excel or CSV file and generate a FASTA file of mutations.

    Args:
        file_path (str): Path to the input file (Excel or CSV).
        wt_fasta_path (str): Path to the WT sequence FASTA file.
        dataset_name (str): Name of the dataset to be used for output files.
        fitness_column (str): Name of the column containing fitness values.
        cutoff_value (float): Base cutoff value for binary classification.
        output_dir (str): Directory to save output files.
        sheet_name (str, optional): Name of the sheet if Excel file. Defaults to None.
        cutoff_rule (str, optional): Rule for cutoff ('greater_than' or 'less_than'). Defaults to 'greater_than'.
        cutoff_percentiles (list, optional): List of percentiles to calculate additional cutoff values. Defaults to None.
        AA_shift (int, optional): Amino acid position shift. Defaults to None.
        drop_columns (bool, optional): Drop unnecessary columns. Defaults to False.

    Returns:
        tuple: Filtered DataFrame and list of fractions above cutoff values.
    """
    # Read the input file
    if file_path.endswith('.xlsx'):
        dataframe = pd.read_excel(file_path, sheet_name=sheet_name) if isinstance(sheet_name, str) else pd.read_excel(file_path)
    elif file_path.endswith('.csv'):
        dataframe = pd.read_csv(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide an Excel (.xlsx) or CSV (.csv) file.")

    # Filter and process the data
    filtered_df = dataframe.dropna(subset=[fitness_column]).copy()
    wt_sequence = next(SeqIO.parse(wt_fasta_path, 'fasta')).seq
    
    # Generate FASTA file of mutations
    generate_mutation_fasta(filtered_df, wt_sequence, dataset_name, output_dir, AA_shift)

    # Process fitness data
    filtered_df['fitness'] = filtered_df[fitness_column]
    filtered_df['fitness_scaled'] = (filtered_df[fitness_column] - filtered_df[fitness_column].min()) / (filtered_df[fitness_column].max() - filtered_df[fitness_column].min())
    
    # Calculate cutoffs including the base cutoff value
    cutoff_percentiles = cutoff_percentiles or []
    all_cutoffs = [cutoff_value] + [np.percentile(filtered_df[fitness_column], p) for p in cutoff_percentiles]
    cutoff_labels = [''] + [f'_{p}p' for p in cutoff_percentiles]
    
    # Apply cutoffs and calculate fractions
    total_values = len(filtered_df)
    fractions_above_cutoff = []
    number_above_cutoff = []

    for i, (cutoff, label) in enumerate(zip(all_cutoffs, cutoff_labels)):
        column_name = f'fitness_binary{label}'
        if cutoff_rule == 'greater_than':
            filtered_df[column_name] = (filtered_df[fitness_column] > cutoff).astype(int)
        elif cutoff_rule == 'less_than':
            filtered_df[column_name] = (filtered_df[fitness_column] < cutoff).astype(int)
        elif cutoff_rule == 'custom':
            if cutoff_function is None:
                raise ValueError("Custom function must be provided for 'custom' cutoff rule.")
            filtered_df[column_name] = cutoff_function(filtered_df, fitness_column, cutoff).astype(int)
        else:
            raise ValueError("Invalid cutoff rule. Please specify 'greater_than', 'less_than', or 'custom'.")
        
        number = filtered_df[column_name].sum()
        fraction = (filtered_df[column_name].sum() / total_values) if total_values > 0 else 0
        number_above_cutoff.append(number)
        fractions_above_cutoff.append(fraction)


    # Print cutoff values and number and fractions above cutoff
    print(f'Cutoff values: {all_cutoffs}')
    print(f'Number above cutoff: {number_above_cutoff}')
    print(f'Fractions above cutoff: {fractions_above_cutoff}')

    # Drop unnecessary columns if specified
    if drop_columns:
        columns_to_keep = ['variant', fitness_column, 'fitness', 'fitness_scaled'] + [f'fitness_binary{label}' for label in cutoff_labels]
        filtered_df = filtered_df[columns_to_keep]

    # Save processed data
    os.makedirs(output_dir, exist_ok=True)
    filtered_df.to_csv(os.path.join(output_dir, f'{dataset_name}_labels.csv'), index=False)

    return filtered_df, fractions_above_cutoff

def generate_mutation_fasta(df, wt_sequence, dataset_name, output_dir, AA_shift):
    """
    Generate a FASTA file of mutations based on the wild-type sequence.

    Args:
        df (pd.DataFrame): DataFrame containing variant information.
        wt_sequence (str): Wild-type sequence.
        dataset_name (str): Name of the dataset.
        output_dir (str): Directory to save output files.
        AA_shift (int): Amino acid position shift.
    """
    # Create output directory and file name
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{dataset_name}.fasta')

    # Write the filtered dataframe to a FASTA file
    with open(output_file, 'w') as f:
        for variant in df['variant']:
            # Add the wild-type sequence to the FASTA file if 'WT' is present
            if 'WT' in variant:
                f.write(f'>{variant}\n{wt_sequence}\n')
            # For other variants, generate the mutated sequence
            else:
                position = int(variant[1:-1]) - 1 if AA_shift is None else int(variant[1:-1]) - AA_shift
                wt_aa, mutated_aa = variant[0], variant[-1]
                # Check if the wild-type amino acid matches the sequence
                if wt_sequence[position] == wt_aa:
                    sequence = wt_sequence[:position] + mutated_aa + wt_sequence[position+1:]
                    f.write(f'>{variant}\n{sequence}\n')
                else:
                    print(f'Error: WT amino acid at position {position} is not {wt_aa}')

def plot_mutations_per_position(df):
    """
    Plot a histogram of mutations per position.

    Args:
        df (pd.DataFrame): DataFrame containing variant information.
    """

    # Filter out rows with missing or "WT" mutations
    df_filtered = df.dropna(subset=["variant"]).query('variant != "WT"')
    mutations_per_position = {}

    # Iterate over the rows of the DataFrame and increment the count of mutations at each position
    for mutation_str in df_filtered["variant"]:
        pos = int(mutation_str[1:-1])
        mutations_per_position[pos] = mutations_per_position.get(pos, 0) + 1

    # Plot a histogram of the number of mutations per position
    plt.bar(mutations_per_position.keys(), mutations_per_position.values())
    plt.xlabel('Number of mutations')
    plt.ylabel('Number of positions')
    plt.title('Mutations per position')
    plt.show()

def plot_histogram_of_readout(df, column_name, cutoff=None):
    """
    Plot a histogram of readout values for all mutants.

    Args:
        df (pd.DataFrame): DataFrame containing readout values.
        column_name (str): Name of the column to plot.
        cutoff (float, optional): Cutoff value to display on the plot. Defaults to None.
    """
    # Plot histogram of readout values for all mutants
    fig, ax = plt.subplots()
    ax.hist(df[column_name].values, bins=100)
    ax.set_xlabel(column_name)
    ax.set_ylabel('Number of mutants')
    ax.set_title(f'{column_name} distribution across mutants')

    # Add vertical line to indicate WT value
    if "WT" in df["variant"].unique():
        wt_val = df.loc[df["variant"] == 'WT', column_name].values[0]
        ax.axvline(wt_val, color='red', linestyle='--', label='WT')
    # Add vertical line to indicate cutoff value
    if cutoff:
        ax.axvline(cutoff, color='black', linestyle='--', label='cutoff')
    ax.legend()
    plt.show()

def markin_custom_cutoff(markin_df, fitness_column, cutoff):
    """
    Apply a custom cutoff to the markin dataframe based on fitness and p-value.

    Parameters:
    - markin_df (DataFrame): The input DataFrame containing fitness and p-values.
    - fitness_column (str): The column name of the fitness values.
    - cutoff (float): The cutoff value for the fitness column.

    Returns:
    - Series: A binary Series ('fitness_binary') with values 0 or 1 based on cutoff and p-value.
    """
    
    # Apply the cutoff on the fitness column and p-value
    fitness_binary = ((markin_df[fitness_column] > cutoff) & (markin_df['kcatOverKM_cMUP_p-value'] < 0.01))

    return fitness_binary

def preprocess_cas12f(input_file, wt_fasta_output, preprocessed_output_file):
    """
    Bespoke function to preprocess the DMS data for AsCas12f.

    Args:
        input_file (str): Path to the input Excel file.
        wt_fasta_output (str): Path to save the wild-type sequence FASTA file.
        preprocessed_output_file (str): Path to save the preprocessed Excel file.

    Returns:
        str: Path to the preprocessed output file.
    """
    
    # Read the input file
    df = pd.read_excel(input_file)

    # Extract substitution and position from variant column
    df['substitution'] = df['variant'].str[0]
    df['position'] = df['variant'].str[1:].astype(int)

    # Identify wild-type rows and create WT sequence
    wt_rows = df[df['mean'] == 1].rename(columns={'substitution': 'WT'})
    wt_sequence = ''.join(wt_rows['WT'])

    # Write WT sequence to FASTA file
    with open(wt_fasta_output, 'w') as f:
        f.write(f'>AsCas12f\n{wt_sequence}\n')

    # Process non-WT rows
    df = df[df['mean'] != 1].rename(columns={'variant': 'variant_raw', 'mean': 'avg_fitness'})

    # Merge WT information
    df = df.merge(wt_rows[['position', 'WT']], on='position', how='left')

    # Create new variant column
    df['variant'] = df['WT'] + df['position'].astype(str) + df['substitution']

    # Clean up the dataframe
    df = df.drop(columns=['No'])
    df = df[(df['substitution'] != '*') & (~df['avg_fitness'].isna())]

    # Add WT row
    wt_row = pd.DataFrame({
        'variant_raw': ['WT'],
        'rep1': [1.0],
        'rep2': [1.0],
        'avg_fitness': [1.0],
        'substitution': [np.nan],
        'position': [np.nan],
        'WT': ['WT'],
        'variant': ['WT']
    })
    df = pd.concat([df, wt_row], ignore_index=True)

    # Save preprocessed data
    df.to_excel(preprocessed_output_file, index=False)

def preprocess_zikv_E(input_file, wt_fasta_output, preprocessed_output_file):
    """
    Preprocess the Zika Envelope DMS data.

    Args:
        input_file (str): Path to the input Excel file.
        wt_fasta_output (str): Path to save the wild-type sequence FASTA file.
        preprocessed_output_file (str): Path to save the preprocessed Excel file.

    Returns:
        str: Path to the preprocessed output file.
    """
    
    # Read the input file
    df = pd.read_excel(input_file, sheet_name='mutational effects')

    # Rename columns
    df = df.rename(columns={'mutation': 'variant'})

    # Extract wild-type sequence
    wt_sequence = ''.join(df[df['wildtype'] == df['mutant']]['wildtype'].values)

    # Write WT sequence to FASTA file
    with open(wt_fasta_output, 'w') as f:
        f.write(f'>Zikv_E\n{wt_sequence}\n')

    # Filter out wild-type rows
    df = df[df['wildtype'] != df['mutant']]

    # Add WT row
    wt_row = pd.DataFrame({
        'variant': ['WT'],
        'site': [np.nan],
        'wildtype': [np.nan],
        'mutant': [np.nan],
        'effect': [1.0],
        'log2effect': [0.0]
    })
    df = pd.concat([df, wt_row], ignore_index=True)

    # Save preprocessed data, with Cas12f as the sheet name
    df.to_excel(preprocessed_output_file, index=False)

def preprocess_cov2_S(input_file, preprocessed_output_file):
    """
    Preprocess the SARS-CoV-2 S protein DMS data.

    Args:
        input_file (str): Path to the input CSV file.
        preprocessed_output_file (str): Path to save the preprocessed CSV file.

    Returns:
        str: Path to the preprocessed output file.
    """
    
    # Read the input file
    df = pd.read_csv(input_file)

    # Create variant column
    df['variant'] = df['wildtype'] + df['site'].astype(str) + df['mutation']

    # Drop unnecessary columns
    df = df.drop(columns=['site_total_escape', 'site_max_escape', 'condition'])

    # Group by variant and calculate mean fitness
    df_averaged = df.groupby(['variant', 'site', 'wildtype', 'mutation'])['mut_escape'].mean().reset_index()

    # Sort by site
    df_averaged = df_averaged.sort_values(by=['site'])

    # Extract wild-type sequence
    wt_sequence = ''.join(df_averaged['wildtype'].drop_duplicates().sort_index())

    # Save preprocessed data
    df_averaged.to_csv(preprocessed_output_file, index=False)
