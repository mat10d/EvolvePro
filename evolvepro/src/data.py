import os
from typing import List, Dict, Any, Tuple, Union
from Bio import SeqIO
import pandas as pd
import numpy as np
import torch

def load_dms_data(dataset_name: str, model_name: str, embeddings_path: str, labels_path: str, 
                  embeddings_file_type: str, embeddings_type_pt: str = 'both') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load DMS data from files and align embeddings with labels.

    Args:
        dataset_name (str): Name of the dataset.
        model_name (str): Name of the model used for embeddings.
        embeddings_path (str): Path to the embeddings file.
        labels_path (str): Path to the labels file.
        embeddings_file_type (str): File type of embeddings ('csv' or 'pt').
        embeddings_type_pt (str, optional): Type of embeddings to use if 'pt' file ('average', 'mutated', or 'both'). Defaults to 'both'.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Aligned embeddings and labels DataFrames.
    """
    # Generate file paths
    labels_file = os.path.join(labels_path, f'{dataset_name}_labels.csv')
    embeddings_file = os.path.join(embeddings_path, f'{dataset_name}_{model_name}.{embeddings_file_type}')

    # Load labels
    labels = pd.read_csv(labels_file)
    
    # Process embeddings based on file type
    if embeddings_file_type == "csv":
        embeddings = pd.read_csv(embeddings_file, index_col=0)
    elif embeddings_file_type == "pt":
        embeddings = torch.load(embeddings_file)
        embeddings = process_pt_embeddings(embeddings, embeddings_type_pt)
        if embeddings is None:
            return None, None
        embeddings = pd.DataFrame.from_dict(embeddings, orient='index')
    else:
        print("Invalid file type. Please choose either 'csv' or 'pt'")
        return None, None

    # Align embeddings with labels
    labels = labels[labels['activity'].notna()]
    embeddings = embeddings[embeddings.index.isin(labels['variant'])]

    # Sort labels and embeddings by variant
    labels = labels.sort_values(by=['variant']).reset_index(drop=True)
    embeddings = embeddings.loc[labels['variant']]

    # Check if embeddings and labels are aligned
    if labels['variant'].tolist() == embeddings.index.tolist():
        print('Embeddings and labels are aligned')
        return embeddings, labels
    else:
        print('Embeddings and labels are not aligned')
        return None, None

def process_pt_embeddings(embeddings: Dict, embeddings_type_pt: str) -> Dict:
    """
    Process embeddings from .pt file based on the specified type.

    Args:
        embeddings (Dict): Dictionary containing embeddings.
        embeddings_type_pt (str): Type of embeddings to use ('average', 'mutated', or 'both').

    Returns:
        Dict: Processed embeddings dictionary.
    """
    # Get average or mutated embeddings
    if embeddings_type_pt == 'average':
        return {key: value['average'].numpy() for key, value in embeddings.items()}
    elif embeddings_type_pt == 'mutated':
        return {key: value['mutated'].numpy() for key, value in embeddings.items()}
    # Concatenate average and mutated embeddings
    elif embeddings_type_pt == 'both':
        return {key: np.concatenate((value['average'].numpy(), value['mutated'].numpy())) for key, value in embeddings.items()}
    else:
        print("Invalid embeddings_type_pt. Please choose 'average', 'mutated', or 'both'")
        return None

def load_experimental_embeddings(base_path: str, embeddings_file_name: str, rename_WT: bool = False) -> pd.DataFrame:
    """
    Load experimental embeddings from file.

    Args:
        base_path (str): Base path to the data directory.
        embeddings_file_name (str): Name of the embeddings file.

    Returns:
        pd.DataFrame: Experimental embeddings.
    """
    file_path = os.path.join(base_path, embeddings_file_name)
    embeddings = pd.read_csv(file_path, index_col=0)

    # Rename 'WT Wild-type sequence' to 'WT'
    if rename_WT:
        embeddings = embeddings.rename(index={'WT Wild-type sequence': 'WT'})

    return embeddings

def load_experimental_data(base_path: str, round_file_name: str, wt_fasta_path: str, single_mutant: bool = True) -> pd.DataFrame:
    """
    Load experimental data from file and process variants.

    Args:
        base_path (str): Base path to the data directory.
        round_file_name (str): Name of the round file.
        wt_fasta_path (str): Path to the wild-type FASTA file.
        single_mutant (bool, optional): Flag for single mutant processing. Defaults to True.

    Returns:
        pd.DataFrame: Processed experimental data.
    """
    # Load experimental data
    file_path = os.path.join(base_path, round_file_name)
    df = pd.read_excel(file_path)

    # Load wild-type sequence
    WT_sequence = str(SeqIO.read(wt_fasta_path, "fasta").seq)

    # Process variants
    if single_mutant:
        df['updated_variant'] = df['Variant'].apply(lambda x: process_variant(x, WT_sequence))
    else:
        df.rename(columns={'Variant': 'updated_variant'}, inplace=True)

    return df

def process_variant(variant: str, WT_sequence: str) -> str:
    """
    Process a single variant.

    Args:
        variant (str): Variant string.
        WT_sequence (str): Wild-type sequence.

    Returns:
        str: Processed variant string.
    """
    # Check if variant is WT
    if variant == 'WT':
        return variant
    
    # Extract position and amino acids
    position = int(variant[:-1])
    wt_aa = WT_sequence[position - 1]
    return wt_aa + variant

def create_iteration_dataframes(df_list: List[pd.DataFrame], expected_variants: List[str]) -> Tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
    """
    Create training and testing dataframes for iterative learning.

    Args:
        df_list (List[pd.DataFrame]): List of DataFrames containing experimental data from each round.
        expected_variants (List[str]): List of all expected variant names.

    Returns:
        Tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]: 
            - iteration DataFrame: Contains variant and iteration information for training.
            - labels DataFrame: Contains variant, activity, and iteration information for testing.
            Returns (None, None) if duplicates are found.
    """
    processed_dfs = []

    # Process each round's data
    for round_num, df in enumerate(df_list, start=1):
        df_copy = df.copy()
        
        # Set iteration for WT in first round, exclude WT from subsequent rounds
        if round_num == 1:
            df_copy.loc[df_copy['updated_variant'] == 'WT', 'iteration'] = 0
        else:
            df_copy = df_copy[df_copy['updated_variant'] != 'WT']
        
        df_copy.loc[df_copy['updated_variant'] != 'WT', 'iteration'] = round_num
        df_copy['iteration'] = df_copy['iteration'].astype(float)
        df_copy.rename(columns={'updated_variant': 'variant'}, inplace=True)
        
        processed_dfs.append(df_copy)

    # Combine all processed dataframes
    combined_df = pd.concat(processed_dfs, ignore_index=True)

    # Check for duplicates
    if has_duplicates(combined_df):
        return None, None

    # Create iter_train dataframe
    iteration = combined_df[['variant', 'iteration']]

    # Create iter_test dataframe
    labels = combined_df[['variant', 'activity', 'iteration']]

    # Add a activity_binary and activity_scaled column to labels
    labels['activity_binary'] = labels['activity'].apply(lambda x: 1 if x >= 1 else 0)
    labels['activity_scaled'] = labels['activity'].apply(lambda x: (x - labels['activity'].min()) / (labels['activity'].max() - labels['activity'].min()))

    # Add missing variants to iter_test
    labels = add_missing_variants(labels, expected_variants)

    # Reorder iter_test based on expected variants
    labels = labels.set_index('variant').reindex(expected_variants, fill_value=np.nan).reset_index()
    labels.rename(columns={'index': 'variant'}, inplace=True)
    
    return iteration, labels

def has_duplicates(df: pd.DataFrame) -> bool:
    """
    Check for duplicates in the 'variant' column of the dataframe.

    Args:
        df (pd.DataFrame): DataFrame to check for duplicates.

    Returns:
        bool: True if duplicates are found, False otherwise.
    """
    # Find duplicates in the 'variant' column
    duplicates = df[df.duplicated(subset=['variant'], keep=False)]
    
    # Print duplicates if found
    if not duplicates.empty:
        print("Duplicates found in variant column:")
        print(duplicates)
        print("Exiting.")
        return True
    return False

def add_missing_variants(df: pd.DataFrame, expected_variants: List[str]) -> pd.DataFrame:
    """
    Add missing variants to the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame to add missing variants to.
        expected_variants (List[str]): List of all expected variant names.

    Returns:
        pd.DataFrame: DataFrame with missing variants added.
    """
    missing_variants = set(expected_variants) - set(df['variant'])
    missing_df = pd.DataFrame({
        'variant': list(missing_variants),
        'activity': np.nan,
        'activity_binary': np.nan,
        'activity_scaled': np.nan,
        'iteration': np.nan
    })
    return pd.concat([df, missing_df], ignore_index=True)