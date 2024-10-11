import os
from typing import List, Dict, Any, Tuple, Union
import pandas as pd
import numpy as np
import torch

def load_dms_data(dataset_name: str, model_name: str, embeddings_path: str, labels_path: str, embeddings_file_type: str, embeddings_type_pt: str = 'both') -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Construct the file paths
    labels_file = os.path.join(labels_path, f'{dataset_name}_labels.csv')
    
    embeddings_file_name = f'{dataset_name}_{model_name}.{embeddings_file_type}'
    embeddings_file = os.path.join(embeddings_path, embeddings_file_name)

    # Read labels
    labels = pd.read_csv(labels_file)

    # Read embeddings
    if embeddings_file_type == "csv":
        embeddings = pd.read_csv(embeddings_file, index_col=0)
    elif embeddings_file_type == "pt":
        embeddings = torch.load(embeddings_file)
        if embeddings_type_pt == 'average':
            embeddings = {key: value['average'].numpy() for key, value in embeddings.items()}
        elif embeddings_type_pt == 'mutated':
            embeddings = {key: value['mutated'].numpy() for key, value in embeddings.items()}
        elif embeddings_type_pt == 'both':
            embeddings = {key: np.concatenate((value['average'].numpy(), value['mutated'].numpy())) for key, value in embeddings.items()}
        else:
            print("Invalid embeddings_type_pt. Please choose 'average', 'mutated', or 'both'")
            return None, None
        embeddings = pd.DataFrame.from_dict(embeddings, orient='index')
    else:
        print("Invalid file type. Please choose either 'csv' or 'pt'")
        return None, None

    # Filter labels and embeddings
    labels = labels[labels['fitness'].notna()]
    embeddings = embeddings[embeddings.index.isin(labels['variant'])]

    # Sort the labels and embeddings by variant
    labels = labels.sort_values(by=['variant']).reset_index(drop=True)
    embeddings = embeddings.loc[labels['variant']]

    # Check if the labels and embeddings are aligned
    if labels['variant'].tolist() == embeddings.index.tolist():
        print('Embeddings and labels are aligned')
        return embeddings, labels
    else:
        print('Embeddings and labels are not aligned')
        return None, None

def read_experimental_data(base_path: str, round_file_name: str, WT_sequence: str, single_mutant: bool = True) -> pd.DataFrame:
    file_path = base_path + '/rounds/' + round_file_name
    df = pd.read_excel(file_path)

    # Iterate through the 'Variant' column and update the values based on t7_sequence
    if single_mutant:
        updated_variants = []
        for _, row in df.iterrows():
            variant = row['Variant']
            if variant == 'WT':
                updated_variants.append(variant)
            else:
                position = int(variant[:-1])
                wt_aa = WT_sequence[position - 1]
                updated_variant = wt_aa + variant
                updated_variants.append(updated_variant)
        
        df['updated_variant'] = updated_variants  # Add the updated variants to the DataFrame
    else:
        df.rename(columns={'Variant': 'updated_variant'}, inplace=True)

    return df

def create_dataframes(df_list: List[pd.DataFrame], expected_index: List[str]) -> Tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
    # First dataframe
    dfs = []  # List to store modified dataframes
    
    for i, df in enumerate(df_list, start=1):
        # Create a copy of the dataframe
        df_copy = df_list[i - 1].copy()
        # If the variant is WT, and i is equal to 1 assign iteration number 0
        if i == 1:
            df_copy.loc[df_copy['updated_variant'] == 'WT', 'iteration'] = 0
        else:
            df_copy = df_copy[df_copy['updated_variant'] != 'WT']
        df_copy.loc[df_copy['updated_variant'] != 'WT', 'iteration'] = i
        df_copy['iteration'] = df_copy['iteration'].astype(int)
        df_copy.rename(columns={'updated_variant': 'variant'}, inplace=True)  # Rename the column
    
        dfs.append(df_copy)

    df1 = pd.concat(dfs, ignore_index=True)
    df2 = pd.concat(dfs, ignore_index=True)

    # Check for duplicates in the 'variant' column of df1 or df2
    df1_duplicates = df1[df1.duplicated(subset=['variant'], keep=False)]
    df2_duplicates = df2[df2.duplicated(subset=['variant'], keep=False)]

    if not df1_duplicates.empty or not df2_duplicates.empty:
        print("Duplicates found in variant column:")
        if not df1_duplicates.empty:
            print("Duplicates in df1:")
            print(df1_duplicates)
        if not df2_duplicates.empty:
            print("Duplicates in df2:")
            print(df2_duplicates)
        print("Exiting.")
        return None, None

    df1 = df1[['variant', 'iteration']]
    df2 = df2[['variant', 'fitness', 'iteration']]

    expected_index_blank = [variant for variant in expected_index if variant not in df2['variant'].tolist()]
    # make a df_external that has a column 'variant' with all the variants in expected_index
    df_external = pd.DataFrame({'variant': expected_index_blank})
    df_external['fitness'] = np.nan  
    df_external['iteration'] = 1001 
    df2 = df2.append(df_external, ignore_index=True)
    # order df2 by expected_index
    df2 = df2.set_index('variant').reindex(expected_index, fill_value=np.nan).reset_index()
    # rename column 'index' to 'variant'
    df2 = df2.rename(columns={'index': 'variant'})
    
    return df1, df2