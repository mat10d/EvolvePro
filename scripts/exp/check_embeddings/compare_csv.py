import os
import pandas as pd
import numpy as np

def compare_single_csv(file_name, dir1, dir2):
    file1_path = os.path.join(dir1, file_name)
    file2_path = os.path.join(dir2, file_name)

    # Check if the files exist in both directories
    if not os.path.exists(file1_path):
        print(f"File {file_name} not found in {dir1}")
        return
    if not os.path.exists(file2_path):
        print(f"File {file_name} not found in {dir2}")
        return

    # Read the CSV files
    df1 = pd.read_csv(file1_path, index_col=0)
    df2 = pd.read_csv(file2_path, index_col=0)

    # Sort both dataframes by index
    df1_sorted = df1.sort_index()
    df2_sorted = df2.sort_index()

    # Compare the sorted dataframes
    if df1_sorted.equals(df2_sorted):
        print(f"{file_name}: Identical after sorting")
    else:
        print(f"{file_name}: Different after sorting")
        
        # Check if indices are the same
        if set(df1_sorted.index) != set(df2_sorted.index):
            print("  Indices are different")
            print(f"  Indices only in file1: {set(df1_sorted.index) - set(df2_sorted.index)}")
            print(f"  Indices only in file2: {set(df2_sorted.index) - set(df1_sorted.index)}")
        else:
            # If indices are the same, find where values differ
            diff_mask = ~(df1_sorted == df2_sorted)
            if diff_mask.any().any():
                diff_indices = diff_mask.any(axis=1)
                print(f"  Differences found at indices: {diff_indices[diff_indices].index.tolist()}")
                
                # Print a few examples of differences
                for idx in diff_indices[diff_indices].index[:5]:  # Show up to 5 examples
                    print(f"    Index {idx}:")
                    print(f"      File1: {df1_sorted.loc[idx].tolist()}")
                    print(f"      File2: {df2_sorted.loc[idx].tolist()}")
            else:
                print("  Dataframes have the same indices and values, but may differ in precision")

    # Check if the dataframes are close (for floating point differences)
    if np.allclose(df1_sorted, df2_sorted, equal_nan=True):
        print("  Dataframes are numerically close (may have small floating point differences)")
    else:
        print("  Dataframes have significant numerical differences")

# Usage
dir1 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm'
dir2 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/csvs'
csv_file_name = 't7_pol_esm2_t48_15B_UR50D.csv'  

compare_single_csv(csv_file_name, dir1, dir2)