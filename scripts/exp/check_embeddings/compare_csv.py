import os
import pandas as pd
import numpy as np

def compare_csv_files(dir1, dir2):
    # Get all CSV files in the first directory
    csv_files = [f for f in os.listdir(dir1) if f.endswith('.csv')]

    for csv_file in csv_files:
        file1_path = os.path.join(dir1, csv_file)
        file2_path = os.path.join(dir2, csv_file)

        # Check if the corresponding file exists in the second directory
        if not os.path.exists(file2_path):
            print(f"File {csv_file} not found in {dir2}")
            continue

        # Read the CSV files
        df1 = pd.read_csv(file1_path, index_col=0)
        df2 = pd.read_csv(file2_path, index_col=0)

        # Sort both dataframes by index
        df1_sorted = df1.sort_index()
        df2_sorted = df2.sort_index()

        # Compare the sorted dataframes
        if df1_sorted.equals(df2_sorted):
            print(f"{csv_file}: Identical after sorting")
        else:
            print(f"{csv_file}: Different after sorting")
            
            # Check if indices are the same
            if set(df1_sorted.index) != set(df2_sorted.index):
                print("  Indices are different")
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
        if np.allclose(df1_sorted, df2_sorted):
            print("  Dataframes are numerically close (may have small floating point differences)")
        else:
            print("  Dataframes have significant numerical differences")

# Usage
dir1 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm'
dir2 = '/orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/csvs'

compare_csv_files(dir1, dir2)