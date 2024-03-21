import os
import torch
import pandas as pd

# Specify the parent directory containing the study subfolders
parent_directory = '/om/group/abugoot/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/'

# Specify the output directory
output_directory = '/om/group/abugoot/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/csvs/'

# Specify subfolders we want to concatenate
concatenate_folders = ["brenan", "jones", "stiffler", "haddox", "doud", "giacomelli", "kelsic", "lee", "markin", "cas12f", "cov2_S", "zikv_E"]

# Iterate over each study subfolder in the parent directory
for study_folder_name in os.listdir(parent_directory):
    if study_folder_name in concatenate_folders:
        study_folder_path = os.path.join(parent_directory, study_folder_name)

        # Skip non-directory items
        print(study_folder_path)
        if not os.path.isdir(study_folder_path):
            continue

        # Iterate over the model subfolders within the study subfolder
        for model_folder_name in os.listdir(study_folder_path):
            model_folder_path = os.path.join(study_folder_path, model_folder_name)

            # Skip non-directory items
            if not os.path.isdir(model_folder_path):
                continue

            # Create a list to store DataFrames
            dataframes = []

            # List all .pt files in the model subfolder
            files = []
            for r, d, f in os.walk(model_folder_path):
                for file in f:
                    if file.endswith('.pt'):
                        files.append(os.path.join(r, file))

            # Iterate over each file in the model subfolder
            for file_path in files:
                file_data = torch.load(file_path)  # Load the file data
                label = file_data['label']  # Extract the label
                representations = file_data['mean_representations']  # Extract the mean representations

                # Extract the single key-value pair from representations
                key, tensor = representations.popitem()

                row_name = label  # Unique row name with label and model name
                row_data = tensor.tolist()  # Convert the tensor to a list
                new_df = pd.DataFrame([row_data], index=[row_name])

                # Append the new DataFrame to the list
                dataframes.append(new_df)

            # Concatenate all DataFrames in the list
            if dataframes:
                concatenated_df = pd.concat(dataframes)
                print("Shape of concatenated DataFrame:", concatenated_df.shape)
                # Save the concatenated DataFrame as a CSV file with the study and model names
                output_filename = f"{study_folder_name}_{model_folder_name}.csv"
                output_path = os.path.join(output_directory, output_filename)
                concatenated_df.to_csv(output_path)
