import os
import pandas as pd
from sklearn.decomposition import PCA

# Function to perform PCA on embeddings DataFrame
def pca_embeddings(embeddings_df, n_components=10):
    # Store the row indices
    indices = embeddings_df.index
    
    # Perform PCA on the embeddings
    pca = PCA(n_components=n_components)
    embeddings_pca = pca.fit_transform(embeddings_df)
    
    # Get the explained variance ratio for each principal component
    explained_variance_ratio = pca.explained_variance_ratio_
    
    # Convert embeddings to a DataFrame and reattach the row indices
    embeddings_pca_df = pd.DataFrame(embeddings_pca, index=indices, columns=[f'PCA {i}' for i in range(1, n_components + 1)])
    
    return embeddings_pca_df, explained_variance_ratio

# Directory containing the CSV files
csv_directory = '/orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/csvs/'

# Directory to save the PCA results
pca_output_directory = '/orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/extract/esm/results_means/pca/'

# Number of components for PCA
n_components = 10

# Iterate over each CSV file in the directory
for csv_file in os.listdir(csv_directory):
    csv_path = os.path.join(csv_directory, csv_file)

    # Print file running for
    print(f"Processing file: {csv_file}")
    
    # Read the CSV file into a DataFrame
    embeddings_df = pd.read_csv(csv_path, index_col=0)
    
    # Perform PCA on the embeddings DataFrame
    embeddings_pca_df, explained_variance_ratio = pca_embeddings(embeddings_df, n_components)
    
    # Save the PCA embeddings to a new CSV file
    output_filename = f"pca_{csv_file}"
    output_path = os.path.join(pca_output_directory, output_filename)
    embeddings_pca_df.to_csv(output_path)
    
    # Save explained variance ratio to a CSV file for each principal component
    explained_variance_df = pd.DataFrame({"Explained Variance Ratio": explained_variance_ratio})
    explained_variance_filename = f"{output_filename}_pct_explained.csv"
    explained_variance_path = os.path.join(pca_output_directory, explained_variance_filename)
    explained_variance_df.to_csv(explained_variance_path, index=False)
