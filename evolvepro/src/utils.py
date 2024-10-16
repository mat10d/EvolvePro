from sklearn.decomposition import PCA
import pandas as pd

# Perform PCA on the embeddings
def pca_embeddings(embeddings_df, n_components=10):
    # Store the row indices
    indices = embeddings_df.index
    
    # Perform PCA on the embeddings
    pca = PCA(n_components=n_components)
    embeddings_pca = pca.fit_transform(embeddings_df)
    # Get the embeddings for the top n_components
    embeddings_pca = embeddings_pca[:, :n_components]
    
    # Convert embeddings to a dataframe and reattach the row indices
    embeddings_pca_df = pd.DataFrame(embeddings_pca, index=indices, columns=[f'PCA {i}' for i in range(1, n_components + 1)])
    
    return embeddings_pca_df

