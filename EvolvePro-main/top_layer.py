import pandas as pd
import numpy as np
from sklearn import linear_model, preprocessing
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy.spatial.distance import cdist
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
import xgboost
from sklearn.decomposition import PCA
import warnings
import random
import time
import os
import sys
import argparse
import torch
from sklearn.cluster import KMeans
from Bio import SeqIO

# Ignore FutureWarnings and ConvergenceWarnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn.neural_network")
pd.options.mode.chained_assignment = None  # default='warn'

# Function to read in the data
def read_data(dataset_name, base_path, file_type, embeddings_type='both', experimental = False):
    # Construct the file paths
    """
    Read in data from files based on the specified parameters.
    
    Args:
    dataset_name (str): Name of the dataset
    base_path (str): Base path for the data files
    file_type (str): Type of files to read ('csvs' or 'pts')
    embeddings_type (str): Type of embeddings to use ('average', 'mutated', or 'both')
    experimental (bool): Whether this is experimental data
    
    Returns:
    tuple: embeddings, labels, and hierarchy data (if not experimental)
    """

    if file_type == "csvs":
        labels_file = os.path.join(base_path, 'labels', dataset_name.split('_')[0] + '_labels.csv')
        print(labels_file)
        hie_file = os.path.join(base_path, 'hie_temp', dataset_name.split('_')[0] + '.csv')
        embeddings_file = os.path.join(base_path, 'csvs', dataset_name + '.csv')
        print(embeddings_file)
        # Read in mean embeddings across all rounds
        embeddings = pd.read_csv(embeddings_file, index_col=0)
    elif file_type == "pts":
        labels_file = os.path.join(base_path, 'labels', dataset_name.split('_')[-1] + '_labels.csv')
        print(labels_file)
        hie_file = os.path.join(base_path, 'hie_temp', dataset_name.split('_')[-1] + '.csv')
        embeddings_file = os.path.join(base_path, 'pts', dataset_name + '.pt')
        print(embeddings_file)
        # Read in pytorch tensor of embeddings
        embeddings = torch.load(embeddings_file)
        # Convert embeddings to a dataframe
        if embeddings_type == 'average':
            embeddings = {key: value['average'].numpy() for key, value in embeddings.items()}
        elif embeddings_type == 'mutated':
            embeddings = {key: value['mutated'].numpy() for key, value in embeddings.items()}
        elif embeddings_type == 'both':
            embeddings = {key: torch.cat((value['average'], value['mutated'])).numpy() for key, value in embeddings.items()}
        else:
            print("Invalid embeddings_type. Please choose 'average', 'mutated', or 'both'")
            return None, None

        # Convert embeddings dictionary to a dataframe
        embeddings = pd.DataFrame.from_dict(embeddings, orient='index')
    else:
        print("Invalid file type. Please choose either 'csvs' or 'pts'")
        return None, None

    # if not experimental
    if not experimental:
        # Read in labels
        labels = pd.read_csv(labels_file)

        # Read in hierarchy
        hie_data = pd.read_csv(hie_file)

        # Filter out rows where fitness is NaN
        labels = labels[labels['fitness'].notna()]

        # Filter out rows in embeddings where row names are not in labels variant column
        embeddings = embeddings[embeddings.index.isin(labels['variant'])]

        # Align labels by variant
        labels = labels.sort_values(by=['variant'])

        # Align embeddings by row name
        embeddings = embeddings.sort_index()

        # Confirm that labels and embeddings are aligned, reset index
        labels = labels.reset_index(drop=True)

        # Get the variants in labels and embeddings, convert to list
        label_variants = labels['variant'].tolist()
        embedding_variants = embeddings.index.tolist()

        # Check if embedding row names and label variants are identical
        if label_variants == embedding_variants:
            print('Embeddings and labels are aligned')

        # return embeddings and labels
        return embeddings, labels, hie_data

    else:
        return embeddings


# Active learning function for one iteration
def top_layer(iter_train, iter_test, embeddings_pd, labels_pd, measured_var, regression_type='ridge', top_n=None, final_round=10):

    """
    Perform one iteration of active learning.
    
    Args:
    iter_train (list): Iterations to use for training
    iter_test (int): Iteration to use for testing
    embeddings_pd (pd.DataFrame): Embeddings data
    labels_pd (pd.DataFrame): Labels data
    measured_var (str): Variable to predict
    regression_type (str): Type of regression model to use
    top_n (int): Number of top predictions to return
    final_round (int): Final round number
    
    Returns:
    tuple: Test predictions and all predictions
    """

    # Get the variants in labels and embeddings, convert to list
    label_variants = labels_pd['variant'].tolist()
    embedding_variants = embeddings_pd.index.tolist()

    # Check if embedding row names and label variants are identical
    if label_variants == embedding_variants:
        print('Embeddings and labels are aligned')

    # reset the indices of embeddings_pd and labels_pd
    embeddings_pd = embeddings_pd.reset_index(drop=True)
    labels_pd = labels_pd.reset_index(drop=True)

    # save column 'iteration' in the labels dataframe
    iteration = labels_pd['iteration']

    # save labels
    labels = labels_pd

    # save mean embeddings as numpy array
    a = embeddings_pd

    # subset a, y to only include the rows where iteration = iter_train and iter_test
    idx_train = iteration[iteration.isin(iter_train)].index.to_numpy()
    idx_test = iteration[iteration.isin([iter_test])].index.to_numpy()

    # subset a to only include the rows where iteration = iter_train and iter_test
    X_train = a.loc[idx_train, :]
    X_test = a.loc[idx_test, :]

    y_train = labels[iteration.isin(iter_train)][measured_var]

    y_test = labels[iteration.isin([iter_test])][measured_var]

    # fit
    if regression_type == 'ridge':
        model = linear_model.RidgeCV()
    elif regression_type == 'lasso':
        model = linear_model.LassoCV(max_iter=100000,tol=1e-3)
    elif regression_type == 'elasticnet':
        model = linear_model.ElasticNetCV(max_iter=100000,tol=1e-3)
    elif regression_type == 'linear':
        model = linear_model.LinearRegression()
    elif regression_type == 'neuralnet':
        model = MLPRegressor(hidden_layer_sizes=(5), max_iter=1000, activation='relu', solver='adam', alpha=0.001,
                             batch_size='auto', learning_rate='constant', learning_rate_init=0.001, power_t=0.5,
                             momentum=0.9, nesterovs_momentum=True, shuffle=True, random_state=1, tol=0.0001,
                             verbose=False, warm_start=False, early_stopping=False, validation_fraction=0.1, beta_1=0.9,
                             beta_2=0.999, epsilon=1e-08)
    elif regression_type == 'randomforest':
        model = RandomForestRegressor(n_estimators=100, criterion='friedman_mse', max_depth=None, min_samples_split=2,
                                      min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_features='auto',
                                      max_leaf_nodes=None, min_impurity_decrease=0.0, bootstrap=True, oob_score=False,
                                      n_jobs=None, random_state=1, verbose=0, warm_start=False, ccp_alpha=0.0,
                                      max_samples=None)
    elif regression_type == 'gradientboosting':
        model = xgboost.XGBRegressor(objective='reg:squarederror', colsample_bytree=0.3, learning_rate=0.1,
                                     max_depth=5, alpha=10, n_estimators=10)

    model.fit(X_train, y_train)

    # make predictions on train data
    y_pred_train = model.predict(X_train)
    y_std_train = np.zeros(len(y_pred_train))
    # make predictions on test data
    # NOTE: can work on alternate 2-n round strategies here
    y_pred_test = model.predict(X_test)
    y_std_test = np.zeros(len(y_pred_test))

    # combine predicted and actual thermostability values with sequence IDs into a new dataframe
    df_train = pd.DataFrame({'variant': labels.variant[idx_train], 'y_pred': y_pred_train, 'y_actual': y_train})
    df_test = pd.DataFrame({'variant': labels.variant[idx_test], 'y_pred': y_pred_test, 'y_actual': y_test})
    
    # sort df_test by y_pred
    df_test = df_test.sort_values(by=['y_pred'], ascending=False)

    df_all = pd.concat([df_train, df_test])

    # sort df_all by y_pred
    df_all = df_all.sort_values(by=['y_pred'], ascending=False)

    return df_test, df_all

def read_experimental_data(base_path, round_file_name, WT_sequence, single_mutant=True):
    """
    Read experimental data from an Excel file and process variant information.
    
    Args:
    base_path (str): Base path for the data file
    round_file_name (str): Name of the Excel file
    WT_sequence (str): Wild-type sequence
    single_mutant (bool): Whether to process single mutants
    
    Returns:
    pd.DataFrame: Processed experimental data
    """

    file_path = base_path + round_file_name
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

def create_dataframes(df_list, expected_index):
    """
    Create and process dataframes from a list of experimental data.
    
    Args:
    df_list (list): List of dataframes with experimental data
    expected_index (list): Expected index for the final dataframe
    
    Returns:
    tuple: Two processed dataframes
    """
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

# import  data
dataset_name = 'hc_esm2_t48_15B_UR50D'
base_path = './results_means/'
file_type = 'csvs'
embeddings_type = 'average'
experimental = True
embeddings = read_data(dataset_name, base_path, file_type, embeddings_type, experimental)
# replace WT Wild-type sequence index in embeddings with 'WT'
embeddings = embeddings.rename(index={'WT Wild-type sequence': 'WT'})
print(embeddings)


experiment_path = './experimental_data/'
round_file_name = 'Antibody_Round1.xlsx'
protein_seq = 'MQVQLVESGGGVVQPGRSLRLSCAASGFTFSNYAMYWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRTEDTAVYYCASGSDYGDYLLVYWGQGTLVTVSSSTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK'
experimental_data = read_experimental_data(experiment_path, round_file_name, protein_seq)
print(experimental_data)
df_list = [experimental_data]
iterations_one, labels_one = create_dataframes(df_list, embeddings.index)
print (iterations_one)
print(labels_one)




iteration_old = iterations_one
embeddings_pd = embeddings
labels_pd = labels_one
measured_var = 'fitness'
regression_type = 'randomforest'
num_mutants_per_round = 16
final_round = 16

df_test, df_all = top_layer(
    iter_train=iteration_old['iteration'].unique().tolist(),
    iter_test=1001,
    embeddings_pd=embeddings_pd,
    labels_pd=labels_pd,
    measured_var=measured_var,
    regression_type=regression_type,
    top_n=None,
    final_round=final_round
)

#df_all.to_csv('prediction/round1_all_new.csv', index=False)
df_all = df_all[df_all['y_actual'].isna()][['variant', 'y_pred']]
df_all = df_all.sort_values(by=['y_pred'], ascending=False)
df_all.to_csv('prediction/round1_prediction.csv', index=False)

