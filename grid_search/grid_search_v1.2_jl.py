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
from sklearn_extra.cluster import KMedoids

# Ignore FutureWarnings and SettingWithCopyWarnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn.neural_network")
pd.options.mode.chained_assignment = None  # default='warn'

# Create the parser for the command line arguments
def create_parser():
    parser = argparse.ArgumentParser(description="Run experiments with different combinations of grid search variables.")
    parser.add_argument("--dataset_name", type=str, help="Name of the esm embeddings csv file of the dataset dataset")
    parser.add_argument("--base_path", type=str, help="Base path of the dataset")
    parser.add_argument("--num_simulations", type=int, help="Number of simulations for each parameter combination. Example: 3, 10")
    parser.add_argument("--num_iterations", type=int, nargs="+", help="List of number of iterations. Example: 3 5 10 (must be greater than 1)")
    parser.add_argument("--measured_var", type=str, nargs="+", help="Fitness type to train on. Options: fitness fitness_scaled")
    parser.add_argument("--learning_strategies", type=str, nargs="+", help="Type of learning strategy. Options: random top5bottom5 top10 dist")
    parser.add_argument("--num_mutants_per_round", type=int, nargs="+", help="Number of mutants per round. Example: 8 10 16 32 128")
    parser.add_argument("--first_round_strategies", type=str, nargs="+", help="Type of first round strategy. Options: random diverse_medoids representative_hie")
    parser.add_argument("--embedding_types", type=str, nargs="+", help="Types of embeddings to train on. Options: embeddings embeddings_norm embeddings_pca")
    parser.add_argument("--regression_types", type=str, nargs="+", help="Regression types. Options: ridge lasso elasticnet linear neuralnet randomforest gradientboosting")
    parser.add_argument("--file_type", type=str, help="Type of file to read. Options: csvs pts")
    parser.add_argument("--embeddings_type_pt", type=str, help="Type of pytorch embeddings to read. Options: average mutated both")
    return parser

# Function to read in the data
def read_data(dataset_name, base_path, file_type, first_round_strategies, embeddings_type='both'):
    # Construct the file paths
    if file_type == "csvs":
        labels_file = os.path.join(base_path, 'labels', dataset_name.split('_')[0] + '_labels.csv')
        hie_file = os.path.join(base_path, 'hie_temp', dataset_name.split('_')[0] + '.csv')
        embeddings_file = os.path.join(base_path, 'csvs', dataset_name + '.csv')
        # Read in mean embeddings across all rounds
        embeddings = pd.read_csv(embeddings_file, index_col=0)
    elif file_type == "pts":
        labels_file = os.path.join(base_path, 'labels', dataset_name.split('_')[-1] + '_labels.csv')
        hie_file = os.path.join(base_path, 'hie_temp', dataset_name.split('_')[-1] + '.csv')
        embeddings_file = os.path.join(base_path, 'pts', dataset_name + '.pt')
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

    # Read in labels
    labels = pd.read_csv(labels_file)

    # Read in hie
    if first_round_strategies == "representative_hie":
        hie_data = pd.read_csv(hie_file)
    else:
        hie_data = pd.DataFrame()

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

# Function to scale the embeddings in the dataframe
def scale_embeddings(embeddings_df):
    
    scaler = StandardScaler()
    scaled_embeddings = scaler.fit_transform(embeddings_df)

    scaled_embeddings_df = pd.DataFrame(scaled_embeddings)

    return scaled_embeddings_df

# Perform PCA on the embeddings
def pca_embeddings(embeddings_df, labels_df, dataset_name, n_components=8):
    # Perform PCA on the embeddings
    pca = PCA(n_components=50)
    embeddings_pca = pca.fit_transform(embeddings_df)
    # Get the embeddings for the top n_components
    embeddings_pca = embeddings_pca[:, :n_components]
    # Convert embeddings to a dataframe
    embeddings_pca_df = pd.DataFrame(embeddings_pca, columns=[f'PCA {i}' for i in range(1, n_components + 1)])

    return embeddings_pca_df

# Function for selecting mutants in the first round
def first_round(labels, embeddings, hie_data, num_mutants_per_round, first_round_strategy='random', random_seed=None):

    # Filter out 'WT' variant from labels
    variants_without_WT = labels.variant[labels.variant != 'WT']

    # Perform random first round search strategy
    if first_round_strategy == 'random':
        # Set random seed
        if random_seed is not None:
            np.random.seed(random_seed)  # Use NumPy's random seed for consistent randomization
        random_mutants = np.random.choice(variants_without_WT, size=num_mutants_per_round, replace=False)
        iteration_one_ids = random_mutants

    elif first_round_strategy == 'diverse_medoids':
        # Set random seed
        if random_seed is not None:
            np.random.seed(random_seed)  # Use NumPy's random seed for consistent randomization
        num_clusters = num_mutants_per_round
        
        # Perform PCA with 100 dimensions
        pca = PCA(n_components=100)
        pca_embeddings = pca.fit_transform(embeddings)
        pca_embeddings_reduced = pca_embeddings[:, :2]
        
        # Perform K-medoids clustering on PCA embeddings
        clusters = KMedoids(n_clusters=num_clusters, metric='euclidean', random_state=random_seed).fit(pca_embeddings_reduced)
        cluster_medoids = clusters.medoid_indices_
        cluster_labels = clusters.labels_

        # Select one medoid per cluster
        selected_mutants = []
        for cluster_idx in range(num_clusters):
            cluster_variants = variants_without_WT[cluster_labels == cluster_idx]
            if len(cluster_variants) > 0:
                # Select the medoid variant from the cluster
                selected_medoid = cluster_variants[cluster_medoids[cluster_idx]]
                selected_mutants.append(selected_medoid)

        iteration_one_ids = selected_mutants

    elif first_round_strategy == 'representative_hie':
        # iteration_one_ids should be the 0th column of the hie_data DataFrame
        iteration_one_ids = hie_data.iloc[:, 0].tolist()
    else:
        print("Invalid first round search strategy.")
        return None, None

    # Create DataFrame for the first round
    iteration_one = pd.DataFrame({'variant': iteration_one_ids, 'iteration': 1})
    WT = pd.DataFrame({'variant': 'WT', 'iteration': 0}, index=[0])
    iteration_one = iteration_one.append(WT)

    # Merge with labels DataFrame and fill null values with 1001
    labels_one = pd.merge(labels, iteration_one, on='variant', how='left')
    labels_one.iteration[labels_one.iteration.isnull()] = 1001

    return labels_one, iteration_one

# Active learning function for one iteration
def top_layer(iter_train, iter_test, embeddings_pd, labels_pd, measured_var, regression_type='ridge', top_n=None, final_round=10):
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
    y_train_fitness_scaled = labels[iteration.isin(iter_train)]['fitness_scaled']
    y_train_fitness_binary = labels[iteration.isin(iter_train)]['fitness_binary']

    y_test = labels[iteration.isin([iter_test])][measured_var]
    y_test_fitness_scaled = labels[iteration.isin([iter_test])]['fitness_scaled']
    y_test_fitness_binary = labels[iteration.isin([iter_test])]['fitness_binary']

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
                                      min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_features=1,
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

    # calculate metrics
    train_error = mean_squared_error(y_train, y_pred_train)
    test_error = mean_squared_error(y_test, y_pred_test)
    # compute train and test r^2
    train_r_squared = r2_score(y_train, y_pred_train)
    test_r_squared = r2_score(y_test, y_pred_test)
    if regression_type == 'linear' or regression_type == 'neuralnet' or regression_type == 'randomforest' or regression_type == 'gradientboosting':
        alpha = 0
    else:
        alpha = model.alpha_
    dist_metric_train = cdist(X_train, X_test, metric='euclidean').min(axis=1)
    dist_metric_test = cdist(X_test, X_train, metric='euclidean').min(axis=1)

    # combine predicted and actual thermostability values with sequence IDs into a new dataframe
    df_train = pd.DataFrame({'variant': labels.variant[idx_train], 'y_pred': y_pred_train, 'y_actual': y_train, 
                             'y_actual_scaled': y_train_fitness_scaled, 'y_actual_binary': y_train_fitness_binary,
                             'dist_metric': dist_metric_train, 'std_predictions': y_std_train})
    df_test = pd.DataFrame({'variant': labels.variant[idx_test], 'y_pred': y_pred_test, 'y_actual': y_test, 
                            'y_actual_scaled': y_test_fitness_scaled, 'y_actual_binary': y_test_fitness_binary,
                            'dist_metric': dist_metric_test, 'std_predictions': y_std_test})
    df_all = pd.concat([df_train, df_test])

    df_sorted_all = df_all.sort_values('y_pred', ascending=False).reset_index(drop=True)

    # Calculate additional metrics
    median_fitness_scaled = df_sorted_all.loc[:final_round, 'y_actual_scaled'].median()
    top_fitness_scaled = df_sorted_all.loc[:final_round, 'y_actual_scaled'].max()
    top_variant = df_sorted_all.loc[df_sorted_all['y_actual_scaled'] == top_fitness_scaled, 'variant'].values[0]

    #get spearman's rank correlation to actual data
    spearman_corr = df_sorted_all.loc[:final_round, ['y_pred', 'y_actual']].corr(method='spearman').iloc[0,1]

    fitness_binary_percentage = df_sorted_all.loc[:final_round, 'y_actual_binary'].mean()

    return train_error, test_error, train_r_squared, test_r_squared, alpha, median_fitness_scaled, top_fitness_scaled,top_variant, fitness_binary_percentage,spearman_corr, df_test

# Function to run n simulations of directed evolution
def directed_evolution_simulation(labels, embeddings, hie_data, num_simulations, num_iterations, num_mutants_per_round=10, measured_var = 'fitness', regression_type='ridge', learning_strategy='top10', top_n=None, final_round = 10, first_round_strategy = 'random'):
    output_list = []

    if first_round_strategy == 'random' or first_round_strategy == 'diverse_medoids':
        for i in range(1,num_simulations+1):
            labels_one, iteration_one = first_round(labels, embeddings, hie_data, num_mutants_per_round, first_round_strategy=first_round_strategy, random_seed=i)
            
            first_round_labels = labels_one[labels_one['iteration'] == 1]["variant"]
            first_round_top_fitness_scaled = labels_one[labels_one['iteration'] == 1]["fitness_scaled"].max()
            first_round_top_variant = labels_one.loc[(labels_one['iteration'] == 1) & (labels_one["fitness_scaled"] == first_round_top_fitness_scaled), "variant"].values[0]
            first_round_fitness_binary_percentage = labels_one[labels_one['iteration'] == 1]["fitness_binary"].mean()
            first_round_median_fitness_scaled = labels_one[labels_one['iteration'] == 1]["fitness_scaled"].median()

            num_mutants_per_round_list = []
            first_round_strategy_list = []
            learning_strategy_list = []
            regression_type_list = []
            simulation_list =[]
            round_list = []
            test_error_list = []
            train_error_list = []
            train_r_squared_list = []
            test_r_squared_list = []
            alpha_list = []
            median_fitness_scaled_list = []
            top_fitness_scaled_list = []
            fitness_binary_percentage_list = []
            labels_list = []
            top_variant_list = []
            spearman_corr_list = []
            
            num_mutants_per_round_list.append(num_mutants_per_round)
            first_round_strategy_list.append(first_round_strategy)
            learning_strategy_list.append(learning_strategy)
            regression_type_list.append(regression_type)
            simulation_list.append(i)
            round_list.append(1)
            test_error_list.append("None")
            train_error_list.append("None")
            train_r_squared_list.append("None")
            test_r_squared_list.append("None")
            alpha_list.append("None")
            spearman_corr_list.append("None")
            median_fitness_scaled_list.append(first_round_median_fitness_scaled)
            top_fitness_scaled_list.append(first_round_top_fitness_scaled)
            top_variant_list.append(first_round_top_variant)
            fitness_binary_percentage_list.append(first_round_fitness_binary_percentage)
            labels_list.append(",".join(first_round_labels))

            labels_new = labels_one
            iteration_new = iteration_one

            for j in range(2, num_iterations + 1):
                iteration_old = iteration_new

                train_error, test_error, train_r_squared, test_r_squared, alpha, median_fitness_scaled, top_fitness_scaled, top_variant, fitness_binary_percentage, spearman_corr, df_test_new = top_layer(
                    iter_train=iteration_old['iteration'].unique().tolist(), iter_test=1001,
                    embeddings_pd=embeddings, labels_pd=labels_new,
                    measured_var=measured_var, regression_type=regression_type, top_n=top_n, final_round=final_round)
            
                num_mutants_per_round_list.append(num_mutants_per_round)
                first_round_strategy_list.append(first_round_strategy)
                learning_strategy_list.append(learning_strategy)
                regression_type_list.append(regression_type)
                simulation_list.append(i)
                round_list.append(j)
                test_error_list.append(test_error)
                train_error_list.append(train_error)
                train_r_squared_list.append(train_r_squared)
                test_r_squared_list.append(test_r_squared)
                alpha_list.append(alpha)
                median_fitness_scaled_list.append(median_fitness_scaled)
                top_fitness_scaled_list.append(top_fitness_scaled)
                top_variant_list.append(top_variant)
                fitness_binary_percentage_list.append(fitness_binary_percentage)
                spearman_corr_list.append(spearman_corr)

                # NOTE: work on alternate 2-n round strategies here
                if learning_strategy == 'dist':
                    iteration_new_ids = df_test_new.sort_values(by='dist_metric', ascending=False).head(num_mutants_per_round).variant
                elif learning_strategy == 'random':
                    iteration_new_ids = random.sample(list(df_test_new.variant), num_mutants_per_round)
                elif learning_strategy == 'top5bottom5':
                    iteration_new_ids = df_test_new.sort_values(by='y_pred', ascending=False).head(int(num_mutants_per_round/2)).variant
                    iteration_new_ids.append(df_test_new.sort_values(by='y_pred', ascending=False).tail(int(num_mutants_per_round/2)).variant)
                elif learning_strategy == 'top10':
                    iteration_new_ids = df_test_new.sort_values(by='y_pred', ascending=False).head(num_mutants_per_round).variant
                
                labels_list.append(",".join(iteration_new_ids))

                iteration_new = pd.DataFrame({'variant': iteration_new_ids, 'iteration': j})
                iteration_new = iteration_new.append(iteration_old)
                labels_new = pd.merge(labels, iteration_new, on='variant', how='left')
                labels_new.iteration[labels_new.iteration.isnull()] = 1001

            df_metrics = pd.DataFrame({'simulation_num':simulation_list,'round_num': round_list, 'num_mutants_per_round': num_mutants_per_round_list, 'first_round_strategy': first_round_strategy_list, 'learning_strategy': learning_strategy_list, 'regression_type': regression_type_list,
                                    'test_error': test_error_list, 'train_error': train_error_list,
                                    'train_r_squared': train_r_squared_list, 'test_r_squared': test_r_squared_list,
                                    'alpha': alpha_list, 'median_fitness_scaled': median_fitness_scaled_list,
                                    'top_fitness_scaled': top_fitness_scaled_list,
                                    'fitness_binary_percentage': fitness_binary_percentage_list, 'labels': labels_list, "top_variant": top_variant_list, "spearman_corr": spearman_corr_list})
            
            output_list.append(df_metrics)        
    
    output_table = pd.concat(output_list)
    return output_table

# Function to run the experiment with different combinations of parameters
def grid_search(dataset_name, base_path, num_simulations, num_iterations, measured_var, learning_strategies,
                   num_mutants_per_round, first_round_strategies, embedding_types, regression_types, file_type, embeddings_type_pt=None):
    
    # read in dataset
    embeddings, labels, hie_data = read_data(dataset_name, base_path, file_type, embeddings_type_pt)

    # scale embeddings
    embeddings_norm = scale_embeddings(embeddings)

    # generate embeddings_pca
    embeddings_pca = pca_embeddings(embeddings, labels, dataset_name, n_components=8)

    # save the embeddings in a list    
    embeddings_list = {
        'embeddings': embeddings,
        'embeddings_norm': embeddings_norm,
        'embeddings_pca': embeddings_pca
    }

    # save the labels in a list
    output_results = {}

    # Initialize the total combinations count
    total_combinations = 0 

    # Print the total number of combinations
    for strategy in learning_strategies:
        for var in measured_var:
            for iterations in num_iterations:
                for mutants_per_round in num_mutants_per_round:
                    if mutants_per_round == 128 and iterations != 3:
                        continue  # Skip other iterations when mutants_per_round is 128
                    for embedding_type in embedding_types:
                        for regression_type in regression_types:
                            for first_round_strategy in first_round_strategies:
                                total_combinations += 1

    # Print the corrected total_combinations count
    print(f"Total combinations: {total_combinations}")

    # Initialize the combination count
    output_list = []
    combination_count = 0

    start_time = time.time()

    for strategy in learning_strategies:
        output_results[strategy] = {}
        for var in measured_var:
            output_results[strategy][var] = {}
            for iterations in num_iterations:
                output_results[strategy][var][iterations] = {}
                for mutants_per_round in num_mutants_per_round:
                    output_results[strategy][var][iterations][mutants_per_round] = {}
                    for embedding_type in embedding_types:
                        output_results[strategy][var][iterations][mutants_per_round][embedding_type] = {}
                        for regression_type in regression_types:
                            output_results[strategy][var][iterations][mutants_per_round][embedding_type][regression_type] = {}
                            for first_round_strategy in first_round_strategies:
                                output_results[strategy][var][iterations][mutants_per_round][embedding_type][regression_type][first_round_strategy] = {}
                                if mutants_per_round == 128 and iterations != 3:
                                    continue  # Skip other iterations when mutants_per_round is 128
                                combination_count += 1
                                # print overall progress
                                print(
                                    f"Progress: {combination_count}/{total_combinations} "
                                    f"({(combination_count/total_combinations)*100:.2f}%)"
                                )

                                # run simulations for current combination of parameters
                                output_table = directed_evolution_simulation(
                                    labels=labels,
                                    embeddings=embeddings_list[embedding_type],
                                    num_simulations=num_simulations,
                                    hie_data=hie_data,
                                    num_iterations=iterations,
                                    num_mutants_per_round=mutants_per_round,
                                    measured_var=var,
                                    regression_type=regression_type,
                                    learning_strategy=strategy,
                                    final_round=mutants_per_round,
                                    first_round_strategy=first_round_strategy  
                                
                                )
                                output_list.append(output_table)


    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Total execution time: {execution_time:.2f} seconds")

    #concat the outputlist into a dataframe
    df_results = pd.concat(output_list)

    # save the dataframe to a csv file using the dataset_name
    if embeddings_type_pt == None:
        df_results.to_csv(f"{base_path}/results/{dataset_name}_results.csv", index=False)
    else:
        df_results.to_csv(f"{base_path}/results/{dataset_name}_{embeddings_type_pt}_results.csv", index=False)

def main():
    parser = create_parser()
    args = parser.parse_args()
    grid_search(
        args.dataset_name, args.base_path, args.num_simulations, args.num_iterations,
        args.measured_var, args.learning_strategies, args.num_mutants_per_round, args.first_round_strategies,
        args.embedding_types, args.regression_types, args.file_type, args.embeddings_type_pt
    )
 
if __name__ == "__main__":
    main()
