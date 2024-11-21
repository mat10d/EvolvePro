import os
import random
import time
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple

from evolvepro.src.data import load_dms_data, load_experimental_embeddings, load_experimental_data, create_iteration_dataframes
from evolvepro.src.utils import pca_embeddings
from evolvepro.src.model import first_round, top_layer

# Function to run the directed evolution simulation
def directed_evolution_simulation(
    labels: pd.DataFrame, 
    embeddings: pd.DataFrame, 
    num_simulations: int, 
    num_iterations: int, 
    num_mutants_per_round: int = 10, 
    measured_var: str = 'activity', 
    regression_type: str = 'ridge', 
    learning_strategy: str = 'topn', 
    top_n: int = None, 
    final_round: int = 10, 
    first_round_strategy: str = 'random', 
    embedding_type: str = None, 
    explicit_variants: Optional[List[str]] = None) -> pd.DataFrame:

    """
    Run the directed evolution simulation.

    Args:
    labels (pd.DataFrame): DataFrame of labels.
    embeddings (pd.DataFrame): DataFrame of embeddings.
    num_simulations (int): Number of simulations to run.
    num_iterations (int): Number of iterations to run.
    num_mutants_per_round (int): Number of mutants to select per round.
    measured_var (str): Measured variable.
    regression_type (str): Type of regression model.
    learning_strategy (str): Learning strategy.
    top_n (int): Number of top variants to consider.
    final_round (int): Number of final round mutants.
    first_round_strategy (str): First round strategy.
    embedding_type (str): Type of embeddings.
    explicit_variants (list): List of explicit variants.

    Returns:
    pd.DataFrame: DataFrame of simulation results.
    """

    # Initialize the output list of metrics
    output_list = []

    for i in range(1, num_simulations + 1):

        # Initialize the variables
        iteration_old = None
        num_mutants_per_round_list = []
        first_round_strategy_list = []
        measured_var_list = []
        learning_strategy_list = []
        regression_type_list = []
        embedding_type_list = []
        simulation_list =[]
        round_list = []
        test_error_list = []
        train_error_list = []
        train_r_squared_list = []
        test_r_squared_list = []
        alpha_list = []
        median_activity_scaled_list = []
        top_variant_list = []
        top_final_round_variants_list = []
        top_activity_scaled_list = []
        spearman_corr_list = []
        activity_binary_percentage_list = []
        
        # Initialize the list of variants for each round
        this_round_variants_list = []
        next_round_variants_list = []

        j = 0    
        while j <= num_iterations:
            # Perform mutant selection for the first round
            if j == 0:
                labels_new, iteration_new, this_round_variants = first_round(
                    labels, 
                    embeddings, 
                    explicit_variants=explicit_variants,
                    num_mutants_per_round=num_mutants_per_round, 
                    first_round_strategy=first_round_strategy, 
                    embedding_type=embedding_type, 
                    random_seed=i
                )

                # Append the results to the output list
                num_mutants_per_round_list.append(num_mutants_per_round)
                first_round_strategy_list.append(first_round_strategy)
                measured_var_list.append(measured_var)
                learning_strategy_list.append(learning_strategy)
                regression_type_list.append(regression_type)
                embedding_type_list.append(embedding_type)                
                simulation_list.append(i)
                round_list.append(j)
                # Append None values for the metrics for the first round
                test_error_list.append("None")
                train_error_list.append("None")
                train_r_squared_list.append("None")
                test_r_squared_list.append("None")
                alpha_list.append("None")
                median_activity_scaled_list.append("None")
                top_activity_scaled_list.append("None")
                top_variant_list.append("None")
                top_final_round_variants_list.append("None")
                activity_binary_percentage_list.append("None")
                spearman_corr_list.append("None")
                # Append the variants for the first round, round 0 will have None
                this_round_variants_list.append("None")
                next_round_variants_list.append(",".join(this_round_variants))

                j += 1

            else:
                # Perform mutant selection for the subsequent rounds
                iteration_old = iteration_new
                print("iterations considered", iteration_old)
                
                train_error, test_error, train_r_squared, test_r_squared, alpha, median_activity_scaled, top_activity_scaled, top_variant, top_final_round_variants, activity_binary_percentage, spearman_corr, df_test_new, this_round_variants = top_layer(
                    iter_train=iteration_old['iteration'].unique().tolist(), iter_test=None,
                    embeddings_pd=embeddings, labels_pd=labels_new,
                    measured_var=measured_var, regression_type=regression_type, top_n=top_n, final_round=final_round)
                # Perform mutant selection for the next round based on the results of the current round
                if learning_strategy == 'dist':
                    iteration_new_ids = df_test_new.sort_values(by='dist_metric', ascending=False).head(num_mutants_per_round).variant
                elif learning_strategy == 'random':
                    iteration_new_ids = random.sample(list(df_test_new.variant), num_mutants_per_round)
                elif learning_strategy == 'topn2bottomn2':
                    top_half = df_test_new.sort_values(by='y_pred', ascending=False).head(int(num_mutants_per_round / 2)).variant
                    bottom_half = df_test_new.sort_values(by='y_pred', ascending=False).tail(int(num_mutants_per_round / 2)).variant
                    iteration_new_ids = pd.concat([top_half, bottom_half])
                elif learning_strategy == 'topn':
                    iteration_new_ids = df_test_new.sort_values(by='y_pred', ascending=False).head(num_mutants_per_round).variant

                iteration_new = pd.DataFrame({'variant': iteration_new_ids, 'iteration': j})
                iteration_new = pd.concat([iteration_new, iteration_old], ignore_index=True)
                labels_new = pd.merge(labels, iteration_new, on='variant', how='left')

                num_mutants_per_round_list.append(num_mutants_per_round)
                first_round_strategy_list.append(first_round_strategy)
                measured_var_list.append(measured_var)
                learning_strategy_list.append(learning_strategy)
                regression_type_list.append(regression_type)
                embedding_type_list.append(embedding_type)
                simulation_list.append(i)
                round_list.append(j)
                test_error_list.append(test_error)
                train_error_list.append(train_error)
                train_r_squared_list.append(train_r_squared)
                test_r_squared_list.append(test_r_squared)
                alpha_list.append(alpha)
                median_activity_scaled_list.append(median_activity_scaled)
                top_activity_scaled_list.append(top_activity_scaled)
                top_variant_list.append(top_variant)
                top_final_round_variants_list.append(top_final_round_variants)
                activity_binary_percentage_list.append(activity_binary_percentage)
                spearman_corr_list.append(spearman_corr)

                this_round_variants_list.append(",".join(iteration_old.variant))
                next_round_variants_list.append(",".join(iteration_new_ids))

                j += 1

            df_metrics = pd.DataFrame({
                'simulation_num': simulation_list, 
                'round_num': round_list, 
                'num_mutants_per_round': num_mutants_per_round_list, 
                'first_round_strategy': first_round_strategy_list, 
                'measured_var': measured_var_list, 
                'learning_strategy': learning_strategy_list, 
                'regression_type': regression_type_list,
                'embedding_type': embedding_type_list,
                'test_error': test_error_list, 
                'train_error': train_error_list,
                'train_r_squared': train_r_squared_list, 
                'test_r_squared': test_r_squared_list,
                'alpha': alpha_list, 
                "spearman_corr": spearman_corr_list,
                'median_activity_scaled': median_activity_scaled_list, 
                'top_activity_scaled': top_activity_scaled_list, 
                'activity_binary_percentage': activity_binary_percentage_list, 
                "top_variant": top_variant_list, 
                "top_final_round_variants": top_final_round_variants_list, 
                "this_round_variants": this_round_variants_list, 
                "next_round_variants": next_round_variants_list
            })

        output_list.append(df_metrics)


    output_table = pd.concat(output_list)
    return output_table

# Function to run the experiment with different combinations of parameters 
def grid_search(
    dataset_name: str,
    experiment_name: str,
    model_name: str,
    embeddings_path: str,
    labels_path: str,
    num_simulations: int,
    num_iterations: List[int],
    measured_var: List[str],
    learning_strategies: List[str],
    num_mutants_per_round: List[int],
    num_final_round_mutants: int,
    first_round_strategies: List[str],
    embedding_types: List[str],
    pca_components: List[int],
    regression_types: List[str],
    embeddings_file_type: str,
    output_dir: str,
    embeddings_type_pt: Optional[str] = None) -> None:

    """
    Run the experiment with different combinations of parameters.

    Args:
    dataset_name (str): Name of the dataset.
    experiment_name (str): Name of the experiment.
    model_name (str): Name of the model.
    embeddings_path (str): Path to the embeddings file.
    labels_path (str): Path to the labels file.
    num_simulations (int): Number of simulations to run.
    num_iterations (List[int]): List of number of iterations to run.
    measured_var (List[str]): List of measured variables.
    learning_strategies (List[str]): List of learning strategies.
    num_mutants_per_round (List[int]): List of number of mutants to select per round.
    num_final_round_mutants (int): Number of final round mutants.
    first_round_strategies (List[str]): List of first round strategies.
    embedding_types (List[str]): List of embedding types.
    pca_components (List[int]): List of PCA components.
    regression_types (List[str]): List of regression types.
    embeddings_file_type (str): Type of embeddings file.
    output_dir (str): Directory to save output files.
    embeddings_type_pt (str): Type of embeddings file (PyTorch).

    Returns:
    None
    """

    # Load dataset
    embeddings, labels = load_dms_data(dataset_name, model_name, embeddings_path, labels_path, embeddings_file_type, embeddings_type_pt)
    
    if embeddings is None or labels is None:
        print("Failed to load data. Exiting.")
        return

    # Generate pca components only if pca_components is not None
    if pca_components is not None:
        embeddings_pca = {
            f'embeddings_pca_{n}': pca_embeddings(embeddings, n_components=n)
            for n in pca_components
        }
        embeddings_list = {
            'embeddings': embeddings,
            **embeddings_pca
        }
    else:
        embeddings_list = {
            'embeddings': embeddings
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
                                combination_count += 1
                                # print overall progress

                                # run simulations for current combination of parameters
                                output_table = directed_evolution_simulation(
                                    labels=labels,
                                    embeddings=embeddings_list[embedding_type],
                                    num_simulations=num_simulations,
                                    num_iterations=iterations,
                                    num_mutants_per_round=mutants_per_round,
                                    measured_var=var,
                                    regression_type=regression_type,
                                    learning_strategy=strategy,
                                    final_round=num_final_round_mutants,
                                    first_round_strategy=first_round_strategy,
                                    embedding_type=embedding_type
                                )
                                print(
                                    f"Progress: {combination_count}/{total_combinations} "
                                    f"({(combination_count/total_combinations)*100:.2f}%)"
                                )

                                output_list.append(output_table)


    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Total execution time: {execution_time:.2f} seconds")

    #concat the outputlist into a dataframe
    df_results = pd.concat(output_list)

    # make the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # save the dataframe to a csv file using the dataset_name
    if embeddings_type_pt == None:
        df_results.to_csv(f"{output_dir}/{dataset_name}_{model_name}_{experiment_name}.csv", index=False)
    else:
        df_results.to_csv(f"{output_dir}/{dataset_name}_{model_name}_{experiment_name}_{embeddings_type_pt}.csv", index=False)

def evolve_experimental(
    protein_name : str,
    round_name : str,
    embeddings_base_path : str,
    embeddings_file_name : str,
    round_base_path : str,
    round_file_names : List[str],
    wt_fasta_path : str,
    rename_WT : bool = False,
    number_of_variants : int = 12,
    output_dir : str = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_results/'
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: 

    """
    Perform one round of directed evolution for a protein.

    Args:
    protein_name (str): Name of the protein.
    round_name (str): Name of the current round (e.g., 'Round1').
    embeddings_base_path (str): Base path for embeddings file.
    embeddings_file_name (str): Name of the embeddings file.
    round_base_path (str): Base path for round data files.
    round_file_names (list): List of round file names.
    wt_fasta_path (str): Path to the wild-type FASTA file.
    rename_WT (bool): Whether to rename the wild-type.
    number_of_variants (int): Number of top variants to display.
    output_dir (str): Directory to save output files.

    Returns:
    tuple: (this_round_variants, df_test, df_sorted_all)
    """
    
    print(f"Processing {protein_name} - {round_name}")
    
    # Load embeddings
    embeddings = load_experimental_embeddings(embeddings_base_path, embeddings_file_name, rename_WT)
    print(f"Embeddings loaded: {embeddings.shape}")
    
    # Load experimental data
    all_experimental_data = []
    for round_file_name in round_file_names:
        experimental_data = load_experimental_data(round_base_path, round_file_name, wt_fasta_path)
        all_experimental_data.append(experimental_data)
        print(f"Loaded experimental data for {round_file_name}: {experimental_data.shape}")
    
    # Create iteration dataframes
    iteration, labels = create_iteration_dataframes(all_experimental_data, embeddings.index)
    print(f"iteration shape: {iteration.shape}")
    print(f"Labels shape: {labels.shape}")
    
    # Perform top layer analysis
    this_round_variants, df_test, df_sorted_all = top_layer(
        iter_train=iteration['iteration'].unique().tolist(),
        iter_test=None,
        embeddings_pd=embeddings,
        labels_pd=labels,
        measured_var='activity',
        regression_type='randomforest',
        experimental=True
    )
    
    # Print results
    print(f"\nTested variants in this round: {len(this_round_variants)}")
    print(this_round_variants)
    print(f"\nTop {number_of_variants} variants predicted by the model:")
    print(df_test.sort_values(by=['y_pred'], ascending=False).head(number_of_variants))
    
    # Save results if an output_dir is provided
    if output_dir is not None:
        output_dir = os.path.join(output_dir, protein_name, round_name)
        os.makedirs(output_dir, exist_ok=True)
        iteration.to_csv(os.path.join(output_dir, 'iteration.csv'))
        this_round_variants.to_csv(os.path.join(output_dir, 'this_round_variants.csv'))
        df_test = df_test.sort_values(by=['y_pred'], ascending=False)
        df_test.to_csv(os.path.join(output_dir, 'df_test.csv'))
        df_sorted_all.to_csv(os.path.join(output_dir, 'df_sorted_all.csv'))
        print(f"\nData saved to {output_dir}")
    
    return this_round_variants, df_test, df_sorted_all

def evolve_experimental_multi(
    protein_name: str,
    round_name: str,
    embeddings_base_path: str,
    embeddings_file_names: List[str],
    round_base_path: str,
    round_file_names_single: List[str],
    round_file_names_multi: List[str],
    wt_fasta_path: str,
    rename_WT: bool = False,
    number_of_variants: int = 12,
    output_dir: str = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_results/'
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Perform one round of directed evolution for a protein with multi-mutant support.
    """
    print(f"Processing {protein_name} - {round_name}")
    
    # Load and concatenate multiple embedding files
    embeddings_list = []
    for i, file_name in enumerate(embeddings_file_names):
        embedding = load_experimental_embeddings(embeddings_base_path, file_name, rename_WT)
        if i > 0:  # If not the first file
            embedding = embedding[~embedding.index.isin(['WT', 'WT Wild-type sequence'])]
        embeddings_list.append(embedding)
    
    embeddings = pd.concat(embeddings_list)
    print(f"Embeddings loaded: {embeddings.shape}")
    
    # Load experimental data
    all_experimental_data = []
    for round_file_name in round_file_names_single:
        experimental_data = load_experimental_data(round_base_path, round_file_name, wt_fasta_path, single_mutant=True)
        all_experimental_data.append(experimental_data)
        print(f"Loaded experimental data for {round_file_name}: {experimental_data.shape}")

    for round_file_name in round_file_names_multi:
        experimental_data = load_experimental_data(round_base_path, round_file_name, wt_fasta_path, single_mutant=False)
        all_experimental_data.append(experimental_data)
        print(f"Loaded experimental data for {round_file_name}: {experimental_data.shape}")
    
    # Create iteration dataframes
    iteration, labels = create_iteration_dataframes(all_experimental_data, embeddings.index)
    
    # Perform top layer analysis
    this_round_variants, df_test, df_sorted_all = top_layer(
        iter_train=iteration['iteration'].unique().tolist(),
        iter_test=None,
        embeddings_pd=embeddings,
        labels_pd=labels,
        measured_var='activity',
        regression_type='randomforest',
        experimental=True
    )
    
    # Print results
    print(f"\nTested variants in this round: {len(this_round_variants)}")
    print(this_round_variants)
    print(f"\nTop {number_of_variants} variants predicted by the model:")
    print(df_test.sort_values(by=['y_pred'], ascending=False).head(number_of_variants)[['variant', 'y_pred']])
    
    # Save results if an output_dir is provided
    if output_dir is not None:
        output_dir = os.path.join(output_dir, protein_name, round_name)
        os.makedirs(output_dir, exist_ok=True)
        iteration.to_csv(os.path.join(output_dir, 'iteration.csv'))
        this_round_variants.to_csv(os.path.join(output_dir, 'this_round_variants.csv'))
        df_test = df_test.sort_values(by=['y_pred'], ascending=False)
        df_test.to_csv(os.path.join(output_dir, 'df_test.csv'))
        df_sorted_all.to_csv(os.path.join(output_dir, 'df_sorted_all.csv'))
        print(f"\nData saved to {output_dir}")
    
    return this_round_variants, df_test, df_sorted_all