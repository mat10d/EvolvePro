import os
import random
import time
import pandas as pd
from typing import List, Dict, Any, Optional
from evolvepro.src.data import load_dms_data
from evolvepro.src.utils import pca_embeddings
from evolvepro.src.model import first_round, top_layer

def directed_evolution_simulation(
    labels: pd.DataFrame, 
    embeddings: pd.DataFrame, 
    num_simulations: int, 
    num_iterations: int, 
    num_mutants_per_round: int = 10, 
    measured_var: str = 'fitness', 
    regression_type: str = 'ridge', 
    learning_strategy: str = 'top10', 
    top_n: int = None, 
    final_round: int = 10, 
    first_round_strategy: str = 'random', 
    embedding_type: str = None, 
    explicit_variants: Optional[List[str]] = None) -> pd.DataFrame:

    output_list = []

    for i in range(1, num_simulations + 1):
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
        median_fitness_scaled_list = []
        top_variant_list = []
        top_final_round_variants_list = []
        top_fitness_scaled_list = []
        spearman_corr_list = []
        fitness_binary_percentage_list = []
        
        this_round_variants_list = []
        next_round_variants_list = []

        j = 0    
        while j <= num_iterations:
            # Perform mutant selection for the current round
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
                num_mutants_per_round_list.append(num_mutants_per_round)
                first_round_strategy_list.append(first_round_strategy)
                measured_var_list.append(measured_var)
                learning_strategy_list.append(learning_strategy)
                regression_type_list.append(regression_type)
                embedding_type_list.append(embedding_type)                
                simulation_list.append(i)
                round_list.append(j)
                test_error_list.append("None")
                train_error_list.append("None")
                train_r_squared_list.append("None")
                test_r_squared_list.append("None")
                alpha_list.append("None")
                median_fitness_scaled_list.append("None")
                top_fitness_scaled_list.append("None")
                top_variant_list.append("None")
                top_final_round_variants_list.append("None")
                fitness_binary_percentage_list.append("None")
                spearman_corr_list.append("None")

                this_round_variants_list.append("None")
                next_round_variants_list.append(",".join(this_round_variants))

                j += 1

            else:
                iteration_old = iteration_new
                print("iterations considered", iteration_old)

                train_error, test_error, train_r_squared, test_r_squared, alpha, median_fitness_scaled, top_fitness_scaled, top_variant, top_final_round_variants, fitness_binary_percentage, spearman_corr, df_test_new, this_round_variants = top_layer(
                    iter_train=iteration_old['iteration'].unique().tolist(), iter_test=None,
                    embeddings_pd=embeddings, labels_pd=labels_new,
                    measured_var=measured_var, regression_type=regression_type, top_n=top_n, final_round=final_round)
                # Perform mutant selection for the next round based on the results of the current round
                if learning_strategy == 'dist':
                    iteration_new_ids = df_test_new.sort_values(by='dist_metric', ascending=False).head(num_mutants_per_round).variant
                elif learning_strategy == 'random':
                    iteration_new_ids = random.sample(list(df_test_new.variant), num_mutants_per_round)
                elif learning_strategy == 'top5bottom5':
                    iteration_new_ids = df_test_new.sort_values(by='y_pred', ascending=False).head(int(num_mutants_per_round / 2)).variant
                    iteration_new_ids.append(df_test_new.sort_values(by='y_pred', ascending=False).tail(int(num_mutants_per_round / 2)).variant)
                elif learning_strategy == 'top10':
                    iteration_new_ids = df_test_new.sort_values(by='y_pred', ascending=False).head(num_mutants_per_round).variant

                iteration_new = pd.DataFrame({'variant': iteration_new_ids, 'iteration': j})
                iteration_new = iteration_new.append(iteration_old)
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
                median_fitness_scaled_list.append(median_fitness_scaled)
                top_fitness_scaled_list.append(top_fitness_scaled)
                top_variant_list.append(top_variant)
                top_final_round_variants_list.append(top_final_round_variants)
                fitness_binary_percentage_list.append(fitness_binary_percentage)
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
                'median_fitness_scaled': median_fitness_scaled_list, 
                'top_fitness_scaled': top_fitness_scaled_list, 
                'fitness_binary_percentage': fitness_binary_percentage_list, 
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

    # Load dataset
    embeddings, labels = load_dms_data(dataset_name, model_name, embeddings_path, labels_path, embeddings_file_type, embeddings_type_pt)
    
    if embeddings is None or labels is None:
        print("Failed to load data. Exiting.")
        return

    # Generate PCA embeddings for each n_components
    embeddings_pca = {
        f'embeddings_pca_{n}': pca_embeddings(embeddings, n_components=n)
        for n in pca_components
    }

    # Save the embeddings in a dictionary
    embeddings_list = {
        'embeddings': embeddings,
        **embeddings_pca
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