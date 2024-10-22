import pandas as pd
import os
from evolvepro.src.plot import *

# Initialize paths and datasets
base_dir = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro"
data_dir = os.path.join(base_dir, "output/dms_results/")
output_dir = os.path.join(base_dir, "output/dms_plots/")
datasets = [
    "cov2_S", "cas12f", "zikv_E", "kelsic", "brenan", "stiffler", "markin",
    "giacomelli", "jones", "haddox", "doud", "lee"
]

# One shot comparison analysis
esm2_15B_one_shot = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_one_shot",
    group_columns=["num_mutants_per_round", "round_num"],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

esm2_15B_one_shot = apply_labels(
    esm2_15B_one_shot,
    column='label',
    prefix='pre-training: ',
    suffix=' mutants',
    value_column='num_mutants_per_round'
)

esm2_15B_optimal = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_optimal",
    group_columns=["num_mutants_per_round", "round_num"],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

esm2_15B_optimal = filter_dataframe(esm2_15B_optimal, conditions={"round_num": [5, 10]})
esm2_15B_optimal = apply_labels(
    esm2_15B_optimal,
    column='label',
    prefix='EVOLVEPro: ',
    suffix=' rounds',
    value_column='round_num'
)

# Load external comparison data
background = load_external_data(
    os.path.join(base_dir, "output/dms/background_df.csv"),
    label="background",
    rename_columns={'defined': 'fitness_binary_percentage_mean'}
)
efficient_evolution = load_external_data(
    os.path.join(base_dir, "output/dms_results/external_data/efficient-evolution.csv"),
    label="efficient-evolution default"
)
efficient_evolution_wider = load_external_data(
    os.path.join(base_dir, "output/dms_results/external_data/efficient-evolution_wider.csv"),
    label="efficient-evolution expanded"
)

# Combine and plot one shot comparison
dataframes = [esm2_15B_one_shot, esm2_15B_optimal, background, efficient_evolution, efficient_evolution_wider]
one_shot = concatenate_dataframes(
    dataframes,
    output_dir=output_dir,
    output_file="one_shot.csv"
)
plot_comparison(
    one_shot,
    title="Comparison to one shot strategies",
    output_dir=output_dir,
    output_file="one_shot"
)

# Base model comparison
models = [
    ("esm2_t48_15B_UR50D", "esm2_15B_optimal"),
    ("esm2_t36_3B_UR50D", "esm2_3B_optimal"),
    ("esm2_t33_650M_UR50D", "esm2_650M_optimal"),
    ("esm1b_t33_650M_UR50S", "esm1b_650M_optimal"),
    ("esm1v_t33_650M_UR90S_1", "esm1v_650M_optimal"),
    ("one_hot_encoded", "one_hot_encoded_optimal"),
    ("integer_encoded", "integer_encoded_optimal"),
    ("proteinbert", "proteinbert_optimal"),
    ("prot_t5", "prot_t5_optimal"),
    ("ankh_large", "ankh_large_optimal"),
    ("ankh_base", "ankh_base_optimal"),
    ("unirep", "unirep_optimal")
]

processed_dfs = []
for model, experiment in models:
    df = read_dms_data(
        directory=data_dir,
        datasets=datasets,
        model=model,
        experiment=experiment,
        group_columns=["num_mutants_per_round", "round_num"],
        aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
    )
    df = apply_labels(df, column='label', value_column='model')
    processed_dfs.append(df)

base_model = concatenate_dataframes(
    processed_dfs,
    output_dir=output_dir,
    output_file="base_model.csv"
)

base_model_5 = filter_dataframe(
    base_model,
    conditions={"round_num": [5, 10]},
    output_dir=output_dir,
    output_file="base_model_5.csv"
)

base_model_10 = filter_dataframe(
    base_model,
    conditions={"round_num": [5, 10]},
    output_dir=output_dir,
    output_file="base_model_10.csv"
)

plot_comparison(
    base_model_5,
    title="Base model comparison, 5 rounds of evolution",
    output_dir=output_dir,
    output_file="base_model_5"
)

plot_comparison(
    base_model_10,
    title="Base model comparison, 10 rounds of evolution",
    output_dir=output_dir,
    output_file="base_model_10"
)

# Grid search comparison
esm2_15B_grid = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_full",
    group_columns=[
        "num_mutants_per_round", "round_num", "first_round_strategy",
        "measured_var", "learning_strategy", "regression_type", "embedding_type"
    ],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

esm2_15B_grid_5 = filter_dataframe(esm2_15B_grid, conditions={"round_num": [5]})
esm2_15B_grid_10 = filter_dataframe(esm2_15B_grid, conditions={"round_num": [10]})

# Plot grid search results
for strategy in ['first_round_strategy', 'measured_var', 'learning_strategy', 'regression_type', 'embedding_type']:
    for rounds, df in [('5', esm2_15B_grid_5), ('10', esm2_15B_grid_10)]:
        plot_grid_search_bar(
            df=df,
            variable='fitness_binary_percentage_mean',
            strategy_column=strategy,
            title=f'{strategy} ({rounds} rounds)',
            output_dir=output_dir,
            output_file=f'{strategy}_{rounds}'
        )

# Grid search heatmaps
for rounds, df in [('5', esm2_15B_grid_5), ('10', esm2_15B_grid_10)]:
    plot_grid_search_heatmap(
        df=df,
        variable='fitness_binary_percentage_mean',
        strategy_columns=['regression_type', 'learning_strategy'],
        output_dir=output_dir,
        output_file=f'regression_type_learning_strategy_{rounds}'
    )

# Dimensionality reduction comparison
esm2_15B_grid_pca = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_pca",
    group_columns=["num_mutants_per_round", "round_num", "embedding_type"],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

esm2_15B_grid_pca_5 = filter_dataframe(esm2_15B_grid_pca, conditions={"round_num": [5]})
esm2_15B_grid_pca_10 = filter_dataframe(esm2_15B_grid_pca, conditions={"round_num": [10]})

for rounds, df in [('5', esm2_15B_grid_pca_5), ('10', esm2_15B_grid_pca_10)]:
    plot_grid_search_bar(
        df=df,
        variable='fitness_binary_percentage_mean',
        strategy_column='embedding_type',
        title=f'embedding_type ({rounds} rounds)',
        output_dir=output_dir,
        output_file=f'embedding_type_pca_{rounds}'
    )

# Across rounds comparison
esm2_15B_optimal = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_optimal",
    group_columns=["num_mutants_per_round", "round_num"],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

save_dataframe(
    esm2_15B_optimal,
    output_dir=output_dir,
    output_file="esm2_15B_optimal.csv"
)

plot_by_round(
    esm2_15B_optimal,
    variable='fitness_binary_percentage_mean',
    output_dir=output_dir,
    output_file='esm2_15B_optimal'
)

# Number of variants comparison
esm2_15B_num_variants = read_dms_data(
    directory=data_dir,
    datasets=datasets,
    model="esm2_t48_15B_UR50D",
    experiment="esm2_15B_num_variants",
    group_columns=["num_mutants_per_round", "round_num"],
    aggregate_columns=["median_fitness_scaled", "top_fitness_scaled", "fitness_binary_percentage"]
)

plot_by_round_split(
    esm2_15B_num_variants,
    variable='fitness_binary_percentage_mean',
    split_variable='num_mutants_per_round',
    output_dir=output_dir,
    output_file='esm2_15B_optimal'
)