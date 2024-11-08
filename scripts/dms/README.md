# Deep Mutational Scanning (DMS) Workflow for EVOLVEpro

This directory contains scripts for running the DMS workflow of EVOLVEpro, which optimizes a few-shot model on deep mutational scanning datasets. The source functions used here are in `evolvepro/src/evolve.py`, which calls the underlying model, in `evolvepro/src/model.py`.

## Main Script: dms_main.py

The main script `dms_main.py` is used to run experiments with different combinations of grid search variables.

### Usage

For dms, use the evolvepro environment:

```bash
conda activate evolvepro
python dms_main.py [arguments]
```

### Arguments

- `--dataset_name`: Name of the dataset
- `--experiment_name`: Name of the experiment
- `--model_name`: Name of the model used for embeddings
- `--embeddings_path`: Path to the embeddings file
- `--labels_path`: Path to the labels file
- `--num_simulations`: Number of simulations for each parameter combination
- `--num_iterations`: List of number of iterations (must be greater than 1)
- `--measured_var`: Fitness type to train on (options: fitness, fitness_scaled)
- `--learning_strategies`: Type of learning strategy (options: random, top5bottom5, top10, dist)
- `--num_mutants_per_round`: Number of mutants per round
- `--num_final_round_mutants`: Number of mutants in final round
- `--first_round_strategies`: Type of first round strategy (options: random, diverse_medoids, representative_hie)
- `--embedding_types`: Types of embeddings to train on (options: embeddings, embeddings_pca)
- `--pca_components`: Number of PCA components to use
- `--regression_types`: Regression types (options: ridge, lasso, elasticnet, linear, neuralnet, randomforest, gradientboosting, knn, gp)
- `--embeddings_file_type`: Type of embeddings file to read (options: csv, pt)
- `--output_dir`: Output directory for grid search results
- `--embeddings_type_pt`: Type of embeddings to use for .pt files (options: average, mutated, both)

## Example SLURM Scripts

We provide two example SLURM scripts for different use cases:

- `esm_2_15B_full.sh`: Performs a comprehensive grid search across multiple datasets and parameters.
- `esm2_15B.sh`: Runs simulations with a specific set of optimal parameters across multiple datasets.

## Running the Scripts

1. Ensure you have the necessary dependencies installed and the correct conda environment activated.
2. Submit the SLURM job using:
   ```bash
   sbatch esm_2_15B_full.sh
   ```
   or
   ```bash
   sbatch esm2_15B.sh
   ```
3. The results will be saved in the specified output directory.

Note: These scripts make use of GNU Parallel to run multiple datasets in parallel. Make sure you have it installed and loaded in your environment.