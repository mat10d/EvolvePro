#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=7-00:00:00
#SBATCH --job-name=esm2_15B_full_grid_search_distributed
#SBATCH -n 1
#SBATCH -N 1   
#SBATCH --cpus-per-task=12  # Set the number of CPU cores
#SBATCH --mem=100gb  
#SBATCH -p abugoot
#SBATCH --output /orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/grid_search/out/esm2_15B_full_grid_search_distributed-%j.out 

source ~/.bashrc
conda activate embeddings
module load openmind8/gnu-parallel/20240222

names=("brenan" "stiffler" "doud" "haddox" "giacomelli" "jones" "kelsic" "lee" "markin" "zikv_E" "cas12f" "cov2_S")
datasets=()

for name in "${names[@]}"; do
    datasets+=("${name}_esm2_t48_15B_UR50D")
done

# Function to run grid search for a given dataset
function run_grid_search() {
    dataset_name=$1

    echo "Running ${dataset_name} dataset:" > out/${dataset_name}-average-esm2_15B_full_grid_search_distributed.out
    python3 -u grid_search.py \
        --dataset_name ${dataset_name} \
        --experiment_name "esm2_15B_full_grid_search_distributed" \
        --base_path ../extract/esm/results_means \
        --num_simulations 3 \
        --num_iterations 10 \
        --measured_var "fitness" "fitness_scaled" \
        --learning_strategies "random" "top5bottom5" "top10" "dist" \
        --num_mutants_per_round 16 \
        --num_final_round_mutants 16 \
        --first_round_strategies "random" "diverse_medoids" \
        --embedding_types "embeddings" "embeddings_pca" \
        --regression_types "ridge" "lasso" "elasticnet" "linear" "neuralnet" "randomforest" "gradientboosting" \
        --file_type "csvs" \
        --embeddings_type_pt "average" \
        >> out/${dataset_name}-average-esm2_15B_full_grid_search_distributed.out
    echo "Done running ${dataset_name} dataset:" >> out/${dataset_name}-average-esm2_15B_full_grid_search_distributed.out
}

# Export the function so it's available to GNU Parallel
export -f run_grid_search

# Use GNU Parallel to run the grid search function in parallel for each dataset
parallel -j12 run_grid_search ::: "${datasets[@]}"
