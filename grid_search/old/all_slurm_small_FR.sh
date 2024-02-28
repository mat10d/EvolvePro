#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=all_slurm_small_FR
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=4gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/Github/directed_evolution/grid_search/out/all_slurm_small_FR-%j.out 

source ~/.bashrc
conda activate embeddings

#datasets=("brenan_esm2_t33_650M_UR50D" "stiffler_esm2_t33_650M_UR50D" "doud_esm2_t33_650M_UR50D" "haddox_esm2_t33_650M_UR50D" "giacomelli_esm2_t33_650M_UR50D" "jones_esm2_t33_650M_UR50D" "kelsic_esm2_t33_650M_UR50D" "lee_esm2_t33_650M_UR50D" "markin_esm2_t33_650M_UR50D")
datasets=("markin_esm2_t33_650M_UR50D")
num_simulations=3
num_iterations=(2 3 4 5 6 7 8 9 10 11)
measured_var="fitness"
learning_strategies="top10"
num_mutants_per_round=16
first_round_strategies=("random" "diverse_medoids" "representative_hie")
embedding_types="embeddings"
regression_types="randomforest"
file_type="csvs"

# Function to run grid search for a given dataset
function run_grid_search() {
    dataset_name=$1

    echo "Running ${dataset_name} dataset:" > out/${dataset_name}-${first_round_strategies}-hc_small_FR.out
    python3 -u grid_search.py \
        --dataset_name ${dataset_name} \
        --base_path ../esm-extract/results_means \
        --num_simulations ${num_simulations} \
        --num_iterations ${num_iterations[*]} \
        --measured_var ${measured_var} \
        --learning_strategies ${learning_strategies} \
        --num_mutants_per_round ${num_mutants_per_round} \
        --first_round_strategies ${first_round_strategies[*]} \
        --embedding_types ${embedding_types} \
        --regression_types ${regression_types} \
        --file_type ${file_type} \
        >> out/${dataset_name}-${first_round_strategies}-hc_small_FR.out
    echo "Done running ${dataset_name} dataset:" >> out/${dataset_name}-${first_round_strategies}-hc_small_FR.out
}

# Loop over datasets
for dataset in "${datasets[@]}"; do
    run_grid_search ${dataset}
done
