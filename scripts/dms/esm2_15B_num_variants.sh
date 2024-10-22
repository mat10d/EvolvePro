#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=1:00:00 
#SBATCH --job-name=esm2_15B_num_variants
#SBATCH -n 12
#SBATCH -N 1   
#SBATCH --cpus-per-task=5  
#SBATCH --mem=50gb  
#SBATCH --output out/esm2_15B_num_variants-%j.out 

source ~/.bashrc
conda activate embeddings
module load openmind8/gnu-parallel/20240222

datasets=("brenan" "doud" "haddox" "giacomelli" "jones" "lee" "zikv_E" "cas12f")

# Function to run dms_main for a given dataset
run_dms_main_500() {
    dataset_name=$1
    output_file="out/${dataset_name}-esm2_15B_num_variants.out"

    echo "Running ${dataset_name} dataset:" > ${output_file}
    python3 -u dms_main.py \
        --dataset_name ${dataset_name} \
        --experiment_name "esm2_15B_num_variants" \
        --model_name "esm2_t48_15B_UR50D" \
        --embeddings_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm" \
        --labels_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms" \
        --num_simulations 3 \
        --num_iterations 10 \
        --measured_var "fitness" \
        --learning_strategies "top10" \
        --num_mutants_per_round 10 20 30 40 50 100 200 500 \
        --num_final_round_mutants 16 \
        --first_round_strategies "random" \
        --embedding_types "embeddings" \
        --regression_types "randomforest" \
        --embeddings_file_type "csv" \
        --output_dir "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms_results" \
        >> ${output_file} 2>&1
    echo "Done running ${dataset_name} dataset:" >> ${output_file}
}

# Export the function so it's available to GNU Parallel
export -f run_dms_main_500

# Use GNU Parallel to run the dms_main function in parallel for each dataset
parallel -j9 run_dms_main_500 ::: "${datasets[@]}"

datasets=("stiffler")

# Function to run dms_main for a given dataset
run_dms_main_200() {
    dataset_name=$1
    output_file="out/${dataset_name}-esm2_15B_num_variants.out"

    echo "Running ${dataset_name} dataset:" > ${output_file}
    python3 -u dms_main.py \
        --dataset_name ${dataset_name} \
        --experiment_name "esm2_15B_num_variants" \
        --model_name "esm2_t48_15B_UR50D" \
        --embeddings_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm" \
        --labels_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms" \
        --num_simulations 10 \
        --num_iterations 1 \
        --measured_var "fitness" \
        --learning_strategies "top10" \
        --num_mutants_per_round 10 20 30 40 50 100 200 \
        --num_final_round_mutants 16 \
        --first_round_strategies "random" \
        --embedding_types "embeddings" \
        --regression_types "randomforest" \
        --embeddings_file_type "csv" \
        --output_dir "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms_results" \
        >> ${output_file} 2>&1
    echo "Done running ${dataset_name} dataset:" >> ${output_file}
}

# Export the function so it's available to GNU Parallel
export -f run_dms_main_200

# Use GNU Parallel to run the dms_main function in parallel for each dataset
parallel -j1 run_dms_main_200 ::: "${datasets[@]}"

datasets=("kelsic" "markin" "cov2_S")

# Function to run dms_main for a given dataset
run_dms_main_100() {
    dataset_name=$1
    output_file="out/${dataset_name}-esm2_15B_num_variants.out"

    echo "Running ${dataset_name} dataset:" > ${output_file}
    python3 -u dms_main.py \
        --dataset_name ${dataset_name} \
        --experiment_name "esm2_15B_num_variants" \
        --model_name "esm2_t48_15B_UR50D" \
        --embeddings_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm" \
        --labels_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms" \
        --num_simulations 10 \
        --num_iterations 1 \
        --measured_var "fitness" \
        --learning_strategies "top10" \
        --num_mutants_per_round 10 20 30 40 50 100 \
        --num_final_round_mutants 16 \
        --first_round_strategies "random" \
        --embedding_types "embeddings" \
        --regression_types "randomforest" \
        --embeddings_file_type "csv" \
        --output_dir "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms_results" \
        >> ${output_file} 2>&1
    echo "Done running ${dataset_name} dataset:" >> ${output_file}
}

# Export the function so it's available to GNU Parallel
export -f run_dms_main_100

# Use GNU Parallel to run the dms_main function in parallel for each dataset
parallel -j3 run_dms_main_100 ::: "${datasets[@]}"

