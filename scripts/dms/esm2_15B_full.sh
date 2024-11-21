#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=7-00:00:00 
#SBATCH -p abugoot
#SBATCH --job-name=esm2_15B_full
#SBATCH -n 12
#SBATCH -N 1   
#SBATCH --cpus-per-task=16
#SBATCH --mem=80gb  
#SBATCH --output out/esm2_15B_full-%j.out 

source ~/.bashrc
conda activate evolvepro
module load openmind8/gnu-parallel/20240222

datasets=("brenan" "stiffler" "doud" "haddox" "giacomelli" "jones" "kelsic" "lee" "markin" "zikv_E" "cas12f" "cov2_S")

# Function to run dms_main for a given dataset
run_dms_main() {
    dataset_name=$1
    output_file="out/${dataset_name}-esm2_15B_full.out"

    echo "Running ${dataset_name} dataset:" > ${output_file}
    python3 -u dms_main.py \
        --dataset_name ${dataset_name} \
        --experiment_name "esm2_15B_full" \
        --model_name "esm2_t48_15B_UR50D" \
        --embeddings_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/plm/esm" \
        --labels_path "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms" \
        --num_simulations 3 \
        --num_iterations 10 \
        --measured_var "activity" "activity_scaled" \
        --learning_strategies "random" "topn2bottomn2" "topn" "dist" \
        --num_mutants_per_round 16 \
        --num_final_round_mutants 16 \
        --first_round_strategies "random" "diverse_medoids" \
        --embedding_types "embeddings" "embeddings_pca_10" \
        --pca_components 10 \
        --regression_types "ridge" "lasso" "elasticnet" "linear" "neuralnet" "randomforest" "gradientboosting" "knn" "gp" \
        --embeddings_file_type "csv" \
        --output_dir "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/dms_results" \
        >> ${output_file} 2>&1
    echo "Done running ${dataset_name} dataset:" >> ${output_file}
}

# Export the function so it's available to GNU Parallel
export -f run_dms_main

# Use GNU Parallel to run the dms_main function in parallel for each dataset
parallel -j12 run_dms_main ::: "${datasets[@]}"