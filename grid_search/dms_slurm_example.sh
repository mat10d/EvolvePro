#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=rounds_esm1b
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=4gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/simulate/out/brenan-%j.out 

source ~/.bashrc
conda activate embeddings

# 1) Run full grid search on esm1b mean embeddings
echo "Running Brenan dataset:" > out/brenan_esm1b-hc.out
python3 -u grid_search.py \
    --dataset_name brenan_esm1b_t33_650M_UR50S \
    --base_path ../esm-extract/results_means \
    --num_simulations 3 \
    --num_iterations 3 5 10 \
    --measured_var fitness fitness_scaled \
    --learning_strategies random top5bottom5 top10 dist \
    --num_mutants_per_round 8 10 16 32 128 \
    --embedding_types embeddings embeddings_norm embeddings_pca \
    --regression_types ridge lasso elasticnet linear neuralnet randomforest gradientboosting \
    >> out/brenan_esm1b-hc.out
echo "Done running Brenan dataset:" >> out/brenan_esm1b-hc.out

# 2) Run full grid search on esm2 mean embeddings
echo "Running Brenan dataset:" > out/brenan_esm2-hc.out
python3 -u grid_search.py \
    --dataset_name brenan_esm2_t33_650M_UR50D \
    --base_path ../esm-extract/results_means \
    --num_simulations 3 \
    --num_iterations 3 5 10 \
    --measured_var fitness fitness_scaled \
    --learning_strategies random top5bottom5 top10 dist \
    --num_mutants_per_round 8 10 16 32 128 \
    --embedding_types embeddings embeddings_norm embeddings_pca \
    --regression_types ridge lasso elasticnet linear neuralnet randomforest gradientboosting \
    >> out/brenan_esm2-hc.out
echo "Done running Brenan dataset:" >> out/brenan_esm2-hc.out

# 3) Run small grid search on esm1b mean embeddings to see how performance changes over rounds
echo "Running Brenan dataset:" > out/brenan_esm1b-hc_small.out
python3 -u grid_search.py \
    --dataset_name brenan_esm1b_t33_650M_UR50S \
    --base_path ../esm-extract/results_means \
    --num_simulations 10 \
    --num_iterations 2 3 4 5 6 7 8 9 10 11 \
    --measured_var fitness \
    --learning_strategies top10 \
    --num_mutants_per_round 16 \
    --embedding_types embeddings \
    --regression_types randomforest \
    >> out/brenan_esm1b-hc_small.out
echo "Done running Brenan dataset:" >> out/brenan_esm1b-hc_small.out

# 4) Run small grid search on esm2 mean embeddings to see how performance changes over rounds
echo "Running Brenan dataset:" > out/brenan_esm2-hc_small.out
python3 -u grid_search.py \
    --dataset_name brenan_esm2_t33_650M_UR50D \
    --base_path ../esm-extract/results_means \
    --num_simulations 10 \
    --num_iterations 2 3 4 5 6 7 8 9 10 11 \
    --measured_var fitness \
    --learning_strategies top10 \
    --num_mutants_per_round 16 \
    --embedding_types embeddings \
    --regression_types randomforest \
    >> out/brenan_esm2-hc_small.out
echo "Done running Brenan dataset:" >> out/brenan_esm2-hc_small.out