#!/bin/bash

# 1) Run full grid search on esm1b mean embeddings
python3 -u grid_search.py \
    --dataset_name brenan_esm1b_t33_650M_UR50S \
    --base_path ../esm-extract/results_means \
    --num_simulations 3 \
    --num_iterations 3 5 10 \
    --measured_var fitness fitness_scaled \
    --learning_strategies random top5bottom5 top10 dist \
    --num_mutants_per_round 8 10 16 32 128 \
    --embedding_types embeddings embeddings_norm embeddings_pca \
    --regression_types ridge lasso elasticnet linear neuralnet randomforest gradientboosting

# 2) Run full grid search on esm2 mean embeddings
python3 -u grid_search.py \
    --dataset_name brenan_esm2_t33_650M_UR50D \
    --base_path ../esm-extract/results_means \
    --num_simulations 3 \
    --num_iterations 3 5 10 \
    --measured_var fitness fitness_scaled \
    --learning_strategies random top5bottom5 top10 dist \
    --num_mutants_per_round 8 10 16 32 128 \
    --embedding_types embeddings embeddings_norm embeddings_pca \
    --regression_types ridge lasso elasticnet linear neuralnet randomforest gradientboosting

# 3) Run small grid search on esm1b mean embeddings to see how performance changes over rounds
python3 grid_search.py \
    --dataset_name brenan_esm1b_t33_650M_UR50S \
    --base_path ../esm-extract/results_means \
    --num_simulations 10 \
    --num_iterations 2 3 4 5 6 7 8 9 10 11 \
    --measured_var fitness \
    --learning_strategies top10 \
    --num_mutants_per_round 16 \
    --embedding_types embeddings \
    --regression_types randomforest

# 4) Run small grid search on esm2 mean embeddings to see how performance changes over rounds
python3 grid_search.py \
    --dataset_name brenan_esm2_t33_650M_UR50D \
    --base_path ../esm-extract/results_means \
    --num_simulations 10 \
    --num_iterations 2 3 4 \
    --measured_var fitness \
    --learning_strategies top10 \
    --num_mutants_per_round 16 \
    --embedding_types embeddings \
    --regression_types randomforest