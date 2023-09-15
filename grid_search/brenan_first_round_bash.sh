#!/bin/bash

# 1) Run full grid search on esm1b mean embeddings
echo "Running Brenan dataset:" > out/brenan_FR_esm2-hc.out
python3 -u grid_search.py \
    --dataset_name brenan_esm2_t33_650M_UR50D \
    --base_path ../esm-extract/results_means \
    --num_simulations 3 \
    --num_iterations 2 3 \
    --measured_var fitness \
    --learning_strategies top10 \
    --num_mutants_per_round 16 \
    --first_round_strategies random diverse_medoids representative_hie \
    --embedding_types embeddings \
    --regression_types randomforest \
    --file_type csvs \
    >> out/brenan_FR_esm2-hc.out
echo "Done running Brenan dataset:" >> out/brenan_FR_esm2-hc.out
