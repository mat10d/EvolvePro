#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=esm2_15B_full_grid_search
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=20gb  
#SBATCH --output /orcd/archive/abugoot/001/Projects/Matteo/Github/directed_evolution/grid_search/out/esm2_15B_full_grid_search-%j.out 

source ~/.bashrc
conda activate embeddings

names=("zikv_E" "cas12f")
# names=("brenan" "stiffler" "doud" "haddox" "giacomelli" "jones" "kelsic" "lee" "markin" "zikv_E" "cas12f" "cov2_S")

datasets=()

for name in "${names[@]}"; do
    datasets+=("${name}_esm2_t48_15B_UR50D")
done

experiment_name="esm2_15B_full_grid_search"
num_simulations=3
measured_var=("fitness" "fitness_scaled")

echo ${measured_var[@]}
echo ${measured_var[*]}

echo "${measured_var[@]}"
echo "${measured_var[*]}"

export experiment_name="esm2_15B_full_grid_search"
export num_simulations=3
export measured_var=("fitness" "fitness_scaled")

echo ${measured_var[@]}
echo ${measured_var[*]}

echo "${measured_var[@]}"
echo "${measured_var[*]}"

