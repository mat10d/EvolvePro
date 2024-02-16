#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=18:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=concatenate_new
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=10gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/Github/directed_evolution/esm-extract/out/concatenate_new-%j.out 

source ~/.bashrc
conda activate embeddings

python3 concatenate.py