#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=4:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=concatenate
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=10gb  

source ~/.bashrc
conda activate embeddings

python3 concatenate.py