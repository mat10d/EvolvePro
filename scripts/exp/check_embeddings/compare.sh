#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=12:00:00 
#SBATCH --job-name=compare
#SBATCH -n 1
#SBATCH -N 1   
#SBATCH --cpus-per-task=12  
#SBATCH --mem=100gb  
#SBATCH --output out/compare-%j.out 

source ~/.bashrc
conda activate embeddings

# python compare_fastas.py
python compare_embeddings.py
