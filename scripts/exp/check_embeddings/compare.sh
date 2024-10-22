#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=1:00:00 
#SBATCH --job-name=scratch
#SBATCH -n 12
#SBATCH -N 1   
#SBATCH --cpus-per-task=5  
#SBATCH --mem=100gb  
#SBATCH --output scratch-%j.out 

source ~/.bashrc
conda activate embeddings

python compare_fasta.py
python compare_tensor.py

# python compare_csv.py
