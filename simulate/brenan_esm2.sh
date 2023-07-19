#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=brenan_esm2
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
#SBATCH --constraint=high-capacity    
#SBATCH --mem=10gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/simulate/out/brenan_esm2-%j.out 

source ~/.bashrc
conda activate embeddings

python3 -u simulate.py brenan_esm2_t33_650M_UR50D > /om/group/abugoot/Projects/Matteo/simulate/out/brenan_esm2-hc.out
