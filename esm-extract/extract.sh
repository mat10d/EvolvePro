#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=12:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=means
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
#SBATCH --constraint=high-capacity    
#SBATCH --mem=10gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/esm-extract/out/means-%j.out 

source ~/.bashrc
conda activate embeddings

study_names=("brenan")

model_names=("esm1b_t33_650M_UR50S" "esm2_t33_650M_UR50D")
fasta_path="/om/group/abugoot/Projects/Matteo/esm-extract/fasta/"
results_path="/om/group/abugoot/Projects/Matteo/esm-extract/results_means/"

repr_layers=33
toks_per_batch=3000

for model_name in "${model_names[@]}"; do
  for study in "${study_names[@]}"; do
    command="python3 extract.py ${model_name} ${fasta_path}${study}.fasta ${results_path}${study}/${model_name} --repr_layers ${repr_layers} --toks_per_batch ${toks_per_batch} --include mean"
    echo "Running command: ${command}"
    eval "${command}"
  done
done