#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
#SBATCH -x node[110]
#SBATCH --job-name=esm
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=1  
#SBATCH --constraint=high-capacity    
#SBATCH --mem=200gb  
#SBATCH --output out/esm-%j.out 

source ~/.bashrc
conda activate embeddings

cd /orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro

study_names=("brenan" "jones" "stiffler" "haddox" "doud" "giacomelli" "kelsic" "lee" "markin" "cas12f" "cov2_S" "zikv_E")
model_names=("esm2_t48_15B_UR50D")
fasta_path="output/dms/"
results_path="output/plm/esm/"

repr_layers=48
toks_per_batch=512

mkdir -p ${results_path}

for model_name in "${model_names[@]}"; do
  for study in "${study_names[@]}"; do
    command="python3 evolvepro/plm/esm/extract.py ${model_name} ${fasta_path}${study}.fasta ${results_path}${study}/${model_name} --toks_per_batch ${toks_per_batch} --include mean --concatenate_dir ${results_path}"
    echo "Running command: ${command}"
    eval "${command}"
  done
done