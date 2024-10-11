#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=24:00:00 
#SBATCH -x node[110]
#SBATCH --job-name=proteinbert
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=1  
#SBATCH --constraint=high-capacity    
#SBATCH --mem=200gb  
#SBATCH --output out/proteinbert-%j.out 

source ~/.bashrc
conda activate embeddings

cd /orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro

study_names=("brenan" "jones" "stiffler" "haddox" "doud" "giacomelli" "kelsic" "lee" "markin" "cas12f" "cov2_S" "zikv_E")

fasta_path="output/dms/"
results_path="output/plm/proteinbert/"
model_names=("proteinbert")

mkdir -p ${results_path}

for model_name in "${model_names[@]}"; do
  for study in "${study_names[@]}"; do
    command="python evolvepro/plm/proteinbert/extract.py --input ${fasta_path}${study}.fasta --output ${results_path}${study}_${model_name}.csv"
    echo "Running command: ${command}"
    eval "${command}"
  done
done