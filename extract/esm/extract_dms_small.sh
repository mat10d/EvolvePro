#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=24:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=means
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
#SBATCH --constraint=high-capacity    
#SBATCH --mem=100gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/Github/directed_evolution/extract/esm/out/means-%j.out 

source ~/.bashrc
conda activate embeddings

cd  /om/group/abugoot/Projects/Matteo/Github

study_names=("brenan" "jones" "stiffler" "haddox" "doud" "giacomelli" "kelsic" "lee" "markin" "cas12f" "cov2_S" "zikv_E")

# model_names=("/om/group/abugoot/Projects/Matteo/Github/directed_evolution/.cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt" "/om/group/abugoot/Projects/Matteo/Github/directed_evolution/.cache/torch/hub/checkpoints/esm1v_t33_650M_UR90S_1.pt" "/om/group/abugoot/Projects/Matteo/Github/directed_evolution/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S.pt")
model_names=("esm1v_t33_650M_UR90S_1")

fasta_path="directed_evolution/data_processing/output/"
results_path="directed_evolution/extract/esm/results_means/"

repr_layers=33
toks_per_batch=2000

for model_name in "${model_names[@]}"; do
  for study in "${study_names[@]}"; do
    model_part=$(basename "${model_name}")
    model_part="${model_part%%.*}"  # Removing extension
    
    command="python3 directed_evolution/extract/esm/extract.py ${model_name} ${fasta_path}${study}.fasta ${results_path}${study}/${model_part} --repr_layers ${repr_layers} --toks_per_batch ${toks_per_batch} --include mean"
    echo "Running command: ${command}"
    eval "${command}"
  done
done