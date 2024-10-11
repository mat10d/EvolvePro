#!/bin/bash
# Configuration values for SLURM job submission.
#SBATCH --time=2:00:00 
#SBATCH --job-name=one-hot
#SBATCH -n 1 
#SBATCH -N 1   
#SBATCH --cpus-per-task=1  
#SBATCH --mem=100gb  
#SBATCH --output out/one-hot-%j.out 

source ~/.bashrc
conda activate evolvepro
cd /orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro

study_names=("brenan" "jones" "stiffler" "haddox" "doud" "giacomelli" "kelsic" "lee" "markin" "cas12f" "cov2_S" "zikv_E")
encoding_methods=("one_hot" "integer")
fasta_path="output/dms/"
results_path="output/plm/one-hot/"

for method in "${encoding_methods[@]}"; do
  for study in "${study_names[@]}"; do
    command="python3 evolvepro/plm/one-hot/extract.py ${fasta_path}${study}.fasta --method ${method} --results_path ${results_path}"
    echo "Running command: ${command}"
    eval "${command}"
  done
done