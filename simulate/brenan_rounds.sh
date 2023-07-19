#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --time=48:00:00 
##SBATCH -x node[110]
#SBATCH --job-name=rounds_esm1b
#SBATCH -n 1 
#SBATCH -N 1   
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1  
##SBATCH --constraint=high-capacity    
#SBATCH --mem=4gb  
#SBATCH --output /om/group/abugoot/Projects/Matteo/simulate/out/rounds_esm1b-%j.out 

source ~/.bashrc
conda activate embeddings

echo "Running Brenan dataset:" > /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out
python3 -u simulate_rounds.py brenan_esm1b_t33_650M_UR50S >> /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out
echo "Done running Brenan dataset:" >> /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out

echo "Running Brenan dataset:" >> /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out
python3 -u simulate_rounds.py brenan_esm2_t33_650M_UR50D >> /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out
echo "Done running Brenan dataset:" >> /om/group/abugoot/Projects/Matteo/simulate/out/brenan_rounds-hc.out