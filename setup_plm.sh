#!/bin/bash

# Print that the script is running:
echo "Setting up PLM environment..."

# Create and activate the conda environment
conda env create --prefix /orcd/archive/abugoot/001/Projects/Matteo/Github/evolvepro_envs/plm -f plm_environment.yml
conda activate plm

# Run the clone script
python evolvepro/plm/clone_git_plms.py

echo "PLM environment setup complete!"