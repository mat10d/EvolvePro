#!/bin/bash

# Print that the script is running:
echo "Setting up PLM environment..."

# Create and activate the conda environment
conda env create -f plm_environment.yml
conda activate plm

# Run the clone script
python evolvepro/plm/clone_git_plms.py

echo "PLM environment setup complete!"