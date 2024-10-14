# EVOLVEpro Code Repository

This repository contains the code for EVOLVEpro, a model for in silico directed evolution of protein activities using few-shot active learning.

## Process Overview

The EVOLVEpro workflow consists of four main steps:

1. **Process**: Generate and clean FASTA and CSV files
2. **PLM**: Extract protein language model (PLM) embeddings for all variants
3. **Run EVOLVEpro**: Apply the model to either DMS or experimental data
4. **Visualize**: Prepare outputs and visualizations

## Step-by-Step Description

### 1. Process

Generate and clean FASTA and CSV files containing protein variant sequences and their corresponding activity data.

For detailed instructions, see the [Process README](scripts/process/README.md).

### 2. PLM

Extract protein language model embeddings for all variants using various PLM models.

For detailed instructions, see the [PLM README](scripts/plm/README.md).

### 3. Run EVOLVEpro

Apply the EVOLVEpro model. Choose the appropriate workflow based on your data type:
- Deep mutational scanning (DMS) data: See the [DMS README](scripts/dms/README.md)
- Experimental data: See the [Experimental README](scripts/experimental/README.md)

### 4. Visualize

Prepare outputs and create visualizations to interpret the results of the EVOLVEpro process.

For detailed instructions, see the [Visualize README](scripts/visualize/README.md).

## Getting Started

[Include setup instructions, dependencies, and other relevant information]

## Overview of approach for *de novo* improvement of a protein's activity

## Citation

If you use this code in your research, please cite our paper:
[Include citation information]

## License

[Include license information]