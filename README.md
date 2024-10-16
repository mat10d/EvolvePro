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

Apply the EVOLVEpro model to optimize protein activity. There are two main workflows:

#### DMS Workflow
Use this workflow to optimize a few-shot model on a deep mutational scanning (DMS) dataset, where fitness values are known for a large number of variants.

For detailed instructions, see the [DMS README](scripts/dms/README.md).

#### Experimental Workflow
Use this workflow for iterative experimental optimization of protein activity.

For detailed instructions, see the [Experimental README](scripts/exp/README.md).

### 4. Visualize

Prepare outputs and create visualizations to interpret the results of the EVOLVEpro process.

For detailed instructions, see the [Plot README](scripts/plot/README.md).

## Getting Started

[Include setup instructions, dependencies, and other relevant information]

## Colab Tutorial

For a step-by-step guide on using EVOLVEpro to improve a protein's activity, see our Google Colab tutorial:

[Include link to Google Colab notebook here]

## Citation

If you use this code in your research, please cite our paper:
[Include citation information]

## License

[Include license information]