# EVOLVEpro
This is the offical codebase to reproduce and use EVOLVEpro

# Execution
For using EVOLVEpro to optimize any protein of interest, download the files in the EVOLVEpro-Main folder and follow the instructions below. 
There are two anaconda environments (esm2.yml and embedding.yml) to download and three executable files to run in sequences for evolving any protein sequences by putting in the folder a fasta file carrying the WT protein sequences. 
1. Upload protein fasta files to the folder.
2. Run extract_15B.sh by changing the name of the input to the target fasta file and this will generate raw embeddings (ESM2-15B) for the target protein.
3. Run concatenate.sh which automatically aggregates the embeddings for each mutant to a single csv file in the results folder.
4. Run toplayer.sh which takes in the experimental results and predicts the top 10 mutations for subsequent rounds testing.
5. Clone and test the predicted mutations and feed the data back in by appending the next round results to the experimental file by re-running toplayer.sh to generate more mutants.

# Further finetuning of the model and reproducing the data in Jiang et al 2024

The order of analysis is specified below. This is a rough repo, but overall describes key steps undertaken in the analysis. We include the Brenan MAPK1 DMS experiment as a reference.

### data_processing:

`data_processing_VEP.ipynb` contains the basic code for processing a file from the Livesey & Marsch Mol Sys Bio manuscript (Using deep mutational scanning to benchmark variant effect predictors and identify disease mutations)

This relies on the `Source.xlsx` file that contains information regarding all the deep mutational scanning experiments, and a fasta file of the WT AA sequence of the protein of interest (`brenan_WT.fasta`).

Example outputs are shown--essentially a csv of the tab in the `Source.xlsx` table with some additional variables (`brenan_labels.csv`) and a fasta file (`brenan.fasta`) with all of the subsitutions tested in the paper. There is a specific format for data read in which is compatible with the `Source.xlsx`

### esm-extract:

`extract.sh` is a basic OpenMind compatible bash file for running the `extract.py` file released with esm. It relies on the fasta file to mean embeddings of mutants. The output format is a .pth file for each mutant, where each file is named after the substitution.

To make the data more workable for downstream tasks, `concatenate.sh` calls on `concatenate.py` to generate a single csv of mean esm embeddings. The results of this are saved in `results_means/csvs`

This output is not here due to size constraints, but can be found at: https://www.dropbox.com/scl/fo/kkizj4qt6s00n6zbq9zbs/h?rlkey=dehfdmy7dhbhauqqp0y3fjmnh&dl=0

### grid_search

`grid_search.py` is the general script for doing ML assisted targeted evolution approaches across a number of parameters in a full grid search including embedding scaling, fitness metric, number of rounds, number of mutants per rounds, top layer type, and esm version:

* dataset_name: Name of the esm embeddings csv file of the dataset you want to use for the simulations.
* base_path: Base path of the dataset, which should contain the necessary ESM embeddings and labels CSV files.
* num_simulations: Number of simulations for each parameter combination.
* num_iterations: List of integers representing the number of iterations for the simulations. Must be greater than 1.
* measured_var: List of strings indicating the fitness type to train on. Choose from: fitness, fitness_scaled.
* learning_strategies: List of strings representing the type of learning strategy to use. Choose from: random, top5bottom5, top10, dist.
* num_mutants_per_round: List of integers representing the number of mutants per round. Must be a positive integer.
* embedding_types: List of strings representing the types of embeddings to train on. Choose from: embeddings, embeddings_pca.
* regression_types: List of strings representing the regression types for the ML models. Choose from: ridge, lasso, elasticnet, linear, neuralnet, randomforest, gradientboosting.

The `brenan_[].sh` files are simple OpenMind/command line compatible bash files for running the above simulations, with examples for a large scale grid search (1,2), and surveying improvement across rounds (3,4).

### top-layer-metrics

This directory contains notebooks that assist with visualizations of simulation outputs. Still rough/needs to be improved.

### top-layer-metrics

A yml file that temporarily contains all packages required to work with this repo is also provided.
