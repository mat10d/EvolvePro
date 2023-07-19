# directed_evolution

The order of analysis is specified below. This is a rough repo, but overall describes key steps undertaken in the analysis. We include the Brenan MAPK1 DMS experiment as a reference.

### data_processing:

`data_processing_VEP.ipynb` contains the basic code for processing a file from the Livesey & Marsch Mol Sys Bio manuscript (Using deep mutational scanning to benchmark variant effect predictors and identify disease mutations)

This relies on the `Source.xlsx` file that contains information regarding all the deep mutational scanning experiments, and a fasta file of the WT AA sequence of the protein of interest.

Example outputs are shown--essentially a csv of the tab in the `Source.xlsx` table with some additional variables (`brenan_labels.csv`) and a fasta file (`brenan.fasta`) with all of the subsitutions tested in the paper. There is a specific format for data read in which is compatible with the `Source.xlsx`

### esm-extract:

`extract.sh` is a basic OpenMind compatible bash file for running the `extract.py` file released with esm. It relies on the fasta file to mean embeddings of mutants. The output format is a .pth file for each mutant, where each file is named after the substitution.

To make the data more workable for downstream tasks, `concatenate.sh` calls on `concatenate.py` to generate a single csv of mean esm embeddings. The results of this are saved in `results_means/csvs`

This output is not here due to size constraints, but can be found at: https://www.dropbox.com/scl/fo/kkizj4qt6s00n6zbq9zbs/h?rlkey=dehfdmy7dhbhauqqp0y3fjmnh&dl=0

### simulate

`simulate.py` is the general script for simulating 3 iterations of this approach across a number of parameters including embedding scaling, fitness metric, number of rounds, number of mutants per rounds, top layer type, and esm version.

A very similar script `simulate_rounds.py`exists to do this across 2-n rounds to see how final output metrics change as we add more rounds. Simulate rounds was run with a fix set of other parameters with the exception of number of mutants per round.

The `brenan_[].sh` files are simple OpenMind compatible bash files for running the above simulations.

### top-layer-metrics

This directory contains notebooks that assist with visualizations of simulation outputs.
