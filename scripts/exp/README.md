# Experimental Workflow for EVOLVEpro

This directory contains scripts for running the experimental workflow of EVOLVEpro, which iteratively optimizes protein activity through experimental rounds of evolution. The source functions used here are in `evolvepro/src/evolve.py`, which calls the underlying model.

For exp, use the evolvepro environment:

```bash
conda activate evolvepro
```

## Evolution Types

1. Single Mutant Evolution: Explores individual amino acid substitutions
2. Multi-Mutant Evolution: Explores combinations of mutations based on previous rounds

## Input Files

1. FASTA file: Contains the wild-type protein sequence
2. PLM embeddings: CSV file(s) containing embeddings for all mutants of interest
3. Round data: Excel files containing fitness measurements for each round of evolution

## Key Parameters

- `protein_name`: Name of the protein being evolved
- `round_name`: Identifier for the current round of evolution
- `number_of_variants`: Number of variants to predict for the next round
- `rename_WT`: Boolean to indicate if the wild-type sequence should be renamed in the output

## Usage

Create a Python script (e.g., `t7_pol.py`) with the following structure:

```python
from evolvepro.src.evolve import evolve_experimental, evolve_experimental_multi

protein_name = 't7_pol'
embeddings_base_path = '/path/to/embeddings'
embeddings_file_name = 'embeddings_file.csv'
round_base_path = '/path/to/round/data'
wt_fasta_path = "/path/to/wildtype/fasta"
number_of_variants = 12
output_dir = '/path/to/output/directory'

# Single variant
round_name = 'Round2'
round_file_names = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx']
rename_WT = True

evolve_experimental(
    protein_name,
    round_name,
    embeddings_base_path,
    embeddings_file_name,
    round_base_path,
    round_file_names,
    wt_fasta_path,
    rename_WT,
    number_of_variants,
    output_dir
)

# Multivariant
embeddings_file_name_2nd = 'embeddings_2nd_file.csv'
embeddings_file_name_3rd = 'embeddings_3rd_file.csv'

round_name = 'Round6'
round_file_names_single = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx', 'T7_pol_Round3.xlsx', 'T7_pol_Round4.xlsx']
round_file_names_multi = ['T7_pol_Round5.xlsx']
rename_WT = True

evolve_experimental_multi(
    protein_name,
    round_name,
    embeddings_base_path,
    [embeddings_file_name, embeddings_file_name_2nd, embeddings_file_name_3rd],
    round_base_path,
    round_file_names_single,
    round_file_names_multi,
    wt_fasta_path,
    rename_WT,
    number_of_variants,
    output_dir
)
```

Detailed results will be saved in the specified output directory, for each round, and the specified top `number_of_variants` to assess for the following round will be returned.