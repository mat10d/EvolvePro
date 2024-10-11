# Process

This directory contains scripts to prepare files for use in EvolvePro. There are two main scripts for processing mutation data: one for experimental data and another for deep mutational scanning (DMS) data.

## 1. Experimental Mutation Processing

Source functions in `evolvepro/process/exp_process.py`.

This script shows examples of how to generate single amino acid mutants and n-mutant variants for proteins to improve. The output target is a FASTA file of all single AA substitutions relative to the WT sequence.

### How to Run:

To generate a wild-type FASTA file and create single amino acid mutants:

```python
generate_wt('MNTINIAKNDFSDIELAAIPFNTLADHYGERLAREQLALEHE...', 'output_path', 'dataset_WT.fasta')
generate_single_aa_mutants('output_path/dataset_WT.fasta', 'output_path/dataset.fasta')
```

To generate n-mutant combinations:

```python
generate_n_mutant_combinations('output_path/dataset_WT.fasta', 'beneficial_mutations.xlsx', 3, 'output_path/dataset_3rd.fasta', threshold=1)
```

## 2. DMS Data Processing

Source functions in `evolvepro/process/dms_process.py`.

This script processes deep mutational scanning (DMS) data for various proteins.

### How to Run:

Many of these DMS datasets were taken from Livesey & Marsh (Mol Syst Biol.) and can be accessed through [Figshare](https://figshare.com/articles/dataset/Raw_variant_effect_predictions_and_DMS_data_for_benchmarking_variant_effect_predictors_/12369359/1?file=22798430). A high fitness cutoff was set manually for each of these datasets. The output target is a FASTA file of all single AA substitutions *in the DMS dataset*, as well as a file that contains fitness measurements in a specific format. 

```python
process_dataset(
    file_path='data/dms/fitness/Source.xlsx',
    dataset_name='brenan',
    wt_fasta_path='data/dms/wt_fasta/brenan_WT.fasta',
    fitness_column='DMS_SCH',
    cutoff_value=2.5,
    output_dir='output/dms',
    sheet_name='MAPK1',
    cutoff_rule='greater_than',
    cutoff_percentiles=[90, 95]
)
```

To create visualization plots:

```python
plot_mutations_per_position(processed_df)
plot_histogram_of_readout(processed_df, 'DMS_SCH', 2.5)
```

Adjust the input files and parameters as needed for your specific datasets, such as in `preprocess_cas12f`.