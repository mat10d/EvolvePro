# Protein Language Model (PLM) Embeddings Generation

This directory contains SLURM-compatible files to generate protein language model embeddings from various models. These files run `evolvepro/plm/[protein_language_model]/extract.py` to generate embeddings in CSV format from a FASTA file containing all single amino acid mutants of interest.

## How to Run

The usage of `extract.py` files may vary depending on the model. Below is a general example of how to use the scripts:

```bash
python plm/[model_name]/extract.py [model_location] [fasta_file] [output_dir] [additional_options]
```

### Example:

To generate a CSV file of ESM-2 (15B parameter model) embeddings of the last hidden layer from a FASTA file of single amino acid mutants:

```bash
python plm/esm/extract.py esm2_t48_15B_UR50D output/dms/brenan.fasta output/plm/esm/brenan --toks_per_batch 512 --include mean --concatenate_dir output/plm/esm/
```

In this example:
- `esm2_t48_15B_UR50D` is the model name
- `output/dms/brenan.fasta` is the input FASTA file
- `output/plm/esm/brenan` is the output directory for the embeddings
- `--toks_per_batch 512` sets the number of tokens per batch
- `--include mean` specifies to include mean pooled embeddings
- `--concatenate_dir output/plm/esm/` sets the directory for concatenated results

Note: Adjust `--toks_per_batch` based on your GPU memory capacity.

## Model-Specific Instructions

Different protein language models may require specific commands or have unique parameters. Refer to the `extract.py` files in each model's subdirectory and/or the SLURM-compatible files in this directory for model-specific instructions. We recommend using a GPU powered HPC to extract protein language model embeddings, but these SLURM-compatible files can be easily adapted into a command line call.