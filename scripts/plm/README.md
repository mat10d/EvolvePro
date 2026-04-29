# Protein Language Model (PLM) Embeddings Generation

This directory contains legacy SLURM-compatible files for generating protein language model embeddings from FASTA inputs. They call `evolvepro/plm/[protein_language_model]/extract.py` directly and remain compatibility scaffolding during the CLI migration.

## Target workspace workflow

For active-learning workspaces, prefer the high-level CLI:

```bash
evolvepro init my-protein --set-default --yes
cd my-protein
# edit evolvepro.toml: [embedding] backend/model/device as needed
evolvepro evolve
# fill activity values in rounds/round_000/candidates.csv
evolvepro evolve --device cuda:0
```

When `inputs/embeddings.csv` exists, `evolvepro evolve` reuses it and records raw/semantic hashes. When embeddings are missing, `evolvepro evolve` can prepare them through the configured backend. Generated FASTA, embedding CSV, metadata, backend/device metadata, and CAS snapshots are recorded in `round.lock` and `manifest.json`.

GPU-first defaults apply to PLM extraction:

- `device = "auto"` prefers CUDA.
- Missing CUDA fails unless CPU fallback is explicitly enabled.
- CPU extraction is compatibility/test mode:

```bash
evolvepro evolve --allow-cpu-embedding --device cpu
```

## Legacy direct extraction

Use these scripts only when reproducing old workflows or preparing standalone embedding CSV files outside the workspace state machine. The direct extractor usage varies by model. A general shape is:

```bash
python evolvepro/plm/[model_name]/extract.py [model_location] [fasta_file] [output_dir] [additional_options]
```

Example ESM-2 15B direct extraction:

```bash
python evolvepro/plm/esm/extract.py esm2_t48_15B_UR50D output/dms/brenan.fasta output/plm/esm/brenan --toks_per_batch 512 --include mean --concatenate_dir output/plm/esm/
```

Arguments in this example:

- `esm2_t48_15B_UR50D`: model name
- `output/dms/brenan.fasta`: input FASTA
- `output/plm/esm/brenan`: per-record output directory
- `--toks_per_batch 512`: token batch size; tune for GPU memory
- `--include mean`: mean-pooled embeddings
- `--concatenate_dir output/plm/esm/`: concatenated CSV output directory

Direct extraction output should be normalized to the workspace CSV contract before target provenance:

```text
variant,feature_0,feature_1,...
```

The current target backend bridge covers ESM-style extraction. Other PLM shell scripts remain legacy compatibility entry points until equivalent backend adapters exist.
