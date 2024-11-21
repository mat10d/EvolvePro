# Plot

This directory contains scripts for plotting functions for visualizing EVOLVEpro results. Examples are contained in `dms.py` and `exp.py`, and are based on source functions in `evolvepro/src/evolve.py`.

For plot, use the evolvepro environment:

```bash
conda activate evolvepro
```

## DMS Benchmark Visualization (dms.py)

### Data Processing
- `read_dms_data()`: Loads and processes DMS benchmark results, averaging across simulations of the same strategy
- `filter_dataframe()`: Filters data based on specified conditions

### Performance Plots
- `plot_comparison()`: Bar plots comparing model/strategy performance
- `plot_grid_search_bar()`: Bar plots of hyperparameter performance
- `plot_grid_search_heatmap()`: Heatmaps showing parameter interactions
- `plot_by_round()`: Line plots of evolution progress
- `plot_by_round_split()`: Evolution progress split by conditions

## Experimental Evolution Visualization (exp.py)
- `read_exp_data()`: Loads experimental round data from Excel files
- `plot_variants_by_iteration()`: Tracks activity progression across rounds

All functions save outputs as PNG/PDF files and processed data as CSVs.