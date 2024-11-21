import pandas as pd
import os
from evolvepro.src.plot import read_exp_data, plot_variants_by_iteration
import matplotlib.pyplot as plt

output_dir = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_plots/'

# T7_pol
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/t7_pol/rounds'
round_file_names_single = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx', 'T7_pol_Round3.xlsx', 'T7_pol_Round4.xlsx']
round_file_names_multi = ['T7_pol_Round5.xlsx']
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/t7_pol_WT.fasta"

df = read_exp_data(round_base_path, round_file_names_single, wt_fasta_path, round_file_names_multi)

plot_variants_by_iteration(df, activity_column='activity', output_dir=output_dir, output_file="t7_pol")

# mlv
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/mlv/rounds'
round_file_names_single = ['mlv_Round4_all.xlsx', 'mlv_Round5.xlsx', 'mlv_Round6.xlsx']
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/mlv_WT.fasta"

df = read_exp_data(round_base_path, round_file_names_single, wt_fasta_path)

plot_variants_by_iteration(df, activity_column='activity', output_dir=output_dir, output_file="mlv")

# bxb1
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/bxb1/rounds'
round_file_names_single = ['bxb1_Round1.xlsx', 'bxb1_Round4.xlsx', 'bxb1_Round5.xlsx', 'bxb1_Round6.xlsx', 'bxb1_Round7.xlsx', 'bxb1_Round8.xlsx']
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/bxb1_WT.fasta"

df = read_exp_data(round_base_path, round_file_names_single, wt_fasta_path)

plot_variants_by_iteration(df, activity_column='activity', output_dir=output_dir, output_file="bxb1")


