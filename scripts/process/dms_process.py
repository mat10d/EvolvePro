import os
import pandas as pd
from evolvepro.src.process import process_dataset, plot_mutations_per_position, plot_histogram_of_readout, preprocess_cas12f, preprocess_zikv_E, preprocess_cov2_S, markin_custom_cutoff

# Find the project root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))

# Brenan
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'brenan'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'brenan_WT.fasta')
activity_column = 'DMS_SCH'
cutoff_value = 2.5
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'MAPK1'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
brenan_df, brenan_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(brenan_df)
plot_histogram_of_readout(brenan_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Giacomelli
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'giacomelli'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'giacomelli_WT.fasta')
activity_column = 'DMS_null_etoposide'
cutoff_value = 1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'P53'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
giacomelli_df, giacomelli_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(giacomelli_df)
plot_histogram_of_readout(giacomelli_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Jones
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'jones'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'jones_WT.fasta')
activity_column = 'DMS_0.625'
cutoff_value = 2.8
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'ADRB2'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
jones_df, jones_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(jones_df)
plot_histogram_of_readout(jones_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Kelsic
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'kelsic'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'kelsic_WT.fasta')
activity_column = 'DMS_rich'
cutoff_value = 0.98
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'infA'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
kelsic_df, kelsic_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(kelsic_df)
plot_histogram_of_readout(kelsic_df, activity_column, cutoff_value)

# Stiffler
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'stiffler'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'stiffler_WT.fasta')
activity_column = 'DMS_amp_2500_(b)'
cutoff_value = 0.01
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'bla'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
stiffler_df, stiffler_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(stiffler_df)
plot_histogram_of_readout(stiffler_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Haddox
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'haddox'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'haddox_WT.fasta')
activity_column = 'DMS'
cutoff_value = 0.1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'env'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
haddox_df, haddox_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(haddox_df)
plot_histogram_of_readout(haddox_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Doud
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'doud'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'doud_WT.fasta')
activity_column = 'DMS'
cutoff_value = 0.1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'HA-H1N1'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
doud_df, doud_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(doud_df)
plot_histogram_of_readout(doud_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Lee
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'Source.xlsx')
dataset_name = 'lee'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'lee_WT.fasta')
activity_column = 'DMS'
cutoff_value = 0.1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'HA-H3N2'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
lee_df, lee_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(lee_df)
plot_histogram_of_readout(lee_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Markin
file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'abf8761_markin_data-s1.csv')
dataset_name = 'markin'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'markin_WT.fasta')
activity_column = 'kcatOverKM_cMUP_M-1s-1'
cutoff_value = 1400000
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = None
cutoff_rule = 'custom'
cutoff_percentiles = [90, 95]
cutoff_function = markin_custom_cutoff
AA_shift = 20
drop_columns = True

# Process the dataset
markin_df, markin_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    cutoff_function=cutoff_function,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(markin_df)
plot_histogram_of_readout(markin_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Cas12f
input_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'DMS_AsCas12f.xlsx')
wt_fasta_output = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'cas12f_WT.fasta')
preprocessed_output_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'DMS_AsCas12f_preprocessed.xlsx')

# Preprocess the dataset
preprocess_cas12f(input_file, wt_fasta_output, preprocessed_output_file)

file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'DMS_AsCas12f_preprocessed.xlsx')
dataset_name = 'cas12f'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'cas12f_WT.fasta')
activity_column = 'avg_activity'
cutoff_value = 1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'Sheet1'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
cas12f_df, cas12f_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(cas12f_df)
plot_histogram_of_readout(cas12f_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Zikv_E
input_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'jvi.01291-19-sd003.xlsx')
wt_fasta_output = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'zikv_E_WT.fasta')
preprocessed_output_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'zikv_E_preprocessed.xlsx')

# Preprocess the dataset
preprocess_zikv_E(input_file, wt_fasta_output, preprocessed_output_file)

file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'zikv_E_preprocessed.xlsx')
dataset_name = 'zikv_E'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'zikv_E_WT.fasta')
activity_column = 'effect'
cutoff_value = 1
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = 'Sheet1'
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
zikv_E_df, zikv_E_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(zikv_E_df)
plot_histogram_of_readout(zikv_E_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Cov2_S
input_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'all_samples_raw_data--sarscov2.csv')
preprocessed_output_file = os.path.join(project_root, 'data' , 'dms', 'activity', 'sarscov2_preprocessed.csv')

# Preprocess the dataset (WT fasta file is not used in this case)
preprocess_cov2_S(input_file, preprocessed_output_file)

file_path = os.path.join(project_root, 'data' , 'dms', 'activity', 'sarscov2_preprocessed.csv')
dataset_name = 'cov2_S'
wt_fasta_path = os.path.join(project_root, 'data' , 'dms', 'wt_fasta', 'cov2_S_WT.fasta')
activity_column = 'mut_escape'
cutoff_value = 0.05
output_dir = os.path.join(project_root, 'output' , 'dms')
sheet_name = None
cutoff_rule = 'greater_than'
cutoff_percentiles = [90, 95]
AA_shift = None
drop_columns = True

# Process the dataset
cov2_S_df, cov2_S_frac = process_dataset(
    file_path=file_path,
    dataset_name=dataset_name,
    wt_fasta_path=wt_fasta_path,
    activity_column=activity_column,
    cutoff_value=cutoff_value,
    output_dir=output_dir,
    sheet_name=sheet_name,
    cutoff_rule=cutoff_rule,
    cutoff_percentiles=cutoff_percentiles,
    AA_shift=AA_shift,
    drop_columns=drop_columns
)

# Print results
plot_mutations_per_position(cov2_S_df)
plot_histogram_of_readout(cov2_S_df, activity_column, cutoff_value)
print(f"Processing complete for dataset: {dataset_name}")

# Save background activity fractions
background_df = pd.DataFrame({
    'brenan': brenan_frac,
    'giacomelli': giacomelli_frac,
    'jones': jones_frac,
    'kelsic': kelsic_frac,
    'stiffler': stiffler_frac,
    'haddox': haddox_frac,
    'doud': doud_frac,
    'lee': lee_frac,
    'markin': markin_frac,
    'cas12f': cas12f_frac,
    'zikv_E': zikv_E_frac,
    'cov2_S': cov2_S_frac
}, index=['defined', '90th percentile', '95th percentile'])

background_df = background_df.T
background_df.to_csv(os.path.join(project_root, 'output' , 'dms', 'background_df.csv'), index_label='dataset')