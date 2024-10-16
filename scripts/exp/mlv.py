from evolvepro.src.evolve import evolve_experimental

protein_name = 'mlv'
embeddings_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/mlv/esm'
embeddings_file_name = 'mlv_esm2_t48_15B_UR50D.csv'
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/mlv/rounds'
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/mlv_WT.fasta"
number_of_variants = 12
output_dir = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_results/'

# Round 1
round_name = 'Round1'
round_file_names = ['mlv_Round1.xlsx']
rename_WT = False

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

# Round 2
round_name = 'Round2'
round_file_names = ['mlv_Round1.xlsx', 'mlv_Round2.xlsx']
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

# Round 3
round_name = 'Round3'
round_file_names = ['mlv_Round1.xlsx', 'mlv_Round2.xlsx', 'mlv_Round3.xlsx']
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

# Round 4
round_name = 'Round4'
round_file_names = ['mlv_Round4_all.xlsx']
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

# Round 5
round_name = 'Round5'
round_file_names = ['mlv_Round4_all.xlsx', 'mlv_Round5.xlsx']
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

# Round 6
round_name = 'Round6'
round_file_names = ['mlv_Round4_all.xlsx', 'mlv_Round5.xlsx', 'mlv_Round6.xlsx']
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

