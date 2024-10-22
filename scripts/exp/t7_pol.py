from evolvepro.src.evolve import evolve_experimental, evolve_experimental_multi

protein_name = 't7_pol'
embeddings_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/t7_pol/esm'
embeddings_file_name = 't7_pol_esm2_t48_15B_UR50D.csv'
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/t7_pol/rounds'
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/t7_pol_WT.fasta"
number_of_variants = 12
output_dir = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_results/'

# Round 1
round_name = 'Round1'
round_file_names = ['T7_pol_Round1.xlsx']
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

# Round 2
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

# Round 3
round_name = 'Round3'
round_file_names = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx', 'T7_pol_Round3.xlsx']
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
round_file_names = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx', 'T7_pol_Round3.xlsx', 'T7_pol_Round4.xlsx']
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

# Multi-mutant embeddings
embeddings_file_name_2nd = 't7_pol_2nd_esm2_t48_15B_UR50D.csv'
embeddings_file_name_3rd = 't7_pol_3rd_esm2_t48_15B_UR50D.csv'

# Round 5 -- the first round of multi-mutant evolution does not consider the previous round data
round_name = 'Round5'
round_file_names_single = ['T7_pol_Round1.xlsx', 'T7_pol_Round2.xlsx', 'T7_pol_Round3.xlsx', 'T7_pol_Round4.xlsx']
round_file_names_multi = []
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

# Round 6
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