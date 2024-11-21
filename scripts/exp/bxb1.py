from evolvepro.src.evolve import evolve_experimental

protein_name = 'bxb1'
embeddings_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/bxb1/esm'
embeddings_file_name = 'bxb1_esm2_t48_15B_UR50D.csv'
round_base_path = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/exp_data/bxb1/rounds'
wt_fasta_path = "/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/data/exp/wt_fasta/bxb1_WT.fasta"
number_of_variants = 12
output_dir = '/orcd/archive/abugoot/001/Projects/Matteo/Github/EvolvePro/output/exp_results/'

# Round 1
round_name = 'Round1'
round_file_names = ['bxb1_Round1.xlsx']
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
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round2.xlsx']
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
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round2.xlsx', 'bxb1_Round3.xlsx']
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

# Round 4 -- round 3 was excluded due to concerns surrounding activity measurements
round_name = 'Round4'
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round2.xlsx', 'bxb1_Round4.xlsx']
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

# Round 5 -- rounds 2 and 3 were excluded due to concerns surrounding activity measurements (for all subsequent rounds)
round_name = 'Round5'
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round4.xlsx', 'bxb1_Round5.xlsx']
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
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round4.xlsx', 'bxb1_Round5.xlsx', 'bxb1_Round6.xlsx']
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

# Round 7
round_name = 'Round7'
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round4.xlsx', 'bxb1_Round5.xlsx', 'bxb1_Round6.xlsx', 'bxb1_Round7.xlsx']
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

# Round 8
round_name = 'Round8'
round_file_names = ['bxb1_Round1.xlsx', 'bxb1_Round4.xlsx', 'bxb1_Round5.xlsx', 'bxb1_Round6.xlsx', 'bxb1_Round7.xlsx', 'bxb1_Round8.xlsx']
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