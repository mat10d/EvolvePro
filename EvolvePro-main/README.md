

# User Instruction
There are two anaconda environments (esm2.yml and embedding.yml) to download and three executable files to run in sequences for evolving any protein sequences by putting in the folder a fasta file carrying the WT protein sequences. 
1. Upload protein fasta files to the folder.
2. Run extract_15B.sh by changing the name of the input to the target fasta file and this will generate raw embeddings (ESM2-15B) for the target protein.
3. Run concatenate.sh which automatically aggregates the embeddings for each mutant to a single csv file in the results folder.
4. Run toplayer.sh which takes in the experimental results and predicts the top 10 mutations for subsequent rounds testing.
5. Clone and test the predicted mutations and feed the data back in by appending the next round results to the experimental file by re-running toplayer.sh to generate more mutants.




