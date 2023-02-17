# analyze pdf seqs that david sent me


# %%

#%%
from collections import defaultdict

# Open the FASTA file for reading
with open('./pdb_human_seqs.fa', 'r') as f:
    # Initialize a dictionary to store sequences and accessions
    seq_dict = defaultdict(list)

    # Initialize a sequence variable to store the current sequence being read
    current_seq = ''

    # Read the file line by line
    for line in f:

        line = line.strip()
        if line.startswith('>'):
            # Extract the accession number from the header line
            accession = line[1:]
            # Add the previous sequence to the dictionary with the accession number as the key
            if current_seq:
                seq_dict[current_seq].append(accession)
                current_seq = ''
        else:
            # Add the sequence to the current sequence variable
            current_seq += line

    # Add the last sequence to the dictionary with the last accession number as the key
    if current_seq:
        seq_dict[current_seq].append(accession)

# Output the results as a table
with open('pdb_human_seqs_unique.tsv', 'w') as ofile:
    ofile.write("Sequence\tAccessions\n")
    for seq, accessions in seq_dict.items():
        ofile.write(f"{seq}\t{', '.join(accessions)}" + "\n")

#%%
# Output the results as a fasta file
with open('pdb_human_seqs_unique.fasta', 'w') as ofile:
    for seq, accessions in seq_dict.items():
        ofile.write('>' + '|'.join(accessions) + '\n')
        ofile.write(seq + '\n')



#%%